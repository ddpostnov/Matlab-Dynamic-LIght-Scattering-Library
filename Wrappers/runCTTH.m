%runCTTH  Bolus transit landmarks (baseline, arrival, upslope, peak, recirculation)
%
%   runCTTH(s,fNames) loads every *b_I_r.mat file in fNames (bolus-intensity
%   datasets that store traces in the same SOURCE/RESULTS structures as the
%   contrast "_K_" files), and—per spatial ROI and, optionally, per pixel—
%   derives model-free transit-time landmarks from each intensity trace:
%
%        tBaselineB   end of the common pre-bolus period (FOV reference)
%        t0B          dye arrival at the ROI/pixel
%        tUpslopeB    moment of steepest sustained rise
%        tPeakB       first-pass concentration maximum
%        tReCircB     end of first-pass washout / onset of recirculation
%
%   No curve fitting is performed. Two artefacts are handled explicitly:
%   (i) motion spikes that would corrupt the steepest-derivative estimate,
%   and (ii) a "bump" on the rise (e.g. partial bolus failure) that would be
%   mistaken for the peak. The intensity "values" reported alongside each
%   landmark time are read from the filtered trace. All window sizes live in s.
%
%   DETECTION (per trace)
%     1. r  = hampel(raw,hampWin,hampSig)          % remove motion spikes
%        y  = sgolayfilt( medfilt1(r,medTime), sgOrder , sgFrame )
%     2. dy = gradient(y)                           % fine, per-frame derivative
%        dyc= movmean(dy,slopeWin)                  % sustained slope (spike-robust)
%     3. tPeak    = first maximum of y with prominence > promFrac*(peak-baseline)
%                   ( skips shoulders / partial-bolus bumps; minor spikes too )
%     4. tUpslope = argmax(dyc) restricted to BEFORE tPeak (leading edge only),
%                   so neither a post-peak recovery nor a transient spike wins.
%     5. t0  : scanning from tUpslope back toward baseline, the first of
%              (a) first local minimum of y,
%              (b) first sample below the baseline median,
%              (c) start of the contiguous run with fine dy > minStep
%              ( minStep = 1 = smallest change resolvable in uint16 ).
%        tBaseline = min(t0) over the FOV minus one frame.
%     6. tReCirc  : earliest of the first prominent local minimum after tPeak,
%                   the first sample where the sustained decline slows below
%                   minStep, or the end of the observation.
%
%   INPUTS
%     s        parameter struct with fields
%                • calcData  "regions" | "all"   ("all" adds per-pixel maps)
%                • medSpace  2-D median window for per-pixel despeckling   (def 3)
%                • medTime   temporal median window                        (def 3)
%                • sgOrder   Savitzky-Golay polynomial order               (def 2)
%                • sgFrame   Savitzky-Golay frame length (odd)             (def 11)
%                • hampWin   Hampel half-window (0 disables despiking)     (def 3)
%                • hampSig   Hampel outlier threshold (n*MAD)              (def 3)
%                • slopeWin  window for the sustained-slope derivative     (def 9)
%                • promFrac  peak/trough prominence as fraction of amplitude(def 0.2)
%                • minStep   smallest detectable change (uint16)           (def 1)
%                • baseWin   baseline window (frames) for the median   (def 10% of frames)
%                • parforCTTHPixels  logical - run the per-pixel loop  (def true)
%                            in parallel.  A WORKER BOUND on the parfor, not a
%                            branch: false runs the identical loop body serially in
%                            the client and starts no pool.  NOTE the parfor is over
%                            image COLUMNS inside a serial loop over ROWS, so a
%                            parallel run enters the pool once per row.
%     fNames   cell array of *b_I_r.mat paths - the RESULTS member of each bolus
%              product.  The SOURCE cube and the SETTINGS are named from it.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUT FILES (side-effects, SOURCE left untouched)
%       *b_I_r.mat   RESULTS  – sMetrics/dvsMetrics extended with t*B/v*B
%                               columns and (calcData="all") imgT*B/imgV*B maps
%       *b_I_s.mat   SETTINGS – field settings.ctthCalculation added
%
%   EXAMPLE
%     s.calcData="all";                 % per-region + per-pixel
%     s.medSpace=3; s.medTime=3;         % despeckle / despike
%     s.sgOrder=2;  s.sgFrame=11;        % smoothing
%     s.hampWin=3;  s.hampSig=3;         % motion-spike rejection
%     s.slopeWin=9; s.promFrac=0.2;      % robust upslope / peak
%     s.minStep=1;                       % uint16 resolution
%     s.parforCTTHPixels=true;           % false: no parallel pool (slower per-pixel pass)
%     files=dir(fullfile(dataRoot,'*b_I_r.mat'));
%     runCTTH(s,fullfile({files.folder}',{files.name}'));
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

function runCTTH(s,fNames)

if ~all( cellfun(@(f) isempty(f) || contains(f,'b_I_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "b_I_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_b_I_r.mat'').']);
end

% ---- defaults (every window size lives in s, nothing hardcoded below) ----
if ~isfield(s,'calcData'); s.calcData="regions"; end   % "regions" | "all"
if ~isfield(s,'medSpace'); s.medSpace=3;  end          % 2-D median window (pixels)
if ~isfield(s,'medTime');  s.medTime =3;  end          % temporal median window
if ~isfield(s,'sgOrder');  s.sgOrder =2;  end          % Savitzky-Golay order
if ~isfield(s,'sgFrame');  s.sgFrame =11; end          % Savitzky-Golay frame
if ~isfield(s,'hampWin');  s.hampWin =3;  end          % Hampel half-window (0=off)
if ~isfield(s,'hampSig');  s.hampSig =3;  end          % Hampel n*MAD threshold
if ~isfield(s,'slopeWin'); s.slopeWin=9;  end          % sustained-slope window
if ~isfield(s,'promFrac'); s.promFrac=0.2;end          % prominence / amplitude
if ~isfield(s,'minStep');  s.minStep =1;  end          % min detectable change (uint16)
% Parallelism is optional: a BOUND on the per-pixel parfor (Inf workers or 0), never
% a branch - parfor(...,0) runs the identical body serially in the client and starts
% no pool.  Default true, so a settings file carrying no such field is unaffected.
if ~isfield(s,'parforCTTHPixels') || isempty(s.parforCTTHPixels)
    s.parforCTTHPixels=true;
end
s.sgFrame=s.sgFrame+1-mod(s.sgFrame,2);                % force odd
s.sgFrame=max(s.sgFrame,s.sgOrder+3-mod(s.sgOrder,2)); % keep frame > order, odd

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  Cancel is only
% checked between files (a hook inside the per-pixel parfor would broadcast oddly).
rep=reportOpen(s,'CTTH',fNames);

for fidx=1:numel(fNames)
     if reportCancelled(rep), break; end        % cooperative cancel between files
     if ~isempty(fNames{fidx})

    s.fName=fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results source settings
    load(getProductPath(s.fName,'d'))
    load(getProductPath(s.fName,'s'));
    load(s.fName);

    dt=results.time(2)-results.time(1);
    if ~isfield(s,'baseWin') || isempty(s.baseWin)
        bw=max(5,round(0.1*size(results.sData,1)));
    else
        bw=s.baseWin;
    end

    % ----------------------------- regions -----------------------------
    results.sMetrics=appendBolusMetrics(results.sMetrics,results.sData,results.time,s,bw,dt);

    % ------------------------- vessel segments -------------------------
    if isfield(results,"dvsData")
        results.dvsMetrics=appendBolusMetrics(results.dvsMetrics,results.dvsData,results.time,s,bw,dt);
    end

    % --------------------------- pixels (opt) --------------------------
    if strcmp(s.calcData,"all")
        data=single(source.data);
        for k=1:size(data,3)                         % 2-D median across frames
            data(:,:,k)=medfilt2(data(:,:,k),[s.medSpace s.medSpace],'symmetric');
        end
        mask=results.cMask>0;
        [R,C,~]=size(data);
        t=source.time(:);
        imgT=nan(R,C,4,'single');
        imgV=nan(R,C,4,'single');
        nwPix=0; if s.parforCTTHPixels, nwPix=Inf; end   % worker bound, not a branch
        for i=1:R
            rowData=squeeze(data(i,:,:));            % C x T
            rowMask=mask(i,:);
            rowT=nan(C,4); rowV=nan(C,4);
            parfor (j=1:C, nwPix)
                if rowMask(j)
                    [idx,y]=getBolusTimes(rowData(j,:),s,bw);
                    if all(~isnan(idx))
                        rowT(j,:)=reshape(t(idx),1,4);
                        rowV(j,:)=reshape(y(idx),1,4);
                    end
                end
            end
            imgT(i,:,:)=reshape(rowT,1,C,4);
            imgV(i,:,:)=reshape(rowV,1,C,4);
        end
        results.imgT0B       =imgT(:,:,1);  results.imgV0B       =imgV(:,:,1);
        results.imgTUpslopeB =imgT(:,:,2);  results.imgVUpslopeB =imgV(:,:,2);
        results.imgTPeakB    =imgT(:,:,3);  results.imgVPeakB    =imgV(:,:,3);
        results.imgTReCircB  =imgT(:,:,4);  results.imgVReCircB  =imgV(:,:,4);
        results.imgTBaselineB=min(imgT(:,:,1),[],'all','omitnan')-dt;
    end

    settings.ctthCalculation=reportSettings(s);
    reportWriting(rep);
    save(s.fName,'results','-v7.3','-nocompression');
    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
    reportSaved(rep);
end
end
reportClose(rep);
end


% ====================================================================== %
function M=appendBolusMetrics(M,data,time,s,bw,dt)
%appendBolusMetrics  landmark times/values for a [frames x ROI] matrix
N=size(data,2);
idxAll=nan(N,4);
valAll=nan(N,4);
for ii=1:N
    [idx,y]=getBolusTimes(data(:,ii),s,bw);
    if all(~isnan(idx))
        idxAll(ii,:)=idx;
        valAll(ii,:)=reshape(y(idx),1,4);
    end
end
tvec=time(:);
tt=nan(N,4);
ok=all(~isnan(idxAll),2);
tt(ok,:)=tvec(idxAll(ok,:));
tBaseline=min(tt(:,1),[],'omitnan')-dt;

M.('t0B')       =tt(:,1);
M.('tUpslopeB') =tt(:,2);
M.('tPeakB')    =tt(:,3);
M.('tReCircB')  =tt(:,4);
M.('tBaselineB')=tBaseline.*ones(height(M),1);
M.('v0B')       =valAll(:,1);
M.('vUpslopeB') =valAll(:,2);
M.('vPeakB')    =valAll(:,3);
M.('vReCircB')  =valAll(:,4);
end


% ====================================================================== %
function [idx,y]=getBolusTimes(r,s,bw)
%getBolusTimes  [t0 tUpslope tPeak tReCirc] indices + filtered trace y
r=double(r(:));
T=numel(r);
idx=nan(1,4);
y=nan(T,1);
if T<s.sgFrame || any(isnan(r)) || any(isinf(r))
    return
end

% despike (motion artefacts), then median + Savitzky-Golay
if s.hampWin>0
    r=hampel(r,s.hampWin,s.hampSig);
end
y  =sgolayfilt(medfilt1(r,s.medTime,'truncate'),s.sgOrder,s.sgFrame);
dy =gradient(y);                             % fine, per-frame derivative
dyc=movmean(dy,s.slopeWin);                  % sustained slope (spike-robust)

% amplitude / prominence reference
baseMed=median(y(1:min(max(bw,1),T)));
prom=max(s.promFrac*(max(y)-baseMed),s.minStep);

% (3) tPeak : first sufficiently-prominent maximum (skips bumps & spikes)
[~,locP]=findpeaks(y,'MinPeakProminence',prom);
if isempty(locP); [~,iPk]=max(y); else; iPk=locP(1); end

% (4) tUpslope : steepest sustained rise on the leading edge (before the peak)
[~,iUp]=max(dyc(1:max(iPk-1,1)));

% (5) onset t0, searching from the upslope back toward baseline
i0=[];
for k=iUp-1:-1:2                             % (a) first local minimum
    if y(k)<=y(k-1) && y(k)<=y(k+1); i0=k; break; end
end
if isempty(i0)                               % (b) first sample below baseline
    i0=find(y(1:iUp)<baseMed,1,'last');
end
if isempty(i0)                               % (c) start of supra-threshold rise
    k=iUp;
    while k>1 && dy(k-1)>s.minStep; k=k-1; end
    i0=k;
end
if isempty(i0); i0=1; end

% (6) tReCirc : earliest of prominent trough / sustained decline < minStep / end
iMin=[];
if iPk<T-1
    [~,locM]=findpeaks(-y(iPk:T),'MinPeakProminence',prom);
    if ~isempty(locM); iMin=iPk+locM(1)-1; end
end
iDec=find(dyc(iPk+1:T)>-s.minStep,1,'first');
if ~isempty(iDec); iDec=iDec+iPk; end
iRc=min([iMin,iDec,T]);

idx=[i0,iUp,iPk,iRc];
end