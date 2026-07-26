%registerRetinaLSCI - motion registration of raw retinal LSCI (.rls) recordings
%
%   [regTrustVector,s] = registerRetinaLSCI(rlsFileName, Name,Value,...)
%
% DESCRIPTION
%   Removes (mostly slow, translational) motion artefacts from a raw laser
%   speckle recording of the human retina / optic disc and writes a new,
%   registered .rls file (e.g. myfile.rls -> myfile_reg.rls).
%
%   In the raw frames the retinal features are too weak to register directly
%   and even a single spatial-contrast frame is very noisy.  Because the motion
%   is slow relative to the 100 fps frame rate, a short block of contrast frames
%   can be averaged without smearing the vasculature.  The function therefore:
%     1) streams the raw data in batches and computes spatial contrast (getK);
%     2) averages every 'avgFrames' contrast frames into a clean "keyframe" in
%        which the optic-disc vessels are visible;
%     3) enhances the keyframes (vessel emphasis) and registers them to a common
%        template built from the roughly-aligned keyframes (phase correlation,
%        translation or rigid) - a couple of refinement passes;
%     4) fills failed keyframes from their neighbours, smooths and interpolates
%        the per-keyframe transform to every frame;
%     5) streams the raw data again, applies the per-frame transform and writes
%        the registered .rls (the header is copied verbatim, so frame count,
%        dtype, timestamps and sampling are preserved).
%
%   The raw data is only RELOCATED, never smoothed: translation uses integer-
%   pixel shifts and rigid uses nearest-neighbour resampling, so the speckle
%   values - and therefore the spatial contrast - are preserved everywhere
%   except the vacated edge.  (Sub-pixel motion is thus rounded to whole pixels;
%   this is the correct trade-off for LSCI, where interpolation would corrupt the
%   speckle statistics.)
%
% INPUT
%   rlsFileName   path to the raw .rls file.
%
% OPTIONAL Name-Value
%   'procType'      'gpu' (default) | 'cpu'  - contrast computation.
%   'kernelSize'    spatial-contrast kernel (default 5).
%   'avgFrames'     contrast frames averaged per keyframe (default 25 -> 0.25 s
%                   at 100 fps).  Larger = cleaner keyframes but more motion
%                   blur within a block.
%   'transformType' 'translation' (default) | 'rigid' (adds rotation).
%   'refPasses'     template-refinement passes (default 2).
%   'rotationLimit' deg; a rigid keyframe rotating more than this is treated as a
%                   failed registration (default 10).
%   'maxShift'      px; keyframe shift larger than this is treated as failed
%                   ([] -> auto = 0.40*min(sizeY,sizeX)).
%   'smoothSpan'    keyframes; moving-median smoothing of the shift series before
%                   interpolation (odd, default 3; 1 = off).
%   'trustCorrThr'  min correlation of a warped keyframe to the template for the
%                   registration to be trusted (default 0.35).
%   'intraMotionThr' px; within-block motion (shift between the two halves of a
%                   keyframe) above this flags the block as rapid motion
%                   (default 4).  Catches fast motion that blurs the keyframe
%                   average but averages out between block centres.
%   'stepMotionThr' px; a jump larger than this between consecutive keyframes
%                   flags rapid motion / a saccade (default 8).
%   'trustMargin'   keyframes flagged around every bad keyframe (default 1).
%   'regMaxDim'     keyframes are down-sampled so their larger side is at most
%                   this many px before registration (default 512) - faster and
%                   less sensitive to residual speckle.
%   'outputSuffix'  appended to the file name (default '_reg').
%   'writeFile'     write the registered .rls (default true).  false = estimate
%                   only (returns transforms + trust, writes nothing).
%   'preview'       save a *_reg_preview.png diagnostic (default true).
%   'memoryCoef'    fraction of free RAM used for batching (default 0.7).
%
% OUTPUT
%   regTrustVector  sizeT x 1 logical.  1 = low motion and a confident
%                   registration; 0 = strong motion / failed registration and
%                   the frames around it.  Use it to weight or discard frames
%                   in the subsequent contrast analysis.
%   s               struct with the estimated transforms and diagnostics:
%                     .outputFileName, .transformType, .avgFrames, .size [Y X T]
%                     .kfCenters   keyframe centre frame indices
%                     .kfShift     [nKf x 2] estimated [x y] shift per keyframe
%                     .kfTheta     [nKf x 1] rotation (deg, 0 if translation)
%                     .regCorr     [nKf x 1] template correlation per keyframe
%                     .stepMotion  [nKf x 1] inter-keyframe drift (px)
%                     .intraMotion [nKf x 1] within-block motion (px)
%                     .trustKf     [nKf x 1] per-keyframe trust
%                     .template    enhanced registration template
%                     .params      the resolved parameter set
%
% DEPENDS ON
%   readRLS, getK, enhanceForRegistration, and MATLAB's Image Processing Toolbox
%   (imregtform/imregconfig, imwarp, imtranslate).
%
% See also: getContrastFromRLS, getK, readRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 23-July-2026

function [regTrustVector,s] = registerRetinaLSCI(rlsFileName, varargin)

p = inputParser;
addRequired(p,'rlsFileName',@(x) ischar(x)||isstring(x));
addParameter(p,'procType','gpu',@(x) any(validatestring(x,{'cpu','gpu'})));
addParameter(p,'kernelSize',5,@(x) isnumeric(x)&&isscalar(x)&&x>=3);
addParameter(p,'avgFrames',25,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'transformType','translation',@(x) any(validatestring(x,{'translation','rigid'})));
addParameter(p,'refPasses',2,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'rotationLimit',10,@(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'maxShift',[],@(x) isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'smoothSpan',3,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'trustCorrThr',0.35,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'intraMotionThr',4,@(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'stepMotionThr',8,@(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'regMaxDim',512,@(x) isnumeric(x)&&isscalar(x)&&x>=64);
addParameter(p,'trustMargin',1,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'outputSuffix','_reg',@(x) ischar(x)||isstring(x));
addParameter(p,'writeFile',true,@(x) islogical(x)||ismember(x,[0 1]));
addParameter(p,'preview',true,@(x) islogical(x)||ismember(x,[0 1]));
addParameter(p,'memoryCoef',0.7,@(x) isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
parse(p,rlsFileName,varargin{:});
r = p.Results;
rlsFileName = char(r.rlsFileName);
avg = round(r.avgFrames);

procType = r.procType;
if strcmpi(procType,'gpu') && gpuDeviceCount==0
    warning('registerRetinaLSCI:noGPU','No GPU found - falling back to CPU.');
    procType = 'cpu';
end
[regOpt,regMet] = imregconfig('monomodal'); regOpt.MaximumIterations = 300;

% ---- metadata --------------------------------------------------------------
[~,~,meta] = readRLS(rlsFileName,'KeepOpen',true,'FramesToRead',1);
fclose(meta.fId);
sizeY = meta.sizeY; sizeX = meta.sizeX; sizeT = meta.sizeT;
dataType = meta.dataType; dataSize = meta.dataSize;
nKf = max(1,floor(sizeT/avg));
maxShift = r.maxShift; if isempty(maxShift), maxShift = 0.40*min(sizeY,sizeX); end

% ---- batch size (whole number of keyframes) --------------------------------
[~,mem] = memory; memAvail = mem.PhysicalMemory.Available;
batchSize = max(avg,floor(memAvail*r.memoryCoef/(sizeX*sizeY*(dataSize+16))));
batchSize = min(max(avg,floor(batchSize/avg)*avg),floor(sizeT/avg)*avg);
if batchSize<avg, batchSize = avg; end

% ===========================================================================
% PASS A - spatial contrast -> averaged keyframes
% ===========================================================================
% Each keyframe is accumulated as two halves so the within-block motion (fast
% motion that blurs the average but averages out between block centres) can be
% measured as the shift between the first and second half.
kfSumA = zeros(sizeY,sizeX,nKf,'single'); kfCntA = zeros(nKf,1);
kfSumB = zeros(sizeY,sizeX,nKf,'single'); kfCntB = zeros(nKf,1);
processed = 0; stream = []; tA = tic;
fprintf('registerRetinaLSCI: %d frames, %d keyframes\n',sizeT,nKf);
while processed<sizeT
    nb = min(batchSize,sizeT-processed);
    if isempty(stream)
        [raw,~,stream] = readRLS(rlsFileName,'KeepOpen',true,'FramesToRead',nb);
    else
        [raw,~,stream] = readRLS(stream,'FramesToRead',nb);
    end
    K = getK(raw,'spatial','kernelSize',r.kernelSize,'procType',procType);
    firstKf = floor(processed/avg)+1;
    lastKf  = min(ceil((processed+nb)/avg),nKf);
    for k = firstKf:lastKf
        glo = (k-1)*avg+1;
        ghi = k*avg; if k==nKf, ghi = sizeT; end        % last keyframe absorbs the remainder
        mid = glo+ceil((ghi-glo+1)/2)-1;                % split the block into two halves
        flo = max(glo-processed,1);     fhi = min(mid-processed,nb);   % first half
        if fhi>=flo, kfSumA(:,:,k)=kfSumA(:,:,k)+sum(K(:,:,flo:fhi),3); kfCntA(k)=kfCntA(k)+(fhi-flo+1); end
        flo = max((mid+1)-processed,1); fhi = min(ghi-processed,nb);   % second half
        if fhi>=flo, kfSumB(:,:,k)=kfSumB(:,:,k)+sum(K(:,:,flo:fhi),3); kfCntB(k)=kfCntB(k)+(fhi-flo+1); end
    end
    processed = processed+nb;
    fprintf('\r  contrast + keyframes:  %5.1f%%  [%.0fs]',100*processed/sizeT,toc(tA));
end
fprintf('\n');
if isfield(stream,'fId') && ~isempty(fopen(stream.fId)), fclose(stream.fId); end
clear raw K

kfCnt = kfCntA+kfCntB;
kf  = (kfSumA+kfSumB) ./ reshape(max(kfCnt,1),1,1,nKf);  % averaged keyframes
kfA = kfSumA ./ reshape(max(kfCntA,1),1,1,nKf);
kfB = kfSumB ./ reshape(max(kfCntB,1),1,1,nKf);
clear kfSumA kfSumB
intraMotion = zeros(nKf,1);                              % within-block motion (px)
for k = 1:nKf
    if kfCntA(k)>=2 && kfCntB(k)>=2
        [tx,ty] = regOneRaw(enhanceKeyframe(kfB(:,:,k)),enhanceKeyframe(kfA(:,:,k)),r.regMaxDim,regOpt,regMet);
        intraMotion(k) = hypot(tx,ty);
    end
    fprintf('\r  within-block motion:   %5.1f%%',100*k/nKf);
end
fprintf('\n');
clear kfA kfB
kfCenters = ((1:nKf)'-1)*avg + (avg+1)/2;
kfCenters(nKf) = ((nKf-1)*avg+1 + sizeT)/2;             % true centre of the last (fuller) block

% ===========================================================================
% Keyframe enhancement + registration to a common template
% ===========================================================================
kfE = zeros(sizeY,sizeX,nKf,'single');
for k = 1:nKf, kfE(:,:,k) = enhanceKeyframe(kf(:,:,k)); end

[kfShift,kfTheta,regCorr,template] = estimateKeyframeTransforms( ...
    kfE,r.transformType,round(r.refPasses),r.rotationLimit,maxShift,r.regMaxDim,regOpt,regMet);

% ---- trust + robust shift series -------------------------------------------
stepMotion = zeros(nKf,1);
if nKf>1
    stepMotion(2:end) = sqrt(sum(diff(kfShift,1,1).^2,2));
    stepMotion(1) = stepMotion(2);
end
badReg = regCorr(:) < r.trustCorrThr | isnan(kfShift(:,1));
% failed registration: drop the estimate and interpolate it from good neighbours
kfShiftFix = kfShift; kfThetaFix = kfTheta;
kfShiftFix(badReg,:) = NaN; kfThetaFix(badReg) = NaN;
kfShiftFix = fillmissing(kfShiftFix,'linear',1,'EndValues','nearest');
kfThetaFix = fillmissing(kfThetaFix,'linear','EndValues','nearest');
if all(isnan(kfShiftFix(:))), kfShiftFix(:) = 0; end
if all(isnan(kfThetaFix)),    kfThetaFix(:) = 0; end
if r.smoothSpan>1 && nKf>2
    kfShiftFix = smoothdata(kfShiftFix,1,'movmedian',round(r.smoothSpan));
    kfThetaFix = smoothdata(kfThetaFix,1,'movmedian',round(r.smoothSpan));
end

% Only rapid, LARGE motion is distrusted: a motion-blurred block (large shift
% between its two halves) or a large jump between consecutive blocks.  Slow
% drift and cardiac pulsation are normal and stay trusted.
motionBad = intraMotion > r.intraMotionThr | stepMotion > r.stepMotionThr;
trustKf = ~(badReg | motionBad);
if r.trustMargin>0                                      % also distrust the neighbours
    trustKf = ~logical(movmax(double(~trustKf),2*round(r.trustMargin)+1));
end

% ---- interpolate the transform to every frame ------------------------------
gAll = (1:sizeT)';
perX = interpClamp(kfCenters,kfShiftFix(:,1),gAll);
perY = interpClamp(kfCenters,kfShiftFix(:,2),gAll);
perTheta = interpClamp(kfCenters,kfThetaFix,gAll);
regTrustVector = false(sizeT,1);
for k = 1:nKf
    glo = (k-1)*avg+1;
    ghi = k*avg; if k==nKf, ghi = sizeT; end
    regTrustVector(glo:ghi) = trustKf(k);
end

% ===========================================================================
% PASS B - apply the transform to the raw data and write the registered .rls
% ===========================================================================
[pp,nn,ee] = fileparts(rlsFileName);
outName = fullfile(pp,[nn,char(r.outputSuffix),ee]);
isRigid = strcmpi(r.transformType,'rigid');
Rout = imref2d([sizeY sizeX]);
if r.writeFile
    fi = fopen(rlsFileName,'r'); hdr = fread(fi,30*1024,'*uint8'); fclose(fi);
    fo = fopen(outName,'w');
    if fo==-1, error('registerRetinaLSCI:write','Cannot open %s for writing.',outName); end
    cleanup = onCleanup(@() fclose(fo));
    fwrite(fo,hdr,'uint8');
    processed = 0; stream = []; tB = tic;
    fprintf('registerRetinaLSCI: writing %s\n',outName);
    while processed<sizeT
        nb = min(batchSize,sizeT-processed);
        if isempty(stream)
            [raw,ts,stream] = readRLS(rlsFileName,'KeepOpen',true,'FramesToRead',nb);
        else
            [raw,ts,stream] = readRLS(stream,'FramesToRead',nb);
        end
        for f = 1:nb
            g = processed+f;
            if isRigid          % nearest-neighbour: relocates speckle, no smoothing of raw data
                w = imwarp(raw(:,:,f),rigidtform2d(perTheta(g),[perX(g) perY(g)]), ...
                    'nearest','OutputView',Rout,'FillValues',0);
            else                % pure integer-pixel shift: raw speckle is only relocated, never smoothed
                w = shiftIntFrame(raw(:,:,f),perX(g),perY(g));
            end
            fwrite(fo,ts(f),'uint64');
            fwrite(fo,cast(w,dataType),dataType);
        end
        processed = processed+nb;
        fprintf('\r  applying + writing:    %5.1f%%  [%.0fs]',100*processed/sizeT,toc(tB));
    end
    fprintf('\n');
    if isfield(stream,'fId') && ~isempty(fopen(stream.fId)), fclose(stream.fId); end
    clear cleanup
else
    outName = '';
end

% ---- pack output + preview -------------------------------------------------
s = struct('outputFileName',outName,'transformType',r.transformType, ...
    'avgFrames',avg,'size',[sizeY sizeX sizeT],'kfCenters',kfCenters, ...
    'kfShift',kfShift,'kfTheta',kfTheta,'regCorr',regCorr,'stepMotion',stepMotion, ...
    'intraMotion',intraMotion,'trustKf',trustKf,'template',template,'params',r);

if r.preview
    try savePreview(fullfile(pp,[nn,char(r.outputSuffix),'_preview.png']), ...
            template,kfCenters,kfShiftFix,regCorr,trustKf,sizeT,avg,nKf); catch, end %#ok<CTCH>
end
fprintf('registerRetinaLSCI: done (%.1f%% of frames trusted).\n',100*mean(regTrustVector));
end

% ===========================================================================
function e = enhanceKeyframe(img)
%enhanceKeyframe  Emphasise the (dark, low-contrast) vessels for registration.
%   Thin wrapper over the shared primitive enhanceForRegistration: no mask
%   (whole-frame percentile clip) plus the Gaussian speckle-suppression tail
%   (sigma = 1).
e = enhanceForRegistration(img,'smoothSigma',1);
end

% ===========================================================================
function [kfShift,kfTheta,regCorr,template] = estimateKeyframeTransforms(kfE,type,nPass,rotLim,maxShift,regMaxDim,opt,met)
%estimateKeyframeTransforms  Register enhanced keyframes to a common template.
%   Intensity-based imregtform on a down-sampled copy - phase correlation
%   (imregcorr) is unreliable here because it locks onto the fixed camera
%   vignette / residual speckle rather than the moving vasculature.
[sy,sx,nKf] = size(kfE);
kfShift = zeros(nKf,2); kfTheta = zeros(nKf,1); regCorr = zeros(nKf,1);
if nKf==1, template = kfE; regCorr = 1; return; end

fx = max(1,max(sy,sx)/regMaxDim);                       % down-sample factor
if fx>1
    kfD = zeros(round(sy/fx),round(sx/fx),nKf,'single');
    for k=1:nKf, kfD(:,:,k) = imresize(kfE(:,:,k),1/fx); end
else
    kfD = kfE;
end
[dy,dx,~] = size(kfD); Rout = imref2d([dy dx]);
isRigid = strcmpi(type,'rigid');

template = kfD(:,:,round(nKf/2));                        % start from the middle keyframe
aligned = zeros(dy,dx,nKf,'single');
for pass = 1:nPass
    for k = 1:nKf
        [tx,ty,th,ok] = regOne(kfD(:,:,k),template,type,rotLim,maxShift/fx,opt,met);
        kfShift(k,:) = [tx ty]*fx; kfTheta(k) = th;     % translation back to full resolution
        if ok
            if isRigid
                aligned(:,:,k) = imwarp(kfD(:,:,k),rigidtform2d(th,[tx ty]),'OutputView',Rout,'FillValues',0);
            else
                aligned(:,:,k) = imtranslate(kfD(:,:,k),[tx ty],'OutputView','same','FillValues',0);
            end
        else
            aligned(:,:,k) = kfD(:,:,k);
            kfShift(k,:) = NaN;
        end
        fprintf('\r  registering keyframes: pass %d/%d  %5.1f%%',pass,nPass,100*k/nKf);
    end
    template = median(aligned,3);                        % robust template for the next pass
end
fprintf('\n');
m = template>0;                                          % correlation over the valid region
for k = 1:nKf
    a = aligned(:,:,k);
    if nnz(m)>10, regCorr(k) = corrSafe(a(m),template(m)); else, regCorr(k) = 0; end
end
template = imresize(template,[sy sx]);                   % full-resolution template for the preview
end

% ===========================================================================
function [tx,ty,th,ok] = regOne(moving,fixed,type,rotLim,maxShift,opt,met)
%regOne  Single intensity-based registration with sanity checks.  A diverged or
%   out-of-range result is treated as a failed keyframe (filled from neighbours).
tx = 0; ty = 0; th = 0; ok = true;
ws = warning('off','all');
try
    tform = imregtform(moving,fixed,type,opt,met);
    [tx,ty] = transformPointsForward(tform,0,0);
    if isprop(tform,'RotationAngle'), th = tform.RotationAngle; end
catch
    ok = false;
end
warning(ws);
if ~ok || ~isfinite(tx) || ~isfinite(ty) || hypot(tx,ty)>maxShift || abs(th)>rotLim
    tx = 0; ty = 0; th = 0; ok = false;
end
end

% ===========================================================================
function [tx,ty] = regOneRaw(moving,fixed,regMaxDim,opt,met)
%regOneRaw  Rough translational shift between two images (down-sampled).
f = max(1,max(size(moving))/regMaxDim);
if f>1, moving = imresize(moving,1/f); fixed = imresize(fixed,1/f); end
ws = warning('off','all');
try
    [tx,ty] = transformPointsForward(imregtform(moving,fixed,'translation',opt,met),0,0);
    tx = tx*f; ty = ty*f;
    if ~isfinite(tx) || ~isfinite(ty), tx = 0; ty = 0; end
catch
    tx = 0; ty = 0;
end
warning(ws);
end

% ===========================================================================
function c = corrSafe(a,b)
a = double(a(:)); b = double(b(:));
a = a-mean(a); b = b-mean(b);
d = sqrt(sum(a.^2)*sum(b.^2));
if d>0, c = sum(a.*b)/d; else, c = 0; end
end

% ===========================================================================
function out = shiftIntFrame(img,dx,dy)
%shiftIntFrame  Rigid integer-pixel translation of a raw frame (NO interpolation).
%   Content moves by round([dx dy]) in [x=col, y=row]; the raw speckle values are
%   only relocated (never averaged/smoothed) so the spatial contrast is preserved
%   everywhere except the vacated edge, which is filled with 0.
dx = round(dx); dy = round(dy);
[H,W] = size(img);
out = zeros(H,W,'like',img);
r = max(1,1-dy):min(H,H-dy);                            % source rows kept in frame
c = max(1,1-dx):min(W,W-dx);                            % source cols kept in frame
if ~isempty(r) && ~isempty(c)
    out(r+dy,c+dx) = img(r,c);
end
end

% ===========================================================================
function v = interpClamp(x,y,xi)
%interpClamp  Linear interpolation that holds the end values (no extrapolation).
if numel(x)<2, v = repmat(y(1),size(xi)); return; end
v = interp1(x,y,xi,'linear');
v(xi<x(1)) = y(1); v(xi>x(end)) = y(end);
end

% ===========================================================================
function savePreview(fname,template,kfC,kfShift,regCorr,trustKf,sizeT,avg,nKf)
f = figure('Visible','off','Color','w','Position',[100 100 1100 450]);
subplot(1,2,1); imagesc(template); axis image off; colormap(gca,'gray');
title('Registration template');
subplot(1,2,2); hold on;
plot(kfC,kfShift(:,1),'-o','MarkerSize',3,'DisplayName','shift x');
plot(kfC,kfShift(:,2),'-o','MarkerSize',3,'DisplayName','shift y');
yl = ylim;
for k = 1:nKf                                            % shade distrusted blocks
    if ~trustKf(k)
        glo = (k-1)*avg+0.5; ghi = k*avg+0.5; if k==nKf, ghi = sizeT+0.5; end
        patch([glo ghi ghi glo],[yl(1) yl(1) yl(2) yl(2)],[1 0.4 0.4], ...
            'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
    end
end
plot(kfC,regCorr*diff(yl)+yl(1),':k','DisplayName','template corr (scaled)');
xlabel('frame'); ylabel('shift (px)'); legend('Location','best'); grid on;
title(sprintf('Estimated motion (%.0f%% trusted)',100*mean(repelem(trustKf,avg,1))));
exportgraphics(f,fname); close(f);
end
