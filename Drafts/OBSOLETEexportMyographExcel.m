%exportMyographExcel  Export myograph results to a multi-sheet Excel workbook
%
%   outPath = exportMyographExcel(sources,xlsxPath) writes one .xlsx summarising
%   one or more *_myograph.mat result files (SOURCES: a path, a cellstr of paths,
%   or a loaded struct with fields intervals/s/meta).  Sheets:
%     settings    - s flattened + file, frame rate, pixel size, crop, timestamp, code version
%     intervals   - name, start, end, duration, nFrames, nY, measured fraction
%     vasomotion markers - long: file, interval, Y, the band-suffixed vasomotion scalars
%                   (ampMean/ampStd/ampSkew per band, the VB frequency-shape set
%                   fCent*/fSprd*/shapePeak/nPeak*, and the flare/silence clustering
%                   dur*/amp* per band; scalars.VB/CB.*)
%     ampPctScalar - long: file, interval, Y, band, percentile, amp - band-envelope
%                   amplitude percentiles scalars.VB/CB.ampPct at the s.pcts LEVELS
%     propagation - file, interval, method, speed, unit, direction, CI, p, R2, confidence, flags
%     spectra     - long: file, interval, Y, f, ampMean, ampStd (decimated in f/Y)
%     ampPctSpectra - file, interval, f, one column per pctCenters bin - the VB
%                   percentile-resolved spectra fVectors.VB.ampMeanPct averaged over Y (decimated in f)
%     diameter    - median-over-Y (and optionally per-Y) diameter time course (decimated in t),
%                   for the default analysed measure (wall centre) of the three measured
%
%   The 1,048,576 x 16,384 sheet limits are respected by decimating the spectra
%   (in f) and the diameter (in t), and the decimation factor is written into the
%   sheet.  Large arrays (spectrum.amp/.phase, timeVectors, rData, idxL/idxR/mask) are never
%   exported in full.
%
%   opts (name/value or struct): fDecim (auto), tDecim (auto), perY (false),
%   maxRows (1e6).
%
%   DEPENDS ON  loadMyographResults, myographMeasureIndex; base MATLAB
%   writetable/writecell.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

function outPath = exportMyographExcel(sources,xlsxPath,opts)

if nargin<3||isempty(opts), opts=struct(); end
def.fDecim=[]; def.tDecim=[]; def.perY=false; def.maxRows=1e6;
def.exportFullSpectra=false; def.fullSpectraFile=[];   % full per-Y x f goes to a SEPARATE file
fn=fieldnames(def); for i=1:numel(fn), if ~isfield(opts,fn{i}), opts.(fn{i})=def.(fn{i}); end, end
recs=loadSources(sources);
if isempty(recs), error('exportMyographExcel:noData','No loadable results.'); end
if exist(xlsxPath,'file'), delete(xlsxPath); end     % start clean (avoid stale sheets)

writeSettings(recs,xlsxPath);
writeIntervals(recs,xlsxPath);
writeVsmMarkers(recs,xlsxPath,opts);
writeAmpPctScalar(recs,xlsxPath,opts);
writePropagation(recs,xlsxPath);
writeSpectra(recs,xlsxPath,opts);
writeAmpPctSpectra(recs,xlsxPath,opts);
writeDiameter(recs,xlsxPath,opts);
outPath=xlsxPath;
end

% =====================================================================
function recs=loadSources(sources)
recs=struct('label',{},'intervals',{},'s',{},'meta',{});
if ischar(sources)||isstring(sources), sources={char(sources)}; end
if isstruct(sources)
    for i=1:numel(sources)
        recs(end+1)=oneRec(labelOf(sources(i)),sources(i).intervals,getf(sources(i),'s'),getf(sources(i),'meta')); %#ok<AGROW>
    end
    return;
end
for i=1:numel(sources)
    [iv,s,meta]=loadMyographResults(sources{i});
    [~,nn]=fileparts(sources{i});
    recs(end+1)=oneRec(nn,iv,s,meta); %#ok<AGROW>
end
end
function r=oneRec(lab,iv,s,meta), r=struct('label',lab,'intervals',iv,'s',s,'meta',meta); end
function v=getf(S,f), if isfield(S,f), v=S.(f); else, v=struct(); end, end
function l=labelOf(S), if isfield(S,'meta')&&isfield(S.meta,'fName'), [~,l]=fileparts(S.meta.fName); else, l='data'; end, end

% =====================================================================
function writeSettings(recs,xlsx)
rows={};
rows(end+1,:)={'field','value'}; %#ok<*AGROW>
for r=1:numel(recs)
    m=recs(r).meta; s=recs(r).s;
    rows(end+1,:)={sprintf('=== file %d ===',r),recs(r).label};
    rows(end+1,:)={'fName',getStr(m,'fName')};
    rows(end+1,:)={'frameRate',num2str(getNum(m,'frameRate'))};
    rows(end+1,:)={'pixelSize_um_per_px',pxStr(m)};
    rows(end+1,:)={'timeCrop',mat2strSafe(getf2(m,'timeCrop'))};
    rows(end+1,:)={'rowRange',mat2strSafe(getf2(m,'rowRange'))};
    rows(end+1,:)={'codeVersion',getStr(m,'codeVersion')};
    rows(end+1,:)={'created',getStr(m,'createdTimestamp')};
    rows(end+1,:)={'formatVersion',num2str(getNum(m,'formatVersion'))};
    if isstruct(s)
        fn=fieldnames(s);
        for i=1:numel(fn)
            v=s.(fn{i});
            if isnumeric(v)||islogical(v), rows(end+1,:)={fn{i},mat2strSafe(v)};
            elseif ischar(v), rows(end+1,:)={fn{i},v}; end
        end
    end
end
writecell(rows,xlsx,'Sheet','settings');
end

% =====================================================================
function writeIntervals(recs,xlsx)
T=table();
for r=1:numel(recs)
    iv=recs(r).intervals;
    for k=1:numel(iv)
        t=iv(k).time; D=iv(k).diameter;
        mf=NaN; if ~isempty(iv(k).mask), mf=mean(iv(k).mask(:)); end
        T=[T; cell2table({recs(r).label,iv(k).name,minv(t),maxv(t),spanv(t), ...
            numel(t),cols(D),mf}, 'VariableNames', ...
            {'file','interval','start_s','end_s','duration_s','nFrames','nY','measuredFraction'})];
    end
end
writetable(T,xlsx,'Sheet','intervals');
end

% =====================================================================
function writeVsmMarkers(recs,xlsx,opts)
%WRITEVSMMARKERS  long-format per-Y vasomotion scalars -> sheet 'vasomotion markers' (band-branched tree, section 2).
%   Columns are the per-Y scalars.VB/CB.* markers, VB/CB-suffixed for the flat table:
%   envelope moments ampMean/ampStd/ampSkew (both bands), the VB-only frequency-shape /
%   multiplicity set fCent*/fSprd*/shapePeak/nPeak*, and the flare/silence clustering
%   dur*/amp* (both bands; CB from its own independent Otsu mask).  The low-value markers
%   retired by the redesign are dropped (see DATA-MODEL section 9); the per-Y percentile
%   VECTOR scalars.*.ampPct is not a single-value marker, so it is exported separately to
%   the 'ampPctScalar' sheet rather than here.
smk={ 'ampMeanVB','VB','ampMean';   'ampStdVB','VB','ampStd';   'ampSkewVB','VB','ampSkew'; ...
      'ampMeanCB','CB','ampMean';   'ampStdCB','CB','ampStd';   'ampSkewCB','CB','ampSkew'; ...
      'fCentMeanVB','VB','fCentMean';   'fCentStdVB','VB','fCentStd'; ...
      'fSprdMeanVB','VB','fSprdMean';   'fSprdStdVB','VB','fSprdStd'; ...
      'shapePeakVB','VB','shapePeak';   'nPeakMeanVB','VB','nPeakMean';   'nPeakStdVB','VB','nPeakStd'; ...
      'durFlareMeanVB','VB','durFlareMean';       'durFlareStdVB','VB','durFlareStd'; ...
      'durSilenceMeanVB','VB','durSilenceMean';   'durSilenceStdVB','VB','durSilenceStd'; ...
      'ampFlareMeanVB','VB','ampFlareMean';       'ampFlareStdVB','VB','ampFlareStd'; ...
      'ampSilenceMeanVB','VB','ampSilenceMean';   'ampSilenceStdVB','VB','ampSilenceStd'; ...
      'durFlareMeanCB','CB','durFlareMean';       'durFlareStdCB','CB','durFlareStd'; ...
      'durSilenceMeanCB','CB','durSilenceMean';   'durSilenceStdCB','CB','durSilenceStd'; ...
      'ampFlareMeanCB','CB','ampFlareMean';       'ampFlareStdCB','CB','ampFlareStd'; ...
      'ampSilenceMeanCB','CB','ampSilenceMean';   'ampSilenceStdCB','CB','ampSilenceStd' };
C={}; hdr=[{'file','interval','Y'},smk(:,1)']; C(1,:)=hdr;
for r=1:numel(recs)
    iv=recs(r).intervals;
    for k=1:numel(iv)
        v=iv(k).vasomotion; if isempty(v)||~isstruct(v)||~isfield(v,'scalars'), continue; end
        nY=numel(v.scalars.VB.ampMean);
        for yy=1:nY
            row={recs(r).label,iv(k).name,yy};
            for b=1:size(smk,1), row{end+1}=scalarMk(v,smk{b,2},smk{b,3},yy); end
            C(end+1,:)=row;
            if size(C,1)>opts.maxRows, break; end
        end
    end
end
writecell(C,xlsx,'Sheet','vasomotion markers');
end

function y=scalarMk(v,band,name,yy)
%SCALARMK  one per-Y band scalar v.scalars.<band>.<name> (NaN if absent or too short).
y=NaN;
if isfield(v,'scalars')&&isfield(v.scalars,band)&&isfield(v.scalars.(band),name)
    d=double(v.scalars.(band).(name)); if numel(d)>=yy, y=d(yy); end
end
end

% =====================================================================
function writeAmpPctScalar(recs,xlsx,opts)
%WRITEAMPPCTSCALAR  long-format per-Y band-envelope amplitude percentiles -> 'ampPctScalar'.
%   scalars.VB/CB.ampPct is [nY x numel(s.pcts)] = prctile(a(t),s.pcts): the band-mean
%   envelope amplitude a(t) evaluated at each percentile LEVEL in s.pcts (0..100 %).  One
%   long row per (file,interval,Y,band,percentile).  NOTE these percentile LEVELS index
%   ampPct - NOT the pctCenters bin CENTRES that key the fVector percentile spectra
%   (sheet 'ampPctSpectra').  s.pcts read from rec.s (fallback 0:10:100) like vFRof.
C={}; C(1,:)={'file','interval','Y','band','percentile','amp'};
for r=1:numel(recs)
    iv=recs(r).intervals; pcts=pctsOf(recs(r));
    for k=1:numel(iv)
        v=iv(k).vasomotion; if isempty(v)||~isstruct(v)||~isfield(v,'scalars'), continue; end
        for bc={'VB','CB'}
            band=bc{1};
            if ~isfield(v.scalars,band)||~isfield(v.scalars.(band),'ampPct'), continue; end
            P=double(v.scalars.(band).ampPct);              % [nY x numel(s.pcts)]
            np=min(numel(pcts),size(P,2));
            for yy=1:size(P,1)
                for j=1:np
                    C(end+1,:)={recs(r).label,iv(k).name,yy,band,pcts(j),P(yy,j)}; %#ok<AGROW>
                    if size(C,1)>opts.maxRows, writecell(C,xlsx,'Sheet','ampPctScalar'); return; end
                end
            end
        end
    end
end
writecell(C,xlsx,'Sheet','ampPctScalar');
end

% =====================================================================
function writeAmpPctSpectra(recs,xlsx,opts)
%WRITEAMPPCTSPECTRA  VB percentile-resolved amplitude spectra, averaged over Y -> 'ampPctSpectra'.
%   fVectors.VB.ampMeanPct is [nY x nF x nBin] (mean |W| per VB-envelope amplitude bin);
%   averaged over Y -> [nF x nBin], written as f_Hz + one column per bin, each column
%   labelled by its pctCenters bin CENTRE (v.pctCenters).  Decimated in f with the same
%   fDecim logic as writeSpectra.  NOTE the bin CENTRES here differ from the percentile
%   LEVELS s.pcts that index the scalar ampPct (sheet 'ampPctScalar').
C={}; haveHdr=false; nBinH=0;
for r=1:numel(recs)
    iv=recs(r).intervals;
    for k=1:numel(iv)
        v=iv(k).vasomotion;
        if isempty(v)||~isstruct(v)||~isfield(v,'fVectors')||~isfield(v.fVectors,'VB') ...
                ||~isfield(v.fVectors.VB,'ampMeanPct')||isempty(v.fVectors.VB.ampMeanPct), continue; end
        f=double(v.f(:)); nf=numel(f);
        M=squeeze(mean(double(v.fVectors.VB.ampMeanPct),1,'omitnan'));   % [nF x nBin] mean over Y
        if isvector(M), M=M(:); end
        if ~haveHdr
            nBinH=size(M,2); pc=double(v.pctCenters(:)); labs=cell(1,nBinH);
            for b=1:nBinH, if b<=numel(pc), labs{b}=sprintf('pct_%g',pc(b)); else, labs{b}=sprintf('bin%d',b); end, end
            C(1,:)=[{'file','interval','f_Hz'},labs,{'fDecim'}]; haveHdr=true;
        end
        fdec=opts.fDecim; if isempty(fdec), fdec=max(1,ceil(nf/60)); end
        for j=1:fdec:nf
            C(end+1,:)=[{recs(r).label,iv(k).name,f(j)},num2cell(fitRow(M(j,:),nBinH)),{fdec}]; %#ok<AGROW>
            if size(C,1)>opts.maxRows, writecell(C,xlsx,'Sheet','ampPctSpectra'); return; end
        end
    end
end
if ~haveHdr, C={'file','interval','f_Hz'}; end     % header-only sheet when no percentile spectra present
writecell(C,xlsx,'Sheet','ampPctSpectra');
end
function pcts=pctsOf(rec)
%PCTSOF  percentile LEVELS s.pcts for this record (fallback 0:10:100), like vFRof for vFR.
pcts=0:10:100;
if isfield(rec,'s')&&isstruct(rec.s)&&isfield(rec.s,'pcts')&&~isempty(rec.s.pcts), pcts=double(rec.s.pcts); end
end
function y=fitRow(x,n)
%FITROW  1xN row: the first n elements of x, NaN-padded if x is shorter (fixed sheet width).
x=double(x(:)'); y=nan(1,n); m=min(n,numel(x)); y(1:m)=x(1:m);
end

% =====================================================================
function D=analysedDiameter(ivk)
%ANALYSEDDIAMETER  the wall-centre diameter of one interval as [frames x nY]
%   getMyographDiameter measures outer / wall centre / lumen and stacks them in the
%   3rd dimension; the exported trace is the default analysed one.  A result written
%   before the three-measure change is already [frames x nY] and comes back unchanged.
if isfield(ivk,'measures'), meas=ivk.measures; else, meas={}; end
k=min(myographMeasureIndex('',meas),size(ivk.diameter,3));
D=double(ivk.diameter(:,:,k));
end

% =====================================================================
function writePropagation(recs,xlsx)
C={}; C(1,:)={'file','interval','method','speed','unit','direction','speedCI_lo','speedCI_hi', ...
    'confidence','confidenceLevel','medianCorr','R2','pValue','rowFraction','totalLagSamples', ...
    'belowResolution','domFreq','nRows','flags','confidenceText'};
for r=1:numel(recs)
    iv=recs(r).intervals;
    for k=1:numel(iv)
        p=iv(k).prop; if isempty(p)||~isstruct(p), continue; end
        m=getf(p,'metrics');
        C(end+1,:)={recs(r).label,iv(k).name,gs(p,'method'),gn(p,'speed'),gs(p,'speedUnit'),gs(p,'direction'), ...
            sci(p,1),sci(p,2),gn(p,'confidence'),gs(p,'confidenceLevel'),mn(m,'medianCorr'),gn(p,'R2'),gn(p,'pValue'), ...
            mn(m,'rowFraction'),mn(m,'totalLagSamples'),gn(p,'belowResolution'),gn(p,'domFreq'),gn(p,'nRows'), ...
            strjoin(gc(p,'qualityFlags'),';'),gs(p,'confidenceText')}; %#ok<AGROW>
    end
end
writecell(C,xlsx,'Sheet','propagation');
end

% =====================================================================
function writeSpectra(recs,xlsx,opts)
%WRITESPECTRA  REFINED spectra: per-interval spectrum averaged over Y (+ scalar band
%   metrics), decimated in f.  The full per-Y x f spectra are NOT put here - they go
%   to a separate file only when opts.exportFullSpectra / opts.fullSpectraFile is set.
C={}; C(1,:)={'file','interval','f_Hz','ampMean_avgY','ampStd_avgY','peakF_vFR','vFR_power','fDecim'};
for r=1:numel(recs)
    iv=recs(r).intervals; vFR=vFRof(recs(r));
    for k=1:numel(iv)
        v=iv(k).vasomotion; if isempty(v)||~isstruct(v)||~isfield(v,'fVectors')||~isfield(v.fVectors,'ampMean'), continue; end
        f=double(v.f(:)); nf=numel(f);
        sm=mean(double(v.fVectors.ampMean),1,'omitnan'); sm=sm(:);
        ss=mean(double(v.fVectors.ampStd),1,'omitnan');  ss=ss(:);
        inV=f>=vFR(1)&f<=vFR(2); ff=f(inV);
        if any(inV), [~,im]=max(sm(inV)); pkF=ff(im); vpow=trapz(f(inV),sm(inV)); else, pkF=NaN; vpow=NaN; end
        fdec=opts.fDecim; if isempty(fdec), fdec=max(1,ceil(nf/60)); end
        for j=1:fdec:nf
            C(end+1,:)={recs(r).label,iv(k).name,f(j),sm(j),ss(j),pkF,vpow,fdec}; %#ok<AGROW>
        end
    end
end
writecell(C,xlsx,'Sheet','spectra');
if opts.exportFullSpectra || ~isempty(opts.fullSpectraFile)
    ff=opts.fullSpectraFile;
    if isempty(ff), [pp,nn]=fileparts(xlsx); ff=fullfile(pp,[nn '_fullSpectra.xlsx']); end
    writeFullSpectra(recs,ff,opts);
end
end

% =====================================================================
function writeFullSpectra(recs,ff,opts)
%WRITEFULLSPECTRA  full per-Y x f ampMean/ampStd to a SEPARATE workbook (decimated, bounded)
if exist(ff,'file'), delete(ff); end
C={}; C(1,:)={'file','interval','Y','f_Hz','ampMean','ampStd','fDecim','yDecim'};
for r=1:numel(recs)
    iv=recs(r).intervals;
    for k=1:numel(iv)
        v=iv(k).vasomotion; if isempty(v)||~isstruct(v)||~isfield(v,'fVectors')||~isfield(v.fVectors,'ampMean'), continue; end
        f=double(v.f(:)); nf=numel(f); nY=size(v.fVectors.ampMean,1);
        fdec=opts.fDecim; if isempty(fdec), fdec=max(1,ceil(nf/60)); end
        ydec=max(1,ceil(nY/200));
        for yy=1:ydec:nY
            for j=1:fdec:nf
                C(end+1,:)={recs(r).label,iv(k).name,yy,f(j),dbl2(v.fVectors.ampMean,yy,j),dbl2(v.fVectors.ampStd,yy,j),fdec,ydec}; %#ok<AGROW>
                if size(C,1)>opts.maxRows, writecell(C,ff,'Sheet','fullSpectra'); return; end
            end
        end
    end
end
writecell(C,ff,'Sheet','fullSpectra');
end
function vFR=vFRof(rec)
vFR=[0.05 0.25];
if isfield(rec,'s')&&isstruct(rec.s)&&isfield(rec.s,'vFR')&&numel(rec.s.vFR)==2, vFR=rec.s.vFR; end
end

% =====================================================================
function writeDiameter(recs,xlsx,opts)
C={}; C(1,:)={'file','interval','t_s','medianDiameter_px','tDecim'};
for r=1:numel(recs)
    iv=recs(r).intervals;
    for k=1:numel(iv)
        t=double(iv(k).time(:)); D=analysedDiameter(iv(k));
        med=median(D,2,'omitnan'); nT=numel(t);
        tdec=opts.tDecim; if isempty(tdec), tdec=max(1,ceil(nT/5000)); end
        ti=1:tdec:nT;
        for j=ti
            C(end+1,:)={recs(r).label,iv(k).name,t(j),med(j),tdec};
            if size(C,1)>opts.maxRows, writecell(C,xlsx,'Sheet','diameter'); return; end
        end
    end
end
writecell(C,xlsx,'Sheet','diameter');
end

% =====================================================================
function y=dbl2(M,i,j), M=double(M); if size(M,1)>=i && size(M,2)>=j, y=M(i,j); else, y=NaN; end, end
function y=gn(p,f), if isfield(p,f)&&~isempty(p.(f)), y=double(p.(f)); y=y(1); else, y=NaN; end, end
function y=gs(p,f), if isfield(p,f)&&(ischar(p.(f))||isstring(p.(f))), y=char(p.(f)); else, y=''; end, end
function y=gc(p,f), if isfield(p,f)&&iscell(p.(f)), y=p.(f); else, y={}; end, end
function y=sci(p,i), y=NaN; if isfield(p,'speedCI')&&numel(p.speedCI)>=i, y=p.speedCI(i); end, end
function y=mn(m,f), if isstruct(m)&&isfield(m,f)&&~isempty(m.(f)), y=double(m.(f)); y=y(1); else, y=NaN; end, end
function s=getStr(m,f), if isfield(m,f)&&(ischar(m.(f))||isstring(m.(f))), s=char(m.(f)); else, s=''; end, end
function n=getNum(m,f), if isfield(m,f)&&~isempty(m.(f))&&isnumeric(m.(f)), n=m.(f); n=n(1); else, n=NaN; end, end
function v=getf2(m,f), if isfield(m,f), v=m.(f); else, v=[]; end, end
function s=pxStr(m), if isfield(m,'pixelSize')&&~isempty(m.pixelSize)&&m.pixelSize>0, s=num2str(m.pixelSize); else, s='(uncalibrated - px)'; end, end
function s=mat2strSafe(v), if isempty(v), s=''; elseif isnumeric(v)||islogical(v), s=mat2str(double(v)); else, s=char(string(v)); end, end
function y=minv(t), t=double(t); if isempty(t), y=NaN; else, y=min(t); end, end
function y=maxv(t), t=double(t); if isempty(t), y=NaN; else, y=max(t); end, end
function y=spanv(t), t=double(t); if isempty(t), y=NaN; else, y=max(t)-min(t); end, end
function y=cols(D), y=size(D,2); end
