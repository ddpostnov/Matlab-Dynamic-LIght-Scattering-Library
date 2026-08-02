%getMyographPropagation  Direction and speed of diameter changes travelling along a vessel
%
%   prop = getMyographPropagation(s,sig,time,mask,valid) estimates whether the
%   constriction/dilation of a myograph vessel travels ALONG the vessel (the Y /
%   row axis), and if so how fast and in which direction.  The vessel is straight
%   and near-vertical, so the row index is a proxy for distance along the vessel
%   and a travelling wave shows up as a lag of each location's diameter
%   fluctuation relative to the others.
%
%   METHOD - one estimator, the LAG AT MAXIMUM CROSS-CORRELATION.  Each location's
%   drift-removed diameter is cross-correlated against a robust reference built
%   from the well-correlated locations; the lag at the correlation maximum is
%   regressed on the row index (slope -> 1/speed, sign -> direction).  There is no
%   phase-delay estimate, no event-onset estimate and no vasomotion input: the
%   speed comes from the diameter and nothing else.
%
%   The result is NOT thresholded: a speed is ALWAYS returned, accompanied by
%   EXPLAINABLE confidence metrics (median cross-correlation between locations,
%   the R^2 and p of the lag-vs-row fit, the fraction of locations used, and
%   whether the total lag is below the sampling resolution) plus a plain-language
%   summary, so the user decides whether to trust it.  The guard against
%   over-claiming is the surrogate: shuffling the row labels destroys the lag
%   ordering, R^2 collapses and the confidence goes to zero.
%
%   INPUTS
%     s        parameter struct (defaults filled if missing/empty):
%                • vFR           vasomotion band [lo hi] Hz (dominant-freq search)
%                • pixelSize     µm per px; [] or 0 -> report in px (px/s)
%                • rowRange      [lo hi] rows that were measured (others excluded)
%                • detrendSec    per-row high-pass window removing bulk drift, s
%                • propMinMeasured min per-row measured fraction to keep a row (mask)
%                • propMinCoh    min |correlation| to accept a location into the fit
%                • propArtifactK MAD factor rejecting constant/erratic artifact rows
%                • propMaxLagFrac lag search limit, fraction of the dominant period
%                • propMinRows   minimum locations for an estimate
%                • propNShuffle  surrogate iterations for the p-value
%                • propResolutionSamples total lag below this many samples is flagged
%                                as below the sampling resolution (magnitude unreliable)
%     sig      [nT x nY] diameter of ONE measure, one column per location along Y.
%     time     [nT x 1] frame times, s.
%     mask     [nT x nY] 1 = measured, 0 = interpolated, or [] (all measured).
%     valid    [nT x 1] false = a wall was off-FOV that frame, or [] (all valid).
%              Off-FOV frames carry an edge-clamped diameter, so they are removed
%              and interpolated over before anything else happens.
%
%   OUTPUT
%     prop  struct with fields
%             • method           'max-correlation-lag'
%             • speed / speedUnit  propagation speed (px/s or µm/s) - ALWAYS finite
%             • direction        'upward' | 'downward'  (increasing row = downward)
%             • speedCI          95% CI of the speed (upper may be Inf near the floor)
%             • R2 / pValue      lag-vs-row fit quality and surrogate significance
%             • confidence       0..1 overall confidence (explained by .metrics/.confidenceText)
%             • confidenceLevel  'low' | 'medium' | 'high' (a label, NOT a gate)
%             • confidenceText   one plain sentence the user can act on
%             • metrics          struct: medianCorr, R2, pValue, rowFraction,
%                                totalLagSamples, belowResolution, nRows
%             • belowResolution  true if the total lag is < propResolutionSamples
%             • domFreq          dominant frequency used to size the lag search, Hz
%             • qualityFlags     cellstr diagnostics
%             • lagByRow         [row lag_s] used for the fit (diagnostics)
%             • diag             struct with map/fit data for Explore
%           An empty/too-weak signal returns a low-confidence result with the
%           reason in confidenceText, never an error.
%
%   Direction convention: increasing row index = downward in the image.
%
%   DEPENDS ON  base MATLAB + Statistics Toolbox (robustfit/prctile/mad).
%   Self-contained; no Signal Processing Toolbox (a local normXcorr is used).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

function prop = getMyographPropagation(s,sig,time,mask,valid)

if nargin<4, mask=[]; end
if nargin<5, valid=[]; end

% ---- default parameters (filled if missing or empty) ----
def.vFR=[0.05 0.25];       def.pixelSize=[];          def.rowRange=[1 Inf];
def.detrendSec=30;         def.propMinMeasured=0.5;   def.propMinCoh=0.3;
def.propArtifactK=4;       def.propMaxLagFrac=0.5;    def.propMinRows=20;
def.propNShuffle=200;      def.propResolutionSamples=1.0;
fnd=fieldnames(def);
for i=1:numel(fnd)
    if ~isfield(s,fnd{i}) || isempty(s.(fnd{i})), s.(fnd{i})=def.(fnd{i}); end
end
wS=warning('off','stats:statrobustfit:IterationLimit'); wC=onCleanup(@()warning(wS)); %#ok<NASGU> (benign, on many permutations)

% ---- the signal: the drift-removed diameter, off-FOV frames interpolated over ----
% This used to take the whole interval struct.  Fail loudly rather than reporting
% "no usable signal" if an un-migrated caller still hands one over.
if ~isempty(sig) && ~(isnumeric(sig)||islogical(sig))
    error('getMyographPropagation:signalNotArray', ...
        ['getMyographPropagation now takes plain arrays: ' ...
         'prop = getMyographPropagation(s,sig,time,mask,valid), with sig [nT x nY] ' ...
         'for ONE diameter measure.']);
end
prop=emptyProp();
if isempty(sig) || size(sig,2)<3
    prop.qualityFlags={'no-signal'}; prop.confidenceText='No usable diameter signal.'; return;
end
sig=double(sig); tt=double(time(:));
% off-FOV frames (a wall dilated out of view) carry a garbage edge-clamped diameter -
% exclude them: NaN then interpolate in time so the series stays uniform for the
% cross-correlation (vasomotion is slow, so short gaps fill safely).
if ~isempty(valid)
    inv=~logical(valid(:));
    if any(inv) && ~all(inv)
        sig(inv,:)=NaN;
        sig=fillmissing(sig,'linear',1,'EndValues','nearest');
    end
end
nT=size(sig,1); nY=size(sig,2); fs=1/mean(diff(tt),'omitnan');

% ---- per-location fluctuation + location-quality mask ----
[Xf,keep,rowStd]=rowFluctuations(s,sig,fs,mask,nY);
prop.diag=struct('Xf',Xf,'tt',tt,'y',(1:nY)','rowStd',rowStd(:),'keep0',keep(:));
dbg = isfield(s,'propDebug') && ~isempty(s.propDebug) && s.propDebug;
usable0=nnz(keep);
if dbg, fprintf('[prop] nT=%d nY=%d keep0=%d maskMed=%.2f\n',nT,nY,usable0,maskMed(mask,nY)); end
if usable0<s.propMinRows
    prop.nRows=usable0; prop.qualityFlags={'too-few-rows'};
    prop.confidenceText='Too few measured locations along the vessel to look for propagation.'; return;
end

% ---- dominant oscillation frequency (sizes the lag search) + robust reference ----
ref0=mean(Xf(:,keep),2);
domFreq=dominantFreqFromRef(ref0-mean(ref0),fs,s.vFR);
f0=refineF0(ref0,fs,s.vFR,domFreq); prop.domFreq=f0; y=(1:nY)';
[ref,keep]=robustReference(Xf,keep);
if dbg, fprintf('[prop] after reference: keep=%d (f0=%.3f)\n',nnz(keep),f0); end
if nnz(keep)<s.propMinRows
    prop.nRows=nnz(keep); prop.qualityFlags={'too-few-coherent-rows'};
    prop.confidenceText='The locations do not correlate with each other, so there is no coherent oscillation to track.'; return;
end
maxLag=max(2,round(s.propMaxLagFrac/max(f0,eps)*fs));

% ---- the estimator: lag at maximum cross-correlation of each location vs the reference ----
[lag,coh]=rowLags(Xf,ref,keep,maxLag,fs);
uidx=find(keep & isfinite(lag) & coh>s.propMinCoh & abs(lag)<(maxLag/fs)*0.98); uidx=uidx(:);
prop.nRows=numel(uidx);
if dbg, fprintf('[prop] final use=%d medianCorr=%.2f\n',numel(uidx),median(coh(uidx))); end
if prop.nRows<s.propMinRows
    prop.qualityFlags={'too-few-coherent-rows'};
    prop.confidenceText='Too few well-correlated locations to fit a propagation speed.'; return;
end
yv=y(uidx); lv=lag(uidx); lv=lv(:);
[bL,stL]=robustfit(yv,lv); slope=bL(2); se=stL.se(2);
yhat=bL(1)+bL(2)*yv;
R2=1-sum((lv-yhat).^2)/max(sum((lv-mean(lv)).^2),eps);
rangeY=max(yv)-min(yv); totalLagSamples=abs(slope)*rangeY*fs;
pVal=permP(yv,lv,slope,s.propNShuffle);
[speed,unit]=toSpeed(slope,s.pixelSize);
speedCI=speedCIfromSlope(slope,se,s.pixelSize);
direction=ternary(slope<0,'downward','upward');

% ---- explainable confidence metrics (NO thresholding of the speed) ----
medianCorr=median(coh(uidx));
rowFraction=prop.nRows/max(usable0,1);
belowResolution=totalLagSamples<s.propResolutionSamples;
[conf,level,txt]=confidenceReport(medianCorr,R2,pVal,rowFraction,belowResolution);

prop.method='max-correlation-lag';
prop.speed=speed; prop.speedUnit=unit; prop.direction=direction; prop.speedCI=speedCI;
prop.slope=slope; prop.R2=R2; prop.pValue=pVal; prop.lagByRow=[yv, lv];
prop.confidence=conf; prop.confidenceLevel=level; prop.confidenceText=txt;
prop.belowResolution=belowResolution;
prop.metrics=struct('medianCorr',medianCorr,'R2',R2,'pValue',pVal,'rowFraction',rowFraction, ...
    'totalLagSamples',totalLagSamples,'belowResolution',belowResolution,'nRows',prop.nRows);
prop.diag=struct('ref',ref,'tt',tt,'y',y,'uidx',uidx,'lag',lag(:),'coh',coh(:), ...
    'Xf',Xf,'f0',f0,'fitLag',[bL(1) bL(2)],'maxLag',maxLag);

% ---- informational quality flags (do not change the reported speed) ----
if belowResolution,           prop.qualityFlags=[prop.qualityFlags,{'lag-below-resolution'}]; end
if R2<0.15,                   prop.qualityFlags=[prop.qualityFlags,{'low-R2'}]; end
if medianCorr<0.4,            prop.qualityFlags=[prop.qualityFlags,{'weak-location-correlation'}]; end
end

% =====================================================================
function [Xf,keep,rowStd]=rowFluctuations(s,sig,fs,mask,nY)
%ROWFLUCTUATIONS  per-location drift-removed fluctuation + initial quality mask
win=max(3,round(s.detrendSec*fs));                % high-pass removes bulk drift, keeps vasomotion
trend=movmean(sig,win,1,'omitnan','Endpoints','shrink');
Xf=(sig-trend)./max(median(sig,1,'omitnan'),eps);
Xf(~isfinite(Xf))=0;
rowStd=std(Xf,0,1);
keep=true(1,nY);
r0=max(1,round(s.rowRange(1))); r1=min(nY,round(s.rowRange(2)));
keep([1:r0-1, r1+1:nY])=false;
if ~isempty(mask) && size(mask,2)==nY
    keep=keep & (mean(double(mask),1)>=s.propMinMeasured);
end
med=median(rowStd(keep)); madv=mad(rowStd(keep),1)+eps;   % drop constant/erratic (artifact) rows
keep=keep & rowStd> max(0.15*med, med-s.propArtifactK*1.4826*madv) ...
          & rowStd< med+s.propArtifactK*1.4826*madv;
end

% =====================================================================
function [ref,keep]=robustReference(Xf,keep)
%ROBUSTREFERENCE  reference waveform from the coherent, clean locations (two-pass)
ref0=mean(Xf(:,keep),2);
c=corrTo(Xf,ref0);
core=keep & c>max(0.5,prctile(c(keep),50));
if nnz(core)<10, core=keep; end
ref=mean(Xf(:,core),2);
c2=corrTo(Xf,ref);
keep=keep & c2>0.3;
end

% =====================================================================
function c=corrTo(Xf,ref)
r=ref-mean(ref); rn=r/max(norm(r),eps);
X=Xf-mean(Xf,1); Xn=X./max(sqrt(sum(X.^2,1)),eps);
c=(rn'*Xn); c=c(:)';
end

% =====================================================================
function [lagS,coh]=rowLags(Xf,ref,keep,maxLag,fs)
%ROWLAGS  lag (s) & value of the maximum cross-correlation of each location vs the reference
nY=size(Xf,2); lagS=nan(1,nY); coh=zeros(1,nY);
r=ref-mean(ref);
for c=find(keep)
    x=Xf(:,c)-mean(Xf(:,c));
    [cc,ll]=normXcorr(x,r,maxLag);
    [mx,ip]=max(cc); coh(c)=mx;
    if ip>1 && ip<numel(cc)
        d=(cc(ip-1)-cc(ip+1))/(2*(cc(ip-1)-2*cc(ip)+cc(ip+1))+eps);
        lagS(c)=(ll(ip)+max(min(d,0.5),-0.5))/fs;
    else
        lagS(c)=ll(ip)/fs;
    end
end
end

% =====================================================================
function [c,lags]=normXcorr(x,y,maxLag)
%NORMXCORR  normalised cross-correlation for lags -maxLag..maxLag (no Signal Toolbox)
x=x(:); y=y(:); n=numel(x); lags=(-maxLag:maxLag)';
den=sqrt(sum(x.^2)*sum(y.^2))+eps; c=zeros(numel(lags),1);
for k=1:numel(lags)
    L=lags(k);
    if L>=0, a=x(1+L:n); b=y(1:n-L); else, a=x(1:n+L); b=y(1-L:n); end
    c(k)=sum(a.*b)/den;
end
end

% =====================================================================
function f0=refineF0(x,fs,vFR,f0guess)
%REFINEF0  dominant frequency in the band, parabolic-interpolated to sub-bin
x=x-mean(x); n=numel(x); P=abs(fft(x)).^2; fax=(0:n-1)'*(fs/n);
idx=find(fax>=vFR(1)&fax<=vFR(2));
if isempty(idx), f0=f0guess; return; end
[~,k]=max(P(idx)); kk=idx(k);
if kk>1 && kk<numel(fax)
    a=log(P(kk-1)+eps); b=log(P(kk)+eps); c=log(P(kk+1)+eps);
    d=0.5*(a-c)/(a-2*b+c+eps); d=max(min(d,0.5),-0.5);
    f0=(kk-1+d)*(fs/n);
else
    f0=fax(kk);
end
if ~isfinite(f0) || f0<vFR(1) || f0>vFR(2), f0=f0guess; end
end

% =====================================================================
function f0=dominantFreqFromRef(x,fs,vFR)
%DOMINANTFREQFROMREF  strongest in-band FFT bin of the mean fluctuation
if nargin<3, vFR=[0.03 0.3]; end
n=numel(x); P=abs(fft(x)).^2; fax=(0:n-1)'*(fs/n);
in=fax>=vFR(1)&fax<=vFR(2); P(~in)=0; [~,im]=max(P); f0=fax(im);
if isempty(f0)||f0<=0, f0=mean(vFR); end
end

% =====================================================================
function [speed,unit]=toSpeed(slope,pixelSize)
%TOSPEED  slope (s/px) -> speed (px/s or µm/s); always finite
speed=1/max(abs(slope),eps);
if ~isempty(pixelSize) && pixelSize>0, speed=speed*pixelSize; unit='µm/s'; else, unit='px/s'; end
end

% =====================================================================
function ci=speedCIfromSlope(slope,se,pixelSize)
%SPEEDCIFROMSLOPE  95% CI on the speed from the slope standard error (upper may be Inf)
lo=abs(slope)-1.96*se; hi=abs(slope)+1.96*se;
sLo=1/max(hi,eps); sHi=1/max(lo,eps); if lo<=0, sHi=Inf; end
if ~isempty(pixelSize) && pixelSize>0, sLo=sLo*pixelSize; if isfinite(sHi), sHi=sHi*pixelSize; end, end
ci=[sLo sHi];
end

% =====================================================================
function pVal=permP(yv,lag,slope,nsh)
%PERMP  permutation null for the lag-vs-row slope (row-label shuffle)
%   THE over-claiming guard: shuffling the row labels keeps every waveform and every
%   lag but destroys their ORDER along the vessel, which is the only thing a
%   travelling wave adds.  A real wave survives it; noise does not.
if nsh<=0 || numel(yv)<10, pVal=NaN; return; end
yv=yv(:); lag=lag(:); n=numel(yv); sl=nan(nsh,1);
for is=1:nsh
    bb=robustfit(yv(randperm(n)),lag); sl(is)=bb(2);
end
pVal=(1+sum(abs(sl)>=abs(slope)))/(nsh+1);
end

% =====================================================================
function [conf,level,txt]=confidenceReport(medianCorr,R2,pVal,rowFraction,belowRes)
%CONFIDENCEREPORT  transparent 0..1 confidence from its explainable components
%   Every component is something the user can see in the diagnostics: how well the
%   locations oscillate together, how cleanly their lags line up along the vessel,
%   how unlikely that ordering is by chance, and how much of the vessel contributed.
%   The lag ordering (R2) is the propagation-specific signal, so it counts twice.
corrScore = max(0,min(1,medianCorr));
fitScore  = max(0,min(1,R2/0.50));
if ~isfinite(pVal), sigScore=0.5; elseif pVal<0.01, sigScore=1; elseif pVal<0.05, sigScore=0.6; else, sigScore=0.15; end
rowScore  = max(0.3,min(1,rowFraction));
conf = geomean([fitScore fitScore corrScore sigScore rowScore]);
if belowRes, conf=min(conf*0.5,0.30); end     % below the sampling resolution -> speed untrustworthy (low)
if conf>=0.6, level='high'; elseif conf>=0.35, level='medium'; else, level='low'; end

% ---- one plain sentence, naming what carried or spoiled the estimate ----
if R2>=0.50
    fitTxt=sprintf('the lags line up cleanly along the vessel (R2=%.2f)',R2);
elseif R2>=0.15
    fitTxt=sprintf('the lags increase along the vessel with some scatter (R2=%.2f)',R2);
else
    fitTxt=sprintf('the lags are poorly ordered along the vessel (R2=%.2f)',R2);
end
if medianCorr>=0.4
    corrTxt=sprintf('the locations oscillate together (correlation %.2f)',medianCorr);
else
    corrTxt=sprintf('the locations correlate only weakly (correlation %.2f)',medianCorr);
end
if ~isfinite(pVal)
    sigTxt='and no significance test was run';
elseif pVal<0.05
    sigTxt=sprintf('and that ordering is unlikely by chance (p=%.3f)',pVal);
else
    sigTxt=sprintf('and that ordering is no better than chance (p=%.3f)',pVal);
end
txt=sprintf('%s confidence: %s, %s %s',upper1(level),fitTxt,corrTxt,sigTxt);
if belowRes
    txt=[txt '; the whole vessel lags by less than one frame, so the direction may hold but the speed does not'];
end
if rowFraction<0.5
    txt=[txt sprintf('; only %.0f%% of the measured locations could be used',100*rowFraction)];
end
txt=[txt '.'];
end

% =====================================================================
function prop=emptyProp()
prop=struct('method','max-correlation-lag','domFreq',NaN, ...
    'speed',NaN,'speedUnit','px/s','direction','','speedCI',[NaN NaN],'slope',NaN, ...
    'R2',NaN,'pValue',NaN,'nRows',0,'confidence',0,'confidenceLevel','low', ...
    'confidenceText','','belowResolution',true,'metrics',struct(), ...
    'qualityFlags',{{}},'lagByRow',[],'diag',struct());
end

% =====================================================================
function m=maskMed(mask,nY)
m=NaN;
if ~isempty(mask) && size(mask,2)==nY
    m=median(mean(double(mask),1));
end
end
function s=upper1(s), if ~isempty(s), s(1)=upper(s(1)); end, end
function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
