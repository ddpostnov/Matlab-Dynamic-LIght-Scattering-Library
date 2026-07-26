%getMyographPropagation  Direction and speed of diameter changes travelling along a vessel
%
%   prop = getMyographPropagation(s,interval) estimates whether the
%   constriction/dilation of a myograph vessel travels ALONG the vessel
%   (the Y / row axis) within one interval, and if so how fast and in which
%   direction.  The vessel is straight and near-vertical, so the row index is a
%   proxy for distance along the vessel and a travelling wave shows up as a lag
%   of each location's diameter fluctuation relative to the others.
%
%   METHOD.  The reported speed/direction come from the LAG AT MAXIMUM
%   CROSS-CORRELATION of each location's (detrended, drift-removed) diameter
%   against a robust reference built from the well-correlated locations,
%   regressed on the row index (slope -> 1/speed, sign -> direction).  A wavelet-
%   free phase-delay estimate and an event-onset estimate are computed as
%   independent cross-checks.
%
%   The result is NOT thresholded: a speed is ALWAYS returned, accompanied by
%   EXPLAINABLE confidence metrics (median cross-correlation between locations,
%   the R^2 and p of the lag-vs-row fit, the fraction of locations used, the
%   agreement between methods, and whether the total lag is below the sampling
%   resolution) plus a one-line summary, so the user decides whether to trust it
%   ("a propagation speed of X was detected with confidence C").
%
%   INPUTS
%     s        parameter struct (defaults filled if missing/empty):
%                • vFR           vasomotion band [lo hi] Hz (dominant-freq search)
%                • pixelSize     µm per px; [] or 0 -> report in px (px/s)
%                • rowRange      [lo hi] rows that were measured (others excluded)
%                • propSignal    'diameter' (default; lag on the drift-removed
%                                diameter) | 'rData' (band-limited reconstruction)
%                • detrendSec    per-row high-pass window for 'diameter', s
%                • propMinMeasured min per-row measured fraction to keep a row (mask)
%                • propMinCoh    min |correlation| to accept a location into the fit
%                • propArtifactK MAD factor rejecting constant/erratic artifact rows
%                • propMaxLagFrac lag search limit, fraction of the dominant period
%                • propMinRows   minimum locations for an estimate
%                • propNShuffle  surrogate iterations for the p-value
%                • propResolutionSamples total lag below this many samples is flagged
%                                as below the sampling resolution (magnitude unreliable)
%     interval struct with fields time, diameter, mask (optional), vasomotion (optional).
%
%   OUTPUT
%     prop  struct with fields
%             • method           'max-correlation-lag'
%             • signalSource     'diameter' | 'rData'
%             • speed / speedUnit  propagation speed (px/s or µm/s) - ALWAYS finite
%             • direction        'upward' | 'downward'  (increasing row = downward)
%             • speedCI          95% CI of the speed (upper may be Inf near the floor)
%             • R2 / pValue      lag-vs-row fit quality and surrogate significance
%             • confidence       0..1 overall confidence (explained by .metrics/.confidenceText)
%             • confidenceLevel  'low' | 'medium' | 'high' (a label, NOT a gate)
%             • confidenceText   one-line human-readable justification
%             • metrics          struct: medianCorr, R2, pValue, rowFraction,
%                                totalLagSamples, belowResolution, methodAgreement,
%                                speedRatioLagPhase, vasomotionBandRatio, nRows
%             • belowResolution  true if the total lag is < propResolutionSamples
%             • phase / event    cross-check estimates (speed/direction)
%             • qualityFlags     cellstr diagnostics
%             • lagByRow         [row lag_s] used for the fit (diagnostics)
%             • diag             struct with map/fit data for Explore
%           An empty/too-weak interval returns a low-confidence result with the
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
% Last revision: 21-July-2026

function prop = getMyographPropagation(s,interval)

% ---- default parameters (filled if missing or empty) ----
def.vFR=[0.05 0.25];       def.pixelSize=[];         def.rowRange=[1 Inf];
def.propSignal='diameter'; def.detrendSec=30;        def.propMethods={'lag','event'};
def.propMinMeasured=0.5;   def.propMinCoh=0.3;        def.propArtifactK=4;
def.propMaxLagFrac=0.5;    def.propMinRows=20;        def.propNShuffle=200;
def.propResolutionSamples=1.0;
fnd=fieldnames(def);
for i=1:numel(fnd)
    if ~isfield(s,fnd{i}) || isempty(s.(fnd{i})), s.(fnd{i})=def.(fnd{i}); end
end
methods=s.propMethods; if ischar(methods), methods={methods}; end
wS=warning('off','stats:statrobustfit:IterationLimit'); wC=onCleanup(@()warning(wS)); %#ok<NASGU> (benign, on many permutations)

% ---- pick the input signal and its time base ----
[sig,tt,src]=selectSignal(s,interval);
prop=emptyProp(src);
if isempty(sig) || size(sig,2)<3
    prop.qualityFlags={'no-signal'}; prop.confidenceText='no usable diameter signal'; return;
end
tt=tt(:); nT=size(sig,1); nY=size(sig,2); fs=1/mean(diff(tt),'omitnan');

% ---- per-location fluctuation + location-quality mask ----
[Xf,keep,rowStd]=rowFluctuations(s,sig,src,fs,interval,nT,nY);
prop.diag=struct('Xf',Xf,'tt',tt,'y',(1:nY)','rowStd',rowStd(:),'keep0',keep(:));
dbg = isfield(s,'propDebug') && ~isempty(s.propDebug) && s.propDebug;
usable0=nnz(keep);
if dbg, fprintf('[prop] src=%s nT=%d nY=%d keep0=%d maskMed=%.2f\n',src,nT,nY,usable0,maskMed(interval,nY)); end
if usable0<s.propMinRows
    prop.nRows=usable0; prop.qualityFlags={'too-few-rows'};
    prop.confidenceText='too few measured locations'; return;
end

% ---- dominant vasomotion frequency + robust reference ----
domFreq=dominantFreq(s,interval,Xf,keep,fs);
f0=refineF0(mean(Xf(:,keep),2),fs,s.vFR,domFreq); prop.domFreq=f0; y=(1:nY)'; w=2*pi*f0;
[ref,keep]=robustReference(Xf,keep,rowStd);
if dbg, fprintf('[prop] after reference: keep=%d (f0=%.3f)\n',nnz(keep),f0); end
if nnz(keep)<s.propMinRows
    prop.nRows=nnz(keep); prop.qualityFlags={'too-few-coherent-rows'};
    prop.confidenceText='locations do not correlate (no coherent oscillation)'; return;
end
maxLag=max(2,round(s.propMaxLagFrac/max(f0,eps)*fs));

% ---- PRIMARY: lag at maximum cross-correlation of each location vs the reference ----
[lag,coh]=rowLags(Xf,ref,keep,maxLag,fs);
uidx=find(keep & isfinite(lag) & coh>s.propMinCoh & abs(lag)<(maxLag/fs)*0.98); uidx=uidx(:);
prop.nRows=numel(uidx);
if dbg, fprintf('[prop] final use=%d medianCorr=%.2f\n',numel(uidx),median(coh(uidx))); end
if prop.nRows<s.propMinRows
    prop.qualityFlags={'too-few-coherent-rows'};
    prop.confidenceText='too few well-correlated locations'; return;
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

% ---- CROSS-CHECK 1: phase-delay of the dominant oscillation ----
[phi,amp]=rowPhaseDelay(Xf,keep,tt,f0);
speedPh=NaN; dirPh='';
pu=intersect(uidx,find(isfinite(phi(:))));
if numel(pu)>=s.propMinRows
    [ys,o]=sort(y(pu)); pv=phi(pu); phu=unwrap(pv(o(:)));
    bP=robustfit(ys,phu); slopeP=bP(2)/w;
    speedPh=toSpeed(slopeP,s.pixelSize); dirPh=ternary(slopeP<0,'downward','upward');
end

% ---- explainable confidence metrics (NO thresholding of the speed) ----
medianCorr=median(coh(uidx));
rowFraction=prop.nRows/max(usable0,1);
belowResolution=totalLagSamples<s.propResolutionSamples;
[methodAgree,speedRatio]=agreeLagPhase(direction,dirPh,speed,speedPh);
bandR=vasoBandRatio(interval);
[conf,level,txt]=confidenceReport(medianCorr,R2,pVal,rowFraction,methodAgree,speedRatio,belowResolution,bandR);

prop.method='max-correlation-lag'; prop.signalSource=src;
prop.speed=speed; prop.speedUnit=unit; prop.direction=direction; prop.speedCI=speedCI;
prop.slope=slope; prop.R2=R2; prop.pValue=pVal; prop.lagByRow=[yv, lv];
prop.confidence=conf; prop.confidenceLevel=level; prop.confidenceText=txt;
prop.belowResolution=belowResolution;
prop.metrics=struct('medianCorr',medianCorr,'R2',R2,'pValue',pVal,'rowFraction',rowFraction, ...
    'totalLagSamples',totalLagSamples,'belowResolution',belowResolution,'methodAgreement',methodAgree, ...
    'speedRatioLagPhase',speedRatio,'vasomotionBandRatio',bandR,'nRows',prop.nRows);
prop.phase=struct('speed',speedPh,'direction',dirPh,'unit',unit);
prop.diag=struct('ref',ref,'tt',tt,'y',y,'uidx',uidx,'lag',lag(:),'coh',coh(:), ...
    'phi',phi(:),'amp',amp(:),'Xf',Xf,'f0',f0,'fitLag',[bL(1) bL(2)],'maxLag',maxLag);

% ---- CROSS-CHECK 2: event-onset regression ----
if ismember('event',methods)
    prop.event=eventEstimate(s,interval,Xf,ref,keep,fs,y);
end

% ---- informational quality flags (do not change the reported speed) ----
if belowResolution,           prop.qualityFlags=[prop.qualityFlags,{'lag-below-resolution'}]; end
if R2<0.15,                   prop.qualityFlags=[prop.qualityFlags,{'low-R2'}]; end
if medianCorr<0.4,            prop.qualityFlags=[prop.qualityFlags,{'weak-location-correlation'}]; end
if isfinite(bandR)&&bandR<3,  prop.qualityFlags=[prop.qualityFlags,{'weak-vasomotion'}]; end
end

% =====================================================================
function [sig,tt,src]=selectSignal(s,interval)
%SELECTSIGNAL  choose the diameter (default, for the lag method) or rData
%   rData is the band-limited reconstruction, stored at
%   vasomotion.timeVectors.VB.rData on the vasomotion.timeWT time base (present only
%   when 'reconstruction' was requested in s.segVsmReturn).
useR = strcmpi(s.propSignal,'rData') && isfield(interval,'vasomotion') && ~isempty(interval.vasomotion) ...
    && isstruct(interval.vasomotion) && isfield(interval.vasomotion,'timeVectors') ...
    && isfield(interval.vasomotion.timeVectors,'VB') && isfield(interval.vasomotion.timeVectors.VB,'rData') ...
    && ~isempty(interval.vasomotion.timeVectors.VB.rData);
if useR
    sig=double(interval.vasomotion.timeVectors.VB.rData); tt=double(interval.vasomotion.timeWT(:)); src='rData';
elseif isfield(interval,'diameter') && ~isempty(interval.diameter)
    sig=double(interval.diameter); tt=double(interval.time(:)); src='diameter';
    % off-FOV frames (a wall dilated out of view) carry a garbage edge-clamped
    % diameter - exclude them: NaN then interpolate in time so the series stays
    % uniform for the lag/cross-correlation (vasomotion is slow, so short gaps fill safely).
    if isfield(interval,'valid') && ~isempty(interval.valid)
        inv=~logical(interval.valid(:));
        if any(inv) && ~all(inv)
            sig(inv,:)=NaN;
            sig=fillmissing(sig,'linear',1,'EndValues','nearest');
        end
    end
else
    sig=[]; tt=[]; src='';
end
end

% =====================================================================
function [Xf,keep,rowStd]=rowFluctuations(s,sig,src,fs,interval,nT,nY)
%ROWFLUCTUATIONS  per-location drift-removed fluctuation + initial quality mask
if strcmp(src,'rData')
    Xf=sig-mean(sig,1,'omitnan');                 % rData is already band-limited
else
    win=max(3,round(s.detrendSec*fs));            % high-pass removes bulk drift, keeps vasomotion
    trend=movmean(sig,win,1,'omitnan','Endpoints','shrink');
    Xf=(sig-trend)./max(median(sig,1,'omitnan'),eps);
end
Xf(~isfinite(Xf))=0;
rowStd=std(Xf,0,1);
keep=true(1,nY);
r0=max(1,round(s.rowRange(1))); r1=min(nY,round(s.rowRange(2)));
keep([1:r0-1, r1+1:nY])=false;
if isfield(interval,'mask') && ~isempty(interval.mask) && size(interval.mask,2)==nY
    keep=keep & (mean(double(interval.mask),1)>=s.propMinMeasured);
end
med=median(rowStd(keep)); madv=mad(rowStd(keep),1)+eps;   % drop constant/erratic (artifact) rows
keep=keep & rowStd> max(0.15*med, med-s.propArtifactK*1.4826*madv) ...
          & rowStd< med+s.propArtifactK*1.4826*madv;
end

% =====================================================================
function [ref,keep]=robustReference(Xf,keep,rowStd) %#ok<INUSD>
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
function [phi,amp]=rowPhaseDelay(Xf,keep,tt,f0)
%ROWPHASEDELAY  phase & amplitude of the dominant oscillation per location (harmonic fit)
tt=tt(:); w=2*pi*f0; C=[cos(w*tt) sin(w*tt)]; M=(C'*C)\C';
nY=size(Xf,2); phi=nan(1,nY); amp=zeros(1,nY);
for c=find(keep)
    co=M*Xf(:,c); phi(c)=atan2(-co(1),co(2)); amp(c)=hypot(co(1),co(2));
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
function [agree,ratio]=agreeLagPhase(dirL,dirP,spL,spP)
if isempty(dirP)||~isfinite(spP), agree=NaN; ratio=NaN; return; end
agree=double(strcmp(dirL,dirP));
ratio=max(spL,spP)/max(min(spL,spP),eps);
end

% =====================================================================
function r=vasoBandRatio(interval)
%VASOBANDRATIO  VB-vs-CB activity ratio - the retired band-ratio marker rebuilt inline
%   (redesign decision 6): median band-mean amplitude in the vasomotion band divided by
%   the same in the control band.  A ratio < 3 flags weak vasomotion (see caller).
r=NaN;
if ~(isfield(interval,'vasomotion')&&isstruct(interval.vasomotion)), return; end
v=interval.vasomotion;
if ~(isfield(v,'scalars')&&isfield(v.scalars,'VB')&&isfield(v.scalars,'CB')), return; end
r = median(double(v.scalars.VB.ampMean),'omitnan') / ...
    median(double(v.scalars.CB.ampMean),'omitnan');
end

% =====================================================================
function [conf,level,txt]=confidenceReport(medianCorr,R2,pVal,rowFraction,methodAgree,speedRatio,belowRes,bandR)
%CONFIDENCEREPORT  transparent 0..1 confidence from its explainable components
%   Direction is usually robust; the speed MAGNITUDE is the uncertain part, so a
%   large lag-vs-phase speed ratio lowers the confidence and is called out.
corrScore = max(0,min(1,medianCorr));
fitScore  = max(0,min(1,R2/0.50));                         % R2 (lag ordering) is the propagation-specific signal
if ~isfinite(pVal), sigScore=0.5; elseif pVal<0.01, sigScore=1; elseif pVal<0.05, sigScore=0.6; else, sigScore=0.15; end
if isnan(methodAgree)                                      % no phase cross-check available
    agreeScore=0.7;
elseif methodAgree<=0                                      % directions disagree
    agreeScore=0.3;
elseif ~isfinite(speedRatio) || speedRatio<1.5            % agree, magnitudes close
    agreeScore=1;
elseif speedRatio<3
    agreeScore=0.75;                                       % agree on direction, magnitude loose
else
    agreeScore=0.5;                                        % agree on direction, magnitude uncertain
end
rowScore  = max(0.3,min(1,rowFraction));
conf = geomean([fitScore fitScore corrScore sigScore agreeScore rowScore]);  % lag-ordering (R2) weighted x2
if belowRes, conf=min(conf*0.5,0.30); end     % below the sampling resolution -> speed untrustworthy (low)
if conf>=0.6, level='high'; elseif conf>=0.35, level='medium'; else, level='low'; end
% ---- plain-language justification, naming the weak factors ----
bits={sprintf('corr=%.2f',medianCorr), sprintf('R2=%.2f',R2), sprintf('p=%.3f',pVal)};
if ~isnan(methodAgree), bits{end+1}=ternary(methodAgree>0,'direction agrees','direction DISAGREES'); end
notes={};
if medianCorr<0.4, notes{end+1}='locations weakly correlated'; end
if R2<0.15, notes{end+1}='lags poorly ordered along Y'; end
if isfinite(pVal)&&pVal>=0.05, notes{end+1}='not significant'; end
if isfinite(speedRatio)&&speedRatio>=2, notes{end+1}=sprintf('speed magnitude uncertain (lag/phase differ %.1fx)',speedRatio); end
if belowRes, notes{end+1}='total lag below sampling resolution (magnitude unreliable)'; end
if isfinite(bandR)&&bandR<3, notes{end+1}='weak vasomotion'; end
txt=sprintf('%s confidence (%s)',upper1(level),strjoin(bits,', '));
if ~isempty(notes), txt=[txt ' - ' strjoin(notes,'; ')]; end
end

% =====================================================================
function ev=eventEstimate(s,interval,Xf,ref,keep,fs,y) %#ok<INUSL>
%EVENTESTIMATE  independent speed from per-location constriction-onset times
ev=struct('speed',NaN,'direction','','nEvents',0,'speedUnit','px/s');
r=-(ref-mean(ref)); thr=std(r);
f0=dominantFreqFromRef(-r,fs);
minSep=round(0.5/max(f0,0.03)*fs); w=max(3,minSep);
[~,loc]=localMaxima(r,thr,minSep);
if numel(loc)<2, return; end
idx=find(keep); speeds=[]; dirs=[];
for e=1:numel(loc)
    lo=max(1,loc(e)-w); hi=min(size(Xf,1),loc(e)+w); tOn=nan(1,numel(idx));
    for k=1:numel(idx)
        seg=-(Xf(lo:hi,idx(k))-mean(Xf(lo:hi,idx(k)))); [~,ip]=max(seg); tOn(k)=(lo+ip-1)/fs;
    end
    good=isfinite(tOn);
    if nnz(good)>=s.propMinRows
        pp=polyfit(y(idx(good)),tOn(good)',1);
        if abs(pp(1))>eps, speeds(end+1)=1/abs(pp(1)); dirs(end+1)=sign(pp(1)); end %#ok<AGROW>
    end
end
ev.nEvents=numel(speeds);
if ev.nEvents>=1
    sp=median(speeds); ev.speed=sp; ev.direction=ternary(sign(median(dirs))<0,'downward','upward');
    if ~isempty(s.pixelSize)&&s.pixelSize>0, ev.speed=sp*s.pixelSize; ev.speedUnit='µm/s'; end
    ev.speedDist=speeds;
end
end

% =====================================================================
function pVal=permP(yv,lag,slope,nsh)
%PERMP  permutation null for the lag-vs-row slope (row-label shuffle)
if nsh<=0 || numel(yv)<10, pVal=NaN; return; end
yv=yv(:); lag=lag(:); n=numel(yv); sl=nan(nsh,1);
for is=1:nsh
    bb=robustfit(yv(randperm(n)),lag); sl(is)=bb(2);
end
pVal=(1+sum(abs(sl)>=abs(slope)))/(nsh+1);
end

% =====================================================================
function domFreq=dominantFreq(s,interval,Xf,keep,fs)
if isfield(interval,'vasomotion') && isstruct(interval.vasomotion) ...
        && isfield(interval.vasomotion,'fVectors') && isfield(interval.vasomotion.fVectors,'ampMean') ...
        && ~isempty(interval.vasomotion.fVectors.ampMean)
    domFreq=dominantFromSpct(interval.vasomotion,s); return;
end
ref=mean(Xf(:,keep),2); domFreq=dominantFreqFromRef(ref-mean(ref),fs,s.vFR);
end
function f0=dominantFromSpct(V,s)
sm=mean(double(V.fVectors.ampMean),1,'omitnan'); f=double(V.f(:));
in=f>=s.vFR(1)&f<=s.vFR(2); [~,im]=max(sm(in)); ff=f(in); f0=ff(im);
if isempty(f0)||~isfinite(f0), f0=mean(s.vFR); end
end
function f0=dominantFreqFromRef(x,fs,vFR)
if nargin<3, vFR=[0.03 0.3]; end
n=numel(x); P=abs(fft(x)).^2; fax=(0:n-1)'*(fs/n);
in=fax>=vFR(1)&fax<=vFR(2); P(~in)=0; [~,im]=max(P); f0=fax(im);
if isempty(f0)||f0<=0, f0=mean(vFR); end
end

% =====================================================================
function prop=emptyProp(src)
prop=struct('method','max-correlation-lag','signalSource',src,'domFreq',NaN, ...
    'speed',NaN,'speedUnit','px/s','direction','','speedCI',[NaN NaN],'slope',NaN, ...
    'R2',NaN,'pValue',NaN,'nRows',0,'confidence',0,'confidenceLevel','low', ...
    'confidenceText','','belowResolution',true,'metrics',struct(), ...
    'phase',struct('speed',NaN,'direction',''),'event',struct('speed',NaN,'nEvents',0), ...
    'qualityFlags',{{}},'lagByRow',[],'diag',struct());
end

% =====================================================================
function m=maskMed(interval,nY)
m=NaN;
if isfield(interval,'mask')&&~isempty(interval.mask)&&size(interval.mask,2)==nY
    m=median(mean(double(interval.mask),1));
end
end
function [mx,loc]=localMaxima(x,thr,minSep)
%LOCALMAXIMA  simple peak picker (avoids Signal Processing Toolbox findpeaks)
x=x(:); isPk=[false; x(2:end-1)>x(1:end-2) & x(2:end-1)>x(3:end); false] & x>thr;
loc=find(isPk); mx=x(loc); if isempty(loc), return; end
[~,ord]=sort(mx,'descend'); loc=loc(ord); keepp=true(size(loc));
for i=1:numel(loc)
    if ~keepp(i), continue; end
    keepp(abs(loc-loc(i))<minSep & (1:numel(loc))'>i)=false;
end
loc=sort(loc(keepp)); mx=x(loc);
end
function s=upper1(s), if ~isempty(s), s(1)=upper(s(1)); end, end
function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
