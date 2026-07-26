%getVasomotionMetrics  Shared wavelet vasomotion core (basis setup + per-series analysis)
%
%   getVasomotionMetrics centralises the continuous-wavelet vasomotion "science"
%   shared by runVasomotion (LSCI per-segment and per-pixel paths) and
%   getMyographVasomotion.  It replaces the previously-duplicated local helpers
%   vsmSetup and vsmSpectrumMetrics (which were byte-identical in both files).
%   Two call modes are dispatched by argument count:
%
%   SETUP  (nargin 2)
%       layout = getVasomotionMetrics(time, s)
%     Builds the analytic-Morlet cwtfilterbank for the sampling implied by `time`
%     and runs one throwaway wt on the time vector to recover the wavelet
%     frequency vector fwt and cone of influence coi.  Returns a small `layout`
%     struct built ONCE per time base and reused for every series on that base:
%       .fb .fwt .coi          wavelet bank, frequencies (Hz), cone of influence
%       .f                     fwt augmented with the s.vFR/s.cFR band edges (desc)
%       .idxsVFR .idxsCFR      logical masks selecting those bands within .f
%       .fs                    sampling frequency, Hz
%       .timeVSM               cone-of-influence-trimmed analysed time window
%       .want                  per-level compute flags from s.segVsmReturn (see below)
%       .nrmStat .nPeakProm    per-series invariants hoisted here so every series/pixel
%                              reuses them: the (x-c)/c centring-statistic handle and
%                              the findpeaks prominence fraction for the VB peak count
%
%   ANALYSIS  (nargin 3)
%       m = getVasomotionMetrics(series, layout, s)    % series = [nT x 1] trace
%     The science for ONE vector series, MODULAR via s.segVsmReturn: only the
%     requested analysis levels are computed and returned.  The series is normalised
%     internally into the relative fluctuation (x-c)/c with centre c per
%     s.normalisation/normsize (per-column and therefore identical to the callers'
%     old whole-matrix normalisation), transformed with the prebuilt bank, optionally
%     reconstructed (icwt over s.vFR) and rescaled with the series' own centre, then
%     reduced to the per-frequency |CWT| spectrum and passed through the metric core.
%
%   ANALYSIS levels/tokens (s.segVsmReturn; a cell subset of the names below).  Default
%   (absent/empty) is the COMPLETE set of COARSE levels {'bands','moments','series',
%   'clustering','reconstruction','spectrum'}; callers narrow it as needed.  The bag `m`
%   is FLAT and band-qualified (vb/cb prefixes, spct for per-frequency vectors); the
%   callers map it into the band-branched results tree (scalars/fVectors/timeVectors/
%   spectrum x VB/CB).  Bands are VB = vasomotion band (s.vFR) and CB = control band (s.cFR).
%       'bands'          scalar temporal reductions of the Group-A band series (VB & CB).
%                          The COARSE 'bands' token is shorthand for the FIVE FINE band
%                          tokens below; request a fine token directly to compute only that
%                          subset - the lean per-pixel path uses this to skip the expensive
%                          peak count.  All envelope reductions use 'omitnan' as before.
%          'bandsAmp'      vbAmpMean/vbAmpStd cbAmpMean/cbAmpStd   (VB & CB envelope mean/std)
%          'bandsSkew'     vbAmpSkew cbAmpSkew                     (VB & CB envelope skewness)
%          'bandsShape'    fCentMean/fCentStd fSprdMean/fSprdStd shapePeak  (VB instantaneous
%                            frequency-shape centroid group; shapePeak=fCentMean/fSprdMean)
%          'bandsPct'      vbAmpPct cbAmpPct                       (VB & CB envelope percentiles)
%          'bandsPeak'     nPeakMean/nPeakStd    (VB peak multiplicity; VB only.  EXPENSIVE:
%                            findpeaks per timepoint via peakCountSeries)
%       'moments'        per-frequency temporal moments spctMean/spctStd/spctSkew and the
%                          VB amplitude-percentile-resolved MEAN spectrum vbSpctPct
%                          [nFreq x numel(s.pcts)-1].  (Was the 'spectrum' level.)
%       'series'         the Group-A time series (axis = timeVSM): ts (VB envelope
%                          a_vFR(t)), tsc (CB envelope a_cFR(t)), and the VB-only
%                          fCent/fSprd/nPeak instantaneous frequency-shape series.
%       'clustering'     independent per-band Otsu-elbow flare/silence segmentation:
%                          VB from ts -> flare, vbSpctFlare/vbSpctSilence,
%                          vbDurFlareMean/Std vbDurSilenceMean/Std,
%                          vbAmpFlareMean/Std vbAmpSilenceMean/Std;
%                          CB from tsc (SAME procedure, independent mask) -> cbFlare,
%                          cbSpctFlare/cbSpctSilence and the cb* duration/amplitude
%                          scalars.  cb amplitude scalars are CB amplitude inside CB's
%                          OWN clusters (unlike the retired adFCFR/adSCFR).
%       'reconstruction' band-limited reconstruction rData (icwt over s.vFR, rescaled).
%       'spectrum'       the AMPLITUDE m.amp and PHASE m.phase of the time-frequency grid
%                          [nFreq x nDec] (was the complex 'complex' level m.WT).  When a
%                          decimated grid is requested (s.tgtFS set) the RECOVERED amplitude
%                          and phase are decimated - amplitude by a moving average, phase by
%                          the circular mean of the unit phasors - NOT the complex
%                          coefficients, so fast-oscillation amplitude is not cancelled by
%                          phase rotation within the averaging window.  At full resolution
%                          (s.tgtFS empty) m.amp=abs / m.phase=angle of the interpolated
%                          complex grid.
%     s.segVsmReturn is the sole selector; when it is absent the complete set is used.
%
%     Returns a FLAT bag `m` holding only the selected levels' quantities, plus a
%     logical m.valid (false when the normalised series is non-finite, so callers
%     keep their preallocated defaults).  m.rData (reconstruction) is returned even
%     when invalid (rescaling the zero reconstruction by the NaN centre reproduces the
%     old post-loop NaN column).
%
%   All retained quantities use the identical expression as the previous core and are
%   bit-identical to it (VB/CB envelope mean/std/skewness/percentiles, the per-frequency
%   moments, the VB percentile/flare/silence MEAN spectra, the VB flare/silence duration
%   mean+std and VB amplitude-in-cluster, the fCent/fSprd/nPeak frequency-shape metrics,
%   the independent CB clustering, the reconstruction and the VB flare mask).  The
%   amplitude/phase grid m.amp/m.phase REPLACES the previous decimated complex m.WT: at
%   full resolution (s.tgtFS empty) it is a mechanical split - m.amp=abs, m.phase=angle of
%   the interpolated complex grid - but when decimated (s.tgtFS set) the amplitude and phase
%   are averaged SEPARATELY (see 'spectrum' above), a deliberate new baseline that preserves
%   fast-oscillation amplitude.  (Still dropped from the older core: kurtosis
%   kurtVFR/kurtCFR/spctKRT, bandRatio/flatness/peakRatio, the std pages of the
%   percentile/flare/silence spectra, the statFLR/statSLC count & median, and the derived
%   occupancy/burstRate/medAmpFlare.)
%
%   INPUTS
%     time    [nT x 1] frame-time vector (SETUP) setting the sampling frequency.
%     series  one raw trace (ANALYSIS, vector) - normalised internally.
%     layout  struct from SETUP mode.
%     s       parameter struct (vFR, cFR, wFR, wVPO, normalisation, normsize,
%             tgtFS, pcts, otsuMaxN, otsuElbow, nPeakProm, segVsmReturn; see
%             runVasomotion for the full description).  s.nPeakProm is the findpeaks
%             MinPeakProminence as a fraction of the per-time VB band max-min range
%             (default 0.10 = 10%); a finite scalar >= 0.
%
%   OUTPUTS
%     layout  (SETUP)    wavelet basis + frequency grid + .want for one time base.
%     m       (ANALYSIS) flat marker bag for one series (selected levels only).
%
%   DEPENDS ON
%     MATLAB Wavelet Toolbox (cwtfilterbank/wt/icwt), Image Processing Toolbox
%     (multithresh/bwareaopen) and Statistics/Signal Toolboxes
%     (skewness/prctile/findpeaks).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

function out = getVasomotionMetrics(arg1,arg2,arg3)
if nargin==2
    out=vsmSetup(arg1,arg2);          % SETUP:    (time, s)        -> layout
elseif nargin==3
    out=vsmAnalyze(arg1,arg2,arg3);   % ANALYSIS: (series, layout, s) -> m
else
    error('getVasomotionMetrics:nargin', ...
        ['Use layout=getVasomotionMetrics(time,s) for SETUP or ' ...
         'm=getVasomotionMetrics(series,layout,s) for ANALYSIS.']);
end
end

% =====================================================================
function layout=vsmSetup(tvec,s)
%vsmSetup  Wavelet filter bank and frequency grid for one time base.
%   Lifted verbatim from runVasomotion.  Builds the analytic-Morlet cwtfilterbank
%   for the sampling implied by tvec so each analysed trace uses its own sampling
%   frequency, and returns the wavelet frequencies fwt, cone-of-influence coi, the
%   vFR/cFR-augmented frequency grid f with its band indices, the sampling
%   frequency fs and the cone-of-influence-trimmed analysed time window timeVSM.
%   layout.want resolves s.segVsmReturn once so every series/pixel on this time base
%   computes the same set of levels; layout.nrmStat/nPeakProm hoist the per-series
%   invariants (centring-statistic handle, findpeaks prominence fraction) here so
%   every series/pixel reuses them instead of rebuilding them per call.
fs=1./mean(diff(tvec));
fb=cwtfilterbank('SignalLength',numel(tvec), ...
    'SamplingFrequency',fs, ...
    'Wavelet','amor', ...
    'FrequencyLimits',s.wFR, ...
    'VoicesPerOctave',s.wVPO);
[~,fwt,coi]=wt(fb,single(tvec));
f=sort(unique(cat(1,fwt(:),s.vFR',s.cFR')),'descend');
idxsVFR=f>=s.vFR(1) & f<=s.vFR(2);
idxsCFR=f>=s.cFR(1) & f<=s.cFR(2);
timeVSM=tvec(coi<0.05);
layout.fb=fb; layout.fwt=fwt; layout.coi=coi; layout.f=f;
layout.idxsVFR=idxsVFR; layout.idxsCFR=idxsCFR; layout.fs=fs;
layout.timeVSM=timeVSM;
layout.want=vsmParseReturn(s);
%per-series invariants built ONCE (bit-identical to building them per call): the
%findpeaks prominence fraction and the centring-statistic handle for (x-c)/c.
layout.nPeakProm=resolveNPeakProm(s);
layout.nrmStat=makeNrmStat(s);
end

% =====================================================================
function want=vsmParseReturn(s)
%vsmParseReturn  Resolve s.segVsmReturn into per-level / per-band compute logicals.
%   Accepted tokens:
%     COARSE levels : 'bands' 'moments' 'series' 'clustering' 'reconstruction' 'spectrum'.
%     FINE band tokens (a caller can request a SUBSET of the band scalars without paying
%       for the rest - the lean per-pixel path uses these to avoid the expensive peak count):
%         'bandsAmp'   VB+CB ampMean/ampStd
%         'bandsSkew'  VB+CB ampSkew
%         'bandsShape' fCentMean/Std fSprdMean/Std shapePeak  (VB frequency-shape group)
%         'bandsPct'   VB+CB ampPct percentiles
%         'bandsPeak'  VB nPeakMean/Std          (EXPENSIVE: findpeaks per timepoint)
%   The COARSE 'bands' token is shorthand for ALL FIVE fine tokens (full per-segment output).
%   Returns a struct with a logical field per fine band token and per other level, plus the
%   convenience summaries .anyBand (any band scalar requested) and .anyMetric (any band /
%   moments / series / clustering requested - i.e. the magnitude metric core is needed).
%     * s.segVsmReturn present -> use it (validated cell subset of the tokens above).
%     * s.segVsmReturn absent   -> documented default (the complete set of coarse levels).
fine  ={'bandsAmp','bandsSkew','bandsShape','bandsPct','bandsPeak'};
other ={'moments','series','clustering','reconstruction','spectrum'};
coarse=[{'bands'},other];
valid =[coarse,fine];
if isfield(s,'segVsmReturn') && ~isempty(s.segVsmReturn)
    sel=s.segVsmReturn;
    if ischar(sel)||isstring(sel), sel=cellstr(sel); end
    if ~iscellstr(sel)
        error('getVasomotionMetrics:segVsmReturn', ...
            's.segVsmReturn must be a cell array of level names.');
    end
    bad=~ismember(sel,valid);
    if any(bad)
        error('getVasomotionMetrics:segVsmReturn', ...
            'Unknown s.segVsmReturn level(s): %s. Valid levels: %s.', ...
            strjoin(sel(bad),', '),strjoin(valid,', '));
    end
else
    sel=coarse;
end
if ismember('bands',sel), sel=[sel(:)',fine]; end   %coarse 'bands' -> all fine band tokens
for i=1:numel(fine),  want.(fine{i}) =ismember(fine{i}, sel); end
for i=1:numel(other), want.(other{i})=ismember(other{i},sel); end
want.anyBand  =want.bandsAmp||want.bandsSkew||want.bandsShape||want.bandsPct||want.bandsPeak;
want.anyMetric=want.anyBand||want.moments||want.series||want.clustering;
end

% =====================================================================
function p=resolveNPeakProm(s)
%resolveNPeakProm  Default/validate s.nPeakProm (findpeaks MinPeakProminence fraction).
%   Fraction of the per-time VB band max-min range used as the peak-prominence
%   threshold for the N(t) multiplicity series.  Default 0.10; a finite scalar >= 0.
if ~isfield(s,'nPeakProm') || isempty(s.nPeakProm)
    p=0.10;
    return
end
p=s.nPeakProm;
if ~(isnumeric(p)&&isscalar(p)&&isfinite(p)&&p>=0)
    error('getVasomotionMetrics:nPeakProm', ...
        's.nPeakProm must be a finite scalar >= 0 (fraction of the per-time VB band max-min range).');
end
end

% =====================================================================
function m=vsmAnalyze(X,layout,s)
%vsmAnalyze  Per-series vasomotion analysis for one vector series.
want=layout.want;
nPeakProm=layout.nPeakProm;   %hoisted into SETUP (same value for every series)

% ---- normalise internally into the relative fluctuation (x-c)/c ----
% mean/median/movmean/movmedian all act per column, so this per-series centring is
% bit-identical to the callers' pre-loop whole-matrix normalisation.  layout.nrmStat
% is the centring-statistic handle built once in SETUP.
X=X(:);
nrm=layout.nrmStat(X,1);
sub=single((X-nrm)./nrm);
m.valid=all(isfinite(sub));

fb=layout.fb; fwt=layout.fwt; coi=layout.coi; f=layout.f;
vFR=s.vFR; cFR=s.cFR;

% ---- which stages do the selected levels require? ----
wantRecon    = want.reconstruction;
wantSpectrum = want.spectrum;    %amplitude/phase grid (was the complex m.WT)
wantMetrics  = want.anyMetric;   %any band scalar / moments / series / clustering
needWT       = wantRecon || wantSpectrum || wantMetrics;

% Transform once and reuse for the reconstruction, the amp/phase grid and the metrics.
wttsFull=[];
if m.valid && needWT, wttsFull=wt(fb,sub); end

% ---- band-limited reconstruction, rescaled with the series' own centre ----
% (was the callers' post-loop rescale).  Returned even for a non-finite series:
% rescaling the zero reconstruction by the NaN centre reproduces today's NaN column.
if wantRecon
    rDataRaw=zeros(numel(layout.timeVSM),1);   %double accumulator (matches old preallocation)
    if m.valid
        tmp=icwt(wttsFull,'amor',fwt,vFR);
        rDataRaw(:)=tmp(coi<0.05);
    end
    nrmResc=nrm;
    if size(nrmResc,1)>1, nrmResc=nrmResc(coi<0.05,:); end   %moving centre -> align with recon samples
    m.rData=rDataRaw.*nrmResc + nrmResc;
end
if ~m.valid, return; end
if ~needWT, return; end   %nothing selected needs the transform

wtts=wttsFull(:,coi<0.05);          %cone-of-influence-trimmed complex coefficients

% ---- (spectrum) amplitude/phase grid m.amp/m.phase (was the complex m.WT) ----
% interp1 on the COMPLEX coefficients first (identical to the old complex path), then
% split into amplitude and phase.  When a decimated grid is requested (s.tgtFS set) the
% RECOVERED amplitude and phase are decimated - amplitude by a moving average, phase by
% the circular mean of the unit phasors - NOT the complex coefficients, so fast-oscillation
% amplitude is not cancelled by phase rotation within the averaging window.  Same subsample
% stride as the old complex path, so the decimated time axis (timeDWT) is unchanged.
if wantSpectrum
    wttsC=interp1(fwt,wtts,f);                                  %full-res complex grid on f
    amp=abs(wttsC); phase=angle(wttsC);
    if ~isempty(s.tgtFS)
        avgN=floor(layout.fs./s.tgtFS/2)*2+1;
        amp  =subsample(movmean(amp,avgN,2,'Endpoints','discard'),avgN);
        u    =subsample(movmean(exp(1i*phase),avgN,2,'Endpoints','discard'),avgN); %circular mean
        phase=angle(u);
    end
    m.amp=amp; m.phase=phase;                                   %[nFreq x nDec] (full-res: nDec=nT)
end

% ---- metric core (bands/moments/series/clustering) on the interpolated |CWT| ----
if wantMetrics
    wttsMag=interp1(fwt,abs(wtts),f);
    mm=vsmSpectrumMetrics(wttsMag,f,layout.idxsVFR,layout.idxsCFR, ...
        vFR,cFR,s.pcts,layout.fs,s.otsuMaxN,s.otsuElbow,nPeakProm,want);
    fnm=fieldnames(mm);
    for ii=1:numel(fnm), m.(fnm{ii})=mm.(fnm{ii}); end
end
end

% =====================================================================
function nrmStat=makeNrmStat(s)
%makeNrmStat  Centring-statistic handle for the relative fluctuation (mirrors callers).
%   'mean'/'median' use a single global value per series; 'mmean'/'mmedian' use a
%   moving statistic with a centred window of length s.normsize (inf or 0 -> global).
if ~isfield(s,'normalisation') || isempty(s.normalisation)
    s.normalisation='median';
end
if ~isfield(s,'normsize') || isempty(s.normsize)
    s.normsize=inf;
end
switch s.normalisation
    case 'mean'
        nrmStat=@(x,dim) mean(x,dim,'omitmissing');
    case 'median'
        nrmStat=@(x,dim) median(x,dim,'omitmissing');
    case 'mmean'
        if isinf(s.normsize) || s.normsize==0
            nrmStat=@(x,dim) mean(x,dim,'omitmissing');
        elseif s.normsize<3
            error('s.normsize must be >=3 (or inf/0 for a global statistic) when using ''mmean''/''mmedian''.');
        else
            nrmStat=@(x,dim) movmean(x,s.normsize,dim,'omitmissing');
        end
    case 'mmedian'
        if isinf(s.normsize) || s.normsize==0
            nrmStat=@(x,dim) median(x,dim,'omitmissing');
        elseif s.normsize<3
            error('s.normsize must be >=3 (or inf/0 for a global statistic) when using ''mmean''/''mmedian''.');
        else
            nrmStat=@(x,dim) movmedian(x,s.normsize,dim,'omitmissing');
        end
    otherwise
        error('s.normalisation must be ''mean'', ''median'', ''mmean'' or ''mmedian''.');
end
end

% =====================================================================
function m = vsmSpectrumMetrics(wtts,f,idxsVFR,idxsCFR,vFR,cFR,pcts,fs,maxN,otsuElbow,nPeakProm,want)
%vsmSpectrumMetrics  Vasomotion metrics from an interpolated |CWT| magnitude spectrum.
%   Shared metric core (modular via `want`): computes only the requested band sub-parts
%   (bandsAmp/Skew/Shape/Pct/Peak) and levels (moments/series/clustering).  Every retained
%   expression is unchanged from the previous core; the band amplitude envelopes ts (VB) /
%   tsc (CB) and the VB frequency-shape series fc/sprd/nPeak are shared Group-A intermediates
%   that the scalars reduce and the 'series' level stores.  Each intermediate is computed
%   ONLY when a requested output needs it - in particular the expensive nPeak (findpeaks per
%   timepoint) is skipped unless 'bandsPeak' or 'series' is requested.  Emits a FLAT,
%   band-qualified bag.
m=struct();

% ---- which shared intermediates does the requested output need? ----
needTs    = want.bandsAmp || want.bandsSkew || want.bandsPct || want.series || want.clustering || want.moments; %VB envelope ts
needTsc   = want.bandsAmp || want.bandsSkew || want.bandsPct || want.series || want.clustering;                 %CB envelope tsc
needVpcts = want.bandsPct || want.moments;                                                                      %VB envelope pct edges
needFcSprd= want.bandsShape || want.series;                                                                     %VB centroid/spread fc/sprd
needNPeak = want.bandsPeak  || want.series;                                                                     %VB peak count (EXPENSIVE)

% ---- (moments) per-frequency temporal moments of |CWT| ----
if want.moments
    m.spctMean=mean(wtts,2)';
    m.spctStd =std(wtts,0,2)';
    m.spctSkew=skewness(wtts,0,2)';
end

% ---- Group-A base: band-integrated amplitude envelopes a_B(t) (trapz form) ----
% ts = a_vFR(t), tsc = a_cFR(t); identical expression to the previous core.
if needTs || needFcSprd || needNPeak
    AbV=wtts(idxsVFR,:); fV=f(idxsVFR);                     %VB band slice (also feeds fc/sprd/nPeak)
end
if needTs
    ts=-trapz(fV,AbV,1); ts=ts(:)./(vFR(2)-vFR(1));
end
if needVpcts, vpcts=prctile(ts,pcts); end                  %shared by vbAmpPct + vbSpctPct binning
if needTsc
    AbC=wtts(idxsCFR,:); fC=f(idxsCFR);
    tsc=-trapz(fC,AbC,1); tsc=tsc(:)./(cFR(2)-cFR(1));
end

% ---- Group-A base: VB instantaneous frequency-shape series (VB band only) ----
% centroid fc(t), spread sprd(t) and peak count nPeak(t) on the |W| band spectrum, each
% computed only when needed.  Exception guard / NaN policy: a column whose band integral is
% <=0 or non-finite (den) -> fc=sprd=NaN; a column with no peaks / a degenerate or errored
% findpeaks -> nPeak=NaN.  The descending-f trapz sign cancels in the fc/sprd ratios.
if needFcSprd
    den=trapz(fV,AbV,1);                                    %[1 x nT] (<=0 on the descending grid)
    fc=trapz(fV,fV.*AbV,1)./den;                            %centroid; sign cancels
    sprd=sqrt(trapz(fV,(fV-fc).^2.*AbV,1)./den);            %spread (Hz)
    bad=~(-den>0);                                          %band integral <=0 or non-finite -> degenerate
    fc(bad)=NaN; sprd(bad)=NaN;
    fc=fc(:); sprd=sprd(:);
end
if needNPeak
    nPeak=peakCountSeries(AbV,nPeakProm);                   %[nT x 1], NaN where no-peak/degenerate
end

% ---- (series) store the Group-A time series (axis = timeVSM) ----
if want.series
    m.ts=ts;        %VB envelope a_vFR(t)
    m.tsc=tsc;      %CB envelope a_cFR(t)
    m.fCent=fc;     %VB instantaneous centroid
    m.fSprd=sprd;   %VB instantaneous spread
    m.nPeak=nPeak;  %VB instantaneous peak count
end

% ---- (bandsAmp) VB/CB envelope mean & std ----
if want.bandsAmp
    m.vbAmpMean=mean(ts);   m.vbAmpStd=std(ts);            %== old adVFR
    m.cbAmpMean=mean(tsc);  m.cbAmpStd=std(tsc);           %== old adCFR
end

% ---- (bandsSkew) VB/CB envelope skewness ----
if want.bandsSkew
    m.vbAmpSkew=skewness(ts,0);                            %== old skewVFR
    m.cbAmpSkew=skewness(tsc,0);                           %== old skewCFR
end

% ---- (bandsPct) VB/CB envelope percentiles ----
if want.bandsPct
    m.vbAmpPct=vpcts;                                      %== old pctsVFR
    m.cbAmpPct=prctile(tsc,pcts);                          %== old pctsCFR
end

% ---- (bandsShape) VB frequency-shape centroid group ----
if want.bandsShape
    m.fCentMean=mean(fc,'omitnan');    m.fCentStd=std(fc,0,'omitnan');
    m.fSprdMean=mean(sprd,'omitnan');  m.fSprdStd=std(sprd,0,'omitnan');
    m.shapePeak=m.fCentMean./m.fSprdMean;                 %peak sharpness (how well-defined the rhythm)
end

% ---- (bandsPeak) VB peak multiplicity ----
if want.bandsPeak
    m.nPeakMean=mean(nPeak,'omitnan'); m.nPeakStd=std(nPeak,0,'omitnan');
end

% ---- (moments) VB amplitude-percentile-resolved MEAN spectrum ----
if want.moments
    nB=numel(pcts)-1;
    pctSpc=zeros(numel(f),nB);
    for ii=1:nB
        sel=ts>=vpcts(ii) & ts<(vpcts(ii+1)+eps);
        pctSpc(:,ii)=mean(wtts(:,sel),2);                  %== old spctPCTS(:,:,1) (mean page only)
    end
    m.vbSpctPct=pctSpc;
end

% ---- (clustering) independent per-band Otsu-elbow flare/silence segmentation ----
if want.clustering
    cV=clusterEnvelope(ts, wtts,fs,vFR,maxN,otsuElbow);    %VB (retained: identical to old VB path)
    cC=clusterEnvelope(tsc,wtts,fs,vFR,maxN,otsuElbow);    %CB (SAME procedure, independent mask)
    % VB (bit-identical to the pre-redesign VB markers)
    m.flare=cV.mask;
    m.vbSpctFlare=cV.spctFlare;          m.vbSpctSilence=cV.spctSilence;   %== old spctFLR/SLC(:,1)
    m.vbDurFlareMean=cV.durFlare(1);     m.vbDurFlareStd=cV.durFlare(2);   %== old statFLR(2:3)
    m.vbDurSilenceMean=cV.durSilence(1); m.vbDurSilenceStd=cV.durSilence(2); %== old statSLC(2:3)
    m.vbAmpFlareMean=cV.ampFlare(1);     m.vbAmpFlareStd=cV.ampFlare(2);   %== old adFVFR
    m.vbAmpSilenceMean=cV.ampSilence(1); m.vbAmpSilenceStd=cV.ampSilence(2); %== old adSVFR
    % CB (independent CB mask -> CB amplitude inside CB's own clusters; new semantics)
    m.cbFlare=cC.mask;
    m.cbSpctFlare=cC.spctFlare;          m.cbSpctSilence=cC.spctSilence;
    m.cbDurFlareMean=cC.durFlare(1);     m.cbDurFlareStd=cC.durFlare(2);
    m.cbDurSilenceMean=cC.durSilence(1); m.cbDurSilenceStd=cC.durSilence(2);
    m.cbAmpFlareMean=cC.ampFlare(1);     m.cbAmpFlareStd=cC.ampFlare(2);
    m.cbAmpSilenceMean=cC.ampSilence(1); m.cbAmpSilenceStd=cC.ampSilence(2);
end
end

% =====================================================================
function nPeak=peakCountSeries(Ab,prom)
%peakCountSeries  Per-column peak multiplicity N(t) with the DATA-MODEL NaN policy.
%   N(t)=numel(findpeaks(Ab(:,t),'MinPeakProminence',prom*(max-min))) per column;
%   a column with no peaks, or a degenerate/errored findpeaks, or a non-finite band
%   spectrum, yields NaN (so the omitnan scalar reductions exclude that timepoint).
nT=size(Ab,2);
nPeak=nan(nT,1);
for tt=1:nT
    col=Ab(:,tt);
    if ~all(isfinite(col)), continue; end         %non-finite column -> NaN
    rng=max(col)-min(col);
    try
        pk=findpeaks(col,'MinPeakProminence',prom*rng);
    catch
        continue                                   %degenerate/errored -> NaN
    end
    if ~isempty(pk), nPeak(tt)=numel(pk); end      %no peaks -> NaN (preset)
end
end

% =====================================================================
function c = clusterEnvelope(env,wtts,fs,vFR,maxN,otsuElbow)
%clusterEnvelope  Otsu-elbow flare/silence segmentation of one band envelope + stats.
%   The segmentation is IDENTICAL to the historical VB path (elbow on the multithresh
%   sequence, threshold on the first Otsu level, minimum-duration morphological open
%   with minT derived from vFR) and is applied independently to whichever band envelope
%   `env` is passed (VB: ts, CB: tsc).  Returns the flare mask plus per-cluster MEAN
%   spectra and duration/amplitude mean+std (the count & median are dropped).
nF=size(wtts,1);

% ---- elbow-selected Otsu threshold -> raw flare mask ----
optN=zeros(maxN,1);
for n=1:maxN, [~,optN(n)]=multithresh(env,n); end
optN=diff(optN);
optN=find(optN<otsuElbow,1,'first');
if isempty(optN)
    mask=false(size(env));
else
    minT=ceil(fs./min(vFR)); %minimum flare duration
    thr=multithresh(env,optN);
    mask=env>thr(1);
    tmp=padarray(mask,[minT,0],"replicate");
    tmp=bwareaopen(tmp,minT)>0;
    mask=tmp(minT+1:end-minT);
end

% ---- per-cluster MEAN spectra + duration/amplitude mean+std ----
spctFlare=zeros(nF,1); spctSilence=zeros(nF,1);
durFlare=[NaN,NaN];    durSilence=[NaN,NaN];      %[mean std] (s)
ampFlare=[0,0];        ampSilence=[0,0];          %[mean std]
if any(mask)
    spctFlare=mean(wtts(:,mask),2);
    ampFlare=[mean(env(mask)),std(env(mask))];
end
if any(~mask)
    spctSilence=mean(wtts(:,~mask),2);
    ampSilence=[mean(env(~mask)),std(env(~mask))];
end
if any(mask) && any(~mask)
    toggleIdx=find(diff(mask(:)'));
    internalLengths=diff(toggleIdx);
    internalStates=mask(toggleIdx(1:end-1)+1);
    durFLR=internalLengths(internalStates==1);
    durSLC=internalLengths(internalStates==0);
    durFlare=[mean(durFLR),std(durFLR)]/fs;
    durSilence=[mean(durSLC),std(durSLC)]/fs;
end

c.mask=mask;
c.spctFlare=spctFlare; c.spctSilence=spctSilence;
c.durFlare=durFlare;   c.durSilence=durSilence;
c.ampFlare=ampFlare;   c.ampSilence=ampSilence;
end

% =====================================================================
function Y=subsample(X,n)
%subsample  Keep every n-th column (the decimation stride shared by amp & phase).
%   The same stride the old complex path used on the movmean('Endpoints','discard') grid,
%   so the decimated amplitude and phase land on the identical timeDWT axis.
Y=X(:,1:n:end);
end
