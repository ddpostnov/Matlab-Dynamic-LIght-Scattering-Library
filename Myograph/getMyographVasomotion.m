%getMyographVasomotion  Wavelet vasomotion analysis of myograph diameter time series
%
%   vasomotion = getMyographVasomotion(s,data,time) runs the same wavelet
%   vasomotion analysis as getVasomotion, but on in-memory myograph data and
%   without the per-pixel machinery.  data is [nFrames x nY] (each column a
%   diameter time series, e.g. one Y position along the vessel) or [nFrames x 1]
%   (a single 1-D trace); the analysis is performed along frames (dim 1).  time is
%   the [nFrames x 1] frame-time vector (s) used to set the sampling frequency.
%
%   The parameters are the SAME as getVasomotion (vFR, cFR, wFR, wVPO,
%   normalisation, normsize, tgtFS, pcts, otsuMaxN, otsuElbow, nPeakProm,
%   segVsmReturn); the LSCI-only fields ppxVsmReturn/ppxSegmentAveraging/vsmSignals
%   are ignored (Myograph has one signal and no per-pixel path).  s.segVsmReturn
%   (cell subset of {'bands','moments','series','clustering','reconstruction',
%   'spectrum'}, default all six) selects which analysis levels are computed.  The
%   metric maths matches getVasomotion because both call the shared core
%   getVasomotionMetrics, and both assemble the result with assembleVasomotionTree.
%
%   OUTPUT
%     vasomotion  ONE <VSM> tree - identical in shape to a single
%          getVasomotion RESULTS.vasomotion.<sig> sub-tree, with the segment
%          dimension replaced by the nY diameter positions.  It carries its OWN axes
%          at the root (a single unit, so nothing is hoisted to a container), then
%          the four band-branched shape branches gated by s.segVsmReturn:
%            .f          [nF x 1]  frequency grid, Hz (descending)
%            .timeWT     [nT x 1]  cone-of-influence-trimmed analysed time, s
%            .timeDWT    [nD x 1]  decimated time base for spectrum.amp/phase, s
%            .pctCenters [1 x nB]  amplitude-percentile bin centres
%            .scalars.VB / .scalars.CB  [nY x 1] ampMean/ampStd/ampSkew, ampPct
%                        ('bands'); VB ONLY fCentMean/fCentStd, fSprdMean/fSprdStd,
%                        shapePeak, nPeakMean/nPeakStd ('bands'); durFlareMean/Std,
%                        durSilenceMean/Std, ampFlareMean/Std, ampSilenceMean/Std
%                        ('clustering'; CB from its own independent Otsu mask)
%            .fVectors   [nY x nF] ampMean/ampStd/ampSkew ('moments');
%                        .fVectors.VB.ampMeanPct [nY x nF x nB] ('moments');
%                        .fVectors.VB/.CB .ampFlare/.ampSilence [nY x nF] ('clustering')
%            .timeVectors.VB [nT x nY] amp/fCent/fSprd/nPeak ('series'), rData
%                        ('reconstruction'), maskFlare ('clustering');
%            .timeVectors.CB [nT x nY] amp ('series'), maskFlare ('clustering')
%            .spectrum.amp / .spectrum.phase  [nY x nF x nD] amplitude & phase of the
%                        decimated CWT ('spectrum'); full-res (s.tgtFS empty) amp=|CWT|/
%                        phase=angle(CWT), decimated they are averaged separately (a
%                        deliberate new baseline, not abs/angle of a decimated complex grid)
%          Each per-Y field has nY rows/slices (the token `vsm` is fully retired:
%          the field, this variable, and every Myograph consumer read `vasomotion`).
%
%   DEPENDS ON
%     Basic functions/getVasomotionMetrics (shared wavelet vasomotion core; the
%     metric maths lives there - Myograph is no longer self-contained) and
%     assembleVasomotionTree (shared <VSM> tree assembler), the MATLAB Wavelet
%     Toolbox (cwtfilterbank/wt/icwt) and Statistics/Signal Toolboxes.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

function vasomotion = getMyographVasomotion(s,data,time)

% ---- parameter defaults (same names as getVasomotion) ----
if ~isfield(s,'normalisation')  || isempty(s.normalisation),  s.normalisation='median'; end
if ~isfield(s,'normsize')       || isempty(s.normsize),       s.normsize=inf; end
if ~isfield(s,'tgtFS'), s.tgtFS=[]; end
% output levels selected by s.segVsmReturn (resolved into layout.want by
% getVasomotionMetrics)

% (the centring statistic for the relative fluctuation is applied per series
%  inside getVasomotionMetrics, selected by s.normalisation/normsize)

if s.vFR(1)<s.wFR(1) || s.vFR(2)>s.wFR(2) || s.cFR(1)<s.wFR(1) || s.cFR(2)>s.wFR(2)
    error('s.vFR and s.cFR must lie within s.wFR.');
end

% ---- data / timing ----
if isvector(data), data=data(:); end
data=single(data);
time=time(:);
nY=size(data,2);

% ---- wavelet basis + frequency grid (shared SETUP; built once for this time base) ----
layout=getVasomotionMetrics(time,s);
fs=layout.fs; f=layout.f; timeVSM=layout.timeVSM;
want=layout.want;   %which analysis levels s.segVsmReturn selects
nT=numel(timeVSM); nf=numel(f);
npc=numel(s.pcts); nB=npc-1;

if want.spectrum && ~isempty(s.tgtFS) && s.tgtFS>fs
    error('Target sampling frequency cannot be higher than source sampling frequency')
end

% ---- decimated time base for spectrum.amp/phase (same even-window decimation as the core) ----
if ~isempty(s.tgtFS)
    avgN=floor(fs./s.tgtFS/2)*2+1;
    timeDWT=movmean(timeVSM,avgN,'Endpoints','discard'); timeDWT=timeDWT(1:avgN:end);
else
    timeDWT=timeVSM;
end
nD=numel(timeDWT);

% ---- preallocate (nY "segments"); nY is small, so preallocate every level and let
%      the tree assembler drop the unselected branches - the parfor writes are gated ----
%Zero-default (invalid columns keep 0): band amp mean/std and the flare/silence
%duration+amplitude clustering scalars.
[vbAmpMean,vbAmpStd,cbAmpMean,cbAmpStd, ...
 vbDurFlareMean,vbDurFlareStd,vbDurSilenceMean,vbDurSilenceStd, ...
 vbAmpFlareMean,vbAmpFlareStd,vbAmpSilenceMean,vbAmpSilenceStd, ...
 cbDurFlareMean,cbDurFlareStd,cbDurSilenceMean,cbDurSilenceStd, ...
 cbAmpFlareMean,cbAmpFlareStd,cbAmpSilenceMean,cbAmpSilenceStd]=deal(zeros(nY,1,'single'));
%NaN-default (undefined for a non-finite series): skewness, percentiles and the VB
%frequency-shape/multiplicity scalars (matches the old marker layout).
[vbAmpSkew,cbAmpSkew,fCentMean,fCentStd,fSprdMean,fSprdStd, ...
 shapePeak,nPeakMean,nPeakStd]=deal(nan(nY,1,'single'));
vbAmpPct=nan(nY,npc,'single'); cbAmpPct=nan(nY,npc,'single');
[spctMean,spctStd,spctSkew,vbSpctFlare,vbSpctSilence, ...
 cbSpctFlare,cbSpctSilence]=deal(zeros(nY,nf,'single'));
vbSpctPct=zeros(nY,nf,nB,'single');
[tsVB,tsCB,fCentTV,fSprdTV,nPeakTV,flareVB,flareCB]=deal(zeros(nT,nY,'single'));
rData=zeros(nT,nY,class(data));
[sAmp,sPhase]=deal(zeros(nY,nf,nD,'single'));   %amplitude & phase grids (was one complex WT)

% ---- broadcast level flags for the parfor; the science runs per column in getVasomotionMetrics.
%      Band scalars are gated per FINE token (bandsAmp/Skew/Shape/Pct/Peak) mirroring the core. ----
wantBandsAmp=want.bandsAmp; wantBandsSkew=want.bandsSkew; wantBandsShape=want.bandsShape;
wantBandsPct=want.bandsPct; wantBandsPeak=want.bandsPeak;
wantMoments=want.moments; wantSeries=want.series;
wantClustering=want.clustering; wantRecon=want.reconstruction; keepSpectrum=want.spectrum;

parfor i=1:nY
    m=getVasomotionMetrics(data(:,i),layout,s);
    if wantRecon
        rData(:,i)=m.rData;   %already rescaled per series (NaN for a non-finite series, as before)
    end
    if m.valid
        if wantMoments
            spctMean(i,:)=m.spctMean; spctStd(i,:)=m.spctStd; spctSkew(i,:)=m.spctSkew;
            vbSpctPct(i,:,:)=m.vbSpctPct;
        end
        if wantBandsAmp
            vbAmpMean(i)=m.vbAmpMean; vbAmpStd(i)=m.vbAmpStd;
            cbAmpMean(i)=m.cbAmpMean; cbAmpStd(i)=m.cbAmpStd;
        end
        if wantBandsSkew
            vbAmpSkew(i)=m.vbAmpSkew; cbAmpSkew(i)=m.cbAmpSkew;
        end
        if wantBandsPct
            vbAmpPct(i,:)=m.vbAmpPct; cbAmpPct(i,:)=m.cbAmpPct;
        end
        if wantBandsShape
            fCentMean(i)=m.fCentMean; fCentStd(i)=m.fCentStd; fSprdMean(i)=m.fSprdMean; fSprdStd(i)=m.fSprdStd;
            shapePeak(i)=m.shapePeak;
        end
        if wantBandsPeak
            nPeakMean(i)=m.nPeakMean; nPeakStd(i)=m.nPeakStd;
        end
        if wantSeries
            tsVB(:,i)=m.ts; tsCB(:,i)=m.tsc; fCentTV(:,i)=m.fCent; fSprdTV(:,i)=m.fSprd; nPeakTV(:,i)=m.nPeak;
        end
        if wantClustering
            vbSpctFlare(i,:)=m.vbSpctFlare; vbSpctSilence(i,:)=m.vbSpctSilence;
            cbSpctFlare(i,:)=m.cbSpctFlare; cbSpctSilence(i,:)=m.cbSpctSilence;
            vbDurFlareMean(i)=m.vbDurFlareMean; vbDurFlareStd(i)=m.vbDurFlareStd;
            vbDurSilenceMean(i)=m.vbDurSilenceMean; vbDurSilenceStd(i)=m.vbDurSilenceStd;
            vbAmpFlareMean(i)=m.vbAmpFlareMean; vbAmpFlareStd(i)=m.vbAmpFlareStd;
            vbAmpSilenceMean(i)=m.vbAmpSilenceMean; vbAmpSilenceStd(i)=m.vbAmpSilenceStd;
            cbDurFlareMean(i)=m.cbDurFlareMean; cbDurFlareStd(i)=m.cbDurFlareStd;
            cbDurSilenceMean(i)=m.cbDurSilenceMean; cbDurSilenceStd(i)=m.cbDurSilenceStd;
            cbAmpFlareMean(i)=m.cbAmpFlareMean; cbAmpFlareStd(i)=m.cbAmpFlareStd;
            cbAmpSilenceMean(i)=m.cbAmpSilenceMean; cbAmpSilenceStd(i)=m.cbAmpSilenceStd;
            flareVB(:,i)=m.flare; flareCB(:,i)=m.cbFlare;
        end
        if keepSpectrum
            sAmp(i,:,:)=m.amp; sPhase(i,:,:)=m.phase;
        end
    end
end

% ---- assemble the <VSM> tree (identity shapes: accumulators already in tree shape) ----
acc=struct('vbAmpMean',vbAmpMean,'vbAmpStd',vbAmpStd,'vbAmpSkew',vbAmpSkew,'vbAmpPct',vbAmpPct, ...
    'cbAmpMean',cbAmpMean,'cbAmpStd',cbAmpStd,'cbAmpSkew',cbAmpSkew,'cbAmpPct',cbAmpPct, ...
    'fCentMean',fCentMean,'fCentStd',fCentStd,'fSprdMean',fSprdMean,'fSprdStd',fSprdStd, ...
    'shapePeak',shapePeak,'nPeakMean',nPeakMean,'nPeakStd',nPeakStd, ...
    'vbDurFlareMean',vbDurFlareMean,'vbDurFlareStd',vbDurFlareStd, ...
    'vbDurSilenceMean',vbDurSilenceMean,'vbDurSilenceStd',vbDurSilenceStd, ...
    'vbAmpFlareMean',vbAmpFlareMean,'vbAmpFlareStd',vbAmpFlareStd, ...
    'vbAmpSilenceMean',vbAmpSilenceMean,'vbAmpSilenceStd',vbAmpSilenceStd, ...
    'cbDurFlareMean',cbDurFlareMean,'cbDurFlareStd',cbDurFlareStd, ...
    'cbDurSilenceMean',cbDurSilenceMean,'cbDurSilenceStd',cbDurSilenceStd, ...
    'cbAmpFlareMean',cbAmpFlareMean,'cbAmpFlareStd',cbAmpFlareStd, ...
    'cbAmpSilenceMean',cbAmpSilenceMean,'cbAmpSilenceStd',cbAmpSilenceStd, ...
    'spctMean',spctMean,'spctStd',spctStd,'spctSkew',spctSkew,'vbSpctPct',vbSpctPct, ...
    'vbSpctFlare',vbSpctFlare,'vbSpctSilence',vbSpctSilence, ...
    'cbSpctFlare',cbSpctFlare,'cbSpctSilence',cbSpctSilence, ...
    'tsVB',tsVB,'tsCB',tsCB,'fCentTV',fCentTV,'fSprdTV',fSprdTV,'nPeakTV',nPeakTV, ...
    'rData',rData,'flareVB',flareVB,'flareCB',flareCB,'amp',sAmp,'phase',sPhase);
[shp.sc,shp.scPct,shp.fv,shp.pc,shp.tv,shp.spec]=deal(@(A)A);   %accumulators already in tree shape
vasomotion=assembleVasomotionTree(acc,want,shp);

% ---- own axes at the tree root (single unit; nothing hoisted to a container) ----
vasomotion.f=f;
vasomotion.timeWT=timeVSM;
vasomotion.timeDWT=timeDWT;
vasomotion.pctCenters=(s.pcts(1:end-1)+s.pcts(2:end))./2;
end
