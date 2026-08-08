%getNVCMetrics - Model-free per-epoch neurovascular metrics (setup + per-epoch analysis)
%
%   getNVCMetrics centralises the model-free neurovascular-coupling "science" that
%   runNVC applies to every segment of every stimulus epoch.  It builds everything
%   that depends only on the epoch time base ONCE, and reuses it for every trace on
%   that base.  Two call modes are dispatched by argument count:
%
%   SETUP  (nargin 2)
%       layout = getNVCMetrics(time, s)
%     Resolves the baseline, response and finale windows, the stimulus boxcar, the
%     median-filter length IN SAMPLES, the threshold rules and s.segNvcReturn.
%     Returns a small `layout` built once per (signal, time base).
%
%   ANALYSIS  (nargin 3)
%       m = getNVCMetrics(series, layout, s)   % series = ONE EPOCH of one segment
%     Reduces one epoch of one trace to the flat, prefix-free, all-single scalar bag
%     of the per-trial NVC data model.  The caller applies the ns / nd / ng branch
%     prefix, exactly as runPulsatility applies ps / pd.
%
%   THERE IS NO MODEL HERE.  Every number below is read off the trace.  The optional
%   response fit is fitNVC's job and runs beside this core, never inside it.
%
%   THE FILTER IS FOR THE THRESHOLD CROSSINGS, AND FOR NOTHING ELSE.  `TRise`, `TDec`
%   and `TDecF` are level crossings, where one noisy sample would otherwise end a
%   response, and a median filter is the right instrument for exactly that.  Everything
%   else is read off the RAW trace, and the reason is that PEAK AND TTP MUST NAME THE
%   SAME SAMPLE OF THE DATA.  A peak read off the filtered trace is the median of a
%   window - some neighbouring sample's value, reported at this sample's time - so the
%   pair (`Peak`, `TTP`) would describe a point the recording does not contain.  The
%   areas are integrated on the raw trace for the same reason and one more: `trapz`
%   averages the noise out on its own, and filtering it would distort the area at the
%   window edges.
%
%   THE PER-EPOCH PEAK IS THEREFORE A BIASED ESTIMATOR, AND THAT IS WHAT `SNR` IS FOR.
%   A max over the ~80-sample response window of ONE epoch sits about 2.4 baseline SDs
%   above the baseline even when no stimulus is in the window at all - measured on the
%   reference recording, against 3.96 for a real epoch.  The bias is a property of the
%   max, not of the noise level, and what removes it is AGGREGATION, not smoothing: the
%   epoch aggregate is a MEDIAN over epochs, and `SNR` travels beside every amplitude so
%   a reader can see how far above that null each epoch actually sat.  Filtering the
%   trace first would make the number smaller without making it a measurement of
%   anything - it would just be a different sample, reported at the wrong time.
%
%   THE MEDIAN-FILTER LENGTH IS IN SECONDS.  s.nvcMedFiltSec is converted here to an
%   ODD SAMPLE COUNT using this signal's own dt.  The library analyses signals on two
%   very different clocks in the same step - the segment traces at 4 Hz and the guided
%   per-region trace at 100 Hz - so a length the user set in samples would mean two
%   different physical filters in one run.  layout.medRealisedSec reports the length
%   actually used, which is the requested one rounded to the nearest odd multiple of
%   dt.
%
%   POLARITY IS RESOLVED ONCE PER SEGMENT, NEVER PER EPOCH.  Constricting vessels and
%   the negative surround response are real, so the sign of the response is a setting
%   (s.nvcPolarity 'auto' | 'positive' | 'negative').  With 'auto' the producer
%   resolves the sign from the SEGMENT's epoch-average and passes it in
%   s.nvcPolarityResolved (+1/-1) for every epoch of that segment.  A sign resolved
%   per epoch would flip on noise in a non-responding segment, which makes `Peak` the
%   maximum of |noise| - exactly the null distribution the bias rule above exists to
%   keep out.  A lone ANALYSIS call with no resolved sign falls back to resolving from
%   the one epoch it has; that is a convenience for a single-epoch caller and a test,
%   and it is NOT what the producer should do.  An explicit 'positive'/'negative'
%   always wins over both.  The sign actually used is reported as `PeakSign`.
%
%   A THRESHOLD CROSSING MUST HOLD.  Every crossing rule requires s.nvcHoldSamples
%   CONSECUTIVE samples at or beyond the level, so a one-sample noise dip is not a
%   response ending.  Unlike the filter length this one IS a sample count, by design:
%   it is a run-length rule about consecutive samples of the trace in front of it.
%
%   A THRESHOLD NEVER CROSSED LEAVES ITS TIME NaN AND SETS ITS FLAG.  `NoRise` and
%   `NoDec` say why, and nothing falls back to a window edge: on real data 1.0 % of
%   epochs never cross the decay threshold, and an epoch that never came back must not
%   be indistinguishable from one that came back on the last sample.
%
%   WHAT TDec MEASURES, AND WHAT IT DOES NOT.  `TDec` is the first time after the peak
%   that the trace holds at or below `Fin + FinSd` - the finale level plus its own
%   noise.  IT IS A WIDTH AT A DATA-DEFINED LEVEL, NOT "THE END OF THE RESPONSE".  On
%   real data that level sits at 43 % of the peak height and TDec lands within 0.25 s
%   of where the population response crosses it, so a `Dur` of about 3.5 s for a 5 s
%   stimulus whose peak is at 4.25 s is the right answer to the question the rule
%   asks.  Its one conditioning weakness is that `Fin + FinSd` FLOATS with the finale
%   noise from epoch to epoch, so `TDecF` - the same width taken at a fixed fraction
%   of the peak height above the finale - is emitted BESIDE it.  Neither replaces the
%   other; they are two levels of the same measurement.
%
%   THE TIMING METRICS ARE AGGREGATION INPUTS, NOT PER-EPOCH READINGS.  Measured on
%   real 4 Hz data, the fraction of each metric's variance that belongs to the SEGMENT
%   rather than to the trial is: Bl 0.95, Fin 0.99, Peak 0.61, AUCblk 0.49, AUCf 0.42,
%   PeakRel 0.39, Dur 0.31, AUCb 0.30, AUCblkRel 0.25, TRise 0.18, TTP 0.13, TDec 0.06.
%   `TTP`, `TRise`, `TDec` and `Dur` from ONE epoch are therefore dominated by trial
%   noise, and neither a longer filter nor fitting the epoch changes that - it is the
%   SNR of one repetition, not a fault in the definitions.  They are stored per epoch
%   because that is where they are computed, and they are meant to be aggregated:
%   across epochs for a per-segment number, across segments for a per-epoch number.
%   The amplitude metrics are readable per epoch; the timing metrics are not.
%
%   EVERY TIME IN THE BAG IS MEASURED FROM THE STIMULUS MARK, `TTP`, `TRise`, `TDec`
%   and `TDecF` alike, so there is one origin in the bag and `Dur` is their plain
%   difference.  The windows in the SETTINGS are on the epoch clock, which is the
%   clock the protocol is written in; the metrics are on the stimulus clock, which is
%   the clock the physiology is read in.
%
%   THIS CORE NEVER CONVERTS ANYTHING.  The NVC step runs strictly AFTER the BFI
%   conversion, so every trace reaching here is already a flow index: no 1/K^2 appears
%   in this file and nothing branches on a trace's format.  The core takes a [nT x 1]
%   of numbers and a layout.  (dvsDiameter is the one input that is not a flow index -
%   it is a length in pixels - and the core does not care, because every metric here
%   is defined on an arbitrary signal.  The caller's `nd` prefix is what records the
%   difference, which is why an nd* and an ns* column must never be pooled.)
%
%   THE SCALARS
%     Window statistics
%       Bl        mean of the RAW trace over the baseline window - the normaliser
%       BlSd      SD over the baseline window - the noise floor AND the rise level
%       BlSlope   least-squares slope over the baseline window, units/s - the drift
%                 that makes AUCb and AUCf differ (they correlate only 0.52 on real
%                 data, so subtracting the finale is not the same measurement as
%                 subtracting the baseline)
%       Fin FinSd the same two over the end-of-epoch finale window
%     Amplitude, on the RAW trace
%       Peak      extremum of (trace - Bl) over the response window, in the resolved
%                 polarity's direction, SIGNED and in the trace's own units.  It is a
%                 sample the recording actually contains, and `TTP` is when it happened
%       PeakRel   Peak / Bl.  NOT redundant with Peak: flow-index units are arbitrary,
%                 so only the relative number compares between recordings
%       TTP       time of that peak, from the stimulus mark
%       SNR       peak height in baseline SDs, measured in the resolved direction -
%                 the per-epoch trust number, without which nothing downstream can
%                 weight or reject an epoch.  Positive when the trace moved the way
%                 the polarity says it should
%       PeakSign  +1 / -1, the polarity actually used
%     Timing, on the FILTERED trace
%       TRise     walking BACK from the peak, the last time the trace held at or below
%                 Bl + BlSd for nvcHoldSamples
%       TDec      walking FORWARD from the peak, the first time it holds at or below
%                 Fin + FinSd (see above)
%       TDecF     the same, at the fixed level Fin + nvcDecayFrac*(Bl + Peak - Fin)
%       Dur       TDec - TRise, a WIDTH at the Fin + FinSd level
%       NoRise    1 when the rise threshold was never held before the peak
%       NoDec     1 when the decay threshold was never held after it
%     Areas, on the RAW trace
%       AUCb      trapz(trace - Bl)  over the response window
%       AUCf      trapz(trace - Fin) over the response window
%       AUCblk    over [TRise TDec]: the total area under the trace MINUS the
%                 trapezoid under the chord joining its two endpoints - the EXTRA
%                 volume the response delivered over what a straight
%                 baseline-to-finale transition would have.  Best conditioned of the
%                 block family
%       AUCblkRel EXACTLY total/chord: how much flow passed through against how much
%                 was supposed to.  Nothing is emitted beside it - total/chord - 1 is
%                 an affine transform of it and carries no information, and a column
%                 that carries none is a column to keep in step with forever.  A zero
%                 chord gives NaN, never Inf
%   EVERY SCALAR NAME STARTS WITH A CAPITAL, so the caller's prefixed name is a plain
%   concatenation (ns + Peak = nsPeak) and there is exactly one spelling of each
%   quantity between this bag, the saved tree and the metrics tables.
%
%   ANALYSIS LEVELS  (s.segNvcReturn; a cell subset of the three names below).
%   Default (absent or empty) is the COMPLETE set.  The per-pixel path reuses SETUP
%   with s.segNvcReturn := s.ppxNvcReturn.
%       'markers'  the window statistics and the amplitude block: Bl BlSd BlSlope Fin
%                  FinSd Peak PeakRel TTP SNR PeakSign.  This is the per-pixel set:
%                  amplitude is spatially structured, per-epoch timing is not
%       'timing'   TRise TDec TDecF Dur NoRise NoDec
%       'areas'    AUCb AUCf AUCblk AUCblkRel
%   THE LEVELS GATE WHAT IS EMITTED, NOT WHAT IS COMPUTED.  Every quantity here is
%   O(nT) arithmetic on one epoch and none of it is worth branching around; `AUCblk`
%   needs the timing crossings whether or not the timing block was asked for, so a
%   level system that skipped work would have to carry that dependency.  What the
%   levels are for is the OUTPUT: layout.scalarNames is the contract the producer
%   builds its [nSeg x nEp] arrays from, and the per-pixel path must not allocate a
%   [Y x X x nEp] array for a metric nobody asked for.
%
%   AN INVALID TRACE IS NaN IN EVERY FIELD.  m.valid = all(isfinite(series)) and the
%   right length; the bag is built NaN-complete from layout.scalarNames BEFORE
%   anything is computed, so that rule is enforced in one place and cannot be broken
%   by a branch that forgets a field.  `NoRise` and `NoDec` are single 0/1 rather than
%   logical for exactly this reason - a logical cannot be NaN, and a flag that reads
%   `false` for an epoch that was never measured is a lie the tree cannot correct.
%
% Syntax:
%    layout = getNVCMetrics(time, s)              % SETUP    - once per time base
%    m      = getNVCMetrics(series, layout, s)    % ANALYSIS - one epoch of one trace
%
% Inputs:
%    time    - [nT x 1] EPOCH time base, seconds from the epoch start (SETUP).
%    series  - [nT x 1] one epoch of one segment's trace on that base (ANALYSIS).
%    layout  - the struct SETUP returned.
%    s       - parameter struct.  Fields:
%                epochStimStartSec  when the stimulus starts within the epoch, s
%                stimDurationSec    how long the stimulus lasts, s (REQUIRED)
%                epochBaselineSec   [t1 t2] baseline window on the epoch clock, s
%                epochFinaleSec     [-t1 0] finale window, back from the epoch's LAST
%                                   SAMPLE, s (REQUIRED)
%                nvcMedFiltSec      median-filter length in SECONDS (default 1.25)
%                nvcHoldSamples     consecutive samples a crossing must hold
%                                   (default 2)
%                nvcDecayFrac       fraction-of-peak level for TDecF (default 0.1)
%                nvcPolarity        'auto' (default) | 'positive' | 'negative'
%                nvcPolarityResolved  (ANALYSIS, 'auto' only) +1/-1 resolved by the
%                                   producer from this SEGMENT's epoch-average
%                segNvcReturn       cell subset of {'markers','timing','areas'}
%                                   (default: all three); the per-pixel caller passes
%                                   s.ppxNvcReturn through this field
%
% Outputs:
%    layout  - (SETUP)    windows, boxcar, filter length, threshold rules and .want
%                         for one time base, plus .scalarNames, the ORDERED list of
%                         scalars ANALYSIS emits (the caller needs no other knowledge
%                         of the field set).
%    m       - (ANALYSIS) flat, prefix-free, all-single bag for one epoch: exactly
%                         layout.scalarNames, plus the logical m.valid.
%
% Example:
%    s.epochStimStartSec=10; s.stimDurationSec=5;
%    s.epochBaselineSec=[0 10]; s.epochFinaleSec=[-5 0];
%    layout = getNVCMetrics(results.nvc.time, s);
%    m      = getNVCMetrics(epochTrace, layout, s);
%    fprintf('peak %+.1f %% at %.2f s, SNR %.1f\n', 100*m.PeakRel, m.TTP, m.SNR);
%
% See also: fitNVC, getPulsatilityMetrics, getVasomotionMetrics, movmedian, movsum
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 05-August-2026

%------------- BEGIN CODE --------------
function out = getNVCMetrics(arg1,arg2,arg3)
if nargin==2
    out=nvcMetricsSetup(arg1,arg2);          % SETUP:    (time, s)           -> layout
elseif nargin==3
    out=nvcMetricsAnalyze(arg1,arg2,arg3);   % ANALYSIS: (series, layout, s) -> m
else
    error('getNVCMetrics:nargin', ...
        ['Use layout=getNVCMetrics(time,s) for SETUP or ' ...
         'm=getNVCMetrics(series,layout,s) for ANALYSIS.']);
end
end

% =====================================================================
function L=nvcMetricsSetup(time,s)
%nvcMetricsSetup  Everything that depends only on the epoch time base and the settings.
%   The per-epoch call is then a handful of window reductions and two crossing
%   searches, which is what makes 1940 segments x 20 epochs affordable without a fit.
time=double(time(:));
L.time=time;
L.nT  =numel(time);
if L.nT<4
    error('getNVCMetrics:shortEpoch', ...
        'An epoch of %d samples is too short to measure.',L.nT);
end
L.dt=median(diff(time));
if ~(isfinite(L.dt) && L.dt>0)
    error('getNVCMetrics:badTimeBase', ...
        'The epoch time base must increase; its median step is %g s.',L.dt);
end

% ---- stimulus geometry, read on the EPOCH clock ----------------------------
% epochStimStartSec and epochBaselineSec are measured from the START of the epoch.
% Anchoring to time(1) rather than assuming it is 0 costs one addition and cannot be
% silently wrong for a producer that hands over a shifted clock.
tStim=requireScalar(s,'epochStimStartSec','when the stimulus starts within the epoch');
D    =requireScalar(s,'stimDurationSec','how long the stimulus lasts');
if ~(D>0)
    error('getNVCMetrics:stimDuration', ...
        's.stimDurationSec must be positive, not %g.',D);
end
L.tStim=time(1)+tStim;
L.D    =D;
L.u    =double(time>=L.tStim & time<=L.tStim+D);   % the boxcar, for the producer

% ---- the three windows -----------------------------------------------------
bl=requireWindow(s,'epochBaselineSec','the pre-stimulus baseline window [t1 t2]');
L.blIdx=time>=time(1)+min(bl) & time<=time(1)+max(bl);
if ~any(L.blIdx)
    error('getNVCMetrics:emptyBaseline', ...
        'No sample of the epoch falls inside the baseline window [%g %g] s.', ...
        min(bl),max(bl));
end

% THE FINALE IS ANCHORED TO THE LAST SAMPLE, not to a nominal epoch duration: the
% core is handed a time base, not a protocol, and time(end) is the only end it can
% see.  s.epochFinaleSec is written [-t1 0], so both ends are added to time(end).
fin=requireWindow(s,'epochFinaleSec','the end-of-epoch finale window [-t1 0]');
L.finIdx=time>=time(end)+min(fin) & time<=time(end)+max(fin);
if ~any(L.finIdx)
    error('getNVCMetrics:emptyFinale', ...
        ['No sample of the epoch falls inside the finale window [%g %g] s, ' ...
         'measured back from its last sample.'],min(fin),max(fin));
end

% The response window is everything from the stimulus mark to the end of the epoch -
% the finale included, because a late excursion is a response the peak search must be
% able to see.  Stored as indices so the per-epoch call does no find().
L.respIdx=find(time>=L.tStim);
if numel(L.respIdx)<3
    error('getNVCMetrics:shortResponse', ...
        ['Only %d samples follow the stimulus start (%g s into a %g s epoch); ' ...
         'there is no response to measure.'],numel(L.respIdx),tStim,time(end)-time(1));
end

% ---- the median filter, converted from SECONDS to an odd sample count ------
L.medSec=optScalar(s,'nvcMedFiltSec',1.25,'the median-filter length in seconds');
if ~(L.medSec>0)
    error('getNVCMetrics:medFilt', ...
        's.nvcMedFiltSec must be positive (seconds), not %g.',L.medSec);
end
L.medN=oddSamples(L.medSec,L.dt,L.nT);
L.medRealisedSec=L.medN*L.dt;

% ---- the threshold rules ---------------------------------------------------
L.holdN=optCount(s,'nvcHoldSamples',2, ...
    'consecutive samples a threshold crossing must hold');
L.holdN=min(L.holdN,L.nT);
L.decayFrac=optScalar(s,'nvcDecayFrac',0.1,'the fraction-of-peak decay level');
if ~(L.decayFrac>0 && L.decayFrac<1)
    error('getNVCMetrics:decayFrac', ...
        's.nvcDecayFrac must lie strictly between 0 and 1, not %g.',L.decayFrac);
end

% ---- polarity and the emitted contract -------------------------------------
L.polarity   =resolvePolarityMode(s);
L.want       =nvcParseReturn(s);
L.scalarNames=nvcScalarNames(L.want);
end

% =====================================================================
function m=nvcMetricsAnalyze(series,L,s)
%nvcMetricsAnalyze  The model-free metrics of ONE epoch of ONE trace.
%
%   THE BAG IS BUILT NaN-COMPLETE FIRST.  Every name the layout promises is written as
%   NaN before anything is computed, so "invalid trace -> NaN everywhere" is enforced
%   in one place and the caller can index the bag by layout.scalarNames without ever
%   testing which branch ran.
y=double(series(:));
m=struct('valid',numel(y)==L.nT && all(isfinite(y)));
for k=1:numel(L.scalarNames)
    m.(L.scalarNames{k})=single(NaN);
end
if ~m.valid, return, end

% ---- the two traces: filtered for amplitude and timing, raw for the areas --
yf=movmedian(y,L.medN);

% ---- window statistics, on the RAW trace -----------------------------------
yb=y(L.blIdx);
yn=y(L.finIdx);
Bl   =mean(yb);   BlSd =std(yb);
Fin  =mean(yn);   FinSd=std(yn);
m=put(m,'Bl',Bl);        m=put(m,'BlSd',BlSd);
m=put(m,'Fin',Fin);      m=put(m,'FinSd',FinSd);
m=put(m,'BlSlope',lsSlope(L.time(L.blIdx),yb));

% ---- polarity: the segment's sign, mirrored into every rule below ----------
sg=resolvePolarity(L,s,yf,Bl);
m=put(m,'PeakSign',sg);

% ---- the peak, ON THE RAW TRACE --------------------------------------------
% rr is the change from baseline mirrored into the resolved direction, so every
% comparison below is written once, for a response that rises.  The peak is taken here
% and not on yf because `Peak` and `TTP` must name the SAME SAMPLE OF THE DATA - see the
% header.  `SNR` is then exactly the quantity the reference recording's null was
% measured with: (max over the response window - baseline mean) / baseline SD.
rr=sg*(y-Bl);
[pkH,kr]=max(rr(L.respIdx));       % pkH = peak HEIGHT in the resolved direction
kPk=L.respIdx(kr);                 % ...and its index in the full epoch
Peak=sg*pkH;                       % ...signed, back in the trace's own units
m=put(m,'Peak',Peak);
m=put(m,'TTP',L.time(kPk)-L.tStim);
if Bl~=0,   m=put(m,'PeakRel',Peak/Bl); end
if BlSd~=0, m=put(m,'SNR',pkH/BlSd);    end   % a flat baseline has no noise floor

% ---- the crossings, ON THE FILTERED TRACE ----------------------------------
% This is the one place the filter is used, and the reason it exists: a level crossing
% must not be decided by a single noisy sample.  Both rules are "the level, HELD for
% holdN consecutive samples", and a held run is found with movsum on the logical rather
% than a loop: a trailing window of holdN samples sums to holdN exactly when all of them
% are below the level, and movsum's shrinking end windows can never reach holdN, which
% is the right answer - a run with fewer than holdN samples left in the epoch is not a
% held crossing.
rm=sg*(yf-Bl);
iRise=lastHeldBefore(rm<=BlSd,L.holdN,kPk);
iDec =firstHeldAfter (sg*(yf-Fin)<=FinSd,L.holdN,kPk);
m=put(m,'NoRise',isempty(iRise));
m=put(m,'NoDec' ,isempty(iDec));
if ~isempty(iRise), m=put(m,'TRise',L.time(iRise)-L.tStim); end
if ~isempty(iDec),  m=put(m,'TDec' ,L.time(iDec) -L.tStim); end
if ~isempty(iRise) && ~isempty(iDec)
    m=put(m,'Dur',L.time(iDec)-L.time(iRise));
end

% TDecF is the same width at a FIXED level: a fraction of the peak height ABOVE the
% finale, so it does not float with the finale noise the way Fin + FinSd does.  A peak
% that never rose above the finale in the resolved direction has no such height, and a
% fraction of it would not be a level.
hF=sg*(Bl+Peak-Fin);
if hF>0
    iDecF=firstHeldAfter(sg*(yf-Fin)<=L.decayFrac*hF,L.holdN,kPk);
    if ~isempty(iDecF), m=put(m,'TDecF',L.time(iDecF)-L.tStim); end
end

% ---- the areas, ON THE RAW TRACE -------------------------------------------
tr=L.time(L.respIdx);
yr=y(L.respIdx);
m=put(m,'AUCb',trapz(tr,yr-Bl));
m=put(m,'AUCf',trapz(tr,yr-Fin));

% The block pair is defined on [TRise TDec], so it exists only when both crossings do.
% Its LIMITS come from the filtered trace and its VALUES do not: the chord endpoints
% are the raw samples at those two times, which keeps AUCblk a raw-trace quantity
% throughout and keeps the "areas off the raw trace" rule whole.
if ~isempty(iRise) && ~isempty(iDec) && iDec>iRise
    tb=L.time(iRise:iDec);
    yblk=y(iRise:iDec);
    total=trapz(tb,yblk);
    chord=0.5*(yblk(1)+yblk(end))*(tb(end)-tb(1));
    m=put(m,'AUCblk',total-chord);
    if chord~=0, m=put(m,'AUCblkRel',total/chord); end
end
end

% =====================================================================
function m=put(m,name,v)
%put  Write one scalar IF the layout promised it.
%   The NaN-complete bag already carries exactly the requested names, so isfield is
%   the level test and there is no second list of what each level contains.
if isfield(m,name), m.(name)=single(v); end
end

% =====================================================================
function i=lastHeldBefore(below,holdN,kPk)
%lastHeldBefore  Last index at or before kPk ending a run of holdN samples `below`.
%   The trailing window [holdN-1 0] sums the sample and the holdN-1 before it, so the
%   sum reaches holdN exactly at the END of a held run - which is "the last time the
%   trace was still at baseline", the quantity TRise names.
held=movsum(double(below(1:kPk)),[holdN-1 0])>=holdN;
i=find(held,1,'last');
end

% =====================================================================
function i=firstHeldAfter(below,holdN,kPk)
%firstHeldAfter  First index at or after kPk starting a run of holdN samples `below`.
%   The leading window [0 holdN-1] sums the sample and the holdN-1 after it, so the
%   sum reaches holdN exactly at the START of a held run - which is "the first time the
%   trace had come back", the quantity TDec names.
held=movsum(double(below(kPk:end)),[0 holdN-1])>=holdN;
i=find(held,1,'first');
if ~isempty(i), i=i+kPk-1; end
end

% =====================================================================
function b=lsSlope(t,y)
%lsSlope  Least-squares slope of y against t, in units per second.
%   Written out rather than through polyfit so a one-sample window is a NaN instead of
%   a rank-deficiency warning printed once per segment per epoch.
t=t(:)-mean(t(:));
d=sum(t.^2);
if d==0, b=NaN; return, end
b=sum(t.*(y(:)-mean(y(:))))/d;
end

% =====================================================================
function sg=resolvePolarity(L,s,yf,Bl)
%resolvePolarity  The sign every threshold in this epoch is mirrored by.
%   'positive'/'negative' are an explicit choice and win over everything.  'auto' takes
%   the sign the producer resolved from this SEGMENT's epoch-average; only when there
%   is none does it fall back to this one epoch, which is the single-epoch convenience
%   documented in the header and NOT what the producer should rely on.
switch L.polarity
    case 'positive', sg=1;  return
    case 'negative', sg=-1; return
end
if isfield(s,'nvcPolarityResolved') && ~isempty(s.nvcPolarityResolved)
    sg=sign(double(s.nvcPolarityResolved(1)));
    if sg==0
        error('getNVCMetrics:polarityResolved', ...
            's.nvcPolarityResolved must be +1 or -1, not 0.');
    end
    return
end
sg=sign(mean(yf(L.respIdx)-Bl));
if sg==0, sg=1; end
end

% =====================================================================
function mode=resolvePolarityMode(s)
%resolvePolarityMode  Default/validate s.nvcPolarity.
mode='auto';
if isfield(s,'nvcPolarity') && ~isempty(s.nvcPolarity)
    mode=char(string(s.nvcPolarity));
end
if ~any(strcmp(mode,{'auto','positive','negative'}))
    error('getNVCMetrics:nvcPolarity', ...
        's.nvcPolarity must be ''auto'', ''positive'' or ''negative'', not ''%s''.',mode);
end
end

% =====================================================================
function n=oddSamples(sec,dt,nT)
%oddSamples  A length in SECONDS as an odd, in-range sample count for movmedian.
%   Odd so the filter is centred and introduces no half-sample shift into TTP; capped
%   at the epoch length so a filter longer than the data is a shorter filter rather
%   than an error.
n=max(1,round(sec/dt));
if mod(n,2)==0, n=n+1; end
if n>nT
    n=nT;
    if mod(n,2)==0, n=n-1; end
end
n=max(1,n);
end

% =====================================================================
function want=nvcParseReturn(s)
%nvcParseReturn  Resolve s.segNvcReturn into per-level emit logicals.
%   Absent OR empty gives the documented default, the complete set.  The per-pixel
%   caller passes s.ppxNvcReturn through this field, exactly as the pulsatility core
%   takes s.ppxPulsReturn through s.segPulsReturn.
levels={'markers','timing','areas'};
if isfield(s,'segNvcReturn') && ~isempty(s.segNvcReturn)
    sel=s.segNvcReturn;
    if ischar(sel)||isstring(sel), sel=cellstr(sel); end
    if ~iscellstr(sel)
        error('getNVCMetrics:segNvcReturn', ...
            's.segNvcReturn must be a cell array of level names.');
    end
    bad=~ismember(sel,levels);
    if any(bad)
        error('getNVCMetrics:segNvcReturn', ...
            'Unknown s.segNvcReturn level(s): %s. Valid levels: %s.', ...
            strjoin(sel(bad),', '),strjoin(levels,', '));
    end
else
    sel=levels;
end
for i=1:numel(levels), want.(levels{i})=ismember(levels{i},sel); end
end

% =====================================================================
function names=nvcScalarNames(want)
%nvcScalarNames  The ORDERED scalar contract for one layout.
%   ANALYSIS emits exactly these names and no others, so the producer builds its
%   [nSeg x nEp] arrays from this list and never has to know which levels ran.  The
%   list lives here, beside the code that fills it, and nowhere else.
names={};
if want.markers
    names=[names,{'Bl','BlSd','BlSlope','Fin','FinSd', ...
        'Peak','PeakRel','TTP','SNR','PeakSign'}];
end
if want.timing
    names=[names,{'TRise','TDec','TDecF','Dur','NoRise','NoDec'}];
end
if want.areas
    names=[names,{'AUCb','AUCf','AUCblk','AUCblkRel'}];
end
end

% =====================================================================
function v=requireScalar(s,name,what)
%requireScalar  A required numeric setting, named in the error rather than found
%   missing three lines later inside an arithmetic expression.
if ~isfield(s,name) || isempty(s.(name))
    error('getNVCMetrics:missingSetting', ...
        's.%s is required - %s, in seconds.',name,what);
end
v=double(s.(name));
if ~(isscalar(v)&&isfinite(v))
    error('getNVCMetrics:badSetting', ...
        's.%s must be a finite scalar (%s).',name,what);
end
end

% =====================================================================
function w=requireWindow(s,name,what)
%requireWindow  A required two-element window setting, in seconds.
if ~isfield(s,name) || numel(s.(name))~=2
    error('getNVCMetrics:missingWindow', ...
        's.%s must be %s, in seconds.',name,what);
end
w=double(s.(name)(:))';
if any(~isfinite(w))
    error('getNVCMetrics:badWindow', ...
        's.%s must be finite (%s).',name,what);
end
end

% =====================================================================
function v=optScalar(s,name,dflt,what)
%optScalar  An optional finite numeric setting with a documented default.
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=double(s.(name)); end
if ~(isscalar(v)&&isfinite(v))
    error('getNVCMetrics:badSetting', ...
        's.%s must be a finite scalar (%s).',name,what);
end
end

% =====================================================================
function n=optCount(s,name,dflt,what)
%optCount  An optional positive-integer setting with a documented default.
n=optScalar(s,name,dflt,what);
if ~(n>=1 && n==round(n))
    error('getNVCMetrics:badSetting', ...
        's.%s must be a positive integer (%s), not %g.',name,what,n);
end
end
