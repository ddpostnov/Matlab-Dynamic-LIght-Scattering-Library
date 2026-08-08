%getNVCMetrics - Model-free neurovascular markers of ONE epoch of ONE trace
%
%   getNVCMetrics is the whole of the neurovascular "science" runNVC applies to every
%   segment of every stimulus epoch.  THERE IS NO MODEL HERE and there is no fit
%   beside it either: the author's ruling of 06-Aug-2026 is that the step measures the
%   trace and reports what it measured.  Everything that depends only on the epoch time
%   base is built ONCE and reused for every trace on that base.  Two call modes are
%   dispatched by argument count:
%
%   SETUP  (nargin 2)
%       layout = getNVCMetrics(time, s)
%     Resolves the baseline window, the finale window, the stimulus mark and duration,
%     the peak window, the cumulative-area levels and the ordered scalar contract.  It
%     also publishes the four STIMULUS-clock anchors the confidence core grades its
%     ramps on, so the geometry of an epoch is resolved in exactly ONE place and one
%     layout serves getNVCMetrics, getNVCConfidence and the pages alike.
%
%   ANALYSIS  (nargin 3)
%       m = getNVCMetrics(series, layout, s)   % series = ONE EPOCH of one segment
%     Reduces one epoch of one trace to the flat, prefix-free, all-single scalar bag.
%     The caller applies the ns / nd / ng branch prefix, exactly as runPulsatility
%     applies ps / pd.
%
%   IT IS A PURE FUNCTION OF (series, layout, s), AND THAT IS NEW.  Nothing is averaged
%   before it is measured, and that now includes the POLARITY: the sign is sign(Auc) of
%   the epoch in hand, so the core needs no outside knowledge of the segment it came
%   from.  The old per-segment epoch-average that resolved a sign is gone, and so are
%   s.nvcPolarity, s.nvcPolarityResolved and the `PeakSign` output - sign(Auc) is
%   recoverable from a stored marker.  A sign taken from an AREA does not flip on noise
%   the way a sign taken from a maximum does, and it guarantees Cb(end) = |Auc| > 0, so
%   the cumulative construction below is well posed for constrictors too.
%
%   MIRROR ONLY WHAT NEEDS A DIRECTION.  `Pk`, `PkRel` and the cumulative `Cb` the times
%   are read off.  NOT `Bl`, `BlStd`, `Fn`, `FnStd`, `St`, `Auc` or the relatives - a
%   mean and an area carry their own sign, and mirroring them would make a constricting
%   vessel's `St` positive, which is a lie.  `Pk` is reported SIGNED, so a constrictor's
%   peak is negative.
%
%   THE TIMES ARE QUANTILES OF ONE CUMULATIVE CURVE, ANCHORED AT 50 %.  Author's rule of
%   06-Aug-2026: the baseline-referenced running integral is accumulated from the start
%   of the cut, the 50 % level is found scanning forward, and every other level is found
%   by walking OUTWARD from there - left for levels below 50 %, right for levels above.
%   The reason is that at 10 % of the delivered area a pre-stimulus noise bump is large
%   against the level, and a forward search from t(1) would place T10 before the mark
%   because of it.  THE 50 % LEVEL IS ALWAYS COMPUTED, whether or not 50 is in
%   s.nvcAreaPcts; it is the anchor, not a marker that happens to be requested.
%
%   Neither search can fail, which is a property rather than a case to guard: cumtrapz
%   makes Cb(1) = 0 and every level is > 0, so a left search always finds its bracket;
%   and i50 is the FIRST index to reach half the maximum, so the maximum's own index is
%   at or after it and a right search always finds one too.  The only NaN case is
%   cbMax <= 0 - a trace that never rose - and it leaves the levels and areas standing.
%
%   EVERY TIME IS LINEARLY INTERPOLATED between the samples that bracket it, never
%   rounded to a sample.  Author's rule for the whole module: at 4 Hz a sample is 0.25 s,
%   which is 14 % of T10's median, and the sample-versus-next-sample convention is
%   otherwise a systematic bias of exactly one sample.  A flat bracket reports the
%   sample the search landed on instead of dividing by zero.
%
%   cbMax, NOT Auc, IS THE TOTAL the levels are fractions of.  They differ wherever the
%   trace dips below baseline late in the epoch - 29 % of segments of the reference
%   recording, median difference 0.3 %.
%
%   THIS CORE CONVERTS NOTHING AND FILTERS NOTHING.  The step runs strictly AFTER the
%   BFI conversion, so no 1/K^2 appears here and nothing branches on a trace's format.
%   No marker reads a filtered trace either: the median filter went with the threshold
%   crossings it existed to serve, and s.nvcMedFiltSec / s.nvcHoldSamples are retired.
%   dvsDiameter is a length in pixels rather than a flow index and the core does not
%   care - every marker here is defined on an arbitrary signal, and the caller's `nd`
%   prefix is what records the difference, which is why an nd* and an ns* column must
%   never be pooled.
%
%   THE 14 SCALARS
%     'levels'
%       Bl      mean over the baseline window - the normaliser
%       BlStd   SD over the baseline window - the noise floor
%       Fn      mean over the finale window
%       FnStd   SD over the finale window.  Fn - Bl is the un-returned tail, which is
%               64 % of Auc on the reference recording because its 30 s ISI is shorter
%               than the response
%       Ep      mean over the WHOLE epoch - the TRUE mean, not the Bl + Auc/T
%               reconstruction (author's choice, 06-Aug-2026).  It exists for confidence
%               factor fEp and it is stored because a factor computed from something the
%               product does not contain cannot be checked afterwards
%     'amplitudes'
%       St      mean(y - Bl) over the stimulus window [tStim, tStim+D]
%       StRel   St / Bl
%       Pk      extremum of (y - Bl) over [tStim, tStim+D+g] in the epoch's own
%               direction, SIGNED and in the trace's own units
%       PkRel   Pk / Bl
%       Auc     trapz(t, y - Bl) over the WHOLE epoch.  The baseline stretch contributes
%               ~0 because Bl is its own mean (0.2 % of Auc, measured), so the whole
%               epoch is both simpler and free
%       AucRel  Auc / Bl, in SECONDS: the time integral of the fractional flow change,
%               i.e. the extra perfusion expressed as seconds of resting flow.  The one
%               marker that carries amplitude and persistence in a single number, and
%               the most valid of the three relatives (0.95)
%     'times'    one per s.nvcAreaPcts entry, named from the level
%       T10 T50 T90   the times by which 10 / 50 / 90 % of the accumulated increase has
%               been delivered, on the STIMULUS clock, so they may be negative.  A
%               negative T10 means the increase had already passed 10 % before the mark:
%               drift or motion, not a response, and confidence factor fT10 is what acts
%               on it
%   The relatives are exact derivations of a stored marker and are stored anyway, by the
%   author's decision: dividing by Bl removes the scale difference between two contrast
%   estimators, which lifts validity from 0.86 to 0.95 for all three.
%
%   EVERY SCALAR NAME STARTS WITH A CAPITAL, so the caller's prefixed name is a plain
%   concatenation (ns + StRel = nsStRel) and there is exactly one spelling of each
%   quantity between this bag, the saved tree and the metrics tables.
%
%   ONE LOWERCASE FIELD TRAVELS WITH THE BAG AND IS NEVER A COLUMN.  `m.tPeak` is the
%   time of the epoch's GLOBAL maximum, which confidence factor fPeak grades - and which
%   is how fPeak can say the peak fell OUTSIDE the peak window while `Pk` reports the
%   in-window one.  No peak time is stored anywhere (author's decision: TTP scored
%   reliability 0.06, two near-equal humps), so it is spelled lowercase, exactly as
%   `m.valid` is: a lowercase name cannot become a prefixed tree column.
%   layout.blockNames names it beside the scalars, and it is the ONLY difference between
%   layout.scalarNames (what reaches the tree) and layout.blockNames (what the producer
%   stacks into the [nUnit x nMetric x nEp] block the confidence core reads).
%
%   ANALYSIS LEVELS  (s.segNvcReturn; a cell subset of the three names above).  Default
%   (absent or empty) is the COMPLETE set.  The per-pixel path reuses SETUP with
%   s.segNvcReturn := s.ppxNvcReturn.
%   THE LEVELS GATE WHAT IS EMITTED, NEVER WHAT IS COMPUTED.  Auc sets the polarity and
%   the confidence core needs every marker, so all of them are always computed; what the
%   levels are for is the per-pixel path, which must not allocate a [Y x X x nEp] array
%   for a metric nobody asked for.  A producer that wants the tree gated and the
%   confidence computed calls SETUP TWICE - once with the complete set to measure with,
%   once with the user's subset to name the tree columns.
%
%   AN INVALID TRACE IS NaN IN EVERY FIELD.  m.valid = all(isfinite(series)) and the
%   right length; the bag is built NaN-complete from layout.blockNames BEFORE anything is
%   computed, so that rule is enforced in one place and cannot be broken by a branch that
%   forgets a field.
%
%   AN INTEGER PERCENTILE, OR A NAMED ERROR.  The scalar name is derived from the level,
%   and T12.5 is not a field name.  s.nvcAreaPcts must be integers strictly inside
%   (0,100) and unique; they are sorted on the way out so the name order is stable.
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
%                nvcPeakGraceSec    how far past the stimulus offset Pk may be found, s
%                                   (default 2; validated - the peak is inside [0, D+2]
%                                   for 95.4 % of segments)
%                nvcAreaPcts        cumulative-area levels the times are read at,
%                                   integers in (0,100) (default [10 50 90])
%                segNvcReturn       cell subset of {'levels','amplitudes','times'}
%                                   (default: all three); the per-pixel caller passes
%                                   s.ppxNvcReturn through this field
%
% Outputs:
%    layout  - (SETUP)    windows, boxcar, levels and the four stimulus-clock anchors
%                         for one time base, plus .scalarNames (the ORDERED scalars that
%                         reach the tree) and .blockNames (those plus the internal
%                         tPeak, which is what ANALYSIS returns as numbers).
%    m       - (ANALYSIS) flat, prefix-free, all-single bag for one epoch: exactly
%                         layout.blockNames, plus the logical m.valid.
%
% Example:
%    s.epochStimStartSec=10; s.stimDurationSec=5;
%    s.epochBaselineSec=[0 10]; s.epochFinaleSec=[-5 0];
%    layout = getNVCMetrics(results.nvc.time, s);
%    m      = getNVCMetrics(epochTrace, layout, s);
%    fprintf('%+.1f %% at the stimulus, half the area by %.2f s\n', ...
%        100*m.StRel, m.T50);
%
% See also: getNVCConfidence, getNVCEpochTrust, getPulsatilityMetrics,
%           getVasomotionMetrics, cumtrapz, trapz
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 06-August-2026

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
%   The per-epoch call is then a handful of window reductions, one cumtrapz and one
%   outward search per level, which is what makes 1940 segments x 20 epochs affordable.
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

% THE STIMULUS MUST FIT INSIDE THE EPOCH.  A protocol whose stimulus outlasts its own
% epoch has no finale to measure and no peak window to clamp; it is a settings mistake
% and it is named here rather than found as an empty window three reductions later.
if L.tStim<time(1) || L.tStim+D>time(end)
    error('getNVCMetrics:stimWindow', ...
        ['The stimulus window [%g %g] s does not fit inside the epoch [%g %g] s ' ...
         '(s.epochStimStartSec %g, s.stimDurationSec %g).'], ...
        L.tStim,L.tStim+D,time(1),time(end),tStim,D);
end
L.stIdx=time>=L.tStim & time<=L.tStim+D;
if ~any(L.stIdx)
    error('getNVCMetrics:stimWindow', ...
        'No sample of the epoch falls inside the %g s stimulus window.',D);
end

% ---- the two windows -------------------------------------------------------
bl=requireWindow(s,'epochBaselineSec','the pre-stimulus baseline window [t1 t2]');
L.blIdx=time>=time(1)+min(bl) & time<=time(1)+max(bl);
if ~any(L.blIdx)
    error('getNVCMetrics:emptyBaseline', ...
        'No sample of the epoch falls inside the baseline window [%g %g] s.', ...
        min(bl),max(bl));
end

% THE FINALE IS ANCHORED TO THE LAST SAMPLE, not to a nominal epoch duration: the core
% is handed a time base, not a protocol, and time(end) is the only end it can see.
% s.epochFinaleSec is written [-t1 0], so both ends are added to time(end).
fin=requireWindow(s,'epochFinaleSec','the end-of-epoch finale window [-t1 0]');
L.finIdx=time>=time(end)+min(fin) & time<=time(end)+max(fin);
if ~any(L.finIdx)
    error('getNVCMetrics:emptyFinale', ...
        ['No sample of the epoch falls inside the finale window [%g %g] s, ' ...
         'measured back from its last sample.'],min(fin),max(fin));
end

% ---- the peak window, CLAMPED to the last sample ---------------------------
% The grace period is how far past the stimulus offset the peak may be found; a grace
% that would run off the end of the epoch is a shorter window, not an error.
L.g=optScalar(s,'nvcPeakGraceSec',2,'how far past the stimulus offset Pk may be found');
if ~(L.g>=0)
    error('getNVCMetrics:peakGrace', ...
        's.nvcPeakGraceSec must be zero or positive (seconds), not %g.',L.g);
end
L.pkIdx=time>=L.tStim & time<=min(L.tStim+D+L.g,time(end));

% ---- the four STIMULUS-clock anchors, resolved HERE and nowhere else -------
% getNVCConfidence grades three ramps on these and must never re-derive them from the
% settings: one layout, one geometry.  Fs is the finale window's nominal START, not the
% first sample inside it, because the geometry belongs to the protocol.
L.aMark   =0;                                   % the stimulus mark
L.aStimEnd=D;                                   % the stimulus offset
L.aFinale =(time(end)+min(fin))-L.tStim;        % the finale window's start
L.aLast   =time(end)-L.tStim;                   % the epoch's last sample
L.span    =time(end)-time(1);                   % the epoch span, for fLoss

% ---- the cumulative-area levels and the emitted contract -------------------
L.pcts       =nvcParsePcts(s);
L.want       =nvcParseReturn(s);
L.scalarNames=nvcScalarNames(L.want,L.pcts);
L.blockNames =[L.scalarNames,{'tPeak'}];
end

% =====================================================================
function m=nvcMetricsAnalyze(series,L,~)
%nvcMetricsAnalyze  The model-free markers of ONE epoch of ONE trace.
%
%   THE BAG IS BUILT NaN-COMPLETE FIRST.  Every name the layout promises is written as
%   NaN before anything is computed, so "invalid trace -> NaN everywhere" is enforced in
%   one place and the caller can index the bag by layout.blockNames without ever testing
%   which branch ran.
%   (s is accepted for signature symmetry with the other cores but UNUSED: every
%   setting is resolved at SETUP, and every direction is resolved from the epoch in
%   hand.  That is what makes ANALYSIS a pure function of the trace and the geometry -
%   the first time this core has been one - and it is why there is no polarity field
%   left for a producer to get wrong.)
y=double(series(:));
m=struct('valid',numel(y)==L.nT && all(isfinite(y)));
for k=1:numel(L.blockNames)
    m.(L.blockNames{k})=single(NaN);
end
if ~m.valid, return, end

t=L.time;

% ---- levels ----------------------------------------------------------------
yb=y(L.blIdx);
yn=y(L.finIdx);
Bl=mean(yb);  BlStd=std(yb);
Fn=mean(yn);  FnStd=std(yn);
m=put(m,'Bl',Bl);      m=put(m,'BlStd',BlStd);
m=put(m,'Fn',Fn);      m=put(m,'FnStd',FnStd);
m=put(m,'Ep',mean(y));                 % the TRUE epoch mean, for confidence factor fEp

% ---- the area, WHICH IS ALSO THE POLARITY ---------------------------------
% sign(Auc) of the epoch in hand, and nothing else: one epoch of one trace, no outside
% knowledge.  An all-zero difference has no direction, so it is read as positive - the
% arbitrary choice a flat trace forces, and the markers it feeds are all zero anyway.
d  =y-Bl;
Auc=trapz(t,d);
sg =sign(Auc); if sg==0, sg=1; end

% ---- amplitudes ------------------------------------------------------------
% Bl, Auc and their relatives are NOT mirrored - a mean and an area carry their own
% sign.  Pk is, because a maximum needs a direction, and it is reported signed.
m=put(m,'St' ,mean(d(L.stIdx)));
m=put(m,'Auc',Auc);
r=sg*d;                                          % the change, in this epoch's direction
Pk=sg*max(r(L.pkIdx));
m=put(m,'Pk',Pk);
if Bl~=0
    m=put(m,'StRel' ,mean(d(L.stIdx))/Bl);
    m=put(m,'PkRel' ,Pk/Bl);
    m=put(m,'AucRel',Auc/Bl);
end

% ---- the peak's TIME: internal, for confidence factor fPeak ---------------
% The GLOBAL maximum of the epoch, not the in-window one, which is exactly how fPeak
% can score a peak that fell outside the window while Pk reports the in-window value.
[~,kGlobal]=max(r);
m.tPeak=single(t(kGlobal)-L.tStim);

% ---- the times: quantiles of the mirrored cumulative ----------------------
Cb=cumtrapz(t,r);
cbMax=max(Cb);
if cbMax>0
    i50=find(Cb>=0.5*cbMax,1,'first');           % THE ANCHOR - always, whether or not
    for k=1:numel(L.pcts)                        % 50 is a requested level
        p=L.pcts(k);
        lev=0.01*p*cbMax;
        if p<50
            % walk LEFT from the anchor: the last sample at or below the level.  Cb(1)
            % is 0 and lev > 0, so this always finds one, and it is strictly before i50
            % because Cb(i50) >= 0.5*cbMax > lev - so kLo+1 exists.
            kLo=find(Cb(1:i50)<=lev,1,'last');
        else
            % walk RIGHT from the anchor: the first sample at or above the level.  The
            % maximum's own index is at or after i50, so this always finds one; and
            % i50 >= 2 because Cb(1) = 0 < 0.5*cbMax - so kLo = k-1 >= 1.
            kHi=i50-1+find(Cb(i50:end)>=lev,1,'first');
            kLo=kHi-1;
        end
        m=put(m,sprintf('T%d',p),interpAt(t,Cb,kLo,lev)-L.tStim);
    end
end
end

% =====================================================================
function tc=interpAt(t,Cb,kLo,lev)
%interpAt  The time Cb passes `lev` between samples kLo and kLo+1, linearly.
%   A flat bracket has no crossing to interpolate, so it reports the sample the search
%   landed on rather than dividing by zero.
dC=Cb(kLo+1)-Cb(kLo);
if dC==0, tc=t(kLo); return, end
tc=t(kLo)+(lev-Cb(kLo))/dC*(t(kLo+1)-t(kLo));
end

% =====================================================================
function m=put(m,name,v)
%put  Write one scalar IF the layout promised it.
%   The NaN-complete bag already carries exactly the requested names, so isfield is the
%   level test and there is no second list of what each level contains.
if isfield(m,name), m.(name)=single(v); end
end

% =====================================================================
function pcts=nvcParsePcts(s)
%nvcParsePcts  Default/validate s.nvcAreaPcts - the cumulative-area levels.
%   Integers strictly inside (0,100), unique, sorted on the way out so the derived
%   scalar names come in a stable order.  The name IS the level (T10), which is why a
%   fractional level is a named error rather than a field called T12.5.
pcts=[10 50 90];
if isfield(s,'nvcAreaPcts') && ~isempty(s.nvcAreaPcts)
    pcts=double(s.nvcAreaPcts(:))';
end
if ~all(isfinite(pcts))
    error('getNVCMetrics:nvcAreaPcts','s.nvcAreaPcts must be finite.');
end
if ~all(pcts==round(pcts))
    error('getNVCMetrics:nvcAreaPcts', ...
        ['s.nvcAreaPcts must be whole percentages - the marker is named from the ' ...
         'level, and T%g is not a field name.'],pcts(find(pcts~=round(pcts),1)));
end
if any(pcts<=0 | pcts>=100)
    error('getNVCMetrics:nvcAreaPcts', ...
        's.nvcAreaPcts must lie strictly between 0 and 100; %g does not.', ...
        pcts(find(pcts<=0|pcts>=100,1)));
end
if numel(unique(pcts))~=numel(pcts)
    error('getNVCMetrics:nvcAreaPcts', ...
        's.nvcAreaPcts must not repeat a level; each one names one scalar.');
end
pcts=sort(pcts);
end

% =====================================================================
function want=nvcParseReturn(s)
%nvcParseReturn  Resolve s.segNvcReturn into per-level emit logicals.
%   Absent OR empty gives the documented default, the complete set.  The per-pixel
%   caller passes s.ppxNvcReturn through this field, exactly as the pulsatility core
%   takes s.ppxPulsReturn through s.segPulsReturn.
levels={'levels','amplitudes','times'};
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
function names=nvcScalarNames(want,pcts)
%nvcScalarNames  The ORDERED scalar contract for one layout.
%   ANALYSIS emits exactly these names (plus the lowercase tPeak and valid), so the
%   producer builds its [nSeg x nEp] arrays from this list and never has to know which
%   levels ran.  The list lives here, beside the code that fills it, and nowhere else.
names={};
if want.levels
    names=[names,{'Bl','BlStd','Fn','FnStd','Ep'}];
end
if want.amplitudes
    names=[names,{'St','StRel','Pk','PkRel','Auc','AucRel'}];
end
if want.times
    names=[names,arrayfun(@(p)sprintf('T%d',p),pcts,'UniformOutput',false)];
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
