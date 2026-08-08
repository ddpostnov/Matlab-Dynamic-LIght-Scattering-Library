%getBolusMetrics - Model-free bolus-transit markers of ONE trace, referenced to the input
%
%   The twenty-one markers of one unit's bolus curve, built against getNVCMetrics so the
%   two model-free steps in this library read the same way.
%
%   TWO CALL MODES, dispatched on whether the second argument is a struct:
%
%   SETUP     L = getBolusMetrics(time, inputSeries, s)
%       Resolves the baseline window, the plateau window, the level list, the slope
%       window - and MEASURES THE INPUT FUNCTION, whose moments every referenced marker
%       is a difference against.  The input is geometry here in exactly the sense the
%       stimulus is geometry in getNVCMetrics: one layout, one origin, resolved once.
%   ANALYSIS  m = getBolusMetrics(series, L, s)
%       Reduces one trace to the flat, prefix-free, all-single bag.
%
%   s.segCtthReturn SAYS WHAT IS KEPT AND NEVER WHAT IS MEASURED.  Building L.in needs
%   the input's plateau step and its 50 % time, so measuring the input through the
%   CALLER'S OWN selection made a return set without 'amplitudes' or without 'times'
%   throw "Unrecognized field name Step" from inside SETUP - naming an internal field
%   instead of the setting that was wrong, and refusing a request that is perfectly
%   meaningful, since {'times'} alone is the absolute times and nothing else.  The input
%   is geometry, and geometry is not gated by a display setting: it is ALWAYS measured
%   with the complete block set, and the selection filters only L.scalarNames - what the
%   caller gets back and, through it, what reaches the tree.  Every documented return set
%   therefore works, and the one extra trace measurement costs nothing.
%
%   WHY THE MARKERS LOOK LIKE THE NVC ONES AND WHY TWO OF THEM DO NOT.  The NVC response
%   returns to baseline, so the monotone quantity is the trace's CUMULATIVE INTEGRAL and
%   the times are its quantiles.  A bolus of a NON-CLEARING intravascular tracer does the
%   opposite: it never returns, and the curve ITSELF is the cumulative - C(t)/C(inf) is
%   the normalised input convolved with the transit-time distribution.  So the times are
%   read off the curve directly, with the same anchor-at-50-%-and-search-outward rule and
%   the same linear interpolation, for the same reason.  Everything else follows: the
%   plateau step replaces the area as the normaliser, and a blood-volume marker appears
%   where a first-pass area would have been.
%
%   THE INPUT IS NOT A DELTA FUNCTION AND EVERY ABSOLUTE TIME KNOWS IT.  A ~4 s infusion
%   convolves every curve with a ~4 s rectangle, which adds the rectangle's own mean
%   (T/2 = 2 s) to every centroid and its own variance (T^2/12 = 1.33 s^2) to every
%   width.  MEASURED on the reference recording: the arterial input's centroid is 1.949 s
%   against a whole-field centroid of 2.623 s, so 74 % of the naive "mean transit time"
%   is the infusion.  Only DIFFERENCES against the input are transit times, which is why
%   Mtt, Cth and Delay exist and why m1 and vr are lowercase.
%
%   FIVE LOWERCASE FIELDS TRAVEL WITH THE BAG AND CAN NEVER BE COLUMNS.  m.m1 and m.vr
%   are the trace's own centroid and variance - the naive, infusion-dominated ones - and
%   m.riseClean, m.orderOk and m.tailSlope are three shape checks getBolusConfidence
%   grades.  All five are spelled lowercase for the reason getNVCMetrics spells m.tPeak
%   lowercase: the naming grammar's capital-first rule is what stops an unstorable
%   quantity becoming a prefixed tree column.  An absolute centroid stored beside a
%   corrected one is a number somebody will eventually plot, and a confidence input
%   stored as a marker is a marker nobody can interpret.  L.scalarNames is what reaches
%   the tree; L.blockNames is those five plus that, and it is what getBolusConfidence is
%   handed.
%
%   THE MOMENTS ARE INTEGRALS, NOT DERIVATIVES.  See getBolusMoments.  Nothing here
%   differentiates the trace except the two slope markers, and those say their window in
%   SECONDS so the same setting holds at 4 Hz and at 100 Hz.
%
%   NO FILTERING.  The step this replaced hampel-despiked, median-filtered and
%   Savitzky-Golay smoothed every trace before it measured anything.  Measured on the
%   reference recording, the single-pixel baseline noise is 13.2 counts against a plateau
%   step of 1199, and 100 noise realisations move T10/T50/T90 by 0.016/0.017/0.026 s - so
%   the filtering bought nothing that the integrals do not already buy, and it is what
%   made the per-pixel pass expensive.  The two slope markers carry their own window
%   because a derivative genuinely needs one.
%
%   THE MARKERS
%     'levels'
%       Bl      mean over the baseline window - the pre-bolus level
%       BlStd   SD over the baseline window - the noise floor every factor is scaled on
%       Fn      mean over the plateau window
%       FnStd   SD over the plateau window
%       Pk      max(y) - Bl, in the trace's own units
%     'amplitudes'
%       Step    Fn - Bl, the PLATEAU STEP.  For a non-clearing intravascular tracer the
%               long-time level is the systemic concentration times the local blood
%               volume, so this is the volume-like quantity - and it is what the first
%               pass would have measured if the first pass separated
%       PkRel   Pk / Step, the overshoot.  Dimensionless and comparable BETWEEN
%               recordings: 1.33 for the whole field of the reference recording, 1.41
%               for its arterial input
%       BvRel   Step / Step of the input - relative blood volume, the honest replacement
%               for rCBV.  It carries an unknown per-pixel optical gain (depth,
%               scattering, out-of-plane haze) and is a RELATIVE map, never an absolute
%               volume - and until background removal has run it reads a wall at 0.42
%               against a lumen's 0.45, which is scatter rather than blood
%     'times'   absolute, on the recording clock - NOT comparable between recordings
%       T0      onset, the 25-75 % rise line back-extrapolated to the baseline.  A LINE
%               FIT, not a crossing, which is why it survives where the old three-way
%               fallback (first local minimum, else first sample below the baseline
%               median, else the start of a supra-threshold run) did not: those three
%               rules find three different things and the trace decides which
%       T<p>    one per s.ctthLevelPcts entry, the time the curve reaches p % of its own
%               plateau step, anchored at 50 % and searched outward
%       TPeak   the time of the maximum.  IT NAMES THE SAME SAMPLE AS Pk
%     'transit' input-referenced - THE ONLY MARKERS COMPARABLE BETWEEN RECORDINGS
%       Delay   T50 - T50 of the input.  The transit delay
%       Mtt     m1 - m1 of the input.  The first-moment mean transit time, with the
%               infusion divided out because it is in both terms
%       Cth     sqrt(vr - vr of the input).  The width of the transit-time distribution,
%               the H in CTTH.  NaN WHEN THE DIFFERENCE IS NEGATIVE, which is what a
%               curve narrower than its own input means: nothing was resolved.  It is
%               additionally gated for the whole recording by getBolusCthFloor
%     'shape'
%       EdgeWid T90 - T10, the leading-edge width.  Immune to a short record because the
%               edge is over long before the record ends - but it is NOT a corrected
%               width: subtracting the input's edge in quadrature over-reads the true
%               kernel sd by 50-86 %, measured, so it compares units WITHIN one
%               recording and nothing else
%       RiseUp  steepest rise on the leading edge, as a fraction of Step per second
%       FallDn  steepest fall after the peak, same units, negative.  It measures
%               redistribution to the systemic level and NOT washout - the tracer does
%               not clear
%
%   WHAT IS NOT IN THE SET, AND WHY - each of these was a candidate in the brief:
%     FIRST-PASS AREA / rCBV.  Author's decision, 07-Aug-2026: the recirculation overlaps
%       a 4 s infusion, the first pass never separates, and an area with no tail to cut
%       is undefined.  MEASURED: the field mean settles at 75 % of its own peak and the
%       only trough after the peak has a prominence of 0.9 %.  BvRel is what replaces it.
%     FWHM.  There is no half maximum on the falling side to find - the curve settles
%       ABOVE half its peak.  EdgeWid is the width that does exist.
%     RECIRCULATION ONSET.  The trough that would define it is 0.9 % prominent on the
%       field mean and is not detectable at all per pixel: a trough after the peak exists
%       in 100.0 % of segments at a median time that IS the peak, i.e. findpeaks returning
%       the first noise wiggle.  A marker whose confidence is zero everywhere is worse
%       than no marker.
%
%   AN INVALID TRACE IS NaN IN EVERY FIELD, and the bag is built NaN-complete from
%   L.blockNames before anything is computed, so that rule lives in one place.
%
% Syntax:
%    L = getBolusMetrics(time, inputSeries, s)     % SETUP    - once per recording
%    m = getBolusMetrics(series, L, s)             % ANALYSIS - one trace
%
% Inputs:
%    time        - [nT x 1] recording clock, seconds, starting at 0.
%    inputSeries - [nT x 1] the arterial input function on that clock.
%    series      - [nT x 1] one unit's trace.
%    L           - the struct SETUP returned.
%    s           - parameter struct.  Fields:
%                    ctthBaselineSec  [t1 t2] pre-bolus window; [] derives it from the
%                                     input's own onset minus ctthGuardSec
%                    ctthGuardSec     margin between the baseline window and the onset,
%                                     s (default 0.5)
%                    ctthPlateauSec   width of the end window the plateau is read over,
%                                     s (default 1)
%                    ctthLevelPcts    the levels, integers in (0,100)
%                                     (default [10 25 50 75 90])
%                    ctthSlopeSec     window the two slope markers are measured over, s
%                                     (default 0.1)
%                    segCtthReturn    cell subset of {'levels','amplitudes','times',
%                                     'transit','shape'} (default all five)
%
% Outputs:
%    L - (SETUP)    windows, levels, and .in - the input's Step, T50, m1 and vr, which
%                   is what every referenced marker is a difference against.
%    m - (ANALYSIS) flat all-single bag: exactly L.blockNames, plus the logical m.valid.
%
% Example:
%    L = getBolusMetrics(results.time, aif, s);
%    m = getBolusMetrics(results.sData(:,k), L, s);
%    fprintf('arrives %.2f s after the artery, MTT %.2f s\n', m.Delay, m.Mtt);
%
% See also: getBolusMoments, getBolusInput, getBolusConfidence, getBolusCthFloor,
%           getNVCMetrics, runCTTH
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function out = getBolusMetrics(arg1,arg2,arg3)
if nargin~=3
    error('getBolusMetrics:nargin', ...
        ['Use L=getBolusMetrics(time,inputSeries,s) for SETUP or ' ...
         'm=getBolusMetrics(series,L,s) for ANALYSIS.']);
end
% THE DISPATCH IS ON THE SECOND ARGUMENT'S TYPE, not on nargin, because both modes take
% three arguments.  A layout is a struct and an input trace is not, so it cannot be
% ambiguous - and it is named here rather than discovered by a caller.
if isstruct(arg2)
    out=btAnalyze(arg1,arg2,arg3);
else
    out=btSetup(arg1,arg2,arg3);
end
end

% =====================================================================
function L=btSetup(time,inputSeries,s)
%btSetup  Everything that depends only on the time base, the settings and the INPUT.
time=double(time(:));
L.time=time; L.nT=numel(time);
if L.nT<20
    error('getBolusMetrics:shortRecord', ...
        'A bolus of %d samples is too short to measure.',L.nT);
end
L.dt=median(diff(time));
L.span=time(end)-time(1);

L.plateauSec=optScalar(s,'ctthPlateauSec',1,'the plateau window width');
L.guardSec  =optScalar(s,'ctthGuardSec'  ,0.5,'the baseline guard');
L.slopeSec  =optScalar(s,'ctthSlopeSec'  ,0.1,'the slope window');
L.endIdx=time>=time(end)-L.plateauSec;
if ~any(L.endIdx)
    error('getBolusMetrics:emptyPlateau', ...
        'No sample falls inside the %g s plateau window.',L.plateauSec);
end

% ---- the baseline window ----------------------------------------------------
% DERIVED FROM THE INPUT, not from a fraction of the frame count.  The pre-bolus period
% is a property of the RECORDING - one window for every unit - and the input is the
% earliest thing in it, so the input's own onset is what bounds it.  The old rule (10 %
% of the frames) happens to work on the reference recording and is a coincidence of its
% span.
u=double(inputSeries(:));
if isfield(s,'ctthBaselineSec') && ~isempty(s.ctthBaselineSec)
    bw=double(s.ctthBaselineSec(:))';
    L.blIdx=time>=time(1)+min(bw) & time<=time(1)+max(bw);
else
    uEnd=mean(u(L.endIdx)); u0=median(u(time<=time(1)+0.2*L.span));
    k=find(u-u0>=0.10*(uEnd-u0),1,'first');
    if isempty(k), k=round(0.1*L.nT); end
    L.blIdx=time<=time(k)-L.guardSec;
end
if nnz(L.blIdx)<5
    error('getBolusMetrics:emptyBaseline', ...
        ['The pre-bolus window holds %d samples. The bolus span starts too close to ' ...
         'the injection - the entry step must keep more of the recording before it.'], ...
        nnz(L.blIdx));
end
% AND IT MUST ACTUALLY BE PRE-BOLUS.  A count of samples is not the test: on a span that
% starts mid-rise the onset search anchors on an already-rising median, lands late, and
% hands back a "baseline" window that is a piece of the upslope - non-empty, plausible,
% and wrong in every marker that divides by Bl.  MEASURED: a record cut at 2.45 s, i.e.
% after the bolus had begun, passed the count test and scored 1.00 on every confidence
% factor.  So the window is checked against the input it was derived from.
uBl=mean(u(L.blIdx)); uEndL=mean(u(L.endIdx));
L.blRise=(max(u(L.blIdx))-uBl)/max(uEndL-uBl,eps);
if L.blRise>0.05
    error('getBolusMetrics:baselineNotFlat', ...
        ['The pre-bolus window already carries %.0f %% of the bolus. The recording ' ...
         'starts after the injection - the entry step must keep more of the ' ...
         'recording before it.'],100*L.blRise);
end

L.pcts=btParsePcts(s);
L.want=btParseReturn(s);
[L.scalarNames,L.blockNames]=btNames(L.want,L.pcts);

% ---- MEASURE THE INPUT with a bare layout, then hang it off the real one ----
% The input's own referenced markers are zero by construction, so it is analysed against
% a layout whose .in is neutral.  That is one code path measuring both, which is what
% stops the input and the units being measured two slightly different ways.
%
% AND WITH THE COMPLETE BLOCK SET, whatever the caller asked to keep.  Step and the 50 %
% time are read off mu below to build L.in, so a selection reaching this measurement
% would decide whether the ORIGIN can be computed - see the header.
L0=L; L0.want=btAllBlocks();
[L0.scalarNames,L0.blockNames]=btNames(L0.want,L.pcts);
L0.in=struct('Step',1,'T50',0,'m1',0,'vr',0);
mu=btAnalyze(u,L0,s);
L.in=struct('Step',double(mu.Step),'T50',double(mu.(sprintf('T%d',nearest50(L.pcts)))), ...
            'm1',double(mu.m1),'vr',double(mu.vr));
L.inputTrace=u;
end

% =====================================================================
function m=btAnalyze(series,L,s) %#ok<INUSD>
%btAnalyze  The markers of ONE trace.  A pure function of (series, layout).
y=double(series(:));
m=struct('valid',numel(y)==L.nT && all(isfinite(y)));
for k=1:numel(L.blockNames), m.(L.blockNames{k})=single(NaN); end
if ~m.valid, return, end
t=L.time;

% ---- levels ----------------------------------------------------------------
yb=y(L.blIdx); ye=y(L.endIdx);
Bl=mean(yb); Fn=mean(ye);
m=put(m,'Bl',Bl); m=put(m,'BlStd',std(yb));
m=put(m,'Fn',Fn); m=put(m,'FnStd',std(ye));
d=y-Bl;
[pkv,ipk]=max(d);
m=put(m,'Pk',pkv);
m=put(m,'TPeak',t(ipk));                 % NAMES THE SAME SAMPLE AS Pk

Step=Fn-Bl;
m=put(m,'Step',Step);
if Step>0
    m=put(m,'PkRel',pkv/Step);
    m=put(m,'BvRel',Step/L.in.Step);
end

% ---- the moments, by integration (see getBolusMoments) ---------------------
o=getBolusMoments(t,d,L.plateauSec);
m.m1=single(o.m1); m.vr=single(o.var);
m=put(m,'Mtt',o.m1-L.in.m1);
dv=o.var-L.in.vr;
% A NEGATIVE DIFFERENCE IS NaN AND NOT A SMALL NUMBER.  It says the curve is narrower
% than the input it is referenced to, i.e. the record resolved no width at all - which is
% the honest answer on a 13.6 s span and must not be rounded up to zero.
if dv>0, m=put(m,'Cth',sqrt(dv)); end

% ---- the times -------------------------------------------------------------
Tq=nan(1,numel(L.pcts));
if Step>0
    i50=find(d>=0.5*Step,1,'first');
    if ~isempty(i50) && i50>=2
        for k=1:numel(L.pcts)
            p=L.pcts(k); lev=0.01*p*Step;
            if p<50
                kLo=find(d(1:i50)<=lev,1,'last');
            else
                kHi=i50-1+find(d(i50:end)>=lev,1,'first');
                kLo=kHi-1;
            end
            if isempty(kLo)||kLo<1||kLo>=L.nT, continue, end
            Tq(k)=interpAt(t,d,kLo,lev);
        end
    end
end
for k=1:numel(L.pcts), m=put(m,sprintf('T%d',L.pcts(k)),Tq(k)); end
p50=nearest50(L.pcts);
m=put(m,'Delay',Tq(L.pcts==p50)-L.in.T50);

% ---- onset: the 25-75 % rise line, back-extrapolated ----------------------
% A LINE FIT RATHER THAN A CROSSING, so a single noisy sample near the foot cannot move
% it.  The old step's three-way fallback tried three different rules in order and
% reported whichever fired; this reports one quantity.
lo=0.25*Step; hi=0.75*Step;
kA=find(d>=lo,1,'first'); kB=find(d>=hi,1,'first');
if ~isempty(kA)&&~isempty(kB)&&kB>kA+1
    pf=polyfit(t(kA:kB),d(kA:kB),1);
    if pf(1)>0, m=put(m,'T0',-pf(2)/pf(1)); end
end

% ---- shape ------------------------------------------------------------------
if isfield(m,'EdgeWid') && all(isfinite(Tq([1 end])))
    m=put(m,'EdgeWid',Tq(end)-Tq(1));
end
w=max(3,2*round(0.5*L.slopeSec/L.dt)+1);
ds=movmean(gradient(d,t),w);
if Step>0
    m=put(m,'RiseUp',max(ds(1:max(ipk,2)))/Step);
    if ipk<L.nT, m=put(m,'FallDn',min(ds(ipk:end))/Step); end
end

% ---- the three lowercase shape checks getBolusConfidence grades ------------
% riseClean: the share of the LEADING EDGE whose local slope keeps up with the rise's own
% mean slope.  A clean bolus rise is steady; motion, a partial injection and a shoulder
% all show up here and nowhere else in the set.  It is measured between T10 and the peak,
% so it says nothing about the baseline, where a flat trace wanders by construction.
%
% IT COUNTS SLOW SAMPLES, NOT REVERSALS, and the difference was measured rather than
% argued.  Counting the samples where the smoothed curve goes back DOWN separated a
% 1 s flat shoulder from an untouched rise by 0.014 (0.948 against 0.962) - nothing.
% Counting the samples whose slope falls below a quarter of the rise's mean separates the
% same pair by 0.218 (0.702 against 0.920), i.e. 16x better, because a shoulder does not
% reverse a rise, it PAUSES one.
if isfinite(Tq(1)) && ipk>2
    kA=max(2,find(t>=Tq(1),1,'first'));
    if ipk>kA
        sg=ds(kA:ipk);
        m.riseClean=single(1-mean(sg<0.25*mean(sg)));
    end
end
% orderOk: the share of adjacent requested levels whose times increase.  With one level
% requested there is no pair to check and the answer is 1 - nothing was found wrong.
if numel(Tq)<2
    m.orderOk=single(1);
elseif all(isfinite(Tq))
    m.orderOk=single(mean(diff(Tq)>0));
else
    m.orderOk=single(0);
end
% tailSlope: how much of its own plateau step the unit would still move over another
% record length.  IT IS A SLOPE AND NOT A SPREAD, and that distinction was measured
% rather than assumed: the first version of this used the SD inside the plateau window,
% which on a segment average is drift but on a single PIXEL is photon noise - so it read
% 0.000-0.496 per pixel against 0.639-0.934 per segment and turned the settling check
% into a second, worse-scaled copy of the noise check.  A least-squares slope averages
% the noise out and leaves the trend, which is the thing that actually invalidates a
% plateau-normalised marker.  It reads 0.092 per segment and 0.095 per pixel, and the
% agreement between the two unit sizes is the evidence that it measures drift.
if Step>0
    kT=t>=t(end)-3*L.plateauSec;
    tt=t(kT); dd=d(kT); tm=mean(tt); vt=sum((tt-tm).^2);
    if vt>0
        m.tailSlope=single(abs(sum((tt-tm).*(dd-mean(dd)))/vt)*L.span/Step);
    end
end
end

% =====================================================================
function tc=interpAt(t,y,kLo,lev)
%interpAt  The time y passes lev between kLo and kLo+1, linearly.  A flat bracket has no
%   crossing to interpolate and reports the sample the search landed on.
dC=y(kLo+1)-y(kLo);
if dC==0, tc=t(kLo); return, end
tc=t(kLo)+(lev-y(kLo))/dC*(t(kLo+1)-t(kLo));
end

% =====================================================================
function p=nearest50(pcts)
%nearest50  The level Delay and the input's anchor are read at.  50 when it is requested,
%   otherwise the closest one - so a user who asks for [20 40 60 80] still gets a Delay,
%   and the settings file records which level it was taken at.
[~,k]=min(abs(pcts-50)); p=pcts(k);
end

% =====================================================================
function m=put(m,name,v)
if isfield(m,name), m.(name)=single(v); end
end

% =====================================================================
function pcts=btParsePcts(s)
pcts=[10 25 50 75 90];
if isfield(s,'ctthLevelPcts') && ~isempty(s.ctthLevelPcts)
    pcts=double(s.ctthLevelPcts(:))';
end
if ~all(isfinite(pcts)) || ~all(pcts==round(pcts))
    error('getBolusMetrics:ctthLevelPcts', ...
        ['s.ctthLevelPcts must be whole percentages - the marker is named from the ' ...
         'level, and T12.5 is not a field name.']);
end
if any(pcts<=0|pcts>=100)
    error('getBolusMetrics:ctthLevelPcts', ...
        's.ctthLevelPcts must lie strictly between 0 and 100.');
end
if numel(unique(pcts))~=numel(pcts)
    error('getBolusMetrics:ctthLevelPcts','s.ctthLevelPcts must not repeat a level.');
end
pcts=sort(pcts);
end

% =====================================================================
function blocks=btBlockList()
%btBlockList  The five return blocks, in the order btScalarNames writes them.  ONE list,
%   because a second spelling of it is what lets a block exist in the parser and not in
%   the names, or the other way about.
blocks={'levels','amplitudes','times','transit','shape'};
end

% =====================================================================
function want=btParseReturn(s)
levels=btBlockList();
sel=levels;
if isfield(s,'segCtthReturn') && ~isempty(s.segCtthReturn)
    sel=s.segCtthReturn;
    if ischar(sel)||isstring(sel), sel=cellstr(sel); end
    bad=~ismember(sel,levels);
    if any(bad)
        error('getBolusMetrics:segCtthReturn', ...
            'Unknown level(s): %s.',strjoin(sel(bad),', '));
    end
end
for i=1:numel(levels), want.(levels{i})=ismember(levels{i},sel); end
end

% =====================================================================
function want=btAllBlocks()
%btAllBlocks  Every block, which is what the INPUT is measured with however little of it
%   the caller keeps.
levels=btBlockList();
for i=1:numel(levels), want.(levels{i})=true; end
end

% =====================================================================
function [scalarNames,blockNames]=btNames(want,pcts)
%btNames  What is STORED, and the bag the analysis builds - those names plus the five
%   lowercase fields that travel with it and can never be columns.  The pair is built
%   here so the two layouts of one SETUP cannot disagree about the second list.
scalarNames=btScalarNames(want,pcts);
blockNames =[scalarNames,{'m1','vr','riseClean','orderOk','tailSlope'}];
end

% =====================================================================
function names=btScalarNames(want,pcts)
names={};
if want.levels,     names=[names,{'Bl','BlStd','Fn','FnStd','Pk'}]; end
if want.amplitudes, names=[names,{'Step','PkRel','BvRel'}]; end
if want.times
    names=[names,{'T0'},arrayfun(@(p)sprintf('T%d',p),pcts,'UniformOutput',false), ...
        {'TPeak'}];
end
if want.transit,    names=[names,{'Delay','Mtt','Cth'}]; end
if want.shape,      names=[names,{'EdgeWid','RiseUp','FallDn'}]; end
end

% =====================================================================
function v=optScalar(s,name,dflt,what)
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=double(s.(name)); end
if ~(isscalar(v)&&isfinite(v))
    error('getBolusMetrics:badSetting', ...
        's.%s must be a finite scalar (%s).',name,what);
end
end
