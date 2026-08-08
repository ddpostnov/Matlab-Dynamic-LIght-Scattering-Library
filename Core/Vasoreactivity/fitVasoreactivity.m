%fitVasoreactivity - Drug / gas-challenge vasoreactivity response fit (setup + per-trace fit)
%
%   fitVasoreactivity centralises the vasoreactivity "science" that
%   runFitVasoreactivity applies to every segment trace and every masked pixel of a
%   whole-recording contrast product.  It builds everything that depends only on the
%   recording clock ONCE - the decimated time base, the injection anchor, the
%   baseline window, the bounds and the multi-start grid - and reuses it for every
%   trace on that base.  Two call modes are dispatched by argument count:
%
%   SETUP  (nargin 2)
%       layout = fitVasoreactivity(time, s)
%     Builds the decimated time base implied by `time` and s.tgtFS, places the
%     injection and the baseline window on it, and resolves s.segVrcReturn once.
%
%   ANALYSIS  (nargin 3)
%       m = fitVasoreactivity(series, layout, s)   % series = [nTraw x 1], FULL rate
%     The science for ONE trace.  The trace is handed in at the RECORDING's rate and
%     decimated here, so a caller never has to decimate a cube itself.  Returns a
%     FLAT, prefix-free bag `m` (the caller applies the rs/rd branch prefix, exactly
%     as runPulsatility applies ps/pd); every numeric field is SINGLE.
%
%   THIS STEP MODELS PHARMACOKINETICS DRIVING A QUASI-STATIC VASCULAR GAIN.  It does
%   NOT model vascular dynamics, and nothing it reports is comparable with anything
%   fitNVC reports.  Decompose the measurement:
%
%       dose -> [PK: concentration at the effect site] -> [PD: dose to tone]
%            -> [vascular dynamics] -> BFI(t)
%
%   For neurovascular coupling the input is a boxcar you control exactly, the PD is
%   near-linear, and the VASCULAR FILTER dominates the observed shape - which is why
%   fitNVC fits a kernel and convolves it.  For a drug injection this inverts: the PK
%   dominates and the vascular filter is invisible.  The numbers settle it.  The NVC
%   filter has tauS ~ 0.8 s, tauF ~ 0.4 s and a ring period of ~4 s; acetazolamide
%   CBF responses have t_max ~ 11-15 min and t_washout ~ 35-49 min (Fahlstrom et al.
%   2024, Magn Reson Imaging 110:35-42).  That is two to three orders of magnitude of
%   separation: on a drug timescale the second-order system is fully settled at every
%   instant and contributes a static gain, not a shape.
%
%   THEREFORE tauS AND tauF ARE UNIDENTIFIABLE FROM A SLOW DRUG CURVE.  Fitting them
%   anyway converges and returns numbers that mean nothing.  The rate constants here
%   are named RateFast / RateSlow and are reported in 1/s rather than as time
%   constants, deliberately, so that nobody sets a "tau" from this file beside a
%   "tau" from fitNVC.  They are not the same quantity and the numbers look
%   tantalisingly comparable.  The decision rule for which core to reach for is only:
%   IS THE INPUT FAST OR SLOW COMPARED WITH ~5 s?
%
%   THE MODEL  (s.vrcModel = 'bateman', the only model in v1)
%     A general bi-exponential (absorption/elimination) response on a drifting
%     baseline, with the bracket NORMALISED TO UNIT PEAK:
%           tp   = log(r2/r1)/(r2-r1)                  % the bracket's own peak time
%           p    = exp(-r1*tp) - exp(-r2*tp)           % its peak value
%           y(t) = b0 + b1*(t-tRef) + (A/p)*( exp(-r1*x) - exp(-r2*x) ),  x = t-t0-tInj
%     zero before the onset, with the constraint r2 > r1 > 0.
%
%     SIX free parameters: b0 b1 A t0 r1 (r2-r1).  A IS SIGNED - positive is
%     dilation, negative is constriction - and is NOT bounded positive.
%
%     THE NORMALISATION IS WHY A MEANS ANYTHING.  As written in the textbooks the
%     peak of exp(-r1*t)-exp(-r2*t) is not 1: it depends on r1 and r2, so A and the
%     rate constants are multiplicatively entangled and the optimiser walks a valley
%     trading one against the other.  This is the same trap fitNVC hits with its
%     second-order kernel (whose area is tauF) and it has the same fix.  Divided by
%     p, A IS the peak change from the baseline in the trace's own units - readable,
%     decorrelated from the rates, and directly comparable across recordings.
%
%     r2 > r1 IS FITTED AS A BOX, by taking (r2-r1) as the parameter rather than r2.
%     A box constraint cannot express an inequality between two parameters, and the
%     normalised model is exactly INVARIANT under swapping r1 and r2 (both the
%     bracket and p change sign, so their ratio does not) - so left free, the two
%     rates are identified only as an unordered pair and the reported columns would
%     flip between traces.  Fitting the gap makes RateSlow the elimination rate and
%     RateFast the absorption rate in every row.
%
%     THE DEGENERATE LIMIT IS A BRANCH, NOT AN ERROR.  As r2 -> r1 the bracket
%     collapses, and both it and p go to zero; the ratio does not.  Its limit is
%           g(x) = e*r1*x*exp(-r1*x)                   % also peak 1, at x = 1/r1
%     which is implemented, selected by a band RELATIVE to r1 (never an == test)
%     because the optimiser walks across that boundary constantly and a 0/0 there
%     poisons the whole trace.  This mirrors the critical-damping branch of fitNVC.
%
%     THE DRIFT TERM b1 IS NOT OPTIONAL.  Over 45 minutes, anaesthesia depth,
%     temperature and window clarity move an LSCI baseline by amounts comparable with
%     the response itself, and without b1 that drift is absorbed into the response.
%     Fitting drift blind is still a compromise: a VEHICLE or SHAM recording is the
%     only way to identify it properly, and a study that wants trustworthy amplitudes
%     should record one.
%
%   n = 1 PER ANIMAL, AND THAT SHAPES THE WHOLE DESIGN.  Neurovascular coupling gives
%   30-60 epochs to average and bootstrap; a drug injection gives ONE curve, so the
%   epoch-bootstrap confidence interval is simply not available.  Three consequences
%   are designed in rather than discovered later:
%     - the parameter count is held at SIX.  This is why the Lindquist-Wager
%       inverse-logit is excluded here even more firmly than in fitNVC;
%     - the uncertainty on the headline amplitude is derived from the residuals with
%       an AR(1) correction (GainSE, CVRSE), and the residual autocorrelation Rho and
%       the effective sample size NEff are stored beside them, so an over-optimistic
%       interval is VISIBLE rather than implied.  A long recording has thousands of
%       samples and tens of independent ones; NEff says which;
%     - AICc is stored per recording, so nested models can be compared on one animal.
%
%   TWO CONFOUNDS THAT HAVE NO NVC COUNTERPART, and neither is corrected here:
%     SYSTEMIC.  Drugs change mean arterial pressure, heart rate and PaCO2, and
%     cerebral autoregulation then acts to OPPOSE the resulting flow change.  The
%     observed curve is drug-on-vessel PLUS autoregulation-on-systemic-pressure.  This
%     is negligible in neurovascular coupling and potentially dominant here.  The
%     wrapper can carry a measured pressure trace alongside the fit (s.mapTrace) so
%     the confound is visible in the product and on the report page; it is NOT a
%     regressor, and reading a fitted amplitude as drug-on-vessel is an assumption.
%     ASYMMETRY.  Constrictors are not negative dilators - different receptor
%     kinetics, different smooth-muscle mechanics, hysteresis in the pressure-diameter
%     relation.  Do not pool a dilation and a constriction protocol, and do not assume
%     one kernel with a sign flip; fit them separately.
%
%   THE STEAL DIAGNOSTIC IS MANDATORY AND MODEL-FREE.
%   A SUM OF REAL EXPONENTIALS CANNOT UNDERSHOOT: it has one extremum and is monotone
%   thereafter.  An undershoot needs a complex-conjugate pole pair, which this model
%   does not have - so a recording with genuine cerebrovascular steal cannot be
%   described by it, and the fit does not fail loudly, it fails QUIETLY: it absorbs the
%   sub-baseline tail into the drift term and reports a plausible R2 with a materially
%   wrong amplitude.  Steal is therefore read off the RAW decimated trace and stored
%   beside TMax and TWash:
%       Steal = 1 when the trace, WITH THE BASELINE-WINDOW LINE REMOVED, stays below
%               -max(vrcStealK*sd, vrcStealFrac*|Peak|) for at least vrcStealSec
%               consecutive seconds AFTER the peak; sd is the spread of the baseline
%               window ABOUT that line
%   THE LINE REMOVAL IS NOT A REFINEMENT.  Measured against the baseline MEAN, as the
%   plain definition would have it, a 30-45 minute recording ends well below its own
%   baseline through drift alone - the same drift b1 exists to model - and the flag
%   fires on EVERY long recording and means nothing.  The line through the baseline
%   window is the only drift estimate available without a fit, and the vrcStealFrac
%   floor is what absorbs the error in extrapolating it across half an hour.
%   A FIT WITH Steal = 1 IS A FIT NOT TO TRUST.  That is the documented failure mode,
%   not a hypothetical, and the flag is the only thing that catches it - R2 will not.
%
%   MODELS DELIBERATELY NOT IN v1, with the decisions already taken so a later
%   session does not re-open them:
%     'fahlstrom'   the literature-exact acetazolamide form y = A*(2*exp(-r1*t) -
%                   exp(-r2*t)), three parameters, baseline structural.  IT IS TO BE
%                   FITTED IN SECONDS with the published bounds converted (the paper's
%                   r1 in [0.005 0.04] and r2 in [0 0.8] are per MINUTE), so the
%                   library stays uniformly in seconds.  Its A is the BASELINE, not
%                   the peak - it is not this file's A, and two quantities with one
%                   letter is exactly how a comparison goes wrong.
%     'expConv'     y = b0 + b1*(t-tRef) + G*(x * (1/tau)exp(-t/tau)), for a MEASURED
%                   input (PetCO2, drug concentration, MAP): tau is the dynamic
%                   component - the fitness of smooth-muscle tone control - and G the
%                   steady-state reactivity.  Scientifically the cleanest of the four
%                   and worth pushing protocols toward.  IT NEEDS A REAL NUMERICAL
%                   conv(): fitNVC has none to lend, because a boxcar input lets it
%                   convolve in closed form through the step response, and that trick
%                   does not generalise.
%     'secondOrder' fitNVC's damped-oscillator kernel driven by the fitted Bateman
%                   curve, for a reproducible sub-baseline dip.  Needs the IMPULSE
%                   response h, and fitNVC's private stepResponse returns the STEP
%                   response H - so it needs a sibling, and fitNVC re-pointed at it in
%                   the same session, or there are two spellings of one kernel.
%
%   ANALYSIS LEVELS  (s.segVrcReturn; a cell subset of the three names below).
%   Default (absent or empty) is the COMPLETE set.  The per-pixel path reuses SETUP
%   with s.segVrcReturn := s.ppxVrcReturn.
%       'markers'         model-free scalars read straight off the DECIMATED RAW
%                         trace.  Always computable, no fit:
%                           Peak      largest excursion from the baseline mean after
%                                     the injection, SIGNED (so a constrictor reports
%                                     a negative peak - this differs from fitNVC,
%                                     where a downward response is not a first-class
%                                     case and the peak is the maximum)
%                           PeakRel   the same as a fraction of the baseline mean
%                           TMax      time of that peak, from the INJECTION
%                           AUC       integral of the change over s.vrcAucSec
%                           TWash     when the trace first returns to within
%                                     s.vrcWashFrac of the peak, from the injection
%                           BlMean BlStd  baseline-window mean and standard deviation,
%                                     as measured - BlStd therefore INCLUDES any drift
%                                     across the window, which is what a reader wants
%                                     to see.  The steal test below uses the spread
%                                     about the baseline LINE instead, for the reason
%                                     given there
%                           Steal StealSec  the diagnostic above, and how long the
%                                     trace actually stayed below the threshold
%       'model'           fitted parameters, derived physiology and fit quality
%                         (runs the fit) - see the scalar list below
%       'reconstruction'  fData [nT x 1] on the DECIMATED base, the fitted curve
%
%   THE MARKERS ARE COMPUTED ON THE RAW TRACE, NOT ON THE RECONSTRUCTION - the same
%   rule fitNVC follows and for the same reason: read off the fit they would be no
%   check at all.  Peak / AUC / TWash are ALSO re-derived from the fitted curve and
%   stored with a Fit suffix, so the two are one measurement of two curves.
%
%   SCALARS EMITTED AT THE 'model' LEVEL
%     Gain Onset Baseline Drift            A, t0 (from the injection), b0, b1
%     RateFast RateSlow                    r2, r1 in 1/s - PK rates, NOT vascular taus
%     TMaxModel                            t0 + tp, the drift-free peak time from the
%                                          injection; the number comparable with the
%                                          literature's t_max.  TMax (marker) is read
%                                          off the raw trace and therefore includes
%                                          the drift, which over 45 min can move it
%     K                                    Gain/tp, the rate of rise of the first phase
%     CVR                                  Gain/Baseline, the relative change - the
%                                          reported reactivity, and the model's twin
%                                          of the model-free PeakRel
%     R2 RMSE AICc                         fit quality, on the residual
%     Rho NEff GainSE CVRSE                AR(1) residual correlation, the effective
%                                          sample size it implies, and the corrected
%                                          standard error of the two headline
%                                          amplitudes.  A 95% interval is 1.96*SE
%     NIter Converged StartsAgree          the optimiser's own diagnostics
%     PeakFit AUCFit TWashFit              the markers re-read off the fitted curve
%   EVERY SCALAR NAME STARTS WITH A CAPITAL, so the caller's prefixed name is a plain
%   concatenation (rs + Peak = rsPeak) and there is exactly one spelling of each
%   quantity between this bag, the saved tree and the metrics tables.
%   Steal and Converged are 1/0 SINGLES rather than logicals, so every scalar the
%   caller assembles into a map or a table column has one class.
%
%   AN INVALID TRACE IS NaN IN EVERY FIELD.  m.valid = all(isfinite(series)); an
%   invalid trace skips the decimation, the fit and every emitted scalar, and fData is
%   NaN.  Same clean rule as fitNVC - there is no legacy golden here.
%
%   DECIMATE BEFORE FITTING (s.tgtFS, default 1 Hz).  A 45-minute recording is tens of
%   thousands of samples per segment and a per-pixel cube that will not fit in memory,
%   and the response has a t_max of MINUTES, so 1 Hz is generous.  The decimation is
%   blockDecimate, a non-overlapping BLOCK MEAN, so it anti-aliases rather than merely
%   sub-sampling and the derived parameters do not move with the target rate; the
%   cardiac and respiratory content of an LSCI trace would otherwise fold straight down
%   into the band the drug response lives in.  s.tgtFS empty, or at or above the
%   recording's own rate, keeps every sample.  THE TRACE ARRIVES AT THE RECORDING'S OWN
%   RATE and is decimated inside ANALYSIS, so a per-pixel caller never has to decimate
%   a cube of its own.
%
%   THE FIT IS MULTI-START, AND THAT IS NOT OPTIONAL.  The objective is non-convex; a
%   single start returns quiet garbage on a minority of traces and nothing in the
%   output says so.  StartsAgree records how many start points reached the winning
%   solution - a trace where 1 start in 16 wins is a trace to distrust.  The grid is
%   Roberts' deterministic R_d sequence (quasiUniform) over a PHYSIOLOGICAL sub-box of
%   the bounds, seeded with two literature start points, so no rng state leaks between
%   traces or workers and a fit is pinnable.
%
%   TWO OF THE EMITTED SCALARS DESCRIBE THE OPTIMISER, NOT THE VESSEL.  NIter and
%   StartsAgree count what the trust region did on the way to the answer, and BLAS
%   threading differs between a parfor worker and the client, so the same optimum is
%   reached in a different number of steps.  A test that compares those two across a
%   parallel split is asserting something that is not true.
%
% Syntax:
%    layout = fitVasoreactivity(time, s)             % SETUP    - once per time base
%    m      = fitVasoreactivity(series, layout, s)   % ANALYSIS - one trace
%
% Inputs:
%    time    - [nTraw x 1] recording clock, seconds, at the product's own rate (SETUP).
%    series  - [nTraw x 1] one trace on that FULL-rate base (ANALYSIS).
%    layout  - the struct SETUP returned.
%    s       - parameter struct.  Fields:
%                injectionSec   when the drug was given, seconds from the start of
%                               this product's time base (REQUIRED; the wrapper
%                               resolves the per-file list to a scalar here)
%                baselineSec    [t1 t2] pre-injection baseline window, seconds on the
%                               same clock; must END before injectionSec (REQUIRED)
%                vrcModel       'bateman' (default; the only model in v1)
%                tgtFS          target rate before fitting, Hz (default 1)
%                vrcAucSec      AUC horizon after the injection, s (default 2700)
%                vrcStealK      steal threshold in baseline SDs (default 2)
%                vrcStealFrac   steal threshold FLOOR, as a fraction of the peak
%                               (default 0.15); the threshold is the larger of the two
%                vrcStealSec    how long it must stay below, s (default 60)
%                vrcWashFrac    washout level as a fraction of the peak (default 0.1)
%                vrcStarts      number of multi-start points (default 16)
%                segVrcReturn   cell subset of {'markers','model','reconstruction'}
%                               (default: all three); the per-pixel caller passes
%                               s.ppxVrcReturn through this field
%
% Outputs:
%    layout  - (SETUP)    the decimated time base, the injection anchor, the baseline
%                         and response windows, the bounds, the start grid and .want,
%                         plus .scalarNames, the ORDERED list of scalars ANALYSIS
%                         emits (the caller needs no other knowledge of the field set).
%    m       - (ANALYSIS) flat, prefix-free, all-single bag for one trace: exactly
%                         layout.scalarNames, plus fData when 'reconstruction' was
%                         asked for, plus the logical m.valid.
%
% Example:
%    s.injectionSec=600; s.baselineSec=[60 540]; s.tgtFS=1;
%    layout = fitVasoreactivity(results.time, s);
%    m      = fitVasoreactivity(results.sData(:,7), layout, s);
%    fprintf('CVR = %.1f %%, t_max = %.1f min, steal = %d\n', ...
%        100*m.CVR, m.TMaxModel/60, m.Steal);
%
% See also: runFitVasoreactivity, getPulsatilityMetrics, getVasomotionMetrics,
%           blockDecimate, quasiUniform, lsqcurvefit
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

%------------- BEGIN CODE --------------
function out = fitVasoreactivity(arg1,arg2,arg3)
if nargin==2
    out=vrcSetup(arg1,arg2);            % SETUP:    (time, s)             -> layout
elseif nargin==3
    out=vrcAnalyze(arg1,arg2,arg3);     % ANALYSIS: (series, layout, s)   -> m
else
    error('fitVasoreactivity:nargin', ...
        ['Use layout=fitVasoreactivity(time,s) for SETUP or ' ...
         'm=fitVasoreactivity(series,layout,s) for ANALYSIS.']);
end
end

% =====================================================================
function L=vrcSetup(time,s)
%vrcSetup  Everything that depends only on the recording clock and the settings.
%   The per-trace call is then a block mean and a bounded least squares, which is what
%   makes even a 45-minute recording affordable.
tRaw=double(time(:));
L.nTraw=numel(tRaw);
if L.nTraw<8
    error('fitVasoreactivity:shortRecording', ...
        'A recording of %d samples is too short to fit a drug response.',L.nTraw);
end
L.fsRaw=1/median(diff(tRaw));
L.t1   =tRaw(1);

% ---- the model and the levels ----------------------------------------------
L.model=resolveModel(s);
L.want =vrcParseReturn(s);

% ---- the decimation --------------------------------------------------------
% An odd window and a stride of the same length: a non-overlapping BLOCK MEAN, which
% low-passes before it sub-samples.  A plain stride would alias the cardiac and
% respiratory content straight into the band the drug response lives in.
L.avgN=blockDecimate('window',L.fsRaw,targetFS(s));
L.time=blockDecimate(tRaw,L.avgN);
L.nT=numel(L.time);
if L.nT<8
    error('fitVasoreactivity:overDecimated', ...
        ['Decimating to %g Hz leaves only %d samples; lower s.tgtFS or leave it ' ...
         'empty to keep the recording''s own rate.'],targetFS(s),L.nT);
end
L.dt=median(diff(L.time));
L.fs=1/L.dt;

% ---- the injection anchor and the baseline window ---------------------------
% Both are read on the PRODUCT's clock, anchored to time(1) rather than assuming the
% base starts at 0: it costs one addition and cannot be silently wrong.  For every
% product this step accepts, time(1) is the recording start.
inj=requireScalar(s,'injectionSec','when the drug was given');
L.tInj=L.t1+inj;
L.inj =inj;
if L.tInj<=L.time(1) || L.tInj>=L.time(end)
    error('fitVasoreactivity:injectionOutsideRecording', ...
        ['The injection at %g s falls outside the recording, which runs from %g to ' ...
         '%g s.'],inj,L.time(1)-L.t1,L.time(end)-L.t1);
end

if ~isfield(s,'baselineSec') || numel(s.baselineSec)~=2
    error('fitVasoreactivity:noBaseline', ...
        ['s.baselineSec must be the two-element pre-injection baseline window ' ...
         '[t1 t2], in seconds from the recording start.']);
end
bl=double(s.baselineSec(:))';
if max(bl)>inj
    error('fitVasoreactivity:baselineAfterInjection', ...
        ['The baseline window [%g %g] s does not end before the injection at %g s, ' ...
         'so the "baseline" already contains the response.'],min(bl),max(bl),inj);
end
L.bl   =[min(bl) max(bl)];
L.blIdx=L.time>=L.t1+L.bl(1) & L.time<=L.t1+L.bl(2);
if sum(L.blIdx)<2
    error('fitVasoreactivity:emptyBaseline', ...
        ['Fewer than two samples fall inside the baseline window [%g %g] s at %g Hz; ' ...
         'widen the window or raise s.tgtFS.'],L.bl(1),L.bl(2),L.fs);
end
% tRef is the baseline-window CENTRE, so the fitted b0 is the baseline LEVEL rather
% than the level at an arbitrary recording edge extrapolated back along the drift.
L.tRef=L.t1+mean(L.bl);

% The response window - everything from the injection onward - is where the markers
% look for the peak.  Stored as indices so the per-trace call does no find().
L.respIdx=find(L.time>=L.tInj);
if numel(L.respIdx)<3
    error('fitVasoreactivity:shortResponse', ...
        ['Only %d samples follow the injection at %g s in a %g s recording; there ' ...
         'is no response to measure.'],numel(L.respIdx),inj,L.time(end)-L.t1);
end

% ---- marker settings -------------------------------------------------------
L.aucSec   =optScalar(s,'vrcAucSec',   2700);  % 45 min, the ACZ washout horizon
L.stealK   =optScalar(s,'vrcStealK',   2);
L.stealFrac=optScalar(s,'vrcStealFrac',0.15);
L.stealSec =optScalar(s,'vrcStealSec', 60);
L.washFrac =optScalar(s,'vrcWashFrac', 0.1);
if ~(L.washFrac>0 && L.washFrac<1)
    error('fitVasoreactivity:vrcWashFrac', ...
        's.vrcWashFrac must be a fraction of the peak strictly between 0 and 1, not %g.', ...
        L.washFrac);
end

% ---- parameter layout, bounds and the multi-start grid ---------------------
L=vrcParameterSetup(L,s);

% ---- the scalar contract ---------------------------------------------------
L.scalarNames=vrcScalarNames(L.want);

% ---- the optimiser ---------------------------------------------------------
L.opts=optimoptions('lsqcurvefit','Display','off', ...
    'Algorithm','trust-region-reflective', ...
    'MaxIterations',200,'MaxFunctionEvaluations',4000, ...
    'FunctionTolerance',1e-10,'StepTolerance',1e-12,'OptimalityTolerance',1e-8);
end

% =====================================================================
function L=vrcParameterSetup(L,s)
%vrcParameterSetup  The parameter vector, its bounds and the start grid.
%   Parameter order:  [b0 b1 A t0 r1 dr],  r2 = r1 + dr  (see the header on why the
%   GAP is the parameter).  Index constants are stored so nothing downstream counts
%   positions by hand.
L.iB0=1; L.iB1=2; L.iA=3; L.iT0=4; L.iR1=5; L.iDR=6;
L.nP =6;

% t0 is the delay from the injection to the response onset.  Negative is allowed: an
% injection mark a minute early is an ordinary protocol error, and clipping it at zero
% would push the error into the rate constants instead.  The range is generous because
% a systemically delivered drug can take minutes to reach the effect site.
t0Lo=-60;  t0Hi=600;
% RATE BOUNDS ARE WIDE, START BOUNDS ARE PHYSIOLOGICAL.  r1 in [1e-5 1] 1/s covers
% elimination from a second to most of a day; dr keeps r2 strictly above r1.
r1Lo=1e-5; r1Hi=1;
drLo=1e-6; drHi=5;
% The start box: elimination from ~50 s to ~9 h, absorption from ~5 s to ~3 h.
startLo=[3e-5 1e-4];  startHi=[2e-2 2e-1];      % [r1 dr], drawn log-uniformly

L.lb=[-Inf -Inf -Inf t0Lo r1Lo drLo];
L.ub=[ Inf  Inf  Inf t0Hi r1Hi drHi];

% ---- the start grid --------------------------------------------------------
% THE SHAPE COLUMNS ONLY.  b0, b1 and A are data-driven - they are the baseline, the
% drift and the largest excursion of the trace in front of us - and cannot be known at
% setup, so they are left NaN here and filled per trace in ANALYSIS.  A NaN that
% reached the optimiser would be an obvious fault rather than a quiet one.
nStarts=16;
if isfield(s,'vrcStarts') && ~isempty(s.vrcStarts)
    nStarts=double(s.vrcStarts);
    if ~(isscalar(nStarts)&&isfinite(nStarts)&&nStarts>=1)
        error('fitVasoreactivity:vrcStarts','s.vrcStarts must be a positive integer scalar.');
    end
    nStarts=max(1,round(nStarts));
end
L.nStarts=nStarts;

% Row 1 is the acetazolamide start point of Fahlstrom et al. 2024, CONVERTED TO 1/s
% (their t is in minutes: r1 = 0.007/min, r2 = 0.33/min).  Row 2 is a faster gas
% challenge - 10 min elimination, 1 min absorption, t_max about 2.6 min - so the grid
% is seeded at both ends of the protocols this step is meant to serve.
seed=[0  0.007/60  (0.33-0.007)/60;             % [t0 r1 dr]
      30 1/600     1/60-1/600     ];
nSeed=size(seed,1);
S=nan(nStarts,L.nP);
S(1:min(nSeed,nStarts),L.iT0:L.nP)=seed(1:min(nSeed,nStarts),:);

if nStarts>nSeed
    nQ=nStarts-nSeed;
    U=quasiUniform(nQ,3);                       % t0 and the two rate columns
    Q=zeros(nQ,3);
    Q(:,1)=0+U(:,1)*300;                        % t0, uniform over the plausible delays
    for k=1:2                                   % r1 and dr, log-uniform
        lo=startLo(k); hi=startHi(k);
        Q(:,k+1)=exp(log(lo)+U(:,k+1)*(log(hi)-log(lo)));
    end
    % the start box is a sub-box of the bounds; clip anyway, so a future edit to
    % either cannot put a start outside its own bound
    Q=min(max(Q,[t0Lo r1Lo drLo]),[t0Hi r1Hi drHi]);
    S(nSeed+1:end,L.iT0:L.nP)=Q;
end
L.starts=S;
end

% =====================================================================
function m=vrcAnalyze(series,L,~)
%vrcAnalyze  The decimation, the markers and the fit for ONE trace.
%   (s is accepted for signature symmetry with the pulsatility, vasomotion and NVC
%   cores but unused: every per-trace parameter is baked into `layout` at SETUP.)
%
%   THE BAG IS BUILT NaN-COMPLETE FIRST.  Every name the layout promises is written as
%   NaN before anything is computed, so the "invalid trace -> NaN in every field" rule
%   is enforced in one place and cannot be broken by a branch that forgets a field, and
%   so the caller can index the bag by layout.scalarNames unconditionally.
series=double(series(:));
m=struct('valid',numel(series)==L.nTraw && all(isfinite(series)));
for k=1:numel(L.scalarNames)
    m.(L.scalarNames{k})=single(NaN);
end
if L.want.reconstruction
    m.fData=nan(L.nT,1,'single');
end
if ~m.valid, return, end

% the trace arrives at the recording's own rate and is decimated HERE, so the caller
% never has to decimate a cube of its own (see the header)
y=blockDecimate(series,L.avgN);

% ---- model-free markers, ON THE DECIMATED RAW TRACE (see the header) --------
if L.want.markers
    mk=vrcMarkers(y,L);
    m.Peak  =single(mk.peak);    m.PeakRel =single(mk.peakRel);
    m.TMax  =single(mk.tmax);    m.AUC     =single(mk.auc);
    m.TWash =single(mk.twash);
    m.BlMean=single(mk.blMean);  m.BlStd   =single(mk.blStd);
    m.Steal =single(mk.steal);   m.StealSec=single(mk.stealSec);
end

if ~(L.want.model || L.want.reconstruction), return, end

% ---- the fit ---------------------------------------------------------------
P0=fillDataStarts(L.starts,y,L);
[p,q]=multiStartFit(y,L,P0);
fData=vrcEval(p,L);

if L.want.reconstruction
    m.fData=single(fData);
end
if ~L.want.model, return, end

b0=p(L.iB0); A=p(L.iA); t0=p(L.iT0);
r1=p(L.iR1); r2=r1+p(L.iDR);
m.Gain =single(A);  m.Onset=single(t0);
m.Drift=single(p(L.iB1)); m.Baseline=single(b0);
m.RateFast=single(r2); m.RateSlow=single(r1);

% tp is the bracket's own peak time - the same closed form the normalisation uses, and
% the same degenerate branch, so the reported peak time and the fitted curve cannot
% disagree about where the peak is
tp=peakTime(r1,r2);
m.TMaxModel=single(t0+tp);
if tp>0,   m.K  =single(A/tp); end
if b0~=0,  m.CVR=single(A/b0); end

m.R2=single(q.r2); m.RMSE=single(q.rmse); m.AICc=single(q.aicc);
m.Rho=single(q.rho); m.NEff=single(q.nEff);
m.GainSE=single(q.gainSE);
if b0~=0, m.CVRSE=single(q.gainSE/abs(b0)); end
m.NIter=single(q.nIter); m.Converged=single(q.converged);
m.StartsAgree=single(q.nAgree);

% the same marker routine on the fitted curve, so a model-free marker and its Fit twin
% are one measurement of two different curves and any difference between them is a
% statement about the fit rather than about two definitions
mf=vrcMarkers(fData,L);
m.PeakFit=single(mf.peak); m.AUCFit=single(mf.auc); m.TWashFit=single(mf.twash);
end

% =====================================================================
function P0=fillDataStarts(P0,y,L)
%fillDataStarts  The data-driven start columns: baseline, drift and amplitude.
%   These are read off the trace rather than guessed, which is why they are not in the
%   setup grid.  The amplitude is the largest excursion from the baseline in the
%   response window, which carries the right SIGN as well as the right order of
%   magnitude - a constrictor and a dilator start from their own data.
bl=mean(y(L.blIdx),'omitnan');
r =y-bl;
[~,ip]=max(abs(r(L.respIdx)));
A0=r(L.respIdx(ip));
if A0==0, A0=std(y); end
if A0==0, A0=1; end
half=floor(L.nT/2);
b1=(mean(y(half+1:end))-mean(y(1:half)))/max(L.dt*half,eps);

P0(:,L.iB0)=bl;
P0(:,L.iB1)=b1;
P0(:,L.iA )=A0;
end

% =====================================================================
function [pBest,q]=multiStartFit(y,L,P0)
%multiStartFit  Bounded least squares from many start points, best SSE wins.
%   THE MULTI-START IS A PLAIN LOOP, not Global Optimization Toolbox: the library
%   depends on Optimization Toolbox and nothing more, and a loop over deterministic
%   start points is reproducible in a way MultiStart's random sets are not.
%
%   A start point that throws is skipped rather than fatal - the optimiser can hit a
%   parameter combination the model cannot evaluate, and losing one of sixteen starts
%   is not a reason to lose the trace.
nStart=size(P0,1);
fun=@(p,~) vrcEval(p,L);

pBest=nan(1,L.nP); best=Inf; bestOut=[]; bestFlag=-99; bestJ=[];
rn=inf(nStart,1);
for k=1:nStart
    try
        [pk,resnorm,~,exitflag,out,~,Jk]= ...
            lsqcurvefit(fun,P0(k,:),L.time,y,L.lb,L.ub,L.opts);
    catch
        continue
    end
    if ~isfinite(resnorm), continue, end
    rn(k)=resnorm;
    if resnorm<best
        best=resnorm; pBest=pk; bestOut=out; bestFlag=exitflag; bestJ=Jk;
    end
end

if ~isfinite(best)
    % every start failed - report the failure rather than a number nobody can trace
    q=struct('r2',NaN,'rmse',NaN,'aicc',NaN,'rho',NaN,'nEff',NaN,'gainSE',NaN, ...
        'nIter',NaN,'converged',0,'nAgree',0);
    return
end

% how many starts landed on the winner (see the header: a lone winner is a warning)
q.nAgree=sum(rn<=best*(1+1e-4)+1e-12);

yhat=vrcEval(pBest,L);
res=y-yhat;
sse=sum(res.^2);
sst=sum((y-mean(y)).^2);
n=L.nT;
kPar=L.nP+1;                          % +1 for the residual variance (Gaussian LS)
q.r2  =1-sse/sst;
q.rmse=sqrt(sse/n);
% THE SAME AICc CONVENTION AS fitNVC, deliberately, so rsAICc and nsAICc are
% comparable numbers rather than two dialects of one name.
if n-kPar-1>0
    q.aicc=n*log(sse/n)+2*kPar+2*kPar*(kPar+1)/(n-kPar-1);
else
    q.aicc=NaN;                       % too few samples for the correction to exist
end
q.nIter=bestOut.iterations;
q.converged=double(bestFlag>0);

% ---- the AR(1)-corrected uncertainty on the amplitude ----------------------
% WITH ONE CURVE PER ANIMAL THERE IS NO BOOTSTRAP, so the only interval available is
% the one the residuals imply - and a physiological residual is strongly
% autocorrelated, so the textbook standard error is wildly optimistic.  Correcting by
% the AR(1) variance inflation (1+rho)/(1-rho) and STORING rho and the effective
% sample size makes that optimism visible instead of implicit.
[q.rho,q.nEff]=residualAR1(res,n);
q.gainSE=NaN;
if ~isempty(bestJ)
    J=full(bestJ);
    if size(J,1)==n && size(J,2)==L.nP && n>kPar && all(isfinite(J(:)))
        % THE COLUMNS ARE EQUILIBRATED BEFORE THE INVERSE.  A baseline near 100, a
        % drift per second near 1e-3 and a rate near 1e-4 put ten orders of magnitude
        % between the columns of J, so J'*J is numerically singular for a reason that
        % has nothing to do with identifiability - and an rcond test on the unscaled
        % matrix would reject every good fit.  Scaling to unit column norm and undoing
        % it on the covariance is exact.
        sc=sqrt(sum(J.^2,1)); sc(sc==0|~isfinite(sc))=1;
        M=(J./sc)'*(J./sc);
        if rcond(M)>1e-12
            C=(sse/(n-L.nP))*((M\eye(L.nP))./(sc'*sc));
            v=C(L.iA,L.iA)*(1+q.rho)/(1-q.rho);
            if v>0, q.gainSE=sqrt(v); end
        end
    end
end
end

% =====================================================================
function [rho,nEff]=residualAR1(res,n)
%residualAR1  Lag-1 autocorrelation of the residual, and the sample size it leaves.
%   rho is CLIPPED to [0, 0.99]: a negative lag-1 correlation would DEFLATE the
%   interval, which is not what this correction is for, and rho -> 1 sends it to
%   infinity for a trace that merely drifts.
rho=0; nEff=n;
if n<3, return, end
r=res-mean(res);
d=sum(r.^2);
if d<=0, return, end
rho=sum(r(1:end-1).*r(2:end))/d;
if ~isfinite(rho), rho=0; end
rho=max(0,min(0.99,rho));
nEff=n*(1-rho)/(1+rho);
end

% =====================================================================
function y=vrcEval(p,L)
%vrcEval  The observation model  y = b0 + b1*(t-tRef) + A*g(t-tInj-t0),  g of unit peak.
x=L.time-L.tInj-p(L.iT0);
y=p(L.iB0)+p(L.iB1)*(L.time-L.tRef)+p(L.iA)*bateNorm(x,p(L.iR1),p(L.iR1)+p(L.iDR));
end

% =====================================================================
function g=bateNorm(x,r1,r2)
%bateNorm  The bi-exponential bracket NORMALISED TO UNIT PEAK, zero for x <= 0.
%   Dividing by the bracket's own peak is what makes A the peak change in the trace's
%   own units and decorrelates it from the two rates (see the header).
g=zeros(size(x));
pos=x>0;
if ~any(pos), return, end
t=x(pos);
dr=r2-r1;
% THE DEGENERATE POINT IS A BAND, NOT AN EQUALITY.  As r2 -> r1 the bracket and its
% peak both go to zero and their ratio tends to e*r1*t*exp(-r1*t); an == test never
% fires and a 0/0 does.  The band is RELATIVE to r1, the natural scale of dr, so it
% means the same thing for a 1 s and a 3 h time constant.
if dr>1e-6*max(r1,realmin)
    tp=log(r2/r1)/dr;
    pk=exp(-r1*tp)-exp(-r2*tp);
    g(pos)=(exp(-r1*t)-exp(-r2*t))/pk;
else
    g(pos)=exp(1)*r1*t.*exp(-r1*t);     % the limit, also of unit peak (at t = 1/r1)
end
end

% =====================================================================
function tp=peakTime(r1,r2)
%peakTime  Where the unit-peak bracket peaks, in seconds after the onset.
%   ONE definition, shared by bateNorm and by TMaxModel, so the reported peak time and
%   the fitted curve cannot disagree about where the peak is.
dr=r2-r1;
if dr>1e-6*max(r1,realmin)
    tp=log(r2/r1)/dr;
else
    tp=1/r1;
end
end

% =====================================================================
function mk=vrcMarkers(y,L)
%vrcMarkers  The model-free markers of ONE curve - the raw decimated trace, or the
%   fitted one.  ONE routine for both, so a marker and its Fit twin are the same
%   measurement of two curves.
mk=struct('peak',NaN,'peakRel',NaN,'tmax',NaN,'auc',NaN,'twash',NaN, ...
    'blMean',NaN,'blStd',NaN,'steal',NaN,'stealSec',NaN);

bl=y(L.blIdx);
mk.blMean=mean(bl,'omitnan');
mk.blStd =std(bl,0,'omitnan');

r =y-mk.blMean;                       % change from the baseline level
rr=r(L.respIdx);                      % ...over the response window
tw=L.time(L.respIdx);

% THE LARGEST EXCURSION IN EITHER DIRECTION, not the maximum.  A constrictor is a
% first-class case for a drug challenge (it is not one for a whisker stimulus, which
% is why fitNVC takes max), and a peak defined as the maximum would report the noise
% ceiling of a purely constricting trace.
[~,k]=max(abs(rr));
pk=rr(k);
if ~isfinite(pk), return, end
mk.peak=pk;
mk.tmax=tw(k)-L.tInj;                 % measured from the INJECTION
if mk.blMean~=0, mk.peakRel=pk/mk.blMean; end

% AUC over a FIXED horizon, so two recordings of different length are comparable
aIdx=tw<=L.tInj+L.aucSec;
if sum(aIdx)>=2, mk.auc=trapz(tw(aIdx),rr(aIdx)); end

% washout: the first sample after the peak back within a fraction of the peak.  NOT a
% closed form - the equation is transcendental - so it is read off the curve, which is
% also what makes it the same measurement on the raw trace and on the fit.  A curve
% that never comes back (a strong drift, or a steal) leaves it NaN, which is itself
% the answer.
if pk~=0
    j=find(abs(rr(k:end))<=L.washFrac*abs(pk),1,'first');
    if ~isempty(j), mk.twash=tw(k+j-1)-L.tInj; end
end

% ---- the steal diagnostic, computed for EVERY trace (see the header) -------
% THE BASELINE-WINDOW LINE IS REMOVED FIRST, and that is not a refinement - without it
% the flag is worthless.  Over half an hour an LSCI baseline drifts by an amount
% comparable with the response itself, so a trace measured against the baseline MEAN
% ends the recording well below it with no steal anywhere, and the flag fires on
% EVERY long recording.  The line through the baseline window is the only drift
% estimate available without a fit, which is what keeps this marker model-free.
[c0,c1,sdr]=baselineLine(y,L);
zz=y(L.respIdx)-(c0+c1*(L.time(L.respIdx)-L.tRef));
% THE THRESHOLD IS THE LARGER OF TWO SCALES, and each guards a different mistake.
% vrcStealK*sd catches a dip that is large against the baseline noise, which is what
% matters on a quiet recording with a small response; vrcStealFrac*|Peak| stops a
% clean, strongly-responding recording being flagged by the error in extrapolating a
% few minutes of drift estimate across half an hour.
below=zz<-max(L.stealK*sdr,L.stealFrac*abs(pk));
below(1:k)=false;                     % only AFTER the peak
mk.stealSec=longestRun(below)*L.dt;
mk.steal=double(mk.stealSec>=L.stealSec);
end

% =====================================================================
function [c0,c1,sd]=baselineLine(y,L)
%baselineLine  The least-squares line through the baseline window, about tRef, and the
%   spread about it.  sd is the residual spread, NOT std(baseline): a drifting baseline
%   would otherwise inflate its own noise estimate and raise the steal threshold for
%   exactly the recordings that need it lowest.
tb=L.time(L.blIdx)-L.tRef;
yb=y(L.blIdx);
X=[ones(numel(tb),1) tb];
c=X\yb;
c0=c(1); c1=c(2);
sd=std(yb-X*c,0);
if ~isfinite(sd) || sd<0, sd=0; end
end

% =====================================================================
function n=longestRun(tf)
%longestRun  The longest run of consecutive true values (0 when there is none).
%   Written as the gap between the rising and falling edges of the padded run rather
%   than as a loop, because the per-pixel path calls it once per pixel per curve.
tf=logical(tf(:));
n=0;
if ~any(tf), return, end
d=diff([0;tf;0]);
n=max(find(d==-1)-find(d==1));
end

% =====================================================================
function tgt=targetFS(s)
%targetFS  s.tgtFS with its default; EMPTY means "keep every sample", not "1 Hz".
tgt=1;
if isfield(s,'tgtFS'), tgt=double(s.tgtFS); end
end

% =====================================================================
function names=vrcScalarNames(want)
%vrcScalarNames  The ORDERED scalar contract for one layout.
%   ANALYSIS emits exactly these names and no others, so the wrapper can build its
%   result columns from this list and never has to know which model ran.  The list
%   lives here, beside the code that fills it, and nowhere else.  A second model would
%   add its own parameter names to the 'model' block; the markers do not move.
names={};
if want.markers
    names=[names,{'Peak','PeakRel','TMax','AUC','TWash','BlMean','BlStd', ...
        'Steal','StealSec'}];
end
if want.model
    names=[names,{'Gain','Onset','Baseline','Drift','RateFast','RateSlow', ...
        'TMaxModel','K','CVR'}];
    names=[names,{'R2','RMSE','AICc','Rho','NEff','GainSE','CVRSE', ...
        'NIter','Converged','StartsAgree'}];
    names=[names,{'PeakFit','AUCFit','TWashFit'}];
end
end

% =====================================================================
function want=vrcParseReturn(s)
%vrcParseReturn  Resolve s.segVrcReturn into per-level compute logicals.
%   Absent OR empty gives the documented default, the complete set.  The per-pixel
%   caller passes s.ppxVrcReturn through this field, exactly as the NVC core takes
%   s.ppxNvcReturn through s.segNvcReturn.
levels={'markers','model','reconstruction'};
if isfield(s,'segVrcReturn') && ~isempty(s.segVrcReturn)
    sel=s.segVrcReturn;
    if ischar(sel)||isstring(sel), sel=cellstr(sel); end
    if ~iscellstr(sel)
        error('fitVasoreactivity:segVrcReturn', ...
            's.segVrcReturn must be a cell array of level names.');
    end
    bad=~ismember(sel,levels);
    if any(bad)
        error('fitVasoreactivity:segVrcReturn', ...
            'Unknown s.segVrcReturn level(s): %s. Valid levels: %s.', ...
            strjoin(sel(bad),', '),strjoin(levels,', '));
    end
else
    sel=levels;
end
for i=1:numel(levels), want.(levels{i})=ismember(levels{i},sel); end
end

% =====================================================================
function model=resolveModel(s)
%resolveModel  Default/validate s.vrcModel.
%   'bateman' is the only model in v1.  The three that are not here are named in the
%   error rather than silently ignored, because someone WILL type 'fahlstrom' expecting
%   it to work, and a wrong answer is worse than a refusal.
model='bateman';
if isfield(s,'vrcModel') && ~isempty(s.vrcModel)
    model=char(string(s.vrcModel));
end
if ~strcmp(model,'bateman')
    error('fitVasoreactivity:vrcModel', ...
        ['s.vrcModel must be ''bateman'', not ''%s''.  ''fahlstrom'', ''expConv'' ' ...
         'and ''secondOrder'' are designed but not built - see the header.'],model);
end
end

% =====================================================================
function v=requireScalar(s,name,what)
%requireScalar  A required numeric setting, named in the error rather than found
%   missing three lines later inside an arithmetic expression.
if ~isfield(s,name) || isempty(s.(name))
    error('fitVasoreactivity:missingSetting', ...
        's.%s is required - %s, in seconds.',name,what);
end
v=double(s.(name));
if ~(isscalar(v)&&isfinite(v))
    error('fitVasoreactivity:badSetting','s.%s must be a finite scalar (%s).',name,what);
end
end

% =====================================================================
function v=optScalar(s,name,dflt)
%optScalar  An optional numeric setting with a documented default.
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=double(s.(name)); end
if ~(isscalar(v)&&isfinite(v))
    error('fitVasoreactivity:badSetting','s.%s must be a finite scalar.',name);
end
end
