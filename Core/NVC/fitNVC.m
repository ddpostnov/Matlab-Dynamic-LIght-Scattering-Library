%fitNVC - Stimulus-locked neurovascular-coupling response fit (setup + per-trace fit)
%
%   fitNVC centralises the neurovascular-coupling "science" that runFitNVC applies
%   to every segment trace and every masked pixel of an epoch-averaged product.  It
%   builds everything that depends only on the epoch time base ONCE, and reuses it
%   for every trace on that base.  Two call modes are dispatched by argument count:
%
%   SETUP  (nargin 2)
%       layout = fitNVC(time, s)
%     Builds the stimulus boxcar, the baseline window, the parameter bounds and the
%     multi-start grid implied by `time` and the stimulus geometry in s, and
%     resolves s.segNvcReturn once.  Returns a small `layout` reused per trace.
%
%   ANALYSIS  (nargin 3)
%       m = fitNVC(series, layout, s)     % series = [nT x 1] one epoch-averaged trace
%     The science for ONE trace.  Returns a FLAT, prefix-free bag `m` (the caller
%     applies the ns/nd branch prefix, exactly as runPulsatility applies ps/pd);
%     every numeric field is SINGLE.
%
%   THE MODEL IS FITTED AS A KERNEL AND CONVOLVED WITH THE STIMULUS, NEVER FITTED
%   AS A SHAPE.  This is the single most important property of this core.  If the
%   observed shape were fitted directly, the time constants would absorb the
%   stimulus duration and stop being properties of the vasculature - two protocols
%   with 5 s and 20 s of whisker stimulation would report different "tissue" time
%   constants from the same animal.  Here the free parameters describe the impulse
%   response h; what is compared against the data is (h * u), the response of that
%   kernel to THIS protocol's stimulus.  testFitNVC asserts the invariance directly.
%
%   THE FLOW MODEL  (s.nvcModel = 'secondOrder', the default)
%     In the Buxton/Friston hemodynamic model the flow-inducing signal obeys
%           fdot = s
%           sdot = eps*u(t) - s/tauS - (f-1)/tauF
%     where tauF is the AUTOREGULATORY FEEDBACK time constant.  Friston et al. (2000)
%     justify that feedback term by citing post-stimulus undershoots in rCBF, and
%     describe the resulting flow as a very dampened oscillator; recent bottom-up
%     arteriolar models reach the same place, the undershoot being delayed myogenic
%     feedback over-correcting after the dilation.  Eliminating s gives a second-order
%     linear system in x = f-1,
%           xddot + (1/tauS)*xdot + (1/tauF)*x = eps*u(t)
%     which is the only low-parameter model of the FLOW response with a mechanistic
%     undershoot.  The Balloon and Windkessel models are NOT alternatives: they
%     consume CBF as an input and cannot be fitted to a flow trace.
%
%     Writing sigma = 1/(2*tauS) and wd = sqrt(1/tauF - sigma^2), the impulse
%     response has the three damping branches
%           underdamped   tauF <  4*tauS^2   h = (1/wd)*exp(-sigma*t)*sin(wd*t)
%           critical      tauF == 4*tauS^2   h = t*exp(-sigma*t)
%           overdamped    tauF >  4*tauS^2   h = (1/(2*wr))*(exp(-(sigma-wr)*t)
%                                                          -exp(-(sigma+wr)*t))
%     ALL THREE ARE IMPLEMENTED, selected by a RELATIVE band around the critical
%     point rather than an exact == test, because the optimiser walks across that
%     boundary constantly and a branch that returns NaN there poisons the fit.
%
%   THE KERNEL IS NORMALISED TO UNIT AREA, which is the one place this departs from
%   the formulae above, and it changes no shape.  The impulse responses written
%   above have area tauF (second order) rather than 1, so A and tauF would be
%   multiplicatively entangled and the optimiser would trade one against the other
%   along a valley.  Dividing by tauF is exactly the "fold eps into A" the model
%   requires anyway, and it makes A readable: A is the plateau a very long stimulus
%   would reach, in the trace's own units.  Every derived quantity (zeta, omega0,
%   the ring period, TTP, the undershoot ratio) is unchanged by it.
%
%   OBSERVATION MODEL
%       u(t) = 1 for t in [tStim, tStim+D], else 0        % D = s.stimDurationSec
%       y(t) = b0 + b1*(t-tRef) + A*(h*u)(t-t0)
%     SIX free parameters: b0 b1 A t0 tauS tauF.  eps is not separately identifiable
%     from A and is folded into it - a seventh parameter here makes the fit wander.
%     THE DRIFT TERM b1 IS NOT OPTIONAL: without it drift is absorbed into the
%     undershoot, and differences in baseline drift are then reported as genotype
%     differences in undershoot depth.
%
%     Because u is a BOXCAR, (h*u)(t) = H(t-tStim) - H(t-tStim-D) where H is the
%     kernel's step response, and every kernel here has a closed-form H (the
%     second-order step response, the regularised incomplete gamma, 1-exp).  The
%     convolution is therefore exact rather than a discrete sum - no grid, no
%     quadrature error, and O(nT) per model evaluation.  A non-boxcar stimulus
%     would need a numerical conv() here instead; nothing else would change.
%
%   OPTIONAL PRE-DIP  (s.nvcDip, default false)
%     Some LSCI traces show a small decrease before the rise.  This is NOT the
%     classic "initial dip" (an oxygenation phenomenon); it is state-dependent,
%     contested, and in LSCI is readily produced by CBV contamination of the speckle
%     contrast, by stimulus-locked motion, or by ROI contamination from the surround
%     negative response.  So it is NOT in the default model.  With s.nvcDip true a
%     second, faster, first-order CONSTRICTOR pathway is added - the dilator/
%     constrictor structure the interneuron / 20-HETE literature supports:
%           hc(t) = (1/tauC)*exp(-t/tauC)
%           y(t)  = b0 + b1*(t-tRef) + A*(h*u)(t-t0) - Ac*(hc*u)(t-t0)
%     BOTH models are fitted in the same call and BOTH AICc values are stored, so
%     the comparison needs no second run.  NOTHING IS AUTO-SELECTED: the reported
%     Gain/TauS/TauF/Zeta are ALWAYS the base model's, so those columns mean the
%     same thing whether or not the dip was fitted, and the dip fit is reported by
%     DipGain / DipTau / DipR2 / AICcDip beside them.  The author decides.
%
%   ALTERNATIVE  (s.nvcModel = 'doubleGamma')
%     A descriptive fallback, fitted as an impulse response and convolved exactly as
%     above:
%       h(t) = t^(a1-1)*b1^a1*exp(-b1*t)/gamma(a1) - c*t^(a2-1)*b2^a2*exp(-b2*t)/gamma(a2)
%     Free: b0 b1drift A t0 a1 beta1 a2 beta2 c (nine).  DO NOT SEED IT WITH THE SPM
%     CANONICAL VALUES (alpha 6/16, beta 1/1, c 1/6) - those are human BOLD.  The
%     start points here come from RODENT FLOW: awake about 2 s latency / 0.6 s width,
%     urethane-anaesthetised about 4 s / 2.5 s (Martin et al. 2006,
%     doi:10.1016/j.neuroimage.2006.02.021).  Its A is not comparable with the
%     second-order model's: the difference-of-gammas has area 1-c, not 1.
%
%   OUT OF SCOPE FOR v1, DELIBERATELY: the Lindquist-Wager inverse-logit model (7
%   parameters).  It is the better instrument for GROUP comparisons of amplitude /
%   latency / width, but it is too many parameters for a per-pixel fit.  It belongs
%   in a later session if the author wants it - do not add it here by accident.
%
%   ANALYSIS LEVELS  (s.segNvcReturn; a cell subset of the three names below).
%   Default (absent or empty) is the COMPLETE set.  The per-pixel path reuses SETUP
%   with s.segNvcReturn := s.ppxNvcReturn.
%       'markers'         model-free scalars read straight off the trace.  Always
%                         computable, no fit:
%                           Peak       max change from the baseline mean, response window
%                           PeakRel    the same as a fraction of the baseline mean
%                           TTP        time of that peak, from the stimulus start
%                           FWHM       width of the response at half the peak change
%                           Undershoot most negative change AFTER the peak
%                           TTU        time of that minimum, from the stimulus start
%                           URatio     -Undershoot/Peak (positive for a real undershoot)
%                           AUC        integral of the change over the response window
%                           BlMean BlStd  baseline-window mean and standard deviation
%       'model'           fitted parameters, derived physiology and fit quality
%                         (runs the fit) - see the scalar list below
%       'reconstruction'  fData [nT x 1], the fitted curve (runs the fit)
%
%   THE MARKERS ARE COMPUTED ON THE RAW TRACE, NOT ON THE RECONSTRUCTION.  This
%   differs deliberately from getPulsatilityMetrics, which computes its markers on
%   fData to reproduce a pre-refactor behaviour.  Here the whole point of the
%   model-free markers is that they are an INDEPENDENT check on the fit: read off
%   the fit, they would check nothing.  The same five markers are ALSO re-derived
%   from the fitted curve and stored with a Fit suffix, so the two are directly
%   comparable in one metrics table.
%
%   SCALARS EMITTED AT THE 'model' LEVEL
%     Gain Onset Drift Baseline            A, t0, b1, b0
%     TauS TauF Zeta Omega0 RingPeriod Damping        ('secondOrder' only)
%     A1 Beta1 A2 Beta2 CRatio                        ('doubleGamma' only)
%     R2 RMSE AICc NIter Converged StartsAgree        fit quality
%     PeakFit TTPFit FWHMFit UndershootFit TTUFit     markers re-read off the fit
%     DipGain DipTau DipR2 AICcDip                    (s.nvcDip only)
%   EVERY SCALAR NAME STARTS WITH A CAPITAL, so the caller's prefixed name is a plain
%   concatenation (ns + Peak = nsPeak) and there is exactly one spelling of each
%   quantity between this bag, the saved tree and the metrics tables.
%   with  Zeta = sqrt(tauF)/(2*tauS),  Omega0 = 1/sqrt(tauF),  RingPeriod = 2*pi/wd
%   (NaN when there is no ringing to have a period), and Damping = sigma = 1/(2*tauS),
%   the exponential decay rate in 1/s.
%
%   ZETA IS THE HEADLINE NUMBER, and it is why this model was chosen over a
%   descriptive one.  The undershoot-to-peak ratio equals exp(-pi*zeta/sqrt(1-zeta^2))
%   in both the impulse and the step limit, which makes zeta far more robust to
%   differences in stimulus duration than raw undershoot depth: two groups stimulated
%   for different lengths of time can be compared on zeta and cannot be compared on
%   Undershoot.
%
%   AN INVALID TRACE IS NaN IN EVERY FIELD.  m.valid = all(isfinite(series)); an
%   invalid trace skips the fit and every emitted scalar, and fData, is NaN.  No
%   field falls back to time(1) the way the pulsatility core does for bit-identity
%   reasons - there is no legacy golden here, so the clean rule is used.
%
%   THE FIT IS MULTI-START, AND THAT IS NOT OPTIONAL EITHER.  The objective is
%   non-convex with many local minima (Lindquist & Wager needed simulated annealing
%   for a comparable model).  A single start returns quiet garbage on a minority of
%   traces and nothing in the output says so.  StartsAgree records how many start
%   points reached the winning solution: a trace where 1 start in 20 wins is a trace
%   to distrust.  The start grid is a deterministic low-discrepancy sequence, so the
%   same trace always gives the same answer.
%
%   TWO OF THE EMITTED SCALARS DESCRIBE THE OPTIMISER, NOT THE VESSEL.  NIter and
%   StartsAgree count what the trust region did on the way to the answer, and BLAS
%   threading differs between a parfor worker and the client, so the same optimum is
%   reached in a different number of steps and the tie test between start points can
%   break differently.  Every parameter, marker and quality number is reproducible;
%   those two are diagnostics, and a test that compares them across a parallel split
%   is asserting something that is not true.
%
% Syntax:
%    layout = fitNVC(time, s)              % SETUP    - once per time base
%    m      = fitNVC(series, layout, s)    % ANALYSIS - one trace
%
% Inputs:
%    time    - [nT x 1] epoch time base, seconds from the epoch start (SETUP).
%    series  - [nT x 1] one epoch-averaged trace on that base (ANALYSIS).
%    layout  - the struct SETUP returned.
%    s       - parameter struct.  Fields:
%                epochStimStartSec  when the stimulus starts within the epoch, s
%                stimDurationSec    how long the stimulus lasts, s (REQUIRED)
%                epochBaselineSec   [t1 t2] pre-stimulus baseline window, s
%                nvcModel           'secondOrder' (default) | 'doubleGamma'
%                nvcDip             false (default) | true
%                segNvcReturn       cell subset of {'markers','model',
%                                   'reconstruction'} (default: all three); the
%                                   per-pixel caller passes s.ppxNvcReturn through
%                                   this field
%                nvcStarts          number of multi-start points (default 16)
%                nvcWeights         [nT x 1] per-timepoint weights for a weighted
%                                   least squares (optional; e.g. 1./SEM.^2 across
%                                   the accepted epochs).  Absent = unweighted.
%
% Outputs:
%    layout  - (SETUP)    stimulus boxcar, windows, bounds, start grid and .want for
%                         one time base, plus .scalarNames, the ORDERED list of
%                         scalars ANALYSIS emits (the caller needs no other
%                         knowledge of the field set).
%    m       - (ANALYSIS) flat, prefix-free, all-single bag for one trace: exactly
%                         layout.scalarNames, plus fData when 'reconstruction' was
%                         asked for, plus the logical m.valid.
%
% Example:
%    s.epochStimStartSec=10; s.stimDurationSec=5; s.epochBaselineSec=[0 10];
%    layout = fitNVC(results.time, s);
%    m      = fitNVC(results.sData(:,7), layout, s);
%    fprintf('zeta = %.2f, undershoot/peak = %.2f\n', m.Zeta, m.URatio);
%
% See also: runFitNVC, getPulsatilityMetrics, getVasomotionMetrics, runExternalCycle,
%           lsqcurvefit
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

%------------- BEGIN CODE --------------
function out = fitNVC(arg1,arg2,arg3)
if nargin==2
    out=nvcSetup(arg1,arg2);            % SETUP:    (time, s)             -> layout
elseif nargin==3
    out=nvcAnalyze(arg1,arg2,arg3);     % ANALYSIS: (series, layout, s)   -> m
else
    error('fitNVC:nargin', ...
        ['Use layout=fitNVC(time,s) for SETUP or ' ...
         'm=fitNVC(series,layout,s) for ANALYSIS.']);
end
end

% =====================================================================
function L=nvcSetup(time,s)
%nvcSetup  Everything that depends only on the epoch time base and the settings.
%   The per-trace call is then a bounded least-squares and nothing else, which is
%   what makes the per-pixel path affordable at all.
time=double(time(:));
L.time=time;
L.nT  =numel(time);
if L.nT<4
    error('fitNVC:shortEpoch','An epoch of %d samples is too short to fit.',L.nT);
end
L.dt=median(diff(time));

% ---- the model, the dip and the levels -------------------------------------
L.model=resolveModel(s);
L.dip  =resolveFlag(s,'nvcDip',false);
L.want =nvcParseReturn(s);

% ---- stimulus geometry -----------------------------------------------------
% Read on the epoch clock: epochStimStartSec and epochBaselineSec are measured from
% the START OF THE EPOCH, and the epoch time base starts at time(1) (0 for every
% product runExternalCycle writes).  Anchoring to time(1) rather than assuming 0
% costs one addition and cannot be silently wrong.
tStim=requireScalar(s,'epochStimStartSec','when the stimulus starts within the epoch');
D    =requireScalar(s,'stimDurationSec',  'how long the stimulus lasts');
if ~(D>0)
    error('fitNVC:stimDuration','s.stimDurationSec must be positive, not %g.',D);
end
L.tStim=time(1)+tStim;
L.D    =D;
L.u    =double(time>=L.tStim & time<=L.tStim+D);     % the boxcar, [nT x 1]

if ~isfield(s,'epochBaselineSec') || numel(s.epochBaselineSec)~=2
    error('fitNVC:noBaseline', ...
        ['s.epochBaselineSec must be the two-element pre-stimulus baseline window ' ...
         '[t1 t2], in seconds from the epoch start.']);
end
bl=double(s.epochBaselineSec(:))';
L.blIdx=time>=time(1)+min(bl) & time<=time(1)+max(bl);
if ~any(L.blIdx)
    error('fitNVC:emptyBaseline', ...
        'No sample of the epoch falls inside the baseline window [%g %g] s.',bl(1),bl(2));
end
% tRef is the baseline-window CENTRE, so the fitted b0 is the baseline level rather
% than the level at an arbitrary epoch edge extrapolated back along the drift.
L.tRef=time(1)+mean(bl);

% The response window - everything from the stimulus onward - is where the markers
% look for the peak.  Stored as indices so the per-trace call does no find().
L.respIdx=find(time>=L.tStim);
if numel(L.respIdx)<3
    error('fitNVC:shortResponse', ...
        ['Only %d samples follow the stimulus start (%g s into a %g s epoch); ' ...
         'there is no response to measure.'],numel(L.respIdx),tStim,time(end)-time(1));
end

% ---- weights ---------------------------------------------------------------
% Weighted least squares is done by scaling BOTH the model and the data by sqrt(w),
% which is the standard reduction of a weighted problem to an unweighted one.
L.sw=ones(L.nT,1);
if isfield(s,'nvcWeights') && ~isempty(s.nvcWeights)
    w=double(s.nvcWeights(:));
    if numel(w)~=L.nT
        error('fitNVC:weights','s.nvcWeights must have %d elements, one per sample.',L.nT);
    end
    if any(~isfinite(w)) || any(w<0)
        error('fitNVC:weights','s.nvcWeights must be finite and non-negative.');
    end
    L.sw=sqrt(w);
end
L.weighted=any(L.sw~=1);

% ---- parameter layout, bounds and the multi-start grid ---------------------
L=nvcParameterSetup(L,s);

% ---- the scalar contract ---------------------------------------------------
L.scalarNames=nvcScalarNames(L.model,L.dip,L.want);

% ---- the optimiser ---------------------------------------------------------
% Trust-region-reflective is lsqcurvefit's bounded algorithm and the default; the
% tolerances are loose enough that a converged trace costs tens of iterations, not
% hundreds, because this runs once per start point per trace.
L.opts=optimoptions('lsqcurvefit','Display','off', ...
    'Algorithm','trust-region-reflective', ...
    'MaxIterations',200,'MaxFunctionEvaluations',4000, ...
    'FunctionTolerance',1e-10,'StepTolerance',1e-10,'OptimalityTolerance',1e-8);
end

% =====================================================================
function L=nvcParameterSetup(L,s)
%nvcParameterSetup  The parameter vector, its bounds and the start grid.
%   Parameter order, for both kernels:  [b0 b1 A t0 <shape...> (Ac tauC)]
%   The first four are common; the shape block is what the kernel adds; the dip pair
%   is appended when it is fitted.  Index constants are stored so nothing downstream
%   counts positions by hand.
L.iB0=1; L.iB1=2; L.iA=3; L.iT0=4;

% t0 is the response onset relative to the stimulus.  Negative is allowed: a trace
% whose baseline window is slightly mis-specified, or a genuinely fast region, can
% start before the nominal stimulus mark, and clipping that at zero would push the
% error into tauS instead.
t0Lo=-2; t0Hi=10;

switch L.model
    case 'secondOrder'
        % tauS is the signal-decay constant, tauF the autoregulatory feedback
        % constant.  The bounds are wide enough to hold every damping regime and
        % tight enough that the optimiser cannot wander into a flat corner.
        shapeLo =[0.05 0.02];   shapeHi =[20 50];      % tauS, tauF
        startLo =[0.20 0.50];   startHi =[8  25];      % where the starts are drawn
        shapeLog=[true true];                          % drawn log-uniformly
        L.shapeNames={'tauS','tauF'};
        seed=[0.5 2 4];                                % [t0 tauS tauF], zeta = 0.5
    case 'doubleGamma'
        % a>1 keeps a mode in the gamma (a<1 is singular at the origin); c is the
        % undershoot FRACTION, so c>=1 would invert the net response, which A
        % already describes.
        shapeLo =[1.01 0.05 1.01 0.05 0];
        shapeHi =[40   20   40   20   1];
        startLo =[2    0.30 2    0.15 0];
        startHi =[20   8    20   4    0.6];
        shapeLog=[true true true true false];
        L.shapeNames={'a1','beta1','a2','beta2','c'};
        % RODENT FLOW seeds, not the SPM canonical human-BOLD ones (see the header).
        % Gamma with shape a and rate b peaks at (a-1)/b and has width sqrt(a)/b, so
        % a latency/width pair inverts to a = m^2, b = m/width, m = (r+sqrt(r^2+4))/2,
        % r = latency/width.
        seed=[0    4.33 0.832 5.83 0.483 0.20;        % urethane: 4 s / 2.5 s
              0   13.03 6.02 13.03 2.41  0.20];       % awake:    2 s / 0.6 s
end

nShape=numel(shapeLo);
L.iShape=L.iT0+(1:nShape);
L.nPBase=L.iT0+nShape;

lb=[-Inf -Inf -Inf t0Lo shapeLo];
ub=[ Inf  Inf  Inf t0Hi shapeHi];

% The dip pair is appended, never interleaved, so the base parameter positions are
% the same in both fits and one evaluator serves both.
L.iAc=L.nPBase+1; L.iTauC=L.nPBase+2;
lbDip=[lb -Inf 0.05];
ubDip=[ub  Inf 20  ];

L.lb=lb; L.ub=ub;
L.lbDip=lbDip; L.ubDip=ubDip;
L.nPDip=L.nPBase+2;

% ---- the start grid --------------------------------------------------------
% THE SHAPE COLUMNS ONLY.  b0, b1, A and Ac are data-driven - they are the
% baseline, the drift and the peak of the trace in front of us - and cannot be
% known at setup, so they are left NaN here and filled per trace in ANALYSIS.  A
% NaN that reached the optimiser would be an obvious fault rather than a quiet one.
nStarts=16;
if isfield(s,'nvcStarts') && ~isempty(s.nvcStarts)
    nStarts=double(s.nvcStarts);
    if ~(isscalar(nStarts)&&isfinite(nStarts)&&nStarts>=1)
        error('fitNVC:nvcStarts','s.nvcStarts must be a positive integer scalar.');
    end
    nStarts=max(1,round(nStarts));
end
L.nStarts=nStarts;

nSeed=size(seed,1);
S=nan(nStarts,L.nPBase);
S(1:min(nSeed,nStarts),L.iT0:L.nPBase)=seed(1:min(nSeed,nStarts),:);

% The remaining rows come from Roberts' R_d sequence - a deterministic additive
% recurrence with the generalised golden ratio.  It spreads points far more evenly
% over the box than rand does, and, being deterministic, the same trace always
% gives the same answer (a fit that changes between runs cannot be pinned by a test).
if nStarts>nSeed
    nQ=nStarts-nSeed;
    U=quasiUniform(nQ,1+nShape);
    span=[t0Lo t0Hi; [startLo(:) startHi(:)]];
    isLog=[false shapeLog];
    Q=zeros(nQ,1+nShape);
    for k=1:1+nShape
        lo=span(k,1); hi=span(k,2);
        if isLog(k), Q(:,k)=exp(log(lo)+U(:,k)*(log(hi)-log(lo)));
        else,        Q(:,k)=lo+U(:,k)*(hi-lo);
        end
    end
    % the start box is a physiological sub-box of the bounds; clip anyway, so a
    % future edit to either cannot put a start outside its own bound
    Q=min(max(Q,[t0Lo shapeLo]),[t0Hi shapeHi]);
    S(nSeed+1:end,L.iT0:L.nPBase)=Q;
end
L.starts=S;
end

% =====================================================================
function U=quasiUniform(n,dim)
%quasiUniform  Roberts' R_d low-discrepancy sequence on the unit cube.
%   phi is the real root of x^(dim+1) = x + 1 (the generalised golden ratio), found
%   by the fixed-point iteration x <- (1+x)^(1/(dim+1)), which converges in a few
%   steps for every dim used here.  No random number generator is involved, so no
%   rng state leaks between traces or between workers.
phi=2;
for k=1:32, phi=(1+phi)^(1/(dim+1)); end
alpha=phi.^-(1:dim);
U=mod(0.5+(1:n)'*alpha,1);
end

% =====================================================================
function m=nvcAnalyze(series,L,~)
%nvcAnalyze  The fit and the markers for ONE trace.
%   (s is accepted for signature symmetry with the pulsatility and vasomotion cores
%   but unused: every per-trace parameter is baked into `layout` at SETUP.)
%
%   THE BAG IS BUILT NaN-COMPLETE FIRST.  Every name the layout promises is written
%   as NaN before anything is computed, so the "invalid trace -> NaN in every field"
%   rule is enforced in one place and cannot be broken by a branch that forgets a
%   field, and so the caller can index the bag by layout.scalarNames unconditionally.
series=double(series(:));
m=struct('valid',all(isfinite(series)) && numel(series)==L.nT);
for k=1:numel(L.scalarNames)
    m.(L.scalarNames{k})=single(NaN);
end
if L.want.reconstruction
    m.fData=nan(L.nT,1,'single');
end
if ~m.valid, return, end

% ---- model-free markers, ON THE RAW TRACE (see the header) -----------------
if L.want.markers
    mk=nvcMarkers(series,L);
    m.Peak=single(mk.peak);         m.PeakRel  =single(mk.peakRel);
    m.TTP =single(mk.ttp);          m.FWHM     =single(mk.fwhm);
    m.Undershoot=single(mk.under);  m.TTU      =single(mk.ttu);
    m.URatio    =single(mk.uRatio); m.AUC      =single(mk.auc);
    m.BlMean    =single(mk.blMean); m.BlStd    =single(mk.blStd);
end

if ~(L.want.model || L.want.reconstruction), return, end

% ---- the fit ---------------------------------------------------------------
P0=fillDataStarts(L.starts,series,L);
[p,q]=multiStartFit(series,L,P0,L.lb,L.ub,false);
fData=nvcEval(p,L,false);

if L.want.reconstruction
    m.fData=single(fData);
end
if ~L.want.model, return, end

m.Gain=single(p(L.iA)); m.Onset=single(p(L.iT0));
m.Drift=single(p(L.iB1)); m.Baseline=single(p(L.iB0));

switch L.model
    case 'secondOrder'
        tauS=p(L.iShape(1)); tauF=p(L.iShape(2));
        sigma=1/(2*tauS); d=1/tauF-sigma^2;         % d = wd^2
        m.TauS=single(tauS); m.TauF=single(tauF);
        m.Zeta   =single(sqrt(tauF)/(2*tauS));      % damping RATIO (dimensionless)
        m.Omega0 =single(1/sqrt(tauF));             % natural frequency, rad/s
        m.Damping=single(sigma);                    % decay RATE, 1/s
        % a critically- or over-damped response does not ring, so it has no ring
        % period - NaN rather than Inf, which would look like a very slow ring
        if d>0, m.RingPeriod=single(2*pi/sqrt(d)); end
    case 'doubleGamma'
        m.A1=single(p(L.iShape(1))); m.Beta1 =single(p(L.iShape(2)));
        m.A2=single(p(L.iShape(3))); m.Beta2 =single(p(L.iShape(4)));
        m.CRatio=single(p(L.iShape(5)));
end

m.R2=single(q.r2); m.RMSE=single(q.rmse); m.AICc=single(q.aicc);
m.NIter=single(q.nIter); m.Converged=single(q.converged);
m.StartsAgree=single(q.nAgree);

% the same marker routine on the fitted curve, so the model-free and the modelled
% markers are the same measurement of two different curves and can be compared
mf=nvcMarkers(fData,L);
m.PeakFit=single(mf.peak); m.TTPFit=single(mf.ttp); m.FWHMFit=single(mf.fwhm);
m.UndershootFit=single(mf.under); m.TTUFit=single(mf.ttu);

% ---- the dip model, fitted BESIDE the base one and never instead of it -----
if L.dip
    P0d=[P0 nan(size(P0,1),2)];
    P0d=fillDataStarts(P0d,series,L);
    [pd,qd]=multiStartFit(series,L,P0d,L.lbDip,L.ubDip,true);
    m.DipGain=single(pd(L.iAc)); m.DipTau=single(pd(L.iTauC));
    m.DipR2  =single(qd.r2);     m.AICcDip=single(qd.aicc);
end
end

% =====================================================================
function P0=fillDataStarts(P0,series,L)
%fillDataStarts  The data-driven start columns: baseline, drift and amplitude.
%   These are read off the trace rather than guessed, which is why they are not in
%   the setup grid.  The baseline is the mean of the baseline window; the drift is
%   the slope across the two halves of the epoch; the amplitude is the largest
%   excursion from the baseline in the response window, which is the right sign and
%   the right order of magnitude for a positive or a negative responder alike.
bl=mean(series(L.blIdx));
r =series-bl;
[~,ip]=max(abs(r(L.respIdx)));
A0=r(L.respIdx(ip));
if A0==0, A0=std(series); end
if A0==0, A0=1; end
half=floor(L.nT/2);
b1=(mean(series(half+1:end))-mean(series(1:half)))/max(L.dt*half,eps);

P0(:,L.iB0)=bl;
P0(:,L.iB1)=b1;
P0(:,L.iA )=A0;
if size(P0,2)>=L.iTauC
    % the constrictor starts as a fifth of the dilator and faster than it, which is
    % what makes it a PRE-dip rather than a second copy of the main response
    P0(:,L.iAc  )=0.2*A0;
    P0(:,L.iTauC)=0.5;
end
end

% =====================================================================
function [pBest,q]=multiStartFit(y,L,P0,lb,ub,useDip)
%multiStartFit  Bounded least squares from many start points, best SSE wins.
%   THE MULTI-START IS A PLAIN LOOP, not Global Optimization Toolbox: the library
%   depends on Optimization Toolbox and nothing more, and a loop over deterministic
%   start points is reproducible in a way MultiStart's random sets are not.
%
%   A start point that throws is skipped rather than fatal - the optimiser can hit a
%   parameter combination the kernel cannot evaluate, and losing one of sixteen
%   starts is not a reason to lose the trace.
nStart=size(P0,1);
sw=L.sw;
yw=sw.*y;
fun=@(p,~) sw.*nvcEval(p,L,useDip);

pBest=nan(1,size(P0,2)); best=Inf; bestOut=[]; bestFlag=-99;
rn=inf(nStart,1);
for k=1:nStart
    try
        [pk,resnorm,~,exitflag,out]=lsqcurvefit(fun,P0(k,:),L.time,yw,lb,ub,L.opts);
    catch
        continue
    end
    if ~isfinite(resnorm), continue, end
    rn(k)=resnorm;
    if resnorm<best
        best=resnorm; pBest=pk; bestOut=out; bestFlag=exitflag;
    end
end

if ~isfinite(best)
    % every start failed - report the failure rather than a number nobody can trace
    q=struct('r2',NaN,'rmse',NaN,'aicc',NaN,'nIter',NaN,'converged',0,'nAgree',0);
    return
end

% how many starts landed on the winner (see the header: a lone winner is a warning)
q.nAgree=sum(rn<=best*(1+1e-4)+1e-12);

% R2, RMSE and AICc are computed on the UNWEIGHTED residual, which is the number a
% reader expects to see, and the same convention is used for the base and the dip
% model so their AICc values are comparable.  With no weights the two coincide.
yhat=nvcEval(pBest,L,useDip);
res=y-yhat;
sse=sum(res.^2);
sst=sum((y-mean(y)).^2);
n=L.nT;
kPar=numel(pBest)+1;                 % +1 for the residual variance (Gaussian LS)
q.r2  =1-sse/sst;
q.rmse=sqrt(sse/n);
if n-kPar-1>0
    q.aicc=n*log(sse/n)+2*kPar+2*kPar*(kPar+1)/(n-kPar-1);
else
    q.aicc=NaN;                      % too few samples for the correction to exist
end
q.nIter=bestOut.iterations;
q.converged=double(bestFlag>0);
end

% =====================================================================
function y=nvcEval(p,L,useDip)
%nvcEval  The observation model  y = b0 + b1*(t-tRef) + A*(h*u)(t-t0) [- Ac*(hc*u)].
%   The convolution with the boxcar is done in closed form through the kernel's STEP
%   response (see the header): (h*u)(x) = H(x-tStim) - H(x-tStim-D).
x=L.time-p(L.iT0);                                   % the clock the kernel sees
g=stepResponse(x-L.tStim,L,p)-stepResponse(x-L.tStim-L.D,L,p);
y=p(L.iB0)+p(L.iB1)*(L.time-L.tRef)+p(L.iA)*g;
if useDip
    tauC=p(L.iTauC);
    gc=expStep(x-L.tStim,tauC)-expStep(x-L.tStim-L.D,tauC);
    y=y-p(L.iAc)*gc;
end
end

% =====================================================================
function H=stepResponse(tau,L,p)
%stepResponse  The kernel's UNIT-AREA step response H(tau), zero for tau < 0.
%   H is the integral of the impulse response from 0 to tau; convolving a boxcar
%   with h is the difference of two of these, which is why no discrete convolution
%   appears anywhere in this file.
H=zeros(size(tau));
pos=tau>0;
if ~any(pos), return, end
t=tau(pos);

switch L.model
    case 'secondOrder'
        tauS=p(L.iShape(1)); tauF=p(L.iShape(2));
        sigma=1/(2*tauS);
        d=1/tauF-sigma^2;                            % d = wd^2, the damping discriminant
        % THE CRITICAL POINT IS A BAND, NOT AN EQUALITY.  d is a difference of two
        % nearly-equal numbers near critical damping, so an == test never fires and
        % a 0/0 does; the band is RELATIVE to sigma^2, the natural scale of d, so it
        % means the same thing for a 0.05 s and a 20 s time constant.
        band=1e-8*max(sigma^2,realmin);
        if d>band
            wd=sqrt(d);                              % underdamped: it rings
            H(pos)=1-exp(-sigma*t).*(cos(wd*t)+(sigma/wd)*sin(wd*t));
        elseif d<-band
            % overdamped.  Written as two decaying exponentials rather than
            % exp(-sigma*t)*cosh(wr*t): wr < sigma strictly, so both exponents below
            % are negative and bounded, whereas cosh(wr*t) alone overflows to Inf on
            % a long epoch and turns the whole trace into NaN.
            wr=sqrt(-d); kk=sigma/wr;
            H(pos)=1-0.5*((1+kk)*exp(-(sigma-wr)*t)+(1-kk)*exp(-(sigma+wr)*t));
        else
            H(pos)=1-exp(-sigma*t).*(1+sigma*t);     % critically damped
        end
    case 'doubleGamma'
        a1=p(L.iShape(1)); b1=p(L.iShape(2));
        a2=p(L.iShape(3)); b2=p(L.iShape(4));
        c =p(L.iShape(5));
        % the integral of a gamma pdf IS the lower regularised incomplete gamma, so
        % this branch is closed form too
        H(pos)=gammainc(b1*t,a1)-c*gammainc(b2*t,a2);
end
end

% =====================================================================
function H=expStep(tau,tauC)
%expStep  Step response of the first-order constrictor pathway hc = exp(-t/tauC)/tauC.
H=zeros(size(tau));
pos=tau>0;
H(pos)=1-exp(-tau(pos)/tauC);
end

% =====================================================================
function mk=nvcMarkers(y,L)
%nvcMarkers  The model-free markers of ONE curve - the raw trace, or the fitted one.
%   ONE routine for both, so a model-free marker and its Fit twin are the same
%   measurement of two different curves and any difference between them is a
%   statement about the fit rather than about the two definitions.
mk=struct('peak',NaN,'peakRel',NaN,'ttp',NaN,'fwhm',NaN,'under',NaN, ...
    'ttu',NaN,'uRatio',NaN,'auc',NaN,'blMean',NaN,'blStd',NaN);

bl=y(L.blIdx);
mk.blMean=mean(bl,'omitnan');
mk.blStd =std(bl,0,'omitnan');

r =y-mk.blMean;                       % change from the baseline level
rr=r(L.respIdx);                      % ...over the response window
tw=L.time(L.respIdx);

[pk,k]=max(rr);
% A CURVE THAT IS NOT FINITE HAS NO MARKERS.  The raw trace cannot reach here
% non-finite (m.valid gates it), but the FITTED curve can: every start point of a
% pathological trace may fail, and max() of an all-NaN vector still returns index 1,
% which would put a real-looking time into TTPFit.
if ~isfinite(pk), return, end
mk.peak=pk;
mk.ttp =tw(k)-L.tStim;                % measured from the STIMULUS, not the epoch start
if mk.blMean~=0, mk.peakRel=pk/mk.blMean; end
mk.auc=trapz(tw,rr);

% ---- FWHM, by linear interpolation of the half-maximum crossings -----------
% Only meaningful for a positive excursion: for a region that does not respond, or
% responds downward, "the width at half the peak" is not a width of anything, and a
% number invented there would be indistinguishable from a real one downstream.
if pk>0
    half=pk/2;
    iL=find(rr(1:k)<half,1,'last');
    iR=find(rr(k:end)<half,1,'first');
    tL=NaN; tR=NaN;
    if ~isempty(iL), tL=crossAt(tw(iL),rr(iL),tw(iL+1),rr(iL+1),half); end
    if ~isempty(iR)
        j=k+iR-1;
        tR=crossAt(tw(j-1),rr(j-1),tw(j),rr(j),half);
    end
    mk.fwhm=tR-tL;                    % NaN when the curve never comes back down
end

% ---- undershoot: the deepest point AFTER the peak --------------------------
[us,ku]=min(rr(k:end));
mk.under=us;
mk.ttu  =tw(k+ku-1)-L.tStim;
% signed so that a classic response (peak up, undershoot down) gives a POSITIVE
% ratio directly comparable with the model's exp(-pi*zeta/sqrt(1-zeta^2))
if pk~=0, mk.uRatio=-us/pk; end
end

% =====================================================================
function t=crossAt(t1,y1,t2,y2,level)
%crossAt  Where a straight line through (t1,y1)-(t2,y2) crosses `level`.
if y2==y1, t=t2; return, end
t=t1+(level-y1)*(t2-t1)/(y2-y1);
end

% =====================================================================
function names=nvcScalarNames(model,dip,want)
%nvcScalarNames  The ORDERED scalar contract for one layout.
%   ANALYSIS emits exactly these names and no others, so the wrapper can build its
%   result columns from this list and never has to know which model ran.  The list
%   lives here, beside the code that fills it, and nowhere else.
names={};
if want.markers
    names=[names,{'Peak','PeakRel','TTP','FWHM','Undershoot','TTU','URatio','AUC', ...
        'BlMean','BlStd'}];
end
if want.model
    names=[names,{'Gain','Onset','Drift','Baseline'}];
    switch model
        case 'secondOrder'
            names=[names,{'TauS','TauF','Zeta','Omega0','RingPeriod','Damping'}];
        case 'doubleGamma'
            names=[names,{'A1','Beta1','A2','Beta2','CRatio'}];
    end
    names=[names,{'R2','RMSE','AICc','NIter','Converged','StartsAgree'}];
    names=[names,{'PeakFit','TTPFit','FWHMFit','UndershootFit','TTUFit'}];
    if dip
        names=[names,{'DipGain','DipTau','DipR2','AICcDip'}];
    end
end
end

% =====================================================================
function want=nvcParseReturn(s)
%nvcParseReturn  Resolve s.segNvcReturn into per-level compute logicals.
%   Absent OR empty gives the documented default, the complete set.  The per-pixel
%   caller passes s.ppxNvcReturn through this field, exactly as the pulsatility core
%   takes s.ppxPulsReturn through s.segPulsReturn.
levels={'markers','model','reconstruction'};
if isfield(s,'segNvcReturn') && ~isempty(s.segNvcReturn)
    sel=s.segNvcReturn;
    if ischar(sel)||isstring(sel), sel=cellstr(sel); end
    if ~iscellstr(sel)
        error('fitNVC:segNvcReturn','s.segNvcReturn must be a cell array of level names.');
    end
    bad=~ismember(sel,levels);
    if any(bad)
        error('fitNVC:segNvcReturn', ...
            'Unknown s.segNvcReturn level(s): %s. Valid levels: %s.', ...
            strjoin(sel(bad),', '),strjoin(levels,', '));
    end
else
    sel=levels;
end
for i=1:numel(levels), want.(levels{i})=ismember(levels{i},sel); end
end

% =====================================================================
function model=resolveModel(s)
%resolveModel  Default/validate s.nvcModel.
model='secondOrder';
if isfield(s,'nvcModel') && ~isempty(s.nvcModel)
    model=char(string(s.nvcModel));
end
if ~any(strcmp(model,{'secondOrder','doubleGamma'}))
    error('fitNVC:nvcModel', ...
        's.nvcModel must be ''secondOrder'' or ''doubleGamma'', not ''%s''.',model);
end
end

% =====================================================================
function tf=resolveFlag(s,name,dflt)
tf=dflt;
if isfield(s,name) && ~isempty(s.(name)), tf=logical(s.(name)); end
end

% =====================================================================
function v=requireScalar(s,name,what)
%requireScalar  A required numeric setting, named in the error rather than found
%   missing three lines later inside an arithmetic expression.
if ~isfield(s,name) || isempty(s.(name))
    error('fitNVC:missingSetting', ...
        's.%s is required - %s, in seconds.',name,what);
end
v=double(s.(name));
if ~(isscalar(v)&&isfinite(v))
    error('fitNVC:badSetting','s.%s must be a finite scalar (%s).',name,what);
end
end
