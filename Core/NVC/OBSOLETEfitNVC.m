%fitNVC - Stimulus-locked neurovascular-coupling response fit (setup + per-trace fit)
%
%   fitNVC centralises the MODELLED neurovascular-coupling science that runNVC applies
%   to a segment's stimulus-locked trace: once per segment on the average of its
%   epochs, and - when s.nvcFit says so - once per segment per EPOCH beside it.  It
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
%       m = fitNVC(series, layout, s)     % series = [nT x 1], ONE epoch or an average
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
%     t0 CANNOT BE NEGATIVE, bounds [0 10] s.  A vascular response does not precede
%     its stimulus, so a negative latency is never an answer - and with the floor at
%     -2 s, Onset came out negative in 79.8 % of segments of the reference recording,
%     median -0.40 s.  What t0 was fitting there is the model's RISE-TIME mismatch, not
%     a latency: at a 0.25 s sample interval the data wants a rise faster than any
%     bound allows, and t0 and the fast time constant trade against each other along a
%     valley.  Bounding it does not cure the degeneracy - it makes it VISIBLE, as an
%     Onset railed at 0 with Identified 0 beside it, where -0.40 s looked like a
%     measurement.  Re-parameterising the kernel does not cure it either: the two decay
%     constants as free parameters directly give R2 0.7913 against 0.7991 and raise the
%     bound-hit rate from 51.5 % to 88.2 % (EVIDENCE.md 4).  The degeneracy is in the
%     DATA, and a faster recording is the only thing that removes it.
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
%   latency / width, but seven free parameters against the 120 samples of one stimulus
%   repetition is more freedom than the six here already have - and six of them do not
%   identify at 4 Hz.  It belongs in a later session if the author wants it, on the
%   aggregate fit and on a faster recording - do not add it here by accident.
%
%   ANALYSIS LEVELS  (s.segNvcReturn; a cell subset of the two names below).
%   Default (absent or empty) is BOTH.
%       'model'           fitted parameters, derived physiology and fit quality
%                         (runs the fit) - see the scalar list below
%       'reconstruction'  fData [nT x 1], the fitted curve (runs the fit)
%   THERE IS NO PER-PIXEL FIT AT ANY SETTING, so no level here is a per-pixel set:
%   23.2 M fits is 266 hours on the reference recording.
%
%   THIS CORE HAS NO MODEL-FREE MEASUREMENT OF ITS OWN, AND THAT IS THE POINT.  It
%   carried one - a 'markers' level and a private marker routine - until 05-Aug-2026,
%   and it was a second definition of what getNVCMetrics already measures better: with
%   a polarity resolved once per segment, a hold rule on every crossing, and a finale
%   window to measure the return against.  Two definitions of nsPeak in one tree is the
%   kind of thing nobody notices until it is quoted, so one of them went.
%
%   THE MARKER TWINS ARE getNVCMetrics READ OFF THE FITTED CURVE.  PeakFit, TTPFit,
%   AUCbFit, TRiseFit, TDecFit and DurFit are that core called on fData, with THIS
%   trace's settings and THIS segment's resolved polarity - so a twin and its
%   model-free original are the same measurement of two different curves, and any
%   difference between them is a statement about the fit rather than about two
%   definitions.  It is also why s.epochFinaleSec is required here: TRise, TDec and Dur
%   are not defined without a finale window.
%   THE TWIN SET IS CURATED RATHER THAN INHERITED.  getNVCMetrics emits twenty scalars
%   and six of them are twinned; the rest describe the NOISE (BlSd, FinSd, SNR) or the
%   rule that read it (PeakSign, NoRise, NoDec), and a noise statistic of a smooth
%   fitted curve measures the fit's own smoothness, not the recording.
%
%   SCALARS EMITTED AT THE 'model' LEVEL
%     Gain Onset Drift Baseline            A, t0, b1, b0
%     TauS TauF                                       ('secondOrder' only)
%     Zeta Omega0 RingPeriod Damping       DERIVED from TauS/TauF, and GATED - below
%     A1 Beta1 A2 Beta2 CRatio                        ('doubleGamma' only)
%     R2 RMSE AICc NIter Converged StartsAgree Identified      fit quality
%     PeakFit TTPFit AUCbFit TRiseFit TDecFit DurFit  getNVCMetrics on the fitted curve
%     DipGain DipTau DipR2 AICcDip                    (s.nvcDip only)
%   EVERY SCALAR NAME STARTS WITH A CAPITAL, so the caller's prefixed name is a plain
%   concatenation (ns + Gain = nsGain) and there is exactly one spelling of each
%   quantity between this bag, the saved tree and the metrics tables.  No name here is
%   a name getNVCMetrics emits, because both cores write into results.nvc.esMetrics
%   under the same prefix and runNVC refuses a collision rather than overwriting one.
%   with  Zeta = sqrt(tauF)/(2*tauS),  Omega0 = 1/sqrt(tauF),  RingPeriod = 2*pi/wd
%   (NaN when there is no ringing to have a period), and Damping = sigma = 1/(2*tauS),
%   the exponential decay rate in 1/s.
%
%   IDENTIFIED IS THE NUMBER THAT SAYS WHETHER THE OTHERS MEAN ANYTHING.  It is 1 when
%   every BOUNDED parameter of the base model - t0 and the kernel's shape block - came
%   to rest strictly inside its bounds, and 0 when any of them sat ON one.  It is not a
%   nicety: measured on real 4 Hz data, tauS sits on its lower bound in 56.4 % of
%   segments, which means sigma = 1/(2*tauS) is running to infinity and the kernel is
%   being used as "instant jump plus one slow decay".  A railed parameter reported as a
%   plain number is the single failure this core exists to prevent.  A fit that reached
%   no optimum at all leaves Identified NaN rather than claiming either answer.
%
%   AND IT GATES THE DERIVED BLOCK.  Zeta, Omega0, RingPeriod and Damping are all
%   functions of tauS, so when tauS is a bound they are all bound artefacts - Zeta's
%   median on real data is 3.57 and 94 % of segments read as overdamped, entirely
%   because of the railing.  They are therefore emitted ONLY when Identified is 1.
%   TauS and TauF themselves are ALWAYS emitted: what the optimiser did is a fact, and
%   Identified is what says how to read it.  A single stimulus repetition rails far
%   more often than the average of twenty, so in practice this block is the aggregate
%   fit's - which is the level whose SNR supports it.
%
%   ZETA IS WHY THIS MODEL WAS CHOSEN OVER A DESCRIPTIVE ONE.  The undershoot-to-peak
%   ratio equals exp(-pi*zeta/sqrt(1-zeta^2)) in both the impulse and the step limit,
%   which makes zeta far more robust to differences in stimulus duration than a raw
%   undershoot depth: two groups stimulated for different lengths of time can be
%   compared on zeta and cannot be compared on how deep their undershoots went.  On
%   4 Hz flow data it is usually gated away, and a recording fast enough to resolve the
%   rise is what would give it back (MODEL.md).
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
%                epochFinaleSec     [-t1 0] end-of-epoch finale window, s (REQUIRED -
%                                   the marker twins are getNVCMetrics, and TRise,
%                                   TDec and Dur are measured against the finale)
%                nvcModel           'secondOrder' (default) | 'doubleGamma'
%                nvcDip             false (default) | true
%                segNvcReturn       cell subset of {'model','reconstruction'}
%                                   (default: both)
%                nvcStarts          number of multi-start points (default 16)
%                nvcWeights         [nT x 1] per-timepoint weights for a weighted
%                                   least squares (optional; e.g. 1./SEM.^2 across
%                                   the accepted epochs).  Absent = unweighted.
%              plus, read by getNVCMetrics for the marker twins and defaulted there:
%                nvcMedFiltSec nvcHoldSamples nvcDecayFrac nvcPolarity, and
%                nvcPolarityResolved - the sign the PRODUCER resolved for this
%                segment, which is what keeps a twin and its original mirrored the
%                same way.
%
% Outputs:
%    layout  - (SETUP)    stimulus boxcar, windows, bounds, start grid, the twin
%                         core's own layout and .want for one time base, plus
%                         .scalarNames, the ORDERED list of scalars ANALYSIS emits
%                         (the caller needs no other knowledge of the field set).
%    m       - (ANALYSIS) flat, prefix-free, all-single bag for one trace: exactly
%                         layout.scalarNames, plus fData when 'reconstruction' was
%                         asked for, plus the logical m.valid.
%
% Example:
%    s.epochStimStartSec=10; s.stimDurationSec=5; s.epochBaselineSec=[0 10];
%    s.epochFinaleSec=[-5 0];
%    layout = fitNVC(results.nvc.time, s);
%    m      = fitNVC(epochAverageTrace, layout, s);
%    if m.Identified
%        fprintf('zeta = %.2f, gain = %.3g\n', m.Zeta, m.Gain);
%    else
%        fprintf('gain = %.3g; the time constants sat on a bound\n', m.Gain);
%    end
%
% See also: runNVC, getNVCMetrics, fitVasoreactivity, quasiUniform,
%           getPulsatilityMetrics, getVasomotionMetrics, lsqcurvefit
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 05-August-2026

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
%   The per-trace call is then a bounded least-squares and one pass of the twin core,
%   and nothing else - which is what makes 1938 segments x 20 epochs cost 27 minutes
%   rather than a day.
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
% product runNVC cuts).  Anchoring to time(1) rather than assuming 0
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

% The response window - everything from the stimulus onward - is where the amplitude
% start point is read off the trace.  Stored as indices so no per-trace call does a
% find(), and checked here because a fit with no response window in it is a fit of the
% baseline.
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

% ---- the marker twins' core, on the SAME time base -------------------------
% getNVCMetrics is the library's model-free measurement, and the twins ARE it, read off
% the fitted curve.  Its layout is built here rather than handed in, so the two cores
% cannot end up with different windows for the same trace: same time base, same s, one
% call.  Its LEVELS are the twins' own and are deliberately not inherited -
% s.segNvcReturn is spelled in this core's vocabulary, and 'model' is not one of that
% core's levels.  Only the model level needs the twins, so a caller who wants nothing
% but the reconstruction still needs no finale window.
L.twinNames=nvcTwinNames();
L.twin=[];
if L.want.model
    sTwin=s;
    sTwin.segNvcReturn={'markers','timing','areas'};
    L.twin=getNVCMetrics(time,sTwin);
end

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

% t0 is the response onset relative to the stimulus, and it CANNOT BE NEGATIVE - see
% the header.  A response does not precede its stimulus, and with the floor at -2 s
% that is exactly what 79.8 % of the reference recording's segments reported, because
% t0 was absorbing a rise-time mismatch and calling it a latency.  Railing at 0 with
% Identified 0 beside it is a visible failure; -0.40 s is an invisible one.
t0Lo=0; t0Hi=10;

% HOW CLOSE TO A BOUND COUNTS AS ON IT (the Identified rule).  A distance from a bound
% is relative TO THAT BOUND, so one number means the same thing for a 0.05 s and a 20 s
% time constant.  1 % is far finer than any difference between two time constants a
% reader would act on, and coarse enough to catch an optimiser that stopped just short.
L.boundTol=1e-2;

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
function m=nvcAnalyze(series,L,s)
%nvcAnalyze  The fit, and getNVCMetrics read off the curve it produced, for ONE trace.
%   s is read for the marker twins alone - the polarity the PRODUCER resolved for this
%   segment travels in it, and everything else about the fit is baked into `layout` at
%   SETUP.
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
if ~(L.want.model || L.want.reconstruction), return, end

% ---- the fit ---------------------------------------------------------------
% THE CURVE IS SINGLE FROM THE MOMENT IT EXISTS, which is the one place the storage
% type is load-bearing rather than merely conventional: the twins below are measured on
% THIS curve, so running getNVCMetrics again on a saved fData reproduces PeakFit and
% its five companions exactly, rather than in all but the last bits.  R2, RMSE and
% AICc are the optimiser's own and come from multiStartFit in double.
P0=fillDataStarts(L.starts,series,L);
[p,q]=multiStartFit(series,L,P0,L.lb,L.ub,false);
fData=single(nvcEval(p,L,false));

if L.want.reconstruction
    m.fData=fData;
end
if ~L.want.model, return, end

m.Gain=single(p(L.iA)); m.Onset=single(p(L.iT0));
m.Drift=single(p(L.iB1)); m.Baseline=single(p(L.iB0));

% ---- did the optimiser MEASURE the shape, or did it run into a wall? -------
% Only the BOUNDED parameters can be at a bound - b0, b1 and A are free, and the dip
% pair belongs to the model fitted beside this one, never to it.  A fit that reached no
% optimum at all leaves Identified NaN: it has nothing to test, and either answer would
% be a claim about a fit that does not exist.
bnd=[L.iT0 L.iShape];
identified=false;
if all(isfinite(p(bnd)))
    identified=~anyAtBound(p,L.lb,L.ub,bnd,L.boundTol);
    m.Identified=single(identified);
end

switch L.model
    case 'secondOrder'
        tauS=p(L.iShape(1)); tauF=p(L.iShape(2));
        m.TauS=single(tauS); m.TauF=single(tauF);
        % THE DERIVED BLOCK IS GATED BY Identified (see the header).  All four are
        % functions of tauS, so a railed tauS makes all four of them a description of
        % the bound rather than of the vessel - sigma = 1/(2*tauS) running to infinity
        % - and a NaN is the honest report of that, where a number is not.
        if identified
            sigma=1/(2*tauS); d=1/tauF-sigma^2;     % d = wd^2
            m.Zeta   =single(sqrt(tauF)/(2*tauS));  % damping RATIO (dimensionless)
            m.Omega0 =single(1/sqrt(tauF));         % natural frequency, rad/s
            m.Damping=single(sigma);                % decay RATE, 1/s
            % a critically- or over-damped response does not ring, so it has no ring
            % period - NaN rather than Inf, which would look like a very slow ring
            if d>0, m.RingPeriod=single(2*pi/sqrt(d)); end
        end
    case 'doubleGamma'
        m.A1=single(p(L.iShape(1))); m.Beta1 =single(p(L.iShape(2)));
        m.A2=single(p(L.iShape(3))); m.Beta2 =single(p(L.iShape(4)));
        m.CRatio=single(p(L.iShape(5)));
end

m.R2=single(q.r2); m.RMSE=single(q.rmse); m.AICc=single(q.aicc);
m.NIter=single(q.nIter); m.Converged=single(q.converged);
m.StartsAgree=single(q.nAgree);

% ---- the marker twins: getNVCMetrics, on the FITTED curve ------------------
% ONE core, two curves.  The polarity this segment resolved, the filter length, the
% hold rule and the finale window all arrive in s and in L.twin, so the twin is
% mirrored and windowed exactly as the model-free original was, and the pair can be
% read side by side as a statement about the fit.
mf=getNVCMetrics(fData,L.twin,s);
for k=1:numel(L.twinNames)
    m.([L.twinNames{k} 'Fit'])=mf.(L.twinNames{k});
end

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
function tf=anyAtBound(p,lb,ub,idx,relTol)
%anyAtBound  Did any of the bounded parameters come to rest ON one of its bounds?
%   A DISTANCE FROM A BOUND IS RELATIVE TO THAT BOUND, which is what lets one
%   tolerance mean the same thing for a 0.05 s and a 20 s time constant.  A bound of
%   ZERO has no relative distance to be within, and there the width of the parameter's
%   own range is the only scale the bounds themselves offer - for t0, whose floor is 0,
%   that makes the band 0.1 s, well under one sample of a 4 Hz recording.
tf=false;
for k=idx
    lo=lb(k); hi=ub(k); width=hi-lo;
    if isfinite(lo) && abs(p(k)-lo)<=relTol*boundScale(lo,width), tf=true; return, end
    if isfinite(hi) && abs(p(k)-hi)<=relTol*boundScale(hi,width), tf=true; return, end
end
end

function sc=boundScale(b,width)
sc=abs(b);
if sc==0, sc=abs(width); end
end

% =====================================================================
function names=nvcTwinNames()
%nvcTwinNames  The getNVCMetrics scalars that are twinned onto the fitted curve.
%   CURATED, NOT INHERITED (see the header): the amplitude, its time, the two crossings
%   and the width and area between them are all properties a fitted curve HAS.  The
%   rest of that core's set describes the noise (BlSd, FinSd, SNR) or the rule that
%   read it (PeakSign, NoRise, NoDec), and neither is a property of a smooth curve.
names={'Peak','TTP','AUCb','TRise','TDec','Dur'};
end

% =====================================================================
function names=nvcScalarNames(model,dip,want)
%nvcScalarNames  The ORDERED scalar contract for one layout.
%   ANALYSIS emits exactly these names and no others, so the wrapper can build its
%   result columns from this list and never has to know which model ran.  The list
%   lives here, beside the code that fills it, and nowhere else.
names={};
if want.model
    names=[names,{'Gain','Onset','Drift','Baseline'}];
    switch model
        case 'secondOrder'
            names=[names,{'TauS','TauF','Zeta','Omega0','RingPeriod','Damping'}];
        case 'doubleGamma'
            names=[names,{'A1','Beta1','A2','Beta2','CRatio'}];
    end
    names=[names,{'R2','RMSE','AICc','NIter','Converged','StartsAgree','Identified'}];
    names=[names,strcat(nvcTwinNames(),'Fit')];
    if dip
        names=[names,{'DipGain','DipTau','DipR2','AICcDip'}];
    end
end
end

% =====================================================================
function want=nvcParseReturn(s)
%nvcParseReturn  Resolve s.segNvcReturn into per-level compute logicals.
%   Absent OR empty gives the documented default, both levels.  THE VOCABULARY IS THIS
%   CORE'S, not getNVCMetrics': the two cores read the same field name out of the same
%   settings struct, so a caller that wants the fit says {'model'} and this core hands
%   the other one its own levels itself.
levels={'model','reconstruction'};
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
