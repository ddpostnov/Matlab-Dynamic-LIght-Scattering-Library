%getPulsatilityMetrics  Shared harmonic pulsatility core (fit setup + per-cycle analysis)
%
%   getPulsatilityMetrics centralises the harmonic-fit pulsatility "science"
%   that runPulsatility applies to every segment cycle and every pixel cycle.
%   It builds the multi-harmonic sine model ONCE per time base (compiling the
%   fittype is the one non-trivial per-time-base cost) and reuses it for every
%   waveform on that base.  Two call modes are dispatched by argument count:
%
%   SETUP  (nargin 2)
%       layout = getPulsatilityMetrics(time, s)
%     Builds the fit grids and the harmonic model implied by `time` and s.nHarm,
%     and resolves s.segPulsReturn once.  Returns a small `layout` struct built
%     ONCE per time base and reused for every waveform on that base:
%       .time     [nT x 1] one-cycle time base (= time(:)); the marker-time source
%       .T        cycle duration = time(end)+time(1).  ONE expression for the TWO
%                 sample conventions this library currently writes, not a fallback:
%                 runContrastInternalCycle samples the cycle at its ENDPOINTS
%                 (time(1) is exactly 0, so T is time(end)), runIntensityInternalCycle
%                 at its BIN CENTRES (time(1) is half a bin, so the period is
%                 time(end) plus that half bin).  Both producers are live; on a
%                 25-bin cardiac cycle taking time(end) alone under-states the period
%                 by two per cent and moves symRatio and the timeMin wrap with it,
%                 silently.  The speckle side is unchanged to the bit.
%       .x        [nT*5 x 1] fit abscissa linspace(0,5,nT*5) (five tiled cycles)
%       .xx       [nT x 1]   reconstruction abscissa linspace(0,1,nT) (one cycle)
%       .nHarm    number of harmonics (s.nHarm, default 5)
%       .model    prebuilt fittype  sum_{n=1..nHarm} a_n*sin(2*pi*n*x + b_n) + c
%                 with Trust-Region NLS options (amps>=0, phases in [-pi,pi],
%                 offset>=0); honours an optional s.fitOptions override
%       .want     per-level compute flags from s.segPulsReturn (see below)
%
%   ANALYSIS  (nargin 3)
%       m = getPulsatilityMetrics(series, layout, s)   % series = [nT x 1] cycle
%     The science for ONE averaged-cycle waveform, MODULAR via s.segPulsReturn:
%     only the requested levels' quantities are returned.  The fit is run iff the
%     'model' or 'reconstruction' level is selected; when it runs, the model-free
%     markers are computed on the reconstruction fData (reproducing runPulsatility,
%     where the fit overwrites the cycle before the markers are taken), otherwise
%     on the raw cycle.  Returns a FLAT, prefix-free bag `m` (the caller applies
%     the ps/pd branch prefix); every numeric field is SINGLE.
%
%   ANALYSIS levels (s.segPulsReturn; a cell subset of the THREE names below).
%   Default (absent/empty) is the COMPLETE set {'markers','model','reconstruction'}.
%   The per-pixel path reuses SETUP with s.segPulsReturn := s.ppxPulsReturn.
%       'markers'         model-free scalars, ALWAYS computed on the marker source W
%                         (fData if the fit ran, else series):
%                           PI  = (max-min)/mean    (Gosling pulsatility index)
%                           RI  = (max-min)/max     (Pourcelot resistivity index)
%                           mean std min max         (cycle mean/fluctuation/extrema)
%                           timeMin timeMax          (times of the extrema; timeMin
%                                                     is wrap-shifted by -T when it
%                                                     would otherwise trail timeMax)
%                           symRatio = (T-ascend)/ascend, ascend = timeMax-timeMin
%                                                    (descent/ascent asymmetry; NEW)
%       'model'           harmonic coefficients (the fit is run):
%                           hAmp   [1 x nHarm] amplitudes a_1..a_nHarm (>=0)
%                           hPhase [1 x nHarm] phases b_1..b_nHarm, each per-element
%                                              wrapped into [0,2*pi)
%                           R2     fit R^2 (gof.rsquare)
%                         (the fit constant c is NOT returned: for a full-cycle
%                          harmonic series c == mean(fData), so the marker `mean`
%                          already carries it.)
%       'reconstruction'  fData [nT x 1] model reconstruction feval(fit, xx).
%
%     Returns m.valid = all(isfinite(series)).  For a non-finite (invalid) series
%     the fit is skipped and fData/hAmp/hPhase/R2 are NaN; the markers are still
%     computed on the NaN cycle, so - exactly as in runPulsatility - min/max/mean/
%     PI/RI come out NaN while timeMin/timeMax fall back to time(1) (argmax/argmin
%     of an all-NaN column is index 1).  This keeps the retained anchors bit-
%     identical to the pre-refactor output.
%
%   INPUTS
%     time    [nT x 1] one-cycle time base (SETUP); sets the fit grids and T.
%     series  [nT x 1] one averaged-cycle waveform (ANALYSIS).
%     layout  struct from SETUP mode.
%     s       parameter struct.  Fields:
%               nHarm         number of harmonics (default 5)
%               segPulsReturn cell subset of {'markers','model','reconstruction'}
%                             (default: all three); the per-pixel caller passes
%                             s.ppxPulsReturn through this field
%               fitOptions    (optional) fitoptions override for the model; when
%                             absent the internal Trust-Region preset is used
%                             (amps [0,Inf], phases [-pi,pi], offset [0,Inf];
%                             StartPoint amps 0.9,0.8,.. phases 0.1,0.2,.. c=100)
%
%   OUTPUTS
%     layout  (SETUP)    fit grids + harmonic model + .want for one time base.
%     m       (ANALYSIS) flat, prefix-free, all-single marker bag for one cycle
%                        (selected levels only, plus logical m.valid).
%
%   DEPENDS ON
%     MATLAB Curve Fitting Toolbox (fittype/fitoptions/fit) - ONLY when the 'model' or
%     'reconstruction' level is selected.  A markers-only call needs no toolbox.
%
%   THE CORE IS AGNOSTIC ABOUT WHAT THE WAVEFORM MEASURES, and that is what lets two
%   wrappers share it: runPulsatility applies a ps/pd prefix to what comes back and
%   runIntensityPulsatility a pv one.  Nothing here knows or asks.
%
% See also: runPulsatility, runIntensityPulsatility, getVasomotionMetrics,
%           runContrastInternalCycle, runIntensityInternalCycle
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

function out = getPulsatilityMetrics(arg1,arg2,arg3)
if nargin==2
    out=pulsSetup(arg1,arg2);          % SETUP:    (time, s)             -> layout
elseif nargin==3
    out=pulsAnalyze(arg1,arg2,arg3);   % ANALYSIS: (series, layout, s)   -> m
else
    error('getPulsatilityMetrics:nargin', ...
        ['Use layout=getPulsatilityMetrics(time,s) for SETUP or ' ...
         'm=getPulsatilityMetrics(series,layout,s) for ANALYSIS.']);
end
end

% =====================================================================
function layout=pulsSetup(time,s)
%pulsSetup  Fit grids and harmonic model for one time base.
%   Builds the fit abscissae (five tiled cycles for fitting, one cycle for the
%   reconstruction), the prebuilt fittype/fitoptions for the nHarm-harmonic sine
%   model, and resolves s.segPulsReturn once so every waveform on this time base
%   computes the same set of levels.
nHarm=resolveNHarm(s);
layout.time =time(:);
%THE PERIOD, FOR BOTH SAMPLE CONVENTIONS, IN ONE EXPRESSION.  An endpoint-sampled cycle
%starts at t=0 and its last sample IS the period; a centre-sampled one starts half a bin
%in and its last sample is half a bin short of it.  time(end)+time(1) is exact for both
%and reduces to the old time(end) on the speckle side, where time(1) is exactly zero -
%so no golden moves.  It is NOT a back-compatibility branch: both producers are current
%(runContrastInternalCycle, runIntensityInternalCycle) and 02-cardiac.md records why the
%centres convention is REQUIRED on '_c_I' - with endpoints the first and last bin hold
%the same phase and a mean over the bins double-counts the cardiac foot, which is
%exactly the quantity [D11] needs to equal the mean image.
layout.T    =time(end)+time(1);
layout.x    =linspace(0,5,numel(time)*5)';     % five tiled cycles (fit domain)
layout.xx   =linspace(0,1,numel(time))';       % one cycle (reconstruction domain)
layout.nHarm=nHarm;
layout.model=buildPulsModel(nHarm,s);
layout.want =pulsParseReturn(s);
end

% =====================================================================
function nHarm=resolveNHarm(s)
%resolveNHarm  Default/validate s.nHarm (number of harmonics in the sine model).
%   Default 5 (the pre-refactor five-harmonic equation); a positive integer scalar.
if ~isfield(s,'nHarm') || isempty(s.nHarm)
    nHarm=5;
    return
end
nHarm=s.nHarm;
if ~(isnumeric(nHarm)&&isscalar(nHarm)&&isfinite(nHarm)&&nHarm>=1&&nHarm==round(nHarm))
    error('getPulsatilityMetrics:nHarm', ...
        's.nHarm must be a positive integer scalar (number of harmonics).');
end
nHarm=double(nHarm);
end

% =====================================================================
function model=buildPulsModel(nHarm,s)
%buildPulsModel  Prebuilt fittype for sum_{n=1..nHarm} a_n*sin(2*pi*n*x+b_n)+c.
%   Generalises runPulsatility's frozen five-harmonic preset to nHarm harmonics.
%   For nHarm=5 the equation string, bounds and StartPoint are bit-identical to
%   the old launcher preset, so the fit (and thus every coefficient) reproduces
%   it exactly.  Coefficient ordering assumes single-digit harmonic indices
%   (coeffnames sorts a1..a_nHarm, b1..b_nHarm, c alphabetically for nHarm<=9),
%   the same assumption the old code makes.  An optional s.fitOptions overrides
%   the internal preset for power users.
terms=cell(1,nHarm);
for n=1:nHarm
    terms{n}=sprintf('a%d*sin(%d*pi*x+b%d)',n,2*n,n);   % 2*n*pi*x written as 2,4,6,..
end
eqStr=[strjoin(terms,'+') '+c'];

if isfield(s,'fitOptions') && ~isempty(s.fitOptions)
    opt=s.fitOptions;                                    % power-user override
else
    lb   =[zeros(1,nHarm)  , -pi*ones(1,nHarm), 0  ];    % a_n>=0, b_n>=-pi, c>=0
    ub   =[Inf(1,nHarm)    ,  pi*ones(1,nHarm), Inf];    % a_n free, b_n<=pi, c free
    amps =(10-(1:nHarm))./10;                            % 0.9,0.8,.. ; k/10 == literal
    phase=(1:nHarm)./10;                                 % 0.1,0.2,..
    start=[amps, phase, 100];                            % [a_n, b_n, c] StartPoint
    opt=fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region', ...
        'Display','off','Lower',lb,'Upper',ub,'StartPoint',start);
end
model=fittype(eqStr,'options',opt);
end

% =====================================================================
function want=pulsParseReturn(s)
%pulsParseReturn  Resolve s.segPulsReturn into per-level compute logicals.
%   Returns a struct with logical fields .markers .model .reconstruction naming
%   which analysis levels to compute.
%     * s.segPulsReturn present -> use it (validated cell subset of the three levels).
%     * s.segPulsReturn absent   -> documented default (the complete set, all three).
levels={'markers','model','reconstruction'};
if isfield(s,'segPulsReturn') && ~isempty(s.segPulsReturn)
    sel=s.segPulsReturn;
    if ischar(sel)||isstring(sel), sel=cellstr(sel); end
    if ~iscellstr(sel)
        error('getPulsatilityMetrics:segPulsReturn', ...
            's.segPulsReturn must be a cell array of level names.');
    end
    bad=~ismember(sel,levels);
    if any(bad)
        error('getPulsatilityMetrics:segPulsReturn', ...
            'Unknown s.segPulsReturn level(s): %s. Valid levels: %s.', ...
            strjoin(sel(bad),', '),strjoin(levels,', '));
    end
else
    sel=levels;
end
for i=1:numel(levels), want.(levels{i})=ismember(levels{i},sel); end
end

% =====================================================================
function m=pulsAnalyze(series,layout,~)
%pulsAnalyze  Per-cycle pulsatility analysis for one averaged-cycle waveform.
%   (s is accepted for signature symmetry with the vasomotion core but unused:
%   every per-cycle parameter is baked into `layout` at SETUP.)
%   Fits the nHarm-harmonic model (when 'model'/'reconstruction' is selected) and
%   reduces the cycle to the model-free markers, always computed on the marker
%   source W = fData (if the fit ran) else the raw series - the whole-matrix path
%   runPulsatility runs after its (destructive) fit overwrite, here per waveform.
%   Every emitted numeric field is single; the markers are computed uniformly on W
%   (no valid-branch) so an invalid all-NaN cycle reproduces the old output exactly
%   (extrema/PI/RI/mean NaN; timeMin/timeMax = time(1)).
want=layout.want;
series=series(:);
m.valid=all(isfinite(series));

runFit=want.model || want.reconstruction;

% ---- harmonic fit (skipped for an invalid cycle; NaN coeffs/fData like the old else) ----
if runFit
    if m.valid
        tsRep=double(repmat(series,5,1));               % five tiled cycles, double for fit
        [f,gof]=fit(layout.x,tsRep,layout.model);
        fData=single(feval(f,layout.xx));               % [nT x 1] reconstruction (single)
        coefs=single(coeffvalues(f));                   % [a_1..a_nHarm, b_1..b_nHarm, c]
        if want.model
            m.hAmp=coefs(1:layout.nHarm);
            hPhase=coefs(layout.nHarm+1:2*layout.nHarm);
            hPhase(hPhase<0)=hPhase(hPhase<0)+2*pi;      % per-element wrap into [0,2*pi)
            m.hPhase=hPhase;
            m.R2=single(gof.rsquare);
        end
    else
        fData=nan(numel(layout.xx),1,'single');
        if want.model
            m.hAmp  =nan(1,layout.nHarm,'single');
            m.hPhase=nan(1,layout.nHarm,'single');
            m.R2    =single(NaN);
        end
    end
    if want.reconstruction
        m.fData=fData;
    end
    W=fData;                                             % markers on the reconstruction
else
    W=series;                                            % markers on the raw cycle
end

% ---- model-free markers (uniform on W; identical expressions to runPulsatility) ----
[maxVal,maxIdx]=max(W,[],1);
[minVal,minIdx]=min(W,[],1);
tMax=layout.time(maxIdx);
tMin=layout.time(minIdx);
if tMin>tMax, tMin=tMin-layout.T; end                   % wrap-shift so tMin<=tMax
meanW=mean(W,1,'omitnan');
ascend =tMax-tMin;
descend=layout.T-ascend;

m.mean    =single(meanW);
m.std     =single(std(W,0,1,'omitnan'));
m.min     =single(minVal);
m.max     =single(maxVal);
m.PI      =single((maxVal-minVal)./meanW);
m.RI      =single((maxVal-minVal)./maxVal);
m.timeMin =single(tMin);
m.timeMax =single(tMax);
m.symRatio=single(descend./ascend);
end
