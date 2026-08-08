%meShift - The magnifier: amplify the in-band phase and put the picture back together.
%
%   THE ONE PLACE THE ORDER OF THE PIPELINE IS WRITTEN DOWN, and the order is not
%   arbitrary:
%
%        pyramid -> phase -> TEMPORAL -> GLOBAL -> BLUR -> MASK -> shift -> collapse
%
%   Temporal first, because the bulk motion to be removed is the bulk motion IN THE
%   BAND - the preparation drifts at every frequency and only the part that shares a
%   band with the pulse can be confused with it.  Global before the blur, because a
%   spatial blur of a field that still contains a whole-field offset just spreads the
%   offset.  Mask last, because it is a decision about where to AMPLIFY and must not
%   touch anything that was measured.
%
%   THE GAIN IS 1+alpha, NOT alpha.  Wadhwa et al., SIGGRAPH 2013, Eq. 5: the
%   magnified band is A*exp(i*w*(x + (1+alpha)*delta)).  The band already carries
%   delta; the shift adds alpha*delta on top of it, so what comes out has moved by
%   (1+alpha)*delta.  The BOUND, confusingly, is stated on alpha*delta < lambda/4 -
%   it is the added shift that breaks the feature up, not the total.  If a
%   calibration ever comes out at gain alpha rather than 1+alpha, a DC term has been
%   doubled or dropped somewhere above; chase it, do not fit a correction factor.
%
%   THE SHIFT IS A QUATERNION EXPONENTIAL, real part taken, exactly as the ICCP 2014
%   pseudocode's PhaseShiftCoefficientRealPart.  cos(|p|) times the band, minus the
%   two Riesz components carrying sin(|p|) along the phase's own orientation.  As
%   the phase goes to zero the sin(|p|)/|p| factor goes to one, which is the
%   identity - so alpha = 0 returns the reconstruction, and that is what the
%   calibration measures the raw amplitude on.
%
%   TWO SYNTAXES, BECAUSE THE ANALYSIS IS THE EXPENSIVE HALF AND alpha IS FREE.
%   Everything up to the shift is independent of alpha, so a sweep over alpha - or
%   over which levels are amplified - runs the analysis once and re-shifts.  That is
%   what makes the calibration curve affordable, and it is exact rather than an
%   approximation: nothing above the shift reads alpha.
%
% Syntax:
%    [out, an] = meShift(stack, fs, s)
%    out       = meShift(an, s)
%
% Inputs:
%    stack - [rows columns frames].
%    fs    - frames per second.
%    an    - the analysis returned by a previous call, to re-shift at a new factor.
%    s     - settings.  Fields read:
%            .levels         which pyramid levels are amplified.  Levels outside
%                            this list pass through untouched, which is how a
%                            per-level calibration is made.
%            .alpha          the amplification factor; scalar, or one per level.
%            .nPyrLevels     how many bands to build.  Empty takes max(s.levels).
%            .riesz          '3tap' or 'exact'   read by meRieszPyramid
%            .temporalFilter .passband .filterOrder  read by meTemporal
%            .globalOrder                          read by meGlobal
%            .blurSigma                            read by meBlur
%            .maskMode .mask                       read by meMask
%
% Outputs:
%    out - the magnified stack, single, the size of the input.
%    an  - the analysis, for re-shifting and for inspection:
%          .pyr    the pyramid, bands UNSHIFTED
%          .pc,.ps {1 x nLevels} the in-band phase after global removal and blur -
%                  this is the field the zero-phase assertion is made on
%          .amp    {1 x nLevels} the weights
%          .mask   {1 x nLevels} from meMask
%          .levels the levels that were analysed
%          .fs     frames per second
%
% Example:
%    s = meSettings;  s.passband = [7.9 14.7];  s.alpha = 20;
%    [mag, an] = meShift(stack, 100, s);
%    s.alpha = 50;  mag50 = meShift(an, s);      % no re-analysis
%
% Dependencies: meRieszPyramid, mePhase, meTemporal, meGlobal, meBlur, meMask,
%               meCollapse.
% See also: meRieszPyramid, mePhase, meTemporal, meGlobal, meBlur, meMask,
%           meCollapse, meLinear, meValidate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function [out, an] = meShift(x, varargin)

if isstruct(x)
    an = x;
    s  = varargin{1};
else
    fs = varargin{1};
    s  = varargin{2};
    an = analyse(x, fs, s);
end

levels = intersect(s.levels(:)', an.levels);
alpha  = s.alpha(:)';
if isscalar(alpha)
    alpha = repmat(alpha, 1, numel(levels));
elseif numel(alpha) ~= numel(levels)
    error('meShift:alpha', ...
        'alpha is one number or one per amplified level; %d given for %d levels.', ...
        numel(alpha), numel(levels));
end

pyr = an.pyr;
for i = 1:numel(levels)
    k  = levels(i);
    pc = alpha(i).*an.pc{k};
    ps = alpha(i).*an.ps{k};
    if ~isempty(an.mask{k})
        pc = pc.*an.mask{k};
        ps = ps.*an.mask{k};
    end
    pyr.lap{k} = shiftReal(pyr.lap{k}, pyr.rx{k}, pyr.ry{k}, pc, ps);
end
out = meCollapse(pyr);
end

% =====================================================================
function an = analyse(x, fs, s)
%analyse  Everything that does not depend on alpha.

if ~isa(x,'single') && ~isa(x,'gpuArray'), x = single(x); end

nLev = s.nPyrLevels;
if isempty(nLev), nLev = max(s.levels); end

an        = struct();
an.pyr    = meRieszPyramid(x, nLev, s.levels, s.riesz);
an.levels = an.pyr.levels;
an.fs     = fs;

ph        = mePhase(an.pyr);
an.pc     = cell(1,nLev);
an.ps     = cell(1,nLev);
an.amp    = ph.amp;

for k = an.levels
    pc = meTemporal(ph.cos{k}, fs, s);
    ps = meTemporal(ph.sin{k}, fs, s);
    [pc,ps] = meGlobal(pc, ps, ph.amp{k}, s.globalOrder);
    an.pc{k} = meBlur(pc, ph.amp{k}, s.blurSigma);
    an.ps{k} = meBlur(ps, ph.amp{k}, s.blurSigma);
end

an.mask = meMask(s, an.pyr);
end

% =====================================================================
function L = shiftReal(L, RX, RY, pc, ps)
%shiftReal  Real part of exp(phase) times the quaternion coefficient - the
%   pseudocode's PhaseShiftCoefficientRealPart.  sin(p)/p is one at p = 0, so a
%   zero shift is the identity and alpha = 0 returns the reconstruction exactly.
p  = hypot(pc, ps);
sc = sin(p)./p;
sc(p==0) = 1;
L  = cos(p).*L - (pc.*sc).*RX - (ps.*sc).*RY;
end
%------------- END OF CODE --------------
