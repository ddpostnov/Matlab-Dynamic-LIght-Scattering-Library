%getTauC - Simplified decorrelation time from an autocorrelation stack.
%
% Syntax:
%    [tauC, decorThreshold] = getTauC(g2)
%    [tauC, decorThreshold] = getTauC(g2, 'Name', Value, ...)
%
% Description:
%    Estimates the per-pixel decorrelation time tauC of a normalized
%    intensity autocorrelation stack g2 (see getNormalizedG2) by thresholding,
%    without fitting a model. For every pixel the first lag at which g2 drops
%    to or below a threshold is found and refined by linear interpolation
%    between the two bracketing lags. tauC is a fast initial estimate; use it
%    to seed and order fitDLSI.
%
%    The threshold can be a per-pixel 1/e^2 decay level (default), its
%    spatially smoothed variant, a single level taken from the last frame, or
%    a fixed user value (see 'method').
%
% Inputs:
%    g2 - Normalized intensity autocorrelation, 3-D [Y, X, nLags], with lag 0
%         in g2(:,:,1) (as returned by getNormalizedG2).
%
% Optional Name-Value Pair Arguments:
%    'method'      - (String) Threshold definition:
%                    'decay'         (Default) per-pixel 1/e^2 level:
%                                    (g2_0 - g2_min)/e^2 + g2_min.
%                    'decayFiltered' as 'decay' but on spatially smoothed
%                                    first/last frames (uses 'filterSigma').
%                    'lastFrame'     single level = max of the smoothed last
%                                    frame, applied to every pixel.
%                    'fixed'         constant level given by 'threshold'.
%    'threshold'   - (Numeric) Level used when 'method' is 'fixed'.
%    'filterSigma' - (Numeric) Gaussian sigma for the smoothed methods.
%                    Default: 1.
%    'fps'         - (Numeric) Frame rate, used to convert lag index to time
%                    when 'lags' is not given. Default: 1 (tauC in lag units).
%    'lags'        - (Numeric) Explicit lag-time vector, length size(g2,3).
%                    Overrides 'fps'. Default: [] ((0:nLags-1)/fps).
%
% Outputs:
%    tauC           - Decorrelation-time map [Y, X], in seconds when 'fps' or
%                     'lags' is given (otherwise in lag units). Pixels that
%                     never cross the threshold are set to the longest lag.
%    decorThreshold - The threshold map [Y, X] actually used.
%
% Examples:
%    % 1. Default 1/e^2 decay estimate, tauC in seconds:
%    tauC = getTauC(g2, 'fps', 24000);
%
%    % 2. Fixed threshold at 1.2:
%    tauC = getTauC(g2, 'method', 'fixed', 'threshold', 1.2, 'fps', 24000);
%
% Dependencies: Image Processing Toolbox (only for the smoothed methods).
% See also: getNormalizedG2, fitDLSI
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 22-July-2026

%------------- BEGIN CODE --------------
function [tauC, decorThreshold] = getTauC(g2, varargin)

p = inputParser;
p.KeepUnmatched = false;
addRequired(p, 'g2', @(x) isnumeric(x) && ndims(x)==3);
addParameter(p, 'method', 'decay', @(x) any(validatestring(x, ...
    {'decay','decayFiltered','lastFrame','fixed'})));
addParameter(p, 'threshold', [], @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'filterSigma', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'fps', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'lags', [], @(x) isempty(x) || isvector(x));
parse(p, g2, varargin{:});

sizeY = size(g2,1);
sizeX = size(g2,2);
nLags = size(g2,3);
method      = validatestring(p.Results.method, {'decay','decayFiltered','lastFrame','fixed'});
filterSigma = p.Results.filterSigma;

% Lag-time vector (seconds when fps/lags supplied, otherwise lag units).
if isempty(p.Results.lags)
    lagTimes = (0:nLags-1) / p.Results.fps;
else
    if numel(p.Results.lags) ~= nLags
        error('getTauC:lagsLength', ...
            '''lags'' must have size(g2,3) = %d elements.', nLags);
    end
    lagTimes = p.Results.lags(:).';
end

% --- Threshold map ---------------------------------------------------------
decorThreshold = zeros(sizeY, sizeX);
switch method
    case 'lastFrame'
        frame = imgaussfilt(g2(:,:,end), filterSigma);
        decorThreshold(:) = max(frame(:));
    case 'decay'
        g2min = min(g2, [], 3);
        decorThreshold = (g2(:,:,1) - g2min)./exp(2) + g2min;
    case 'decayFiltered'
        firstF = imgaussfilt(g2(:,:,1),   filterSigma);
        lastF  = imgaussfilt(g2(:,:,end), filterSigma);
        decorThreshold = (firstF - lastF)./exp(2) + lastF;
    case 'fixed'
        if isempty(p.Results.threshold)
            error('getTauC:missingThreshold', ...
                'Provide ''threshold'' when ''method'' is ''fixed''.');
        end
        decorThreshold(:) = p.Results.threshold;
end

% --- First threshold crossing per pixel, with linear interpolation ---------
g2r  = reshape(g2, sizeY*sizeX, nLags);
thr  = decorThreshold(:);
below = g2r <= thr;
cross = any(below, 2);
[~, idx] = max(below, [], 2);   % first crossing lag index (1 for no-cross rows)
idx(~cross) = nLags;            % never crossing -> longest lag

tauC = zeros(sizeY*sizeX, 1);
lagT = lagTimes(:);

% Pixels crossing at the very first lag, and pixels never crossing.
tauC(cross & idx==1) = lagT(1);
tauC(~cross)         = lagT(nLags);

% Pixels crossing after the first lag: interpolate between idx-1 and idx.
mInterp = cross & idx>1;
rows = find(mInterp);
ii   = idx(mInterp);
linR = sub2ind([sizeY*sizeX, nLags], rows, ii);
linL = sub2ind([sizeY*sizeX, nLags], rows, ii-1);
yR = g2r(linR);
yL = g2r(linL);
den = yR - yL;
cR = (thr(mInterp) - yL) ./ den;
cR(den==0) = 0;                 % flat segment -> take the left node
tauC(mInterp) = lagT(ii-1).*(1-cR) + lagT(ii).*cR;

tauC = reshape(tauC, sizeY, sizeX);
end
