%getBolusMoments - First and second moments of a bolus curve, WITHOUT differentiating it
%
%   A bolus curve C(t) is the arterial input convolved with a unit-area kernel k:
%   C = Ca * k.  Under convolution the CUMULANTS ADD, so
%
%       mean(C')  = mean(Ca') + mean(k)
%       var (C')  = var (Ca') + var (k)
%
%   and every transit-time quantity anybody wants is a difference between the tissue
%   curve's moments and the input's.  THAT IS THE WHOLE OF THE MODEL-FREE METHOD, and it
%   is why this function exists: the moments belong to C', the derivative, and
%   differentiating a noisy fluorescence trace to take a moment of it is the one thing
%   that must not happen.
%
%   INTEGRATE BY PARTS INSTEAD.  With y = C - Bl and y(0) = 0,
%
%       INT t*y'(t) dt  = T*y(T) - INT y dt
%       INT t^2*y'(t)dt = T^2*y(T) - 2*INT t*y dt
%
%   so both moments are PLAIN INTEGRALS OF THE CURVE:
%
%       m1  = T - INT y dt / yT
%       m2  = T^2 - 2*INT t*y dt / yT
%       var = m2 - m1^2
%
%   No derivative, no smoothing choice, and the noise is averaged rather than amplified.
%   Verified against a rectangle (mean W/2, sd W/sqrt(12)) and against a rectangle
%   convolved with an exponential (moments add): both to better than 2 %.
%
%   y' MAY BE NEGATIVE AND THE ALGEBRA STILL HOLDS.  A bolus curve overshoots its own
%   plateau, so y' is positive on the rise and negative on the fall; the "distribution"
%   whose moments these are is a signed measure.  Cumulants add for signed measures too,
%   as long as the total masses are non-zero - and yT, the total mass, is the plateau
%   step, which is large and positive.  What DOES break is the interpretation of var as a
%   square: a curve whose input is wider than itself returns var < 0, which is the
%   correct report that the measurement resolved nothing, and the caller must not sqrt it
%   into a small positive number.
%
%   yT IS THE PLATEAU, MEASURED OVER A WINDOW.  Reading y(T) off the last sample makes
%   every moment inherit that sample's noise; the mean over the last wEnd seconds does
%   not.  The small inconsistency this introduces (the integral runs to T, the level is
%   an average ending at T) is worth it and is stated rather than hidden.
%
%   THE CURVE MUST HAVE REACHED ITS PLATEAU BY T.  If it has not, m1 is biased towards T
%   and var with it, and no arithmetic here can tell.  getBolusCthFloor is what calibrates
%   how far that bias reaches on a given recording, and getBolusConfidence's settling
%   factor is what grades it per unit; neither lives here.
%
% Syntax:
%    o = getBolusMoments(t,y,wEnd)
%
% Inputs:
%    t     - [nT x 1] time base, seconds, starting at 0 at the first sample.
%    y     - [nT x 1] BASELINE-REFERENCED curve, so y(1) ~ 0.
%    wEnd  - width of the end window the plateau level is averaged over, seconds.
%
% Outputs:
%    o - struct with .yT (the plateau step), .m1 (first moment, s), .m2 (second raw
%        moment, s^2), .var (central second moment, s^2, MAY BE NEGATIVE) and .sd
%        (sqrt of var, NaN when var < 0).
%
% Example:
%    o = getBolusMoments(results.time, trace-baseline, 1);
%    fprintf('centroid %.2f s, spread %.2f s\n', o.m1, o.sd);
%
% See also: getBolusMetrics, getBolusInput, getBolusCthFloor, getBolusConfidence, trapz
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function o = getBolusMoments(t,y,wEnd)

t=double(t(:)); y=double(y(:));
T=t(end);
o.yT=mean(y(t>=T-wEnd));
o.m1 =T   -  trapz(t,y)     /o.yT;
o.m2 =T^2 - 2*trapz(t,t.*y) /o.yT;
o.var=o.m2-o.m1^2;
% A NEGATIVE VARIANCE IS A RESULT, NOT A ROUNDING ERROR - it says the curve is narrower
% than the input it is being referenced to, i.e. nothing was resolved.  sqrt of it is NaN
% on purpose, so a caller cannot quietly turn "unresolved" into "zero spread".
if o.var>=0, o.sd=sqrt(o.var); else, o.sd=NaN; end
end
