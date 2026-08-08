%meBlur - Smooth the phase where it means something and ignore it where it does not.
%
%   NOT OPTIONAL, and it is the one line of the pipeline that is about the data
%   rather than the mathematics.  Phase is only defined where there is something to
%   have a phase: at a pyramid coefficient of near-zero amplitude - flat tissue,
%   the middle of the lumen - the measured phase is the direction of that pixel's
%   photon noise, and it is uniformly distributed.  A plain Gaussian blur would
%   average that garbage into the wall beside it.  Weighting by the amplitude is
%   what stops it, and it is Eq. 17 of Wadhwa et al., SIGGRAPH 2013:
%
%        blurred = (phase * amplitude  convolved with K)  /  (amplitude  convolved with K)
%
%   Where the neighbourhood carries signal this is close to a plain blur; where it
%   carries none the numerator and the denominator are both small and the answer is
%   the phase of whatever signal was nearby, not the noise underfoot.
%
%   IT IS ALSO WHERE THE SENSITIVITY COMES FROM.  Averaging n independent pixels
%   divides the phase noise by root n while leaving a wall's phase - coherent over
%   the whole wall - alone.  On a shot-noise-limited fluorescence recording that is
%   not a refinement, it is a large part of why anything is visible at all.
%
%   THE SIGMA IS IN PIXELS OF THE LEVEL, not of the frame, so the same number means
%   the same fraction of a wavelength at every level.  Two pixels at level 3 is
%   eight pixels of the original frame and half a wavelength of that level's band.
%
% Syntax:
%    p = meBlur(p, amp, sigma)
%
% Inputs:
%    p     - the phase component, [rows columns frames].
%    amp   - the weight, same size, from mePhase.
%    sigma - standard deviation of the Gaussian, in pixels OF THIS LEVEL.  Zero or
%            negative returns p untouched.
%
% Outputs:
%    p     - the same size and class, blurred.
%
% Example:
%    pc = meBlur(pc, ph.amp{3}, 2);
%
% Dependencies: Image Processing Toolbox (imfilter).
% See also: mePhase, meGlobal, meShift
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function p = meBlur(p, amp, sigma)

if sigma <= 0, return; end

r = max(1, ceil(3*sigma));
g = exp(-((-r:r).^2)./(2*sigma^2));
g = g./sum(g);

num = sepBlur(p.*amp, g);
den = sepBlur(amp,    g);

% A neighbourhood with no signal in it at all has no phase to report.  Guarding on
% a fraction of the level's own mean amplitude, rather than on zero, keeps the ratio
% from amplifying a single stray count into a phase.
floorAmp = cast(1e-6*mean(amp,'all'), 'like', p);
p = num./max(den, floorAmp);
end

% =====================================================================
function y = sepBlur(x, g)
%sepBlur  The separable Gaussian, rows then columns, frames untouched.
y = imfilter(x, g.', 'symmetric', 'same');
y = imfilter(y, g  , 'symmetric', 'same');
end
%------------- END OF CODE --------------
