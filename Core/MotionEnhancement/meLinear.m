%meLinear - Linear Eulerian magnification.  THE CONTROL, NOT A FALLBACK.
%
%   THIS IS HERE TO BE BEATEN.  Wu et al., SIGGRAPH 2012 amplify the temporally
%   band-passed INTENSITY of each Laplacian band.  The method is older, simpler and
%   easier to implement than the phase-based one, and on shot-noise-limited
%   fluorescence it is the wrong choice for one reason, stated in Wadhwa et al.,
%   SIGGRAPH 2013 Table 1 and audible in one word: linear magnification AMPLIFIES
%   noise, phase-based magnification TRANSLATES it.  Multiply a band by alpha and
%   its photon noise is multiplied by alpha too.
%
%   NOBODY SHOULD EVER PROMOTE THIS FUNCTION.  If the phase-based path is
%   disappointing the answer is more averaging before the magnifier, a coarser
%   level, or a learned magnifier - never this.  It exists because the VISIBLE
%   DIFFERENCE between its output and the phase-based output on the same recording
%   is evidence about where the noise floor sits, and because it is the only honest
%   way to show a reader why the harder method was worth writing.
%
%   IT ALSO SEPARATES INTENSITY FROM MOTION.  The lumen of a FITC-filled vessel
%   brightens and dims at the cardiac frequency as haematocrit and flow change.  That
%   is a pure intensity change with no motion in it at all.  This function will
%   render it as motion; the phase-based one is largely blind to it.  Run both and
%   the difference names the confound.
%
%   THE BOUND IS DIFFERENT AND TIGHTER: (1+alpha)*delta < lambda/8, against
%   alpha*delta < lambda/4 for the phase-based path.  No attenuation of the finest
%   levels is applied here; s.levels chooses the bands and everything else passes
%   through, which keeps the comparison against meShift like for like.
%
% Syntax:
%    out = meLinear(stack, fs, s)
%
% Inputs:
%    stack - [rows columns frames].
%    fs    - frames per second.
%    s     - settings.  Fields read: .levels, .alpha, .nPyrLevels,
%            .temporalFilter, .passband, .filterOrder.
%
% Outputs:
%    out   - the magnified stack, single, the size of the input.
%
% Example:
%    s = meSettings;  s.passband = [7.9 14.7];  s.alpha = 20;
%    ctrl = meLinear(stack, 100, s);        % beside meShift(stack,100,s)
%
% Dependencies: meRieszPyramid, meTemporal, meCollapse.
% See also: meShift, meValidate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function out = meLinear(stack, fs, s)

if ~isa(stack,'single') && ~isa(stack,'gpuArray'), stack = single(stack); end

nLev = s.nPyrLevels;
if isempty(nLev), nLev = max(s.levels); end

% No Riesz transform is needed: this path never looks at phase.
pyr    = meRieszPyramid(stack, nLev, []);
levels = intersect(s.levels(:)', 1:nLev);

alpha = s.alpha(:)';
if isscalar(alpha)
    alpha = repmat(alpha, 1, numel(levels));
elseif numel(alpha) ~= numel(levels)
    error('meLinear:alpha', ...
        'alpha is one number or one per amplified level; %d given for %d levels.', ...
        numel(alpha), numel(levels));
end

for i = 1:numel(levels)
    k = levels(i);
    pyr.lap{k} = pyr.lap{k} + alpha(i).*meTemporal(pyr.lap{k}, fs, s);
end
out = meCollapse(pyr);
end
%------------- END OF CODE --------------
