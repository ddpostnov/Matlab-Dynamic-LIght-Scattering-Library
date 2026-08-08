%mePhase - Where the structure is, frame by frame: the accumulated quaternionic phase.
%
%   THE QUANTITY THE WHOLE METHOD RESTS ON.  A Laplacian band and its two Riesz
%   components form a quaternion whose phase says where the structure at that scale
%   sits.  Move the wall a hundredth of a pixel and the intensity barely changes,
%   but the phase changes by omega times that displacement - which at level 3
%   (lambda = 16 px) is a fortieth of a radian, a number, not a difference of two
%   pictures that look the same.
%
%   ACCUMULATED FROM DIFFERENCES, NOT READ OFF DIRECTLY, and that is the unwrapping
%   trick from the ICCP 2014 pseudocode.  An absolute phase is only known modulo
%   2*pi and jumps at the wrap; the difference between two consecutive frames of a
%   small motion never approaches pi, so summing the differences gives a phase that
%   runs smoothly over the whole recording.  Frame one is the origin: what comes out
%   is displacement relative to the first frame, which is exactly what a temporal
%   band-pass wants.
%
%   atan2 RATHER THAN THE PSEUDOCODE'S acos.  They give the same angle, but acos
%   near an argument of one loses half its significant digits - and near one is
%   precisely where a sub-pixel motion lives, so the pseudocode's own form is at its
%   worst exactly where this block needs it at its best.  atan2(|imag|, real) is
%   accurate across the whole range and needs no clamp against a rounding error that
%   pushes the argument past unity.
%
%   THE AMPLITUDE IS A WEIGHT, NOT A MEASUREMENT.  It is the geometric mean of the
%   two frames' quaternion magnitudes, and it exists so that meBlur and meGlobal can
%   ignore the places where phase means nothing.  In flat tissue the band is near
%   zero and its phase is the direction of the noise; weighting by amplitude is what
%   keeps that out of the answer.  Frame one has no pair, so it borrows frame two's
%   weight - it contributes no phase, so nothing else is affected.
%
% Syntax:
%    ph = mePhase(pyr)
%    ph = mePhase(pyr, levels)
%
% Inputs:
%    pyr    - the struct from meRieszPyramid.  The levels asked for must carry a
%             Riesz transform.
%    levels - which levels to compute.  Absent or empty takes pyr.levels.
%
% Outputs:
%    ph - .cos    {1 x nLevels} accumulated phase times cos of its orientation,
%                 [rows columns frames]; empty at the levels not asked for
%         .sin    {1 x nLevels} the same times sin of the orientation
%         .amp    {1 x nLevels} the weight, same size
%         .levels the levels that were computed
%
% Example:
%    pyr = meRieszPyramid(stack, 5, 3:5);
%    ph  = mePhase(pyr);
%    d3  = hypot(ph.cos{3}, ph.sin{3});     % radians of motion at lambda = 16 px
%
% Dependencies: none beyond core MATLAB.
% See also: meRieszPyramid, meTemporal, meGlobal, meBlur, meShift
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function ph = mePhase(pyr, levels)

if nargin<2 || isempty(levels), levels = pyr.levels; end
levels = levels(:)';

nL        = numel(pyr.lap);
ph        = struct();
ph.cos    = cell(1,nL);
ph.sin    = cell(1,nL);
ph.amp    = cell(1,nL);
ph.levels = levels;

for k = levels
    if isempty(pyr.rx{k})
        error('mePhase:noRiesz', ...
            'Level %d has no Riesz transform - name it in meRieszPyramid''s third argument.',k);
    end
    A = pyr.lap{k};  B = pyr.rx{k};  C = pyr.ry{k};
    if size(A,3) < 2
        error('mePhase:frames','Phase is a difference between frames; one frame has none.');
    end

    % q(t) * conj(q(t-1)).  The real part is the dot product of the two
    % quaternions and the imaginary part is their wedge; the pseudocode writes it
    % out this way and so does this.
    qR = A(:,:,2:end).*A(:,:,1:end-1) + B(:,:,2:end).*B(:,:,1:end-1) + ...
         C(:,:,2:end).*C(:,:,1:end-1);
    qX = A(:,:,1:end-1).*B(:,:,2:end) - A(:,:,2:end).*B(:,:,1:end-1);
    qY = A(:,:,1:end-1).*C(:,:,2:end) - A(:,:,2:end).*C(:,:,1:end-1);

    den = hypot(qX,qY);
    th  = atan2(den, qR);

    % th/den is the angle per unit of the imaginary part, and it stays finite as
    % den goes to zero (th itself goes to zero with it).  Only an exact zero, where
    % the orientation is genuinely undefined, needs a guard.
    scl = th./den;
    scl(den==0) = 0;

    dCos = scl.*qX;
    dSin = scl.*qY;
    amp  = sqrt(hypot(den,qR));
    clear qR qX qY den th scl

    z        = zeros([size(A,1) size(A,2) 1], 'like', A);
    ph.cos{k} = cat(3, z, cumsum(dCos,3));
    ph.sin{k} = cat(3, z, cumsum(dSin,3));
    ph.amp{k} = cat(3, amp(:,:,1), amp);
end
end
%------------- END OF CODE --------------
