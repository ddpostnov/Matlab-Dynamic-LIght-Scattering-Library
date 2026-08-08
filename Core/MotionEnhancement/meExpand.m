%meExpand - Upsample one pyramid band to the size of the level below it.
%
%   THE ONE EXPAND OPERATOR OF THE BLOCK, and it is public rather than a
%   subfunction for a single reason: meRieszPyramid uses it to FORM the Laplacian
%   bands and meCollapse uses it to PUT THEM BACK.  A Laplacian pyramid is exactly
%   invertible only because those two operations are the same operation - write it
%   twice and the two copies are one edit away from disagreeing, at which point the
%   round trip stops being exact and every displacement measured downstream is
%   wrong by an amount nobody can see.
%
%   ZEROS ARE INSERTED, THEN BLURRED WITH TWICE THE REDUCE KERNEL.  Half the
%   samples of the upsampled grid are zero, so the binomial [1 4 6 4 1]/16 that
%   reduces would return half the brightness on the way back; [1 4 6 4 1]/8 is the
%   same filter scaled to put it back.  Both the even and the odd output samples
%   collect exactly 8/16 of the kernel, which is why one factor serves both.
%
%   THE BOUNDARY RULE DOES NOT HAVE TO BE RIGHT, IT HAS TO BE THE SAME.  Whatever
%   this does at the frame edge, meRieszPyramid subtracts it and meCollapse adds it
%   back, so the round trip is exact for any consistent choice.  Symmetric padding
%   is the choice because it is the one that does not invent an edge.
%
% Syntax:
%    U = meExpand(L, sz)
%
% Inputs:
%    L  - a pyramid band, [rows columns] or [rows columns frames].
%    sz - [rows columns] of the level being expanded to.  It must be the size the
%         reduce step came from: numel(1:2:sz(1)) by numel(1:2:sz(2)) has to be
%         size(L,1) by size(L,2).
%
% Outputs:
%    U  - the band at sz, same class and same number of frames as L.
%
% Example:
%    U = meExpand(coarse, size(fine,[1 2]));
%    lap = fine - U;
%
% Dependencies: Image Processing Toolbox (imfilter).
% See also: meRieszPyramid, meCollapse
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function U = meExpand(L, sz)

nr = numel(1:2:sz(1));
nc = numel(1:2:sz(2));
if size(L,1)~=nr || size(L,2)~=nc
    error('meExpand:size', ...
        'A band of %dx%d cannot expand to %dx%d - that size came from a different level.', ...
        size(L,1), size(L,2), sz(1), sz(2));
end

% imfilter takes its kernel in double whatever the image is in, and returns the
% image's class.
k = [1 4 6 4 1]./8;

U = zeros([sz(1) sz(2) size(L,3)], 'like', L);
U(1:2:sz(1), 1:2:sz(2), :) = L;
U = imfilter(U, k.', 'symmetric', 'same');
U = imfilter(U, k  , 'symmetric', 'same');
end
%------------- END OF CODE --------------
