%meCollapse - Put a Laplacian pyramid back together into a picture.
%
%   EXACTLY THE INVERSE OF meRieszPyramid, and that is an assertion rather than a
%   hope.  Each band was formed as G_k minus the expanded coarser level, so adding
%   the expanded coarser level back returns G_k identically, whatever the expand
%   operator does at the frame edge.  The round trip is therefore exact to float
%   rounding, and the test asserts it: if it is not exact, the reconstruction is
%   wrong and every displacement measured downstream is wrong with it, silently and
%   by an amount that looks like biology.
%
%   THE RIESZ COMPONENTS TAKE NO PART.  They are a way of READING the phase of a
%   band, not a part of the band.  meShift rewrites lap and leaves rx and ry alone,
%   and this function reads only lap and res - so a pyramid whose bands have been
%   phase-shifted collapses by the same code as one that has not.
%
% Syntax:
%    I = meCollapse(pyr)
%
% Inputs:
%    pyr - the struct from meRieszPyramid.  Only .lap, .res and .size are read, so
%          a pyramid whose bands meShift has rewritten collapses the same way.
%
% Outputs:
%    I   - [rows columns frames], the class the bands were held in.
%
% Example:
%    pyr = meRieszPyramid(single(frame));
%    err = max(abs(meCollapse(pyr) - single(frame)),[],'all');   % ~1e-6
%
% Dependencies: meExpand.
% See also: meRieszPyramid, meShift, meExpand
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function I = meCollapse(pyr)

I = pyr.res;
for k = numel(pyr.lap):-1:1
    I = pyr.lap{k} + meExpand(I, pyr.size{k});
end
end
%------------- END OF CODE --------------
