%quasiUniform - Roberts' R_d low-discrepancy sequence on the unit cube.
%
%   quasiUniform(n,dim) returns n points spread over the dim-dimensional unit cube
%   by an additive recurrence with the generalised golden ratio.  It exists so that
%   a multi-start fit can be SEEDED DETERMINISTICALLY: every fitting core in this
%   library draws its start grid from here, so the same trace always gives the same
%   answer and a fit can be pinned by a test.
%
%   WHY NOT rand.  Two reasons, and both of them have cost time elsewhere.  A random
%   start set makes a fit irreproducible between runs, which means a golden can never
%   be written for it; and it leaks rng state, so a trace fitted on a parfor worker
%   and the same trace fitted in the client can take different paths.  Roberts' R_d
%   sequence also covers a box far more evenly than the same number of random draws,
%   which is the property the multi-start actually needs.
%
%   phi is the real root of x^(dim+1) = x + 1 (the generalised golden ratio), found
%   by the fixed-point iteration x <- (1+x)^(1/(dim+1)); it converges in a handful of
%   steps for every dimension a fit in this library uses.  The 0.5 offset centres the
%   first point rather than starting it on the box edge.
%
% Syntax:
%    U = quasiUniform(n, dim)
%
% Inputs:
%    n    - how many points to draw.
%    dim  - how many dimensions (one per fitted parameter drawn from the grid).
%
% Outputs:
%    U    - [n x dim] points in [0,1).  The caller maps each column onto its own
%           parameter range, linearly or logarithmically as that parameter needs.
%
% Example:
%    U = quasiUniform(12,2);                       % 12 starts over 2 shape parameters
%    tau = exp(log(0.2)+U(:,1)*(log(8)-log(0.2))); % mapped log-uniformly onto [0.2 8]
%
% See also: fitVasoreactivity
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

%------------- BEGIN CODE --------------
function U=quasiUniform(n,dim)

phi=2;
for k=1:32, phi=(1+phi)^(1/(dim+1)); end
alpha=phi.^-(1:dim);
U=mod(0.5+(1:n)'*alpha,1);
end
