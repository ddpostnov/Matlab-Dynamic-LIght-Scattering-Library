%getMetric - Fetch a metrics-table column as a double vector (NaNs if absent)
%
%   v = getMetric(M,name) returns column NAME of the per-segment metrics table
%   M as a double column vector, or a column of NaNs (height(M)x1) when the
%   column is not present.  Pure helper shared by the vascular-tree derivation
%   core (getVascularTree/combinePhi) and the setVascularTree editor GUI.
%
% Syntax:
%    v = getMetric(M,name)
%
% Inputs:
%    M    - table of per-segment metrics (e.g. results.sMetrics).
%    name - char/string column name to fetch.
%
% Outputs:
%    v    - double column vector; NaNs if the column is absent.
%
% See also: getVascularTree, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026

%------------- BEGIN CODE --------------
function v = getMetric(M,name)
% Column as a double vector, or NaNs if the column is absent.
if any(strcmp(name,M.Properties.VariableNames))
    v = double(M.(name));
else
    v = nan(height(M),1);
end
v=v(:);
end
