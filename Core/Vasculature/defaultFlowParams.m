%defaultFlowParams - Default flow-ordering parameters for the vascular tree
%
%   fp = defaultFlowParams() returns the default struct array of parameters that
%   order the vascular hierarchy in getVascularTree.  Each element has fields
%   name / label / role / scope / weight / enabled:
%     role  = 'arrival'     (increases downstream),
%             'pulsatility' (decreases downstream),
%             'caliber'     (U-shaped: decreases downstream on the arterial
%                            side, increases on the venous side), or
%             'tip'         (geometric: rewards a connection near the vessel tip
%                            OPPOSITE its committed vessel link).
%     scope = 'vessel' / 'parench' for the two 'tip' parameters (else '').
%     weight 0 = off.
%   Arrival/pulsatility drive the global order + the parenchyma links; caliber
%   (BFI/diameter) additionally shapes the per-side vessel trees; tip params do
%   not enter the flow potential.  setVascularTree seeds s.flowParams with this
%   (editable live in the GUI); call it to build s.flowParams for a direct
%   getVascularTree invocation.
%
% Syntax:
%    fp = defaultFlowParams()
%
% Outputs:
%    fp - 1x7 struct array (name/label/role/scope/weight/enabled).
%
% See also: getVascularTree, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026

%------------- BEGIN CODE --------------
function fp = defaultFlowParams()
% Parameters that order the tree.  role = 'arrival' (increases downstream),
% 'pulsatility' (decreases downstream), 'caliber' (U-shaped: decreases
% downstream on the arterial side, increases on the venous side), or 'tip'
% (geometric: rewards a connection near the vessel tip OPPOSITE its committed
% vessel link - scope 'vessel' for the vessel tree, 'parench' for the
% parenchyma links).  weight 0 = off.  Arrival/pulsatility drive the global
% order + the parenchyma links; caliber (BFI/diameter) additionally shapes the
% per-side vessel trees; tip params do not enter the flow potential.
fp = struct( ...
    'name',   {'psTimeMax','psTimeMin','psPI','BFI','diameter','tipVessel','tipParench'}, ...
    'label',  {'peak arrival','foot arrival','pulsatility PI','BFI','diameter', ...
               'tip dist (vessel-vessel)','tip dist (vessel-parench)'}, ...
    'role',   {'arrival','arrival','pulsatility','caliber','caliber','tip','tip'}, ...
    'scope',  {'','','','','','vessel','parench'}, ...
    'weight', {  1,        1,        1,        1,     1,        1,          1 }, ...
    'enabled',{ true,     true,     true,     true,  true,    false,       true} );
end
