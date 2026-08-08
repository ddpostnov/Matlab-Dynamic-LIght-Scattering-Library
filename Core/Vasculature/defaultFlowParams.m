%defaultFlowParams - Default flow-ordering parameters for the vascular tree
%
%   fp = defaultFlowParams(fName) returns the default struct array of parameters that
%   order the vascular hierarchy in getVascularTree for the product named fName.  Each
%   element has fields name / label / role / scope / weight / enabled:
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
%   IT TAKES THE FILE NAME BECAUSE THE DRIVING COLUMNS ARE A FACT ABOUT THE PRODUCT.
%   getVascularCues resolves them - psTimeMax/psTimeMin/psPI on a speckle internal cycle,
%   bsDelay/bsMtt on a fluorescence bolus, pvTimeMax/pvTimeMin/pvPI on a fluorescence
%   beat - and this function is ONE code path over that answer, not a branch per product.
%   That is deliberate: an arm written out per product is an arm that can go stale
%   unexercised, and the fluorescence branch reaches three of them through one step.
%
%   AN ANGIOGRAM RETURNS NO ARRIVAL AND NO PULSATILITY PARAMETER AT ALL, because it
%   carries no timing.  getVascularTree then refuses by name rather than inventing a
%   direction, which is the loud half of the pair described in getVesselTypeGuess.
%
% Syntax:
%    fp = defaultFlowParams(fName)
%
% Inputs:
%    fName - the RESULTS member's path ('*_c_BFI_r.mat', '*_b_I_r.mat', ...).
%
% Outputs:
%    fp - 1xN struct array (name/label/role/scope/weight/enabled).  Seven entries on a
%         speckle internal cycle, five on a fluorescence bolus (no pulsatility cue and no
%         flow image), six on a fluorescence beat, three on an angiogram.
%
% See also: getVascularCues, getVascularTree, setVascularTree, getVesselTypeGuess
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function fp = defaultFlowParams(fName)
% Parameters that order the tree.  role = 'arrival' (increases downstream),
% 'pulsatility' (decreases downstream), 'caliber' (U-shaped: decreases
% downstream on the arterial side, increases on the venous side), or 'tip'
% (geometric: rewards a connection near the vessel tip OPPOSITE its committed
% vessel link - scope 'vessel' for the vessel tree, 'parench' for the
% parenchyma links).  weight 0 = off.  Arrival/pulsatility drive the global
% order + the parenchyma links; caliber (BFI/diameter) additionally shapes the
% per-side vessel trees; tip params do not enter the flow potential.

c = getVascularCues(fName);

fp = struct('name',{},'label',{},'role',{},'scope',{},'weight',{},'enabled',{});

% arrival first, in the order the product offers them - it is the ordering signal, and
% on the fluorescence bolus branch it is the only one there is
for k = 1:numel(c.arrival)
    fp(end+1) = one(c.arrival{k}, c.arrivalLabels{k}, 'arrival', '', 1, true); %#ok<AGROW>
end
if ~isempty(c.pulsatility)
    fp(end+1) = one(c.pulsatility, 'pulsatility PI', 'pulsatility', '', 1, true);
end
% caliber is U-shaped and never enters the GLOBAL potential - it shapes the per-side
% trees only, which is why a branch with no flow image still enables the width
for k = 1:numel(c.caliber)
    fp(end+1) = one(c.caliber{k}, c.caliber{k}, 'caliber', '', 1, true); %#ok<AGROW>
end
% the two geometric tip parameters are the same on every product: they are about the
% SHAPE of a segment and know nothing about what was measured on it
fp(end+1) = one('tipVessel', 'tip dist (vessel-vessel)',   'tip', 'vessel',  1, false);
fp(end+1) = one('tipParench','tip dist (vessel-parench)',  'tip', 'parench', 1, true);
end

% =====================================================================
function s = one(name,label,role,scope,weight,enabled)
s = struct('name',name,'label',label,'role',role,'scope',scope, ...
           'weight',weight,'enabled',enabled);
end
%------------- END OF CODE --------------
