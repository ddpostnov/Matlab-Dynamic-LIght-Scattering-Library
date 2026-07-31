%wbTypePresets - Named protocol pipelines a recording TYPE can be built from.
%
%   The Constructor lets a type copy its configuration from another type; these
%   are the same thing for a user who has no other type yet - the standard
%   protocols this library exists to run, each as the list of steps it ticks.
%   Picking one POPULATES the checkboxes and nothing more: every box stays
%   editable afterwards under the ordinary Constructor rules (prerequisites are
%   pulled in, dependants pushed out), so a preset is a starting point, never a
%   mode the type is locked into.
%
%   The lists are the author's protocols, transcribed from the launchers:
%
%     Pulsatility           the cardiac branch: internal cycle -> ... -> pulsatility
%     Vasomotion            the contrast branch: contrast -> ... -> vasomotion
%     NVC                   vasomotion plus the external (stimulus) cycle
%     Vasoreactivity        the contrast branch without a dynamics metric
%     Pulsatility+Vasomotion  BOTH raw entry steps on one recording - the type that
%                           drives two branches at once (spec §6 branchScope 'one':
%                           vasomotion resolves to _t, pulsatility to _c, with no
%                           ambiguity to expose in the UI)
%
%   perAnimal steps (registration, vessel typing) are listed here too, in
%   .animalSteps, because a protocol is not complete without them - the caller
%   applies those to the animal panel rather than to the type matrix.
%
%   Pure data, no graphics: the Constructor reads it, the tests read it, and
%   adding a protocol is a data edit here.
%
% Syntax:
%    names = wbTypePresets('names')
%    p     = wbTypePresets('get', name)
%    tf    = wbTypePresets('exists', name)
%
% Outputs:
%    names - cellstr of the preset names, in menu order.
%    p     - struct with fields:
%              name        char, as listed;
%              steps       cellstr of per-file step ids to tick;
%              animalSteps cellstr of per-animal step ids to tick;
%              note        one line describing what the protocol produces.
%
% See also: wbTypeSelection, wbStepRegistry, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function out = wbTypePresets(action, name)

P = presetTable();
switch action
    case 'names',  out = {P.name};
    case 'exists', out = any(strcmp(char(name), {P.name}));
    case 'get'
        k = find(strcmp(char(name), {P.name}), 1);
        if isempty(k), out = []; else, out = P(k); end
    otherwise
        error('wbTypePresets:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function P = presetTable()
%presetTable  THE definition.  Step ids are wbStepRegistry ids, in pipeline order.
P = struct('name',{},'steps',{},'animalSteps',{},'note',{});

P(end+1) = mk('Pulsatility', ...
    {'internalCycle','setRegions','segmentation','BFI','pulsatility'}, ...
    {'registration','vesselTypes'}, ...
    'cardiac branch: per-segment pulsatility from the internal cycle');

P(end+1) = mk('Vasomotion', ...
    {'contrast','setRegions','segmentation','dynamicSegmentation','BFI','vasomotion'}, ...
    {'registration','vesselTypes'}, ...
    'contrast branch: slow-dynamics vasomotion metrics');

P(end+1) = mk('NVC', ...
    {'contrast','externalCycle','setRegions','segmentation','dynamicSegmentation', ...
     'BFI','vasomotion'}, ...
    {'registration','vesselTypes'}, ...
    'stimulus-locked epochs on the contrast branch, plus vasomotion');

P(end+1) = mk('Vasoreactivity', ...
    {'contrast','setRegions','segmentation','dynamicSegmentation','BFI'}, ...
    {'registration','vesselTypes'}, ...
    'contrast branch to BFI, no dynamics metric');

P(end+1) = mk('Pulsatility+Vasomotion', ...
    {'contrast','internalCycle','setRegions','segmentation','dynamicSegmentation', ...
     'BFI','vasomotion','pulsatility'}, ...
    {'registration','vesselTypes'}, ...
    'both raw entry steps: the recording drives the contrast AND the cardiac branch');
end

% =====================================================================
function p = mk(name, steps, animalSteps, note)
p = struct('name',name,'steps',{steps},'animalSteps',{animalSteps},'note',note);
end
