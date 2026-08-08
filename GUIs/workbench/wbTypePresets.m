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
%     Fluorescence angiogram  the '.cxd' angiogram branch: the vasculature, how dense
%                           it is, and how it oscillates slowly
%     Fluorescence cardiac  the '.cxd' cardiac branch: the averaged heartbeat, then
%                           BOTH of its metrics - the pulsatile plasma volume per
%                           vessel and the wall displacement of the same vessels.
%                           It is the only preset that ticks two metrics on one
%                           product, and that is master plan Q1's answer made real:
%                           one step owns the magnification, the other stays a trace
%                           step, and the protocol asks for both
%     Fluorescence bolus    an injection: where the tracer arrives, when, how spread
%                           out the arrivals are, and which vessel feeds which.  It is
%                           the only fluorescence protocol that names an ANIMAL step,
%                           because the vessel types are painted once per animal
%     Pressure myograph     a vessel video: diameter, intervals, propagation and
%                           vasomotion, all on one '_MYO' product
%     Wire myograph         a LabChart recording: read it, choose the windows,
%                           analyse each interval's channels - the same '_MYO'
%                           product and the same two later steps
%
%   perAnimal steps (registration, vessel typing) are listed here too, in
%   .animalSteps, because a protocol is not complete without them - the caller
%   applies those to the animal panel rather than to the type matrix.  A protocol
%   whose modality has none (the myograph) simply lists none.
%
%   A PRESET MAY NAME A STEP THAT IS NOT THERE.  It is applied per step against the
%   type's own modality-filtered registry, and a step that registry does not offer
%   is skipped - so a protocol can be written whole while the steps behind it are
%   still arriving, and a modality that never gains one is never broken by it.
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
% Last revision: 08-August-2026

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
    {'contrast','setRegions','segmentation','dynamicSegmentation', ...
     'BFI','vasomotion','nvc'}, ...
    {'registration','vesselTypes'}, ...
    'stimulus-locked responses on the contrast branch, plus vasomotion');

P(end+1) = mk('Vasoreactivity', ...
    {'contrast','setRegions','segmentation','dynamicSegmentation','BFI'}, ...
    {'registration','vesselTypes'}, ...
    'contrast branch to BFI, no dynamics metric');

P(end+1) = mk('Pulsatility+Vasomotion', ...
    {'contrast','internalCycle','setRegions','segmentation','dynamicSegmentation', ...
     'BFI','vasomotion','pulsatility'}, ...
    {'registration','vesselTypes'}, ...
    'both raw entry steps: the recording drives the contrast AND the cardiac branch');

% THE FLUORESCENCE ANGIOGRAM PROTOCOL.  The vasculature, how much of the field it fills,
% and how it oscillates slowly - which is everything the angiogram branch supports, and
% it is the only fluorescence branch that supports the last of them: the averaged
% heartbeat is a tenth of a second on a phase axis and an injection is half a minute, so
% neither carries a clock a vasomotion band can be measured on.
%
% BACKGROUND REMOVAL IS IN, unlike the cardiac protocol below.  Its whole benefit is a
% clearer picture for the steps that read vessel SHAPE - the segmentation and the
% vascular density, both of which are here - and what has never been measured is what it
% does to a sub-pixel wall position, which is the cardiac protocol's quantity and not
% this one's.  It needs to be told how many micrometres a pixel is and refuses until it
% has been; that is the intended behaviour of a step whose every length is a real
% distance, not a defect in the protocol.
%
% DYNAMIC SEGMENTATION IS IN, and it is the opposite decision from the cardiac protocol
% for the same measured reason.  Its per-frame diameter is quantised at a quarter of a
% pixel, which is far too coarse for a heartbeat and exactly right for vasomotion, where
% a diameter moves whole pixels over tens of seconds - and it is what puts two of
% Vasomotion's four default signals on the product at all.
%
% WALL MOTION IS OUT.  On a continuous recording it has no averaged beat to work from, so
% it measures a band against an empty band of the same width - and on the reference
% recording that came out at 1.0 times its own comparison, i.e. nothing, while looking
% entirely convincing.  It is a tick away for whoever wants to confirm that.
P(end+1) = mk('Fluorescence angiogram', ...
    {'intensity','setRegionsI','backgroundRemoval','segmentationI', ...
     'dynamicSegmentationI','topology','vasomotionI'}, ...
    {}, ...
    'the vasculature of a fluorescence recording: how dense it is and how it oscillates slowly');

% THE FLUORESCENCE CARDIAC PROTOCOL, AND IT TICKS BOTH CARDIAC METRICS ON PURPOSE
% (master plan section 7, Q1, answered 07-Aug-2026).  Pulsatility measures how much the
% labelled plasma in each vessel rises and falls over one beat; wall motion measures how
% far the walls of the same vessels moved, sub-pixel, against the scrambled-timing copies
% the entry step stored beside the cycle.  They are two instruments on one product and
% neither is the other's approximation - which is exactly why ONE step owns the
% magnification rather than the pulsatility step growing a second copy of it.
%
% THREE STEPS ARE DELIBERATELY OUT.  Dynamic segmentation is out because on a cardiac
% product its distinctive output is a per-bin diameter quantised at a quarter of a pixel,
% and the cardiac width change is smaller than that step - the pulsatility step therefore
% writes no diameter column and wall motion is where that quantity is answered.
% Background removal is out because nothing has yet measured what cleaning a picture does
% to a sub-pixel wall position, and this is the one protocol whose numbers are sub-pixel
% wall positions.  Vessel typing is out because it is not a workbench step on this
% modality yet.  Every one of them is a tick away.
P(end+1) = mk('Fluorescence cardiac', ...
    {'intensityCycle','setRegionsI','segmentationI','pulsatilityI','motionEnhancement'}, ...
    {}, ...
    'the averaged heartbeat of a fluorescence recording: how much each vessel pulses, and how far its walls moved');

% THE FLUORESCENCE BOLUS PROTOCOL, and it is the only one that ends in a hierarchy.  An
% injection is the one recording on this modality that carries a tracer arrival, and an
% arrival is what separates an artery from a vein and gives the flow a direction: it
% spreads over seconds across the field, where a heartbeat's phase is confined inside a
% tenth of a second and wraps.
%
% THE LAST TWO STEPS BELONG TO THIS PROTOCOL AND TO NO OTHER.  On an angiogram there is
% no timing at all, so the vessel types would be hand-painted (the step says so out loud
% rather than reporting every vessel undecided) and the hierarchy REFUSES BY NAME - so
% ticking either into the angiogram protocol above would put a step in a protocol that
% cannot do its job.
%
% THE VESSEL TYPES GO IN .animalSteps AND NOT IN .steps, and that is not a stylistic
% choice: the step is perAnimal, the type matrix holds perFile steps only, and a
% perAnimal id listed among .steps is silently dropped by the tick loop - it is not
% stored in a row and it is not reported as an animal step either.  The hierarchy step
% below it is perFile and DOES require it, but a cascade cannot tick an animal step from
% the constructor's preset path, so a protocol that wants vessel typing has to name it
% here.  Vasoreactivity, Pulsatility and the rest name their perAnimal steps the same way.
%
% ONE TICK STILL PULLS IN THE REST OF THE CHAIN.  The hierarchy requires the vessel types
% and the transit time, the transit time requires the segmentation, and the segmentation
% names the three entry steps - so the list below is what the cascade would build anyway,
% written out so the protocol says what it produces rather than implying it.
P(end+1) = mk('Fluorescence bolus', ...
    {'intensityBolus','setRegionsI','backgroundRemoval','segmentationI', ...
     'topology','ctth','vascularTreeI'}, ...
    {'vesselTypesI'}, ...
    'an injection: where the tracer arrives, when, how spread out the arrivals are, and which vessel feeds which');

% The myograph protocols.  A preset is applied step by step and a step the type's
% registry does not offer is simply skipped, so listing the interval editor here
% before it exists costs nothing and saves the list being edited twice: the crop
% and the pre-set intervals stay OUT, because "optional" is what off by default
% means, and the interval editor stays IN, because "enabled by default" is what on
% means.
P(end+1) = mk('Pressure myograph', ...
    {'myoVideo','myoDiameter','myoIntervals','myoPropagation','myoVasomotion'}, ...
    {}, ...
    'a vessel video: diameter, then how it travels along the vessel and how it oscillates');

% The wire myograph is three steps, and that is the whole protocol: reading the
% LabChart file IS the measurement, so there is nothing between it and choosing the
% windows.  Propagation is absent rather than off - it compares locations ALONG one
% vessel, and a wire myograph records one number per chamber.
P(end+1) = mk('Wire myograph', ...
    {'labChart','myoIntervals','myoVasomotion'}, ...
    {}, ...
    'a LabChart recording: the windows to compare, and how each channel oscillates in them');
end

% =====================================================================
function p = mk(name, steps, animalSteps, note)
p = struct('name',name,'steps',{steps},'animalSteps',{animalSteps},'note',note);
end
