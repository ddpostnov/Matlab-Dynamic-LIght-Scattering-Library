%getVascularCues - which columns and which picture describe THIS product's vessels
%
%   cues = getVascularCues(fName)
%
% DESCRIPTION
%   THE ONE PLACE THE QUESTION IS ASKED, and the sibling of getVesselPolarity.  The two
%   vascular steps - the artery/vein labelling and the parent/daughter hierarchy - read
%   per-segment columns BY NAME, and the names differ per PRODUCT rather than per
%   modality:
%
%     '_c_BFI'  pulse arrival + pulsatility, from the pulsatility step   psTimeMax
%                                                                       psTimeMin  psPI
%     '_c_I'    the same quantities on a fluorescence beat               pvTimeMax
%                                                                       pvTimeMin  pvPI
%     '_b_I'    TRACER arrival, from the transit-time step               bsDelay  bsMtt
%     '_a_I'    an angiogram is one picture of a whole recording, so it carries no
%               timing at all and NOTHING orders flow on it
%
%   A PREFIX SWAP WOULD NOT HAVE BEEN ENOUGH, and that is why this function names
%   COLUMNS.  'ps' -> 'pv' is a re-prefixing of one family and could have been a two-
%   letter argument; the bolus product's markers are not that family at all - a tracer
%   arrival and a mean transit time have no 'bsPI' and no 'bsTimeMin' - so anything that
%   resolved a PREFIX could not describe the branch this library actually has.
%
%   IT NEVER LOOKS AT THE DATA, only at the name.  What it returns is what the product
%   WOULD carry when the step that writes those columns has run; whether the columns are
%   actually there is the caller's presence test, and saying which of them are missing is
%   how getVesselTypeGuess turns a silent "Uncertain" into a sentence.
%
%   AND IT ANSWERS THE OTHER TWO NAME QUESTIONS the same steps used to answer inline:
%   which RESULTS field holds the mean picture the two editors draw on, and which sibling
%   recordings share this one's segmentation so a derived hierarchy can be copied to them.
%   All three are the same question - what does a product of this name carry - and asking
%   it in three places is how 'imgBFI' ended up hard-coded in a step that now serves a
%   product which has no BFI.
%
% INPUT
%   fName     the RESULTS member's path, i.e. the '*_BFI_r.mat' or '*_I_r.mat' the step
%             was handed.  Reading the name is this function's JOB - it is why the two
%             wrappers and the two cores below it can be told never to do it.
%
% OUTPUT
%   cues      struct describing the product:
%     .product        'BFI' | 'I'      the product token
%     .stage          't'|'s'|'c'|'a'|'b'|''   the stage flag, '' when there is none
%     .branch         'contrast'|'cardiac'|'angiogram'|'bolus'|''
%     .image          RESULTS field holding the mean picture ('imgBFI' | 'imgI')
%     .partners       stage letters this product's hierarchy is INHERITED by - the other
%                     stages of the same recording, which share its segmentation.  Empty
%                     unless this is the stage the hierarchy is derived on (see below)
%     .pulsatility    the pulsatility-index column ('psPI'|'pvPI'|''), '' when the
%                     product carries no beat
%     .arrival        cellstr of the arrival/transit columns, in the order they are
%                     offered ({} when the product carries no timing)
%     .arrivalLabels  a human label per .arrival entry, same size
%     .caliber        cellstr of the U-shaped calibre cues ({'BFI','diameter'} on the
%                     speckle branch, {'diameter'} where there is no flow image)
%     .brightness     per-segment brightness column ('BFI' | '')
%     .brightnessStd  its spread ('std(BFI)' | '')
%     .what           one phrase naming what orders flow here, for a message ('' when
%                     nothing does)
%     .writtenBy      the REGISTRY LABEL of the step that writes those columns, so a
%                     message can tell a user which box to tick ('' when there is none)
%
%   THE HIERARCHY IS DERIVED WHERE THE ARRIVAL COLUMNS ARE WRITTEN and inherited
%   everywhere else, which is what .partners encodes: '_c' on the speckle branch (the
%   pulsatility step runs on the internal cycle) and '_b' on the fluorescence one (the
%   transit-time step runs on the bolus).  A file on any other stage of the same
%   recording gets an empty .partners, because it is a recipient rather than a source.
%
% See also: getVesselPolarity, defaultFlowParams, getVesselTypeGuess, getVascularTree,
%           setVesselTypes, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function cues = getVascularCues(fName)

fName = char(fName);
[stage,product] = parseName(fName);

cues = struct('product',product,'stage',stage,'branch',branchOf(stage), ...
    'image','','partners',{{}},'pulsatility','','arrival',{{}},'arrivalLabels',{{}}, ...
    'caliber',{{}},'brightness','','brightnessStd','','what','','writtenBy','');

switch product
    case 'BFI'
        % THE SPECKLE BRANCH, UNCHANGED.  A BFI product carries a flow image, so the
        % calibre cue is a PAIR - flow and width - and the brightness sliders exist.
        cues.image         = 'imgBFI';
        cues.caliber       = {'BFI','diameter'};
        cues.brightness    = 'BFI';
        cues.brightnessStd = 'std(BFI)';
        cues.pulsatility   = 'psPI';
        cues.arrival       = {'psTimeMax','psTimeMin'};
        cues.arrivalLabels = {'peak arrival','foot arrival'};
        cues.what          = 'when the pulse reaches each vessel and how pulsatile it is';
        cues.writtenBy     = 'Pulsatility';
        if strcmp(stage,'c'), cues.partners = {'t','s'}; end

    case 'I'
        % THE FLUORESCENCE BRANCH.  There is NO BFI twin [D10] - fluorescence intensity
        % is already proportional to plasma volume and there is nothing to invert - so
        % the calibre cue is the width alone and there is no per-segment brightness.
        cues.image   = 'imgI';
        cues.caliber = {'diameter'};
        switch stage
            case 'b'
                % THE TRACER IS THE ARRIVAL MEASUREMENT, and it is a different family of
                % markers rather than a re-spelling of the pulse ones: Delay is a
                % difference of 50 % times against the recording's own arterial input and
                % Mtt a difference of first moments, both larger downstream.  Cth and
                % BvRel are deliberately absent - a width is not monotone along a tree,
                % and BvRel is exposed to the halo until background removal has run.
                cues.arrival       = {'bsDelay','bsMtt'};
                cues.arrivalLabels = {'tracer arrival','mean transit time'};
                cues.what          = 'when the tracer reaches each vessel';
                cues.writtenBy     = 'Transit time';
                cues.partners      = {'a','c'};
            case 'c'
                cues.pulsatility   = 'pvPI';
                cues.arrival       = {'pvTimeMax','pvTimeMin'};
                cues.arrivalLabels = {'peak arrival','foot arrival'};
                cues.what          = 'when the pulse reaches each vessel and how pulsatile it is';
                cues.writtenBy     = 'Pulsatility';
            otherwise
                % AN ANGIOGRAM IS ONE PICTURE OF A WHOLE RECORDING.  It has a vasomotion
                % band and a shape and no timing whatever, so nothing here orders flow -
                % which is a fact about the product and not a gap to be filled in later.
        end

    otherwise
        error('getVascularCues:unknownProduct', ...
            ['"%s" is neither a blood-flow-index nor an intensity result, so there is ' ...
             'no way to tell which of its columns describe flow.'], fName);
end
end

%% =========================  HELPERS (private)  ====================== %%
function [stage,product] = parseName(fName)
%parseName  The stage flag and the product token off the tail '_<stage>_<product>_r.mat'.
%   TWO PATTERNS RATHER THAN ONE WITH AN OPTIONAL GROUP: a named token inside a
%   '(?:...)?' is DROPPED by MATLAB's regexp when the group does not participate, so the
%   stageless spelling is matched separately instead of made optional.
stage = ''; product = '';
t = regexp(fName,'_(?<stage>[a-z])_(?<product>BFI|I)_r\.mat$','names','once');
if ~isempty(t), stage = t.stage; product = t.product; return; end
t = regexp(fName,'_(?<product>BFI|I)_r\.mat$','names','once');
if ~isempty(t), product = t.product; end
end

% =====================================================================
function b = branchOf(stage)
%branchOf  Which analysis pipeline a stage flag names.  The same table wbFileModel
%   keeps, restated here because Core may not read a workbench file.
switch stage
    case {'t','s'}, b = 'contrast';
    case 'c',       b = 'cardiac';
    case 'a',       b = 'angiogram';
    case 'b',       b = 'bolus';
    otherwise,      b = '';
end
end
%------------- END OF CODE --------------
