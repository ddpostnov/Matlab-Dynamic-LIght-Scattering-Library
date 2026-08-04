%myographProduct  The one place the myograph _MYO data pair is handled
%
%   A myograph recording - pressure (a video) or wire (a LabChart file) - yields
%   ONE data product per recording, the _MYO PAIR, and it joins the library's
%   ordinary RESULTS / SETTINGS seam:
%
%       Mouse1.avi      ->  Mouse1_MYO_r.mat   Mouse1_MYO_s.mat
%       Rat3.adicht     ->  Rat3_MYO_r.mat     Rat3_MYO_s.mat
%
%   THERE IS NO SOURCE MEMBER, ON PURPOSE.  A myograph recording is not changed by
%   the analysis and re-reading it is fast, so the library describes it and keeps
%   the MEASUREMENT, and copies neither.  What the recording IS lives in
%   results.recording (an identity card, no arrays); what was measured of it lives
%   per window in results.intervals(k).diameter, and for a wire recording in
%   results.channel(i).  The '.avi' or '.adicht' is opened again when a frame is
%   wanted.
%
%   There is no stage flag: every myograph step appends to the SAME pair.  This
%   function is where opening it, writing it and naming it live, so that the three
%   rules below are stated once and cannot drift apart between the steps.
%
%   EXACTLY ONE FUNCTION CREATES THE PAIR - the entry step of the modality
%   (runMyographVideo for the pressure myograph, runLabChart for the wire one).
%   Every other step is strictly load-modify-save: it opens what is there, adds its
%   own results, and writes back.  'open' therefore errors, with one plain sentence,
%   when the pair is not on disk - a step that finds no product has not been
%   given a hard problem to solve, it has been run out of order.
%
%   SAVE REPLACES A FILE WHOLE.  There is no '-append' anywhere in this library: a
%   step that adds a field re-saves the whole variable it changed, and only that
%   one.  Pass [] for the variable a step did not touch and its file is left
%   alone - which is what keeps a settings-only step from rewriting the results.
%
%   SETTINGS NEVER CARRY GRAPHICS.  The workbench hooks that ride in s are closures
%   over a GUI figure; saving one resurrects a ghost window when the file is loaded
%   back.  Every wrapper writes settings.<its own name> = reportSettings(s), and
%   'save' strips the transport hooks again on the way out, so the rule holds even
%   if a caller forgets.
%
%   Syntax:
%      [results,settings] = myographProduct('open',fName)
%      myographProduct('save',fName,results,settings)
%      p        = myographProduct('path',fName,role)
%      ivs      = myographProduct('intervals')
%      chs      = myographProduct('channels')
%
%   INPUTS
%     fName    the recording, named however the caller happens to hold it: the RAW
%              file ('Mouse1.avi', 'Rat3.adicht') or either member of the pair
%              ('Mouse1_MYO_r.mat', '..._s.mat').
%     results  the RESULTS struct, or [] to leave *_MYO_r.mat untouched.
%     settings the SETTINGS struct, or [] to leave *_MYO_s.mat untouched.
%     role     'r' (default) | 's' - which member of the pair to name.
%
%   OUTPUTS
%     results/settings  the two loaded structs ('open').
%     p        the full path of the requested member ('path').
%     ivs      a 0x0 struct array with the interval fields of the result tree (see
%              THE RESULT TREE below for what .frames and .samples mean), the
%              prototype every step grows its results.intervals from ('intervals').
%              Having it here is what keeps the interval a fixed shape: a step that
%              only knows about one branch of it still writes an element the others
%              can be added to.
%     chs      a 0x0 struct array .name .units .fs .data .time .intervals - the
%              prototype of the wire myograph's channel axis ('channels').  It
%              carries the SAMPLES as well as the windows, because a wire recording
%              has no other place for them now that there is no source member -
%              until the windows are chosen, when they move into the windows.
%
%   THE RESULT TREE this product carries (spec claude-docs/myograph-no-source):
%     results.recording          THE RECORDING'S IDENTITY CARD, and no arrays:
%                                .fName .modality .frameRate .fs .nFrames .duration
%                                .measures, plus PRESSURE: .size .rowRange
%                                .pixelSize / WIRE: .channels (the names it read).
%                                There is no whole-recording time base in it, on
%                                purpose: only each window has one, its own
%                                .diameter.time, and nothing is filled with NaN
%                                between the windows.
%     results.timeCrop           [t0 t1] or []
%     results.comments           LabChart comments (wire myograph)
%     results.blocks             LabChart record boundaries [n x 2] (wire myograph)
%     results.intervals(k)       .name .tStart .tEnd .frames .channels
%                                .diameter .samples .propagation .vasomotion
%                                PRESSURE myograph: one vessel, one flat list
%     results.channel(i)         .name .units .fs .data .time .intervals(k)
%                                WIRE myograph: one LabChart file is several
%                                chambers, so the CHANNEL is the outer axis and
%                                results.intervals stays empty.  A window left on
%                                'all channels' appears under every channel and
%                                keeps an empty .channels; one assigned to a chamber
%                                appears under it alone.  .data and .time are the
%                                WHOLE recording only until the windows are chosen;
%                                the intervals step cuts them into .intervals(k)
%                                .samples and empties them - see .samples below.
%     results.meta               formatVersion codeVersion createdTimestamp
%
%   .frames MEANS THE WINDOW'S FRAME RANGE IN THE ORIGINAL RECORDING, [first last],
%   1-based - frames of the '.avi', or samples of the '.adicht'.  It has to: once
%   everything outside the windows is discarded, the original recording is the only
%   thing left for a frame number to be a number of, and finding the window in the
%   video again is exactly what getMyographWallFrame needs.  Empty when the window
%   holds nothing, and empty on a wire recording whose channels do not share one rate
%   - there is no single index that would be true for all of them.
%
%   .samples IS THE WIRE MYOGRAPH'S MEASUREMENT BRANCH, the twin of .diameter: .data
%   and .time for that window on the channel it is stored under, at that channel's
%   own rate.  A pressure window leaves it empty, exactly as a wire window leaves
%   .diameter empty.  It exists because the samples outside the windows are thrown
%   away with everything else, and what is kept has to live where the window is.
%
%   READ THEM THROUGH myographIntervals, which knows about both shapes and hands
%   back one flat list with the channel written onto each element.  Nothing else
%   should have to test which modality it is looking at.
%
%   getProductPath answers the same question for the rest of the library, but only
%   for a name that already carries a role letter.  This one also takes the RAW
%   recording, which is what the entry step is holding before any product exists.
%
% See also: runMyographVideo, runMyographDiameter, runMyographPropagation,
%           runMyographVasomotion, myographIntervals, myographDiameterBranch,
%           reportSettings, wbFileModel, getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function varargout = myographProduct(action,varargin)

switch lower(char(action))
    case 'open'
        [varargout{1:2}]=openProduct(varargin{:});
    case 'save'
        saveProduct(varargin{:});
    case 'path'
        varargout{1}=productPath(varargin{:});
    case 'intervals'
        varargout{1}=emptyIntervals();
    case 'channels'
        varargout{1}=emptyChannels();
    otherwise
        error('myographProduct:badAction', ...
            'Unknown action ''%s'' (expected open / save / path / intervals / channels).', ...
            char(action));
end
end

% =====================================================================
function p=productPath(fName,role)
%productPath  The _MYO_<role>.mat that belongs beside this recording.
%   THE TOKEN LIVES HERE AND NOWHERE ELSE.  Both myograph modalities write the same
%   one on purpose (the interval editor and the vasomotion step are literally the
%   same registry step for each, which needs one input glob), so it is spelled once.
%
%   THE PAIR IS THE WHOLE PRODUCT.  There are two roles because there are two files;
%   'd' is not one of them, and is refused by the same line that refuses any other
%   letter.  The name GRAMMAR is a different thing and is unchanged - wbFileModel
%   still parses 'd' as a legal role letter, because the speckle branch writes one.
if nargin<2 || isempty(role), role='r'; end
role=lower(char(role));
if ~any(strcmp(role,{'r','s'}))
    error('myographProduct:badRole', ...
        ['A myograph product is the pair _MYO_r.mat + _MYO_s.mat, so the role must ' ...
         'be ''r'' or ''s'', not ''%s''.'],role);
end
[fPath,stem]=fileparts(char(fName));
stem=regexprep(stem,'_MYO_[rs]$','');           % either member names the same product
p=fullfile(fPath,[stem '_MYO_' role '.mat']);
end

% =====================================================================
function [results,settings]=openProduct(fName)
%openProduct  Load both members, or say plainly that they are not there.
rName=productPath(fName,'r');
sName=productPath(fName,'s');
if ~isfile(rName) || ~isfile(sName)
    [~,stem]=fileparts(rName);
    stem=regexprep(stem,'_MYO_r$','');
    error('myographProduct:noProduct', ...
        'No myograph data for %s yet - run the Video step on this recording first.',stem);
end
R=load(rName,'results');  results=R.results;
S=load(sName,'settings'); settings=S.settings;
end

% =====================================================================
function saveProduct(fName,results,settings)
%saveProduct  Write the parts that changed, whole; leave the rest alone.
if nargin<2, results=[];  end
if nargin<3, settings=[]; end
if ~isempty(results)
    save(productPath(fName,'r'),'results','-v7.3','-nocompression');
end
if ~isempty(settings)
    settings=stripHooks(settings);
    save(productPath(fName,'s'),'settings','-v7.3','-nocompression');
end
end

% =====================================================================
function settings=stripHooks(settings)
%stripHooks  reportSettings applied to every step's own sub-struct, so the ghost-
%   window rule is enforced HERE as well as in the wrapper that wrote the field.
if ~isstruct(settings), return; end
fn=fieldnames(settings);
for i=1:1:numel(fn)
    if isstruct(settings.(fn{i})) && isscalar(settings.(fn{i}))
        settings.(fn{i})=reportSettings(settings.(fn{i}));
    end
end
end

% =====================================================================
function ivs=emptyIntervals()
%emptyIntervals  The interval prototype: a 0x0 struct array with every branch of
%   the tree, so elements written by different steps concatenate.
ivs=struct('name',{},'tStart',{},'tEnd',{},'frames',{},'channels',{}, ...
    'diameter',{},'samples',{},'propagation',{},'vasomotion',{});
end

% =====================================================================
function chs=emptyChannels()
%emptyChannels  The wire myograph's channel axis, empty.  It is a prototype for the
%   same reason the interval one is: the entry step writes it before any window
%   exists, so every consumer finds the field whether or not the intervals step has
%   been run.
%
%   IT CARRIES THE SAMPLES, UNTIL THE WINDOWS ARE CHOSEN.  A wire recording has no
%   diameter and no frames, so its channels ARE its measurement, and with no source
%   member .data and .time per channel is where they live, each at its own .fs
%   because a rate chosen per channel is a decision about the measurement and not
%   ours to undo.  The intervals step then cuts them into the windows
%   (.intervals(k).samples) and leaves .data and .time empty: what falls outside the
%   analysed windows is discarded, and a whole-recording copy beside the windows
%   would be exactly the copy this product exists not to keep.
chs=struct('name',{},'units',{},'fs',{},'data',{},'time',{},'intervals',{});
end
