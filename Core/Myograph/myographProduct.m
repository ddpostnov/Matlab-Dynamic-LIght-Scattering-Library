%myographProduct  The one place the myograph _MYO data triplet is handled
%
%   A myograph recording - pressure (a video) or wire (a LabChart file) - yields
%   ONE data product per recording, the _MYO triplet, and it joins the library's
%   ordinary SOURCE / RESULTS / SETTINGS seam:
%
%       Mouse1.avi      ->  Mouse1_MYO_d.mat   Mouse1_MYO_r.mat   Mouse1_MYO_s.mat
%       Rat3.adicht     ->  Rat3_MYO_d.mat     Rat3_MYO_r.mat     Rat3_MYO_s.mat
%
%   There is no stage flag: every myograph step appends to the SAME triplet.  This
%   function is where opening it, writing it and naming it live, so that the three
%   rules below are stated once and cannot drift apart between the steps.
%
%   EXACTLY ONE FUNCTION CREATES THE TRIPLET - the entry step of the modality
%   (runMyographVideo for the pressure myograph, runLabChart for the wire one).
%   Every other step is strictly load-modify-save: it opens what is there, adds its
%   own results, and writes back.  'open' therefore errors, with one plain sentence,
%   when the triplet is not on disk - a step that finds no product has not been
%   given a hard problem to solve, it has been run out of order.
%
%   SAVE REPLACES A FILE WHOLE.  There is no '-append' anywhere in this library: a
%   step that adds a field re-saves the whole variable it changed, and only that
%   one.  Pass [] for the variables a step did not touch and their files are left
%   alone - which is what keeps a per-interval analysis from rewriting the SOURCE
%   cube it only read.
%
%   SETTINGS NEVER CARRY GRAPHICS.  The workbench hooks that ride in s are closures
%   over a GUI figure; saving one resurrects a ghost window when the file is loaded
%   back.  Every wrapper writes settings.<its own name> = reportSettings(s), and
%   'save' strips the transport hooks again on the way out, so the rule holds even
%   if a caller forgets.
%
%   Syntax:
%      [source,results,settings] = myographProduct('open',fName)
%      myographProduct('save',fName,source,results,settings)
%      p        = myographProduct('path',fName,role)
%      ivs      = myographProduct('intervals')
%      chs      = myographProduct('channels')
%
%   INPUTS
%     fName    the recording, named however the caller happens to hold it: the RAW
%              file ('Mouse1.avi', 'Rat3.adicht') or any member of the triplet
%              ('Mouse1_MYO_d.mat', '..._r.mat', '..._s.mat').  All four name the
%              same product.
%     source   the SOURCE struct, or [] to leave *_MYO_d.mat untouched.
%     results  the RESULTS struct, or [] to leave *_MYO_r.mat untouched.
%     settings the SETTINGS struct, or [] to leave *_MYO_s.mat untouched.
%     role     'd' (default) | 'r' | 's' - which member of the triplet to name.
%
%   OUTPUTS
%     source/results/settings  the three loaded structs ('open').
%     p        the full path of the requested triplet member ('path').
%     ivs      a 0x0 struct array with the interval fields of the result tree, the
%              prototype every step grows its results.intervals from ('intervals').
%              Having it here is what keeps the interval a fixed shape: a step that
%              only knows about one branch of it still writes an element the others
%              can be added to.
%     chs      a 0x0 struct array .name .intervals - the prototype of the wire
%              myograph's channel axis ('channels').
%
%   THE RESULT TREE this product carries (spec claude-docs/myograph-workbench):
%     results.timeCrop           [t0 t1] or []
%     results.comments           LabChart comments (wire myograph)
%     results.intervals(k)       .name .tStart .tEnd .frames .channels
%                                .diameter .propagation .vasomotion
%                                PRESSURE myograph: one vessel, one flat list
%     results.channel(i)         .name .intervals(k)   WIRE myograph: one LabChart
%                                file is several chambers, so the CHANNEL is the
%                                outer axis and results.intervals stays empty.  A
%                                window left on 'all channels' appears under every
%                                channel and keeps an empty .channels; one assigned
%                                to a chamber appears under it alone.
%     results.meta               formatVersion codeVersion createdTimestamp
%
%   READ THEM THROUGH myographIntervals, which knows about both shapes and hands
%   back one flat list with the channel written onto each element.  Nothing else
%   should have to test which modality it is looking at.
%
% See also: runMyographVideo, runMyographDiameter, runMyographPropagation,
%           runMyographVasomotion, myographIntervals, reportSettings, wbFileModel
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

function varargout = myographProduct(action,varargin)

switch lower(char(action))
    case 'open'
        [varargout{1:3}]=openProduct(varargin{:});
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
if nargin<2 || isempty(role), role='d'; end
role=lower(char(role));
if ~any(strcmp(role,{'d','r','s'}))
    error('myographProduct:badRole','Role must be ''d'', ''r'' or ''s'', not ''%s''.',role);
end
[fPath,stem]=fileparts(char(fName));
stem=regexprep(stem,'_MYO_[drs]$','');          % a triplet member names the same product
p=fullfile(fPath,[stem '_MYO_' role '.mat']);
end

% =====================================================================
function [source,results,settings]=openProduct(fName)
%openProduct  Load the whole triplet, or say plainly that it is not there.
dName=productPath(fName,'d');
rName=productPath(fName,'r');
sName=productPath(fName,'s');
if ~isfile(dName) || ~isfile(rName) || ~isfile(sName)
    [~,stem]=fileparts(dName);
    stem=regexprep(stem,'_MYO_d$','');
    error('myographProduct:noProduct', ...
        'No myograph data for %s yet - run the Video step on this recording first.',stem);
end
D=load(dName,'source');   source=D.source;
R=load(rName,'results');  results=R.results;
S=load(sName,'settings'); settings=S.settings;
end

% =====================================================================
function saveProduct(fName,source,results,settings)
%saveProduct  Write the parts that changed, whole; leave the rest alone.
if nargin<2, source=[];   end
if nargin<3, results=[];  end
if nargin<4, settings=[]; end
if ~isempty(source)
    save(productPath(fName,'d'),'source','-v7.3','-nocompression');
end
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
    'diameter',{},'propagation',{},'vasomotion',{});
end

% =====================================================================
function chs=emptyChannels()
%emptyChannels  The wire myograph's channel axis, empty.  It is a prototype for the
%   same reason the interval one is: the entry step writes it before any window
%   exists, so every consumer finds the field whether or not the intervals step has
%   been run.
chs=struct('name',{},'intervals',{});
end
