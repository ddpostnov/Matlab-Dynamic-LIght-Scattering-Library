%wbFileModel - Parse / compose the workbench's stage-flagged file names.
%
%   The Processing Workbench identifies each recording by a NAME MODEL derived
%   from its on-disk file name.  A processed file name has the shape
%
%       [RoiN_]<stem>[_<stage>][_<product>]_<role>.mat
%
%   where the tokens are the library's stage flags (see 01-pipeline-map.md §4):
%       stage   : a SINGLE flag naming the processing stage -
%                   t = temporal contrast    s = spatial contrast
%                   c = internal/cardiac cycle    e = external/epoch cycle
%                   b = bolus
%       product : K (speckle contrast) | I (intensity) | g (g2/DLSI) | BFI |
%                 MYO (myograph - diameter or force channels)
%       role    : d (SOURCE) | r (RESULTS) | s (SETTINGS)
%
%   RATIONALE for the single flag (author decision, 2026-07-28).  A strictly
%   logical name would encode the contrast BASE of a cycle - '_e_t_K'/'_e_s_K'
%   for an external cycle, and then '_c_t_K'/'_c_s_K' for an internal cycle,
%   since a cycle can in principle be built from a temporal OR a spatial contrast
%   base.  But that base is ALSO recorded in the associated SETTINGS file, and it
%   matters little downstream: a project is not expected to mix contrast bases
%   for its cycles.  So the suffix is kept simple, carrying only the flag the
%   next step needs - '_t_K', '_s_K', '_c_K', '_e_K'.  t and s are the
%   interchangeable "contrast" branch and are not expected to coexist in one run
%   by default.
%
%   A myograph product carries NO stage flag, only the product token: both
%   myograph modalities write one '_MYO' PAIR per recording and every later
%   step appends to it, so there is no stage for a name to carry.  It is still on a
%   BRANCH, though - 'myograph', the one its registry steps declare - because the
%   branch is the pipeline a file belongs to and a flagless name does not make it
%   belong to none.
%
%   Raw recordings (.rls/.cxd/.avi/.mraw/.adicht) carry no flags.  A region crop
%   keeps its leading 'RoiN_' prefix as part of the recording identity.
%
%   This function both DECOMPOSES a path into that model and COMPOSES a sibling
%   file name for a given flag chain + role - the workbench needs the latter to
%   derive a step's expected output/input names without a directory scan.
%
%   MODALITY is a GUESS, and the user owns it.  A file's extension narrows what
%   it can possibly be - a .rls is a speckle recording of some sort, an .avi is a
%   video - but not which: the parser picks the most likely member and the Files
%   tab offers the rest in a dropdown.  wbFileModel('modalities') is the whole
%   vocabulary and wbFileModel('modalities',ext) the subset an extension allows:
%
%       LSCI    laser speckle contrast imaging      HSLSCI  high-speed LSCI
%       DLSI    dynamic light scattering imaging    EPFL    fluorescence
%       HYPER   hyperspectral imaging               WMYO    wire myograph
%       PMYO    pressure myograph
%
%   Most are stubs for modalities the library will grow into; the RULE is what
%   matters - the vocabulary lives here, once, and nothing hardcodes it elsewhere.
%   A wire myograph records FORCE CHANNELS, not video: its container is the
%   LabChart '.adicht' file, and a video is guessed to be the PRESSURE myograph.
%
%   THE MODALITY OF A '_MYO' PRODUCT IS NOT LOAD-BEARING.  One product token
%   serves both myograph modalities - the interval editor and the vasomotion step
%   are literally the same registry step for each - so a '_MYO' name cannot say
%   which of the two wrote it and the guess is simply PMYO.  Nothing may depend on
%   that: a workbench row is a RAW RECORDING whose modality the user set on the
%   Files tab, and anything that genuinely needs a product's modality reads
%   results.recording.modality out of the file.
%
% Syntax:
%    model = wbFileModel(path)                         % decompose one path
%    model = wbFileModel(path, rootFolder, resultsFolder)  % ... knowing the two roots
%    name  = wbFileModel('compose', model, chain, role)% build a sibling .mat name
%    p     = wbFileModel('identity', model)            % [RoiN_]stem (no flags)
%    b     = wbFileModel('branch', stage, product)     % the pipeline a file sits on
%    ms    = wbFileModel('modalities')                 % the whole vocabulary
%    ms    = wbFileModel('modalities', ext)            % those an extension allows
%    es    = wbFileModel('extensions')                 % the raw containers it knows
%
% Inputs:
%    path   - char/string full path (or bare name) of a recording or product.
%    rootFolder    - (optional) where the RECORDINGS are.
%    resultsFolder - (optional) where the RESULTS go.  Empty, or equal to
%                    rootFolder, means the two trees are one.
%    model  - a struct returned by the decompose form.
%    chain  - flag chain to compose, e.g. '_t_K', '_c_BFI', '_e_K' or '' .
%    role   - 'd' | 'r' | 's' (SOURCE/RESULTS/SETTINGS).
%    stage  - a single stage flag: 't' | 's' | 'c' | 'e' | 'b' | '' .
%    product- a product token: 'K' | 'BFI' | 'I' | 'g' | 'MYO' | '' .
%
% Outputs:
%    model - struct with fields: path, folder, rawFolder, resultsFolder, name,
%            ext, modality, roi (double or []), roiPrefix, stem, identity, flags
%            (cellstr, name order), stage, branch, product, role, isRaw,
%            isReference, animal, type, index, expGroup.  branch is derived from
%            stage: t|s -> 'contrast', c -> 'cardiac', e -> 'epoch',
%            b -> 'bolus'; a stage-less product falls back to its product token,
%            so 'MYO' -> 'myograph'.
%
%   THREE FOLDERS, AND A NAME ALONE ONLY SETTLES ONE.  '.folder' is where the file
%   itself is, which is all a name says.  '.rawFolder' is where this recording's
%   raw container is and '.resultsFolder' is where its products go - and when a
%   project keeps its results apart from its recordings (getResultsPath) those are
%   two different trees.  Which they are cannot be read off the name, so the caller
%   supplies the two ROOTS and this function maps the path in each direction:
%
%       m = wbFileModel(path, rootFolder, resultsFolder)
%
%   Give it no roots - every caller that predates the results folder, and every
%   project that never names one - and all three folders are the same, which is why
%   nothing downstream branches on whether a results folder was set.  A consumer
%   looking for products reads .resultsFolder, one looking for the recording reads
%   .rawFolder, and both are .folder until they are not.  wbDiscoverFiles forwards
%   the roots it was given, so every model of a working set carries the answer and
%   no consumer re-derives it.
%
%   LABEL AXES.  'animal' is the SUBJECT id (what getFileNamesList calls the
%   animalIdentifier) - the scope of registration / vessel typing / the reference
%   recording.  'type' is the recording's experimental role, which owns the
%   processing configuration.  'index' is the recording index within an animal.
%   'expGroup' is the EXPERIMENTAL group, a comparison label used by Export and
%   Explore only; processing ignores it.  They are INDEPENDENT per-file labels -
%   a group may span animals and an animal may span groups - so none is ever
%   derived from another, and none has a fixed vocabulary.  All four are stamped
%   by the discovery layer (wbDiscoverFiles via wbTypeModel), not by this name
%   parser, which leaves them empty.
%
% Notes:
%    The stage flag is a single token, but a LEGACY external-cycle file written
%    the stacked way ('_t_e_K', which runExternalCycle currently emits) still
%    parses to the same identity: up to two stage tokens are stripped and the
%    reported stage is the one adjacent to the product (here 'e').  A stem that
%    literally ends in a stage letter before a real '_<product>_<role>' is the
%    inherent (rare) ambiguity of the flat naming scheme.
%
% See also: wbStepRegistry, wbStateEngine, wbDiscoverFiles, wbTypeModel,
%           getFileNamesList
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function out = wbFileModel(varargin)

if nargin>=1 && (ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1})))
    action = char(varargin{1});
    switch action
        case 'compose'
            out = composeName(varargin{2:end});
            return
        case 'identity'
            m = varargin{2};
            out = fullfile(m.folder,[m.roiPrefix m.stem]);
            return
        case 'branch'
            % PROMOTED so a caller that already has a (stage, product) pair - the
            % explorer's pipeline filter parses the stage through wbProducts - can
            % name the branch without re-deriving the rule.  ONE definition.
            out = branchOf(varargin{2:end});
            return
        case 'modalities'
            if numel(varargin)>=2, out = modalitiesFor(varargin{2});
            else,                  out = allModalities(); end
            return
        case 'extensions'
            tbl = extTable();  out = reshape(tbl(:,1),1,[]);   % raw containers, in order
            return
    end
end

% default: decompose a single path, optionally knowing the two roots
out = decompose(varargin{:});
end

% =====================================================================
function model = decompose(pth, rootFolder, resultsFolder)
if nargin<2, rootFolder = ''; end
if nargin<3, resultsFolder = ''; end
pth = char(pth);
[folder,name,ext] = fileparts(pth);

roleSet    = {'d','r','s'};
productSet = {'BFI','K','I','g','MYO'};          % MYO carries no stage flag
stageSet   = {'t','s','c','e','b'};              % one stage slot; t|s|c|e|b

% ---- leading region-crop prefix (part of the identity) --------------
roi = [];  roiPrefix = '';
tk = regexp(name,'^Roi(\d+)_','tokens','once');
core = name;
if ~isempty(tk)
    roi = str2double(tk{1});
    roiPrefix = sprintf('Roi%d_',roi);
    core = regexprep(name,'^Roi\d+_','');
end

role=''; product=''; flags={};
isRaw = ~strcmpi(ext,'.mat');

if ~isRaw
    parts = strsplit(core,'_');
    if ~isempty(parts) && ismember(parts{end},roleSet)
        role = parts{end};  parts(end) = [];
    end
    if ~isempty(parts) && ismember(parts{end},productSet)
        product = parts{end};  parts(end) = [];
        % strip the stage flag, which lives only alongside a product.  Canonically
        % ONE flag ('_e_K'), but tolerate a legacy stacked name ('_t_e_K') by
        % taking up to two tokens; the reported stage is the one nearest product.
        for hop = 1:2
            if ~isempty(parts) && ismember(parts{end},stageSet)
                flags = [parts(end) flags];  parts(end) = []; %#ok<AGROW> (<=2 hops)
            else
                break
            end
        end
    end
    stem = strjoin(parts,'_');
else
    stem = core;
end

stage = '';
if ~isempty(flags), stage = flags{end}; end      % flag adjacent to the product

% THE TWO TREES (see the header).  getResultsPath owns the rule and is asked in
% each direction; with no roots, or a project that never named a results folder, it
% returns the path verbatim and all three folders agree.  It is asked for the NAME
% and not for the room: a model is built once per file on every repaint of the Files
% tab, and making a folder here would leave one empty folder per scanned recording,
% cost a disk stat per file per repaint and let an unreachable results folder throw
% out of an edit box.  The wrapper that writes there makes its own room.
rawOf     = fileparts(getResultsPath(pth, rootFolder, resultsFolder, 'root'));
resultsOf = fileparts(getResultsPath(pth, rootFolder, resultsFolder, 'results', 'name'));

model = struct( ...
    'path',        pth, ...
    'folder',      folder, ...
    'rawFolder',   rawOf, ...
    'resultsFolder', resultsOf, ...
    'name',        [name ext], ...
    'ext',         ext, ...
    'modality',    modalityOf(ext,product), ...
    'roi',         roi, ...
    'roiPrefix',   roiPrefix, ...
    'stem',        stem, ...
    'identity',    fullfile(folder,[roiPrefix stem]), ...
    'flags',       {flags}, ...
    'stage',       stage, ...
    'branch',      branchOf(stage,product), ...
    'product',     product, ...
    'role',        role, ...
    'isRaw',       isRaw, ...
    'isReference', false, ...
    'animal',      '', ...
    'type',        '', ...
    'index',       '', ...
    'expGroup',    '');
end

% =====================================================================
function name = composeName(model,chain,role)
%composeName  Build fullfile for a sibling triplet member of this recording.
%   IN THE RESULTS FOLDER, which is where a triplet member is: this composes an
%   OUTPUT name, and .resultsFolder is .folder itself unless a results folder was
%   set (see the header), so the ordinary project reads exactly as it always did.
if nargin<3 || isempty(role), role='d'; end
chain = char(chain);
if ~isempty(chain) && chain(1)~='_', chain = ['_' chain]; end
name = fullfile(model.resultsFolder,[model.roiPrefix model.stem chain '_' role '.mat']);
end

% =====================================================================
function b = branchOf(stage,product)
if nargin<2, product = ''; end
stage = char(stage); product = char(product);
%branchOf  Which analysis pipeline a file sits on - the key everything per-branch is
%   filtered and keyed by.
%
%   THE STAGE FLAG SAYS IT, EXCEPT WHEN THERE IS NO STAGE FLAG.  A myograph
%   recording writes ONE '_MYO' PAIR - _r and _s, and no source member - and every
%   step of the pipeline appends to it, so the name carries a product token and no
%   flag at all.  The name GRAMMAR is unchanged by that: 'd' is still a legal role
%   letter here, because a speckle product still writes one and an old myograph
%   '_MYO_d.mat' still sitting on disk must still parse.  And the file is
%   nonetheless on the myograph pipeline, which is the branch its registry steps
%   declare.  Left empty, nothing keyed by (recording, pipeline) could find it: the
%   run expansion's resume test asks for a file of the row's branch, found none, and
%   so re-ran a finished myograph chain on every Run.
switch stage
    case {'t','s'}, b = 'contrast';   % temporal / spatial contrast (interchangeable)
    case 'c',       b = 'cardiac';    % internal cycle
    case 'e',       b = 'epoch';      % external / NVC cycle
    case 'b',       b = 'bolus';      % bolus (CTTH)
    otherwise
        if strcmp(product,'MYO'), b = 'myograph'; else, b = ''; end
end
end

% =====================================================================
function ms = allModalities()
%allModalities  The modality vocabulary, in menu order.  THE single definition.
ms = {'LSCI','HSLSCI','DLSI','EPFL','HYPER','WMYO','PMYO'};
end

% =====================================================================
function tbl = extTable()
%extTable  THE extension vocabulary, in one place: each raw container, what it
%   can legitimately be, and the best guess among those.
%
%   The container constrains the physics.  A .rls holds raw speckle frames, so it
%   is one of the speckle modalities; a video can be a pressure myograph (a
%   diameter is measured from it) or a wide-field fluorescence/speckle recording;
%   a .adicht is a LabChart file and nothing but a wire myograph writes one.  A
%   processed .mat can have come from any of them and is left unconstrained.
%
%   The GUESS is the first column's most likely member, not a claim: the Files tab
%   offers the whole allowed list in a dropdown and the user's answer wins.
sp = setdiff(allModalities(),{'WMYO','PMYO'},'stable');   % everything imaging
vid = {'PMYO','EPFL','LSCI'};                             % a video: myograph first
tbl = { '.rls',    {'LSCI','HSLSCI','DLSI'}, 'LSCI'
        '.cxd',    {'EPFL','HYPER'},         'EPFL'       % vendor fluorescence stack
        '.mraw',   sp,                       'DLSI'
        '.cihx',   sp,                       'DLSI'
        '.avi',    vid,                      'PMYO'
        '.mp4',    vid,                      'PMYO'
        '.mov',    vid,                      'PMYO'
        '.mkv',    vid,                      'PMYO'
        '.adicht', {'WMYO'},                 'WMYO' };
end

function k = extRow(ext)
%extRow  The table row of one extension (0 when it is not a raw container).
ext = lower(char(ext));
if ~isempty(ext) && ext(1)~='.', ext = ['.' ext]; end
tbl = extTable();
k = find(strcmp(ext, tbl(:,1)),1);
if isempty(k), k = 0; end
end

% =====================================================================
function ms = modalitiesFor(ext)
%modalitiesFor  What a given file extension can legitimately be.
k = extRow(ext);
if k==0, ms = allModalities(); return; end    % .mat products: any origin
tbl = extTable();
ms = tbl{k,2};
end

% =====================================================================
function m = modalityOf(ext,product)
%modalityOf  The BEST GUESS at a file's modality (the user can override it).
%   A '.mat' is read from its PRODUCT token, and 'MYO' guesses the pressure
%   myograph - see the header: that guess is not load-bearing.
if strcmpi(ext,'.mat')
    switch product
        case 'K',   m = 'LSCI';
        case 'I',   m = 'EPFL';
        case 'g',   m = 'DLSI';
        case 'BFI', m = 'LSCI';   % BFI is produced from _K in the LSCI path
        case 'MYO', m = 'PMYO';   % one token, two modalities: results.recording.modality decides
        otherwise,  m = 'LSCI';
    end
    return
end
k = extRow(ext);
if k==0, m = ''; return; end
tbl = extTable();
m = tbl{k,3};
end
