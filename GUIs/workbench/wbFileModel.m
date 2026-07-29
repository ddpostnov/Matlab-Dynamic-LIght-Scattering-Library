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
%       product : K (speckle contrast) | I (intensity) | g (g2/DLSI) | BFI
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
%   Raw recordings (.rls/.cxd/.avi/.mraw) carry no flags.  A region crop keeps
%   its leading 'RoiN_' prefix as part of the recording identity.
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
%
% Syntax:
%    model = wbFileModel(path)                         % decompose one path
%    name  = wbFileModel('compose', model, chain, role)% build a sibling .mat name
%    p     = wbFileModel('identity', model)            % [RoiN_]stem (no flags)
%    ms    = wbFileModel('modalities')                 % the whole vocabulary
%    ms    = wbFileModel('modalities', ext)            % those an extension allows
%
% Inputs:
%    path   - char/string full path (or bare name) of a recording or product.
%    model  - a struct returned by the decompose form.
%    chain  - flag chain to compose, e.g. '_t_K', '_c_BFI', '_e_K' or '' .
%    role   - 'd' | 'r' | 's' (SOURCE/RESULTS/SETTINGS).
%
% Outputs:
%    model - struct with fields: path, folder, name, ext, modality, roi (double
%            or []), roiPrefix, stem, identity, flags (cellstr, name order),
%            stage, branch, product, role, isRaw, isReference, animal, type,
%            index, expGroup.  branch is derived from stage: t|s -> 'contrast',
%            c -> 'cardiac', e -> 'epoch', b -> 'bolus'.
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
% Last revision: 28-July-2026

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
        case 'modalities'
            if numel(varargin)>=2, out = modalitiesFor(varargin{2});
            else,                  out = allModalities(); end
            return
    end
end

% default: decompose a single path
out = decompose(varargin{1});
end

% =====================================================================
function model = decompose(pth)
pth = char(pth);
[folder,name,ext] = fileparts(pth);

roleSet    = {'d','r','s'};
productSet = {'BFI','K','I','g'};
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

model = struct( ...
    'path',        pth, ...
    'folder',      folder, ...
    'name',        [name ext], ...
    'ext',         ext, ...
    'modality',    modalityOf(ext,product), ...
    'roi',         roi, ...
    'roiPrefix',   roiPrefix, ...
    'stem',        stem, ...
    'identity',    fullfile(folder,[roiPrefix stem]), ...
    'flags',       {flags}, ...
    'stage',       stage, ...
    'branch',      branchOf(stage), ...
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
if nargin<3 || isempty(role), role='d'; end
chain = char(chain);
if ~isempty(chain) && chain(1)~='_', chain = ['_' chain]; end
name = fullfile(model.folder,[model.roiPrefix model.stem chain '_' role '.mat']);
end

% =====================================================================
function b = branchOf(stage)
%branchOf  Map a stage flag to the analysis branch used for column filtering.
switch stage
    case {'t','s'}, b = 'contrast';   % temporal / spatial contrast (interchangeable)
    case 'c',       b = 'cardiac';    % internal cycle
    case 'e',       b = 'epoch';      % external / NVC cycle
    case 'b',       b = 'bolus';      % bolus (CTTH)
    otherwise,      b = '';
end
end

% =====================================================================
function ms = allModalities()
%allModalities  The modality vocabulary, in menu order.  THE single definition.
ms = {'LSCI','HSLSCI','DLSI','EPFL','HYPER','WMYO','PMYO'};
end

% =====================================================================
function ms = modalitiesFor(ext)
%modalitiesFor  What a given file extension can legitimately be.
%   The container constrains the physics: a .rls holds raw speckle frames, so it
%   is one of the speckle modalities; a video can be a myograph or a wide-field
%   fluorescence/speckle recording.  A processed .mat can have come from any of
%   them, so it is left unconstrained.
ext = lower(char(ext));
if ~isempty(ext) && ext(1)~='.', ext = ['.' ext]; end
switch ext
    case '.rls',                        ms = {'LSCI','HSLSCI','DLSI'};
    case {'.mraw','.cihx'},             ms = setdiff(allModalities(),{'WMYO','PMYO'},'stable');
    case {'.avi','.mp4','.mov','.mkv'}, ms = {'WMYO','EPFL','LSCI'};
    case '.cxd',                        ms = {'EPFL','HYPER'};   % vendor fluorescence stack
    otherwise,                          ms = allModalities();    % .mat products: any origin
end
end

% =====================================================================
function m = modalityOf(ext,product)
%modalityOf  The BEST GUESS at a file's modality (the user can override it).
switch lower(ext)
    case '.rls',                     m = 'LSCI';
    case '.cxd',                     m = 'EPFL';
    case {'.avi','.mp4','.mov','.mkv'}, m = 'WMYO';
    case {'.mraw','.cihx'},          m = 'DLSI';
    case '.mat'
        switch product
            case 'K',   m = 'LSCI';
            case 'I',   m = 'EPFL';
            case 'g',   m = 'DLSI';
            case 'BFI', m = 'LSCI';   % BFI is produced from _K in the LSCI path
            otherwise,  m = 'LSCI';
        end
    otherwise,                       m = '';
end
end
