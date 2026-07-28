%wbFileModel - Parse / compose the workbench's stage-flagged file names.
%
%   The Processing Workbench identifies each recording by a NAME MODEL derived
%   from its on-disk file name.  A processed file name has the shape
%
%       [RoiN_]<stem>[_<temporal>][_<modifier>][_<product>]_<role>.mat
%
%   where the tokens are the library's stage flags (see 01-pipeline-map.md §4):
%       temporal : t (temporal contrast) | s (spatial) | c (internal/cardiac cycle)
%       modifier : e (external/epoch average) | b (bolus)          [optional]
%       product  : K (speckle contrast) | I (intensity) | g (g2/DLSI) | BFI
%       role     : d (SOURCE) | r (RESULTS) | s (SETTINGS)
%   Raw recordings (.rls/.cxd/.avi/.mraw) carry no flags.  A region crop keeps
%   its leading 'RoiN_' prefix as part of the recording identity.
%
%   This function both DECOMPOSES a path into that model and COMPOSES a sibling
%   file name for a given flag chain + role - the workbench needs the latter to
%   derive a step's expected output/input names without a directory scan.
%
% Syntax:
%    model = wbFileModel(path)                         % decompose one path
%    name  = wbFileModel('compose', model, chain, role)% build a sibling .mat name
%    p     = wbFileModel('identity', model)            % [RoiN_]stem (no flags)
%
% Inputs:
%    path   - char/string full path (or bare name) of a recording or product.
%    model  - a struct returned by the decompose form.
%    chain  - flag chain to compose, e.g. '_t_K', '_c_BFI', '_t_e_K' or '' .
%    role   - 'd' | 'r' | 's' (SOURCE/RESULTS/SETTINGS).
%
% Outputs:
%    model - struct with fields: path, folder, name, ext, modality, roi (double
%            or []), roiPrefix, stem, identity, flags (cellstr, name order),
%            temporal, modifier, product, role, isRaw, isReference, group.
%
% Notes:
%    Flag tokens are only stripped when a product token is present, so a stem
%    that happens to end in a single flag letter is not over-stripped except in
%    the (rare) pathological case of a stem literally ending '..._t' in front of
%    a real '_K_d'.  That ambiguity is inherent to the flat naming scheme.
%
% See also: wbStepRegistry, wbStateEngine, wbDiscoverFiles, getFileNamesList
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
    end
end

% default: decompose a single path
out = decompose(varargin{1});
end

% =====================================================================
function model = decompose(pth)
pth = char(pth);
[folder,name,ext] = fileparts(pth);

roleSet     = {'d','r','s'};
productSet  = {'BFI','K','I','g'};
temporalSet = {'t','s','c'};
modifierSet = {'e','b'};

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
        % flag chain lives only alongside a product; strip a modifier then a
        % temporal token (name order is <temporal>_<modifier>_<product>)
        if ~isempty(parts) && ismember(parts{end},modifierSet)
            flags = [parts(end) flags];  parts(end) = [];
        end
        if ~isempty(parts) && ismember(parts{end},temporalSet)
            flags = [parts(end) flags];  parts(end) = [];
        end
    end
    stem = strjoin(parts,'_');
else
    stem = core;
end

temporal = firstMember(flags,temporalSet);
modifier = firstMember(flags,modifierSet);

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
    'temporal',    temporal, ...
    'modifier',    modifier, ...
    'product',     product, ...
    'role',        role, ...
    'isRaw',       isRaw, ...
    'isReference', false, ...
    'group',       '');
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
function t = firstMember(flags,set)
t = '';
for k=1:numel(flags)
    if ismember(flags{k},set), t = flags{k}; return; end
end
end

% =====================================================================
function m = modalityOf(ext,product)
switch lower(ext)
    case '.rls',                     m = 'LSCI';
    case '.cxd',                     m = 'CXD';
    case {'.avi','.mp4','.mov','.mkv'}, m = 'Myograph';
    case {'.mraw','.cihx'},          m = 'DLSI';
    case '.mat'
        switch product
            case 'K',   m = 'LSCI';
            case 'I',   m = 'CXD';
            case 'g',   m = 'DLSI';
            case 'BFI', m = 'LSCI';   % BFI is produced from _K in the LSCI path
            otherwise,  m = 'LSCI';
        end
    otherwise,                       m = '';
end
end
