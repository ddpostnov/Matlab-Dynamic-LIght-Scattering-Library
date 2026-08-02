%wbSettingsModel - Settings bag keyed by STEP and by recording TYPE, with presets.
%
%   Holds the workbench's parameters and resolves the concrete settings struct s
%   for a given step - and, since the Constructor tab, for a given recording TYPE,
%   because processing choices are a property of the type, not of the file (spec
%   D4: type configuration is GLOBAL, no per-animal and no per-experimental-group
%   override; a divergent animal is a second type).
%
%   Layering, applied in order - each layer shadows the ones above it:
%     1. the step's own defaults            (step.presets.default, from the launchers)
%     2. shared-key bag                      (a shared key edited once propagates to
%                                             every step that reads it)
%     3. per-step override                   (shadows the bag for one step)
%     4. per-TYPE shared-key bag             (2, scoped to one type)
%     5. per-(TYPE,step) override            (3, scoped to one type)
%     6. per-file override                   (legacy; the Constructor writes none)
%   The bag starts empty, so before the user touches a shared key each step keeps
%   its own launcher default (the launchers legitimately diverge, e.g. a
%   different trustLimitsK per step); the first edit of a shared key is what makes
%   it propagate.  Layers 4-5 are what make one type's edits invisible to another:
%   resolving without a type yields exactly the pre-Constructor answer, and the
%   type layers only ever ADD on top of it.
%
%   KEYS.  Per-step overrides are keyed by STEP ID and the type layer by
%   'type||stepId' - never by field name alone.  Several steps read a field of the
%   same name ('method' lives on both the internal cycle and BFI, trustLimitsK on
%   three steps), so a field-keyed store would silently couple them and their
%   staleness cascades.
%
%   THE TYPE VOCABULARY IS OPEN: a type is any token the user's regexp matched or
%   typed, it becomes a key the first time it is edited, and nothing here
%   enumerates or validates it.  A type's configuration is NOT dropped when no
%   file currently carries that type - re-typing a file restores it.
%
% Syntax:
%    model = wbSettingsModel('new')
%    model = wbSettingsModel('new', initBag)          % seed the shared bag
%    s     = wbSettingsModel('resolve', model, step)
%    s     = wbSettingsModel('resolve', model, step, fileModel)
%    s     = wbSettingsModel('resolve', model, step, type)        % char = a TYPE
%    s     = wbSettingsModel('resolve', model, step, fileModel, type)
%    model = wbSettingsModel('setShared',     model, key, value)
%    model = wbSettingsModel('setStep',       model, stepId, field, value)
%    model = wbSettingsModel('setFile',       model, fileKey, stepId, field, value)
%    model = wbSettingsModel('setTypeShared', model, type, key, value)
%    model = wbSettingsModel('setTypeStep',   model, type, stepId, field, value)
%    model = wbSettingsModel('copyType',      model, srcType, dstType)
%    model = wbSettingsModel('resetType',     model, type)
%    model = wbSettingsModel('renameStepField', model, stepId, oldField, newField)
%    types = wbSettingsModel('types', model)
%            wbSettingsModel('save',  model, path)
%    model = wbSettingsModel('load',  path)
%
% Inputs:
%    step      - a step struct from wbStepRegistry (carries presets/sharedKeys).
%    fileModel - a wbFileModel struct (its .identity keys per-file overrides).
%    type      - a recording type token (char).  Empty = the global answer.
%
% Outputs:
%    model - struct with fields: bag (struct), stepOverrides (struct of structs),
%            fileOverrides (containers.Map key 'identity||stepId' -> struct),
%            typeBag       (containers.Map type -> struct, the per-type shared keys),
%            typeOverrides (containers.Map 'type||stepId' -> struct).
%    s     - the resolved settings struct for a run.
%
% See also: wbStepRegistry, wbStateEngine, wbFileModel, wbTypeSelection,
%           wbInvalidate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbSettingsModel(action, varargin)

switch action
    case 'new',           out = newModel(varargin{:});
    case 'resolve',       out = resolve(varargin{:});
    case 'setShared',     out = setShared(varargin{:});
    case 'setStep',       out = setStep(varargin{:});
    case 'setFile',       out = setFile(varargin{:});
    case 'setTypeShared', out = setTypeShared(varargin{:});
    case 'setTypeStep',   out = setTypeStep(varargin{:});
    case 'copyType',      out = copyType(varargin{:});
    case 'resetType',     out = resetType(varargin{:});
    case 'renameStepField', out = renameStepField(varargin{:});
    case 'types',         out = configuredTypes(varargin{:});
    case 'save',          saveModel(varargin{:});  out = [];
    case 'load',          out = loadModel(varargin{:});
    otherwise, error('wbSettingsModel:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function model = newModel(initBag)
if nargin<1 || isempty(initBag), initBag = struct(); end
model = struct( ...
    'bag',           initBag, ...
    'stepOverrides', struct(), ...
    'fileOverrides', containers.Map('KeyType','char','ValueType','any'), ...
    'typeBag',       containers.Map('KeyType','char','ValueType','any'), ...
    'typeOverrides', containers.Map('KeyType','char','ValueType','any'));
end

% =====================================================================
function s = resolve(model, step, varargin)
%resolve  Layer the six sources for one step, optionally scoped to a file/type.
[fileModel, type] = parseScope(varargin);
model = fillModel(model);

% 1) step defaults
if isfield(step.presets,'default') && isstruct(step.presets.default)
    s = step.presets.default;
else
    s = struct();
end
% 2) shared-key bag (only the keys this step reads)
s = applyShared(s, model.bag, step);
% 3) per-step override
if isfield(model.stepOverrides, step.id)
    s = overlay(s, model.stepOverrides.(step.id));
end
% 4) per-type shared-key bag, 5) per-(type,step) override
if ~isempty(type)
    if isKey(model.typeBag, type)
        s = applyShared(s, model.typeBag(type), step);
    end
    tk = typeKey(type, step.id);
    if isKey(model.typeOverrides, tk)
        s = overlay(s, model.typeOverrides(tk));
    end
end
% 6) per-file override (legacy layer; the Constructor writes none)
if ~isempty(fileModel)
    key = [fileModel.identity '||' step.id];
    if isKey(model.fileOverrides, key)
        s = overlay(s, model.fileOverrides(key));
    end
end
end

% =====================================================================
function [fileModel, type] = parseScope(args)
%parseScope  A char scope argument is a TYPE, a struct one is a file model, so
%   resolve(model,step,fileModel) keeps its old meaning and resolve(model,step,type)
%   reads naturally.  Either order, both optional.
fileModel = []; type = '';
for i = 1:numel(args)
    a = args{i};
    if isempty(a), continue; end
    if ischar(a) || (isstring(a) && isscalar(a)), type = char(a);
    elseif isstruct(a),                          fileModel = a;
    end
end
end

% =====================================================================
function s = applyShared(s, bag, step)
%applyShared  Copy the shared keys this step reads out of one bag.
if ~isstruct(bag), return; end
shared = step.sharedKeys;
for i = 1:numel(shared)
    k = shared{i};
    if isfield(bag, k), s.(k) = bag.(k); end
end
% libraryFolder is shared for every step
if isfield(bag,'libraryFolder'), s.libraryFolder = bag.libraryFolder; end
end

% =====================================================================
function model = setShared(model, key, value)
model = fillModel(model);
model.bag.(key) = value;
end

% =====================================================================
function model = setStep(model, stepId, field, value)
model = fillModel(model);
if ~isfield(model.stepOverrides, stepId)
    model.stepOverrides.(stepId) = struct();
end
model.stepOverrides.(stepId).(field) = value;
end

% =====================================================================
function model = setFile(model, fileKey, stepId, field, value)
model = fillModel(model);
key = [fileKey '||' stepId];
if isKey(model.fileOverrides, key), ov = model.fileOverrides(key); else, ov = struct(); end
ov.(field) = value;
model.fileOverrides(key) = ov;
end

% =====================================================================
function model = setTypeShared(model, type, key, value)
%setTypeShared  A shared key edited for ONE type: it propagates to every step of
%   that type which reads it, and to no other type.
model = fillModel(model);
type = char(type);
if isKey(model.typeBag, type), bag = model.typeBag(type); else, bag = struct(); end
bag.(key) = value;
model.typeBag(type) = bag;
end

% =====================================================================
function model = setTypeStep(model, type, stepId, field, value)
model = fillModel(model);
k = typeKey(type, stepId);
if isKey(model.typeOverrides, k), ov = model.typeOverrides(k); else, ov = struct(); end
ov.(field) = value;
model.typeOverrides(k) = ov;
end

% =====================================================================
function model = copyType(model, src, dst)
%copyType  Build one type's configuration from another's ("BP starts as a BV").
%   dst's own type layers are REPLACED, not merged - the point is to make the two
%   identical and then diverge by hand.  The global layers are untouched.
model = fillModel(model);
src = char(src); dst = char(dst);
if strcmp(src,dst), return; end
model = resetType(model, dst);
if isKey(model.typeBag, src), model.typeBag(dst) = model.typeBag(src); end
k = keys(model.typeOverrides);
for i = 1:numel(k)
    [t, stepId] = splitTypeKey(k{i});
    if strcmp(t, src), model.typeOverrides(typeKey(dst,stepId)) = model.typeOverrides(k{i}); end
end
end

% =====================================================================
function model = resetType(model, type)
%resetType  Drop a type's own layers - it falls back to the launcher defaults
%   (plus whatever the global bag/step layers say), which is the "start over"
%   button next to copyType.
model = fillModel(model);
type = char(type);
if isKey(model.typeBag, type), remove(model.typeBag, type); end
k = keys(model.typeOverrides);
for i = 1:numel(k)
    t = splitTypeKey(k{i});
    if strcmp(t, type), remove(model.typeOverrides, k{i}); end
end
end

% =====================================================================
function model = renameStepField(model, stepId, oldField, newField)
%renameStepField  Move a field from oldField to newField in every layer scoped to ONE
%   step - the per-step override, the per-(type,step) override and the legacy per-file
%   override.  This is what a stored settings layer needs when a WRAPPER renames a
%   parameter: the layers are keyed by step id and field name, so an override recorded
%   against the old name would otherwise keep being applied, under a name that now
%   means something else.
%
%   THE SHARED BAGS ARE DELIBERATELY NOT TOUCHED.  A shared key is shared BY NAME with
%   other steps (see sharedKeys in wbStepRegistry), so renaming it in the bag would
%   rewrite their settings too - and a key that is shared is, by definition, one whose
%   meaning did not diverge.  Only the step-scoped layers are this step's alone.
%
%   Warns when something actually moved, because a resumed session quietly changing
%   which parameter a number lands on is the one thing a migration must not do silently.
model = fillModel(model);
n = 0;
if isfield(model.stepOverrides, stepId)
    [model.stepOverrides.(stepId), moved] = moveField(model.stepOverrides.(stepId), oldField, newField);
    n = n + moved;
end
n = n + renameInMap(model.typeOverrides, stepId, oldField, newField);
n = n + renameInMap(model.fileOverrides, stepId, oldField, newField);
if n > 0
    warning('wbSettingsModel:fieldRenamed', ...
        ['%d saved setting(s) of the %s step were recorded under the old name "%s" ' ...
        'and have been moved to "%s".'], n, stepId, oldField, newField);
end
end

% =====================================================================
function n = renameInMap(m, stepId, oldField, newField)
%renameInMap  The rename over one Map whose keys END in '||stepId'.
%   Both Maps use that shape ('type||stepId' and 'identity||stepId'), and matching the
%   SUFFIX rather than splitting means nothing depends on what the first segment is or
%   how many separators it contains.  containers.Map is a handle, so the caller's Map is
%   updated in place and only the count comes back.
n = 0;
k = keys(m);
tail = ['||' stepId];
for i = 1:numel(k)
    if ~endsWith(k{i}, tail), continue; end
    [ov, moved] = moveField(m(k{i}), oldField, newField);
    m(k{i}) = ov;
    n = n + moved;
end
end

% =====================================================================
function [s, moved] = moveField(s, oldField, newField)
%moveField  oldField -> newField in one override struct; 0/1 for whether it was there.
%   A layer already carrying the NEW name was written after the rename, so its value
%   wins and the stale one is simply dropped.
moved = 0;
if ~isstruct(s) || ~isfield(s, oldField), return; end
if ~isfield(s, newField)
    s.(newField) = s.(oldField);
    moved = 1;
end
s = rmfield(s, oldField);
end

% =====================================================================
function t = configuredTypes(model)
%configuredTypes  Every type carrying its own settings (order of the Map keys).
model = fillModel(model);
t = keys(model.typeBag);
k = keys(model.typeOverrides);
for i = 1:numel(k), t{end+1} = splitTypeKey(k{i}); end %#ok<AGROW>
t = unique(t,'stable');
end

% =====================================================================
function k = typeKey(type, stepId), k = [char(type) '||' char(stepId)]; end
function [t, stepId] = splitTypeKey(k)
p = strsplit(k,'||');
t = p{1}; stepId = '';
if numel(p)>1, stepId = strjoin(p(2:end),'||'); end
end

% =====================================================================
function model = fillModel(model)
%fillModel  Add any layer an older model (a Phase-3 preset) predates.
if ~isfield(model,'typeBag') || ~isa(model.typeBag,'containers.Map')
    model.typeBag = containers.Map('KeyType','char','ValueType','any');
end
if ~isfield(model,'typeOverrides') || ~isa(model.typeOverrides,'containers.Map')
    model.typeOverrides = containers.Map('KeyType','char','ValueType','any');
end
end

% =====================================================================
function v = presetVersion(), v = 2; end

% =====================================================================
function saveModel(model, pth)
%saveModel  Serialise the model to a .mat (containers.Map -> plain struct).
%   wbPreset.schema is the PRESET layout version - 1 (or absent, which is how every
%   preset written before this field reads) is the layout in which runInternalCycle
%   still called its reference-trace limits trustLimitsK/I; 2 is after that rename.
%   A preset carries the same step-scoped override layers a session does, so it needs
%   the same migration on load and therefore the same way of dating itself.
model = fillModel(model);
wbPreset = model;
wbPreset.schema = presetVersion();
% store the Map layers as parallel cell arrays so the preset is Map-free on disk
wbPreset.fileOverrides = [];
wbPreset.foKeys = keys(model.fileOverrides);
wbPreset.foVals = values(model.fileOverrides);
wbPreset.typeBag = [];
wbPreset.tbKeys = keys(model.typeBag);
wbPreset.tbVals = values(model.typeBag);
wbPreset.typeOverrides = [];
wbPreset.toKeys = keys(model.typeOverrides);
wbPreset.toVals = values(model.typeOverrides);
save(pth,'wbPreset','-v7.3','-nocompression');
end

% =====================================================================
function model = loadModel(pth)
S = load(pth,'wbPreset');
p = S.wbPreset;
model = newModel(p.bag);
model.stepOverrides = p.stepOverrides;
if isfield(p,'foKeys')
    for i = 1:numel(p.foKeys)
        model.fileOverrides(p.foKeys{i}) = p.foVals{i};
    end
end
if isfield(p,'tbKeys')
    for i = 1:numel(p.tbKeys), model.typeBag(p.tbKeys{i}) = p.tbVals{i}; end
end
if isfield(p,'toKeys')
    for i = 1:numel(p.toKeys), model.typeOverrides(p.toKeys{i}) = p.toVals{i}; end
end
model = migrateMaskLimits(model, presetOf(p));
end

% =====================================================================
function v = presetOf(p)
if isfield(p,'schema') && isnumeric(p.schema) && isscalar(p.schema), v = p.schema; else, v = 1; end
end

% =====================================================================
function model = migrateMaskLimits(model, version)
%migrateMaskLimits  Preset layout 2: runInternalCycle's reference-trace mask limits.
%   They used to be called trustLimitsK/trustLimitsI - the names the contrast step uses
%   for the mask it PROPAGATES - and are now maskLimitsK/maskLimitsI, while
%   trustLimitsK/trustLimitsI on the internal cycle mean the propagated mask, as they do
%   everywhere else.  An override stored under the old name is therefore MOVED rather
%   than dropped: it is the user's tuning of the pulse-detection area and it still means
%   exactly that, whereas left alone the same number would silently start masking the
%   OUTPUT.  wbSession runs the same migration for a saved session.
if version >= presetVersion(), return; end
model = renameStepField(model, 'internalCycle', 'trustLimitsK', 'maskLimitsK');
model = renameStepField(model, 'internalCycle', 'trustLimitsI', 'maskLimitsI');
end

% =====================================================================
function dst = overlay(dst, src)
if ~isstruct(src), return; end
fn = fieldnames(src);
for i = 1:numel(fn), dst.(fn{i}) = src.(fn{i}); end
end
