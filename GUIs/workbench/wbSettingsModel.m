%wbSettingsModel - Global settings bag, shared-key propagation, overrides, presets.
%
%   Holds the workbench's parameters and resolves the concrete settings struct s
%   for a given step (and optional file).  Layering, applied in order:
%     1. the step's own defaults          (step.presets.default)
%     2. shared-key bag                    (a shared key edited once propagates to
%                                           every step that reads it)
%     3. per-step override                 (shadows the bag for one step)
%     4. per-file override                 (shadows for one file+step)
%   The bag starts empty, so before the user touches a shared key each step keeps
%   its own launcher default (the launchers legitimately diverge, e.g. a
%   different trustLimitsK per step); the first edit of a shared key is what makes
%   it propagate.  Presets serialise the whole model to a .mat.
%
% Syntax:
%    model = wbSettingsModel('new')
%    model = wbSettingsModel('new', initBag)          % seed the shared bag
%    s     = wbSettingsModel('resolve', model, step)
%    s     = wbSettingsModel('resolve', model, step, fileModel)
%    model = wbSettingsModel('setShared', model, key, value)
%    model = wbSettingsModel('setStep',   model, stepId, field, value)
%    model = wbSettingsModel('setFile',   model, fileKey, stepId, field, value)
%            wbSettingsModel('save',  model, path)
%    model = wbSettingsModel('load',  path)
%
% Inputs:
%    step - a step struct from wbStepRegistry (carries presets/sharedKeys).
%    fileModel - a wbFileModel struct (its .identity keys per-file overrides).
%
% Outputs:
%    model - struct with fields: bag (struct), stepOverrides (struct of structs),
%            fileOverrides (containers.Map key 'identity||stepId' -> struct).
%    s     - the resolved settings struct for a run.
%
% See also: wbStepRegistry, wbStateEngine, wbFileModel
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbSettingsModel(action, varargin)

switch action
    case 'new',       out = newModel(varargin{:});
    case 'resolve',   out = resolve(varargin{:});
    case 'setShared', out = setShared(varargin{:});
    case 'setStep',   out = setStep(varargin{:});
    case 'setFile',   out = setFile(varargin{:});
    case 'save',      saveModel(varargin{:});  out = [];
    case 'load',      out = loadModel(varargin{:});
    otherwise, error('wbSettingsModel:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function model = newModel(initBag)
if nargin<1 || isempty(initBag), initBag = struct(); end
model = struct( ...
    'bag',           initBag, ...
    'stepOverrides', struct(), ...
    'fileOverrides', containers.Map('KeyType','char','ValueType','any'));
end

% =====================================================================
function s = resolve(model, step, fileModel)
% 1) step defaults
if isfield(step.presets,'default') && isstruct(step.presets.default)
    s = step.presets.default;
else
    s = struct();
end
% 2) shared-key bag (only the keys this step reads)
shared = step.sharedKeys;
for i = 1:numel(shared)
    k = shared{i};
    if isfield(model.bag, k), s.(k) = model.bag.(k); end
end
% libraryFolder is shared for every step
if isfield(model.bag,'libraryFolder'), s.libraryFolder = model.bag.libraryFolder; end
% 3) per-step override
if isfield(model.stepOverrides, step.id)
    s = overlay(s, model.stepOverrides.(step.id));
end
% 4) per-file override
if nargin>=3 && ~isempty(fileModel)
    key = [fileModel.identity '||' step.id];
    if isKey(model.fileOverrides, key)
        s = overlay(s, model.fileOverrides(key));
    end
end
end

% =====================================================================
function model = setShared(model, key, value)
model.bag.(key) = value;
end

% =====================================================================
function model = setStep(model, stepId, field, value)
if ~isfield(model.stepOverrides, stepId)
    model.stepOverrides.(stepId) = struct();
end
model.stepOverrides.(stepId).(field) = value;
end

% =====================================================================
function model = setFile(model, fileKey, stepId, field, value)
key = [fileKey '||' stepId];
if isKey(model.fileOverrides, key), ov = model.fileOverrides(key); else, ov = struct(); end
ov.(field) = value;
model.fileOverrides(key) = ov;
end

% =====================================================================
function saveModel(model, pth)
%saveModel  Serialise the model to a .mat (containers.Map -> plain struct).
wbPreset = model;
% store fileOverrides as parallel cell arrays so the preset is Map-free on disk
wbPreset.fileOverrides = [];
wbPreset.foKeys = keys(model.fileOverrides);
wbPreset.foVals = values(model.fileOverrides);
save(pth,'wbPreset','-v7.3');
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
end

% =====================================================================
function dst = overlay(dst, src)
if ~isstruct(src), return; end
fn = fieldnames(src);
for i = 1:numel(fn), dst.(fn{i}) = src.(fn{i}); end
end
