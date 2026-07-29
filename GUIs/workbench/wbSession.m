%wbSession - Save / load a Processing-Workbench session sidecar (.mat).
%
%   Round-trips the workbench's in-memory session state to a plain .mat file so a
%   long processing job survives a MATLAB restart AND so the read-only tools can
%   be handed a curated file set.  The session captures
%     * the SOURCE of the file set - root folder, glob, and the label regexps -
%       plus the curated path list itself, so a reload reproduces the Files table
%       with NO mandatory re-scan;
%     * the CURATION - the per-axis hand overrides (animal / type / index /
%       experimental group, and the modality correction, all keyed by path so
%       they survive a re-scan) and one reference RECORDING IDENTITY per animal;
%     * the CONFIGURATION - which steps each recording TYPE runs, which per-animal
%       steps are on, and the per-type settings layers (the Constructor tab's whole
%       output, spec D4);
%     * the discovery grid + animal rows, the settings model (bag + step/file
%       overrides), the preset reference, and the per-cell overlay (which cells
%       the user checked and which an in-session edit pushed stale).
%   Nothing on disk is re-derived here - wbSession only (de)serialises; the caller
%   re-runs wbStateEngine after load to overlay saved state on the fresh disk
%   picture.
%
%   Every containers.Map member is flattened to parallel key/value cell arrays so
%   the file is Map-free on disk (the trick wbSettingsModel uses for its preset),
%   and only plain data is written - no handles, no closures.
%
%   SCHEMA.  wbSessionData.schema names the layout version (currently 4; a file
%   without the field is the unversioned Phase-3 layout, read as 1).  Loading an
%   older session is supported: the fields it predates come back at their empty
%   defaults, so an old sidecar still restores its file set and settings - a
%   schema-2 file simply carries no type configuration, and the Constructor opens
%   with every type unticked.  Schema 4 changed the SHAPE of a typeSel key from
%   'type||stepId' to 'type||branch||stepId', because one raw recording can drive
%   two independent pipelines and each is configured separately; the keys are
%   written verbatim here and upgraded by wbTypeSelection('fromCells',keys,reg),
%   which is the only place that knows what a step's branch is.
%
% Syntax:
%    wbSession('save', path, session)      % write the sidecar
%    session = wbSession('load', path)     % read it back (Maps rebuilt)
%    s0      = wbSession('empty')          % a blank session struct
%
% Inputs:
%    path    - full path of the .mat sidecar.
%    session - struct with fields:
%       schema        double, the layout version (set by 'save').
%       root          char, the scanned root folder ('' if none).
%       glob          char, the scan glob ('*.rls', '*_K_d.mat', ...).
%       patterns      struct of regexps: animal/type/index/expGroup/ref.
%       paths         cellstr, the curated working set in table order.
%       overrides     struct of containers.Map, one per label axis, path->value.
%       modalityOvr   containers.Map path -> hand-corrected modality.
%       animalRef     containers.Map animal -> reference recording IDENTITY.
%       animalRefMan  containers.Map animal -> hand-pinned identity (manual wins).
%       fNames        2-D cell of char paths (the discovery grid).
%       referenceMode logical (column 1 is a reference/template).
%       animalNames   cellstr, one per animal row.
%       modality      char, the active modality driving the columns.
%       rowOrder      double vector, the display order of the flat rows ([]=natural).
%       bag           struct, the shared settings bag (wbSettingsModel).
%       stepOverrides struct of structs, per-step overrides.
%       fileOverrides containers.Map 'identity||stepId' -> struct.
%       typeBag       containers.Map type -> struct (per-type shared keys).
%       typeOverrides containers.Map 'type||stepId' -> struct (per-type settings).
%       typeSel       containers.Map 'type||branch||stepId' -> true (the row set).
%       animalSel     containers.Map stepId -> true (the per-animal steps that are on).
%       checked       containers.Map 'identity||stepId' -> logical true.
%       staleOverlay  containers.Map 'identity||stepId' -> logical true.
%       presetRef     char, path of the last-used preset ('' if none).
%
% Outputs:
%    session - the same struct shape, with every Map reconstructed.
%
% See also: wbSettingsModel, wbStateEngine, wbDiscoverFiles, wbTypeModel,
%           wbTypeSelection, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbSession(action, varargin)

switch action
    case 'save',  saveSession(varargin{:}); out = [];
    case 'load',  out = loadSession(varargin{:});
    case 'empty', out = emptySession();
    otherwise,    error('wbSession:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function v = schemaVersion(), v = 4; end

% =====================================================================
function s = emptySession()
%emptySession  A blank session at the current schema (the load-time defaults).
s = struct();
s.schema        = schemaVersion();
s.root          = '';
s.glob          = '';
s.patterns      = refPatterns();
s.paths         = {};
s.overrides     = wbTypeModel('emptyOverrides');
s.modalityOvr   = charMap();
s.animalRef     = charMap();
s.animalRefMan  = charMap();
s.fNames        = {};
s.referenceMode = false;
s.animalNames   = {};
s.modality      = 'LSCI';
s.rowOrder      = [];
s.bag           = struct();
s.stepOverrides = struct();
s.fileOverrides = anyMap();
s.typeBag       = anyMap();
s.typeOverrides = anyMap();
s.typeSel       = anyMap();
s.animalSel     = anyMap();
s.checked       = anyMap();
s.staleOverlay  = anyMap();
s.presetRef     = '';
end

% =====================================================================
function p = refPatterns()
%refPatterns  The label patterns plus the reference regexp (not a label axis).
p = wbTypeModel('emptyPatterns');
p.ref = '';
end

% =====================================================================
function saveSession(pth, session)
%saveSession  Flatten every Map and write a plain struct as 'wbSessionData'.
session = fillDefaults(session, emptySession());

wbSessionData = struct();
wbSessionData.schema        = schemaVersion();
wbSessionData.root          = session.root;
wbSessionData.glob          = session.glob;
wbSessionData.patterns      = session.patterns;
wbSessionData.paths         = session.paths;
wbSessionData.fNames        = session.fNames;
wbSessionData.referenceMode = session.referenceMode;
wbSessionData.animalNames   = session.animalNames;
wbSessionData.modality      = session.modality;
wbSessionData.rowOrder      = session.rowOrder;
wbSessionData.bag           = session.bag;
wbSessionData.stepOverrides = session.stepOverrides;
wbSessionData.presetRef     = session.presetRef;

% ---- label overrides: one key/value cell-array pair per axis ----------------
ax = wbTypeModel('axes');
wbSessionData.ovrAxes = ax;
wbSessionData.ovrKeys = cell(1,numel(ax));
wbSessionData.ovrVals = cell(1,numel(ax));
for i = 1:numel(ax)
    m = charMap();
    if isfield(session.overrides,ax{i}), m = session.overrides.(ax{i}); end
    [wbSessionData.ovrKeys{i}, wbSessionData.ovrVals{i}] = mapToCells(m);
end

[wbSessionData.modKeys,    wbSessionData.modVals]    = mapToCells(session.modalityOvr);
[wbSessionData.refKeys,    wbSessionData.refVals]    = mapToCells(session.animalRef);
[wbSessionData.refManKeys, wbSessionData.refManVals] = mapToCells(session.animalRefMan);
[wbSessionData.foKeys,     wbSessionData.foVals]     = mapToCells(session.fileOverrides);
[wbSessionData.tbKeys,     wbSessionData.tbVals]     = mapToCells(session.typeBag);
[wbSessionData.toKeys,     wbSessionData.toVals]     = mapToCells(session.typeOverrides);
[wbSessionData.tselKeys,   ~]                        = mapToCells(session.typeSel);
[wbSessionData.aselKeys,   ~]                        = mapToCells(session.animalSel);
[wbSessionData.chkKeys,    ~]                        = mapToCells(session.checked);
[wbSessionData.staleKeys,  ~]                        = mapToCells(session.staleOverlay);

save(pth,'wbSessionData','-v7.3');
end

% =====================================================================
function session = loadSession(pth)
%loadSession  Read the sidecar, rebuild the Maps, default anything it predates.
S = load(pth,'wbSessionData');
p = S.wbSessionData;

session = emptySession();
if isfield(p,'schema'), session.schema = p.schema; else, session.schema = 1; end

plain = {'root','glob','patterns','paths','fNames','referenceMode', ...
         'animalNames','modality','rowOrder','bag','stepOverrides','presetRef'};
for i = 1:numel(plain)
    if isfield(p,plain{i}), session.(plain{i}) = p.(plain{i}); end
end
session.patterns = fillDefaults(session.patterns, refPatterns());

if isfield(p,'ovrAxes') && isfield(p,'ovrKeys')
    for i = 1:numel(p.ovrAxes)
        session.overrides.(p.ovrAxes{i}) = cellsToMap(p.ovrKeys{i}, p.ovrVals{i}, 'char');
    end
end
if isfield(p,'modKeys'),    session.modalityOvr  = cellsToMap(p.modKeys,    p.modVals,    'char'); end
if isfield(p,'refKeys'),    session.animalRef    = cellsToMap(p.refKeys,    p.refVals,    'char'); end
if isfield(p,'refManKeys'), session.animalRefMan = cellsToMap(p.refManKeys, p.refManVals, 'char'); end

session.fileOverrides = cellsToMap(p.foKeys,    p.foVals,              'any');
session.checked       = cellsToMap(p.chkKeys,   trueVals(p.chkKeys),   'any');
session.staleOverlay  = cellsToMap(p.staleKeys, trueVals(p.staleKeys), 'any');

% ---- the Constructor's output (schema 3; a schema-2 file simply has none) ----
if isfield(p,'tbKeys'),   session.typeBag       = cellsToMap(p.tbKeys, p.tbVals, 'any'); end
if isfield(p,'toKeys'),   session.typeOverrides = cellsToMap(p.toKeys, p.toVals, 'any'); end
if isfield(p,'tselKeys'), session.typeSel   = cellsToMap(p.tselKeys, trueVals(p.tselKeys), 'any'); end
if isfield(p,'aselKeys'), session.animalSel = cellsToMap(p.aselKeys, trueVals(p.aselKeys), 'any'); end
end

% =====================================================================
function s = fillDefaults(s, d)
%fillDefaults  Add any field of the default struct the caller did not supply.
if ~isstruct(s), s = d; return; end
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s,f{i}), s.(f{i}) = d.(f{i}); end
end
end

% =====================================================================
function m = charMap(), m = containers.Map('KeyType','char','ValueType','char'); end
function m = anyMap(),  m = containers.Map('KeyType','char','ValueType','any');  end

% =====================================================================
function [k, v] = mapToCells(m)
%mapToCells  containers.Map -> parallel cell arrays (Map-free on disk).
if isa(m,'containers.Map') && m.Count > 0
    k = keys(m);
    v = values(m);
else
    k = {};
    v = {};
end
end

% =====================================================================
function m = cellsToMap(k, v, vtype)
%cellsToMap  Parallel cell arrays -> containers.Map (char keys, any/char values).
if nargin<3, vtype = 'any'; end
m = containers.Map('KeyType','char','ValueType',vtype);
for i = 1:numel(k)
    m(k{i}) = v{i};
end
end

% =====================================================================
function v = trueVals(k)
%trueVals  A logical-true value for each key (checked/stale are boolean sets).
v = num2cell(true(1,numel(k)));
end
