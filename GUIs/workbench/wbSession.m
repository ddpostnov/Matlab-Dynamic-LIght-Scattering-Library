%wbSession - THE durable, resumable Processing-Workbench session (.mat).
%
%   The session is the CONTRACT between the three programs, not a convenience: the
%   workbench writes it, and guiExport / guiExplore read it without the workbench
%   being open (spec §5).  It round-trips to a plain .mat so a long processing job
%   survives a MATLAB restart - or a crash, or a Stop - and always leaves something
%   resumable behind.  The session captures
%     * the SOURCE of the file set - root folder, glob, and the label regexps -
%       plus the curated path list itself, so a reload reproduces the Files table
%       with NO mandatory re-scan;
%     * the CURATION - the RESOLVED per-file labels (animal / type / index /
%       experimental group / modality / reference), plus the per-axis hand
%       overrides that produced them (keyed by path, so they survive a re-scan)
%       and one reference RECORDING IDENTITY per animal;
%     * the CONFIGURATION - which steps each recording TYPE runs, which per-animal
%       steps are on, and the per-type settings layers (the Constructor tab's whole
%       output, spec D4);
%     * the COMPLETION RECORD - per FILE and per STEP, what ran, when, and with
%       which settings fingerprint (see 'completed' below);
%     * the discovery grid + animal rows, the settings model (bag + step/file
%       overrides), the preset reference, and the INVALIDATION overlay (which
%       done cells an in-session settings edit pushed stale).
%   Nothing on disk is re-derived here - wbSession only (de)serialises; the caller
%   re-runs wbStateEngine after load to overlay saved state on the fresh disk
%   picture.
%
%   WHY THE LABELS ARE STORED RESOLVED.  Processing ignores the experimental group
%   entirely (spec §2), but Export and Explore are built on it, and they must not
%   have to re-run the regexps - or know the override rules - to get it.  So the
%   .files list carries the answer, and 'read' is the door they come in through:
%   neither tool reaches into the file layout, and a schema bump stays a change to
%   this one function.
%
%   KEYED BY PATH OR BY RECORDING IDENTITY.  Nine of the members below are - the
%   four .overrides axes and .modalityOvr by PATH, .completed by 'path||stepId',
%   .fileOverrides and .staleOverlay by 'identity||stepId', .animalRef and
%   .animalRefMan by an identity VALUE.  Renaming a recording on disk therefore has
%   to move all of them together or a session stops matching the files it describes;
%   guiWorkbench>adoptRename owns that, and a new keyed field added here has to be
%   added there too.  The SHAPE never changes, so a rename needs no schema bump.
%
%   PLAIN DATA ONLY.  Every containers.Map member is flattened to parallel
%   key/value cell arrays so the file is Map-free on disk (the trick
%   wbSettingsModel uses for its preset), and the whole payload is swept for
%   function handles before it is written (stripHandles, mirroring the wrappers'
%   reportSettings discipline) - a session must never serialise a closure over a
%   dead figure.
%
%   SCHEMA.  wbSessionData.schema names the layout version (currently 7; a file
%   without the field is the unversioned Phase-3 layout, read as 1).  Loading an
%   older session is supported: the fields it predates come back at their empty
%   defaults, so an old sidecar still restores its file set and settings - a
%   schema-2 file simply carries no type configuration, and the Constructor opens
%   with every type unticked.  Schema 4 changed the SHAPE of a typeSel key from
%   'type||stepId' to 'type||branch||stepId', because one raw recording can drive
%   two independent pipelines and each is configured separately; the keys are
%   written verbatim here and upgraded by wbTypeSelection('fromCells',keys,reg),
%   which is the only place that knows what a step's branch is - and which also
%   DROPS a key naming a step the registry no longer has, so a session outlives a
%   retired step without a schema bump (nothing about the layout changed, only
%   which ids are meaningful).  Schema 5 added
%   the resolved .files list and the .completed record.  Schema 6 DROPPED the
%   'checked' set: the per-cell run queue was retired with the per-file matrix
%   (Phase 6) and selection is now wholly the Constructor's typeSel/animalSel, so
%   the field is neither written nor read - a schema-5 sidecar's chkKeys are simply
%   ignored, and everything else about it still loads.  Schema 7 added .reprocess,
%   the Process tab's re-run switch; a schema-6 sidecar reads it as false, which is
%   the default and the safe answer.  Schema 8 is a MIGRATION rather than a new
%   field: runInternalCycle's reference-trace mask limits were renamed
%   trustLimitsK/trustLimitsI -> maskLimitsK/maskLimitsI, so a schema-7 sidecar's
%   step-scoped override of the old name is moved to the new one on load - see
%   migrateMaskLimits.  It is the one kind of change that MUST bump the schema even
%   though the layout is untouched: the same number under the same key would
%   otherwise keep being applied to a completely different parameter.
%
% Syntax:
%    wbSession('save', path, session)      % write the sidecar
%    session = wbSession('load', path)     % read it back (Maps rebuilt)
%    session = wbSession('read', path)     % load + validate (the consumer API)
%    s0      = wbSession('empty')          % a blank session struct
%    v       = wbSession('schema')         % the version this code writes
%    [tf,why]= wbSession('validate', session)
%
% Inputs:
%    path    - full path of the .mat sidecar.
%    session - struct with fields:
%       schema        double, the layout version (set by 'save').
%       root          char, the scanned root folder ('' if none).
%       glob          char, the scan glob ('*.rls', '*_K_d.mat', ...).
%       patterns      struct of regexps: animal/type/index/expGroup/ref.
%       paths         cellstr, the curated working set in table order.
%       files         1xN struct array, ONE ENTRY PER FILE, parallel to .paths:
%                       .path .name .animal .type .index .expGroup .modality
%                       .identity .branch .isRef .use
%                     The labels are RESOLVED (regexp + hand override already
%                     applied) - this is what guiExport / guiExplore read.  .use
%                     is always true: Phase 1 dropped the use flag, a file that
%                     will not be processed is DELETED from the working set, so
%                     the field survives only because spec §5 names it.
%       completed     containers.Map 'path||stepId' -> struct with
%                       .state ('done'|'error'|'skipped'), .when (datestr,
%                       ISO-ish local time), .fingerprint (char hash of the
%                       settings the step ran with), .message (error text or '').
%                     Written per FILE, per STEP, on every state change, so a
%                     killed run still resumes and a re-run can skip what is
%                     already done unless the fingerprint moved.
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
%       staleOverlay  containers.Map 'identity||stepId' -> logical true (the done
%                     cells an in-session settings edit re-opened).
%       presetRef     char, path of the last-used preset ('' if none).
%       reprocess     logical, the Process tab's 'Re-process finished files' tick:
%                     whether the next Run repeats work that is already on disk.
%                     It is remembered rather than treated as a window preference
%                     because it changes what a resumed session DOES - coming back
%                     to a half-finished protocol and silently re-running it would
%                     be the same trap the From-column coupling was.
%
% Outputs:
%    session - the same struct shape, with every Map reconstructed.
%    tf, why - 'validate': whether the struct is a usable session of a schema
%              this code understands, and a one-line reason when it is not.
%
% See also: wbSettingsModel, wbStateEngine, wbDiscoverFiles, wbTypeModel,
%           wbTypeSelection, guiWorkbench, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function [out, why] = wbSession(action, varargin)

why = '';
switch action
    case 'save',     saveSession(varargin{:}); out = [];
    case 'load',     out = loadSession(varargin{:});
    case 'read',     out = readSession(varargin{:});
    case 'empty',    out = emptySession();
    case 'schema',   out = schemaVersion();
    case 'validate', [out, why] = validateSession(varargin{:});
    otherwise,       error('wbSession:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function v = schemaVersion(), v = 8; end

% =====================================================================
function s = emptySession()
%emptySession  A blank session at the current schema (the load-time defaults).
s = struct();
s.schema        = schemaVersion();
s.root          = '';
s.glob          = '';
s.patterns      = refPatterns();
s.paths         = {};
s.files         = emptyFileList();
s.completed     = anyMap();
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
s.staleOverlay  = anyMap();
s.presetRef     = '';
s.reprocess     = false;
end

% =====================================================================
function p = refPatterns()
%refPatterns  The label patterns plus the reference regexp (not a label axis).
p = wbTypeModel('emptyPatterns');
p.ref = '';
end

% =====================================================================
function f = emptyFileList()
%emptyFileList  The 0x0 shape of session.files - THE list the read-only tools
%   consume.  Declared once so an empty session and a full one have the same
%   fields and a consumer can index it without a guard.
f = struct('path',{},'name',{},'animal',{},'type',{},'index',{}, ...
    'expGroup',{},'modality',{},'identity',{},'branch',{},'isRef',{},'use',{});
end

% =====================================================================
function f = fileListFrom(files)
%fileListFrom  Coerce whatever the caller supplied into the documented shape:
%   every field present, char-valued, in .paths order.  A caller that fills only
%   some fields (an older workbench, a hand-built session) still round-trips.
f = emptyFileList();
if isempty(files) || ~isstruct(files), return; end
proto = struct('path','','name','','animal','','type','','index','', ...
    'expGroup','','modality','','identity','','branch','','isRef',false,'use',true);
fn = fieldnames(proto);
for i = 1:numel(files)
    e = proto;
    for k = 1:numel(fn)
        if ~isfield(files(i),fn{k}), continue; end
        v = files(i).(fn{k});
        if islogical(proto.(fn{k})), e.(fn{k}) = logical(firstOr(v,false));
        else,                        e.(fn{k}) = charOf(v);
        end
    end
    f(end+1) = e; %#ok<AGROW>
end
end
function v = firstOr(x, d)
if isempty(x), v = d; else, v = x(1); end
end
function s = charOf(v)
if ischar(v), s = v; elseif isstring(v), s = char(v);
elseif isnumeric(v) || islogical(v), s = num2str(v); else, s = ''; end
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
wbSessionData.files         = fileListFrom(session.files);
wbSessionData.fNames        = session.fNames;
wbSessionData.referenceMode = session.referenceMode;
wbSessionData.animalNames   = session.animalNames;
wbSessionData.modality      = session.modality;
wbSessionData.rowOrder      = session.rowOrder;
wbSessionData.bag           = session.bag;
wbSessionData.stepOverrides = session.stepOverrides;
wbSessionData.presetRef     = session.presetRef;
wbSessionData.reprocess     = logical(session.reprocess);

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
[wbSessionData.staleKeys,  ~]                        = mapToCells(session.staleOverlay);
[wbSessionData.doneKeys,   wbSessionData.doneVals]   = mapToCells(session.completed);

% a session is read by two OTHER programs, long after this figure is gone: no
% closure over it may ever reach the file (spec §5)
wbSessionData = stripHandles(wbSessionData);
save(pth,'wbSessionData','-v7.3','-nocompression');
end

% =====================================================================
function session = loadSession(pth)
%loadSession  Read the sidecar, rebuild the Maps, default anything it predates.
S = load(pth,'wbSessionData');
p = S.wbSessionData;

session = emptySession();
if isfield(p,'schema'), session.schema = p.schema; else, session.schema = 1; end

plain = {'root','glob','patterns','paths','fNames','referenceMode', ...
         'animalNames','modality','rowOrder','bag','stepOverrides','presetRef', ...
         'reprocess'};
for i = 1:numel(plain)
    if isfield(p,plain{i}), session.(plain{i}) = p.(plain{i}); end
end
session.patterns = fillDefaults(session.patterns, refPatterns());
if isfield(p,'files'), session.files = fileListFrom(p.files); end

if isfield(p,'ovrAxes') && isfield(p,'ovrKeys')
    for i = 1:numel(p.ovrAxes)
        session.overrides.(p.ovrAxes{i}) = cellsToMap(p.ovrKeys{i}, p.ovrVals{i}, 'char');
    end
end
if isfield(p,'modKeys'),    session.modalityOvr  = cellsToMap(p.modKeys,    p.modVals,    'char'); end
if isfield(p,'refKeys'),    session.animalRef    = cellsToMap(p.refKeys,    p.refVals,    'char'); end
if isfield(p,'refManKeys'), session.animalRefMan = cellsToMap(p.refManKeys, p.refManVals, 'char'); end

% (schema 5 and older also carried p.chkKeys, the retired per-cell run queue: it is
% read by nothing now, so it is simply left in the file and not mapped back)
if isfield(p,'foKeys'),    session.fileOverrides = cellsToMap(p.foKeys,    p.foVals,              'any'); end
if isfield(p,'staleKeys'), session.staleOverlay  = cellsToMap(p.staleKeys, trueVals(p.staleKeys), 'any'); end

% ---- the Constructor's output (schema 3; a schema-2 file simply has none) ----
if isfield(p,'tbKeys'),   session.typeBag       = cellsToMap(p.tbKeys, p.tbVals, 'any'); end
if isfield(p,'toKeys'),   session.typeOverrides = cellsToMap(p.toKeys, p.toVals, 'any'); end
if isfield(p,'tselKeys'), session.typeSel   = cellsToMap(p.tselKeys, trueVals(p.tselKeys), 'any'); end
if isfield(p,'aselKeys'), session.animalSel = cellsToMap(p.aselKeys, trueVals(p.aselKeys), 'any'); end

% ---- the completion record (schema 5; older sidecars simply have none) -------
if isfield(p,'doneKeys'), session.completed = cellsToMap(p.doneKeys, p.doneVals, 'any'); end

session = migrateMaskLimits(session);
end

% =====================================================================
function session = migrateMaskLimits(session)
%migrateMaskLimits  Schema 8: runInternalCycle's reference-trace mask limits.
%   The step used to call them trustLimitsK/trustLimitsI - the names the contrast step
%   uses for the mask it PROPAGATES, which is why the two steps had to carry different
%   values under one shared key.  They are now maskLimitsK/maskLimitsI, and
%   trustLimitsK/trustLimitsI on the internal cycle mean the propagated mask, exactly as
%   they do on the contrast step.
%
%   An override recorded against the old name is MOVED, not dropped: it is the user's
%   tuning of which pixels the pulse is detected from, and it still means that.  Left
%   alone it would keep being applied - the same number, under the same key, doing a
%   completely different job.  wbSettingsModel owns the walk over the step-scoped
%   layers, and deliberately leaves the SHARED bags alone: there trustLimitsK is shared
%   by name with the contrast and segmentation steps, and there it always meant trust.
if session.schema >= schemaVersion(), return; end
m = struct('stepOverrides',session.stepOverrides, ...
           'fileOverrides',session.fileOverrides, ...
           'typeOverrides',session.typeOverrides);
m = wbSettingsModel('renameStepField', m, 'internalCycle','trustLimitsK','maskLimitsK');
m = wbSettingsModel('renameStepField', m, 'internalCycle','trustLimitsI','maskLimitsI');
session.stepOverrides = m.stepOverrides;
session.fileOverrides = m.fileOverrides;
session.typeOverrides = m.typeOverrides;
session.schema        = schemaVersion();
end

% =====================================================================
function session = readSession(pth)
%readSession  THE consumer door (spec §5).  guiExport / guiExplore call this, not
%   load(), so neither of them ever has to know the file layout - and a schema
%   bump stays a change to this function.  Errors rather than returning half a
%   session: a tool that got past this call may trust every documented field.
if isempty(pth) || ~isfile(char(pth))
    error('wbSession:noFile','No session file at "%s".', char(pth));
end
try
    session = loadSession(char(pth));
catch ME
    error('wbSession:unreadable','"%s" is not a workbench session (%s).', ...
        char(pth), ME.message);
end
[ok, why] = validateSession(session);
if ~ok, error('wbSession:invalid','"%s": %s', char(pth), why); end
end

% =====================================================================
function [tf, why] = validateSession(session)
%validateSession  Is this a usable session of a schema this code understands?
%   Deliberately shallow - it guards the CONTRACT (the fields a consumer is
%   promised, and a version it can read), not the contents, which are the user's.
tf = false; why = '';
if ~isstruct(session) || ~isscalar(session)
    why = 'not a scalar session struct'; return
end
need = {'schema','root','glob','patterns','paths','files','completed', ...
        'overrides','animalRef','typeSel','animalSel','bag'};
miss = need(~isfield(session, need));
if ~isempty(miss)
    why = ['missing field(s): ' strjoin(miss,', ')]; return
end
if ~isnumeric(session.schema) || ~isscalar(session.schema)
    why = 'schema is not a version number'; return
end
if session.schema > schemaVersion()
    why = sprintf('schema %g was written by a newer workbench (this one reads %g)', ...
        session.schema, schemaVersion()); return
end
if ~isempty(session.files) && ~isstruct(session.files)
    why = 'files is not a struct array'; return
end
tf = true;
end

% =====================================================================
function v = stripHandles(v)
%stripHandles  Recursively replace every function handle by its char text.
%   The wrappers already drop their transport hooks before writing a settings file
%   (Core/Reporting > reportSettings), the myograph wrappers included; a session is
%   read by other programs entirely, so it gets the same treatment for the whole
%   payload - one place, no exceptions.
%   The text is kept rather than deleted so a stray handle is visible in the file
%   instead of silently vanishing.
if isa(v,'function_handle')
    v = ['<fcn:' func2str(v) '>'];
elseif iscell(v)
    for i = 1:numel(v), v{i} = stripHandles(v{i}); end
elseif isstruct(v)
    fn = fieldnames(v);
    for k = 1:numel(v)
        for i = 1:numel(fn), v(k).(fn{i}) = stripHandles(v(k).(fn{i})); end
    end
end
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
%trueVals  A logical-true value for each key (the overlays are boolean sets).
v = num2cell(true(1,numel(k)));
end
