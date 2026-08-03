%wbStateEngine - Per-(file,step) state for the Processing Workbench matrix.
%
%   For one FILE (a wbFileModel handle) and the step registry, compute the state
%   of every step from three inputs - mirroring removeProcessedFiles but richer:
%     1. Disk markers.  A step is DONE when its settings gating field is present
%        in the SETTINGS files OF THAT FILE'S OWN PIPELINE.  Each step appends its
%        gating field to the settings and every downstream _s file carries it
%        forward, so the union runs along the product chain ('_t_K_s' ->
%        '_t_BFI_s'), which is what makes detection survive the BFI rename
%        (_K_->_BFI_) and a deleted original.  The REAL gating field per step is
%        used (01 A1).
%     2. Prerequisites.  Every id in step.requires must be done, else unavailable.
%     3. Settings fingerprint.  The current step settings are compared to the
%        stored settings.(gatingField); a mismatch flips a done step to stale.
%
%   A PREREQUISITE IS ALSO SATISFIED BY THIS RUN (opts.plannedIds, author
%   2026-07-31).  The disk is not the only source: a step whose producer is
%   SELECTED TO RUN in the same sequence will have its input by the time it is
%   reached, because the run order is step-major and puts the producer first.
%   Judging it on the disk alone made the workbench say "no input file" about a file
%   it was itself three lines away from writing.  Such a step therefore reads READY,
%   and NOT done - nothing has happened yet - with the reason saying exactly why it
%   is ready: 'input will be produced earlier in this run'.  With no plannedIds the
%   answer is the disk's alone, which is what every caller outside the live
%   workbench (the tests, a bare query) gets.
%
%   PER FILE, NOT PER RECORDING (spec D6/D8, author 2026-07-28).  One raw
%   recording can drive TWO independent pipelines - the raw producers each write
%   their own triplet ('_t_K', '_c_K') and everything later APPENDS to one of them
%   - so unioning the gating fields across everything named after the recording
%   made a half-processed recording read as fully done: a project where only the
%   cardiac product ever got a BFI reported BFI as done for the contrast file too.
%   The union is therefore confined to the queried file's own pipeline: settings
%   files sitting on ANOTHER raw-producer branch are excluded, while a file
%   derived from this one WITHIN the pipeline (the external cycle's '_e_K',
%   which is a new stage flag but not a new pipeline) still counts.  Which
%   branches are independent pipelines is read from the registry
%   (wbTypeSelection('branches')), never listed by name here.
%
%   A model with NO stage flag - a raw '.rls'/'.cxd' recording, the natural row of
%   a freshly scanned project - has no pipeline of its own yet and legitimately
%   stands for the whole recording, so it keeps the full union.
%
%   States: 'unavailable' | 'ready' | 'done' | 'stale' | 'error'.  ('error' is
%   set by the executor at run time, never here.)
%
% Syntax:
%    st = wbStateEngine(model, reg)
%    st = wbStateEngine(model, reg, curSettings)
%    st = wbStateEngine(model, reg, curSettings, opts)
%
% Inputs:
%    model       - a struct from wbFileModel (identifies the FILE, and through
%                  its .branch the pipeline whose settings are read).
%    reg         - the struct array from wbStepRegistry.
%    curSettings - (optional) struct whose fields are step ids, each a resolved
%                  settings struct s (from wbSettingsModel).  Used only for the
%                  staleness fingerprint; omit/[] to skip (done stays done).
%    opts        - (optional) struct with:
%                    .sData       containers.Map name->settings struct, injected in
%                                 place of a disk scan (for pure unit tests).
%                    .plannedIds  cellstr of step ids THIS FILE'S configuration will
%                                 run in the sequence about to be executed.  They
%                                 count as satisfied prerequisites (see above); they
%                                 never make a step read 'done'.  {} or absent = the
%                                 disk alone decides.
%
% Outputs:
%    st - 1xN struct array aligned to reg, with fields:
%           id, applicable (logical), state (char), reason (char).
%         'ready' carries a reason only when it was reached through plannedIds.
%
% See also: wbFileModel, wbStepRegistry, wbTypeSelection, wbInvalidate,
%           removeProcessedFiles
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function st = wbStateEngine(model, reg, curSettings, opts)

if nargin<3, curSettings = []; end
if nargin<4, opts = struct(); end

% ---- gather THIS FILE's accumulated settings fields -------------------------
[doneFields, storedMap] = gatherSettings(model, reg, opts);

% id -> gating field, for prerequisite testing
ids   = {reg.id};
gate  = {reg.gatingField};
gateOf = containers.Map(ids, gate);

% what the disk says has run, and what THIS RUN will add to it before it gets here
onDiskIds  = doneStepIds(reg, gateOf, doneFields);
plannedIds = plannedOf(opts);
availIds   = uniqueStable([onDiskIds, plannedIds]);

st = repmat(struct('id','','applicable',false,'state','','reason',''),1,numel(reg));

for k = 1:numel(reg)
    step = reg(k);
    st(k).id = step.id;
    st(k).applicable = any(strcmp(model.modality, step.modalities));

    if ~st(k).applicable
        st(k).state  = 'unavailable';
        st(k).reason = 'modality';
        continue
    end

    if isempty(step.gatingField)
        % a step that writes no settings field (e.g. export) is done once its
        % output artifact exists on disk
        doneOnDisk = strcmp(step.outKind,'none') && outputExists(model);
    else
        doneOnDisk = ismember(step.gatingField, doneFields);
    end
    % wbPrereqs stays THE single definition of "satisfied" - it is simply asked
    % twice, of two different sets, so the answer can say which one carried it
    requiresDone = wbPrereqs('met', step, availIds);
    requiresNow  = wbPrereqs('met', step, onDiskIds);
    hasInput = inputAvailable(step, model, requiresDone, doneOnDisk);

    if doneOnDisk
        stored = [];
        if isKey(storedMap, step.gatingField), stored = storedMap(step.gatingField); end
        if fingerprintMatches(step, curSettings, stored)
            st(k).state = 'done';
        else
            st(k).state  = 'stale';
            st(k).reason = 'settings changed';
        end
    elseif ~requiresDone
        st(k).state  = 'unavailable';
        st(k).reason = 'prerequisites not met';
    elseif ~hasInput
        st(k).state  = 'unavailable';
        st(k).reason = 'no input file';
    else
        st(k).state = 'ready';
        if ~requiresNow
            st(k).reason = 'input will be produced earlier in this run';
        end
    end
end
end

% =====================================================================
function ids = plannedOf(opts)
%plannedOf  The step ids this run will produce, as supplied by the host.  Coerced
%   rather than trusted: an absent field, [] and a string array all mean the same
%   thing here, and a caller with no plan at all is the normal case.
ids = {};
if ~isstruct(opts) || ~isfield(opts,'plannedIds') || isempty(opts.plannedIds), return; end
v = opts.plannedIds;
if ischar(v), v = {v}; end
if isstring(v), v = cellstr(v); end
if iscell(v), ids = reshape(cellfun(@char, v, 'UniformOutput', false), 1, []); end
end

% =====================================================================
function u = uniqueStable(c)
%uniqueStable  Unique cellstr, order preserved (unique(...,'stable') on a 1xN cell).
u = {};
if isempty(c), return; end
u = unique(c,'stable');
u = reshape(u,1,[]);
end

% =====================================================================
function [doneFields, storedMap] = gatherSettings(model, reg, opts)
%gatherSettings  Union of settings fields across THIS FILE's pipeline (D6/D8).
%   Same recording, and - when the queried file sits on a pipeline of its own -
%   not on a different one: see samePipeline.
doneFields = {};
storedMap  = containers.Map('KeyType','char','ValueType','any');

if isfield(opts,'sData') && ~isempty(opts.sData)
    % injected: containers.Map name -> settings struct
    names = keys(opts.sData);
    for i = 1:numel(names)
        settings = opts.sData(names{i});
        [doneFields, storedMap] = absorb(settings, doneFields, storedMap);
    end
    return
end

if isempty(model.folder) || ~isfolder(model.folder), return; end
rawBr = rawBranches(reg);
d = dir(fullfile(model.folder,'*_s.mat'));
for i = 1:numel(d)
    cand = wbFileModel(fullfile(d(i).folder, d(i).name));
    if ~strcmp(cand.role,'s'), continue; end
    if ~strcmp(cand.identity, model.identity), continue; end   % same recording only
    if ~samePipeline(model, cand, rawBr), continue; end        % same pipeline only
    try
        S = load(fullfile(d(i).folder, d(i).name),'settings');
    catch
        continue
    end
    if ~isfield(S,'settings'), continue; end
    [doneFields, storedMap] = absorb(S.settings, doneFields, storedMap);
end
end

% =====================================================================
function tf = samePipeline(model, cand, rawBranchList)
%samePipeline  Whether a settings file belongs to the queried file's pipeline.
%   A pipeline is what ONE raw producer started, so the only thing that puts a
%   settings file out of scope is sitting on a DIFFERENT raw-producer branch.
%   Everything else - the same branch, a stage the pipeline grew into later (the
%   external cycle's '_e'), or a flagless settings file - is the same pipeline
%   and its gating fields count.  A queried file that is itself flagless (a raw
%   recording) has no pipeline yet and takes them all.
tf = true;
if isempty(model.branch), return; end
if isempty(cand.branch) || strcmp(cand.branch, model.branch), return; end
tf = ~any(strcmp(cand.branch, rawBranchList));
end

% =====================================================================
function brs = rawBranches(reg)
%rawBranches  The branches that are INDEPENDENT pipelines - one per raw producer.
%   Read from the registry (wbTypeSelection owns the raw/derived split), so a
%   third entry step for a new modality is a registry edit and nothing here.
brs = {};
try
    brs = wbTypeSelection('branches', reg);
catch
    % a hand-built registry (tests, a trimmed array) may not carry the fields
    % wbTypeSelection reads; with no raw producers nothing is out of scope.
end
end

% =====================================================================
function [doneFields, storedMap] = absorb(settings, doneFields, storedMap)
%absorb  Fold one settings struct's gating fields into the running union.
if ~isstruct(settings), return; end
fn = fieldnames(settings);
for j = 1:numel(fn)
    if ~ismember(fn{j}, doneFields), doneFields{end+1} = fn{j}; end %#ok<AGROW>
    if ~isKey(storedMap, fn{j}), storedMap(fn{j}) = settings.(fn{j}); end
end
end

% =====================================================================
function ids = doneStepIds(reg, gateOf, doneFields)
%doneStepIds  Which STEPS the recording's settings say have run.  Prerequisites
%   are then judged by wbPrereqs, the single definition of "satisfied" (a step
%   needing ANY entry product must not demand the contrast one - see wbPrereqs).
ids = {};
for k = 1:numel(reg)
    id = reg(k).id;
    if ~isKey(gateOf, id), continue; end
    gf = gateOf(id);
    if ~isempty(gf) && ismember(gf, doneFields), ids{end+1} = id; end %#ok<AGROW>
end
end

% =====================================================================
function tf = inputAvailable(step, model, requiresDone, doneOnDisk)
%inputAvailable  Whether the step will have its input when it is reached.
%   requiresDone already counts the producers this run will run before it (see
%   plannedIds in the header), so a step whose input is a PRODUCT is available as
%   soon as its producer is on disk OR queued ahead of it.  An entry step reading a
%   RAW recording is the one exception that no plan can help: the recording either
%   sits beside the model or it does not.
if doneOnDisk, tf = true; return; end            % already ran => input existed
if isempty(step.requires)
    % entry step (contrast/internalCycle): needs the raw recording on disk
    ext = rawExtFromGlob(step.inGlob);
    if isempty(ext)
        tf = requiresDone;
    else
        tf = isfile(fullfile(model.folder,[model.stem ext]));
    end
else
    tf = requiresDone;                            % the prereq's output IS the input
end
% guided-style steps additionally need the raw recording alongside the product
if tf && step.needsRaw, tf = rawBeside(model); end
end

% =====================================================================
function tf = rawBeside(model)
%rawBeside  Is the recording this product came from still next to it?  The
%   candidate extensions are the ones wbFileModel knows as raw containers, in its
%   order - the modality vocabulary lives THERE and a second copy of it here would
%   silently stop recognising a modality the moment one is added.  It did: the pair
%   this replaced named '.rls' and '.cxd' only, so every needsRaw step of a video
%   modality read 'no input file' however plainly the video sat beside it.
tf = false;
if isempty(model.folder), return; end
exts = wbFileModel('extensions');
for i = 1:numel(exts)
    if isfile(fullfile(model.folder,[model.stem exts{i}])), tf = true; return; end
end
end

% =====================================================================
function tf = outputExists(model)
%outputExists  Whether the recording already has an exported workbook on disk.
%   exportToExcel writes '<input>_r.mat' -> '<input>.xlsx', so any
%   '<identity>*.xlsx' next to the recording means export has run.
tf = false;
if isempty(model.folder) || ~isfolder(model.folder), return; end
d = dir(fullfile(model.folder,[model.roiPrefix model.stem '*.xlsx']));
tf = ~isempty(d);
end

% =====================================================================
function ext = rawExtFromGlob(glob)
%rawExtFromGlob  '.rls'/'.cxd' from a raw inGlob like '*.rls'; '' for .mat globs.
ext = '';
if contains(glob,'.mat'), return; end
[~,~,e] = fileparts(glob);
ext = e;
end

% =====================================================================
function tf = fingerprintMatches(step, curSettings, stored)
%fingerprintMatches  Current step settings vs the stored run-time settings.
%   Compares only the step's own tunable fields present in BOTH structs; missing
%   comparisons (no current snapshot / older run) count as a match.
tf = true;
if isempty(curSettings) || isempty(stored) || ~isstruct(stored), return; end
if ~isfield(curSettings, step.id), return; end
cur = curSettings.(step.id);
if ~isstruct(cur), return; end

fields = stepFields(step);
ignore = {'libraryFolder','fNamesCopyTo','optimizer','metric','categories', ...
          'refFName','fName','stageFcn','cancelFcn'};
for i = 1:numel(fields)
    f = fields{i};
    if ismember(f, ignore), continue; end
    if isfield(cur,f) && isfield(stored,f)
        if ~isequaln(cur.(f), stored.(f)), tf = false; return; end
    end
end
end

% =====================================================================
function fields = stepFields(step)
%stepFields  The tunable field names a step reads (from its settingGroups + shared keys).
fields = {};
if ~isempty(step.settingGroups)
    fields = [step.settingGroups{:,2}];
end
fields = unique([fields, step.sharedKeys]);
end
