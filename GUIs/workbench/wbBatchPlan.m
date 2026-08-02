%wbBatchPlan - Turn an ordered workbench entry list into WRAPPER CALLS.
%
%   Every wrapper in Wrappers/ is a MULTI-FILE routine, and the launchers call it
%   that way: one getFileNamesList cell, one call, the whole list.  The workbench
%   expands its configuration into one work item per (recording, step), which is
%   the right unit for the matrix and for resuming - but it is NOT the unit a
%   wrapper is written for.  Handing a multi-file routine a 1x1 list makes
%   everything it does BETWEEN files unreachable: setRegions can never carry a
%   drawn ROI from one recording of an animal to the next, and an external cycle
%   can never index its per-file stimulus list.
%
%   This module is the translation.  It takes the entry list exactly as
%   guiWorkbench>buildRunOrder produced it and answers, once, the only question
%   the run loop should not have to: WHICH FILES DOES THIS CALL GET.  It owns all
%   file resolution - the input glob, the branch fan-out, the copy targets, the
%   raw partner, the per-animal reference - so wbExecutor is left as a thin loop
%   over calls.
%
%   FOLDING (spec D1/D3).  One wrapper call takes ONE settings struct, and the
%   workbench resolves settings per (type, step), so a call can only span
%   recordings whose resolved s is identical.  Adjacent perFile entries of the
%   same step therefore fold when their settings compare isequaln once the
%   per-call fields (fName, fNamesCopyTo, refFName) and the transport hooks are
%   removed; two types configured differently keep their own calls, two types
%   that happen to resolve the same merge into one.  setRegions has no settings
%   at all, so all of its recordings land in a single call - literally the
%   launcher.  perAnimal entries are never folded: they already ARE the launcher's
%   per-animal call.  Only RUNS OF ADJACENT entries fold, so a batch can never
%   reorder execution, and the list is never re-sorted (the run order is
%   step-major and already contiguous per step).
%
%   WHAT A FOLD IS FOR DEPENDS ON THE STEP'S fileOrder.  For an 'ordered' step the
%   fold IS the call - the whole point, since that is the launcher cell its
%   cross-file work needs.  For an 'independent' one the executor invokes the
%   wrapper once per recording out of the same fold (see wbExecutor), because there
%   is by declaration nothing happening between the files; the fold then remains the
%   unit of settings, of the log header and of the error boundary, and each
%   recording's own files and copy targets are carried on its unit.
%
%   SHAPING (spec §5).  How the files are laid out is the step's own declaration,
%   registry field fanOut:
%     'flat'   - fNames is a COLUMN of every file of every recording in the batch,
%                in batch order, and fNamesCopyTo is the runSegmentation
%                convention: one ROW per source file, its targets across the
%                columns.
%     'animal' - fNames is 2-D, one ROW per animal and one column per file of that
%                animal, ragged rows padded with '', and fNamesCopyTo is the
%                setRegions ELEMENTWISE convention: the same size as fNames, each
%                element a char or a nested cellstr.  This is the shape whose rows
%                a carry-forward loop walks.
%
%   branchScope is unchanged and still resolves PER RECORDING: 'one' contributes
%   the stage-preferred file, 'all' every branch product, 'copy' the contrast file
%   with its siblings routed into the copy list.
%
%   THE WORKING SET IS A FENCE (round-2 item 8).  branchScope says how many of a
%   recording's branch products a step covers; it does NOT say which products the
%   SESSION is about.  The input glob below is a directory listing filtered only by
%   identity, so for an 'all' or 'copy' step every product on disk joined the
%   wrapper's fNames - including a '_c_K' triplet a project processed last month
%   left behind, in a session that only ever configured the contrast branch.  The
%   host now hands the answer down through ctx.admits, built from the very
%   wbTypeSelection rows the derived monitor table is built from (wbProducts), so a
%   file that is not in the table cannot reach a wrapper either.  ABSENT HOOK = NO
%   FENCE: a ctx without .admits (the headless tests, any other caller) resolves
%   exactly as it always did.  When the fence is what emptied the list the caller
%   is told SO, in its own words - an operator must never be told a file is missing
%   while it is sitting in the folder.
%
%   PLANNING AND RESOLVING ARE TWO DIFFERENT MOMENTS (author, 2026-07-31).  Until
%   now this module answered "which files does this call get" for EVERY entry in one
%   pass, before the run started - and a step whose input its own run had not
%   written yet was therefore reported 'input file not found' and marked as an error
%   before anything had been attempted.  The run order is step-major, so contrast
%   runs three lines above the setRegions entry that consumes its '_t_K_d': the
%   ORDER was right all along, only the timing of the question was wrong.
%   So the two halves are now separate actions:
%     'build'   is pure arithmetic over the entry list - the registry, the models
%               and the resolved settings.  IT TOUCHES NO DISK.  It folds and
%               shapes the calls and leaves fNames/copyTo/rawNames/refName unset.
%     'resolve' is the disk half - the input glob, the branch fan-out, the copy
%               targets, the raw partner, the per-animal reference and the layout -
%               and the executor runs it IMMEDIATELY BEFORE each call, when the disk
%               finally reflects what the earlier steps of the same run wrote.
%   A batch whose input is still missing at that moment is the caller's to report;
%   this module only says so (ok=false) and the run carries on with the next call.
%
%   Pure planning either way: it reads the registry, the models and the disk, and it
%   changes no state, marks no cell and calls no wrapper.
%
% Syntax:
%    batches            = wbBatchPlan('build', entries, ctx)
%    [batches, skipped] = wbBatchPlan('build', entries, ctx)
%    [b, ok, reason]    = wbBatchPlan('resolve', b, ctx)
%
% Inputs:
%    entries - struct array from guiWorkbench>buildRunOrder, already step-major
%              (registry step -> animal -> reference-first -> file).  Each entry
%              has at least: stepId, identity, animalIdx, arity, label, animal.
%    b       - ONE batch from 'build', to be resolved against the disk as it is now.
%    ctx     - the wbExecutor context struct.  This module uses .reg, .modelOf,
%              .animalModels, .resolve, .contrastStage and the OPTIONAL .refFile
%              and .admits(identity,path) (see wbExecutor for what each of them
%              promises).
%
% Outputs:
%    batches - 1xN struct array, one element per WRAPPER CALL, in run order:
%       stepId   the step's id
%       step     the registry step struct (carries the wrapper handle)
%       arity    'perFile' | 'perAnimal'
%       fanOut   'flat' | 'animal' - the shape below
%       fileOrder 'independent' | 'ordered' - whether the call may be cut into one
%                invocation per recording (registry fileOrder; see wbExecutor)
%       label    how the call names itself in the log
%       entries  the folded entries, in order (afterDone is owed one call each)
%       heads    one model per ENTRY: the recording a perFile entry is, or the
%                reference a perAnimal entry starts at
%       models   every recording model the call touches, in order (the cells to
%                mark running / done / error)
%       units    1xK cell, ONE PER FOLDED ENTRY, each a struct with .entry .head
%                .models and - after 'resolve' - the .fNames / .copyTo / .rawNames
%                / .refName of that entry ALONE, shaped exactly as the whole batch
%                is.  This is what makes a per-recording invocation possible without
%                a second definition of the layout.
%       s        the shared resolved settings, hooks NOT yet injected
%       fNames   the wrapper's file list, shaped per fanOut       ('resolve' only)
%       rawNames the parallel raw recordings for a needsRaw step, else {}
%       copyTo   the s.fNamesCopyTo list, shaped per fanOut, {} when nothing is
%                inherited
%       refName  the per-animal reference file ('' when the call has none)
%    skipped - 1xM struct array of entries that could not become part of a call.
%       STRUCTURAL PROBLEMS ONLY now ('unknown step', 'no recording'): a missing
%       input is no longer knowable at this point, and is reported by 'resolve'.
%       Fields: entry, stepId, stepLabel, who, reason, mark (true = the caller
%       should mark the cell as an error).
%    ok      - 'resolve': whether every unit of the batch found its input.
%    reason  - 'resolve': why it did not, else ''.  TWO different answers, because
%              they mean two different things to the operator: nothing on disk at
%              all, or a file that is there but is not part of this session.
%
% Notes:
%    * The entry list is authoritative for ORDER.  Nothing here sorts, and folding
%      only ever merges neighbours, so the launcher property the run order exists
%      to give - every recording reaches a level before anything moves past it -
%      survives untouched.
%    * A batch's settings are the FIRST folded entry's.  That is safe precisely
%      because the fold key is the settings themselves; the discarded copies were
%      equal.
%    * 'resolve' IS ALL-OR-NOTHING for one batch, because a call is: the file list a
%      wrapper receives is one object, and half of it is not a launcher cell.  A
%      batch that loses its input therefore reports the whole call, and the executor
%      marks that call's recordings; the other calls of the same step are untouched.
%
% See also: wbExecutor, wbStepRegistry, wbRunRange, wbRefBranch, wbFileModel,
%           guiWorkbench, setRegions, runSegmentation
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function varargout = wbBatchPlan(action, varargin)

switch lower(char(action))
    case 'build'
        [varargout{1:max(1,nargout)}] = buildBatches(varargin{:});
    case 'resolve'
        [varargout{1:max(1,nargout)}] = resolveBatch(varargin{:});
    otherwise
        error('wbBatchPlan:action','Unknown action ''%s''.', char(action));
end
end

% =====================================================================
function [batches, skipped] = buildBatches(entries, ctx)
%buildBatches  Fold the entry list into wrapper CALLS.  No disk access anywhere in
%   here: what a call is - which entries share it, which settings it takes, how its
%   files will be laid out - is decided by the registry and the resolved settings
%   alone, all of which are knowable before the run starts.  WHICH files is asked
%   later, per call, by 'resolve'.
batches = emptyBatch(); skipped = emptySkip();
if isempty(entries), return; end

% ---- 1. one UNIT per entry: its models and its settings -----------------------
U = {}; sk = {};
for i = 1:numel(entries)
    e    = entries(i);
    step = stepById(ctx.reg, e.stepId);
    if isempty(step)
        sk{end+1} = mkSkip(e, e.stepId, e.stepId, whoOf(e), 'unknown step', false); %#ok<AGROW>
        continue
    end
    if strcmp(e.arity,'perAnimal')
        models = ctx.animalModels(e.animalIdx);
    else
        models = ctx.modelOf(e.identity);
    end
    if isempty(models)
        sk{end+1} = mkSkip(e, step.id, step.label, whoOf(e), 'no recording', false); %#ok<AGROW>
        continue
    end
    mHead = models(1);
    s     = ctx.resolve(step, mHead);
    U{end+1} = newUnit(e, step, s, mHead, models); %#ok<AGROW>
end
if ~isempty(sk), skipped = [sk{:}]; end
if isempty(U), return; end

% ---- 2. fold adjacent perFile units of one step with identical settings -------
B = {}; i = 1;
while i <= numel(U)
    j = i;
    if strcmp(U{i}.entry.arity,'perFile')
        while j < numel(U) && foldable(U{i}, U{j+1}), j = j + 1; end
    end
    B{end+1} = shapeBatch(U(i:j)); %#ok<AGROW>
    i = j + 1;
end
batches = [B{:}];
end

% =====================================================================
function [b, ok, reason] = resolveBatch(b, ctx)
%resolveBatch  THE DISK HALF, asked at the last possible moment.  Every unit of the
%   batch contributes its concrete files, and the batch is laid out per fanOut - the
%   same shaping as before, only now run when the steps ahead of this one in the
%   very same run have already written their products.
%
%   Each unit is ALSO shaped on its own, so a step the registry calls file-order
%   'independent' can be invoked once per recording without anyone re-deriving what
%   that recording's file list or copy list should look like.
ok = false; reason = '';
if isempty(b) || isempty(b.units), reason = 'no recording'; return; end
step = b.step;

for k = 1:numel(b.units)
    u      = b.units{k};
    cstage = ctx.contrastStage(u.head);
    [files, copies, raws, found, why] = entryFiles(step, u.models, cstage, ...
        refFileFor(ctx, u.entry, step), ctx);
    if ~found, reason = why; return; end
    u.files = files; u.copies = copies; u.raws = raws;
    [u.fNames, u.copyTo, u.rawNames] = shapeUnits({u}, step);
    u.refName = refOf(b.arity, u.fNames);
    b.units{k} = u;
end

[b.fNames, b.copyTo, b.rawNames] = shapeUnits(b.units, step);
b.refName = refOf(b.arity, b.fNames);
ok = true;
end

% =====================================================================
function u = newUnit(e, step, s, head, models)
%newUnit  One entry's share of a call, before the disk has been consulted.  The file
%   fields are filled in by 'resolve' and are {} until then.
u = struct('entry',e, 'step',step, 's',s, 'head',head, 'models',models, ...
    'key',foldKey(s), 'files',{{}}, 'copies',{{}}, 'raws',{{}}, ...
    'fNames',{{}}, 'copyTo',{{}}, 'rawNames',{{}}, 'refName','');
end

% =====================================================================
function [fNames, copyTo, rawNames] = shapeUnits(units, step)
%shapeUnits  Lay a list of resolved units out the way the step declares (fanOut).
%   ONE definition, used for the whole batch and for a single unit alike, so a
%   per-recording invocation gets exactly the shape the wrapper expects.
switch fanOutOf(step)
    case 'animal'
        [fNames, copyTo, rawNames] = shapeByAnimal(units, step);
    otherwise
        [fNames, copyTo, rawNames] = shapeFlat(units, step);
end
end

% =====================================================================
function r = refOf(arity, fNames)
%refOf  The per-animal reference file of a laid-out list ('' when there is none):
%   the first file, whichever shape it was laid out in.
r = '';
if strcmp(arity,'perAnimal') && ~isempty(fNames), r = fNames{1}; end
end

% =====================================================================
function tf = foldable(a, b)
%foldable  Whether unit b can join unit a's call: the same perFile step, and the
%   same settings once the per-call fields are out of the way.  Nothing else is
%   consulted - a step that must not merge says so by resolving differently.
tf = strcmp(b.entry.arity,'perFile') && strcmp(a.entry.stepId, b.entry.stepId) && ...
     isequaln(a.key, b.key);
end

% =====================================================================
function k = foldKey(s)
%foldKey  The settings a fold compares.  The fields a CALL owns rather than a
%   configuration - the file it is pointed at, the list it inherits onto, the
%   reference it templates on - are dropped, along with the transport hooks, which
%   are closures and would never compare equal anyway.
k = s;
for f = {'fName','fNamesCopyTo','refFName','stageFcn','cancelFcn'}
    if isfield(k,f{1}), k = rmfield(k,f{1}); end
end
if numel(fieldnames(k)) > 1, k = orderfields(k); end
end

% =====================================================================
function b = shapeBatch(units)
%shapeBatch  One run of folded units becomes ONE call.  The files are NOT part of
%   it yet - only what the call is: its step, its arity, its declared layout and
%   order sensitivity, its shared settings, and the recordings it will touch.
step = units{1}.step;
b           = newBatch();
b.stepId    = step.id;
b.step      = step;
b.arity     = units{1}.entry.arity;
b.fanOut    = fanOutOf(step);
b.fileOrder = orderOf(step);
b.s         = units{1}.s;

ents = cell(1,numel(units)); heads = cell(1,numel(units)); mdls = cell(1,numel(units));
for k = 1:numel(units)
    ents{k}  = units{k}.entry;
    heads{k} = units{k}.head;
    mdls{k}  = units{k}.models;
end
b.entries = [ents{:}];
b.heads   = [heads{:}];
b.models  = [mdls{:}];
b.units   = units;
b.label   = batchLabel(units);
end

% =====================================================================
function [fNames, copyTo, rawNames] = shapeFlat(units, step)
%shapeFlat  The default: a COLUMN of every file of every recording, in batch
%   order, with the copy list in runSegmentation's row-per-source-file convention
%   (its copyTargets reads row fidx, so the rows must line up with fNames).
files = {}; copies = {}; raws = {};
for k = 1:numel(units)
    files  = [files,  units{k}.files];   %#ok<AGROW>
    copies = [copies, units{k}.copies];  %#ok<AGROW>
    raws   = [raws,   units{k}.raws];    %#ok<AGROW>
end
fNames   = files(:);
rawNames = {};
if step.needsRaw, rawNames = raws(:); end
copyTo   = rowCopyList(copies);
end

% =====================================================================
function [fNames, copyTo, rawNames] = shapeByAnimal(units, step)
%shapeByAnimal  ROWS = ANIMALS.  A wrapper that carries state across the columns
%   of a row (setRegions and its ROIs) needs the row to BE the animal, which is
%   what the launcher's getFileNamesList(root,'*_t_K_d.mat','[A-Z]+\d+') hands it.
%   Rows are the runs of adjacent units sharing an animal - the run order groups an
%   animal's recordings together, so this never reorders anything - and ragged rows
%   are padded with '' (setRegions skips empty cells).
rows = {}; rowCopies = {}; rowRaws = {};
k = 1;
while k <= numel(units)
    a = units{k}.entry.animalIdx;
    f = {}; c = {}; r = {};
    while k <= numel(units) && isequal(units{k}.entry.animalIdx, a)
        f = [f, units{k}.files];   %#ok<AGROW>
        c = [c, units{k}.copies];  %#ok<AGROW>
        r = [r, units{k}.raws];    %#ok<AGROW>
        k = k + 1;
    end
    rows{end+1} = f; rowCopies{end+1} = c; rowRaws{end+1} = r; %#ok<AGROW>
end

nR = numel(rows);
nC = max(cellfun(@numel, rows));
fNames = repmat({''}, nR, nC);
copyTo = repmat({''}, nR, nC);
rawNames = repmat({''}, nR, nC);
anyCopy = false; anyRaw = false;
for i = 1:nR
    for j = 1:numel(rows{i})
        fNames{i,j} = rows{i}{j};
        t = rowCopies{i}{j};
        if ~isempty(t)
            anyCopy = true;
            if isscalar(t), copyTo{i,j} = t{1}; else, copyTo{i,j} = t(:)'; end
        end
        if j <= numel(rowRaws{i}) && ~isempty(rowRaws{i}{j})
            anyRaw = true; rawNames{i,j} = rowRaws{i}{j};
        end
    end
end
if ~anyCopy, copyTo = {}; end
if ~step.needsRaw || ~anyRaw, rawNames = {}; end
end

% =====================================================================
function copyTo = rowCopyList(copies)
%rowCopyList  The runSegmentation convention: N rows (one per source file) x the
%   widest target list, padded with ''.  {} when nothing is inherited at all, so a
%   step that copies nowhere never grows an s.fNamesCopyTo field.
copyTo = {};
nT = max([0, cellfun(@numel, copies)]);
if nT == 0, return; end
copyTo = repmat({''}, numel(copies), nT);
for i = 1:numel(copies)
    t = copies{i};
    for j = 1:numel(t), copyTo{i,j} = t{j}; end
end
end

% =====================================================================
function [files, copies, raws, ok, why] = entryFiles(step, models, cstage, refPath, ctx)
%entryFiles  ONE entry's contribution to a call: the concrete files it brings, the
%   copy targets each of them carries, and the raw partner of each (needsRaw only).
%   Everything is 1xK and parallel, so the shapers can lay it out either way.
%
%   perAnimal - the animal's members in order, each contributing as many branch
%   products as branchScope asks for ('one' the stage-preferred file, 'all' every
%   product, which is the launcher's branch-wide 'Roi*_K_d.mat' cell).  refPath,
%   when the host resolved one, IS the first file: the branch of the pinned
%   reference recording the step declared through refBranch.  A member with no
%   input on disk means the animal is not ready, so the whole entry waits.
%   perFile  - one recording, whose co-registered branch products are resolved per
%   branchScope: 'one' the stage-preferred file alone, 'all' every branch (the
%   wrapper's own loop covers them), 'copy' the stage-preferred file with the rest
%   returned as its copy targets.
files = {}; copies = {}; raws = {}; ok = false; why = '';
if nargin < 4, refPath = ''; end
if nargin < 5, ctx = struct(); end

if strcmp(step.arity,'perAnimal')
    wide = ~strcmp(scopeOf(step),'one');
    for k = 1:numel(models)
        [p, fenced] = resolveStepInputs(models(k), step, cstage, ctx);  % preferred first
        if k==1 && ~isempty(refPath) && isfile(refPath)
            p = [{refPath}, p(~strcmp(p,refPath))];           % the pinned reference leads
        end
        if isempty(p), why = missingReason(fenced); return; end   % a member has no input
        if ~wide, p = p(1); end
        files  = [files, p];                                  %#ok<AGROW>
        copies = [copies, repmat({{}},1,numel(p))];           %#ok<AGROW>
        raws   = [raws, repmat({rawPathFor(models(k))},1,numel(p))]; %#ok<AGROW>
    end
    ok = true; return
end

[p, fenced] = resolveStepInputs(models(1), step, cstage, ctx);
if isempty(p), why = missingReason(fenced); return; end
switch scopeOf(step)
    case 'all'
        files = p; copies = repmat({{}},1,numel(p));          % every branch, one call
    case 'copy'
        files = p(1); copies = {p(2:end)};                    % run here, the rest inherit
    otherwise
        files = p(1); copies = {{}};
end
raws = repmat({rawPathFor(models(1))},1,numel(files));
ok = true;
end

% =====================================================================
function p = refFileFor(ctx, e, step)
%refFileFor  The host's resolved reference FILE for this entry ('' when there is
%   none).  Only a per-animal step has a reference column, and only the host knows
%   which recording is pinned, so this is a pure look-up through ctx - the branch
%   rule itself lives in wbRefBranch and is never re-derived here.
p = '';
if ~strcmp(step.arity,'perAnimal'), return; end
if ~isfield(ctx,'refFile') || ~isa(ctx.refFile,'function_handle'), return; end
try
    p = char(ctx.refFile(e.animalIdx, step.id));
catch
    p = '';
end
end

% =====================================================================
function [p, fenced] = resolveStepInputs(model, step, cstage, ctx)
%resolveStepInputs  EVERY concrete _d.mat (or raw) file this step could consume for
%   this recording, ORDERED with the stage-preferred branch first.  Located by base
%   name + the step's input glob and filtered to this recording's identity AND to
%   the session's own working set (ctx.admits - see the header); mirrors
%   getFileNamesList scoped to one recording, so it is robust to the BFI rename and
%   the t|s|c|e branches.  Callers take p{1} for a single-branch step and the whole
%   list (or its tail) for the 'all' / 'copy' scopes - see entryFiles.  {} = the
%   step has no input it may use, and 'fenced' says which of the two that was: a
%   folder with nothing in it, or one whose candidates all belong to some other
%   session.
p = {}; fenced = false;
if nargin < 4, ctx = struct(); end

% entry step consuming the raw recording (contrast / internalCycle).  NEVER fenced:
% the raw recording is the working set, not a product of a configuration.
if isempty(step.requires) && ~contains(step.inGlob,'.mat')
    p = rawEntryFile(model, step);
    return
end

tail = regexprep(step.inGlob,'^\*','');                  % '_K_d.mat' | '_BFI_d.mat' | '_c_BFI_d.mat'
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.folder,[base '*' tail]));
if isempty(d), return; end

want = desiredStage(step, cstage);                       % 't'|'s'|'c' or '' (any)
exact = {}; rest = {}; mine = 0;
for i = 1:numel(d)
    fp = fullfile(d(i).folder, d(i).name);
    cm = wbFileModel(fp);
    if ~strcmp(cm.identity, model.identity), continue; end
    mine = mine + 1;
    if ~admitted(ctx, model.identity, fp), continue; end  % not part of this session
    if ~isempty(want) && strcmp(cm.stage, want)
        exact{end+1} = fp; %#ok<AGROW>
    else
        rest{end+1} = fp; %#ok<AGROW>
    end
end
p = [exact, rest];                                       % preferred branch first, then the others
fenced = isempty(p) && mine > 0;
end

% =====================================================================
function p = rawEntryFile(model, step)
%rawEntryFile  THE RAW RECORDING an entry step reads.  The step's inGlob names ONE
%   extension ('*.rls', '*.avi'), but a container is a family - a video recording
%   the user scanned may legitimately be a .mp4 - so the glob is the FALL-BACK and
%   the model's own path is the answer whenever the model IS the raw recording.
%   Which entry step a raw file may drive is settled long before this, by the
%   modality gate (wbStateEngine), so nothing is lost by not re-testing it here;
%   deriving the name from the glob instead sent the workbench looking for a file
%   that was never written.
p = {};
if isfield(model,'isRaw') && model.isRaw && ~isempty(model.path) && isfile(model.path)
    p = {model.path};  return
end
[~,~,ext] = fileparts(step.inGlob);
cand = fullfile(model.folder,[model.stem ext]);
if isfile(cand), p = {cand}; end
end

% =====================================================================
function tf = admitted(ctx, identity, path)
%admitted  THE FENCE (item 8).  The host answers it from the configuration; this
%   module only asks.  A ctx without the hook - every headless caller and every
%   test written before the fence existed - resolves exactly as it always did, and
%   a hook that throws is treated as no answer, because a fence that cannot say
%   what a session contains may never be the reason a file disappears from it.
tf = true;
if ~isstruct(ctx) || ~isfield(ctx,'admits') || ~isa(ctx.admits,'function_handle'), return; end
try
    tf = logical(ctx.admits(identity, path));
catch
    tf = true;
end
end

% =====================================================================
function r = missingReason(fenced)
%missingReason  Why a unit found no input, in the operator's words.  The two cases
%   have to read differently: being told a file is missing while it sits in the
%   folder is the one message that sends someone looking in the wrong place.
if fenced
    r = 'the file it needs is on disk but this session does not use it';
else
    r = 'input file not found';
end
end

% =====================================================================
function st = desiredStage(step, cstage)
%desiredStage  Which stage flag the step's input should carry - the branch it runs
%   on when a base-name glob matches several.  For a 'copy'-scope step this also
%   picks the branch the work is DONE on (the others inherit), and for 'all' it only
%   fixes the order.  '' = no preference.
switch step.branch
    case 'cardiac', st = 'c';
    case 'contrast', st = cstage;
    otherwise
        % A branch-agnostic column can still be anchored to the contrast side: the
        % 'copy' steps (setRegions, segmentation) draw/compute there and hand the
        % result to the cycle branches, and guided needs the contrast product.
        if strcmp(scopeOf(step),'copy') || strcmp(step.id,'guided')
            st = cstage;
        else
            st = '';
        end
end
end

% =====================================================================
function r = rawPathFor(model)
%rawPathFor  The raw recording beside a processed product (guided steps).  The
%   candidate extensions are the ones wbFileModel knows as raw containers, in its
%   order - the modality vocabulary lives there and this list must not become a
%   second copy of it that grows a modality late.
if isfield(model,'isRaw') && model.isRaw && ~isempty(model.path)
    r = model.path; return                               % it IS the raw recording
end
exts = wbFileModel('extensions');
for i = 1:numel(exts)
    cand = fullfile(model.folder,[model.stem exts{i}]);
    if isfile(cand), r = cand; return; end
end
r = fullfile(model.folder,[model.stem '.rls']);          % derived default (may not exist yet)
end

% =====================================================================
function sc = scopeOf(step)
%scopeOf  A step's branchScope, defaulting to 'one' when the field is absent - a
%   hand-built step struct (tests, a trimmed registry) keeps the old single-file
%   behaviour rather than erroring.
sc = 'one';
if isfield(step,'branchScope') && ~isempty(step.branchScope), sc = step.branchScope; end
end

% =====================================================================
function fo = fanOutOf(step)
%fanOutOf  A step's fanOut, defaulting to 'flat' when the field is absent (same
%   tolerance as scopeOf: a hand-built step struct keeps the column shape).
fo = 'flat';
if isfield(step,'fanOut') && ~isempty(step.fanOut), fo = step.fanOut; end
end

% =====================================================================
function fo = orderOf(step)
%orderOf  A step's fileOrder, defaulting to 'independent' when the field is absent.
%   Same tolerance again - a hand-built step struct (tests, a trimmed registry) then
%   gets the per-recording invocation, which is the safer of the two: it can only
%   narrow which cell an error reddens, never widen it.
fo = 'independent';
if isfield(step,'fileOrder') && ~isempty(step.fileOrder), fo = step.fileOrder; end
end

% =====================================================================
function l = batchLabel(units)
%batchLabel  How a call names itself in the run log: the animal for a per-animal
%   call, the recording for a single one, and a count once several were folded -
%   the per-file detail then comes from the wrapper's own lines.
e = units{1}.entry;
if strcmp(e.arity,'perAnimal')
    l = ['ANIMAL ' e.animal]; return
end
if isscalar(units), l = e.label; return; end
l = sprintf('%d recordings', numel(units));
end

% =====================================================================
function w = whoOf(e)
%whoOf  The name one entry answers to in a log line.
if strcmp(e.arity,'perAnimal'), w = ['ANIMAL ' e.animal]; else, w = e.label; end
end

% =====================================================================
function s = stepById(reg, id)
s = reg(strcmp({reg.id}, id));
if ~isempty(s), s = s(1); end
end

% =====================================================================
function b = emptyBatch()
%emptyBatch  A 0x0 batch struct carrying every field, so a run that produced no
%   call still returns something the caller can ask numel() and fieldnames() of.
b = newBatch(); b(1) = [];
end

% =====================================================================
function b = newBatch()
%newBatch  One batch with every field defaulted, so building them in any order
%   never triggers a field mismatch when they are concatenated.  The four file
%   fields stay empty until 'resolve' fills them.
b = struct('stepId','', 'step',[], 'arity','perFile', 'fanOut','flat', ...
    'fileOrder','independent', 'label','', ...
    'entries',[], 'heads',[], 'models',[], 'units',{{}}, 's',struct(), ...
    'fNames',{{}}, 'rawNames',{{}}, 'copyTo',{{}}, 'refName','');
end

% =====================================================================
function s = emptySkip()
%emptySkip  A 0x0 skip struct with every field.
s = struct('entry',{},'stepId',{},'stepLabel',{},'who',{},'reason',{},'mark',{});
end

% =====================================================================
function s = mkSkip(entry, stepId, stepLabel, who, reason, mark)
s = struct('entry',entry,'stepId',stepId,'stepLabel',stepLabel,'who',who, ...
    'reason',reason,'mark',logical(mark));
end
