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
%   with its siblings routed into the copy list.  A recording whose input is not on
%   disk yet cannot be part of a call, so it is reported back through the second
%   output and the rest of the batch runs without it.
%
%   Pure planning: it reads the registry, the models and the disk, and it changes
%   no state, marks no cell and calls no wrapper.
%
% Syntax:
%    batches            = wbBatchPlan('build', entries, ctx)
%    [batches, skipped] = wbBatchPlan('build', entries, ctx)
%
% Inputs:
%    entries - struct array from guiWorkbench>buildRunOrder, already step-major
%              (registry step -> animal -> reference-first -> file).  Each entry
%              has at least: stepId, identity, animalIdx, arity, label, animal.
%    ctx     - the wbExecutor context struct.  This module uses .reg, .modelOf,
%              .animalModels, .resolve, .contrastStage and the OPTIONAL .refFile
%              (see wbExecutor for what each of them promises).
%
% Outputs:
%    batches - 1xN struct array, one element per WRAPPER CALL, in run order:
%       stepId   the step's id
%       step     the registry step struct (carries the wrapper handle)
%       arity    'perFile' | 'perAnimal'
%       fanOut   'flat' | 'animal' - the shape below
%       label    how the call names itself in the log
%       entries  the folded entries, in order (afterDone is owed one call each)
%       heads    one model per ENTRY: the recording a perFile entry is, or the
%                reference a perAnimal entry starts at
%       models   every recording model the call touches, in order (the cells to
%                mark running / done / error)
%       s        the shared resolved settings, hooks NOT yet injected
%       fNames   the wrapper's file list, shaped per fanOut
%       rawNames the parallel raw recordings for a needsRaw step, else {}
%       copyTo   the s.fNamesCopyTo list, shaped per fanOut, {} when nothing is
%                inherited
%       refName  the per-animal reference file ('' when the call has none)
%    skipped - 1xM struct array of entries that could not become part of a call:
%       entry, stepId, stepLabel, who, reason, mark (true = the caller should mark
%       the cell as an error, which is what a missing input has always done).
%
% Notes:
%    * The entry list is authoritative for ORDER.  Nothing here sorts, and folding
%      only ever merges neighbours, so the launcher property the run order exists
%      to give - every recording reaches a level before anything moves past it -
%      survives untouched.
%    * A batch's settings are the FIRST folded entry's.  That is safe precisely
%      because the fold key is the settings themselves; the discarded copies were
%      equal.
%    * An error inside a call skips the remainder of that call, exactly as it does
%      inside a launcher cell (spec D1, accepted).
%
% See also: wbExecutor, wbStepRegistry, wbRunRange, wbRefBranch, wbFileModel,
%           guiWorkbench, setRegions, runSegmentation
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function [batches, skipped] = wbBatchPlan(action, entries, ctx)

switch lower(char(action))
    case 'build'
        [batches, skipped] = buildBatches(entries, ctx);
    otherwise
        error('wbBatchPlan:action','Unknown action ''%s''.', char(action));
end
end

% =====================================================================
function [batches, skipped] = buildBatches(entries, ctx)
%buildBatches  Resolve every entry to its files, then fold and shape the runs.
batches = emptyBatch(); skipped = emptySkip();
if isempty(entries), return; end

% ---- 1. one UNIT per entry: its models, its settings, its concrete files ------
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
    mHead  = models(1);
    s      = ctx.resolve(step, mHead);
    cstage = ctx.contrastStage(mHead);
    [files, copies, raws, ok] = entryFiles(step, models, cstage, refFileFor(ctx, e, step));
    if ~ok
        sk{end+1} = mkSkip(e, step.id, step.label, whoOf(e), 'input file not found', true); %#ok<AGROW>
        continue
    end
    U{end+1} = struct('entry',e, 'step',step, 's',s, 'head',mHead, 'models',models, ...
        'files',{files}, 'copies',{copies}, 'raws',{raws}, 'key',foldKey(s)); %#ok<AGROW>
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
for f = {'fName','fNamesCopyTo','refFName','progressFcn','stageFcn','cancelFcn'}
    if isfield(k,f{1}), k = rmfield(k,f{1}); end
end
if numel(fieldnames(k)) > 1, k = orderfields(k); end
end

% =====================================================================
function b = shapeBatch(units)
%shapeBatch  One run of folded units becomes ONE call, laid out per step.fanOut.
step = units{1}.step;
b        = newBatch();
b.stepId = step.id;
b.step   = step;
b.arity  = units{1}.entry.arity;
b.fanOut = fanOutOf(step);
b.s      = units{1}.s;

ents = cell(1,numel(units)); heads = cell(1,numel(units)); mdls = cell(1,numel(units));
for k = 1:numel(units)
    ents{k}  = units{k}.entry;
    heads{k} = units{k}.head;
    mdls{k}  = units{k}.models;
end
b.entries = [ents{:}];
b.heads   = [heads{:}];
b.models  = [mdls{:}];
b.label   = batchLabel(units);

switch b.fanOut
    case 'animal'
        [b.fNames, b.copyTo, b.rawNames] = shapeByAnimal(units, step);
    otherwise
        [b.fNames, b.copyTo, b.rawNames] = shapeFlat(units, step);
end

b.refName = '';
if strcmp(b.arity,'perAnimal') && ~isempty(b.fNames), b.refName = b.fNames{1}; end
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
function [files, copies, raws, ok] = entryFiles(step, models, cstage, refPath)
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
files = {}; copies = {}; raws = {}; ok = false;
if nargin < 4, refPath = ''; end

if strcmp(step.arity,'perAnimal')
    wide = ~strcmp(scopeOf(step),'one');
    for k = 1:numel(models)
        p = resolveStepInputs(models(k), step, cstage);       % stage-preferred first
        if k==1 && ~isempty(refPath) && isfile(refPath)
            p = [{refPath}, p(~strcmp(p,refPath))];           % the pinned reference leads
        end
        if isempty(p), return; end                            % a member has no input
        if ~wide, p = p(1); end
        files  = [files, p];                                  %#ok<AGROW>
        copies = [copies, repmat({{}},1,numel(p))];           %#ok<AGROW>
        raws   = [raws, repmat({rawPathFor(models(k))},1,numel(p))]; %#ok<AGROW>
    end
    ok = true; return
end

p = resolveStepInputs(models(1), step, cstage);
if isempty(p), return; end
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
function p = resolveStepInputs(model, step, cstage)
%resolveStepInputs  EVERY concrete _d.mat (or raw) file this step could consume for
%   this recording, ORDERED with the stage-preferred branch first.  Located by base
%   name + the step's input glob and filtered to this recording's identity; mirrors
%   getFileNamesList scoped to one recording, so it is robust to the BFI rename and
%   the t|s|c|e branches.  Callers take p{1} for a single-branch step and the whole
%   list (or its tail) for the 'all' / 'copy' scopes - see entryFiles.  {} = the
%   step has no input on disk yet.
p = {};

% entry step consuming the raw recording (contrast / internalCycle)
if isempty(step.requires) && ~contains(step.inGlob,'.mat')
    [~,~,ext] = fileparts(step.inGlob);
    cand = fullfile(model.folder,[model.stem ext]);
    if isfile(cand), p = {cand}; end
    return
end

tail = regexprep(step.inGlob,'^\*','');                  % '_K_d.mat' | '_BFI_d.mat' | '_c_BFI_d.mat'
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.folder,[base '*' tail]));
if isempty(d), return; end

want = desiredStage(step, cstage);                       % 't'|'s'|'c' or '' (any)
exact = {}; rest = {};
for i = 1:numel(d)
    fp = fullfile(d(i).folder, d(i).name);
    cm = wbFileModel(fp);
    if ~strcmp(cm.identity, model.identity), continue; end
    if ~isempty(want) && strcmp(cm.stage, want)
        exact{end+1} = fp; %#ok<AGROW>
    else
        rest{end+1} = fp; %#ok<AGROW>
    end
end
p = [exact, rest];                                       % preferred branch first, then the others
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
%rawPathFor  The raw recording beside a processed product (guided steps).
for ext = {'.rls','.cxd'}
    cand = fullfile(model.folder,[model.stem ext{1}]);
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
%   never triggers a field mismatch when they are concatenated.
b = struct('stepId','', 'step',[], 'arity','perFile', 'fanOut','flat', 'label','', ...
    'entries',[], 'heads',[], 'models',[], 's',struct(), ...
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
