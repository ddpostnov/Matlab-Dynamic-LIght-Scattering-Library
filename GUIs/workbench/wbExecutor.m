%wbExecutor - Serial, ordered run loop for the Processing Workbench.
%
%   Executes an ordered list of checked (file,step) cells by calling the REAL
%   pipeline wrappers - it ORCHESTRATES, it never reimplements any science.  It
%   is the guiMyograph runList pattern generalised to the file x step matrix:
%   serial, on the main thread, streaming progress through the Phase-1 hook seam
%   (s.progressFcn / s.stageFcn / s.cancelFcn), with continue-on-error and
%   cooperative cancel.
%
%   For each entry it:
%     1. looks up the step spec and the recording model,
%     2. resolves the settings struct s (wbSettingsModel, via ctx.resolve),
%     3. resolves the concrete input _d.mat (or raw) file(s) the wrapper consumes
%        for this recording - perAnimal = the animal's files as a column, reference
%        first; perFile = as many of the recording's co-registered BRANCH products
%        ('_t_K', '_c_K', '_e_K', ...) as the step's branchScope asks for, so a
%        workbench row reproduces what the launcher's file list covers,
%     4. injects the progress/stage/cancel hooks bound to this cell, and - for a
%        'copy'-scope step - the derived s.fNamesCopyTo sibling list,
%     5. marks the cell running, calls the wrapper (interactive steps go through
%        ctx.modalGuard so the parent window is parked), then marks it done or
%        error and surfaces its report artifacts + invalidates downstream.
%
%   The executor is decoupled from the figure through a CONTEXT struct of
%   callbacks (ctx), so it can be unit-tested headlessly (drive it with recording
%   callbacks and assert the call order) as well as from the live uifigure.
%
% Syntax:
%    wbExecutor(entries, ctx)
%
% Inputs:
%    entries - struct array from guiWorkbench>buildRunOrder, already sorted
%              animal -> reference-first -> row -> registry-step-order.  Each has
%              at least: stepId, identity, animalIdx, arity, label, animal.
%    ctx     - struct of callbacks / handles the host supplies:
%       .reg              the wbStepRegistry array
%       .modelOf(id)      -> the wbFileModel for a recording identity
%       .animalModels(ai) -> ordered model array for an animal (reference first)
%       .resolve(step,mdl)-> resolved settings s (wbSettingsModel)
%       .contrastStage(m) -> 't' | 's'  (which contrast flag this project uses)
%       .setState(id,step,state,msg)   set a cell running/done/error (''=revert)
%       .progress(id,step,frac,label)  progress-hook sink
%       .log(msg)                      append one line to the log / CW mirror
%       .isCancelled()    -> tf        cooperative cancel flag
%       .modalGuard(fcn)               run fcn with the parent window parked
%       .afterDone(id,step,model)      surface artifacts + invalidate downstream
%
% Notes:
%    * Cancel is checked between cells and, via s.cancelFcn, between files inside
%      a wrapper (the Phase-1 seam only cancels on file boundaries).  A cancelled
%      cell is reverted (not marked done) and the batch stops cleanly.
%    * Any error in one cell is caught -> that cell goes 'error' and the batch
%      CONTINUES with the next cell (never aborts the whole run).
%    * Steps whose prerequisite input file cannot be located are logged and
%      skipped (left for a later run), not errored.
%    * ONE ROW = ONE RECORDING, several files.  A step with branchScope 'all'
%      (splitRegions, BFI, dynamicSegmentation, export) receives every branch of the
%      recording in a single call and its own per-file loop covers them; a step with
%      'copy' (setRegions, segmentation) runs on the contrast branch only and the
%      wrapper propagates the result to the rest through s.fNamesCopyTo.  Without
%      this the executor picked ONE file per recording and the first dir() match won
%      - which quietly ran the interactive region editor on '_c' and never touched
%      '_t' (dir sorts '_c_K_d.mat' before '_t_K_d.mat').
%
% See also: guiWorkbench, wbStepRegistry, wbSettingsModel, wbModalGuard,
%           wbArtifacts, wbInvalidate, guiMyograph
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function wbExecutor(entries, ctx)

if isempty(entries), ctx.log('Nothing checked - nothing to run.'); return; end
ctx.log(sprintf('=== RUN: %d step(s) ===', numel(entries)));

for i = 1:numel(entries)
    if ctx.isCancelled(), ctx.log('Stopped.'); break; end
    e    = entries(i);
    step = stepById(ctx.reg, e.stepId);
    if isempty(step), ctx.log(['skip: unknown step ' e.stepId]); continue; end

    % ---- resolve the recording model(s) + the concrete input file(s) ----------
    if strcmp(e.arity,'perAnimal')
        models = ctx.animalModels(e.animalIdx);
        who    = ['ANIMAL ' e.animal];
    else
        models = ctx.modelOf(e.identity);
        who    = e.label;
    end
    if isempty(models)
        ctx.log(['skip: no recording for ' who ' :: ' step.label]); continue;
    end
    mHead = models(1);                                   % perFile: the file; perAnimal: the reference
    s     = ctx.resolve(step, mHead);
    cstage = ctx.contrastStage(mHead);

    [fNames, rawNames, refName, copyTo, okInput] = buildFNames(step, models, cstage);
    if ~okInput
        ctx.setState(e.identity, step.id, 'error', 'input file not found');
        ctx.log(sprintf('skip %s :: %s - input file not found', who, step.label));
        continue
    end

    % ---- inject the hooks bound to THIS cell ----------------------------------
    id = e.identity; sid = step.id;
    s.progressFcn = @(f,l) ctx.progress(id, sid, f, hookLabel(l));
    s.stageFcn    = @(st,d) ctx.log(sprintf('  [%s] %s', char(st), hookLabel(d)));
    s.cancelFcn   = @() ctx.isCancelled();
    if strcmp(step.id,'vesselTypes') && ~isempty(refName)
        s.refFName = refName;                            % per-animal paint reference (launcher idiom)
    end
    if ~isempty(copyTo)
        s.fNamesCopyTo = copyTo;                         % 'copy' scope: the other branches inherit
        ctx.log(sprintf('  copy to %d sibling(s): %s', numel(copyTo), ...
            strjoin(cellfun(@shortName,copyTo,'UniformOutput',false),', ')));
    end

    % ---- run it (a perAnimal step marks EVERY of the animal's cells) ----------
    markCells(ctx, sid, models, 'running', '');
    ctx.log(sprintf('run %s :: %s', who, step.label));
    call = @() invokeWrapper(step, s, fNames, rawNames);
    try
        if isInteractiveStep(step, s)
            ctx.modalGuard(call);
        else
            call();
        end
    catch ME
        if isCancelErr(ME)
            markCells(ctx, sid, models, '', '');         % revert; may be partial
            ctx.log('Stopped by user.');
            break
        end
        markCells(ctx, sid, models, 'error', ME.message);
        ctx.log(sprintf('ERROR %s :: %s - %s', who, step.label, ME.message));
        continue                                         % continue-on-error
    end

    % cooperative cancel can fire mid-batch without throwing (checked between
    % files inside the wrapper): if it did, do NOT mark this cell done.
    if ctx.isCancelled()
        markCells(ctx, sid, models, '', '');
        ctx.log('Stopped.');
        break
    end

    markCells(ctx, sid, models, 'done', '');
    ctx.afterDone(mHead.identity, sid, mHead);
end

ctx.log('=== RUN complete ===');
end

% =====================================================================
function [fNames, rawNames, refName, copyTo, ok] = buildFNames(step, models, cstage)
%buildFNames  Concrete wrapper inputs for a step.
%   perAnimal: one file per animal member, Nx1 column, reference first (the branch
%   fan-out below is a perFile concern - an animal row already carries one file
%   per recording).
%   perFile:  ONE recording, whose several co-registered branch products ('_t_K',
%   '_c_K', '_e_K', ...) are resolved according to step.branchScope -
%     'one'  -> the stage-preferred file alone (1x1),
%     'all'  -> every branch as an Nx1 column, the wrapper's own loop covers them,
%     'copy' -> the stage-preferred file (1x1) plus the rest returned in copyTo, to
%               be handed to the wrapper as s.fNamesCopyTo.
%   rawNames parallels fNames for needsRaw steps; refName is the reference file
%   (perAnimal vesselTypes paint target).
fNames = {}; rawNames = {}; refName = ''; copyTo = {}; ok = false;

if strcmp(step.arity,'perAnimal')
    n = numel(models);
    paths = cell(n,1); raws = cell(n,1);
    for k = 1:n
        p = resolveStepInputs(models(k), step, cstage);
        if isempty(p), return; end                       % a member has no input -> not ready
        paths{k} = p{1};
        if step.needsRaw, raws{k} = rawPathFor(models(k)); end
    end
    fNames = paths; rawNames = raws; refName = paths{1};
    ok = true; return
end

p = resolveStepInputs(models(1), step, cstage);           % stage-preferred first
if isempty(p), return; end
switch scopeOf(step)
    case 'all'
        fNames = p(:);                                   % every branch, one cell for the lot
    case 'copy'
        fNames = p(1);                                   % run here...
        copyTo = p(2:end);                               % ...the other branches inherit
    otherwise
        fNames = p(1);
end
if step.needsRaw
    rawNames = repmat({rawPathFor(models(1))}, size(fNames));
end
ok = true;
end

% =====================================================================
function p = resolveStepInputs(model, step, cstage)
%resolveStepInputs  EVERY concrete _d.mat (or raw) file this step could consume for
%   this recording, ORDERED with the stage-preferred branch first.  Located by base
%   name + the step's input glob and filtered to this recording's identity; mirrors
%   getFileNamesList scoped to one recording, so it is robust to the BFI rename and
%   the t|s|c|e branches.  Callers take p{1} for a single-branch step and the whole
%   list (or its tail) for the 'all' / 'copy' scopes - see buildFNames.  {} = the
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
function markCells(ctx, sid, models, state, msg)
%markCells  Set the run-state of one step across every recording it touched
%   (one file for perFile; the whole animal for perAnimal), so the matrix reflects
%   the animal together.
for k = 1:numel(models)
    ctx.setState(models(k).identity, sid, state, msg);
end
end

% =====================================================================
function invokeWrapper(step, s, fNames, rawNames)
%invokeWrapper  Dispatch to the seam function with the right arity.
if step.needsRaw
    step.wrapper(s, fNames, rawNames);
else
    step.wrapper(s, fNames);
end
end

% =====================================================================
function tf = isInteractiveStep(step, s)
%isInteractiveStep  Resolve step.interactive (false | true | @(s)tf) for these settings.
it = step.interactive;
if islogical(it)
    tf = it;
elseif isa(it,'function_handle')
    try
        tf = logical(it(s));
    catch
        tf = false;
    end
else
    tf = false;
end
end

% =====================================================================
function tf = isCancelErr(ME)
%isCancelErr  Whether an error is really a user-cancel signalled by a wrapper.
tf = contains(lower(ME.identifier),'cancel') || ...
     contains(lower(ME.message),'stopped by user') || ...
     contains(lower(ME.message),'cancelled');
end

% =====================================================================
function s = stepById(reg, id)
s = reg(strcmp({reg.id}, id));
if ~isempty(s), s = s(1); end
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
function n = shortName(p)
%shortName  Bare file name of a path, for the copy-target log line.
[~,nm,ex] = fileparts(char(p));
n = [nm ex];
end

% =====================================================================
function l = hookLabel(l)
%hookLabel  Coerce a hook detail argument to a short char for the log.
if isnumeric(l), l = num2str(l); elseif isstring(l), l = char(l); elseif ~ischar(l), l = ''; end
end
