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
%        for this recording at this stage - perFile = one file, perGroup = the
%        group's files as a column, reference first (matching the launchers),
%     4. injects the progress/stage/cancel hooks bound to this cell,
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
%              group -> reference-first -> row -> registry-step-order.  Each has
%              at least: stepId, identity, groupIdx, arity, label, group.
%    ctx     - struct of callbacks / handles the host supplies:
%       .reg              the wbStepRegistry array
%       .modelOf(id)      -> the wbFileModel for a recording identity
%       .groupModels(gi)  -> ordered model array for a group (reference first)
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
    if strcmp(e.arity,'perGroup')
        models = ctx.groupModels(e.groupIdx);
        who    = ['GROUP ' e.group];
    else
        models = ctx.modelOf(e.identity);
        who    = e.label;
    end
    if isempty(models)
        ctx.log(['skip: no recording for ' who ' :: ' step.label]); continue;
    end
    mHead = models(1);                                   % perFile: the file; perGroup: the reference
    s     = ctx.resolve(step, mHead);
    cstage = ctx.contrastStage(mHead);

    [fNames, rawNames, refName, okInput] = buildFNames(step, models, cstage);
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
        s.refFName = refName;                            % per-group paint reference (launcher idiom)
    end

    % ---- run it (a perGroup step marks EVERY group member's cell) -------------
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
function [fNames, rawNames, refName, ok] = buildFNames(step, models, cstage)
%buildFNames  Concrete wrapper inputs for a step: perFile = {file}; perGroup =
%   {files} column, reference first.  rawNames parallels fNames for needsRaw
%   steps; refName is the reference file (perGroup vesselTypes paint target).
fNames = {}; rawNames = {}; refName = ''; ok = false;
n = numel(models);
paths = cell(n,1); raws = cell(n,1);
for k = 1:n
    p = resolveStepInput(models(k), step, cstage);
    if isempty(p), return; end                           % a member has no input -> not ready
    paths{k} = p;
    if step.needsRaw, raws{k} = rawPathFor(models(k)); end
end
if strcmp(step.arity,'perGroup')
    fNames = paths;                                      % Nx1 column, reference first
    rawNames = raws;
    refName = paths{1};
else
    fNames = paths(1);                                   % 1x1 - workbench owns the per-file loop
    rawNames = raws(1);
end
ok = true;
end

% =====================================================================
function p = resolveStepInput(model, step, cstage)
%resolveStepInput  The concrete _d.mat (or raw) file this step consumes for this
%   recording, located by base name + the step's input glob, disambiguated by the
%   branch stage when several products coexist.  Mirrors getFileNamesList scoped
%   to one recording, so it is robust to the BFI rename and the t|s|c branches.
p = '';

% entry step consuming the raw recording (contrast / internalCycle)
if isempty(step.requires) && ~contains(step.inGlob,'.mat')
    [~,~,ext] = fileparts(step.inGlob);
    cand = fullfile(model.folder,[model.stem ext]);
    if isfile(cand), p = cand; end
    return
end

tail = regexprep(step.inGlob,'^\*','');                  % '_K_d.mat' | '_BFI_d.mat' | '_c_BFI_d.mat'
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.folder,[base '*' tail]));
if isempty(d), return; end

want = desiredStage(step, cstage);                       % 't'|'s'|'c' or '' (any)
cands = {}; exact = {};
for i = 1:numel(d)
    fp = fullfile(d(i).folder, d(i).name);
    cm = wbFileModel(fp);
    if ~strcmp(cm.identity, model.identity), continue; end
    cands{end+1} = fp; %#ok<AGROW>
    if ~isempty(want) && strcmp(cm.stage, want), exact{end+1} = fp; end %#ok<AGROW>
end
if ~isempty(exact)
    p = exact{1};
elseif ~isempty(cands)
    p = cands{1};                                        % branch-agnostic step: first match
end
end

% =====================================================================
function st = desiredStage(step, cstage)
%desiredStage  Which stage flag the step's input should carry, to disambiguate a
%   base-name glob that matches more than one branch.  '' = no preference.
switch step.branch
    case 'cardiac', st = 'c';
    case 'contrast', st = cstage;
    otherwise
        % segmentation/guided are computed on the contrast side even though their
        % column is branch-agnostic; prefer the contrast stage for them.
        if any(strcmp(step.id,{'segmentation','guided'})), st = cstage; else, st = ''; end
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
%   (one file for perFile; the whole group for perGroup), so the matrix reflects
%   the group together.
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
function l = hookLabel(l)
%hookLabel  Coerce a hook detail argument to a short char for the log.
if isnumeric(l), l = num2str(l); elseif isstring(l), l = char(l); elseif ~ischar(l), l = ''; end
end
