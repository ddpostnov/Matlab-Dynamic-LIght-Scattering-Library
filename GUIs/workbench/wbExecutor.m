%wbExecutor - Serial, ordered run loop for the Processing Workbench.
%
%   Executes the checked work of a run by calling the REAL pipeline wrappers - it
%   ORCHESTRATES, it never reimplements any science.  It is the guiMyograph runList
%   pattern generalised to the file x step matrix: serial, on the main thread,
%   streaming progress through the Phase-1 hook seam (s.progressFcn / s.stageFcn /
%   s.cancelFcn), with continue-on-error and cooperative cancel.
%
%   ONE CALL IS A BATCH, NOT A CELL.  Every wrapper in Wrappers/ is a multi-file
%   routine and the launchers call it with a whole getFileNamesList cell; anything
%   a wrapper does BETWEEN files is unreachable when it is handed one recording at
%   a time.  So the entry list is first translated into WRAPPER CALLS by
%   wbBatchPlan, which folds the entries of one step whose resolved settings agree
%   and shapes their files the way that step declares (registry fanOut).  This
%   module then does one thing per call:
%     1. mark every recording of the call 'running',
%     2. inject the progress/stage/cancel hooks and, when the step inherits its
%        result onto co-registered siblings, s.fNamesCopyTo,
%     3. invoke the wrapper once (an interactive step goes through ctx.modalGuard
%        so the parent window is parked),
%     4. mark every recording 'done', or every recording 'error' with the message,
%     5. call ctx.afterDone ONCE PER FOLDED ENTRY, so per-recording artifact
%        collection, downstream invalidation and the per-column PDF accounting are
%        exactly what they were when a call was a single cell (spec D6).
%
%   The executor is decoupled from the figure through a CONTEXT struct of
%   callbacks (ctx), so it can be unit-tested headlessly (drive it with recording
%   callbacks and assert the calls) as well as from the live uifigure.
%
% Syntax:
%    wbExecutor(entries, ctx)
%
% Inputs:
%    entries - struct array from guiWorkbench>buildRunOrder, already sorted
%              STEP-MAJOR: registry step -> animal -> reference-first -> file.
%              That is the launcher's order, where one cell hands a single step the
%              whole fNames list, so every recording reaches a level before
%              anything moves past it - which is what the cross-file steps
%              (registration, vessel typing, split regions) rely on.  The executor
%              does not re-sort, and batching only ever merges NEIGHBOURS, so the
%              list still runs as given.  Each entry has at least: stepId,
%              identity, animalIdx, arity, label, animal.
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
%       .refFile(ai,step) -> OPTIONAL.  Which FILE of an animal's pinned
%                         reference RECORDING this step takes as its template /
%                         paint target - the registry's refBranch resolved by
%                         wbRefBranch ('_t'/'_s' for registration, '_c' for
%                         vessel typing).  The host owns that resolution because
%                         it owns the pin; supply '' (or omit the field) and the
%                         reference falls back to whatever branch product of the
%                         reference recording the step's input glob finds first.
%
% Notes:
%    * WHICH FILES A CALL GETS IS NOT DECIDED HERE.  wbBatchPlan owns the whole of
%      it - the input glob, the branch fan-out (branchScope), the copy targets, the
%      raw partner, the per-animal reference, and the shape (fanOut: a flat column,
%      or one row per animal for a wrapper that carries state across a row).  This
%      module knows only that a batch has an fNames and a list of recordings.
%    * STATE IS REPORTED PER BATCH (spec D4).  Every recording of a call reads
%      'running' from the moment it starts until it returns, and then 'done' or
%      'error' together; the percent rides on the call's first recording, and the
%      per-file detail comes from the wrappers' own disp / stageFcn lines.  No
%      wrapper was modified to report progress.
%    * Cancel is checked between calls and, via s.cancelFcn, between files inside a
%      wrapper (the Phase-1 seam only cancels on file boundaries).  A cancel that
%      lands mid-call reverts THAT call's recordings - they may be partly processed
%      - and the run stops cleanly.
%    * Any error in one call is caught -> that call's recordings go 'error' and the
%      run CONTINUES with the next call.  The remainder of the failed call is lost,
%      exactly as it is inside a launcher cell (spec D1, accepted).
%    * A recording whose prerequisite input cannot be located never reaches a
%      wrapper: it is reported before the run starts and left for a later one.  It
%      does not take the rest of its call down with it.
%
% See also: wbBatchPlan, guiWorkbench, wbStepRegistry, wbSettingsModel,
%           wbModalGuard, wbArtifacts, wbInvalidate, wbRefBranch, guiMyograph
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function wbExecutor(entries, ctx)

if isempty(entries), ctx.log('Nothing checked - nothing to run.'); return; end

[batches, skipped] = wbBatchPlan('build', entries, ctx);
reportSkipped(ctx, skipped);
ctx.log(sprintf('=== RUN: %d job(s) in %d call(s) ===', numel(entries), numel(batches)));

for i = 1:numel(batches)
    if ctx.isCancelled(), ctx.log('Stopped.'); break; end
    b   = batches(i);
    sid = b.stepId;

    % ---- the hooks, bound to the CALL (its first recording carries the percent) ---
    s = b.s;
    s.progressFcn = @(f,l) ctx.progress(b.heads(1).identity, sid, f, hookLabel(l));
    s.stageFcn    = @(st,d) ctx.log(sprintf('  [%s] %s', char(st), hookLabel(d)));
    s.cancelFcn   = @() ctx.isCancelled();
    if strcmp(sid,'vesselTypes') && ~isempty(b.refName)
        s.refFName = b.refName;                          % per-animal paint reference (launcher idiom)
    end
    if ~isempty(b.copyTo)
        s.fNamesCopyTo = b.copyTo;                       % the other branches inherit the result
        tg = copyTargetList(b.copyTo);
        ctx.log(sprintf('  copy to %d sibling(s): %s', numel(tg), ...
            strjoin(cellfun(@shortName,tg,'UniformOutput',false),', ')));
    end

    % ---- run it (every recording of the call is marked together) ------------------
    markCells(ctx, sid, b.models, 'running', '');
    % the count is of REAL files: a per-animal shape pads its ragged rows with ''
    ctx.log(sprintf('run %s :: %s (%d file(s))', b.label, b.step.label, ...
        nnz(~cellfun(@isempty, b.fNames(:)))));
    call = @() invokeWrapper(b.step, s, b.fNames, b.rawNames);
    try
        if isInteractiveStep(b.step, s)
            ctx.modalGuard(call);
        else
            call();
        end
    catch ME
        if isCancelErr(ME)
            markCells(ctx, sid, b.models, '', '');       % revert; may be partial
            ctx.log('Stopped by user.');
            break
        end
        markCells(ctx, sid, b.models, 'error', ME.message);
        ctx.log(sprintf('ERROR %s :: %s - %s', b.label, b.step.label, ME.message));
        continue                                         % continue-on-error
    end

    % cooperative cancel can fire mid-call without throwing (checked between files
    % inside the wrapper): if it did, do NOT mark this call done.
    if ctx.isCancelled()
        markCells(ctx, sid, b.models, '', '');
        ctx.log('Stopped.');
        break
    end

    markCells(ctx, sid, b.models, 'done', '');
    for k = 1:numel(b.entries)                           % once per FOLDED ENTRY (spec D6)
        ctx.afterDone(b.heads(k).identity, sid, b.heads(k));
    end
end

ctx.log('=== RUN complete ===');
end

% =====================================================================
function reportSkipped(ctx, skipped)
%reportSkipped  Entries that never became part of a call.  A missing input has
%   always marked its own cell, so it still does - but it is settled BEFORE the run
%   starts, because a call now spans several recordings and one absent file may not
%   decide the fate of the others.
for i = 1:numel(skipped)
    sk = skipped(i);
    if sk.mark, ctx.setState(sk.entry.identity, sk.stepId, 'error', sk.reason); end
    ctx.log(sprintf('skip %s :: %s - %s', sk.who, sk.stepLabel, sk.reason));
end
end

% =====================================================================
function markCells(ctx, sid, models, state, msg)
%markCells  Set the run-state of one step across every recording the call touched,
%   so the matrix reflects the whole call together (spec D4).
for k = 1:numel(models)
    ctx.setState(models(k).identity, sid, state, msg);
end
end

% =====================================================================
function t = copyTargetList(copyTo)
%copyTargetList  The copy targets as one flat cellstr, for the log line only.  Both
%   accepted shapes are flattened here: the row-per-source-file list and the
%   elementwise one, whose entries may themselves be nested cellstr.
t = {};
for i = 1:numel(copyTo)
    e = copyTo{i};
    if isempty(e), continue; end
    if ischar(e) || isstring(e), t{end+1} = char(e); else, t = [t, e(:)']; end %#ok<AGROW>
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
