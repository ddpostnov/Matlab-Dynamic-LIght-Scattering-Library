%wbExecutor - Serial, ordered run loop for the Processing Workbench.
%
%   Executes the checked work of a run by calling the REAL pipeline wrappers - it
%   ORCHESTRATES, it never reimplements any science.  It is a run-list loop
%   generalised to the file x step matrix: serial, on the main thread,
%   streaming the wrappers' own narration through the hook seam (s.stageFcn /
%   s.cancelFcn), with continue-on-error and cooperative cancel.
%
%   ONE CALL IS A BATCH, NOT A CELL.  Every wrapper in Wrappers/ is a multi-file
%   routine and the launchers call it with a whole getFileNamesList cell; anything
%   a wrapper does BETWEEN files is unreachable when it is handed one recording at
%   a time.  So the entry list is first translated into WRAPPER CALLS by
%   wbBatchPlan, which folds the entries of one step whose resolved settings agree
%   and shapes their files the way that step declares (registry fanOut).  This
%   module then does one thing per call:
%     1. RESOLVE the call's files, against the disk as it is at that moment,
%     2. mark the recordings about to run 'running',
%     3. inject the stage/cancel hooks and, when the step inherits its result onto
%        co-registered siblings, s.fNamesCopyTo,
%     4. invoke the wrapper (an interactive step goes through ctx.modalGuard so the
%        parent window is parked),
%     5. mark 'done', or 'error' with the message,
%     6. call ctx.afterDone ONCE PER FOLDED ENTRY, so per-recording artifact
%        collection, downstream invalidation and the per-column PDF accounting are
%        exactly what they were when a call was a single cell (spec D6).
%
%   THE WHOLE SLICE IS MARKED 'queued' BEFORE THE FIRST CALL (spec D2).  Pressing
%   Run used to change nothing on screen until the first wrapper had finished a
%   file, so a step queued for RE-PROCESSING went on reading 'done' from disk - the
%   monitor was describing the last run rather than this one.  The plan is known in
%   full the moment wbBatchPlan has built it, so it is stated then: one setState per
%   batch over that batch's recordings, before the loop.  Nothing has to undo it -
%   the host drops every non-error overlay when the run ends, so a batch that was
%   never reached (a Stop, an error upstream) reverts to its disk state on its own.
%
%   THE FILES ARE RESOLVED HERE, NOT BEFORE THE RUN (author, 2026-07-31).  The run
%   order is step-major, so contrast runs before the setRegions entry that consumes
%   its product - but resolving every batch's files in one pass up front asked for
%   that product while it did not exist yet, and setRegions was marked as an error
%   before anything had been attempted.  Step 1 above is that fix: wbBatchPlan
%   'build' now touches no disk, and 'resolve' runs immediately before each call.
%   A batch whose input is genuinely still missing is marked and logged THERE, one
%   line, and the loop carries on with the next call.
%
%   A RECOMPUTED SOURCE INVALIDATES ITS DERIVATIVES (round-2 item 9, spec D9).
%   MATLAB's save() already replaces a file whole - there is no '-append' anywhere
%   in the library and every producer clears its variables per file - so a step
%   that re-runs overwrites its OWN triplet.  What used to survive is everything
%   BELOW it: re-running contrast wrote a fresh '_t_K' while the '_t_BFI' triplet
%   computed from the old one, the external cycle's, and their report pages stayed
%   on disk, went on reading 'done', were pruned out of the next run and were
%   picked up by the reports.  So immediately before a step that writes a NEW
%   triplet is invoked, the products wbProducts derives as downstream of it are
%   removed and NAMED in the log.  Deliberately PER INVOCATION, not per run: a step
%   the run never reaches must never have deleted anything.  A delete that fails is
%   logged and the run carries on - this is housekeeping, it may not take a run
%   down.
%
%   ONE FILE ERRORS, ONE CELL GOES RED (spec D8/D9).  Whether a fold may be cut up
%   is the step's own declaration, registry fileOrder:
%     'independent' - the wrapper simply loops its list, so the executor invokes it
%                     ONCE PER RECORDING out of the fold.  A throw then reddens that
%                     recording alone; the ones already finished stay done and keep
%                     their afterDone, the ones after it still run, and the ones
%                     never reached are left untouched rather than being marked.
%     'ordered'     - the wrapper reads ACROSS its list (a template, a carry-forward,
%                     a per-file stimulus index), so it gets the single call it
%                     always got and a throw marks the whole call.  That is correct
%                     rather than regrettable: the cross-file state is gone with it.
%   The try/catch lives HERE, once (spec D9).  No wrapper has one and no wrapper was
%   asked to change.
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
%       .setState(id,step,state,msg)   set a cell queued/running/done/error
%                         (''=revert).  'id' is ONE identity or a CELLSTR of
%                         them: a call marks every recording it touches in one
%                         go, so the pre-run pass below costs one repaint per
%                         call rather than one per recording (spec D3).
%       .log(msg)                      append one line to the log / CW mirror
%       .isCancelled()    -> tf        cooperative cancel flag
%       .modalGuard(fcn)               run fcn with the parent window parked
%       .afterDone(id,step,model)      surface artifacts + invalidate downstream
%       .admits(id,path)  -> OPTIONAL.  Is this file on disk part of THIS
%                         session?  The host answers it from the configuration
%                         (guiWorkbench > wbProducts 'flags'), so the files a
%                         wrapper is handed are the ones the derived monitor table
%                         lists and nothing else.  Consumed by wbBatchPlan, not
%                         here.  Omit the field and nothing is fenced.
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
%      module knows only that a batch has an fNames and a list of recordings, and
%      that each of its units carries the same thing for one recording.
%    * STATE IS REPORTED PER INVOCATION.  For an 'ordered' step that is the whole
%      call, exactly as before (spec D4); for an 'independent' one it is the
%      recording, which is what makes a single red cell possible.  'running' carries
%      NO percentage - there is no progress axis - and the per-file detail comes
%      from the wrapper, which narrates through Core/Reporting.
%    * THE LOG IS A LIST OF RECORDINGS, NOT A TRANSCRIPT.  A wrapper emits exactly
%      three lines per recording - Starting, Writing results, Finished with the
%      time it took - and every one of them reaches s.stageFcn and so this log.  A
%      60-file column is therefore 3 lines per recording and nothing else, which is
%      what keeps it inside wbLog's 400-line cap.
%    * Cancel is checked between calls and, via s.cancelFcn, between files inside a
%      wrapper (the Phase-1 seam only cancels on file boundaries).  A cancel that
%      lands mid-call reverts THAT invocation's recordings - they may be partly
%      processed - and the run stops cleanly.  Cancel keeps precedence over
%      everything below.
%    * WHERE A FAILED RUN STOPS (spec D7).  An error never aborts a step: the step
%      is attempted on every one of its files whatever happened to the others.  The
%      run then carries on until the expansion reaches the first step whose
%      fileOrder is 'ordered' while a recording is still red, and stops cleanly
%      there - a step that reads across the whole file set must not be handed a set
%      with holes in it.  The log line names what failed and what was not attempted.
%      With nothing red the run goes to the end, as it always did.
%
% See also: wbBatchPlan, guiWorkbench, wbStepRegistry, wbSettingsModel,
%           wbModalGuard, wbArtifacts, wbInvalidate, wbRefBranch
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function wbExecutor(entries, ctx)

if isempty(entries), ctx.log('Nothing checked - nothing to run.'); return; end

[batches, skipped] = wbBatchPlan('build', entries, ctx);
reportSkipped(ctx, skipped);
ctx.log(sprintf('=== RUN: %d job(s) in %d call(s) ===', numel(entries), numel(batches)));
markQueued(ctx, batches);

failed  = {};                    % what has gone red so far, in words (spec D7)
stopped = false;
for i = 1:numel(batches)
    if ctx.isCancelled(), ctx.log('Stopped.'); stopped = true; break; end
    b   = batches(i);
    sid = b.stepId;

    % ---- a cross-file step must not meet a set with holes in it (spec D7) --------
    if ~isempty(failed) && isOrdered(b)
        ctx.log(stopLine(b, batches(i:end), failed));
        stopped = true;
        break
    end

    % ---- WHICH FILES, asked now that the earlier steps have written theirs -------
    [b, okFiles, whyFiles] = wbBatchPlan('resolve', b, ctx);
    if ~okFiles
        markCells(ctx, sid, b.models, 'error', whyFiles);
        ctx.log(sprintf('skip %s - %s: %s', b.step.label, b.label, whyFiles));
        failed{end+1} = sprintf('%s / %s', b.step.label, b.label); %#ok<AGROW>
        continue
    end

    % ---- the hooks, bound to the CALL -------------------------------------------
    %
    % NO TAG ON THE LOG LINE.  A wrapper emits exactly three lines per recording
    % (Starting / Writing results / Finished ... time elapsed), each of them a
    % finished sentence that names the recording, so there is nothing to filter and
    % nothing to label.  The executor knows which step is running without being told
    % - the tag it used to prepend was the WRAPPER's name, which is not a word the
    % operator has ever seen.  The step is named once, in the call header below, and
    % everything the call says is indented under it.
    %
    % NO PROGRESS HOOK.  There is no progress axis any more: a long loop is silent
    % between its two lines, and the monitor cell reads 'running' with no percentage.
    % NO PDF SWITCHES either - a wrapper never assembles a document, so there is
    % nothing to opt out of; flushPdfColumn assembles the column's own, between the
    % steps, out of the images wbArtifacts re-resolves from disk.
    %
    % s.fNamesCopyTo and s.refFName are NOT set here: they belong to one invocation,
    % and an independent step has one per recording.  callSettings adds them.
    s = b.s;
    s.stageFcn  = @(~,d) ctx.log(['  ' hookLabel(d)]);
    s.cancelFcn = @() ctx.isCancelled();

    % ONE HEADER PER CALL, and the wrapper's own per-file banners are the indented
    % lines under it.  The step comes first because that is what the column is; the
    % count is of REAL files (a per-animal shape pads its ragged rows with '').
    ctx.log(sprintf('%s - %s (%d file(s))', b.step.label, b.label, fileCount(b.fNames)));
    logCopyTargets(ctx, b.copyTo);

    if isOrdered(b)
        [outcome, whoFailed] = runBatchWhole(ctx, b, s);
    else
        [outcome, whoFailed] = runBatchPerUnit(ctx, b, s);
    end
    if strcmp(outcome,'cancel'), stopped = true; break; end
    failed = [failed, whoFailed]; %#ok<AGROW>   already marked and logged where it happened
end

if stopped
    ctx.log('=== RUN stopped ===');
else
    ctx.log('=== RUN complete ===');
end
end

% =====================================================================
function [outcome, failedNames] = runBatchWhole(ctx, b, s)
%runBatchWhole  An 'ordered' step: ONE call over the whole fold, marked together.
%   A throw takes the call with it, which is the honest reading - whatever the
%   wrapper was carrying across its files is gone.
outcome = 'ok'; failedNames = {};
sid = b.stepId;
markCells(ctx, sid, b.models, 'running', '');
pruneSuperseded(ctx, b.step, b.models);
[ok, ME] = invokeGuarded(ctx, b.step, callSettings(s, sid, b), b.fNames, b.rawNames);
if ~ok
    if isCancelErr(ME)
        markCells(ctx, sid, b.models, '', '');           % revert; may be partial
        ctx.log('Stopped by user.');
        outcome = 'cancel'; return
    end
    markCells(ctx, sid, b.models, 'error', ME.message);
    ctx.log(errorLine(b.step, b.label, ME));
    outcome = 'error'; failedNames = {sprintf('%s / %s', b.step.label, b.label)};
    return
end
% cooperative cancel can fire mid-call without throwing (checked between files
% inside the wrapper): if it did, do NOT mark this call done.
if ctx.isCancelled()
    markCells(ctx, sid, b.models, '', '');
    ctx.log('Stopped.');
    outcome = 'cancel'; return
end
markCells(ctx, sid, b.models, 'done', '');
finishUnits(ctx, sid, b.units);
end

% =====================================================================
function [outcome, failedNames] = runBatchPerUnit(ctx, b, s)
%runBatchPerUnit  An 'independent' step: one invocation per RECORDING, out of the
%   same fold.  The recordings are marked as they go, so a failure reddens one cell,
%   the ones before it keep their 'done' and their afterDone, the ones after it
%   still run, and the ones a cancel never reached carry no state at all.
outcome = 'ok'; failedNames = {};
sid = b.stepId;
for k = 1:numel(b.units)
    u = b.units{k};
    markCells(ctx, sid, u.models, 'running', '');
    pruneSuperseded(ctx, b.step, u.models);
    [ok, ME] = invokeGuarded(ctx, b.step, callSettings(s, sid, u), u.fNames, u.rawNames);
    if ~ok
        if isCancelErr(ME)
            markCells(ctx, sid, u.models, '', '');       % revert; may be partial
            ctx.log('Stopped by user.');
            outcome = 'cancel'; return
        end
        markCells(ctx, sid, u.models, 'error', ME.message);
        ctx.log(errorLine(b.step, unitName(u), ME));
        outcome = 'error';
        failedNames{end+1} = sprintf('%s / %s', b.step.label, unitName(u)); %#ok<AGROW>
        continue                                          % the rest of the step still runs
    end
    if ctx.isCancelled()
        markCells(ctx, sid, u.models, '', '');
        ctx.log('Stopped.');
        outcome = 'cancel'; return
    end
    markCells(ctx, sid, u.models, 'done', '');
    finishUnits(ctx, sid, b.units(k));
end
end

% =====================================================================
function pruneSuperseded(ctx, step, models)
%pruneSuperseded  A producer re-run cleans up after itself (item 9, spec D9).
%   Only a step that writes a NEW triplet has anything below it; an in-place step
%   appends to a product that is still the same product, so nothing it could touch
%   goes stale.  The STAGE being rewritten is the step's own (wbProducts 'writes',
%   which asks the settings for the contrast producer's t|s), and wbProducts owns
%   the whole rule about which files that makes superseded - this loop only carries
%   out the verdict and says what it did.
if isempty(models) || ~isstruct(step), return; end
if ~isfield(step,'outKind') || ~strcmp(step.outKind,'new'), return; end
for k = 1:numel(models)
    m = models(k);
    try
        stage = wbProducts('writes', step, ctx.contrastStage(m));
        [gone, unsure] = wbProducts('below', ctx.reg, m, step.id, stage);
    catch ME
        ctx.log(sprintf('  could not check for superseded results: %s', ME.message));
        continue
    end
    removeSuperseded(ctx, gone);
    for q = 1:numel(unsure)
        % D9a: it matches a downstream product but its name does not say which
        % result it was derived from, so it is left alone and named rather than
        % deleted on a guess.
        ctx.log(sprintf('  kept %s - its name does not say which result it came from', ...
            shortName(unsure{q})));
    end
end
end

% =====================================================================
function removeSuperseded(ctx, files)
%removeSuperseded  Delete, and NAME what was deleted.  A silent deletion is the one
%   thing a destructive housekeeping pass must not be; a failed one is a log line
%   and never the end of a run.
if isempty(files), return; end
gone = cell(1,0);
for i = 1:numel(files)
    try
        delete(files{i});
        gone{end+1} = shortName(files{i}); %#ok<AGROW>
    catch ME
        ctx.log(sprintf('  could not remove %s: %s', shortName(files{i}), ME.message));
    end
end
if isempty(gone), return; end
ctx.log(sprintf('  removed %d superseded result file(s): %s', numel(gone), strjoin(gone,', ')));
end

% =====================================================================
function finishUnits(ctx, sid, units)
%finishUnits  ctx.afterDone, owed exactly ONCE PER FOLDED ENTRY (spec D6) whether
%   the fold ran as one call or as N - artifact collection, downstream invalidation
%   and the per-column PDF accounting all count entries, not invocations.
for k = 1:numel(units)
    h = units{k}.head;
    ctx.afterDone(h.identity, sid, h);
end
end

% =====================================================================
function [ok, ME] = invokeGuarded(ctx, step, s, fNames, rawNames)
%invokeGuarded  THE try/catch, in the one place it is allowed to be (spec D9).  An
%   interactive step goes through ctx.modalGuard, which parks the parent window; it
%   is taken and released per invocation, so a per-recording interactive step parks
%   the window once per file rather than once per fold.
ok = true; ME = [];
call = @() invokeWrapper(step, s, fNames, rawNames);
try
    if isInteractiveStep(step, s)
        ctx.modalGuard(call);
    else
        call();
    end
catch err
    ok = false; ME = err;
end
end

% =====================================================================
function s = callSettings(s, sid, part)
%callSettings  The two settings fields that belong to ONE invocation rather than to
%   the configuration: the sibling list this call's result is inherited onto, and
%   the per-animal reference file.  'part' is a batch or a single unit - they carry
%   the same four file fields, shaped the same way, which is exactly what lets an
%   independent step be cut up without anyone re-deriving the layout.
if ~isempty(part.copyTo)
    s.fNamesCopyTo = part.copyTo;
elseif isfield(s,'fNamesCopyTo')
    s = rmfield(s,'fNamesCopyTo');                % never leak the fold's list into a unit call
end
if strcmp(sid,'vesselTypes') && ~isempty(part.refName)
    s.refFName = part.refName;                    % per-animal paint reference (launcher idiom)
end
end

% =====================================================================
function reportSkipped(ctx, skipped)
%reportSkipped  Entries that never became a call at all.  Only STRUCTURAL problems
%   reach here now - an unknown step, a recording the host cannot model - since a
%   missing input is not knowable until the moment the call is made.
for i = 1:numel(skipped)
    sk = skipped(i);
    if sk.mark, ctx.setState(sk.entry.identity, sk.stepId, 'error', sk.reason); end
    ctx.log(sprintf('skip %s - %s: %s', sk.stepLabel, sk.who, sk.reason));
end
end

% =====================================================================
function markQueued(ctx, batches)
%markQueued  Say what this run is ABOUT to do, before it does any of it (spec D2).
%   One call per batch, each one carrying that batch's whole recording list, so a
%   200-file project repaints once per call rather than once per recording.  It is
%   deliberately the FIRST thing after the plan exists and before any file is
%   resolved: 'queued' is a statement of intent, not a claim that the input is
%   there, and a batch whose input turns out to be missing is re-marked 'error'
%   when the loop reaches it.
for i = 1:numel(batches)
    markCells(ctx, batches(i).stepId, batches(i).models, 'queued', '');
end
end

% =====================================================================
function markCells(ctx, sid, models, state, msg)
%markCells  Set the run-state of one step across every recording the call touched,
%   so the matrix reflects the whole call together (spec D4).  The identities go
%   over in ONE setState - the host folds them into a single repaint (spec D3).
if isempty(models), return; end
ctx.setState({models.identity}, sid, state, msg);
end

% =====================================================================
function tf = isOrdered(b)
%isOrdered  Whether this call reads ACROSS its file list (registry fileOrder), and
%   so may neither be cut into per-recording invocations nor be handed a set with a
%   failed recording missing from it.  Absent field = 'independent' (wbBatchPlan).
tf = isfield(b,'fileOrder') && strcmp(b.fileOrder,'ordered');
end

% =====================================================================
function l = stopLine(b, remaining, failed)
%stopLine  Why the run is stopping here, in the operator's words: the step it
%   refused to start, what is still red, and what will therefore not be attempted.
%   Everything named is something the user can act on - a recording and a step.
steps = {};
for i = 1:numel(remaining)
    lbl = remaining(i).step.label;
    if ~any(strcmp(lbl, steps)), steps{end+1} = lbl; end %#ok<AGROW>
end
l = sprintf(['Stopping before %s: it works across the whole set of files at once, ' ...
    'and %d did not finish (%s).  Not attempted: %s.  Fix those and run again.'], ...
    b.step.label, numel(failed), strjoin(failed,'; '), strjoin(steps,', '));
end

% =====================================================================
function l = errorLine(step, who, ME)
%errorLine  ONE line worth reading: which step, which recording, the error's own
%   identifier when it has one, what it said, and where it came from.  The file name
%   carries no folder (the folder is the same for the whole run and would push the
%   message off the line) and the frame is the top of the stack, which is the line
%   that actually threw rather than the wrapper that was called.
l = sprintf('ERROR %s - %s: %s', step.label, who, idTag(ME));
l = [l ME.message];
fr = topFrame(ME);
if ~isempty(fr), l = sprintf('%s  [%s]', l, fr); end
end

function t = idTag(ME)
t = '';
if ~isempty(ME.identifier), t = ['(' ME.identifier ') ']; end
end

function fr = topFrame(ME)
%topFrame  The throwing site as 'file:line', or '' when the error carries no stack.
fr = '';
if isempty(ME.stack), return; end
[~,nm,ex] = fileparts(ME.stack(1).file);
if isempty(nm), nm = ME.stack(1).name; ex = ''; end
fr = sprintf('%s%s:%d', nm, ex, ME.stack(1).line);
end

% =====================================================================
function n = unitName(u)
%unitName  What ONE recording of a fold answers to in a log line: the file the
%   wrapper was pointed at, without its folder, falling back to the entry's label.
n = '';
if ~isempty(u.fNames)
    real = u.fNames(~cellfun(@isempty, u.fNames(:)));
    if ~isempty(real), n = shortName(real{1}); end
end
if isempty(n), n = u.entry.label; end
end

% =====================================================================
function n = fileCount(fNames)
%fileCount  How many REAL files a laid-out list holds (a per-animal shape pads its
%   ragged rows with '').
n = nnz(~cellfun(@isempty, fNames(:)));
end

% =====================================================================
function logCopyTargets(ctx, copyTo)
%logCopyTargets  Name the siblings a step's result is inherited onto, once per call.
if isempty(copyTo), return; end
tg = copyTargetList(copyTo);
ctx.log(sprintf('  copy to %d sibling(s): %s', numel(tg), ...
    strjoin(cellfun(@shortName,tg,'UniformOutput',false),', ')));
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
