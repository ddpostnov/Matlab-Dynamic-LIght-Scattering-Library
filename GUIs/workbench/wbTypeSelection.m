%wbTypeSelection - Which steps are selected for which (recording TYPE, BRANCH).
%
%   ONE RAW RECORDING CAN DRIVE TWO INDEPENDENT PIPELINES.  The steps whose input
%   is the RAW recording rather than a .mat product each write a NEW, independent
%   triplet - the contrast step writes '_t_K' (or '_s_K'), the internal cycle
%   '_c_K' - and every later step APPENDS to one of those products.  A type that
%   runs both therefore has TWO pipelines on one recording, not one, so the
%   Constructor asks its question per (type, BRANCH) and this set is keyed that
%   way:
%
%     RAW PRODUCERS     (inGlob is not a '*.mat' glob) - ticked once per TYPE, and
%                       ticking one is what brings its branch ROW into existence.
%                       Stored under their own branch, so the producer tick and
%                       the row it creates are the same fact.
%     DERIVED CONSUMERS (inGlob is a '*.mat' glob)     - ticked per (type, branch)
%                       row.  A row OFFERS a step when the step is branch-agnostic
%                       (registry branch '') or declares that row's branch;
%                       anything else is not offerable there ('why' says so).
%
%   NEITHER LIST IS WRITTEN DOWN.  Both are derived from the registry's inGlob and
%   branch fields, so adding a step - or a third entry step for a new modality -
%   stays a registry edit and nothing here enumerates step ids or branch names.
%
%   THE VOCABULARY IS OPEN.  A type is any token the user's regexp matched or
%   typed; nothing here enumerates, validates or branches on a type's value, and
%   the number of types is unknown until the files are scanned.  A type simply
%   appears as a key when it is first ticked.  Configuration is deliberately NOT
%   pruned when a type - or a row - stops being present: re-ticking the producer,
%   or re-typing a file, restores what it had.
%
%   TWO CASCADES, both driven by the registry's requires graph, both reported so
%   the caller can log WHY a box moved, both confined to ONE row:
%     tick   - every unticked prerequisite of the step is ticked too, so a chain
%              (BFI -> vasomotion) is constructible in one click.  The row's own
%              producer always counts as available - the row exists because that
%              producer runs - so a cardiac row never pulls in a contrast step;
%     untick - every step of that row that transitively requires it is unticked
%              too, since it could no longer run.
%
%   AND ONE EXCLUSION, driven by the registry's conflictsWith.  Two steps that are
%   alternative ways to do one thing cannot both be on: ticking one UNTICKS its
%   conflicts on that row first - taking their own dependants with them, through
%   the same untick cascade - and only then pulls its prerequisites in.  Both
%   halves are reported in 'changed', which is why the caller asks what state each
%   reported step ended in rather than assuming the list is all one direction.
%   Unticking triggers nothing: dropping one alternative does not choose another.
%
%   RECORDING-LEVEL STEPS ARE TICKED ONCE AND INHERITED.  A step whose branchScope
%   is not 'one' does not belong to a single row, because the wrapper covers the
%   whole recording in ONE call: 'copy' (setRegions, segmentation) runs on the
%   contrast-side file and propagates through s.fNamesCopyTo, and 'all'
%   (dynamicSegmentation, BFI) is handed every branch product as one fNames
%   column.  It is therefore impossible for one branch to have its BFI and
%   another not, so offering an independent tick per row would be a lie.  Such a
%   step is stored ONCE, on the ANCHOR row, and shown as INHERITED on the others
%   ('effective') - which is also what keeps the selection summary from counting it
%   twice.  Unticking it drops the dependants of EVERY row, since they all lose it.
%   The anchor is the branch of the FIRST raw producer the registry declares, the
%   same side wbExecutor>desiredStage anchors 'copy' steps to; when that row is not
%   active the ticks move to the rows that are, so nothing is lost either way.
%   Note what falls out of the registry: every branch-agnostic derived step is
%   recording-level, and every row-local ('one') one names its branch.
%
%   PER-ANIMAL STEPS ARE NOT PART OF THIS SET.  registration and vesselTypes span
%   the animal, not the type (spec D3), and live in the Constructor's own animal
%   panel.  When a cascade walks into one of them it is reported separately (the
%   third output of 'tick' / 'untick'), never written into the row set - so a row
%   ticking vascularTree tells the caller "this also needs vesselTypes, which is an
%   animal step" instead of silently inventing a per-row animal step.
%
% Syntax:
%    sel                  = wbTypeSelection('new')
%    [sel,added,animal]   = wbTypeSelection('tick',   sel, reg, type, branch, stepId)
%    [sel,removed,animal] = wbTypeSelection('untick', sel, reg, type, branch, stepId)
%    [sel,ch,an]          = wbTypeSelection('set',    sel, reg, type, branch, stepId, tf)
%    tf                   = wbTypeSelection('isOn',   sel, type, branch, stepId)
%    [tf,inherited]       = wbTypeSelection('effective', sel, reg, type, branch, stepId)
%    tf                   = wbTypeSelection('rowOn',  sel, reg, type, branch)
%    brs                  = wbTypeSelection('rows',   sel, reg, type)
%    ids                  = wbTypeSelection('steps',  sel, reg, type, branch)
%    ids                  = wbTypeSelection('inherited', sel, reg, type, branch)
%    sel                  = wbTypeSelection('copy',   sel, reg, srcType, dstType)
%    sel                  = wbTypeSelection('clear',  sel, reg, type)
%    types                = wbTypeSelection('types',  sel)
%    ids                  = wbTypeSelection('typeSteps',    reg)  % every per-file step
%    ids                  = wbTypeSelection('rawSteps',     reg)  % the top panel
%    ids                  = wbTypeSelection('derivedSteps', reg)  % the row columns
%    ids                  = wbTypeSelection('animalSteps',  reg)  % the animal panel
%    brs                  = wbTypeSelection('branches',     reg)  % the possible rows
%    id                   = wbTypeSelection('producer',  reg, branch)
%    br                   = wbTypeSelection('anchorBranch', reg)
%    tf                   = wbTypeSelection('offers',    reg, branch, stepId)
%    ids                  = wbTypeSelection('rowSteps',  reg, branch)
%    txt                  = wbTypeSelection('why',       reg, branch, stepId)
%    [k,v]                = wbTypeSelection('toCells',   sel)
%    sel                  = wbTypeSelection('fromCells', keys)
%    sel                  = wbTypeSelection('fromCells', keys, reg)   % + upgrade/prune
%
% Inputs:
%    sel    - containers.Map 'type||branch||stepId' -> true (the selected set).
%    reg    - the struct array from wbStepRegistry.
%    type   - any type token (open vocabulary - never validated).
%    branch - a row branch, i.e. one of the branches a raw producer writes.
%    stepId - a registry step id.
%    tf     - logical, the wanted state ('set' is tick/untick in one call).
%    keys   - a saved key list; with reg given, keys of the older 2-part layout
%             ('type||stepId') are upgraded to their branch rows and keys naming
%             a step the registry no longer has are dropped.
%
% Outputs:
%    sel     - the updated map (value semantics: always take the returned one).
%    added / removed - cellstr of the OTHER step ids the cascade moved, in
%              registry order, excluding the step the caller named.  A tick can
%              move a box in EITHER direction (a conflict is unticked), so ask
%              'isOn' about a reported id rather than assuming which way it went.
%    animal  - cellstr of per-animal step ids the cascade reached; the caller
%              applies them to the animal panel and logs them.
%    ids     - step ids in registry (dependency) order.
%    txt     - why a cell is not offerable on a row ('' when it is).
%
% See also: wbStepRegistry, wbPrereqs, wbInvalidate, wbSettingsModel, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function varargout = wbTypeSelection(action, varargin)

switch action
    case 'new',          varargout{1} = newSel();
    case 'isOn',         varargout{1} = isOn(varargin{:});
    case 'effective',    [varargout{1:2}] = effective(varargin{:});
    case 'rowOn',        varargout{1} = rowOn(varargin{:});
    case 'rows',         varargout{1} = activeBranches(varargin{:});
    case 'steps',        varargout{1} = selectedSteps(varargin{:});
    case 'inherited',    varargout{1} = inheritedSteps(varargin{:});
    case 'types',        varargout{1} = typesOf(varargin{:});
    case 'typeSteps',    varargout{1} = idsByArity(varargin{1},'perFile');
    case 'animalSteps',  varargout{1} = idsByArity(varargin{1},'perAnimal');
    case 'rawSteps',     varargout{1} = idsByRole(varargin{1},true);
    case 'derivedSteps', varargout{1} = idsByRole(varargin{1},false);
    case 'branches',     varargout{1} = branchesOf(varargin{1});
    case 'producer',     varargout{1} = producerFor(varargin{:});
    case 'anchorBranch', varargout{1} = anchorBranch(varargin{1});
    case 'offers',       varargout{1} = offers(varargin{:});
    case 'rowSteps',     varargout{1} = rowSteps(varargin{:});
    case 'why',          varargout{1} = whyNot(varargin{:});
    case 'copy',         varargout{1} = copyType(varargin{:});
    case 'clear',        varargout{1} = clearType(varargin{:});
    case 'toCells',      [varargout{1},varargout{2}] = toCells(varargin{:});
    case 'fromCells',    varargout{1} = fromCells(varargin{:});
    case 'tick',         [varargout{1:3}] = setOne(varargin{1:5},true);
    case 'untick',       [varargout{1:3}] = setOne(varargin{1:5},false);
    case 'set',          [varargout{1:3}] = setOne(varargin{:});
    otherwise
        error('wbTypeSelection:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function s = newSel(), s = containers.Map('KeyType','char','ValueType','any'); end

% =====================================================================
function k = cellKey(type, branch, stepId)
k = [char(type) '||' char(branch) '||' char(stepId)];
end
function [type, branch, stepId] = splitKey(k)
%splitKey  A 3-part key; a 2-part one is the pre-Phase-2b layout (no branch).
p = strsplit(k,'||');
type = p{1}; branch = ''; stepId = '';
if numel(p)==2
    stepId = p{2};
elseif numel(p)>=3
    branch = p{2}; stepId = strjoin(p(3:end),'||');
end
end

% =====================================================================
function tf = isOn(sel, type, branch, stepId)
%isOn  The RAW stored state of one cell (see 'effective' for what is displayed).
tf = isa(sel,'containers.Map') && isKey(sel, cellKey(type,branch,stepId));
end

% =====================================================================
function [tf, inh] = effective(sel, reg, type, branch, stepId)
%effective  What the cell SHOWS: its own state, or the anchor row's when this step
%   is recording-level and therefore ticked once for the whole recording.
tf  = isOn(sel, type, branch, stepId);
inh = false;
src = anchorBranch(reg);
if isSharedStep(reg,stepId) && ~isempty(src) && ~strcmp(branch,src) && ...
        offers(reg,branch,stepId) && rowOn(sel,reg,type,src)
    tf  = isOn(sel, type, src, stepId);
    inh = true;
end
end

% =====================================================================
function b = storeBranch(sel, reg, type, branch, stepId)
%storeBranch  WHERE a tick is kept: its own row, or the anchor row when the step
%   covers the whole recording in one call (so there is only ever one tick).
b = char(branch);
if ~isSharedStep(reg,stepId), return; end
src = anchorBranch(reg);
if ~isempty(src) && rowOn(sel,reg,type,src), b = src; end
end

% =====================================================================
function tf = rowOn(sel, reg, type, branch)
%rowOn  A (type,branch) row EXISTS exactly while its raw producer is ticked.
prod = producerFor(reg, branch);
tf   = ~isempty(prod) && isOn(sel, type, branch, prod);
end

% =====================================================================
function brs = activeBranches(sel, reg, type)
%activeBranches  The rows this type currently has, in registry order.
brs = branchesOf(reg);
brs = brs(cellfun(@(b) rowOn(sel,reg,type,b), brs));
end

% =====================================================================
function [sel, changed, animal] = setOne(sel, reg, type, branch, stepId, tf)
%setOne  Tick or untick one cell and run the matching cascade, within ONE row.
changed = {}; animal = {};
type = char(type); branch = char(branch); stepId = char(stepId);
ids  = {reg.id};
k = find(strcmp(stepId, ids),1);
if isempty(k), return; end

% ---- a raw producer is a TYPE-level fact: it creates or removes its row ------
if isRaw(reg(k))
    key = cellKey(type, reg(k).branch, stepId);
    if tf && ~isKey(sel,key)
        sel(key) = true;
    elseif ~tf && isKey(sel,key)
        remove(sel, key);
    end
    sel = syncSharedSteps(sel, reg, type);      % inherited ticks follow the rows
    return
end

if ~offers(reg, branch, stepId), return; end  % a greyed cell cannot be driven
branch = storeBranch(sel, reg, type, branch, stepId);   % one tick per recording

dropped = {};
if tf
    [sel, dropped] = dropConflicts(sel, reg, type, branch, stepId);
    want = ancestors(reg, sel, type, branch, stepId);   % self + what it still needs
else
    want = orphaned(reg, sel, type, branch, stepId);    % self + what would lose its input
end
[sel, changed, animal] = writeWant(sel, reg, type, branch, want, stepId, tf);
changed = [dropped, changed];

if ~tf && isSharedStep(reg, stepId)
    % it covered every branch product in one call, so every OTHER row loses it too
    brs = activeBranches(sel, reg, type);
    for i = 1:numel(brs)
        if strcmp(brs{i}, branch), continue; end
        w2 = orphaned(reg, sel, type, brs{i}, stepId);
        [sel, ch2] = writeWant(sel, reg, type, brs{i}, w2, stepId, false);
        changed = [changed, ch2]; %#ok<AGROW>
    end
end
changed = orderByRegistry(changed, ids);
animal  = orderByRegistry(animal,  ids);
end

% =====================================================================
function [sel, dropped] = dropConflicts(sel, reg, type, branch, stepId)
%dropConflicts  Untick, on THIS row, every step the one being ticked declares it
%   cannot run beside - and whatever those steps were feeding, through the ordinary
%   untick cascade.  Run BEFORE the prerequisite cascade so the incoming step's own
%   chain is never the thing that gets dropped; the registry guarantees the two
%   rules cannot want the same step (wbStepRegistry>validateRegistry).
dropped = {};
k = find(strcmp(char(stepId), {reg.id}),1);
if isempty(k) || ~isfield(reg,'conflictsWith') || isempty(reg(k).conflictsWith), return; end
cw = reshape(reg(k).conflictsWith,1,[]);
for i = 1:numel(cw)
    if ~isOn(sel, type, branch, cw{i}), continue; end
    want = orphaned(reg, sel, type, branch, cw{i});
    [sel, ch] = writeWant(sel, reg, type, branch, want, stepId, false);
    dropped = [dropped, ch]; %#ok<AGROW> (a handful, once per click - ch names cw{i})
end
end

% =====================================================================
function [sel, changed, animal] = writeWant(sel, reg, type, branch, want, stepId, tf)
%writeWant  Apply one cascade's worth of moves to a row, reporting what moved.
changed = {}; animal = {};
ids = {reg.id};
for i = 1:numel(want)
    id = want{i};
    if isAnimalStep(reg, id)
        if ~strcmp(id, stepId), animal{end+1} = id; end %#ok<AGROW>
        continue                                    % never stored in the row set
    end
    if isRaw(reg(strcmp(ids,id)))
        continue                                    % the row's producer, already on
    end
    key = cellKey(type, storeBranch(sel,reg,type,branch,id), id);
    was = isKey(sel, key);
    if tf && ~was
        sel(key) = true;
        if ~strcmp(id, stepId), changed{end+1} = id; end %#ok<AGROW>
    elseif ~tf && was
        remove(sel, key);
        if ~strcmp(id, stepId), changed{end+1} = id; end %#ok<AGROW>
    end
end
end

% =====================================================================
function out = ancestors(reg, sel, type, branch, stepId)
%ancestors  stepId plus the prerequisites it is still MISSING, in registry order.
%   Missing-ness is judged against what this ROW already has, INCLUDING its own raw
%   producer - the row exists because that producer runs - so a step whose
%   requirement is satisfiable several ways (wbPrereqs requiresAny) pulls in
%   nothing extra: a cardiac row never acquires the contrast step it does not run.
ids  = {reg.id};
mask = false(1,numel(reg));
mask(strcmp(ids,stepId)) = true;
changed = true;
while changed
    changed = false;
    for k = find(mask)
        have = [rowHave(sel,reg,type,branch), ids(mask)];   % ticked, or about to be
        need = wbPrereqs('missing', reg(k), have);
        for i = 1:numel(need)
            j = find(strcmp(ids,need{i}),1);
            if isempty(j) || mask(j), continue; end
            if isRaw(reg(j)), continue; end          % the row already has its producer
            if ~offers(reg,branch,need{i}) && ~isAnimalStep(reg,need{i}), continue; end
            mask(j) = true; changed = true;
        end
    end
end
out = ids(mask);
end

% =====================================================================
function out = orphaned(reg, sel, type, branch, stepId)
%orphaned  stepId plus every step OF THIS ROW that would be left without an input.
%   Selection-aware on purpose: a step fed by EITHER entry product survives losing
%   one of them.  (wbInvalidate's staleness closure is broader by design - see
%   wbPrereqs - but a selection must only drop what truly could not run.)
ids  = {reg.id};
gone = {stepId};
changed = true;
while changed
    changed = false;
    for k = 1:numel(reg)
        id = ids{k};
        % per-animal steps are not in this set but are still worth REPORTING when
        % the row stops feeding them - the caller decides what to say about it
        onRow = isOn(sel,type,branch,id) || isAnimalStep(reg,id);
        if any(strcmp(id,gone)) || ~onRow, continue; end
        have = setdiff(rowHave(sel,reg,type,branch), gone, 'stable');
        if ~wbPrereqs('met', reg(k), have)
            gone{end+1} = id; %#ok<AGROW>
            changed = true;                      % a new drop may orphan more
        end
    end
end
out = ids(ismember(ids,gone));
end

% =====================================================================
function have = rowHave(sel, reg, type, branch)
%rowHave  What is available to a step of this row: the row's own raw producer plus
%   every derived step effectively on it (an INHERITED 'copy' step counts - the
%   result really is there, it was just computed on the source branch).
have = {};
prod = producerFor(reg, branch);
if ~isempty(prod), have{end+1} = prod; end
cols = rowSteps(reg, branch);
for i = 1:numel(cols)
    if effective(sel, reg, type, branch, cols{i}), have{end+1} = cols{i}; end %#ok<AGROW>
end
have = orderByRegistry(have, {reg.id});
end

% =====================================================================
function ids = selectedSteps(sel, reg, type, branch)
%selectedSteps  What this row RUNS, in registry order: its raw producer plus the
%   derived steps ticked on it.  An inherited 'copy' step belongs to the source
%   row and is deliberately absent here - that is what keeps it from being counted
%   (or executed) twice.
ids = {};
if ~rowOn(sel, reg, type, branch), return; end
ids{end+1} = producerFor(reg, branch);
cols = rowSteps(reg, branch);
for i = 1:numel(cols)
    if isOn(sel, type, branch, cols{i}), ids{end+1} = cols{i}; end %#ok<AGROW>
end
ids = orderByRegistry(ids, {reg.id});
end

% =====================================================================
function ids = inheritedSteps(sel, reg, type, branch)
%inheritedSteps  The recording-level steps this row shows as inherited from the
%   anchor row - ticked for the whole recording, run once, counted once.
ids = {};
cols = rowSteps(reg, branch);
for i = 1:numel(cols)
    [tf, inh] = effective(sel, reg, type, branch, cols{i});
    if tf && inh, ids{end+1} = cols{i}; end %#ok<AGROW>
end
ids = orderByRegistry(ids, {reg.id});
end

% =====================================================================
function sel = syncSharedSteps(sel, reg, type)
%syncSharedSteps  Keep the single stored tick of a recording-level step on a row
%   that exists.  It lives on the anchor row while that row is active; when the
%   anchor is switched off its ticks move to the rows that remain, and when it
%   comes back they move home again.  Either way the user's choice survives
%   switching a raw producer on and off.
src  = anchorBranch(reg);
if isempty(src), return; end
brs  = activeBranches(sel, reg, type);
cps  = sharedStepIds(reg);
if any(strcmp(src, brs))
    others = brs(~strcmp(brs, src));
    for i = 1:numel(cps)
        id = cps{i};
        on = isOn(sel,type,src,id) || any(cellfun(@(b) isOn(sel,type,b,id), others));
        if on, sel(cellKey(type,src,id)) = true; end
        for b = 1:numel(others)                    % they inherit; they do not store
            k = cellKey(type,others{b},id);
            if isKey(sel,k), remove(sel,k); end
        end
    end
else
    for i = 1:numel(cps)
        id = cps{i};
        if ~isOn(sel,type,src,id), continue; end
        for b = 1:numel(brs)
            if offers(reg,brs{b},id), sel(cellKey(type,brs{b},id)) = true; end
        end
    end
end
end

% =====================================================================
function t = typesOf(sel)
%typesOf  Every type this set carries a selection for (may include types no longer
%   present in the file set - configuration is kept, not pruned).
t = {};
if ~isa(sel,'containers.Map') || sel.Count==0, return; end
k = keys(sel);
for i = 1:numel(k)
    t{end+1} = splitKey(k{i}); %#ok<AGROW>
end
t = unique(t,'stable');
end

% =====================================================================
function sel = copyType(sel, reg, src, dst)
%copyType  Give dst exactly src's rows and their steps (dst's own are replaced).
src = char(src); dst = char(dst);
if strcmp(src,dst), return; end
k = keys(sel);
mine = {};
for i = 1:numel(k)
    [t,b,s] = splitKey(k{i});
    if strcmp(t,src), mine(end+1,:) = {b,s}; end %#ok<AGROW>
end
sel = clearType(sel, reg, dst);
for i = 1:size(mine,1), sel(cellKey(dst,mine{i,1},mine{i,2})) = true; end
sel = syncSharedSteps(sel, reg, dst);
end

% =====================================================================
function sel = clearType(sel, reg, type) %#ok<INUSD>
%clearType  Untick every row of one type (its settings are cleared separately).
type = char(type);
k = keys(sel);
for i = 1:numel(k)
    if strcmp(splitKey(k{i}), type), remove(sel, k{i}); end
end
end

% =====================================================================
function [k, v] = toCells(sel)
%toCells  Map -> parallel cell arrays, so a session file stays Map-free.
if isa(sel,'containers.Map') && sel.Count>0
    k = keys(sel); v = values(sel);
else
    k = {}; v = {};
end
end

% =====================================================================
function sel = fromCells(k, reg)
%fromCells  Rebuild the set from its saved key list, upgrading the pre-Phase-2b
%   2-part layout when the registry is given: a raw producer moves to the row it
%   creates, a branch-specific step to that branch, and a branch-agnostic one to
%   every row the type actually has (it was configured before rows existed, so the
%   only faithful reading is "on both pipelines").
%
%   A KEY NAMING A STEP THE REGISTRY NO LONGER HAS IS DROPPED, QUIETLY.  The
%   registry is the only list of what a step id means, and it changes over the
%   life of a project - a step is added, or one is dropped from the workbench
%   while its wrapper stays in the library (2026-07-31).  A session saved before
%   that carries the id in its selection, and the honest reading of a tick for
%   something that is no longer a step is "nothing": it names no column, so it can
%   never be seen, untangled or unticked, and keeping it would only write it back
%   out again on the next save.  Everything else about the session still loads.
%   This is THE place that knows the registry, which is why the rule lives here
%   and not in wbSession.
sel = newSel();
for i = 1:numel(k)
    if nargin<2 || isempty(reg), sel(k{i}) = true; continue; end
    if numel(strsplit(k{i},'||'))<3, continue; end      % 2-part: upgraded below
    [~,~,stepId] = splitKey(k{i});
    if any(strcmp(stepId,{reg.id})), sel(k{i}) = true; end
end
if nargin<2 || isempty(reg), return; end

legacy = {};
for i = 1:numel(k)
    if numel(strsplit(k{i},'||'))==2, legacy(end+1,:) = mapCells(k{i}); end %#ok<AGROW>
end
for pass = 1:2                                  % producers first: they make the rows
    for i = 1:size(legacy,1)
        t = legacy{i,1}; s = legacy{i,2};
        j = find(strcmp({reg.id},s),1);
        if isempty(j), continue; end
        if isRaw(reg(j))
            if pass==1, sel(cellKey(t,reg(j).branch,s)) = true; end
            continue
        end
        if pass~=2, continue; end
        brs = activeBranches(sel, reg, t);
        for b = 1:numel(brs)
            if offers(reg,brs{b},s), sel(cellKey(t,brs{b},s)) = true; end
        end
    end
end
for t = typesOf(sel), sel = syncSharedSteps(sel, reg, t{1}); end
end
function c = mapCells(key)
p = strsplit(key,'||');
c = {p{1}, p{end}};
end

% =====================================================================
function tf = isRaw(step)
%isRaw  A RAW PRODUCER reads the recording itself, not a .mat product, and so
%   writes a new independent triplet.  Derived from inGlob - never a name list.
tf = ~isempty(step) && ~endsWith(lower(char(step.inGlob)), '.mat');
end

% =====================================================================
function ids = idsByRole(reg, wantRaw)
%idsByRole  The per-file steps of one role, in dependency (registry) order.
keep = arrayfun(@(st) isRaw(st)==wantRaw && strcmp(st.arity,'perFile'), reg);
ids  = reshape({reg(keep).id},1,[]);
end

% =====================================================================
function ids = idsByArity(reg, arity)
%idsByArity  The registry ids of one arity, in dependency (registry) order.
ids = reshape({reg(strcmp({reg.arity},arity)).id},1,[]);
end

% =====================================================================
function brs = branchesOf(reg)
%branchesOf  The branches a recording can be driven down - one per raw producer,
%   in registry order.  A third entry step for a new modality adds its row here.
raw = idsByRole(reg,true);
brs = {};
for i = 1:numel(raw)
    b = reg(strcmp({reg.id},raw{i})).branch;
    if ~isempty(b), brs{end+1} = b; end %#ok<AGROW>
end
brs = unique(brs,'stable');
end

% =====================================================================
function id = producerFor(reg, branch)
%producerFor  The raw producer whose product IS this row ('' when there is none).
id = '';
raw = idsByRole(reg,true);
for i = 1:numel(raw)
    if strcmp(reg(strcmp({reg.id},raw{i})).branch, branch), id = raw{i}; return; end
end
end

% =====================================================================
function b = anchorBranch(reg)
%anchorBranch  The row a recording-level step is ticked on: the first branch the
%   registry declares (its entry step is listed first).  This mirrors
%   wbExecutor>desiredStage, which anchors 'copy' steps to the contrast stage.
b = '';
brs = branchesOf(reg);
if ~isempty(brs), b = brs{1}; end
end

% =====================================================================
function tf = offers(reg, branch, stepId)
%offers  THE row rule: a row offers a step when the step consumes a product (it is
%   not a raw producer) and is either branch-agnostic or declares this branch.
k = find(strcmp({reg.id},char(stepId)),1);
tf = false;
if isempty(k) || isRaw(reg(k)), return; end
tf = isempty(reg(k).branch) || strcmp(reg(k).branch, char(branch));
end

% =====================================================================
function ids = rowSteps(reg, branch)
%rowSteps  The columns of one row: the per-file derived steps it offers.
ids = idsByRole(reg,false);
ids = ids(cellfun(@(id) offers(reg,branch,id), ids));
end

% =====================================================================
function txt = whyNot(reg, branch, stepId)
%whyNot  Why a cell is greyed on this row ('' when it is offerable).
txt = '';
k = find(strcmp({reg.id},char(stepId)),1);
if isempty(k), return; end
if isRaw(reg(k))
    txt = sprintf(['%s reads the raw recording and writes a new product - it is ' ...
        'ticked once per type above, not per branch.'], reg(k).label);
elseif ~offers(reg,branch,stepId)
    txt = sprintf('%s only runs on the %s branch; this row is the %s branch.', ...
        reg(k).label, reg(k).branch, char(branch));
end
end

% =====================================================================
function tf = isSharedStep(reg, stepId)
%isSharedStep  A step the wrapper runs over the WHOLE recording in one call, so it
%   cannot be ticked per branch: branchScope 'copy' (computed on one branch,
%   propagated to the rest) or 'all' (handed every branch product at once).  Only
%   'one' is genuinely per row.
k = find(strcmp({reg.id},char(stepId)),1);
tf = ~isempty(k) && ~strcmp(reg(k).branchScope,'one');
end
function ids = sharedStepIds(reg)
%sharedStepIds  The derived per-file steps ticked once for the whole recording.
ids = idsByRole(reg,false);
ids = ids(cellfun(@(id) isSharedStep(reg,id), ids));
end

% =====================================================================
function tf = isAnimalStep(reg, stepId)
k = find(strcmp({reg.id},char(stepId)),1);
tf = ~isempty(k) && strcmp(reg(k).arity,'perAnimal');
end

% =====================================================================
function out = orderByRegistry(ids, regIds)
%orderByRegistry  Report members in dependency order, without repeats.
ids = unique(ids,'stable');
out = regIds(ismember(regIds, ids));
end
