%wbRunRange - The FROM/TO run range: configured columns, the resume frontier, slicing.
%
%   A RANGE OVER REGISTRY COLUMNS, NOT A NEW SELECTION AXIS (spec D1).  What runs
%   is decided per TYPE on the Constructor tab and expanded into a step-major list
%   of (file, step) work items; this module only decides HOW MUCH of that list to
%   execute.  Nothing here reads a figure, a handle or a settings bag - it is
%   handed the registry and an entry list and answers four questions:
%
%     'columns'  - WHICH COLUMNS ARE CONFIGURED.  Read off the entries the
%                  expansion produced, in registry order.  It is DATA, never a
%                  written-down list: a protocol running three steps offers three
%                  items and one running fifteen offers fifteen, and adding a step
%                  to wbStepRegistry stays a data edit (spec D5).
%     'frontier' - WHERE A RESUME STARTS (spec D3).  Walking the configured columns
%                  in order, the first one that is not FINISHED.  Finished-ness is
%                  the caller's predicate - the GUI passes "every (file,step) this
%                  column plans is already done on disk" - so an errored column and
%                  a never-touched column both stop the walk, which is what makes
%                  'Last valid' mean "carry on from where it broke".  Nothing
%                  finished gives the first column (a full run); EVERYTHING finished
%                  gives the LAST column, so a resume is a no-op rather than a
%                  silent re-run of the whole protocol.  It is the DEFAULT From and
%                  nothing else - see the note below.
%     'slice'    - THE RANGE ITSELF.  The entries whose step sits in [from..to] BY
%                  REGISTRY INDEX, order preserved.  The expansion's own step-major
%                  order is never disturbed - a slice is a filter, not a re-sort.
%     'resolve'  - THE TWO DROPDOWN VALUES -> two real step ids.  'From' may be the
%                  sentinel (see 'lastValid'), which resolves to the frontier; 'To'
%                  is clamped to be at or after 'From' (spec D4), so the pair can
%                  never describe an empty or inverted range.
%
%   A RANGE NO LONGER DECIDES WHETHER FINISHED WORK RE-RUNS (author, 2026-07-31).
%   Until now an explicit From meant "and do it again whatever is on disk", reported
%   by an 'isForced' action - two unrelated requests carried by one gesture, with no
%   way to start at segmentation without also overwriting it and nothing in the
%   window to say so.  That is now a checkbox of its own on the Process tab, so
%   'isForced' is gone and the sentinel means exactly one thing: the DEFAULT From is
%   the frontier.  This module still does no pruning of its own - it never looks at
%   disk - and the caller decides which expansion, pruned or not, it slices.
%
% Syntax:
%    tok       = wbRunRange('lastValid')
%    cols      = wbRunRange('columns',  reg, entries)
%    id        = wbRunRange('frontier', reg, cols, isDone)
%    sel       = wbRunRange('slice',    reg, entries, fromId, toId)
%    [f,t]     = wbRunRange('resolve',  reg, cols, fromSel, toSel, frontierId)
%    lbl       = wbRunRange('label',    cols, stepId)
%    k         = wbRunRange('index',    cols, stepId)
%
% Inputs:
%    reg        - the struct array from wbStepRegistry.
%    entries    - the run-order entry list (needs a .stepId field; buildRunOrder's).
%    cols       - the 1xN struct('id','label') from 'columns'.
%    isDone     - function_handle stepId -> logical, "this column is finished".
%    fromSel    - a step id, '' or the 'lastValid' sentinel (the dropdown value).
%    toSel      - a step id or '' (the dropdown value).
%    frontierId - the resume frontier, as returned by 'frontier'.
%
% Outputs:
%    cols  - configured columns, registry order, struct('id',..,'label',..).
%    id    - a step id ('' when nothing is configured).
%    sel   - the entries inside the range, in their original order.
%    lbl   - a column's registry label ('' when it is not configured).
%    k     - a column's 1-based position in cols (0 when it is not configured).
%
% Example:
%    cols = wbRunRange('columns', reg, entries);
%    id   = wbRunRange('frontier', reg, cols, @(s) ~any(strcmp(s,{pending.stepId})));
%    [f,t]= wbRunRange('resolve', reg, cols, wbRunRange('lastValid'), '', id);
%    run  = wbRunRange('slice', reg, entries, f, t);
%
% See also: guiWorkbench, wbStepRegistry, wbTypeSelection, wbExecutor
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function varargout = wbRunRange(action, varargin)

switch action
    case 'lastValid', varargout{1} = lastValidToken();
    case 'columns',   varargout{1} = columnsOf(varargin{:});
    case 'frontier',  varargout{1} = frontierOf(varargin{:});
    case 'slice',     varargout{1} = sliceOf(varargin{:});
    case 'resolve',   [varargout{1:2}] = resolveOf(varargin{:});
    case 'label',     varargout{1} = labelOf(varargin{:});
    case 'index',     varargout{1} = indexOf(varargin{:});
    otherwise
        error('wbRunRange:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function t = lastValidToken()
%lastValidToken  The From sentinel: "start wherever the frontier turns out to be".
%   It is deliberately not a step id, so it can never collide with one (spec §2's
%   open vocabulary rule applies to steps as well as to types).
t = '__last__';
end

% =====================================================================
function cols = columnsOf(reg, entries)
%columnsOf  The configured columns, in registry order, derived from the entries the
%   expansion produced.  Pass the UNPRUNED expansion: a column that is already
%   finished is still configured, and losing it would take away the only way to
%   ask for it again.
cols = struct('id',{},'label',{});
ids = {};
if ~isempty(entries) && isfield(entries,'stepId'), ids = {entries.stepId}; end
if isempty(ids), return; end
for k = 1:numel(reg)
    if any(strcmp(reg(k).id, ids))
        cols(end+1) = struct('id',reg(k).id,'label',reg(k).label); %#ok<AGROW>
    end
end
end

% =====================================================================
function id = frontierOf(~, cols, isDone)
%frontierOf  The first configured column that is NOT finished (spec D3).  With
%   everything finished it is the LAST column, so 'Last valid' + 'to the end' asks
%   for nothing instead of quietly repeating the whole protocol.
id = '';
if isempty(cols), return; end
for i = 1:numel(cols)
    if ~isDone(cols(i).id), id = cols(i).id; return; end
end
id = cols(end).id;
end

% =====================================================================
function sel = sliceOf(reg, entries, fromId, toId)
%sliceOf  The entries whose step lies in [from..to] by REGISTRY INDEX.  Order is
%   preserved: the expansion decided it (step-major, the launcher's order) and a
%   range must not reorder what it merely narrows.
sel = entries;
if isempty(entries), return; end
a = regIndex(reg, fromId); if isempty(a), a = 1; end
b = regIndex(reg, toId);   if isempty(b), b = numel(reg); end
if b < a, b = a; end
keep = false(1,numel(entries));
for i = 1:numel(entries)
    k = regIndex(reg, entries(i).stepId);
    keep(i) = ~isempty(k) && k>=a && k<=b;
end
sel = entries(keep);
end

% =====================================================================
function [fromId, toId] = resolveOf(~, cols, fromSel, toSel, frontierId)
%resolveOf  Two dropdown values -> two real step ids.  The sentinel becomes the
%   frontier; a value naming a column that is no longer configured falls back to
%   the ends of the range rather than erroring, since the Constructor may have
%   dropped that column since the dropdown was last painted; and To is clamped to
%   be at or after From (spec D4).
fromId = ''; toId = '';
if isempty(cols), return; end
fromSel = char(fromSel); toSel = char(toSel);

if isempty(fromSel) || strcmp(fromSel, lastValidToken())
    fromId = char(frontierId);
    if indexOf(cols, fromId)==0, fromId = cols(1).id; end
elseif indexOf(cols, fromSel)>0
    fromId = fromSel;
else
    fromId = cols(1).id;
end

if indexOf(cols, toSel)>0, toId = toSel; else, toId = cols(end).id; end
if indexOf(cols, toId) < indexOf(cols, fromId), toId = fromId; end
end

% =====================================================================
function lbl = labelOf(cols, stepId)
%labelOf  A column's registry label, for the toolbar and the run headline.
lbl = '';
k = indexOf(cols, stepId);
if k>0, lbl = cols(k).label; end
end

% =====================================================================
function k = indexOf(cols, stepId)
%indexOf  A column's position in the CONFIGURED list (0 = not configured).
k = 0;
if isempty(cols) || isempty(stepId), return; end
j = find(strcmp({cols.id}, char(stepId)),1);
if ~isempty(j), k = j; end
end

% =====================================================================
function k = regIndex(reg, stepId)
%regIndex  A step's position in the REGISTRY - the order a range is measured in.
k = [];
if isempty(stepId), return; end
k = find(strcmp({reg.id}, char(stepId)),1);
end
