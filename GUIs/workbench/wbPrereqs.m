%wbPrereqs - THE rule for "can this step run yet?", in one place.
%
%   A pipeline step declares two kinds of prerequisite, and the difference matters
%   because this library's steps are BRANCH-AGNOSTIC in the middle:
%
%     requires     ALL of these must be done/selected.  A hard chain -
%                  dynamicSegmentation genuinely cannot run before segmentation.
%     requiresAny  AT LEAST ONE of these must be done/selected.  The middle of the
%                  pipeline (regions, segmentation, BFI, registration) consumes
%                  '*_K_r.mat' - ANY branch product of the recording - so it is
%                  ready as soon as SOME entry step has produced one.  The contrast
%                  step is the usual producer, but a purely pulsatile protocol
%                  starts at the internal cycle and never computes a contrast
%                  product at all; writing 'contrast' as a hard requirement there
%                  forced protocols to tick a step they do not run.
%
%   Every consumer of the dependency graph goes through here so the selection
%   cascade (wbTypeSelection), the on-disk gating (wbStateEngine), the staleness
%   cascade (wbInvalidate) and the workbench's look-ahead can never disagree.
%
%   STALENESS IS DELIBERATELY BROADER than readiness: 'all' returns the union of
%   both lists, because if ANY producer of a step re-ran, that step's inputs
%   changed - even a producer it did not happen to use this time.
%
% Syntax:
%    tf  = wbPrereqs('met',     step, availableIds)
%    ids = wbPrereqs('missing', step, availableIds)
%    ids = wbPrereqs('all',     step)
%    ids = wbPrereqs('hard',    step)
%    ids = wbPrereqs('any',     step)
%    txt = wbPrereqs('describe',step)
%
% Inputs:
%    step         - a step struct from wbStepRegistry.
%    availableIds - cellstr of step ids that are already done (or selected).
%
% Outputs:
%    tf  - whether every requirement is satisfied by availableIds.
%    ids - 'missing': what to add to satisfy it (the hard ones that are absent,
%          plus the FIRST of an unsatisfied requiresAny list - its default
%          producer).  'all'/'hard'/'any': the declared lists.
%    txt - a short human-readable requirement line for a tooltip.
%
% See also: wbStepRegistry, wbTypeSelection, wbStateEngine, wbInvalidate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbPrereqs(action, step, availableIds)

if nargin<3, availableIds = {}; end
hard = listOf(step,'requires');
any_ = listOf(step,'requiresAny');

switch action
    case 'hard',     out = hard;
    case 'any',      out = any_;
    case 'all',      out = unique([hard, any_],'stable');
    case 'met',      out = isempty(missing(hard,any_,availableIds));
    case 'missing',  out = missing(hard,any_,availableIds);
    case 'describe', out = describe(hard,any_);
    otherwise
        error('wbPrereqs:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function ids = missing(hard, any_, availableIds)
%missing  What still has to happen before this step can run.
ids = hard(~ismember(hard, availableIds));
if ~isempty(any_) && ~any(ismember(any_, availableIds))
    ids = [ids, any_(1)];        % none of the alternatives yet: take the default one
end
ids = unique(ids,'stable');
end

% =====================================================================
function l = listOf(step, field)
l = {};
if isstruct(step) && isfield(step,field) && ~isempty(step.(field))
    l = reshape(step.(field),1,[]);
end
end

% =====================================================================
function txt = describe(hard, any_)
bits = {};
if ~isempty(hard), bits{end+1} = strjoin(hard,' + '); end
if ~isempty(any_), bits{end+1} = strjoin(any_,' or '); end
txt = strjoin(bits,', and ');
end
