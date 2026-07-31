%wbPool - Parallel-pool lifecycle for one workbench run
%
%   The first parfor in a run auto-starts the default 'Processes' pool and NOTHING
%   in the library ever takes it down.  Measured on the author's machine: 16 workers
%   pin ~1.05 GB each and never give it back, so a finished run sat on ~18 GB until
%   MATLAB was closed, and delete(gcp('nocreate')) recovered all of it in under 5 s.
%   The footprint also RATCHETS - a worker's heap is a high-water mark over the whole
%   batch, so the single largest recording sets the resident cost of every later one.
%
%   Two decisions shape this function and both are deliberate:
%
%   THE WORKBENCH OWNS THE POOL, THE LAUNCHERS DO NOT.  Pool lifecycle is a
%   GUI-boundary concern.  Nothing in Core/ or Wrappers/ knows this function exists;
%   they keep starting a pool implicitly by running a parfor, exactly as a launcher
%   run does, and the workbench closes the brackets around them.  That is why the
%   file lives here and not in Core/.
%
%   ONLY RELEASE A POOL THE RUN ITSELF STARTED.  If the author already had a pool
%   open for his own console work, a workbench run must not close it - so 'open'
%   records whether one existed BEFORE the run and 'close' releases only when it did
%   not.  'open' starts nothing: a run whose every parfor is switched off (see the
%   s.parfor* fields) never creates a pool, and this is then a no-op end to end.
%
%   Pool SIZE is not this function's business.  Sizing is a performance decision with
%   no measurement behind it, and shutting the pool down already recovers 100 % of the
%   leak.  The cost of the shutdown is that the next run pays pool startup again -
%   measured at 18 s, which is noise against a multi-hour batch.
%
% Syntax:
%    tok = wbPool('open')
%    n   = wbPool('close', tok)
%
% Inputs:
%    action - 'open' or 'close'.
%    tok    - 'close' only: the token 'open' returned.  A missing, empty or
%             malformed token is treated as "there WAS a pool before", i.e. hands
%             off - the safe reading, since wrongly closing someone else's pool is
%             the only outcome here that loses work.
%
% Outputs:
%    tok - 'open': a struct recording whether a pool was already alive.
%    n   - 'close': how many workers were released (0 if none).
%
% Notes:
%    NEVER THROWS.  It is called from a run's onCleanup path, where an exception
%    would mask the error that got the run there in the first place.
%
% Example:
%    tok = wbPool('open');
%    ...run the pipeline...
%    n = wbPool('close', tok);
%    if n>0, wbLog(fig, sprintf('released %d parallel workers', n)); end
%
% Dependencies: Parallel Computing Toolbox (absent -> every call is a no-op).
% See also: guiWorkbench, wbExecutor
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 30-July-2026

%------------- BEGIN CODE --------------
function out = wbPool(action, tok)

if nargin<2, tok = []; end
switch lower(char(action))
    case 'open',  out = openPool();
    case 'close', out = closePool(tok);
    otherwise,    error('wbPool: action must be ''open'' or ''close''.');
end
end

% =====================================================================
function tok = openPool()
%openPool  Record the pre-run state.  Starts NOTHING - the token is the only
%   product, and its single job is to answer "was there a pool before this run".
tok = struct('hadPool',true);        % pessimistic default: hands off
try
    tok.hadPool = ~isempty(gcp('nocreate'));
catch
    % No Parallel Computing Toolbox, or a broken cluster profile.  Either way there
    % is nothing this run could have started, so leave hadPool true (hands off).
end
end

% =====================================================================
function n = closePool(tok)
%closePool  Release the pool, but only if THIS run started it.
n = 0;
try
    if isstruct(tok) && isscalar(tok) && isfield(tok,'hadPool') && ~tok.hadPool
        p = gcp('nocreate');
        if ~isempty(p)
            n = p.NumWorkers;
            delete(p);
        end
    end
catch
    % A pool that is already gone, mid-shutdown, or unreachable is not an error
    % worth reporting: the goal was for it not to be there.
    n = 0;
end
end
