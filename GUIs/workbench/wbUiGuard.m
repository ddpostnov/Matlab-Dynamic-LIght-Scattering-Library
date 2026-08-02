%wbUiGuard - One user action at a time, for the workbench window.
%
%   MATLAB QUEUES CALLBACKS, AND A QUEUE IS THE WRONG ANSWER HERE.  A workbench
%   run holds the main thread inside a wrapper for minutes.  Clicks made in that
%   time are not lost: they are queued, and the next drawnow dispatches every one
%   of them at once.  The run path is full of drawnow (that is what makes the
%   monitor tell the truth), so a user who presses Run twice, or presses Run and
%   then Preview order while the window looks frozen, gets both actions - the
%   second one on top of a window whose state the first has already changed.
%
%   A DROPPED CLICK IS THE CORRECT ANSWER, not a deferred one.  Every command in
%   this window is an action on the CURRENT state: Run runs what is configured
%   now, Scan replaces the file list, a tick changes the configuration.  Replaying
%   a click the user made against a picture that no longer exists is not what they
%   asked for - it is what they asked for one minute ago.  Nobody presses a button
%   twice meaning "do it twice"; they press it twice meaning "did that register?".
%   So the second one is DROPPED, silently: it is not an error, it is not worth a
%   log line, and telling the user off for it would be noise.
%
%   The latch lives in the figure's appdata under its OWN key, deliberately NOT
%   inside the 'app' struct: a callback that fired in the middle of another one's
%   getApp/setApp round-trip would write a stale copy back and clobber it.
%
%   NOT EVERY CALLBACK TAKES IT.  Stop must work precisely BECAUSE the latch is
%   held, Exit must be able to interrupt a run, and a resize is not a command.
%   Those three stay unwrapped; the caller decides, this module only latches.
%
% Syntax:
%    h   = wbUiGuard('wrap', fig, fcn)
%    tf  = wbUiGuard('busy', fig)
%    tok = wbUiGuard('hold', fig)
%
% Inputs:
%    action - 'wrap', 'busy' or 'hold'.
%    fig    - the workbench uifigure the latch belongs to (it is per figure).
%    fcn    - 'wrap' only: the callback to protect.  It keeps its own signature;
%             the returned handle forwards whatever MATLAB passes it.
%
% Outputs:
%    h   - 'wrap': a callback handle to assign to a ButtonPushedFcn (or any other
%          callback property).  Fired while the latch is held it returns at once
%          and does nothing; otherwise it takes the latch, calls fcn and releases
%          it through an onCleanup, so an error inside a callback can never leave
%          the window permanently deaf.
%    tf  - 'busy': whether the latch is held right now.
%    tok - 'hold': an onCleanup token holding the latch for as long as the token
%          is alive.  For work that outlives its own callback - a run - where the
%          latch must cover the executor loop and the unwind that follows it.
%          Holding is re-entrant (it nests inside 'wrap'), because the Run button
%          is itself a wrapped callback.
%
% Notes:
%    * NEVER THROWS on the latch itself.  It is read from callbacks and released
%      from an unwind path, including one where the figure has just been deleted.
%    * The latch is not a lock: everything here runs on the main thread, and it
%      only ever discriminates between a callback dispatched from a drawnow inside
%      another callback and one dispatched from an idle window.
%
% Example:
%    b = uibutton(p,'Text','Run','ButtonPushedFcn',wbUiGuard('wrap',fig,@(~,~)runChecked(fig)));
%
% See also: guiWorkbench, wbExecutor, wbModalGuard
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function out = wbUiGuard(action, fig, fcn)

if nargin<3, fcn = []; end
switch lower(char(action))
    case 'wrap', out = @(varargin) fire(fig, fcn, varargin{:});
    case 'busy', out = depthOf(fig) > 0;
    case 'hold', out = takeLatch(fig);
    otherwise,   error('wbUiGuard: action must be ''wrap'', ''busy'' or ''hold''.');
end
end

% =====================================================================
function fire(fig, fcn, varargin)
%fire  The wrapped callback.  A click that arrives while the latch is held is the
%   replay of one made against a window that has since moved on, so it is dropped
%   here and nowhere else - silently, and before anything is taken.
if depthOf(fig) > 0, return; end
rel = takeLatch(fig); %#ok<NASGU>            released on normal OR error exit
fcn(varargin{:});
end

% =====================================================================
function tok = takeLatch(fig)
%takeLatch  Hold the latch until the returned token is cleared.  It counts rather
%   than flips, so a 'hold' inside a 'wrap' - which is exactly what the Run button
%   does - releases to held, not to free, when the inner one goes.
setDepth(fig, depthOf(fig) + 1);
tok = onCleanup(@() setDepth(fig, depthOf(fig) - 1));
end

% =====================================================================
function n = depthOf(fig)
%depthOf  How many holders the latch has.  A window that is gone holds nothing.
n = 0;
try
    if ~isempty(fig) && isvalid(fig) && isappdata(fig, latchKey())
        n = getappdata(fig, latchKey());
    end
catch
end
end

% =====================================================================
function setDepth(fig, n)
%setDepth  Write the latch back, clamped at zero so an unbalanced release (a
%   callback that deleted the figure and came back) can never leave it negative.
try
    if ~isempty(fig) && isvalid(fig)
        setappdata(fig, latchKey(), max(0, n));
    end
catch
end
end

% =====================================================================
function k = latchKey()
%latchKey  The appdata name.  Its own key, NOT a field of 'app' - see the header.
k = 'wbUiGuardDepth';
end
