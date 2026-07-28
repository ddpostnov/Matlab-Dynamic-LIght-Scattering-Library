%wbModalGuard - Park the workbench window around a blocking interactive step.
%
%   The interactive pipeline steps (setRegions, setVesselTypes, setVascularTree,
%   and the conditionally-blocking contrast/registration/externalCycle) open
%   their OWN editor figures and block until the user presses Done.  They cannot
%   be made modal from outside.  wbModalGuard covers the workbench figure with an
%   overlay panel ("finish the pop-up window"), guards its close button, runs the
%   step, and - always, even if the step errors - removes the overlay and
%   restores the close callback on return.
%
%   Because execution is synchronous the run-loop is already parked inside the
%   step call while the editor is open; the overlay only stops stray clicks / a
%   close on the PARENT window.
%
% Syntax:
%    wbModalGuard(fig, fcn)
%    wbModalGuard(fig, fcn, message)
%
% Inputs:
%    fig     - the workbench uifigure to park.
%    fcn     - zero-argument function handle to run while the window is parked
%              (typically the interactive wrapper call).
%    message - optional overlay text (defaults to the standard notice).
%
% Notes:
%    uipanel has no alpha channel, so the overlay is a solid light cover rather
%    than a true semi-transparent scrim; it still visually parks the window and
%    intercepts clicks over its area.  The guarded CloseRequestFcn refuses to
%    close and requests a stop-after-this-step (sets app.cancel when present).
%
% See also: wbExecutor, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function wbModalGuard(fig, fcn, message)

if nargin<3 || isempty(message)
    message = 'Step in progress - finish the pop-up window to continue.';
end
if ~isvalid(fig), fcn(); return; end

% ---- park the window: overlay + guarded close ------------------------------
oldClose = fig.CloseRequestFcn;
overlay = makeOverlay(fig, message);
fig.CloseRequestFcn = @(~,~) guardClose(fig);
restore = onCleanup(@() unpark(fig, overlay, oldClose));
drawnow;                                     % render the overlay before we block

% ---- run the blocking step -------------------------------------------------
fcn();                                       % onCleanup restores on normal OR error exit
end

% =====================================================================
function overlay = makeOverlay(fig, message)
%makeOverlay  A solid cover panel with a centred notice, tagged for cleanup.
pos = fig.Position;
overlay = uipanel(fig,'Tag','wbModalOverlay','BorderType','none', ...
    'BackgroundColor',[0.90 0.90 0.94], ...
    'Position',[1 1 max(1,pos(3)) max(1,pos(4))], ...
    'AutoResizeChildren','off');
gl = uigridlayout(overlay,[3 1],'RowHeight',{'1x','fit','1x'},'Padding',[20 20 20 20]);
box = uigridlayout(gl,[2 1],'RowHeight',{'fit','fit'},'RowSpacing',6);
box.Layout.Row = 2;
uilabel(box,'Text',['[busy]  ' message],'FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','center');
uilabel(box,'Text','The workbench is disabled until the pop-up window is closed.', ...
    'HorizontalAlignment','center','FontAngle','italic');
uistack(overlay,'top');
end

% =====================================================================
function unpark(fig, overlay, oldClose)
%unpark  Remove the overlay and restore the close callback (cleanup path).
if isgraphics(overlay), delete(overlay); end
if isvalid(fig), fig.CloseRequestFcn = oldClose; end
end

% =====================================================================
function guardClose(fig)
%guardClose  Refuse to close mid-step; request a stop after the current step.
try
    app = getappdata(fig,'app');
    if isstruct(app) && isfield(app,'cancel')
        app.cancel = true; setappdata(fig,'app',app);
    end
catch
end
try
    uialert(fig,['A step''s pop-up window is open. Finish it to continue; ' ...
        'the batch will stop after this step.'],'Step in progress');
catch
end
end
