%reportCancelled - Has the host asked this run to stop?
%
%   The cooperative cancel check.  The hook rides in the reporting context like
%   every other sink, so a wrapper never touches s.cancelFcn and never has to ask
%   whether it is there: reportOpen defaulted it to @()false when the caller is a
%   launcher rather than the workbench.
%
% Syntax:
%    tf = reportCancelled(rep)
%
% Inputs:
%    rep - the context from reportOpen.
%
% Outputs:
%    tf - true when the host wants the run stopped at the next file boundary.
%
% Example:
%    if reportCancelled(rep), break; end
%
% See also: reportOpen, reportFile
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function tf = reportCancelled(rep)

tf = rep.cancelFcn();
end
