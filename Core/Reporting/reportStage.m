%reportStage - One line naming the procedure a wrapper has just entered.
%
%   Indented two spaces under the file banner, and carrying NO file name: the
%   banner already said which recording this is, and repeating it on every stage
%   line is what made the old output read as a transcript rather than a list of
%   recordings.  The "<procedure> on <file>" form belongs only where there is no
%   banner - a single-file or interactive context.
%
%   Stage lines are command-window detail.  They are not forwarded to the host
%   log, which already shows which step is running; see reportOpen for the
%   class-to-sink table.
%
% Syntax:
%    reportStage(rep, text)
%
% Inputs:
%    rep  - the context from reportOpen.
%    text - the procedure name ('Speckle contrast', 'Saving').
%
% Outputs:
%    None.
%
% Example:
%    reportStage(rep,'Speckle contrast');     % ->   Speckle contrast
%
% See also: reportOpen, reportFile, reportProgress
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reportStage(rep, text)

rep.emit('stage', ['  ' char(string(text))]);
end
