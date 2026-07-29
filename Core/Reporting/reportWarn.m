%reportWarn - Something did not work, but the run carries on.
%
%   The 'warn' class: unlike a stage line it reaches the host log as well as the
%   command window, and unlike everything except 'error' it survives
%   s.reportLevel='quiet'.  That is the whole point - a quiet batch run still has
%   to say which recordings it could not process.
%
%   USE IT, NOT warning().  MATLAB's warning() writes to stderr only, so a
%   workbench operator watching the log sees nothing and finds out at the end
%   that a segment was skipped.  Reserve warning() for a genuine programming
%   fault the caller should fix; use this for the ordinary "this recording did
%   not cooperate" case, which is data, not a bug.
%
%   It does NOT repeat the file name - the banner above it already named the
%   recording.
%
% Syntax:
%    reportWarn(rep, text)
%
% Inputs:
%    rep  - the context from reportOpen.
%    text - what could not be done, and enough detail to act on it ('segment 42
%           reaches outside the image - skipped').
%
% Outputs:
%    None.
%
% Example:
%    reportWarn(rep, sprintf('segment %d reaches outside the image - skipped',k));
%
% See also: reportOpen, reportStage, reportSaved
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reportWarn(rep, text)

rep.emit('warn', ['  ! ' char(string(text))]);
end
