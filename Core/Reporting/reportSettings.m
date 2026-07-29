%reportSettings - Drop the transport callbacks before s is written to a settings file.
%
%   A settings file is a record of HOW a recording was processed.  The three
%   optional workbench hooks are transport, not parameters: they are function
%   handles bound to a GUI that no longer exists by the time anyone reads the
%   file back, and saving them would make the same processing produce a different
%   *_s.mat depending on who launched it.
%
%   This replaces the 18 private stripHooks copies and strips EXACTLY the fields
%   they stripped - progressFcn, stageFcn, cancelFcn - and no more, so a
%   hook-free call still saves a byte-identical settings struct.
%
% Syntax:
%    s = reportSettings(s)
%
% Inputs:
%    s - the wrapper's parameter struct, with or without the hooks.
%
% Outputs:
%    s - the same struct with the three hook fields removed (a no-op when they
%        were never there).
%
% Example:
%    settings.runContrastFromRLS = reportSettings(s);
%    save(fName,'settings','-v7.3');
%
% See also: reportOpen, reportCancelled
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function s = reportSettings(s)

for h={'progressFcn','stageFcn','cancelFcn'}
    if isfield(s,h{1}), s=rmfield(s,h{1}); end
end
end
