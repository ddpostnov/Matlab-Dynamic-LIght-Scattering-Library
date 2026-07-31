%reportSaved - The closing line of a file, with the time that file took.
%
%       Starting Contrast on Mouse1.rls
%       Writing results
%       Finished Contrast on Mouse1.rls, time elapsed: 12.4 s   <- this one
%
%   It repeats the procedure and the file VERBATIM from the Starting line, so the
%   pair can be read as a block however many recordings, wrappers or workbench
%   batches are interleaved between them.  That repetition is deliberate: a
%   closing line that says only "done" is unreadable in a sixty-file log.
%
%   The elapsed time is measured from the file's own Starting line, so a wrapper
%   never keeps a tic of its own, and it covers everything that file cost -
%   loading, computing and writing alike.
%
%   It does NOT count the files it wrote.  How many .mat files a step leaves
%   behind is a property of the step, not news about this recording.
%
% Syntax:
%    reportSaved(rep)
%
% Inputs:
%    rep - the context from reportOpen, with a current file armed by reportFile.
%
% Outputs:
%    None.
%
% Example:
%    reportWriting(rep);
%    save(dName,'source','-v7.3'); save(rName,'results','-v7.3');
%    save(sName,'settings','-v7.3');
%    reportSaved(rep);   % -> Finished Contrast on Mouse1.rls, time elapsed: 12.4 s
%
% See also: reportOpen, reportFile, reportWriting, reportClose
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function reportSaved(rep)

st = rep.state;
% Since the Starting line when there is one; since the call otherwise, so a
% wrapper that saves without ever announcing a file still reports a sane number.
if isKey(st,'fT0') && ~isempty(st('fT0'))
    secs = toc(st('fT0'));
else
    secs = toc(st('t0'));
end

rep.emit(sprintf('Finished %s on %s, time elapsed: %.1f s', ...
    rep.procLabel, shortName(st('fName')), secs));
end

% =====================================================================
function n = shortName(f)
%shortName  Name and extension - the same rendering the Starting line used, so the
%   two lines of a file are textually identical up to the verb.
[~,b,e] = fileparts(char(f));
n = [b e];
end
