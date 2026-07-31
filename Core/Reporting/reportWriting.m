%reportWriting - The middle of the three lines: the results are going to disk.
%
%       Starting Contrast on Mouse1.rls
%       Writing results                 <- this one
%       Finished Contrast on Mouse1.rls, time elapsed: 12.4 s
%
%   It goes immediately before the FIRST save of the _d/_r/_s triplet, and it is
%   the only thing that separates "still computing" from "still saving" now that
%   there is no progress axis.  That distinction is worth a line because the two
%   failure modes look identical from outside: a wrapper stuck on a 40 GB write
%   and a wrapper stuck in a solver both look like nothing happening, and this
%   line says which one you are watching.
%
%   It carries NO argument and NO file name.  The line above it named the
%   recording; a save is a save.
%
%   Intermediate saves - a report image, a scratch .mat - do NOT get this line.
%   It marks the results, which is what the next step reads back.
%
% Syntax:
%    reportWriting(rep)
%
% Inputs:
%    rep - the context from reportOpen.
%
% Outputs:
%    None.
%
% Example:
%    reportWriting(rep);
%    save(dName,'source','-v7.3');
%    save(rName,'results','-v7.3');
%    save(sName,'settings','-v7.3');
%    reportSaved(rep);
%
% See also: reportOpen, reportFile, reportSaved
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function reportWriting(rep)

rep.emit('Writing results');
end
