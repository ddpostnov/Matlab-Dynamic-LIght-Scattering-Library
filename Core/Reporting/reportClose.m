%reportClose - End of a wrapper call: close the context and hand back the ledger.
%
%   THE PDF IS NOT BUILT HERE.  It used to be, and that put document assembly in
%   the middle of a processing call: a step that had just written sixty pages
%   spent another minute rasterising them before the next step could start, and
%   the wrapper - the layer that is supposed to process one recording at a time -
%   owned a decision about the whole batch.  A document is made BETWEEN steps
%   instead, by whoever is running them: a launcher cell calling makeReportPdf,
%   or the workbench's own per-column assembly.  Both call the same function; a
%   wrapper calls neither.
%
%   What is left is bookkeeping.  The ledger of images this call wrote is returned
%   in page order so the caller can assemble them without going back to disk, and
%   the call timer is read one last time.  Nothing here prints: the three lines
%   per file are the whole text surface, and a per-call summary would be a fourth.
%
%   NOTHING HERE THROWS.
%
% Syntax:
%    out = reportClose(rep)
%
% Inputs:
%    rep - the context from reportOpen.
%
% Outputs:
%    out - struct describing the call:
%          .images  the report images this call wrote, in page order
%          .elapsed seconds since reportOpen
%
% Example:
%    out = reportClose(rep);
%    makeReportPdf(out.images, fullfile(dataFolder,'contrast_report.pdf'));
%
% See also: reportOpen, reportSave, makeReportPdf
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function out = reportClose(rep)

st  = rep.state;
out = struct('images',{st('images')}, 'elapsed',toc(st('t0')));
end
