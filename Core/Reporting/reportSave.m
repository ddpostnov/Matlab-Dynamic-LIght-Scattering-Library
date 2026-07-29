%reportSave - Export a report page beside its recording, then delete the figure.
%
%   The image is written next to the data file it describes, named from that
%   file's stem plus _rep_<tag>.  Keeping it beside the recording is what keeps a
%   report with its data when a folder is copied, and it is how wbArtifacts finds
%   report images at all - it globs on the stem and matches by tail.
%
%   A REPORT IS A BY-PRODUCT AND MUST NEVER KILL A RUN.  A failed export costs a
%   warning line and nothing else, and THE FIGURE IS DELETED ON EVERY PATH - the
%   cleanup object below fires whether the export succeeded, failed or threw
%   something nobody anticipated.  The wrappers this replaces printed a maximised
%   window per file and never closed it, which left sixty windows behind a
%   sixty-file column and made close all look necessary.
%
%   RESOLUTION IS SPECIFIED AS OUTPUT PIXELS, NOT AS DPI.  exportgraphics scales
%   a pixel-sized figure against MATLAB's nominal 96 dpi, so 150 dpi on the
%   1400 x 700 canvas lands at roughly 2100 x 1050 px, about 2.2 MP.  200 dpi is
%   the hard ceiling; test the pixel size, because dpi semantics shift with the
%   figure's Units and pixel dimensions do not.
%
% Syntax:
%    reportSave(rep, fh, tag)
%
% Inputs:
%    rep - the context from reportOpen, with a current file armed by reportFile.
%    fh  - the figure from reportFigure.  It is deleted before this returns.
%    tag - short stage name; the image's tail is _rep_<tag>.jpg.
%
% Outputs:
%    None - the path is appended to rep's ledger for reportClose to assemble.
%
% Example:
%    fh = reportFigure(rep,'segments');
%    ...
%    reportSave(rep,fh,'segments');    % -> Mouse1_t_K_d_rep_segments.jpg
%
% See also: reportFigure, reportOpen, reportClose, wbReportPdf, exportgraphics
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reportSave(rep, fh, tag)

DPI     = 150;                  % 200 is the hard ceiling - never raise this past it
DPI_MAX = 200;

tag = char(string(tag));

% The figure goes, whatever happens below.  The cleanup handle is ANONYMOUS and
% captures fh by value on purpose: an onCleanup holding a nested function handle
% would keep this workspace alive and never fire.
guard = onCleanup(@() delete(fh));

st    = rep.state;
fName = st('fName');
if isempty(fName)
    rep.emit('warn', sprintf('  report image "%s" not written (no current file)', tag));
    return
end

[fPath,stem] = fileparts(char(fName));
imgName      = [stem '_rep_' tag '.jpg'];
outFile      = fullfile(fPath, imgName);

try
    exportgraphics(fh, outFile, 'ContentType','image', 'Resolution', min(DPI,DPI_MAX));
catch ME
    rep.emit('warn', sprintf('  report image "%s" not written (%s)', tag, ME.message));
    return
end

imgs         = st('images');
imgs{end+1}  = outFile;
st('images') = imgs; %#ok<NASGU>  % st is a handle Map, so this lands in rep's ledger
rep.emit('report', ['  report  ' imgName]);
end
