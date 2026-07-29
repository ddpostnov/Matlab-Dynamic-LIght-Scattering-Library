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
%   1400 x 700 canvas lands at 2188 x 1094 px, about 2.4 MP.  200 dpi is the hard
%   ceiling; test the pixel size, because dpi semantics shift with the figure's
%   Units and pixel dimensions do not.
%
%   AND THE PAGE IS ANCHORED, OR THE SIZE WOULD BE A PROPERTY OF THE CONTENT.
%   exportgraphics TRIMS the white margin around whatever was drawn, so the same
%   canvas exports at 1994 x 1072 with a two-tile layout, 1752 x 953 with a lone
%   axes and 1741 x 763 with a subplot pair - three different page shapes from one
%   "canonical" figure, which is the mismatch this module exists to remove.  A
%   hairline frame drawn at the edge just before the export gives the crop
%   something to stop at, so EVERY page comes out at exactly the canvas size.  It
%   has to be added last: tiledlayout deletes axes that already exist, and white
%   or near-white anchors are trimmed away like any other margin.
%
% Syntax:
%    reportSave(rep, fh, tag)
%    reportSave(rep, fh, tag, fName)
%
% Inputs:
%    rep   - the context from reportOpen, with a current file armed by reportFile.
%    fh    - the figure from reportFigure.  It is deleted before this returns.
%    tag   - short stage name; the image's tail is _rep_<tag>.jpg.
%    fName - (optional) name the image from THIS file instead of the current one.
%            For a page that describes a recording the banner did not name: a copy
%            target inheriting a mask (setRegions, runSegmentation) or a report
%            written outside the per-file loop (runRegistration's overlays).
%
% Outputs:
%    None - the path is appended to rep's ledger for reportClose to assemble.
%
% Example:
%    fh = reportFigure(rep,'segments');
%    ...
%    reportSave(rep,fh,'segments');            % -> Mouse1_t_K_d_rep_segments.jpg
%    reportSave(rep,fh,'regions',targetName);  % -> the SIBLING's _rep_regions.jpg
%
% See also: reportFigure, reportOpen, reportClose, wbReportPdf, exportgraphics
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reportSave(rep, fh, tag, fName)

DPI     = 150;                  % 200 is the hard ceiling - never raise this past it
DPI_MAX = 200;

tag = char(string(tag));

% The figure goes, whatever happens below.  The cleanup handle is ANONYMOUS and
% captures fh by value on purpose: an onCleanup holding a nested function handle
% would keep this workspace alive and never fire.
guard = onCleanup(@() delete(fh));

st = rep.state;
if nargin<4 || isempty(fName)
    fName = st('fName');
end
if isempty(fName)
    rep.emit('warn', sprintf('  report image "%s" not written (no current file)', tag));
    return
end

[fPath,stem] = fileparts(char(fName));
imgName      = [stem '_rep_' tag '.jpg'];
outFile      = fullfile(fPath, imgName);

try
    anchorPage(fh);
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

% =====================================================================
function anchorPage(fh)
%anchorPage  The hairline page frame that fixes the exported size (see the header).
%   Drawn at the very edge of the canvas, in the lightest grey that still survives
%   exportgraphics' crop, with HandleVisibility off so it cannot be picked up by a
%   gca in code that runs afterwards.  A failure here is not worth losing the page
%   over - the export just falls back to a trimmed one.
try
    ax = axes(fh,'Units','normalized','Position',[0.0015 0.003 0.997 0.994], ...
        'HandleVisibility','off','HitTest','off','Tag','reportPageFrame', ...
        'XLim',[0 1],'YLim',[0 1],'Color','none','XTick',[],'YTick',[], ...
        'Box','on','XColor',[0.88 0.88 0.88],'YColor',[0.88 0.88 0.88]);
    uistack(ax,'bottom');
catch
end
end
