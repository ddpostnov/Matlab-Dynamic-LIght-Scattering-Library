%reportSave - Export a report page beside its recording, then delete the figure.
%
%   The image is written next to the data file it describes, named from that
%   file's stem plus _rep_<tag>.  Keeping it beside the recording is what keeps a
%   report with its data when a folder is copied, and it is how wbArtifacts finds
%   report images at all - it globs on the stem and matches by tail.
%
%   AND THE PAGE SAYS WHICH RECORDING IT IS.  The same stem is written across the
%   top of the page, so a PDF column assembled from ten steps is a stack of NAMED
%   pages rather than a stack of images.  It is done here, once, because here is
%   the only place that knows the answer: a wrapper writing a page for a SIBLING
%   (fName below) does not have the sibling's name in the banner, and left to
%   themselves the wrappers disagreed - five of ten titled their page, five did
%   not, and one titled one of its two pages and not the other.  A wrapper that
%   wants its page titled now does nothing at all.
%
%   A REPORT IS A BY-PRODUCT AND MUST NEVER KILL A RUN - NOR NARRATE ONE.  A
%   failed export costs nothing at all: no line, no warning, no entry in the
%   ledger, and THE FIGURE IS DELETED ON EVERY PATH - the cleanup object below
%   fires whether the export succeeded, failed or threw something nobody
%   anticipated.  The three lines a wrapper emits per recording are the entire
%   text surface, and a page that did not get written is not one of them.  The
%   wrappers this replaces printed a maximised window per file and never closed
%   it, which left sixty windows behind a sixty-file column and made close all
%   look necessary.
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
%    fName - (optional) name the image from THIS file instead of the current one,
%            and TITLE THE PAGE with it too - the page is about that recording, so
%            it says so.  For a page that describes a recording the banner did not
%            name: a copy target inheriting a mask (setRegions, runSegmentation) or
%            a report written outside the per-file loop (runRegistration's
%            overlays).
%
% Outputs:
%    None - the page is titled with the recording it is named after and the path is
%    appended to rep's ledger, which reportClose hands back for whoever assembles
%    the document between the steps.
%
% Example:
%    fh = reportFigure(rep,'segments');
%    ...
%    reportSave(rep,fh,'segments');            % -> Mouse1_t_K_r_rep_segments.jpg
%    reportSave(rep,fh,'regions',targetName);  % -> the SIBLING's _rep_regions.jpg
%
% See also: reportFigure, reportOpen, reportClose, makeReportPdf, exportgraphics
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

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
    return                      % nothing to name the page after - silently, see above
end

[fPath,stem] = fileparts(char(fName));
outFile      = fullfile(fPath, [stem '_rep_' tag '.jpg']);

try
    titlePage(fh, stem);
    anchorPage(fh);
    exportgraphics(fh, outFile, 'ContentType','image', 'Resolution', min(DPI,DPI_MAX));
catch
    return                      % a page is a by-product; the run does not hear about it
end

imgs         = st('images');
imgs{end+1}  = outFile;
st('images') = imgs; %#ok<NASGU>  % st is a handle Map, so this lands in rep's ledger
end

% =====================================================================
function titlePage(fh, stem)
%titlePage  The recording's name across the top of the page (see the header).
%   Underscores are softened to spaces because the default interpreter reads them
%   as subscripts, which is the form the wrappers used before this moved here.
%   IT GOES ON BEFORE THE ANCHOR: sgtitle attaches to whichever grid the figure
%   holds, and anchorPage's full-canvas axes is a grid of its own, so a title added
%   after it would sit on the frame instead of on the page.  A failure here costs
%   the title and nothing else - the page is still exported, untitled.
try
    sgtitle(fh, strrep(stem,'_',' '));
catch
end
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
