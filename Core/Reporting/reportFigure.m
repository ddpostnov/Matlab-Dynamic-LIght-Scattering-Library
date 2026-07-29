%reportFigure - A canonical, invisible, white figure to draw a report page into.
%
%   THE REPORT FIGURE AND THE INTERACTIVE FIGURE ARE TWO DIFFERENT OBJECTS.  An
%   interactive figure is shown, maximised to whatever monitor the operator has,
%   and sized for a human; reusing it as the report page is what makes a JPEG's
%   pixel size a function of the screen it was drawn on, and pages assembled from
%   a laptop run and a workstation run then do not match.  This one is invisible,
%   fixed in size, and exists only long enough to be exported.
%
%   THE GEOMETRY CONSTANTS LIVE HERE AND NOWHERE ELSE.  1400 x 700 px for the
%   usual two-panel page, 900 x 700 for a single panel, white background, and the
%   font defaults set ON THE FIGURE so that every axes and every text object in
%   the page inherits them.  Every page then matches because every figure did -
%   that is guaranteed at the source, not afterwards in the PDF.
%
%   The figure has no menu bar, no toolbar and no integer handle, so it cannot be
%   picked up by a stray gcf, cannot be closed by close all in some other code,
%   and costs nothing to create.  reportSave deletes it.
%
% Syntax:
%    fh = reportFigure(rep, tag)
%    fh = reportFigure(rep, tag, layout)
%
% Inputs:
%    rep    - the context from reportOpen.
%    tag    - short name of this page ('contrast', 'segments'); it becomes the
%             figure's Name and, in reportSave, the image's _rep_<tag> tail.
%    layout - (optional) 'wide' (default, 1400 x 700) or 'single' (900 x 700).
%
% Outputs:
%    fh - the figure handle.  Draw into it, then hand it to reportSave.
%
% Example:
%    fh = reportFigure(rep,'contrast');
%    t  = tiledlayout(fh,1,2);
%    nexttile(t); imagesc(imgBFI); axis image off; title('BFI');
%    nexttile(t); imagesc(mask);   axis image off; title('Mask');
%    reportSave(rep,fh,'contrast');
%
% See also: reportOpen, reportSave, reportClose, exportgraphics
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function fh = reportFigure(rep, tag, layout)

if nargin<2, tag = ''; end
if nargin<3 || isempty(layout), layout = 'wide'; end
tag = char(string(tag));

% ---- the canonical page (00-analysis Q5) ----------------------------------
W = 1400; H = 700;                       % two panels side by side
if strcmpi(char(layout),'single'), W = 900; end
FONTSIZE = 11;
FONTNAME = 'Helvetica';

fh = figure('Visible','off','Color','w','Units','pixels', ...
    'Position',[100 100 W H], ...
    'MenuBar','none','ToolBar','none','NumberTitle','off','IntegerHandle','off', ...
    'Name',strtrim([rep.procLabel ' ' tag]),'Tag',['report_' tag], ...
    'DefaultAxesFontSize',FONTSIZE,'DefaultTextFontSize',FONTSIZE, ...
    'DefaultAxesFontName',FONTNAME);
end
