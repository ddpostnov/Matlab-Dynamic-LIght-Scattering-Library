%wbReportPdf - Combine report images into ONE PDF, a page per image.
%
%   A 60-file column leaves 60 JPGs in the workbench's link list, which is a pile
%   rather than a report.  This helper turns a list of report images into a single
%   document you can page through: ONE PAGE PER IMAGE, each page fitted to ITS OWN
%   image with no cropping and no common paper size - the pages of a run are as
%   different in shape as the figures the wrappers wrote, and forcing them onto a
%   letter sheet would either letterbox every one of them or cut the wide ones.
%   Each image is drawn into an invisible figure whose axes fill it and whose
%   aspect matches the image, and exportgraphics appends that axes to the file.
%
%   A REPORT IS A BY-PRODUCT AND MUST NEVER KILL A RUN.  Nothing here throws for
%   bad input: a file that is missing, unreadable or not an image at all is SKIPPED
%   and named in the returned outcome list, and a page that fails to export is
%   recorded the same way.  The caller logs the outcome; the batch carries on.  An
%   empty list writes nothing and returns ''.
%
%   It knows nothing about the workbench - no handles, no settings bag, no step
%   registry - so it is callable from a launcher or a test as it stands.  Where the
%   images come from is wbArtifacts' business (the ONE place that knows how a report
%   image is named); this function is handed paths and asked for a document.
%
% Syntax:
%    out = wbReportPdf(files, outFile)
%    out = wbReportPdf(files, outFile, dpi)
%
% Inputs:
%    files   - cellstr (or a single char path) of report images, in page order.
%    outFile - full path of the PDF to write.  Its folder is created if missing;
%              an existing file is REPLACED (the first page is written with
%              'Append',false, every later one appends).
%    dpi     - export resolution, default 150.
%
% Outputs:
%    out - struct describing what happened:
%          .path    the PDF written, '' when nothing was (empty or all-bad list)
%          .pages   how many pages it carries
%          .items   1xN struct('file','ok','reason') - one per INPUT file, in
%                   order, so a skipped file is always named
%          .skipped cellstr of the files that did not make it (the .items
%                   subset with ok==false), for a one-line log message
%
% Example:
%    out = wbReportPdf({'A_t_K_cm.jpg','B_t_K_cm.jpg'}, 'reports/segmentation.pdf');
%    fprintf('%d page(s), %d skipped\n', out.pages, numel(out.skipped));
%
% See also: wbArtifacts, guiWorkbench, wbExecutor, exportgraphics
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function out = wbReportPdf(files, outFile, dpi)

if nargin<3 || isempty(dpi), dpi = 150; end
out = struct('path','','pages',0, ...
    'items',struct('file',{},'ok',{},'reason',{}),'skipped',{{}});

files = asList(files);
if isempty(files), return; end
if nargin<2 || isempty(outFile)
    error('wbReportPdf:noOutput','An output PDF path is required.');
end
outFile = char(outFile);
d = fileparts(outFile);
if ~isempty(d) && ~isfolder(d), mkdir(d); end

% ONE invisible scratch figure for the whole document, deleted by the cleanup
% object whatever happens below - a failed export must not leave a window behind,
% least of all an invisible one nobody can close.  The cleanup handle is
% ANONYMOUS and captures the handle by value on purpose: an onCleanup holding a
% NESTED function handle would keep this workspace alive and never fire (which is
% exactly how the leak was found).
fh = figure('Visible','off','Color','w','Units','pixels', ...
    'Position',[80 80 400 400],'MenuBar','none','ToolBar','none', ...
    'NumberTitle','off','IntegerHandle','off');
guard = onCleanup(@() delete(fh));

for i = 1:numel(files)
    f = files{i};
    [img, why] = readImage(f);
    if isempty(img)
        out.items(end+1) = mkItem(f,false,why);
        continue
    end
    try
        ax = drawPage(fh, img);
        exportgraphics(ax, outFile, ...
            'Append', out.pages>0, 'ContentType','image', 'Resolution', dpi);
        out.pages = out.pages + 1;
        out.items(end+1) = mkItem(f,true,'');
    catch ME
        out.items(end+1) = mkItem(f,false,ME.message);
    end
end

if out.pages>0, out.path = outFile; end
if ~isempty(out.items), out.skipped = {out.items(~[out.items.ok]).file}; end
end

% =====================================================================
function ax = drawPage(fh, img)
%drawPage  One page: the scratch figure RESHAPED to the image, an axes filling it,
%   and 'axis image' inside - so the exported axes IS the image, at its own aspect
%   ratio and whole.  exportgraphics crops to the axes' content, so a page is never
%   padded to a common sheet and never cut to one either.
[h,w,~] = size(img);
clf(fh);
fh.Position = [80 80 max(2,double(w)) max(2,double(h))];
ax = axes('Parent',fh,'Units','normalized','Position',[0 0 1 1]);
image(ax,img); axis(ax,'image'); axis(ax,'off');
end

% =====================================================================
function [img, why] = readImage(f)
%readImage  An image, or empty + the reason it is not one.  A corrupt file, a .mat
%   that wandered into the list and a path that no longer exists all land here and
%   all mean the same thing to the caller: skip it, say so, carry on.
img = []; why = '';
f = char(f);
if isempty(f) || ~isfile(f), why = 'file not found'; return; end
try
    [img, map] = imread(f);
    if ~isempty(map)
        img = ind2rgb(img, map);
    elseif size(img,3)==1
        img = repmat(unitScale(img),1,1,3);      % a grayscale page must not be
    end                                          % re-coloured by the axes colormap
catch ME
    img = []; why = ME.message;
end
if isempty(img) && isempty(why), why = 'not a readable image'; end
end

% =====================================================================
function u = unitScale(x)
%unitScale  Any single-channel class -> double in [0 1], without the Image
%   Processing Toolbox (this module is meant to run anywhere).
if islogical(x)
    u = double(x);
elseif isinteger(x)
    u = double(x)/double(intmax(class(x)));
else
    u = min(1,max(0,double(x)));
end
end

% =====================================================================
function it = mkItem(f, ok, why)
it = struct('file',char(f),'ok',logical(ok),'reason',char(why));
end

% =====================================================================
function c = asList(v)
%asList  One path, a cellstr or a string array -> a 1xN cellstr with the blanks out.
c = {};
if isempty(v), return; end
if ischar(v), c = {v};
elseif isstring(v), c = cellstr(v(:)');
elseif iscell(v), c = cellfun(@char, v(:)', 'UniformOutput', false);
end
c = c(~cellfun(@isempty,c));
end
