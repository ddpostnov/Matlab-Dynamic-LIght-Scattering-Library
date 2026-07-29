%reportClose - End of a wrapper call: assemble the PDF, tidy up, say what happened.
%
%   A wrapper call is the launcher's "step", so this is where its report images
%   become one document you can page through instead of a pile of JPEGs in a data
%   folder.  The assembly itself is wbReportPdf, reused verbatim: it already fits
%   each page to its own image without cropping and never throws, and it knows
%   nothing about the workbench, so a launcher can use it as it stands.
%
%   THE HOST CAN OPT OUT.  The Processing Workbench sets s.reportPdf=false and
%   keeps its own per-column assembly, because a column spans several batched
%   wrapper calls and a column-level PDF is the document its operator wants.  One
%   code path in the wrapper either way.
%
%   DELETING THE SOURCE IMAGES IS OFF BY DEFAULT (s.reportKeepImages).  wbArtifacts
%   re-resolves report images from disk both to paint the result list and to
%   recover artifacts when a session is resumed, so deleting them also empties
%   that list after a reload - a surprise that should be opted into.  When it IS
%   on, deletion happens only after the PDF is confirmed written, and never
%   touches a file wbReportPdf reported as skipped: an image the document does not
%   contain is the one image you must not throw away.
%
%   NOTHING HERE THROWS.  A missing file, a .mat that wandered into the ledger, a
%   corrupt JPEG, a PDF that cannot be written - each costs a line and the run
%   carries on.
%
% Syntax:
%    out = reportClose(rep)
%
% Inputs:
%    rep - the context from reportOpen.
%
% Outputs:
%    out - struct describing the call:
%          .path    the PDF written, '' when none was
%          .pages   how many pages it carries
%          .images  the ledger, in page order
%          .skipped images that did not make it into the PDF
%          .deleted images removed afterwards (empty unless reportKeepImages is
%                   false)
%          .elapsed seconds since reportOpen
%
% Example:
%    out = reportClose(rep);
%    fprintf('%d page(s) in %s\n', out.pages, out.path);
%
% See also: reportOpen, reportSave, wbReportPdf
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function out = reportClose(rep)

st   = rep.state;
imgs = st('images');

% Close the host's bar.  reportFile moves it to "files finished", so without this
% a completed call would sit at (N-1)/N for ever.
rep.progressFcn(1, rep.procLabel);
out  = struct('path','','pages',0,'images',{imgs},'skipped',{{}},'deleted',{{}}, ...
    'elapsed',toc(st('t0')));

% ---- the document ---------------------------------------------------------
if rep.pdf && ~isempty(imgs)
    try
        out = assemble(rep, imgs, out);
    catch ME
        rep.emit('warn', ['  report PDF not written (' ME.message ')']);
    end
end

% ---- optional cleanup, after the PDF is confirmed --------------------------
if ~rep.keepImages && out.pages>0
    for i = 1:numel(imgs)
        f = imgs{i};
        if any(strcmp(f, out.skipped)), continue; end   % never delete what is not in the PDF
        try
            if isfile(f), delete(f); out.deleted{end+1} = f; end
        catch ME
            rep.emit('warn', ['  could not remove ' shortName(f) ' (' ME.message ')']);
        end
    end
end

% ---- one line ------------------------------------------------------------
txt = sprintf('  %s: %d file(s), %d report image(s), %.1f s', ...
    rep.procLabel, rep.nFiles, numel(imgs), out.elapsed);
if out.pages>0
    txt = [txt sprintf(', %d PDF page(s)', out.pages)];
end
if ~isempty(out.skipped)
    txt = [txt sprintf('; skipped %s', strjoin(cellfun(@shortName, out.skipped, ...
        'UniformOutput', false), ', '))];
end
rep.emit('save', txt);
end

% =====================================================================
function out = assemble(rep, imgs, out)
%assemble  The ledger -> one PDF beside the images, via wbReportPdf.
if exist('wbReportPdf','file')~=2
    rep.emit('warn','  report PDF not written (wbReportPdf is not on the path)');
    return
end
fPath   = fileparts(imgs{1});
outFile = fullfile(fPath, [slug(rep.procLabel) '_report.pdf']);
r       = wbReportPdf(imgs, outFile);
out.path    = r.path;
out.pages   = r.pages;
out.skipped = r.skipped;
end

% =====================================================================
function t = slug(label)
%slug  A procedure label -> a file-name fragment ('Speckle contrast' ->
%   'speckle-contrast').
t = lower(char(label));
t = regexprep(t,'[^a-z0-9]+','-');
t = regexprep(t,'^-+|-+$','');
if isempty(t), t = 'report'; end
end

% =====================================================================
function n = shortName(f)
%shortName  Name and extension only - the folder is the one thing every entry in
%   these lists has in common.
[~,b,e] = fileparts(char(f));
n = [b e];
end
