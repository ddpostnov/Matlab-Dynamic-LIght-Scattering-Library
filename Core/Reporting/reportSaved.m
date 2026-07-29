%reportSaved - One line saying what a file's processing left on disk.
%
%   The closing line of a file's block, indented under the banner:
%
%       [3/12] Mouse1_t_K
%         Speckle contrast
%         Saving
%             saved  3 files, 12.4 s
%
%   It replaces the old pair of lines - "Saving the results" followed by "Saving
%   complete" or an unlabelled toc - and it does NOT repeat the file name: the
%   banner three lines up already gave it.
%
%   WHY THIS IS NOT reportStage.  Stage lines are command-window detail and are
%   not forwarded to the host log.  This one is 'save'-class, so it reaches both:
%   together with the banner it gives the workbench log the two lines per file
%   that make it read as a list of recordings rather than a transcript (see the
%   class-to-sink table in reportOpen).
%
%   The elapsed time is measured from the file's own banner, so a wrapper does not
%   have to keep a tic of its own; pass secs explicitly only when the interval you
%   want to report is not "since this file started".
%
% Syntax:
%    reportSaved(rep, nProducts)
%    reportSaved(rep, nProducts, secs)
%
% Inputs:
%    rep       - the context from reportOpen.
%    nProducts - how many files were written (the _d/_r/_s triplet is 3).
%    secs      - (optional) seconds to report instead of the time since the
%                file banner.
%
% Outputs:
%    None.
%
% Example:
%    reportStage(rep,'Saving');
%    save(dName,'source','-v7.3'); save(rName,'results','-v7.3');
%    save(sName,'settings','-v7.3');
%    reportSaved(rep,3);                  % ->     saved  3 files, 12.4 s
%
% See also: reportOpen, reportFile, reportStage, reportClose
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reportSaved(rep, nProducts, secs)

if nargin<2 || isempty(nProducts), nProducts = 0; end

st = rep.state;
if nargin<3 || isempty(secs)
    % Since the banner when there is one; since the call otherwise, so a wrapper
    % that saves without ever announcing a file still reports a sane number.
    if isKey(st,'fT0') && ~isempty(st('fT0'))
        secs = toc(st('fT0'));
    else
        secs = toc(st('t0'));
    end
end

rep.emit('save', sprintf('      saved  %d file(s), %.1f s', nProducts, secs));
end
