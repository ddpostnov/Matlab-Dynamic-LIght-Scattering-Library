%reportFile - Announce that a wrapper has started a new file, and arm the context.
%
%   One line per recording, and only one:
%
%       [3/12] Mouse1_t_K
%
%   It names the RECORDING, not a position: the workbench already shows where a
%   file sits in the matrix, and the command window used to lose the identity
%   entirely behind "Processing file 3 out of 12".  The stem comes from fileparts,
%   never from splitting a path on a backslash.
%
%   It also arms the context: every report image written until the next
%   reportFile call is named from THIS file, so reportSave needs nothing but a
%   tag; the progress throttle restarts so a new file's first tick is not
%   swallowed by the previous file's; and the file's own timer starts, so
%   reportSaved can close the block with an elapsed time and no wrapper needs a
%   bare tic/toc of its own.
%
% Syntax:
%    reportFile(rep, fIdx, fName)
%
% Inputs:
%    rep   - the context from reportOpen.
%    fIdx  - 1-based index of this file within the call's file list.
%    fName - full path of the file being processed.
%
% Outputs:
%    None - the banner is emitted and rep's state Map is updated in place.
%
% Example:
%    reportFile(rep, 3, 'D:\data\Mouse1_t_K_d.mat');   % -> [3/12] Mouse1_t_K_d
%
% See also: reportOpen, reportStage, reportSave, reportClose
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reportFile(rep, fIdx, fName)

if nargin<3, fName = ''; end
fName    = char(fName);
[~,stem] = fileparts(fName);

% st is a containers.Map, so these writes land in rep's state and are visible to the
% caller; the NASGU pragma is Code Analyzer not knowing that a Map is a handle.
st = rep.state;
st('fIdx')      = fIdx;
st('fName')     = fName;
st('fT0')       = tic;       % what reportSaved reports as this file's elapsed time
st('lastTickT') = -inf;      % a new file restarts the progress throttle
st('lastTickF') = -inf; %#ok<NASGU>  Map is a handle, so these writes are visible

if rep.nFiles > 0
    rep.emit('file', sprintf('[%d/%d] %s', fIdx, rep.nFiles, stem));
    % The banner IS the coarse progress, so the host's bar is moved from here and
    % the command window is spared a second rendering of what [3/12] already says.
    % It reports work FINISHED, hence fIdx-1; reportClose closes it at 100 %.
    rep.progressFcn(max(0,(fIdx-1))/rep.nFiles, rep.procLabel);
else
    rep.emit('file', sprintf('[%d] %s', fIdx, stem));
end
end
