%reportFile - Announce that a wrapper has started a new file, and arm the context.
%
%   The first of the three lines a wrapper emits per recording:
%
%       Starting Contrast on Mouse1_t_K_r.mat
%
%   It names the PROCEDURE and the FILE AS IT IS ON DISK - name and extension,
%   from fileparts, never a path and never a stem with the role suffix chopped
%   off.  The procedure is the label given to reportOpen, and the closing line
%   from reportSaved repeats it verbatim, so the two lines of a file always pair
%   up even when several wrapper calls are interleaved in one log.
%
%   It also arms the context: every report image written until the next
%   reportFile call is named from THIS file, so reportSave needs nothing but a
%   tag; and the file's own timer starts, so reportSaved can close the block with
%   an elapsed time and no wrapper needs a bare tic/toc of its own.
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
%    None - the line is emitted and rep's state Map is updated in place.
%
% Example:
%    reportFile(rep, 3, 'D:\data\Mouse1_t_K_r.mat');
%    % -> Starting Contrast on Mouse1_t_K_r.mat
%
% See also: reportOpen, reportWriting, reportSaved, reportSave, reportClose
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function reportFile(rep, fIdx, fName)

if nargin<3, fName = ''; end
fName = char(fName);

% st is a containers.Map, so these writes land in rep's state and are visible to the
% caller; the NASGU pragma is Code Analyzer not knowing that a Map is a handle.
st = rep.state;
st('fIdx')  = fIdx;
st('fName') = fName;
st('fT0')   = tic; %#ok<NASGU>  what reportSaved reports as this file's elapsed time

rep.emit(sprintf('Starting %s on %s', rep.procLabel, shortName(fName)));
end

% =====================================================================
function n = shortName(f)
%shortName  Name and extension - what the user sees on disk.  fileparts, not a
%   split on a backslash, and the extension is kept: 'Mouse1_t_K_r.mat' is a file,
%   'Mouse1_t_K_d' is a guess at one.
[~,b,e] = fileparts(char(f));
n = [b e];
end
