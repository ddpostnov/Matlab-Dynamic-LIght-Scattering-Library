%myographRecordingPath  The recording a myograph product was made from, if it is still there
%
%   [p,wanted] = myographRecordingPath(results,fName) finds the '.avi' (or the
%   '.adicht') that a *_MYO_r.mat was made from, and returns '' when it is not there.
%
%   THE PRODUCT DESCRIBES THE RECORDING AND DOES NOT COPY IT, so more than one step
%   has to find it again - the diameter step to measure it, the interval editor to
%   preview it, the two cutting steps to decide whether they are allowed to throw
%   anything away.  They all have to agree about where to look, or a step would
%   refuse over a file another step is happily reading.  This is where the order is
%   written down:
%
%     1. THE PATH THE CALLER PASSED, when the caller was given one.
%     2. THE PATH THE ENTRY STEP RECORDED (results.recording.fName).  It is tried
%        before the neighbours because it is the only one that knows the original
%        container's extension.
%     3. A RECORDING OF THE RIGHT KIND SITTING BESIDE THE PRODUCT, with the same
%        stem.  A whole folder moved to another disk lands here, which is the
%        ordinary way a set of results and its recordings travel.
%
%   WHICH KIND depends on the modality: a wire myograph was read from a LabChart
%   file and a pressure one from a video, and offering the wrong extensions would
%   find nothing and say the recording was missing when it is sitting there.
%
%   IT DOES NOT THROW.  A results file that travelled without its recording is a
%   normal thing to be handed; every caller has its own sentence to say about it,
%   and WANTED is the name to put in that sentence.
%
%   Syntax:
%      p            = myographRecordingPath(results,fName)
%      [p,wanted]   = myographRecordingPath(results,fName,rawName)
%
%   INPUTS
%     results  the recording's RESULTS struct (myographProduct 'open').
%     fName    the product this is being asked about - any member of the pair.
%     rawName  (optional) a path the caller was handed for this recording.  Tried
%              first, and ignored when it is not a file.
%
%   OUTPUTS
%     p        the full path of the recording, or '' when it was not found.
%     wanted   the recording AS THE OPERATOR WOULD LOOK FOR IT - name and extension,
%              never a path: the name the entry step wrote down, or the product's
%              own stem with the extension the modality implies.  Always filled, so
%              a message can name the file whether or not it exists.
%
%   EXAMPLE
%     [p,wanted] = myographRecordingPath(results,fName);
%     if isempty(p)
%         warning('Keep %s beside the results.',wanted);
%     end
%
% See also: runMyographDiameter, setMyographIntervals, setMyographCrop,
%           getMyographTrace, getMyographWallFrame, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function [p,wanted] = myographRecordingPath(results,fName,rawName)

if nargin<3, rawName=''; end
rec=fieldOr(results,'recording',struct());
exts=extensionsFor(rec);

[fPath,stem]=fileparts(regexprep(char(fName),'_MYO_[drs]\.mat$',''));

recorded=char(string(fieldOr(rec,'fName','')));
wanted=nameOnly(recorded);
if isempty(wanted), wanted=[stem exts{1}]; end

p=char(rawName);
if ~isempty(p) && isfile(p), return; end
p=recorded;
if ~isempty(p) && isfile(p), return; end
for i=1:1:numel(exts)
    cand=fullfile(fPath,[stem exts{i}]);
    if isfile(cand), p=cand; return; end
end
p='';
end

% =====================================================================
function exts=extensionsFor(rec)
%extensionsFor  What kind of file this recording is.  A wire myograph is a LabChart
%   file and a pressure one is a video; the four video containers are the ones
%   VideoReader is asked for elsewhere in this branch, in the same order.
if strcmpi(char(fieldOr(rec,'modality','PMYO')),'WMYO')
    exts={'.adicht'};
else
    exts={'.avi','.mp4','.mov','.mkv'};
end
end

% =====================================================================
function n=nameOnly(f)
%nameOnly  Name and extension - what the operator sees in a folder.
n='';
if isempty(f), return; end
[~,b,e]=fileparts(char(f));
n=[b e];
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
