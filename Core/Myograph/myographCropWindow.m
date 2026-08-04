%myographCropWindow  Whether a new time crop narrows a myograph measurement or clears it
%
%   [ivT,names] = myographCropWindow(results,crop) decides what a moved time crop does
%   to the measurement a myograph product already holds, and it is the ONE place that
%   decision is made.  It answers the crop CLIPPED to the window it falls inside, or
%   NOTHING - and nothing means the windows have to be cleared, because the stretch
%   the crop asks for was never measured.
%
%   THE DECISION IT MAKES, AND WHY IT IS WORTH A FILE OF ITS OWN.  Changing a crop
%   after a diameter exists discards something either way, and which thing depends
%   entirely on where the new crop falls:
%
%     INSIDE what was already measured  -> the windows are NARROWED to it.  The frames
%       outside go and the measurement inside is kept, so tightening a crop after a
%       first look costs the frames it drops rather than the whole run.  That is the
%       ordinary reason to come back to the crop step.
%     ANYWHERE ELSE                     -> the windows are CLEARED.  The stretch asked
%       for was never read, and a window answering only the part of it that happens to
%       exist would under-report without saying so.
%
%   It used to be twelve lines inside setMyographCrop, where nothing could reach it:
%   that wrapper opens the interval editor, which blocks until Done, so a test that
%   drove it would HANG rather than fail.  A decision that throws a measurement away
%   deserves better than being reachable only by hand, and a readyFcn riding in the
%   wrapper's s is not the answer either - editMyographIntervals' header forbids it.
%   So it is a core, and it is tested like one.
%
%   THE COMPARISON IS IN FRAMES, NOT IN SECONDS, and getting that wrong is expensive.
%   A crop is DRAWN, so its edges land between frames: a boundary a hundredth of a
%   second outside the measured window asks for no frame that is missing.  Compared
%   exactly, such a crop falls through to the clear branch and throws a whole
%   measurement away for a sliver holding nothing.  So a crop is outside only when it
%   reaches a WHOLE FRAME past either end - the first point at which a frame it asks
%   for really was never read.
%
%   A WINDOW IS A CONTIGUOUS STRETCH OF FRAMES, so "inside what was measured" means
%   inside ONE of them.  A crop spanning two disjoint windows spans the gap between
%   them as well, and the gap was never read; no crop at all means the whole
%   recording, which is the widest ask there is.  Both answer nothing.
%
%   Syntax:
%      [ivT,names] = myographCropWindow(results,crop)
%
%   INPUTS
%     results  the recording's RESULTS struct (myographProduct 'open').  Only
%              .intervals(k).diameter.time and .recording.fs are read.
%     crop     [t0 t1] the new time crop, ABSOLUTE seconds.  Empty means no crop,
%              which is the whole recording.
%
%   OUTPUTS
%     ivT      [1 x 2] the crop clipped to the window it lies inside, ready for
%              cutMyographIntervals - or [0 x 2] when the windows must be cleared.
%     names    {1 x 1} that window's name, or {}.  A window that never had one is
%              called 'analysed window', which is what it is.
%
%   EXAMPLE
%     [ivT,names] = myographCropWindow(results,[120 180]);
%     if isempty(ivT)
%         results.intervals = myographProduct('intervals');       % cleared
%     else
%         results.intervals = cutMyographIntervals(results,ivT,names);   % narrowed
%     end
%
%   DEPENDS ON
%     nothing.
%
% See also: setMyographCrop, cutMyographIntervals, myographProduct,
%           setMyographIntervals, runMyographDiameter
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function [ivT,names] = myographCropWindow(results,crop)

ivT=zeros(0,2); names={};
if numel(crop)~=2, return; end
ivs=fieldOr(results,'intervals',[]);
if isempty(ivs) || ~isstruct(ivs), return; end

lo=min(crop); hi=max(crop);
for k=1:1:numel(ivs)
    t=double(fieldOr(fieldOr(ivs(k),'diameter',[]),'time',[]));
    if isempty(t), continue; end
    step=frameStep(results,t);
    if lo<=t(1)-step || hi>=t(end)+step, continue; end
    ivT=[double(max(lo,t(1))) double(min(hi,t(end)))];
    nm=char(fieldOr(ivs(k),'name',''));
    if isempty(nm), nm='analysed window'; end
    names={nm};
    return
end
ivT=zeros(0,2); names={};
end

% =====================================================================
function step=frameStep(results,t)
%frameStep  One frame, in seconds.  The recording says so; a window whose recording
%   did not is read off its own times, which are uniform inside a window.  Zero when
%   neither can answer, which makes the comparison above the exact one it was before
%   the tolerance existed - never wider than the data can justify.
step=0;
fs=double(fieldOr(fieldOr(results,'recording',struct()),'fs',[]));
if ~isempty(fs) && isfinite(fs) && fs>0, step=1/fs; return; end
if numel(t)>1, step=median(diff(t)); end
if ~isfinite(step) || step<0, step=0; end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
