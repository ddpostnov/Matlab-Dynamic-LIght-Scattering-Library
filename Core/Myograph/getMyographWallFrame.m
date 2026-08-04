%getMyographWallFrame  One frame of a myograph recording, with its detected walls on it
%
%   F = getMyographWallFrame(results,k,measure,when) re-opens the recording a
%   pressure-myograph product was made from and returns one frame of it, together
%   with the wall positions window K carries for that frame.  It is the "did it find
%   the right edges" picture, and the only function in the myograph branch that
%   opens a video for display.
%
%   THE RECORDING IS RE-OPENED, NEVER COPIED.  A myograph product keeps the
%   measurement and describes the recording; the frames themselves stay in the
%   '.avi'.  Asking for one is therefore a matter of finding it again, which is
%   exactly what results.intervals(k).frames is for: it is the window's frame range
%   in the ORIGINAL recording, so frames(1)+offset names a frame of the file however
%   many times the windows have since been narrowed.
%
%   IT IS NOT AN EXPLORER VARIABLE, ON PURPOSE.  guiExplore reads the results and
%   only the results, and a video frame is not a result - that rule is why the
%   explorer's "Detected walls" view was removed and it has not changed.  What
%   changed is that a step which is ALREADY HOLDING THE RECORDING - the interval
%   editor, deciding where the windows go - is not copying anything by drawing one,
%   so the picture lives here, where such a step can reach it.
%
%   A MISSING RECORDING IS ORDINARY, NOT AN ERROR.  A results file that travelled
%   without its video is a normal thing to be handed, so .frame comes back empty and
%   .note says why, and the caller draws the walls alone or says so on the axes.
%   Nothing throws.
%
%   THE WALLS ARE A LOWER BOUND WHEN .valid IS FALSE - one wall had dilated out of
%   the field of view that frame, so the diameter under it is not a measurement.
%   Callers draw those walls RED; the interval editor already does.
%
%   Syntax:
%      F = getMyographWallFrame(results,k)
%      F = getMyographWallFrame(results,k,measure,when)
%      F = getMyographWallFrame(results,k,measure,when,recording)
%
%   INPUTS
%     results   the recording's RESULTS struct (myographProduct 'open').
%     k         the window, as myographIntervals numbers the flat list.  EMPTY means
%               the window that CONTAINS when, which is what a caller scrubbing
%               through a recording has - a time, and no window index of its own.
%               A time no window contains still returns its FRAME, with no walls and
%               a note saying so: the recording holds every frame whether or not the
%               diameter step read it, and the nearest window's walls drawn over a
%               frame they do not describe would be worse than none.
%     measure   (optional, default 'mid') which of the three measured diameters the
%               walls are of - 'outer', 'mid' or 'inner'.
%     when      (optional) a time in ABSOLUTE seconds.  EMPTY (the default) means the
%               MIDDLE of the window, as the frame most likely to be representative
%               of it.  Inside a window, the nearest stored frame is the one shown.
%     recording (optional) the recording itself: a path, or a VideoReader ALREADY
%               OPEN on it.  Omitted, results.recording.fName is used.  A caller
%               stepping through frames hands over its own reader rather than have
%               one constructed per frame.
%
%   OUTPUT
%     F        a struct, always assigned:
%                .frame  the image, or EMPTY when the recording could not be read
%                .left   [rows x 1] left  wall position, px
%                .right  [rows x 1] right wall position, px
%                .rows   [rows x 1] the IMAGE ROWS those positions belong to - the
%                        measured band, not 1:height, because the rows outside it
%                        carry an interpolated fill and are not a measurement
%                .time   the frame's absolute time, s (NaN when there is none)
%                .valid  false = a wall had left the field of view
%                .note   why .frame is empty, in one sentence, or '' when it is not
%
%   EXAMPLE
%     F = getMyographWallFrame(results,1,'mid',[]);
%     if ~isempty(F.frame), image(ax,F.frame); axis(ax,'image'); end
%     col = [0.1 0.9 0.1]; if ~F.valid, col = [0.95 0.15 0.15]; end
%     hold(ax,'on');
%     plot(ax,F.left,F.rows,'-','Color',col); plot(ax,F.right,F.rows,'-','Color',col);
%
%   DEPENDS ON
%     myographIntervals, myographMeasureIndex (Core/Myograph) and MATLAB's base
%     VideoReader.
%
% See also: editMyographIntervals, getMyographTrace, runMyographDiameter,
%           myographIntervals, myographProduct, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function F = getMyographWallFrame(results,k,measure,when,recording)

if nargin<2, k=[]; end
if nargin<3 || isempty(measure), measure='mid'; end
if nargin<4, when=[]; end
if nargin<5, recording=[]; end

F=struct('frame',[],'left',[],'right',[],'rows',[],'time',NaN,'valid',true,'note','');

rec=fieldOr(results,'recording',struct());
ivs=myographIntervals(results);
k=windowIndex(ivs,k,when);
fIdx=[];
if isempty(k)
    % THE PICTURE STILL COMES BACK.  A time nobody measured is an ordinary thing to
    % scrub past - the recording holds every frame whether or not the diameter step
    % read it - so the frame is fetched by time and the walls are simply absent.
    F.time=double(when);
    if isempty(F.time) || ~isfinite(F.time)
        F.note='This recording has no measured window to show walls from.';
        return
    end
    F.note='No walls were measured at this time.';
else
    [F,fIdx]=wallsOf(F,ivs(k),rec,results,measure,when);
end

[frame,note]=readOneFrame(recording,rec,fIdx,F.time);
F.frame=frame;
if ~isempty(note), F.note=note; end
end

% =====================================================================
function [F,fIdx]=wallsOf(F,iv,rec,results,measure,when)
%wallsOf  One window's wall positions at the wanted time, and the frame of the
%   ORIGINAL recording that row was read from.
fIdx=[];
d=fieldOr(iv,'diameter',[]);
if ~isstruct(d) || ~isfield(d,'wallL') || isempty(fieldOr(d,'time',[]))
    F.time=double(when);
    if isempty(F.time) || ~isfinite(F.time), F.time=NaN; end
    F.note=['The wall positions of this window were not kept, so only the ' ...
            'recording can be shown.'];
    return
end

[~,mName]=myographMeasureIndex(measure,measuresOf(results));
if ~isfield(d.wallL,mName)
    fn=fieldnames(d.wallL);
    if isempty(fn), F.note='This window carries no wall positions.'; return; end
    mName=fn{1};
end

t=double(d.time(:));
fi=rowAt(t,when);
F.time=t(fi);
F.left =double(d.wallL.(mName)(fi,:)).';
F.right=double(d.wallR.(mName)(fi,:)).';
F.rows =measuredRows(rec,numel(F.left)).';
vv=fieldOr(d,'valid',[]);
if ~isempty(vv) && numel(vv)>=fi, F.valid=logical(vv(fi)); end

% THE FRAME OF THE ORIGINAL RECORDING this row was read from.  .frames is the
% window's range in that recording, and the window's rows run contiguously from its
% first frame, so this is the whole of the arithmetic - and it survives every
% narrowing, because neither term is an index into anything that was cut.
fr=double(fieldOr(iv,'frames',[]));
if numel(fr)==2 && all(isfinite(fr)), fIdx=fr(1)+fi-1; end
end

% =====================================================================
function k=windowIndex(ivs,k,when)
%windowIndex  Which window is being asked for.  A number is taken as given; EMPTY
%   means the window that CONTAINS the time - a caller scrubbing a recording has a
%   time and no index of its own.  NO WINDOW CONTAINS IT IS AN ANSWER: the nearest
%   window's walls would be drawn over a frame they do not describe, which is worse
%   than no walls at all.  With no time either, the first measured window is meant.
if isempty(ivs), k=[]; return; end
if ~isempty(k)
    k=round(double(k(1)));
    if k<1 || k>numel(ivs), k=[]; end
    return
end
for i=1:1:numel(ivs)
    t=double(fieldOr(fieldOr(ivs(i),'diameter',[]),'time',[]));
    if isempty(t), continue; end
    if isempty(when) || ~isfinite(when), k=i; return; end
    if when>=t(1) && when<=t(end), k=i; return; end
end
k=[];
end

% =====================================================================
function fi=rowAt(t,when)
%rowAt  Which stored row of the window is wanted.  No time means the MIDDLE of it -
%   the frame most likely to be representative - and a time takes the nearest row,
%   which for a time inside the window is the frame it names.
if isempty(when) || ~isfinite(when)
    fi=max(1,round(numel(t)/2));
    return
end
[~,fi]=min(abs(t-double(when(1))));
end

% =====================================================================
function [frame,note]=readOneFrame(recording,rec,fIdx,tWanted)
%readOneFrame  The picture itself, from the recording named or handed over.
%   BY FRAME NUMBER WHEN THERE IS ONE, because that is exact; by absolute time
%   otherwise, with the clamp the retired explorer view had right - CurrentTime must
%   stay one frame short of Duration or the read runs off the end of the file.
frame=[];
[vr,own,note]=openRecording(recording,rec);
if isempty(vr), return; end
try
    if ~isempty(fIdx) && isfinite(fIdx)
        n=frameCount(vr);
        frame=read(vr,min(max(1,round(fIdx)),n));
    else
        vr.CurrentTime=min(max(tWanted,0),max(0,vr.Duration-1/max(vr.FrameRate,eps)));
        frame=readFrame(vr);
    end
catch
    frame=[];
    note='That frame could not be read from the recording.';
end
% A READER THIS FUNCTION OPENED IS RELEASED HERE.  The caller's own reader is left
% alone: it is scrubbing with it, and closing it would stop the preview.
if own, delete(vr); end
end

% =====================================================================
function [vr,own,note]=openRecording(recording,rec)
%openRecording  A reader on the recording: the one the caller handed over, or one
%   opened on the path it named.  'own' says whether this function opened it, which
%   is what decides whether it may release it - closing a caller's reader would stop
%   the preview it is scrubbing with.
vr=[]; own=false; note='';
if isa(recording,'VideoReader'), vr=recording; return; end

p='';
if ~isempty(recording), p=char(string(recording)); end
if isempty(p) || ~isfile(p), p=char(string(fieldOr(rec,'fName',''))); end
if isempty(p) || ~isfile(p)
    note='The recording is not beside this result - keep the video with it.';
    return
end
try
    vr=VideoReader(p); own=true;
catch
    vr=[];
    note='The recording is there but could not be opened.';
end
end

% =====================================================================
function n=frameCount(vr)
%frameCount  How many frames the file has, however the container answers.
try
    n=vr.NumFrames;
catch
    n=[];
end
if isempty(n) || ~isfinite(n), n=floor(vr.Duration*vr.FrameRate); end
n=max(1,double(n));
end

% =====================================================================
function rows=measuredRows(rec,nY)
%measuredRows  The image rows the stored wall positions belong to.  The arrays hold
%   the measured band only, so column 1 is image row rowRange(1) and not image row 1;
%   drawing them at 1:nY would put the walls in the wrong place on the frame.
%   A BAND THAT IS NOT A BAND is read as the whole image, which is the same fallback
%   the diameter step makes when it measures: an inverted or out-of-range rowRange
%   never selected rows there either, so offsetting by it here would move the walls
%   off the vessel on exactly the recordings whose settings were wrong.
rows=1:nY;
rr=double(fieldOr(rec,'rowRange',[]));
if numel(rr)~=2 || ~all(isfinite(rr(1))) || rr(2)<rr(1), return; end
r0=max(1,round(rr(1)));
rows=r0:r0+nY-1;
end

% =====================================================================
function meas=measuresOf(results)
%measuresOf  The measures, and their order.  results.recording is the one authority
%   on it, written by the diameter step in the same breath as the windows.
meas=reshape(cellstr(fieldOr(fieldOr(results,'recording',struct()),'measures',{})),1,[]);
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
