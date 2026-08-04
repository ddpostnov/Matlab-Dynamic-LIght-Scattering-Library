%getMyographTrace  The curve a myograph recording's intervals are chosen on
%
%   data = getMyographTrace(results,s) assembles everything editMyographIntervals
%   draws for ONE recording, from that recording's RESULTS struct.  There are three
%   cases and the function decides between them from the data itself:
%
%     A DIAMETER HAS BEEN MEASURED (some window carries a .diameter).  The trace is
%     the LINE-AVERAGED diameter of the analysed measure (s.edgeMode), read back
%     from the windows rather than recomputed, so the curve is exactly the one the
%     diameter step wrote.  The detected walls are the video overlay, and the frames
%     where a wall left the field of view are flagged so they can be drawn red.
%
%     NO DIAMETER YET (the pre-set-intervals and time-crop steps, which run before
%     anything has been measured).  The trace is a PER-FRAME BRIGHTNESS PROFILE read
%     from the recording itself - enough to see where the protocol changed, which is
%     what a window is chosen against at that stage - and there is no wall overlay.
%     The video is then the point of the window.
%
%     A WIRE MYOGRAPH (results.recording.modality 'WMYO').  There is no diameter and
%     no video: the recording IS a set of channels, so they are handed over whole,
%     together with the operator's comments and the LabChart record boundaries, and
%     the editor draws whichever of them the operator ticks.  Nothing is averaged
%     and nothing is resampled here - a channel keeps its own sampling rate.  Once
%     the windows have been chosen the channels hold only those windows' samples,
%     and myographChannelSamples is what reads a channel either way.
%
%   THE WINDOWS ARE LAID END TO END, AND THE GAPS ARE DRAWN AS BREAKS.  A recording
%   measured over disjoint windows has no samples between them, so the returned
%   .time jumps; .fs is the authority on the sampling rate and diff(.time) is not,
%   and the editor breaks the drawn line wherever the two disagree.  THAT IS A BREAK
%   IN A LINE, NOT A NaN IN THE DATA - nothing here fills the gaps, because the rule
%   the product is built on is that what falls outside the windows is not stored.
%   Do not "fix" it by interpolating: the stretch genuinely was never measured.
%
%   THE PROFILE IS A COARSE PRE-PASS, NOT A SCRUB.  It is read once, before the
%   window opens, over every Nth frame, with N set so about s.profileSamples points
%   come back.  The alternative - computing it as the operator scrubs - leaves the
%   trace with holes exactly where nobody has looked yet, and an interval cannot be
%   chosen on a curve that only exists where you have already been.  Measured on a
%   2 h 15 fps recording (124371 frames, 480x640): seeking to one frame costs ~25 ms
%   whatever the distance, decoding one frame in sequence costs ~1.5 ms.  So the pass
%   SEEKS when it wants fewer than one frame in sixteen - 39 s measured on that
%   recording at the default 1200 points, and it scales with the number of POINTS,
%   not with the length of the recording - and DECODES STRAIGHT THROUGH when it wants
%   more (a 150 s recording costs ~2 s).  Sixteen is where the two costs cross.
%
%   Syntax:
%      data = getMyographTrace(results,s)
%      data = getMyographTrace(results,s,rawName)
%      data = getMyographTrace(results,s,rawName,comments)
%
%   INPUTS
%     results  the recording's RESULTS struct (myographProduct 'open').
%     s        parameter structure:
%       .edgeMode       (optional, default 'mid') which of the three diameters the
%                       trace shows - the measure that is analysed by default
%       .profileSamples (optional, default 1200) how many points the brightness
%                       profile has.  Only read when no diameter exists yet.
%     rawName  (optional) the recording to read and to preview.  Defaults to
%              results.recording.fName, then to the video sitting beside the product.
%     comments (optional) the recording's LabChart comments.  They are already in
%              results.comments; the argument is kept so a caller that holds an
%              edited set can pass it instead.
%
%   OUTPUTS
%     data     the struct editMyographIntervals takes: .time .trace .traceName
%              .traceUnit .valid .walls .video .fs, plus
%              .channels .comments .blocks for a wire myograph.  Times are ABSOLUTE
%              seconds from the start of the recording, and .fs is the sampling rate
%              OF THE TRACE (the profile is sampled coarsely, so it is not the frame
%              rate) - which is what breaks the drawn line over the gaps a set of
%              disjoint intervals leaves in the time base.
%              .wallFrame(t,vr) rides with the walls when there are any: ask it for a
%              time and it gives back one FRAME of the recording with that frame's
%              own walls on it (getMyographWallFrame).  The concatenated .walls are
%              drawn against the trace, where a position in the concatenation is all
%              that is wanted; a picture needs the frame of the original file, and
%              only the window's .frames can name that.  Pass an open VideoReader as
%              the second argument to reuse it rather than have one opened per frame.
%              .measurement is true when .trace (or .channels) IS the measurement, so
%              a stretch left outside every window is discarded by saving them, and
%              false when the trace is only a look at the recording - the brightness
%              profile, which costs nothing to leave out.  It is what the editor asks
%              the operator about before Done, and it is asked of the DATA rather
%              than of the step, so a crop chosen after a diameter exists is treated
%              like the narrowing it is.
%
%   DEPENDS ON
%     myographMeasureIndex, myographChannelSamples, getMyographWallFrame
%     (Core/Myograph) and MATLAB's base VideoReader.
%
% See also: editMyographIntervals, setMyographIntervals, setMyographPresetIntervals,
%           setMyographCrop, myographProduct, getMyographDiameter,
%           getMyographWallFrame, myographChannelSamples, readLabChart
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function data = getMyographTrace(results,s,rawName,comments)

if nargin<2, s=struct(); end
if nargin<3, rawName=''; end
if nargin<4, comments=[]; end
rec=fieldOr(results,'recording',struct());
rawName=recordingPath(rec,rawName);

% .measurement SAYS WHETHER LEAVING A STRETCH OUT OF THE WINDOWS COSTS ANYTHING.
% The diameter and the channels ARE the measurement, and the product keeps no copy
% of them, so a window drawn around part of them discards the rest; the brightness
% profile is only a look at the recording, which is still there afterwards.  The
% editor asks before Done on the first and not on the second.
data=struct('time',[],'trace',[],'traceName','','traceUnit','','valid',[], ...
    'walls',[],'video',rawName,'fs',[],'measurement',false);

ch=fieldOr(results,'channel',[]);
measured=measuredWindows(results);
% A WIRE MYOGRAPH IS ONE THAT HAS A CHANNEL AXIS.  Testing whether the channels still
% carry whole-recording samples was right while they were only ever written whole;
% the intervals step now moves them into the windows, and a recording re-opened after
% that would be read as a pressure myograph with nothing measured.
if ~isempty(ch) && isstruct(ch)
    data=channelTraces(data,results,rec,comments);
elseif ~isempty(measured)
    data=diameterTrace(data,measured,rec,s);
else
    data=brightnessTrace(data,rec,s,rawName);
end
end

% =====================================================================
function data=channelTraces(data,results,rec,comments)
%channelTraces  The wire myograph, handed over as it was recorded.  The editor
%   draws whichever channels are ticked, so all of them go across; picking one here
%   would decide for the operator, and averaging them would be meaningless - they
%   are different vessels in different chambers, sometimes in different units.
%
%   THE SAMPLES COME THROUGH myographChannelSamples, which reads a channel whichever
%   way the product holds it - whole before the windows were chosen, per window
%   after.  Once they are per window the time base JUMPS across what was discarded,
%   and the editor already draws that as a break in the line rather than a ramp.
data.channels=results.channel;
for i=1:1:numel(data.channels)
    [y,t]=myographChannelSamples(results.channel(i));
    data.channels(i).data=y;
    data.channels(i).time=t;
end
data.blocks=fieldOr(results,'blocks',zeros(0,2));
data.comments=comments;
if isempty(data.comments), data.comments=fieldOr(results,'comments',[]); end
if isempty(data.comments)
    data.comments=struct('time',{},'text',{},'channel',{});
end
data.video='';                         % there is nothing to preview
data.traceName='Recording';
data.fs=fieldOr(rec,'fs',[]);
data.measurement=true;                 % the samples ARE the measurement
end

% =====================================================================
function data=diameterTrace(data,ivs,rec,s)
%diameterTrace  The measured vessel: its averaged diameter, and the walls behind it,
%   with the windows laid end to end in the order they were measured.  That order is
%   the order their .frames tile, so a time here and a frame index there mean the
%   same sample.
meas=measuresOf(rec);
[~,mName]=myographMeasureIndex(fieldOr(s,'edgeMode','mid'),meas);
if ~isfield(ivs(1).diameter,mName) && ~isempty(meas), mName=meas{1}; end

data.time =partsOf(ivs,@(d) double(reshape(d.time,[],1)));
data.trace=partsOf(ivs,@(d) double(reshape(d.(mName),[],1)));
data.traceName=[measureWord(mName) ' diameter'];
data.traceUnit='diameter (px)';
data.valid=logical(partsOf(ivs,@(d) logical(reshape(fieldOr(d,'valid',[]),[],1))));
if numel(data.valid)~=numel(data.time), data.valid=[]; end

% The walls are carried with THE ROWS THEY BELONG TO: they are drawn on the frame at
% their own image rows, and the rows outside the measured band hold a fill rather
% than a detection, so they are not drawn at all.  A recording measured with
% s.keepArrays false kept no walls, and then there is no overlay to draw.
if all(arrayfun(@(iv) isfield(iv.diameter,'wallL'),ivs))
    L=partsOf(ivs,@(d) double(d.wallL.(mName)));
    R=partsOf(ivs,@(d) double(d.wallR.(mName)));
    data.walls=struct('L',L,'R',R,'rows',measuredRows(rec,size(L,2)));
    % AND THE ONE WAY TO ASK FOR A FRAME OF THE RECORDING WITH ITS OWN WALLS ON IT.
    % The concatenated arrays above are drawn against the trace, where a position in
    % the concatenation is all that is wanted; a picture needs the frame of the
    % ORIGINAL file, and only the window's .frames can name that.  Handing over the
    % question rather than the answer keeps the editor from re-deriving it: it has a
    % time, it asks, and what comes back is the frame and the walls of that frame,
    % which therefore cannot disagree.  It is handed the MEASURED windows in the
    % order they were laid end to end, so the window it finds for a time is the one
    % that stretch of the trace came from.
    R=struct('recording',rec,'intervals',ivs);
    data.wallFrame=@(t,vr) getMyographWallFrame(R,[],mName,t,vr);
end
data.fs=fieldOr(rec,'fs',[]);
data.measurement=true;                 % the windows ARE the measurement
end

% =====================================================================
function data=brightnessTrace(data,rec,s,rawName)
%brightnessTrace  Before anything has been measured: how bright each frame is, over
%   the whole recording.  It is not a measurement of the vessel and does not pretend
%   to be one - it is the cheapest curve that shows where the protocol changed.
nWant=fieldOr(s,'profileSamples',1200);
if isempty(nWant) || ~isfinite(nWant) || nWant<2, nWant=1200; end
[t,y,step,fps]=videoProfile(rawName,nWant);
data.time=t; data.trace=y;
data.traceName='Frame brightness (no diameter measured yet)';
data.traceUnit='brightness (a.u.)';
data.fs=fps/max(step,1);            % the sampling of the TRACE, not of the recording
if isempty(data.fs) || ~isfinite(data.fs) || data.fs<=0
    data.fs=fieldOr(rec,'fs',1);
end
end

% =====================================================================
function ivs=measuredWindows(results)
%measuredWindows  The windows that actually carry a measurement, in the order the
%   diameter step wrote them.  A window with no .diameter is a boundary somebody
%   drew before anything was measured, and there is nothing of it to draw.
ivs=[];
defined=fieldOr(results,'intervals',[]);
if isempty(defined) || ~isstruct(defined), return; end
keep=arrayfun(@(iv) isstruct(fieldOr(iv,'diameter',[])) && ...
    ~isempty(fieldOr(iv.diameter,'time',[])),defined);
if any(keep), ivs=defined(keep); end
end

% =====================================================================
function X=partsOf(ivs,get)
%partsOf  One window's array after another, laid end to end along time.
pieces=cell(1,numel(ivs));
for k=1:1:numel(ivs), pieces{k}=get(ivs(k).diameter); end
X=cat(1,pieces{:});
end

% =====================================================================
function meas=measuresOf(rec)
%measuresOf  The measures, and their order.  results.recording is the one authority
%   on it, written by the diameter step in the same breath as the windows.
meas=reshape(cellstr(fieldOr(rec,'measures',{})),1,[]);
end

% =====================================================================
function [t,y,step,fps]=videoProfile(rawName,nWant)
%videoProfile  The mean brightness of every Nth frame, read the cheaper of the two
%   ways (see the header for the measurement behind the threshold).
v=VideoReader(rawName);
fps=v.FrameRate;
n=frameCount(v);
step=max(1,ceil(n/nWant));
idx=(1:step:n)';
y=nan(numel(idx),1);
if step<16
    % Wanting more than one frame in sixteen, decoding straight through is cheaper
    % than seeking to each one - and short recordings always land here.
    j=1; k=0;
    while j<=numel(idx) && hasFrame(v)
        frame=readFrame(v); k=k+1;
        if k==idx(j), y(j)=frameMean(frame); j=j+1; end
    end
    idx=idx(1:j-1); y=y(1:j-1);
else
    for j=1:1:numel(idx)
        try
            y(j)=frameMean(read(v,idx(j)));
        catch
        end
    end
end
t=(double(idx)-1)/max(fps,eps);
end

% =====================================================================
function m=frameMean(frame)
%frameMean  One number per frame.  A colour frame is averaged over its channels too:
%   the profile says how the scene changed, not what colour it was.
m=mean(double(frame(:)),'omitnan');
end

% =====================================================================
function n=frameCount(v)
%frameCount  The number of frames, from the index when the container has one and
%   from the duration when it does not (a variable-frame-rate container may not).
n=[];
try
    n=v.NumFrames;
catch
end
if isempty(n) || ~isfinite(n), n=floor(v.Duration*v.FrameRate); end
n=max(double(n),0);
end

% =====================================================================
function raw=recordingPath(rec,raw)
%recordingPath  THE RECORDING this product was made from: the path the caller passed,
%   or the one the entry step recorded.
if ~isempty(raw) && isfile(raw), return; end
raw=char(fieldOr(rec,'fName',''));
if ~isempty(raw) && isfile(raw), return; end
raw='';
end

% =====================================================================
function rows=measuredRows(rec,nY)
%measuredRows  The image rows the vessel was measured on.  The stored walls already
%   cover those rows and no others, so this only names them for the overlay - and it
%   falls back to counting them when the recording's row band no longer agrees with
%   what is stored, because the array on disk is the authority on its own width.
rows=1:nY;
rowRange=fieldOr(rec,'rowRange',[]);
sz=fieldOr(rec,'size',[]);
if numel(rowRange)~=2 || numel(sz)<1, return; end
r0=max(1,round(rowRange(1)));
r1=min(double(sz(1)),round(rowRange(2)));
if r1>=r0 && (r1-r0+1)==nY, rows=r0:r1; end
end

% =====================================================================
function w=measureWord(name)
%measureWord  The measure, as a person says it rather than as the field is spelled.
switch lower(char(name))
    case 'outer', w='Outer-wall';
    case 'inner', w='Luminal';
    otherwise,    w='Wall-centre';
end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
