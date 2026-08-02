%getMyographTrace  The curve a myograph recording's intervals are chosen on
%
%   data = getMyographTrace(source,s) assembles everything editMyographIntervals
%   draws for ONE recording, from that recording's SOURCE struct.  There are two
%   cases and the function decides between them from the data itself:
%
%     A DIAMETER HAS BEEN MEASURED (source.data is filled).  The trace is the
%     LINE-AVERAGED diameter of the analysed measure (s.edgeMode), over the measured
%     rows only - the rows outside source.rowRange were never measured and carry an
%     interpolated fill, and averaging them in would drag the trace towards the edge
%     of the vessel.  The detected walls are the video overlay, and the frames where a
%     wall left the field of view are flagged so they can be drawn red.
%
%     NO DIAMETER YET (the pre-set-intervals and time-crop steps, which run before
%     anything has been measured).  The trace is a PER-FRAME BRIGHTNESS PROFILE read
%     from the recording itself - enough to see where the protocol changed, which is
%     what a window is chosen against at that stage - and there is no wall overlay.
%     The video is then the point of the window.
%
%     A WIRE MYOGRAPH (source.modality 'WMYO').  There is no diameter and no video:
%     the recording IS a set of channels, so they are handed over whole, together
%     with the operator's comments and the LabChart record boundaries, and the
%     editor draws whichever of them the operator ticks.  Nothing is averaged and
%     nothing is resampled here - a channel keeps its own sampling rate.
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
%      data = getMyographTrace(source,s)
%      data = getMyographTrace(source,s,rawName)
%      data = getMyographTrace(source,s,rawName,comments)
%
%   INPUTS
%     source   the recording's SOURCE struct (myographProduct 'open').
%     s        parameter structure:
%       .edgeMode       (optional, default 'mid') which of the three diameters the
%                       trace shows - the measure that is analysed by default
%       .profileSamples (optional, default 1200) how many points the brightness
%                       profile has.  Only read when no diameter exists yet.
%     rawName  (optional) the recording to read and to preview.  Defaults to
%              source.fName, then to the video sitting beside the product.
%     comments (optional) the recording's LabChart comments (results.comments).
%              They live in RESULTS rather than in SOURCE - they are something the
%              operator wrote, not something that was sampled - so the caller that
%              has both is the one that hands them over.
%
%   OUTPUTS
%     data     the struct editMyographIntervals takes: .time .trace .traceName
%              .traceUnit .valid .walls .video .fs, plus
%              .channels .comments .blocks for a wire myograph.  Times are ABSOLUTE
%              seconds from the start of the recording, and .fs is the sampling rate
%              OF THE TRACE (the profile is sampled coarsely, so it is not the frame
%              rate) - which is what breaks the drawn line over the gaps a set of
%              disjoint intervals leaves in the time base.
%
%   DEPENDS ON
%     myographMeasureIndex (Core/Myograph) and MATLAB's base VideoReader.
%
% See also: editMyographIntervals, setMyographIntervals, setMyographPresetIntervals,
%           setMyographCrop, myographProduct, getMyographDiameter, readLabChart
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function data = getMyographTrace(source,s,rawName,comments)

if nargin<2, s=struct(); end
if nargin<3, rawName=''; end
if nargin<4, comments=[]; end
rawName=recordingPath(source,rawName);

data=struct('time',[],'trace',[],'traceName','','traceUnit','','valid',[], ...
    'walls',[],'video',rawName,'fs',[]);

if ~isempty(fieldOr(source,'channels',[]))
    data=channelTraces(data,source,comments);
elseif ~isempty(fieldOr(source,'data',[]))
    data=diameterTrace(data,source,s);
else
    data=brightnessTrace(data,source,s,rawName);
end
end

% =====================================================================
function data=channelTraces(data,source,comments)
%channelTraces  The wire myograph, handed over as it was recorded.  The editor
%   draws whichever channels are ticked, so all of them go across; picking one here
%   would decide for the operator, and averaging them would be meaningless - they
%   are different vessels in different chambers, sometimes in different units.
data.channels=source.channels;
data.blocks=fieldOr(source,'blocks',zeros(0,2));
data.comments=comments;
if isempty(data.comments)
    data.comments=struct('time',{},'text',{},'channel',{});
end
data.video='';                         % there is nothing to preview
data.traceName='Recording';
data.fs=fieldOr(source,'fs',[]);
end

% =====================================================================
function data=diameterTrace(data,source,s)
%diameterTrace  The measured vessel: its averaged diameter, and the walls behind it.
meas=fieldOr(source,'measures',{});
[k,mName]=myographMeasureIndex(fieldOr(s,'edgeMode','mid'),meas);
k=min(k,size(source.data,3));
rows=measuredRows(fieldOr(source,'rowRange',[1 Inf]),size(source.data,2));

data.time=double(source.time(:));
data.trace=mean(double(source.data(:,rows,k)),2,'omitnan');
data.traceName=[measureWord(mName) ' diameter'];
data.traceUnit='diameter (px)';
data.valid=fieldOr(source,'valid',[]);
% The walls are carried with THE ROWS THEY BELONG TO: they are drawn on the frame at
% their own image rows, and the rows outside the measured band hold a fill rather
% than a detection, so they are not drawn at all.
wl=fieldOr(source,'wallL',[]); wr=fieldOr(source,'wallR',[]);
if ~isempty(wl) && ~isempty(wr)
    data.walls=struct('L',double(wl(:,rows,k)),'R',double(wr(:,rows,k)),'rows',rows);
end
data.fs=fieldOr(source,'fs',[]);
end

% =====================================================================
function data=brightnessTrace(data,source,s,rawName)
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
    data.fs=fieldOr(source,'fs',1);
end
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
function raw=recordingPath(source,raw)
%recordingPath  THE RECORDING this product was made from: the path the caller passed,
%   the one the entry step recorded, or the video sitting beside the product.
if ~isempty(raw) && isfile(raw), return; end
raw=char(fieldOr(source,'fName',''));
if ~isempty(raw) && isfile(raw), return; end
raw='';
end

% =====================================================================
function rows=measuredRows(rowRange,nY)
%measuredRows  The image rows the vessel was measured on.  Rows outside them carry
%   an interpolated fill, so nothing that averages along the vessel may include them.
if numel(rowRange)~=2, rowRange=[1 Inf]; end
r0=max(1,round(rowRange(1)));
r1=min(nY,round(rowRange(2)));
if ~(r1>=r0), rows=1:nY; else, rows=r0:r1; end
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
