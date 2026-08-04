%OBSOLETEcutMyographSource  RETIRED 2026-08-03 - superseded by cutMyographIntervals
%
%   It cut a myograph SOURCE struct, and a myograph product no longer has one: the
%   measurement lives per window in results.intervals(k).diameter, and the span this
%   used to slice is assembled from those windows.  Its replacement is
%   Core/Myograph/cutMyographIntervals, which takes the RESULTS instead and builds
%   the diameter block through the shared myographDiameterBranch.
%
%   Kept for reference only.  Nothing calls it and nothing should.
%
% =====================================================================
%cutMyographSource  Cut a measured myograph recording into its analysis intervals
%
%   intervals = cutMyographSource(source,ivT,names) turns a set of time windows into
%   the results.intervals elements the myograph result tree holds: each window's frame
%   range into source.time, its per-line diameter, its line-AVERAGED traces and the
%   statistics a protocol is compared on.
%
%   THE SHAPE IS runMyographDiameter'S, BECAUSE IT IS THE SAME THING - the diameter of
%   this window - and moving a boundary must not change what a window contains.  So
%   .diameter.lines.<measure> is written here too: a window rebuilt after an edit that
%   carried only the trace would silently drop the measurement the trace was made of,
%   halfway down a pipeline that had it a moment earlier.  Everything is taken over
%   the MEASURED rows only (source.rowRange): the rows outside were never measured and
%   carry an interpolated fill, and averaging them in would drag the trace towards the
%   edge of the vessel.
%
%   (runMyographDiameter can be told not to keep the arrays, with s.keepLines false.
%   That setting is not visible from here - this function is handed a source and a set
%   of boundaries, nothing else - and the wrong way to guess is the one that loses a
%   measurement, so they are always written.)
%
%   EVERY ELEMENT IS REBUILT FROM THE BOUNDARIES AS THEY NOW STAND, and no attempt is
%   made to match a window against a previously stored one.  A window still called
%   'drug' is a DIFFERENT window once a boundary has moved, and nothing outside these
%   arguments can tell which of two same-named windows an older analysis belonged to;
%   rebuilding is the only answer that cannot pair a trace with the wrong window.  For
%   the same reason .propagation and .vasomotion come back EMPTY - they describe a
%   window, so a window that moved has none until those steps run again.
%
%   NO WINDOWS MEANS THE WHOLE RECORDING, as one interval named for what it is - the
%   same rule setRegions applies when no ROI is drawn, so the analyses downstream
%   always have a window to run in.
%
%   A WIRE MYOGRAPH IS CUT BY TIME, NOT BY FRAME.  It has no diameter and its
%   channels may be sampled at different rates, so an interval carries its
%   boundaries and the CHANNEL it was selected on, and each analysis picks its own
%   samples out of its own channel by time.  .frames is still filled when the
%   channels happen to share one base (source.time is then real), and is left empty
%   when they do not - there is no single index range that would be true for all of
%   them, and a wrong one would be worse than none.
%
%   AND IT COMES BACK SPLIT BY CHANNEL, because that is what a wire myograph result
%   is: one LabChart file is several chambers, and a window may belong to one of
%   them.  The second output is the results.channel axis -
%       channels(i).name  channels(i).intervals(k)
%   - with a window assigned to one channel appearing under that channel only, and a
%   window left on 'all channels' appearing under EVERY channel, which is what it
%   means.  Its own .channels stays empty, so re-opening the recording in the editor
%   tells the two apart again instead of finding one window per channel.
%   A PRESSURE MYOGRAPH FILLS THE FIRST OUTPUT AND LEAVES THE SECOND EMPTY; a wire
%   one does the opposite.  Both are always assigned, so the caller stores the pair
%   and never branches on the modality.
%
%   Syntax:
%      intervals             = cutMyographSource(source,ivT,names)
%      [intervals,channels]  = cutMyographSource(source,ivT,names,channels)
%
%   INPUTS
%     source   the recording's SOURCE struct after the diameter step: .time, .data
%              [T x nY x 3], .mask, .valid, .measures, .rowRange.  For a wire
%              myograph: .channels, and .time when the rates agree.
%     ivT      [n x 2] window start/end in ABSOLUTE seconds.  A window is clipped to
%              the analysed span; one that falls entirely outside it is kept (the
%              operator defined it) with no frames and no diameter (there is nothing
%              in it).
%     names    1 x n cellstr of window names.  Missing ones are auto-named.
%     channels (wire myograph) 1 x n cell of cellstr - the channel each window was
%              selected on.  {} means every channel; one name means that channel.
%
%   OUTPUTS
%     intervals  1 x n struct array with the interval fields of the result tree:
%                .name .tStart .tEnd .frames .channels .diameter .propagation
%                .vasomotion (the last two empty).  EMPTY for a wire myograph.
%     channels   1 x m struct array .name .intervals - the results.channel axis.
%                EMPTY for a pressure myograph.
%
%   DEPENDS ON
%     myographProduct (Core/Myograph) for the interval prototype.
%
% See also: setMyographIntervals, runMyographDiameter, myographProduct,
%           editMyographIntervals, runLabChart
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function [intervals,channels] = cutMyographSource(source,ivT,names,chans)

if nargin<2 || isempty(ivT), ivT=zeros(0,2); end
if nargin<3, names={}; end
if nargin<4, chans={}; end
names=cellstr(names);
chans=reshape(chans,1,[]);

channels=emptyChannels();
intervals=cutWindows(source,ivT,names,chans);
if isempty(fieldOr(source,'channels',[])), return; end

% A WIRE MYOGRAPH ANSWERS PER CHANNEL.  The windows were cut once, above; splitting
% them is a matter of which channel each belongs to, so nothing is re-cut and the
% element under two channels is the same element.
channels=splitByChannel(intervals,source.channels);
intervals=myographProduct('intervals');
end

% =====================================================================
function channels=splitByChannel(intervals,ch)
%splitByChannel  The results.channel axis: every channel of the recording, with the
%   windows that are analysed on it.  A channel with no window of its own still
%   appears - it is part of the recording, and an empty list says "nothing was
%   defined here" where a missing element would say "there is no such channel".
channels=emptyChannels();
for i=1:1:numel(ch)
    nm=char(fieldOr(ch(i),'name',sprintf('channel%d',i)));
    keep=false(1,numel(intervals));
    for k=1:1:numel(intervals)
        want=intervals(k).channels;
        keep(k)=isempty(want) || any(strcmpi(strtrim(reshape(cellstr(want),1,[])),strtrim(nm)));
    end
    channels(i).name=nm;
    channels(i).intervals=intervals(keep);
end
end

% =====================================================================
function c=emptyChannels()
c=struct('name',{},'intervals',{});
end

% =====================================================================
function intervals=cutWindows(source,ivT,names,channels)
%cutWindows  The windows themselves, before anything is split - one element per
%   window, in the order the editor handed them over.

intervals=myographProduct('intervals');
t=double(fieldOr(source,'time',[])); t=t(:);
wire=~isempty(fieldOr(source,'channels',[]));
if wire
    rows=[]; meas={};
else
    rows=measuredRows(fieldOr(source,'rowRange',[1 Inf]),size(source.data,2));
    meas=fieldOr(source,'measures',{'outer','mid','inner'});
end

for k=1:1:size(ivT,1)
    if numel(names)>=k && ~isempty(names{k})
        intervals(k).name=char(names{k});
    else
        intervals(k).name=sprintf('interval%d',k);
    end
    intervals(k).tStart=ivT(k,1);
    intervals(k).tEnd=ivT(k,2);
    intervals(k).channels=channelsOf(channels,k);
    idx=find(t>=ivT(k,1) & t<=ivT(k,2));
    if isempty(idx)
        intervals(k).frames=[]; intervals(k).diameter=[];
    else
        intervals(k).frames=[idx(1) idx(end)];
        if wire, intervals(k).diameter=[];
        else,    intervals(k).diameter=diameterBranch(source,idx(1):idx(end),rows,meas);
        end
    end
    intervals(k).propagation=[];
    intervals(k).vasomotion=[];
end

if ~isempty(intervals), return; end
sp=analysedSpan(source,t,wire);
if isempty(sp), return; end
intervals(1).name='whole recording';
intervals(1).tStart=sp(1); intervals(1).tEnd=sp(2);
intervals(1).channels=channelsOf(channels,1);
if wire && isempty(t)
    intervals(1).frames=[]; intervals(1).diameter=[];
elseif wire
    intervals(1).frames=[1 numel(t)]; intervals(1).diameter=[];
else
    intervals(1).frames=[1 numel(t)];
    intervals(1).diameter=diameterBranch(source,1:numel(t),rows,meas);
end
intervals(1).propagation=[]; intervals(1).vasomotion=[];
end

% =====================================================================
function sp=analysedSpan(source,t,wire)
%analysedSpan  The whole recording, as one window.  For a wire myograph the span is
%   the union over the channels: they may differ in length because they may differ
%   in sampling rate, and a window cut to the shortest of them would quietly drop
%   the end of every other one.
sp=[];
if ~isempty(t), sp=[t(1) t(end)]; if ~wire, return; end, end
if ~wire, return; end
lo=inf; hi=-inf;
ch=source.channels;
for i=1:1:numel(ch)
    tt=double(ch(i).time);
    if isempty(tt), continue; end
    lo=min(lo,tt(1)); hi=max(hi,tt(end));
end
if isfinite(lo) && isfinite(hi), sp=[lo hi]; end
end

% =====================================================================
function c=channelsOf(channels,k)
%channelsOf  The channels this window was selected on, always as a 1xN cellstr.
c={};
if numel(channels)>=k && ~isempty(channels{k}), c=reshape(cellstr(channels{k}),1,[]); end
end

% =====================================================================
function d=diameterBranch(source,idx,rows,meas)
%diameterBranch  One window's own diameter: every line's own over the rows that were
%   really measured, the same numbers averaged along the vessel, and the statistics a
%   protocol is compared on.  Same shape runMyographDiameter writes, because it is the
%   same thing - the diameter of this window - and every consumer reads it by that
%   shape.  The field ORDER matches it too, so a window cut here and a window measured
%   there are the same struct and not merely equivalent ones.
d=struct('time',[],'nY',numel(rows),'lines',struct(),'stats',struct());
d.time=double(source.time(idx)); d.time=d.time(:);
valid=fieldOr(source,'valid',[]);
if isempty(valid), valid=true(numel(idx),1); else, valid=logical(valid(idx)); valid=valid(:); end
mask=[];
if ~isempty(fieldOr(source,'mask',[])), mask=source.mask(idx,rows); end
for m=1:1:numel(meas)
    lines=single(source.data(idx,rows,m));      % source.data is already single
    d.lines.(meas{m})=lines;
    trace=mean(double(lines),2,'omitnan');
    d.(meas{m})=trace;
    d.stats.(meas{m})=traceStats(trace,valid,mask);
end
end

% =====================================================================
function st=traceStats(trace,valid,mask)
%traceStats  What a protocol compares: the level, the swing and how much of it was
%   actually seen.  OFF-FOV FRAMES ARE EXCLUDED from the level and the swing - a wall
%   that dilated out of view leaves an edge-clamped lower bound behind, and averaging
%   it in would report a constriction that never happened.  How many such frames
%   there were is not hidden: it is validFraction.
if isempty(valid) || numel(valid)~=numel(trace), valid=true(size(trace)); end
v=trace(valid); v=v(isfinite(v));
st=struct('mean',NaN,'std',NaN,'min',NaN,'max',NaN,'amplitude',NaN, ...
    'measuredFraction',NaN,'validFraction',NaN);
if ~isempty(v)
    st.mean=mean(v);  st.std=std(v);
    st.min=min(v);    st.max=max(v);
    st.amplitude=st.max-st.min;
end
if ~isempty(mask), st.measuredFraction=mean(double(mask(:))); end
if ~isempty(trace), st.validFraction=mean(double(valid)); end
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
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
