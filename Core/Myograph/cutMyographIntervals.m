%cutMyographIntervals  Re-cut a myograph result into a new set of analysis windows
%
%   intervals = cutMyographIntervals(results,ivT,names) turns a set of time windows
%   into the results.intervals elements the myograph result tree holds: each
%   window's frame range in the ORIGINAL RECORDING, its per-line diameter and wall
%   positions, its line-AVERAGED traces and the statistics a protocol is compared on.
%
%   IT CUTS THE RESULT, NOT A RECORDING.  A myograph product has no SOURCE member -
%   the measurement lives per window in results.intervals(k).diameter - so the span
%   this works on is the windows already stored, laid end to end in the order the
%   diameter step measured them.
%
%   WHAT IT LEAVES OUT IS GONE.  The windows it returns are the whole of the
%   measurement from then on: the caller replaces results.intervals with them, and
%   the frames between them are not written anywhere, not as a hole and not as a run
%   of NaN.  That is the point of it and it cannot be undone from the product - only
%   by measuring the recording again.  The two wrappers that call it say so before
%   they save, and refuse when the recording that would answer a re-measure is not
%   beside the results.
%
%   THE SHAPE IS runMyographDiameter'S, BECAUSE IT IS THE SAME THING - the diameter
%   of this window - and moving a boundary must not change what a window contains.
%   Both producers call myographDiameterBranch, so the fields and their ORDER are
%   decided once: a window cut here and a window measured there are the same struct
%   and not merely equivalent ones.
%
%   IT CAN ONLY EVER NARROW.  A boundary outside every stored window has no data
%   behind it - the frames it names were never measured, or were discarded with the
%   windows they sat outside - so such a window comes back with no frames and no
%   diameter, and says so by name.  Widening a window means measuring the recording
%   again.
%
%   .frames IS THE WINDOW'S RANGE IN THE ORIGINAL RECORDING, and it is carried
%   through every re-cut rather than recomputed from a position in the stored data:
%   the stored data is what keeps shrinking, and the recording does not.  A window
%   that reaches across two measured stretches names the frames from the first of
%   them to the last, so its arrays are SHORTER than the range it names - which is
%   the honest reading of a window somebody drew over a gap, and the only one that
%   still points at the right frames of the '.avi'.
%
%   A WINDOW WRITTEN WITHOUT ITS ARRAYS (runMyographDiameter with s.keepArrays
%   false) is narrowed by its TRACE, and the result carries the trace alone, as it
%   did before.  Cutting rows out of a trace and pooling the cut rows of the array
%   are the same arithmetic on the same numbers, so nothing disagrees - what is lost
%   is what was never written, and .stats.<measure>.measuredFraction says so by
%   coming back empty.
%
%   EVERY ELEMENT IS REBUILT FROM THE BOUNDARIES AS THEY NOW STAND, and no attempt
%   is made to match a window against a previously stored one.  A window still
%   called 'drug' is a DIFFERENT window once a boundary has moved, and nothing
%   outside these arguments can tell which of two same-named windows an older
%   analysis belonged to; rebuilding is the only answer that cannot pair a trace
%   with the wrong window.  For the same reason .propagation and .vasomotion come
%   back EMPTY - they describe a window, so a window that moved has none until those
%   steps run again.
%
%   NO WINDOWS MEANS THE WHOLE ANALYSED SPAN, as one interval named for what it is -
%   the same rule setRegions applies when no ROI is drawn, so the analyses
%   downstream always have a window to run in.
%
%   A WIRE MYOGRAPH IS CUT BY TIME, NOT BY FRAME.  It has no diameter and its
%   channels may be sampled at different rates, so an interval carries its
%   boundaries and the CHANNEL it was selected on, and each analysis picks its own
%   samples out of its own channel.  .frames is still filled when the channels share
%   one rate, and is left empty when they do not - there is no single index range
%   that would be true for all of them, and a wrong one would be worse than none.
%
%   AND ITS SAMPLES ARE CUT WITH IT.  runLabChart writes each channel whole because
%   nothing has been chosen yet; this function moves the samples inside each window
%   into that window - results.channel(i).intervals(k).samples, .data and .time at
%   that channel's own rate - and leaves the channel's whole-recording .data and
%   .time EMPTY.  A wire recording's samples are its measurement, so rule 3 applies
%   to them exactly as it applies to a diameter: what falls outside the analysed
%   windows is discarded, and re-reading the '.adicht' is what brings it back.  Read
%   a channel through myographChannelSamples, which knows both shapes.
%
%   AND IT COMES BACK SPLIT BY CHANNEL, because that is what a wire myograph result
%   is: one LabChart file is several chambers, and a window may belong to one of
%   them.  The second output is the results.channel axis -
%       channels(i).name .units .fs .data .time .intervals(k)
%   - with a window assigned to one channel appearing under that channel only, and a
%   window left on 'all channels' appearing under EVERY channel, which is what it
%   means.  Its own .channels stays empty, so re-opening the recording in the editor
%   tells the two apart again instead of finding one window per channel.
%   A PRESSURE MYOGRAPH FILLS THE FIRST OUTPUT AND LEAVES THE SECOND EMPTY; a wire
%   one does the opposite.  Both are always assigned, so the caller stores the pair
%   and never branches on the modality.
%
%   Syntax:
%      intervals                    = cutMyographIntervals(results,ivT,names)
%      [intervals,channels]         = cutMyographIntervals(results,ivT,names,chans)
%      [intervals,channels,tally]   = cutMyographIntervals(results,ivT,names,chans)
%
%   INPUTS
%     results  the recording's RESULTS struct (myographProduct 'open'): .recording,
%              .intervals for a pressure myograph, .channel for a wire one.
%     ivT      [n x 2] window start/end in ABSOLUTE seconds.  A window is clipped to
%              the analysed span; one that falls entirely outside it is kept (the
%              operator defined it) with no frames and no data, and is named in a
%              warning.
%     names    1 x n cellstr of window names.  Missing ones are auto-named.
%     chans    (wire myograph) 1 x n cell of cellstr - the channel each window was
%              selected on.  {} means every channel; one name means that channel.
%
%   OUTPUTS
%     intervals  1 x n struct array with the interval fields of the result tree:
%                .name .tStart .tEnd .frames .channels .diameter .samples
%                .propagation .vasomotion (the last two empty).  EMPTY for a wire
%                myograph.
%     channels   1 x m struct array .name .units .fs .data .time .intervals - the
%                results.channel axis, with .data and .time emptied and the samples
%                moved into the windows.  EMPTY for a pressure myograph.
%     tally      WHAT THIS CUT WOULD THROW AWAY, so the caller can say it and can
%                refuse: .before is how much of the recording the product held,
%                .kept is how much the new windows hold, and .unit is the word for
%                one of them ('frames' for a measured vessel, 'samples' for a wire
%                recording).  Counted on what EXISTS, not on the seconds a window
%                covers - a stretch that was never measured is not lost by being
%                left out - and counted per channel, so a window on every channel
%                counts once on each rather than several times over.
%
%   NOTHING IS WRITTEN HERE.  The cut is returned, so a caller can weigh the tally,
%   refuse, and leave the product exactly as it found it.
%
%   DEPENDS ON
%     myographProduct, myographDiameterBranch, myographChannelSamples (Core/Myograph).
%
% See also: setMyographIntervals, setMyographCrop, myographCropWindow,
%           runMyographDiameter, myographDiameterBranch, myographChannelSamples,
%           myographProduct, editMyographIntervals, runLabChart
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function [intervals,channels,tally] = cutMyographIntervals(results,ivT,names,chans)

if nargin<2 || isempty(ivT), ivT=zeros(0,2); end
if nargin<3, names={}; end
if nargin<4, chans={}; end
names=cellstr(names);
chans=reshape(chans,1,[]);

channels=myographProduct('channels');
ch=fieldOr(results,'channel',[]);
% A WIRE MYOGRAPH IS ONE THAT HAS A CHANNEL AXIS, not one whose channels still carry
% whole-recording samples.  Testing the samples was right while they were only ever
% written whole; now that this function empties them, it would read a recording that
% has been through here once as a pressure myograph the second time.
wire=~isempty(ch) && isstruct(ch);

if wire
    intervals=cutWire(ch,ivT,names,chans);
    % A WIRE MYOGRAPH ANSWERS PER CHANNEL.  The windows were defined once, above;
    % splitting them is a matter of which channel each belongs to, and the samples
    % each one keeps come out of the channel it is stored under.
    [channels,tally]=splitByChannel(intervals,ch);
    intervals=myographProduct('intervals');
    return
end

[intervals,tally]=cutPressure(results,ivT,names,chans);
end

% =====================================================================
function [intervals,tally]=cutPressure(results,ivT,names,chans)
%cutPressure  The measured vessel: the stored windows laid end to end, then cut with
%   the new boundaries.
intervals=myographProduct('intervals');
S=storedSpan(results);
tally=struct('before',numel(S.time),'kept',0,'unit','frames');
cov=false(numel(S.time),1);

for k=1:1:size(ivT,1)
    intervals(k)=blankWindow(nameOf(names,k),ivT(k,1),ivT(k,2),channelsOf(chans,k));
    idx=find(S.time>=ivT(k,1) & S.time<=ivT(k,2));
    if isempty(idx)
        % NOTHING WAS EVER MEASURED HERE.  Said by name rather than returned as a
        % silently empty window: after the recording has been cut to its windows, a
        % boundary outside them is a request to widen, and widening is not something
        % a result can answer.
        warning('cutMyographIntervals:outsideSpan', ...
            ['Window "%s" (%.3g-%.3g s) falls outside every measured window of this ' ...
             'recording, so it has no data behind it.  Measure the recording again ' ...
             'to analyse that stretch.'],intervals(k).name,ivT(k,1),ivT(k,2));
        continue
    end
    intervals(k).frames=frameSpan(S.frame,idx);
    intervals(k).diameter=sliceSpan(S,idx(1):idx(end));
    cov(idx(1):idx(end))=true;
end

if ~isempty(intervals) || isempty(S.time), tally.kept=nnz(cov); return; end
% NO WINDOWS MEANS THE WHOLE ANALYSED SPAN, as one interval named for what it is.
intervals(1)=blankWindow('whole recording',S.time(1),S.time(end),channelsOf(chans,1));
intervals(1).frames=frameSpan(S.frame,1:numel(S.time));
intervals(1).diameter=sliceSpan(S,1:numel(S.time));
tally.kept=numel(S.time);
end

% =====================================================================
function S=storedSpan(results)
%storedSpan  THE ANALYSED SPAN, assembled from the windows that hold it, with the
%   ORIGINAL RECORDING'S FRAME NUMBER carried alongside every row.  That number is
%   the one thing a re-cut cannot re-derive: the stored rows know their absolute
%   time, and only the window they came from knows which frame of the file that was.
%   A window with no diameter contributes nothing: it is a boundary somebody drew
%   before anything was measured.
S=struct('time',[],'frame',[],'nY',0,'measures',{{}},'lines',{{}},'wallL',{{}}, ...
    'wallR',{{}},'measured',[],'valid',[],'hasArrays',true);
ivs=fieldOr(results,'intervals',[]);
if isempty(ivs) || ~isstruct(ivs), return; end

keep=arrayfun(@(iv) isstruct(fieldOr(iv,'diameter',[])),ivs);
ivs=ivs(keep);
if isempty(ivs), return; end

S.measures=measuresOf(results);
S.nY=double(fieldOr(ivs(1).diameter,'nY',0));
S.hasArrays=all(arrayfun(@(iv) isfield(iv.diameter,'lines'),ivs));

nM=numel(S.measures);
[S.lines,S.wallL,S.wallR]=deal(cell(1,nM));
parts=cell(1,numel(ivs));
for k=1:1:numel(ivs), parts{k}=ivs(k).diameter; end

S.time=cellCat(parts,@(d) double(fieldOr(d,'time',[])));
S.frame=frameNumbers(ivs);
S.valid=logical(cellCat(parts,@(d) logical(reshape(fieldOr(d,'valid',[]),[],1))));
if S.hasArrays
    S.measured=logical(cellCat(parts,@(d) logical(fieldOr(d,'measured',[]))));
end
for m=1:1:nM
    nm=S.measures{m};
    if S.hasArrays
        S.lines{m}=cellCat(parts,@(d) d.lines.(nm));
        S.wallL{m}=cellCat(parts,@(d) d.wallL.(nm));
        S.wallR{m}=cellCat(parts,@(d) d.wallR.(nm));
    else
        % THE TRACE IS THE ONLY FORM THE QUANTITY HAS IN THIS FILE, so it is what is
        % narrowed.  myographDiameterBranch pools one column to itself, exactly.
        S.lines{m}=cellCat(parts,@(d) double(reshape(d.(nm),[],1)));
    end
end
end

% =====================================================================
function f=frameNumbers(ivs)
%frameNumbers  The recording's own frame number for every row of the assembled span.
%   A window states its first frame in .frames and its rows are contiguous from
%   there, which is what getMyographDiameter guaranteed when it cut them out of the
%   file's frame axis.  IT IS NOT RE-DERIVED FROM THE TIMES, and it must not be: a
%   window with no frame range is one nothing could say the range of - a wire window
%   whose channels differ in rate is the case - and re-deriving it would put a
%   confident number where the honest answer is none.  Such a window contributes NaN,
%   so a cut landing on it comes back with no frame range rather than a wrong one.
pieces=cell(1,numel(ivs));
for k=1:1:numel(ivs)
    n=numel(double(fieldOr(ivs(k).diameter,'time',[])));
    fr=double(fieldOr(ivs(k),'frames',[]));
    if numel(fr)==2 && isfinite(fr(1)), f0=fr(1); else, f0=NaN; end
    pieces{k}=f0+(0:n-1)';
end
f=cat(1,pieces{:});
end

% =====================================================================
function fr=frameSpan(frames,idx)
%frameSpan  The frame range of a cut window, in the ORIGINAL recording.  Empty when
%   the rows it covers cannot say which frames they were - a window whose range
%   cannot be trusted must not carry one, because the only thing that reads it opens
%   the video with it.
fr=[];
if isempty(frames) || isempty(idx), return; end
f=[frames(idx(1)) frames(idx(end))];
if all(isfinite(f)), fr=f; end
end

% =====================================================================
function d=sliceSpan(S,idx)
%sliceSpan  One window's own block, out of the assembled span.
[lines,wl,wr]=deal(cell(1,numel(S.measures)));
for m=1:1:numel(S.measures)
    lines{m}=S.lines{m}(idx,:);
    if S.hasArrays
        wl{m}=S.wallL{m}(idx,:);
        wr{m}=S.wallR{m}(idx,:);
    end
end
measured=[];
if S.hasArrays && ~isempty(S.measured), measured=S.measured(idx,:); end
d=myographDiameterBranch(S.time(idx),S.nY,S.measures,lines,wl,wr, ...
    measured,S.valid(idx),S.hasArrays);
end

% =====================================================================
function intervals=cutWire(ch,ivT,names,chans)
%cutWire  A wire recording's windows: boundaries and the channel they were selected
%   on, and no diameter.  Their samples and their .frames are filled per channel by
%   splitByChannel, which is where the channel each one is stored under is decided.
intervals=myographProduct('intervals');
for k=1:1:size(ivT,1)
    intervals(k)=blankWindow(nameOf(names,k),ivT(k,1),ivT(k,2),channelsOf(chans,k));
end
if ~isempty(intervals), return; end

% NO WINDOWS MEANS THE WHOLE RECORDING, as one interval named for what it is.
sp=channelSpan(ch);
if isempty(sp), return; end
intervals(1)=blankWindow('whole recording',sp(1),sp(2),channelsOf(chans,1));
end

% =====================================================================
function sp=channelSpan(ch)
%channelSpan  The whole recording, as one window.  The span is the UNION over the
%   channels: they may differ in length because they may differ in sampling rate,
%   and a window cut to the shortest of them would quietly drop the end of every
%   other one.
sp=[]; lo=inf; hi=-inf;
for i=1:1:numel(ch)
    [~,tt]=myographChannelSamples(ch(i));
    if isempty(tt), continue; end
    lo=min(lo,tt(1)); hi=max(hi,tt(end));
end
if isfinite(lo) && isfinite(hi), sp=[lo hi]; end
end

% =====================================================================
function [channels,tally]=splitByChannel(intervals,ch)
%splitByChannel  The results.channel axis: every channel of the recording, with the
%   windows that are analysed on it and THE SAMPLES OF THOSE WINDOWS.  A channel with
%   no window of its own still appears - it is part of the recording, and an empty
%   list says "nothing was defined here" where a missing element would say "there is
%   no such channel".
%
%   THE SAMPLES ARE CUT, NOT COPIED ACROSS.  Each window keeps the samples inside its
%   own boundaries, at that channel's own rate, and the channel's whole-recording
%   .data and .time are emptied: a wire recording's samples are its measurement, and
%   what falls outside the analysed windows is discarded like everything else.  The
%   fields are copied one by one ONTO THE PROTOTYPE rather than by keeping the input
%   element, which is what fixes their ORDER: a channel that came back with its
%   fields in another order would be a different struct holding the same numbers.
channels=myographProduct('channels');
tally=struct('before',0,'kept',0,'unit','samples');
if isempty(intervals), nW=0; else, nW=numel(intervals); end
got=zeros(1,nW);
fps=sharedRate(ch);
fn=reshape(fieldnames(channels),1,[]);
for i=1:1:numel(ch)
    for f=1:1:numel(fn)
        if any(strcmp(fn{f},{'intervals','data','time'})), continue; end
        channels(i).(fn{f})=fieldOr(ch(i),fn{f},[]);
    end
    if isempty(channels(i).name), channels(i).name=sprintf('channel%d',i); end
    channels(i).data=[]; channels(i).time=[];    % the whole-recording copy goes
    nm=char(channels(i).name);
    [y,t]=myographChannelSamples(ch(i));

    keep=false(1,nW);
    for k=1:1:nW
        want=intervals(k).channels;
        keep(k)=isempty(want) || any(strcmpi(strtrim(reshape(cellstr(want),1,[])),strtrim(nm)));
    end
    ivs=intervals(keep);
    kk=find(keep);
    cov=false(numel(t),1);
    for j=1:1:numel(ivs)
        idx=t>=min(ivs(j).tStart,ivs(j).tEnd) & t<=max(ivs(j).tStart,ivs(j).tEnd);
        ivs(j).samples=struct('data',y(idx),'time',t(idx));
        ivs(j).frames=sampleSpan(t(idx),fps);
        got(kk(j))=got(kk(j))+nnz(idx);
        cov=cov | idx;
    end
    channels(i).intervals=ivs;
    tally.before=tally.before+numel(t);
    tally.kept=tally.kept+nnz(cov);
end

% A WINDOW WITH NOTHING BEHIND IT IS NAMED, exactly as the pressure myograph names
% one: after a recording has been cut once, a boundary drawn in a discarded stretch
% is a request to widen, and no result can answer it.
for k=1:1:nW
    if got(k)>0, continue; end
    warning('cutMyographIntervals:outsideSpan', ...
        ['Window "%s" (%.3g-%.3g s) holds no samples of this recording, so it has ' ...
         'no data behind it.  Read the recording again to analyse that stretch.'], ...
        intervals(k).name,intervals(k).tStart,intervals(k).tEnd);
end
end

% =====================================================================
function iv=blankWindow(name,t0,t1,chans)
%blankWindow  One element with its definition filled in and nothing measured yet.
%   THE ANALYSES ARE EMPTY BY CONSTRUCTION: they describe a window, so a window that
%   moved has none until those steps run again.
iv=myographProduct('intervals');
iv(1).name=name;
iv(1).tStart=t0;
iv(1).tEnd=t1;
iv(1).frames=[];
iv(1).channels=chans;
iv(1).diameter=[];
iv(1).samples=[];
iv(1).propagation=[];
iv(1).vasomotion=[];
end

% =====================================================================
function fps=sharedRate(ch)
%sharedRate  The one sampling rate every channel was recorded at, or EMPTY when they
%   were not recorded at one.  It is asked of the RATES alone and not of the sample
%   counts, because the counts stop agreeing the moment the windows are cut and the
%   rate does not: the whole point of .frames is to keep naming the original file.
fps=[];
if isempty(ch), return; end
r=arrayfun(@(c) double(fieldOr(c,'fs',NaN)),ch);
if any(~isfinite(r)) || any(r<=0), return; end
if max(r)-min(r)>1e-9*max(max(r),1), return; end
fps=r(1);
end

% =====================================================================
function fr=sampleSpan(t,fps)
%sampleSpan  A wire window's range in the ORIGINAL recording: the sample numbers of
%   its first and last kept sample, on the rate every channel shares.  Empty when the
%   channels share no rate - there is no single index that would be true for all of
%   them - and empty when the window kept nothing.
fr=[];
if isempty(t) || isempty(fps) || ~isfinite(fps) || fps<=0, return; end
fr=[round(t(1)*fps)+1 round(t(end)*fps)+1];
end

% =====================================================================
function meas=measuresOf(results)
%measuresOf  The measures, and their ORDER, which is the order the cores return them
%   in.  results.recording is the ONE authority on it - the diameter step writes it
%   there in the same breath as the windows - so there is nothing to fall back to and
%   nothing that could disagree.
meas=reshape(cellstr(fieldOr(fieldOr(results,'recording',struct()),'measures',{})),1,[]);
end

% =====================================================================
function X=cellCat(parts,get)
%cellCat  One array per window, laid end to end along time.
X=[];
if isempty(parts), return; end
pieces=cell(1,numel(parts));
for k=1:1:numel(parts), pieces{k}=get(parts{k}); end
X=cat(1,pieces{:});
end

% =====================================================================
function nm=nameOf(names,k)
if numel(names)>=k && ~isempty(names{k}), nm=char(names{k});
else,                                     nm=sprintf('interval%d',k); end
end

% =====================================================================
function c=channelsOf(channels,k)
%channelsOf  The channels this window was selected on, always as a 1xN cellstr.
c={};
if numel(channels)>=k && ~isempty(channels{k}), c=reshape(cellstr(channels{k}),1,[]); end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
