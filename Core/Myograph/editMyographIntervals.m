%editMyographIntervals  Draw, name and adjust a myograph recording's time intervals
%
%   [ivT,names,extra] = editMyographIntervals(data,opts) opens ONE window on ONE
%   recording, lets the operator define the time windows the analyses will run in,
%   and hands them back when Done is pressed.  It is a CORE GUI function - it takes
%   data and returns intervals, opens no file and narrates nothing - so the wrapper
%   that calls it is the one that loops the recordings, exactly the way setRegions
%   calls its own ROI editor.  There is deliberately NO file dropdown and no
%   reference recording: everything is assigned for one recording, then the next
%   window opens.
%
%   THE TABLE IS THE INTERVALS.  A window that has been added is DRAWN on the trace
%   as a shaded band and is NOT draggable there: its start, its end, its name and -
%   for a wire myograph - the channel it is analysed on are all set in the table
%   beside the plot.  That is not a restriction, it is what makes the editor usable
%   on a real recording.  A LabChart file is tens of millions of samples; an
%   interactive band sitting on top of that many points cannot be grabbed by its
%   edge, because every mouse event has to re-render the plot underneath it before
%   the next one arrives.  Typing 512.4 into a cell is also simply more accurate than
%   dragging, which is why the table was there from the start.
%
%   WHAT THE WINDOW HOLDS
%     * a full-width, zoomable trace, DECIMATED TO WHAT IS ON SCREEN.  The drawn line
%       is a min/max envelope of opts.drawPoints buckets across the visible span, so
%       a spike survives but fifteen million samples are never handed to the
%       renderer; zooming in re-decimates and the detail comes back.  IT IS THE
%       DRAWING THAT IS DECIMATED AND NOTHING ELSE - the intervals are cut out of the
%       raw samples, and every analysis downstream reads those;
%     * the intervals already defined, drawn as shaded, labelled, NON-INTERACTIVE
%       bands - the selected one brighter than the rest;
%     * ADD - drag the yellow selection box over a stretch, type a name, press Add.
%       The box then jumps to the RIGHT of the band just made, so a run of intervals
%       is added by repeating one gesture rather than by first dragging the box out
%       from under the band it has just become;
%     * DELETE - removes the interval selected in the table.  It is a plain button,
%       not a mode: a band is not clickable, so there is nothing to arm;
%     * an editable start / end / name TABLE, which is how a boundary is set exactly
%       while the trace is zoomed in on it.  RENAME is that table's name cell;
%     * NUDGE - arrow keys move the interval selected in the table (Ctrl = 5 s at a
%       time, Shift resizes instead of moving).  The nudge does NOT round the
%       boundary: a start typed as 12.35 s stays 12.35 s;
%     * and one hard limit: a band CANNOT LEAVE THE ANALYSED SPAN.  The trace is
%       exactly what was analysed - after a time crop, only the cropped stretch - so
%       a window outside it would select frames that do not exist.  Anything that
%       would take a band past either end is clipped to that end;
%     * a plain-language warning under the table when the intervals overlap or leave
%       gaps.  It WARNS, it does not block - overlapping windows and ignored stretches
%       are both legitimate;
%     * DONE.  Closing the window means Done.
%
%   AND DONE ASKS, WHEN THESE WINDOWS WOULD THROW A MEASUREMENT AWAY.  A myograph
%   product keeps no copy of the recording, so the windows ARE the storage: a stretch
%   left outside every one of them is not written back, and only the recording itself
%   can bring it back.  When data.measurement says the trace is the measurement and
%   some of it falls outside the windows as they stand, Done asks once and says how
%   many frames (or samples) that is.  NOTHING TO LOSE IS SILENT - a window chosen
%   before anything was measured, one that still covers everything, or no windows at
%   all, which means the whole recording.  So is a window that is not on screen: an
%   invisible editor is being driven by a script, and a modal dialog on it would wait
%   for an answer that never comes.  There is no setting that turns the question off.
%
%   PMYO MODE (opts.mode='PMYO') adds, when the data carries it, a VIDEO PREVIEW of
%   the selected interval with the detected walls drawn on it - play / stop, a
%   frame-step spinner and a frame slider.  A frame where a wall left the field of
%   view is drawn RED, in the preview and in the trace, because its diameter is an
%   edge-clamped lower bound rather than a measurement.  THE FRAME AND THE WALLS ON
%   IT COME FROM ONE PLACE - data.wallFrame, which getMyographTrace fills with
%   getMyographWallFrame - so the walls drawn are always the walls of the frame
%   shown; the concatenated data.walls are for the trace, where a position in the
%   windows laid end to end is all that is wanted, and a picture needs the frame of
%   the ORIGINAL recording instead.  With no diameter measured yet (the
%   pre-set-intervals step) the trace is whatever cheap profile the caller passes,
%   the wall overlay is absent, and the video is the point of the window.
%
%   WMYO MODE (opts.mode='WMYO') is the same window with the video column replaced.
%   A wire myograph records channels, not pictures, so what stands where the preview
%   stood is:
%     * a CHANNEL LIST, multi-select, which chooses WHAT IS DRAWN.  It is a view
%       control and nothing else - ticking a channel off does not change any
%       interval - because a four-channel recording is read one or two channels at a
%       time and the analysis must not depend on how the operator was looking at it;
%     * an interval's own CHANNEL, which is the table's fourth column: 'all channels'
%       (the default) or one named channel.  That is what makes a window channel
%       specific - two chambers on one LabChart file get different drug windows -
%       and it is what splits the results (see below);
%     * a COMMENTS TABLE - the marks the operator made while recording, with their
%       time and the channel they were placed on.  Each is also drawn as a vertical
%       line on the trace, so a window can be placed against the event that defines
%       it rather than against the shape of the curve, and clicking a row scrolls
%       the trace to that comment;
%     * the LabChart RECORD BOUNDARIES, shaded, because the time base has real gaps
%       between them - the stretches the operator was stopped.
%   The channels are drawn STACKED, each scaled into its own lane exactly as
%   LabChart draws them: they carry different units, and overlaying millimetres of
%   mercury on degrees celsius says nothing about either.  A band still spans the
%   whole height, because it still selects a stretch of TIME.
%
%   Syntax:
%      [ivT,names,extra] = editMyographIntervals(data,opts)
%
%   INPUTS
%     data   what is drawn.  Times are ABSOLUTE seconds from the start of the
%            recording, never re-based (see the myograph spec), so a band can always
%            be located back in the original footage.
%       .time      [T x 1] the trace's time base (PMYO; empty in WMYO)
%       .trace     [T x 1] the curve itself (PMYO; empty in WMYO)
%       .traceName (optional) what the curve is, for the title
%       .traceUnit (optional) the y-axis label
%       .valid     (optional) [T x 1] logical; false = drawn red (a wall off-FOV)
%       .walls     (optional) struct .L .R (each [T x m]) and .rows - the walls of
%                  the whole trace, at the image rows they were detected on
%       .wallFrame (optional) @(t,vr) -> the frame of the recording at time t with
%                  THAT frame's own walls on it (getMyographWallFrame).  It is what
%                  the preview draws when it is there, so the picture and the walls
%                  on it always describe the same frame; vr is this window's own
%                  VideoReader, handed over so none is opened per frame
%       .measurement (optional, default false) true when .trace or .channels IS the
%                  measurement, so a stretch left outside every window is discarded
%                  by keeping them.  It is what Done asks the operator about
%       .video     (optional) the recording to preview; '' = no preview panel
%       .fs        (optional) sampling rate.  It is the AUTHORITY on the sampling -
%                  with disjoint intervals data.time has gaps, so it cannot be read
%                  off diff(time) - and it is what breaks the drawn line over them.
%       .channels  (WMYO) struct array .name .units .fs .data .time - one element
%                  per channel, each on ITS OWN time base and sampling rate
%       .comments  (WMYO, optional) struct array .time .text .channel
%       .blocks    (WMYO, optional) [nBlk x 2] the LabChart record boundaries
%     opts   how it opens.
%       .mode         'PMYO' (default) | 'WMYO'
%       .intervals    [n x 2] start/end seconds to pre-load, as editable windows
%       .names        1 x n cellstr of their names
%       .channels     (WMYO) 1 x n cell - the channel each pre-loaded interval is
%                     analysed on, as a cellstr.  {} means every channel; a cellstr
%                     naming ONE channel means that channel.  A NAME THIS RECORDING
%                     DOES NOT HAVE is dropped rather than shown - a set of windows
%                     carried over to another recording of the protocol must not
%                     claim a chamber that is not there.
%       .maxIntervals (default Inf) the largest number of windows allowed.  1 is the
%                     time crop: one window, so Add greys out once it exists.
%       .drawPoints   (default 2000) HOW MUCH OF THE SIGNAL IS DRAWN - the number of
%                     min/max buckets the visible span is cut into, per channel.  It
%                     is a DRAWING setting and touches nothing else: the raw samples
%                     are what the intervals are cut out of and what every analysis
%                     reads, and they are never decimated.  Lower it on a slow
%                     machine or a very long recording, raise it if the envelope is
%                     too blocky to place a window against.
%       .title        the window title
%       .readyFcn     (optional) @(fig) called once, with the window built and the
%                     windows pre-loaded, BEFORE the blocking wait.  It is how the
%                     editor is driven without gestures (see the API below) - a hook
%                     that edits and then presses Done.
%       .visible      'on' (default) | 'off'
%
%   OUTPUTS
%     ivT    [n x 2] the intervals in seconds, sorted by start time
%     names  1 x n cellstr of their names, in the same order
%     extra  per-interval side information.  PMYO has none, so it is a 0x0 struct
%            array; the wire myograph returns one element per interval whose
%            .channels is {} for a window analysed on every channel and {'name'} for
%            one that belongs to a single channel.
%
%   DRIVING IT WITHOUT GESTURES.  getappdata(fig,'intervalsAPI') is a struct of
%   handles onto the same callbacks the buttons use - add / rename / editCell /
%   remove / select / press / preview / done, plus setShownChannels /
%   setIntervalChannel / channels / clickComment for the wire myograph - so a test
%   adds, renames, nudges and deletes intervals exactly as a person would and reads
%   back what comes out.  Same arrangement as guiExplore and guiExport.  It is
%   reached through opts.readyFcn, which is the editor's own contract: nothing of
%   this rides in a wrapper's s.
%
%   THE BLOCKING WAIT IS waitfor ON A PRIVATE SENTINEL, never uiwait: the selection
%   box runs its OWN uiwait/uiresume on this figure while it is being dragged, and
%   that inner uiresume also releases an outer uiwait - the editor would return the
%   moment the box moved.  setRegions>editRegions documents the identical trap.
%   THE VIDEO AND THE PLAYBACK TIMER ARE RELEASED ON THE WAY OUT, whatever happens,
%   or the window leaks a timer that keeps ticking at a figure that is gone.
%   NOTHING OUTSIDE THIS WINDOW IS TOUCHED: writing to a component of the caller's
%   parked window would pull it back to the front, over this one.
%
%   EXAMPLE
%     data = getMyographTrace(results,s);
%     [ivT,names] = editMyographIntervals(data, ...
%         struct('title','Intervals - Mouse1.avi','intervals',ivT0,'names',{names0}));
%
%   DEPENDS ON
%     MATLAB's Image Processing Toolbox (drawrectangle, for the selection box only)
%     and base VideoReader.
%
% See also: setMyographIntervals, setMyographPresetIntervals, setMyographCrop,
%           getMyographTrace, readLabChart, myographProduct, setRegions
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function [ivT,names,extra] = editMyographIntervals(data,opts)

if nargin<2, opts=struct(); end
opts=withOptsDefaults(opts);
data=withDataDefaults(data,opts);

fig=buildWindow(data,opts);
% The window, the video and the playback timer are released on the way out however
% this function ends - Done, a closed window, or an error inside a callback.
guard=onCleanup(@() releaseResources(fig));

preloadIntervals(fig,opts);
drawSignal(fig);                 % which draws the bands over it, in that order
refreshTable(fig);
parkSelectionBox(fig);           % out from under the windows that are already there

if ~isempty(opts.readyFcn)
    opts.readyFcn(fig);          % headless: the caller drives the API and presses Done
end

waitfor(fig,'UserData','done');  % NOT uiwait - see the header

[ivT,names]=collect(fig);
extra=collectExtra(fig);
end

%% ===================== defaults ===================================== %%
function d=withDataDefaults(data,opts)
%withDataDefaults  The drawn data, in the one shape the rest of the file assumes.
if ~isstruct(data), error('editMyographIntervals:badData','data must be a struct.'); end
d=data;
d.time=double(fieldOr(d,'time',[])); d.time=d.time(:);
d.trace=double(fieldOr(d,'trace',[])); d.trace=d.trace(:);
if numel(d.trace)~=numel(d.time)
    error('editMyographIntervals:badData','data.time and data.trace must be the same length.');
end
for f={'traceName','traceUnit','video'}
    if ~isfield(d,f{1}) || isempty(d.(f{1})), d.(f{1})=''; else, d.(f{1})=char(d.(f{1})); end
end
for f={'valid','walls','wallFrame'}
    if ~isfield(d,f{1}), d.(f{1})=[]; end
end
% A caller that says nothing is taken to be showing the recording rather than the
% measurement: the question Done asks is one nobody should be asked by accident.
if ~isfield(d,'measurement') || isempty(d.measurement), d.measurement=false; end
d.measurement=logical(d.measurement(1));
d=withChannelDefaults(d,opts);
if ~isempty(d.walls) && ~isfield(d.walls,'rows'), d.walls.rows=1:size(d.walls.L,2); end
if ~isempty(d.valid), d.valid=logical(d.valid(:)); end
if numel(d.valid)~=numel(d.time), d.valid=[]; end
% fs is the authority on the sampling: with disjoint intervals the time base has
% gaps, so the median step is the only sane fallback and diff() is not.
if ~isfield(d,'fs') || isempty(d.fs) || ~isfinite(d.fs)
    if numel(d.time)>1, d.fs=1/median(diff(d.time)); else, d.fs=1; end
end
d.fs=double(d.fs);
end

function d=withChannelDefaults(d,opts)
%withChannelDefaults  The wire myograph's own half of the data: the channels that
%   are drawn and analysed, the comments that are marked on them, and the record
%   boundaries that are shaded.  Only the channels are required, because they are
%   the recording; a file with no comments and one record is ordinary.
if ~isfield(d,'channels'), d.channels=struct('name',{},'units',{},'fs',{},'data',{},'time',{}); end
if ~isfield(d,'comments') || isempty(d.comments)
    d.comments=struct('time',{},'text',{},'channel',{});
end
if ~isfield(d,'blocks') || isempty(d.blocks), d.blocks=zeros(0,2); end
d.blocks=double(d.blocks);
if size(d.blocks,2)~=2, d.blocks=zeros(0,2); end
if ~isWmyo(opts), return; end
if isempty(d.channels)
    error('editMyographIntervals:noChannels', ...
        'A wire myograph recording has to carry at least one channel to show.');
end
for i=1:1:numel(d.channels)
    c=d.channels(i);
    d.channels(i).name=char(fieldOr(c,'name',sprintf('channel%d',i)));
    d.channels(i).units=char(fieldOr(c,'units',''));
    d.channels(i).data=double(fieldOr(c,'data',[])); d.channels(i).data=d.channels(i).data(:);
    d.channels(i).time=double(fieldOr(c,'time',[])); d.channels(i).time=d.channels(i).time(:);
    fs=double(fieldOr(c,'fs',[]));
    if isempty(fs) || ~isfinite(fs) || fs<=0
        tt=d.channels(i).time;
        if numel(tt)>1, fs=1/median(diff(tt)); else, fs=1; end
    end
    d.channels(i).fs=fs;
end
end

function o=withOptsDefaults(opts)
%withOptsDefaults  How the window opens, with everything optional defaulted.
o=opts;
def=struct('mode','PMYO','intervals',[],'names',{{}},'channels',{{}},'maxIntervals',Inf, ...
    'drawPoints',2000,'title','Define intervals','readyFcn',[],'visible','on');
fn=fieldnames(def);
for i=1:1:numel(fn)
    if ~isfield(o,fn{i}) || (isempty(o.(fn{i})) && ~isnumeric(def.(fn{i}))), o.(fn{i})=def.(fn{i}); end
end
if isempty(o.intervals), o.intervals=zeros(0,2); end
o.intervals=double(o.intervals);
if size(o.intervals,2)~=2, o.intervals=zeros(0,2); end
if isempty(o.names), o.names={}; else, o.names=reshape(cellstr(o.names),1,[]); end
if isempty(o.channels), o.channels={}; else, o.channels=reshape(o.channels,1,[]); end
% Stated as what a ceiling IS, not as the ways of being wrong: NaN passes every
% negative test ever written for it (NaN<1 is false, so is NaN>1), and a NaN ceiling
% would compare false against every count and silently mean "unlimited".
okCeiling=isnumeric(o.maxIntervals) && isscalar(o.maxIntervals) && isreal(o.maxIntervals) && ...
    o.maxIntervals>=1 && (isinf(o.maxIntervals) || mod(o.maxIntervals,1)==0);
if ~okCeiling, o.maxIntervals=Inf; end
% Same shape of test for the drawing budget, and for the same reason: it is stated
% as what a valid budget IS.  An infinite one would ask for every sample of a
% fifteen-million-sample channel, which is the thing this setting exists to prevent.
okDraw=isnumeric(o.drawPoints) && isscalar(o.drawPoints) && isreal(o.drawPoints) && ...
    isfinite(o.drawPoints) && o.drawPoints>=50;
if ~okDraw, o.drawPoints=2000; end
o.drawPoints=round(double(o.drawPoints));
o.mode=upper(char(o.mode)); o.title=char(o.title); o.visible=char(o.visible);
if ~any(strcmp(o.mode,{'PMYO','WMYO'}))
    error('editMyographIntervals:unknownMode', ...
        'This editor opens on a pressure myograph or a wire myograph, not on "%s".',o.mode);
end
end

function tf=isWmyo(o), tf=isstruct(o) && isfield(o,'mode') && strcmpi(o.mode,'WMYO'); end

%% ===================== the window =================================== %%
function fig=buildWindow(data,opts)
%buildWindow  The trace on the left, everything that is decided on the right.  The
%   preview panel exists only when there is a recording to preview.

fig=uifigure('Name',opts.title,'Visible',opts.visible, ...
    'Position',figurePosition(),'Color','w');
fig.UserData='editing';                                  % the private wait sentinel
fig.CloseRequestFcn=@(~,~) requestDone(fig);             % closing the window == Done
fig.WindowKeyPressFcn=@(~,e) onKey(fig,e);

wire=isWmyo(opts);
hasVideo=~wire && ~isempty(data.video);

gl=uigridlayout(fig,[1 2],'ColumnWidth',{'2.2x','1x'}, ...
    'Padding',[8 8 8 8],'ColumnSpacing',10,'BackgroundColor','w');

% ---- left: the trace, over the whole height ----------------------------------
lp=uigridlayout(gl,[1 1],'Padding',[0 0 0 0],'BackgroundColor','w');
c.axTrace=uiaxes(lp);

% ---- right: name it, list it, preview it, finish ------------------------------
rows={'fit','fit','1x','fit'};
if hasVideo, rows=[rows,{'1.5x'}]; end
if wire
    % The table is where an interval is now edited - four columns of it - so it gets
    % the room, and the two panels beside it share what is left.
    rows{3}='1.6x'; rows=[rows,{'1x','1x'}];
end
rows=[rows,{'fit'}];
rp=uigridlayout(gl,[numel(rows) 1],'RowHeight',rows,'RowSpacing',7, ...
    'Padding',[0 0 0 0],'BackgroundColor','w');

uilabel(rp,'Text',hintText(opts),'WordWrap','on','FontSize',11);

nr=uigridlayout(rp,[1 3],'ColumnWidth',{'1x','fit','fit'},'Padding',[0 0 0 0], ...
    'ColumnSpacing',5,'BackgroundColor','w');
c.name=uieditfield(nr,'text','Value','','Placeholder','name of the next interval');
c.addBtn=uibutton(nr,'Text','Add','ButtonPushedFcn',@(~,~) onAddPressed(fig));
c.delBtn=uibutton(nr,'Text','Delete','Enable','off', ...
    'ButtonPushedFcn',@(~,~) onDeletePressed(fig), ...
    'Tooltip','remove the interval selected in the table');

c.table=uitable(rp,'ColumnName',tableColumns(data,opts), ...
    'ColumnEditable',true(1,3+wire),'ColumnFormat',tableFormats(data,opts), ...
    'Data',cell(0,3+wire),'CellEditCallback',@(~,e) onTableEdit(fig,e), ...
    'CellSelectionCallback',@(~,e) onRowSelected(fig,e));

c.warn=uilabel(rp,'Text','','WordWrap','on','FontSize',11,'FontColor',[0.75 0.35 0]);

if hasVideo, c=buildPreview(fig,rp,c); end
if wire, c=buildChannelList(fig,rp,c,data); c=buildCommentTable(fig,rp,c,data); end

c.doneBtn=uibutton(rp,'Text','Done','FontWeight','bold', ...
    'BackgroundColor',[0.82 0.92 0.82],'ButtonPushedFcn',@(~,~) requestDone(fig));

app=struct('data',data,'opts',opts,'c',c,'iv',emptyIntervalModel(), ...
    'sel',0,'selRoi',[],'yl',[0 1],'shown',{{}},'sigLayer',gobjects(0), ...
    'bandLayer',gobjects(0),'cache',[],'busy',false, ...
    'preview',struct('vr',[],'times',[],'curIdx',1,'timer',[]));
if wire, app.shown=selectedChannelNames(c); end
setApp(fig,app);

fitSpan(fig);
% RE-DECIMATE ON ZOOM.  The drawn line is an envelope of what is on screen, so it
% has to be rebuilt when the screen changes; without this the operator zooms into a
% flat bucket instead of into the samples.  Guarded, because a release that does not
% let the axes be listened to should cost the detail on deep zoom, not the window.
try
    c.xlimL=addlistener(c.axTrace,'XLim','PostSet',@(~,~) onXLimChanged(fig));
    app=getApp(fig); app.c=c; setApp(fig,app);
catch
end
setappdata(fig,'intervalsAPI',buildAPI(fig));
end

function m=emptyIntervalModel()
%emptyIntervalModel  THE INTERVALS THEMSELVES - a plain struct array, and the only
%   truth in this window.  The bands are drawn FROM it and the table is a view OF
%   it; when the band was the truth, every read had to go and ask a graphics object
%   where it had been dragged to, and a band that is no longer draggable has nothing
%   to be asked.  .channel is '' for a window analysed on every channel.
m=struct('t0',{},'t1',{},'name',{},'channel',{});
end

function cn=tableColumns(~,opts)
if isWmyo(opts), cn={'start (s)','end (s)','name','channel'};
else,            cn={'start (s)','end (s)','name'};
end
end

function cf=tableFormats(data,opts)
%tableFormats  The channel cell is a DROP-DOWN of this recording's own channels plus
%   'all channels', so a window can only ever name a channel that exists.
cf={'bank','bank','char'};
if ~isWmyo(opts), return; end
cf{4}=[{allChannelsWord()},reshape({data.channels.name},1,[])];
end

function w=allChannelsWord(), w='all channels'; end

function c=buildChannelList(fig,parent,c,data)
%buildChannelList  WHICH CHANNELS ARE DRAWN, and nothing else.  It used to set what
%   the selected interval was analysed on as well; that now lives in the table's own
%   channel column, because an analysis that changed with how the operator happened
%   to be looking at the recording was a trap - tick a channel off to read another
%   one and the window you had already defined quietly changed meaning.
c.chPanel=uipanel(parent,'Title','Channels shown','FontWeight','bold','BackgroundColor','w');
g=uigridlayout(c.chPanel,[1 1],'Padding',[6 6 6 6],'BackgroundColor','w');
items=reshape({data.channels.name},1,[]);
c.chList=uilistbox(g,'Items',items,'Multiselect','on','Value',items, ...
    'ValueChangedFcn',@(~,~) onShownChannelsChanged(fig), ...
    'Tooltip','which channels are drawn - an interval''s own channel is set in the table');
end

function c=buildCommentTable(fig,parent,c,data)
%buildCommentTable  The marks the operator made while recording.  They are the
%   reason a window sits where it does - 'drug in' is an event, not a bend in the
%   curve - so they are both a list to click and a line on the trace.
c.comPanel=uipanel(parent,'Title','Comments','FontWeight','bold','BackgroundColor','w');
g=uigridlayout(c.comPanel,[1 1],'Padding',[6 6 6 6],'BackgroundColor','w');
c.comTable=uitable(g,'ColumnName',{'time (s)','comment','channel'}, ...
    'ColumnEditable',[false false false],'ColumnFormat',{'bank','char','char'}, ...
    'Data',commentRows(data.comments), ...
    'CellSelectionCallback',@(~,e) onCommentSelected(fig,e));
end

function D=commentRows(com)
%commentRows  The comments as the table shows them, earliest first.
D=cell(numel(com),3);
for k=1:1:numel(com)
    ch=char(fieldOr(com(k),'channel',''));
    if isempty(ch), ch=allChannelsWord(); end
    D(k,:)={double(com(k).time),char(string(fieldOr(com(k),'text',''))),ch};
end
end

function c=buildPreview(fig,parent,c)
%buildPreview  The recording itself, at the selected interval, with the walls the
%   detector found drawn over it - the one panel that says whether a window holds
%   what the operator thinks it holds.
c.prevPanel=uipanel(parent,'Title','Recording','FontWeight','bold','BackgroundColor','w');
g=uigridlayout(c.prevPanel,[4 1],'RowHeight',{'1x','fit','fit','fit'},'RowSpacing',5, ...
    'Padding',[6 6 6 6],'BackgroundColor','w');
c.axPrev=uiaxes(g); c.axPrev.XTick=[]; c.axPrev.YTick=[];
title(c.axPrev,'select an interval to preview it');
c.slider=uislider(g,'Limits',[1 2],'Value',1, ...
    'ValueChangedFcn',@(o,~) onSlide(fig,o.Value));
br=uigridlayout(g,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',5,'BackgroundColor','w');
uibutton(br,'Text','Play','ButtonPushedFcn',@(~,~) previewPlay(fig));
uibutton(br,'Text','Stop','ButtonPushedFcn',@(~,~) previewStop(fig));
sr=uigridlayout(g,[1 2],'ColumnWidth',{'fit','1x'},'Padding',[0 0 0 0], ...
    'ColumnSpacing',5,'BackgroundColor','w');
uilabel(sr,'Text','Frames per step');
c.step=uispinner(sr,'Limits',[1 200],'Value',2,'Step',1,'RoundFractionalValues','on', ...
    'Tooltip','how many frames playback advances at a time - higher plays faster and skips more');
end

function pos=figurePosition()
%figurePosition  A wide window, centred, that fits on a laptop screen.
sz=get(0,'ScreenSize');
w=min(1500,round(sz(3)*0.9)); h=min(820,round(sz(4)*0.82));
pos=[max(1,round((sz(3)-w)/2)) max(41,round((sz(4)-h)/2)) w h];
end

function t=hintText(opts)
%hintText  What to do, in the order it is done.  One crop window reads differently
%   from a set of intervals, so it says so.
if opts.maxIntervals==1
    t=['Drag the yellow box over the stretch to analyse, then press Add.  Set the ' ...
       'exact start and end in the table.  Only one window can be set here; delete ' ...
       'it to draw another.  Press Done when ready.'];
else
    t=['Drag the yellow box over a stretch, name it and press Add.  An interval is ' ...
       'then adjusted in the TABLE - its start, its end and its name - or nudged ' ...
       'with the arrow keys (Ctrl for 5 s, Shift to resize).  Delete removes the ' ...
       'selected one.  Press Done when ready.'];
end
if isWmyo(opts)
    t=[t '  Set an interval''s channel in the table to analyse it on that channel ' ...
       'alone.  Tick the channel list to choose what is drawn.  Click a comment to ' ...
       'jump to it.'];
end
end

%% ===================== the drawn signal ============================= %%
function fitSpan(fig)
%fitSpan  The axes over the whole recording, once, before anything is drawn.  The
%   x limits are what the decimation is computed against, so they exist first.
app=getApp(fig); ax=app.c.axTrace;
sp=spanOf(app.data,app.opts);
if isWmyo(app.opts), app.yl=[0 1]; else, app.yl=traceYLimits(app.data); end
ylim(ax,app.yl);
setApp(fig,app);
if ~isempty(sp) && sp(2)>sp(1), xlim(ax,sp); end
end

function yl=traceYLimits(d)
%traceYLimits  The pressure myograph's y range, computed from the data ONCE and then
%   frozen - the bands are drawn full height against it, so it may not drift when a
%   band is added.
yl=[0 1];
y=d.trace(isfinite(d.trace));
if isempty(y), return; end
lo=min(y); hi=max(y);
if hi<=lo, hi=lo+1; end
pad=0.05*(hi-lo);
yl=[lo-pad hi+pad];
end

function onXLimChanged(fig)
%onXLimChanged  Zoom or pan: re-decimate to the new view.  The latch is not
%   politeness, it is required - drawSignal sets the very property this listens to.
if ~isvalid(fig), return; end
app=getApp(fig);
if ~isstruct(app) || app.busy, return; end
drawSignal(fig);
end

function drawSignal(fig)
%drawSignal  THE ONE EXPENSIVE LAYER, and the only one decimated.  It is redrawn
%   when the view changes or a channel is ticked; the interval bands are a separate
%   layer on top of it and are not touched here, so adding an interval never pays
%   for redrawing fifteen million samples.
app=getApp(fig);
if app.busy, return; end
app.busy=true; setApp(fig,app);
restore=onCleanup(@() clearBusy(fig));
delete(app.sigLayer(isvalid(app.sigLayer)));
if isWmyo(app.opts), h=drawChannels(fig); else, h=drawTrace(fig); end
app=getApp(fig); app.sigLayer=h; setApp(fig,app);
% THE BANDS GO BACK ON TOP.  New graphics are added at the front of the axes, so a
% re-drawn signal would otherwise sit over the windows it is supposed to sit under -
% and the record shading of a wire recording would wash them out entirely.  Redrawing
% them here rather than at each call site is what keeps that rule in one place.
drawBands(fig);
end

function clearBusy(fig)
if ~isvalid(fig), return; end
app=getApp(fig);
if isstruct(app), app.busy=false; setApp(fig,app); end
end

function h=drawTrace(fig)
%drawTrace  The pressure myograph's curve, BROKEN WHERE THE RECORDING IS.  With
%   disjoint pre-set intervals the time base jumps, and a straight line drawn across
%   an hour that was never measured is the one thing this plot must not show.
app=getApp(fig); d=app.data; ax=app.c.axTrace;
cache=traceCache(fig);
h=gobjects(0);
xl=xlim(ax); nb=viewBuckets(app.opts);
hold(ax,'on');
[x,y]=decimateRuns(d.time,d.trace,cache.runs,xl,nb);
if ~isempty(x)
    h(end+1)=plot(ax,x,y,'-','Color',[0 0.3 0.8],'LineWidth',1);
end
if ~isempty(d.valid) && any(~d.valid)
    yb=d.trace; yb(d.valid)=NaN;
    [xr,yr]=decimateRuns(d.time,yb,cache.runs,xl,nb);
    if ~isempty(xr)
        h(end+1)=plot(ax,xr,yr,'-','Color',[0.85 0.1 0.1],'LineWidth',1.6);
    end
end
hold(ax,'off'); grid(ax,'on'); box(ax,'on');
xlabel(ax,'time (s)');
ylabel(ax,traceUnitOf(d));
title(ax,traceTitle(d));
ylim(ax,app.yl);                               % frozen: the bands are drawn full height
set(h,'PickableParts','none','HitTest','off');
end

function h=drawChannels(fig)
%drawChannels  The wire myograph's trace panel: the selected channels STACKED, one
%   lane each, exactly as LabChart draws them.  Every lane is scaled to its own
%   range - over the WHOLE channel, not over what is on screen, so a lane does not
%   jump while the operator pans - and the numbers are on the axis label rather than
%   on the y ticks, because the question this panel answers is WHEN, not how much.
%
%   The y limits are FIXED at [0 1] however many lanes there are, so a band never
%   has to be re-fitted after a channel is ticked.
app=getApp(fig); d=app.data; ax=app.c.axTrace;
h=gobjects(0);
sel=shownChannelIdx(app);
xl=xlim(ax); nb=viewBuckets(app.opts);
hold(ax,'on');
h=[h shadeBlocks(ax,d.blocks)];
n=numel(sel);
ticks=[]; labels={};
for k=1:1:n
    c=d.channels(sel(k));
    cache=channelCache(fig,sel(k));
    lo=(n-k)/n; hh=1/n;                              % first selected channel on top
    [x,y]=decimateRuns(c.time,c.data,cache.runs,xl,nb);
    if ~isempty(x)
        y=laneScale(y,cache.mn,cache.mx,lo+0.12*hh,lo+0.88*hh);
        h(end+1)=plot(ax,x,y,'-','Color',[0 0.3 0.8],'LineWidth',1); %#ok<AGROW>
    end
    ticks(end+1)=lo+0.5*hh; %#ok<AGROW>
    labels{end+1}=laneLabel(c); %#ok<AGROW>
end
h=[h markComments(ax,d.comments,{d.channels(sel).name},xl)];
hold(ax,'off'); grid(ax,'on'); box(ax,'on');
ylim(ax,[0 1]);
% The lanes are laid out top-down, so the ticks come out descending; an axis wants
% them ascending, and the labels have to follow them.
[ticks,ord]=sort(ticks);
if isempty(ticks), ax.YTick=[]; else, ax.YTick=ticks; ax.YTickLabel=labels(ord); end
xlabel(ax,'time (s)');
ylabel(ax,'');
title(ax,channelsTitle(n,numel(d.channels)));
set(h,'PickableParts','none','HitTest','off');
end

%% ===================== decimation =================================== %%
% THE WHOLE REASON THIS WINDOW IS USABLE.  A LabChart recording is 15 million samples
% per channel; handing that to the renderer costs about a second per redraw, and a
% second per redraw is what made a band impossible to grab by its edge - every mouse
% event had to wait for the one before it.  What is drawn instead is a MIN/MAX
% ENVELOPE of the visible span: the samples are bucketed, and each bucket contributes
% its smallest and its largest value AT THEIR OWN TIMES.  A spike therefore survives
% decimation - which a plain 1:N:end subsample does not - and zooming in re-buckets,
% so the real samples come back as soon as there are few enough of them to draw.

function n=viewBuckets(opts)
%viewBuckets  HOW MANY BUCKETS the visible span is cut into - opts.drawPoints, the
%   one knob on the decimation.  It is a CEILING ON WHAT IS DRAWN and nothing else:
%   the recording, the intervals and everything the analyses read come from the raw
%   samples, which this never touches.  Lower it on a slow machine or a very long
%   recording; raise it if the envelope looks too blocky to place a window against.
n=2000;
if isstruct(opts) && isfield(opts,'drawPoints') && ~isempty(opts.drawPoints)
    n=double(opts.drawPoints);
end
n=max(50,round(n));
end

function [x,y]=decimateRuns(t,sig,runs,xl,nBucket)
%decimateRuns  The visible part of one signal, decimated run by run and joined with
%   NaN, so ONE line object carries the whole channel and the gaps between LabChart
%   records stay gaps.  The budget is shared between the runs THAT CAN BE SEEN, not
%   between all of them: with twenty disjoint pre-set intervals, nineteen of them are
%   off screen and must not each be charged a twentieth of the detail.
x=[]; y=[];
if isempty(t), return; end
vis=zeros(0,2);
for i=1:1:numel(runs)
    k=runs{i};
    [a,b]=visibleRange(t,k(1),k(end),xl);
    if isempty(a), continue; end
    vis(end+1,:)=[a b]; %#ok<AGROW>
end
if isempty(vis), return; end
per=max(32,floor(nBucket/size(vis,1)));
for i=1:1:size(vis,1)
    [xr,yr]=envelopeRange(t,sig,vis(i,1),vis(i,2),per);
    if isempty(xr), continue; end
    x=[x;xr;NaN]; y=[y;yr;NaN]; %#ok<AGROW>
end
end

function [a,b]=visibleRange(t,i0,i1,xl)
%visibleRange  The index range of one run that the current view can see, found by
%   BINARY SEARCH over t itself.  A linear find() over fifteen million samples on
%   every pan frame is the thing this window exists to stop doing, and discretize
%   would be worse still - it takes edges by value, so it would COPY the run.
a=[]; b=[];
lo=xl(1); hi=xl(2);
w=hi-lo; lo=lo-0.02*w; hi=hi+0.02*w;            % a margin, so a pan has something to show
if t(i1)<lo || t(i0)>hi, return; end
a=lowerBound(t,i0,i1,lo);
b=min(lowerBound(t,i0,i1,hi)+1,i1);
if b<a, a=[]; b=[]; end
end

function j=lowerBound(t,i0,i1,v)
%lowerBound  The last index of t(i0:i1) whose value is at or below v, clamped into
%   the run.  Plain binary search on the array in place: no copy, log2(n) reads.
if v<=t(i0), j=i0; return; end
if v>=t(i1), j=i1; return; end
lo=i0; hi=i1;
while hi-lo>1
    mid=floor((lo+hi)/2);
    if t(mid)<=v, lo=mid; else, hi=mid; end
end
j=lo;
end

function [x,y]=envelopeRange(t,sig,a,b,nBucket)
%envelopeRange  One run's visible samples as at most 2*nBucket points: each bucket
%   contributes its smallest and its largest value, AT THEIR OWN TIMES, in time
%   order.  Under the budget the real samples are drawn, so a zoomed-in view is
%   exact rather than an envelope of an envelope.
n=b-a+1;
if n<=2*nBucket
    x=t(a:b); y=sig(a:b); return
end
step=ceil(n/nBucket);
m=floor(n/step)*step;
B=reshape(sig(a:a+m-1),step,[]);
[~,ilo]=min(B,[],1); [~,ihi]=max(B,[],1);
base=a+(0:size(B,2)-1)*step;
p=base+ilo-1; q=base+ihi-1;
first=min(p,q); second=max(p,q);                 % keep the drawn line monotone in x
x=reshape([t(first).';t(second).'],[],1);
y=reshape([sig(first).';sig(second).'],[],1);
if a+m<=b
    x=[x;t(a+m:b)]; y=[y;sig(a+m:b)];            % the tail, drawn as it is
end
end

function cache=traceCache(fig)
%traceCache  The pressure myograph's contiguous runs, found ONCE.  diff() over the
%   whole time base is cheap compared with plotting it, and free compared with doing
%   it again on every pan frame.
app=getApp(fig);
if ~isempty(app.cache), cache=app.cache; return; end
cache=struct('runs',{contiguousRuns(app.data.time,app.data.fs)});
app.cache=cache; setApp(fig,app);
end

function cache=channelCache(fig,i)
%channelCache  One channel's runs and its full-recording range, found ONCE.  The
%   range has to be over the whole channel and not over the view, or the lane would
%   rescale under the operator every time they panned.
app=getApp(fig);
if isempty(app.cache), app.cache=cell(1,numel(app.data.channels)); end
if numel(app.cache)>=i && ~isempty(app.cache{i}), cache=app.cache{i}; return; end
c=app.data.channels(i);
mn=min(c.data); mx=max(c.data);
if ~isfinite(mn) || ~isfinite(mx) || mx<=mn, mn=0; mx=0; end
cache=struct('runs',{contiguousRuns(c.time,c.fs)},'mn',mn,'mx',mx);
app.cache{i}=cache; setApp(fig,app);
end

function y=laneScale(x,mn,mx,lo,hi)
%laneScale  One channel into one lane, against the channel's OWN full range.  A flat
%   channel sits in the middle of its lane rather than collapsing onto its floor,
%   which is what a zero range would otherwise draw.
x=double(x(:));
if mx<=mn
    y=repmat((lo+hi)/2,size(x));
else
    y=lo+(hi-lo)*(x-mn)/(mx-mn);
end
end

function h=shadeBlocks(ax,blocks)
%shadeBlocks  The LabChart records, alternately shaded, as ONE patch with one face
%   per record.  Between two records the operator was stopped, and the time base
%   really does jump there; the shading is what stops that gap being read as a
%   stretch of quiet signal.
h=gobjects(0);
r=2:2:size(blocks,1);
if isempty(r), return; end
V=zeros(4*numel(r),2); F=zeros(numel(r),4);
for i=1:1:numel(r)
    b=blocks(r(i),:);
    V(4*i-3:4*i,:)=[b(1) 0; b(2) 0; b(2) 1; b(1) 1];
    F(i,:)=4*i-3:4*i;
end
h=patch(ax,'Faces',F,'Vertices',V,'FaceColor',[0.92 0.92 0.95], ...
    'EdgeColor','none','FaceAlpha',0.8,'HandleVisibility','off');
end

function h=markComments(ax,com,shown,xl)
%markComments  Each comment where it happened, as ONE line object holding every
%   mark and ONE text call for the labels that fall inside the view.  It used to be
%   an xline apiece: fifty constant-line objects, each re-laying out its own label
%   whenever the axes changed, on a plot that changes on every pan.
%
%   A comment placed on a channel is drawn only while that channel is on screen; one
%   placed on the whole file is always drawn, because it always applies.
h=gobjects(0);
if isempty(com), return; end
t=[]; txt={};
for k=1:1:numel(com)
    ch=char(fieldOr(com(k),'channel',''));
    if ~isempty(ch) && ~any(strcmpi(strtrim(ch),strtrim(shown))), continue; end
    t(end+1)=double(com(k).time); %#ok<AGROW>
    txt{end+1}=shorten(fieldOr(com(k),'text','')); %#ok<AGROW>
end
if isempty(t), return; end
x=reshape([t;t;nan(1,numel(t))],[],1);
y=repmat([0;1;NaN],numel(t),1);
h(end+1)=plot(ax,x,y,'--','Color',[0.85 0.35 0],'LineWidth',0.8);
in=t>=xl(1) & t<=xl(2);
% The labels are the expensive half, so only the ones that can be read are drawn -
% and only while there are few enough of them to tell apart.  text() with a vector
% position returns ONE HANDLE PER LABEL, so they are appended as an array.
if any(in) && nnz(in)<=30
    ht=text(ax,t(in),repmat(0.985,1,nnz(in)),txt(in), ...
        'FontSize',9,'Color',[0.85 0.35 0],'Rotation',0, ...
        'HorizontalAlignment','left','VerticalAlignment','top','Clipping','on');
    h=[h reshape(ht,1,[])];
end
end

function s=laneLabel(c)
u=char(c.units);
if isempty(u), s=char(c.name); else, s=sprintf('%s (%s)',char(c.name),u); end
end

function s=shorten(txt)
s=strtrim(char(string(txt)));
if numel(s)>24, s=[s(1:21) '...']; end
end

function s=channelsTitle(nSel,nAll)
if nSel==0
    s='No channel selected - tick one to draw it';
elseif nSel==nAll
    s='Recording';
else
    s=sprintf('Recording (%d of %d channels shown)',nSel,nAll);
end
end

function idx=shownChannelIdx(app)
%shownChannelIdx  Which channels the list currently holds, as indices into the data,
%   in the data's own order - so the lanes do not reshuffle when the operator ticks
%   one back on.
idx=[];
if ~isWmyo(app.opts) || isempty(app.data.channels), return; end
names={app.data.channels.name};
want=app.shown;
if isempty(want), return; end
idx=find(ismember(lower(strtrim(names)),lower(strtrim(reshape(cellstr(want),1,[])))));
end

function names=selectedChannelNames(c)
%selectedChannelNames  The list's own Value, always as a 1xN cellstr - a uilistbox
%   hands back a bare char when exactly one item is selected.
names={};
if ~isfield(c,'chList') || ~isvalid(c.chList), return; end
v=c.chList.Value;
if isempty(v), return; end
if ischar(v) || isstring(v), names={char(v)}; else, names=reshape(cellstr(v),1,[]); end
end

function sp=spanOf(d,opts)
%spanOf  THE STRETCH A BAND MAY LIVE IN.  For the pressure myograph it is the
%   analysed time base; for the wire myograph it is the union over the channels,
%   which may differ in length because they may differ in sampling rate, so a band
%   is not clipped to whichever channel happens to be shown.
sp=[];
if isWmyo(opts)
    lo=inf; hi=-inf;
    for i=1:1:numel(d.channels)
        t=d.channels(i).time;
        if isempty(t), continue; end
        lo=min(lo,t(1)); hi=max(hi,t(end));
    end
    if isfinite(lo) && isfinite(hi), sp=[lo hi]; end
    return
end
if ~isempty(d.time), sp=[d.time(1) d.time(end)]; end
end

function runs=contiguousRuns(t,fs)
%contiguousRuns  The index blocks the time base is continuous over.  fs is the
%   authority: a gap is a step of more than one and a half samples, which no real
%   recording produces and a concatenation of disjoint intervals always does.
runs={};
n=numel(t);
if n==0, return; end
if n==1, runs={1}; return; end
brk=find(diff(t(:))>1.5/max(fs,eps));
starts=[1;brk+1]; stops=[brk;n];
for i=1:1:numel(starts), runs{end+1}=(starts(i):stops(i))'; end %#ok<AGROW>
end

function s=traceTitle(d)
if isempty(d.traceName), s='Recording'; else, s=char(d.traceName); end
if ~isempty(d.valid) && any(~d.valid)
    s=[s ' (red = a wall left the field of view)'];
end
end
function s=traceUnitOf(d)
if isempty(d.traceUnit), s='signal'; else, s=char(d.traceUnit); end
end

%% ===================== the interval bands =========================== %%
function drawBands(fig)
%drawBands  The intervals as SHADED, LABELLED, NON-INTERACTIVE bands.  They are a
%   layer of their own, drawn from the model and never read back from - so adding,
%   nudging or deleting one costs a handful of patches and not a redraw of the
%   signal underneath.  PickableParts is off on every one of them: a band that
%   swallowed a mouse press would take it away from the zoom and from the selection
%   box, which are the two things the plot is still for.
app=getApp(fig);
delete(app.bandLayer(isvalid(app.bandLayer)));
h=gobjects(0);
ax=app.c.axTrace; yl=app.yl;
hold(ax,'on');
for k=1:1:numel(app.iv)
    isSel=(k==app.sel);
    if isSel, fc=[0.98 0.72 0.10]; fa=0.30; else, fc=[0.20 0.60 1.00]; fa=0.13; end
    x=[app.iv(k).t0 app.iv(k).t1 app.iv(k).t1 app.iv(k).t0];
    h(end+1)=patch(ax,'XData',x,'YData',yl([1 1 2 2]), ...
        'FaceColor',fc,'FaceAlpha',fa,'EdgeColor',fc*0.8,'LineWidth',1); %#ok<AGROW>
    h(end+1)=text(ax,app.iv(k).t0,yl(2),[' ' bandLabel(app.iv(k))], ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontSize',9,'FontWeight',boldIf(isSel),'Color',[0.15 0.15 0.15], ...
        'Clipping','on'); %#ok<AGROW>
end
hold(ax,'off');
set(h,'PickableParts','none','HitTest','off');
app=getApp(fig); app.bandLayer=h; setApp(fig,app);
end

function s=bandLabel(iv)
%bandLabel  What the band says on the plot: its name, and its channel when it has
%   one - a window that belongs to one chamber has to look different from one that
%   belongs to the file.
s=char(iv.name);
if ~isempty(iv.channel), s=[s ' [' char(iv.channel) ']']; end
end

function w=boldIf(tf), if tf, w='bold'; else, w='normal'; end, end

function placeSelectionBox(fig)
%placeSelectionBox  The yellow box a new interval is dragged out with.  It exists
%   from the start rather than being drawn on demand, so adding an interval is one
%   drag and one press instead of a mode to enter first.  IT IS THE ONLY
%   INTERACTIVE OBJECT ON THIS PLOT, which is what makes it responsive: it is the
%   only thing a mouse event has to be hit-tested against.
app=getApp(fig); d=app.data; ax=app.c.axTrace;
sp=spanOf(d,app.opts);
if isempty(sp), return; end
t0=sp(1); span=max(sp(2)-sp(1),eps);
app.selRoi=drawrectangle(ax,'Position',[t0+0.02*span app.yl(1) 0.10*span diff(app.yl)], ...
    'Color',[0.95 0.75 0],'FaceAlpha',0.12,'LineWidth',1,'Label','new');
addlistener(app.selRoi,'ROIMoved',@(o,~) snapFullHeight(o,app.yl));
setApp(fig,app);
end

function parkSelectionBox(fig)
%parkSelectionBox  Put the yellow box to the RIGHT of the windows already defined, so
%   a recording that opens with three intervals does not open with the box sitting
%   under the first of them - the same rule Add applies after every new band.
app=getApp(fig);
if isempty(app.selRoi) || ~isvalid(app.selRoi), placeSelectionBox(fig); app=getApp(fig); end
if isempty(app.iv) || isempty(app.selRoi) || ~isvalid(app.selRoi), return; end
np=app.selRoi.Position;
xl=xlim(app.c.axTrace); w=np(3);
x=max([app.iv.t1]);
if x+w>xl(2), x=max(xl(1),xl(2)-w); end
np(1)=x; app.selRoi.Position=np;
end

function snapFullHeight(roi,yl)
%snapFullHeight  Keep the box spanning the whole axes: it selects a stretch of TIME,
%   and a box that is also dragged in y would say something about the signal value
%   that it does not mean.
p=roi.Position;
if numel(p)==4 && (p(2)~=yl(1) || p(4)~=diff(yl))
    roi.Position=[p(1) yl(1) p(3) diff(yl)];
end
end

%% ===================== the model ==================================== %%
function preloadIntervals(fig,opts)
%preloadIntervals  The intervals the recording already carries.  Over the ceiling (a
%   crop file that somehow holds two windows) only the first are re-offered, so the
%   editor never opens already over budget.
n=size(opts.intervals,1);
if n>opts.maxIntervals, n=opts.maxIntervals; end
app=getApp(fig);
for k=1:1:n
    nm='';
    if numel(opts.names)>=k, nm=char(opts.names{k}); end
    if isempty(nm), nm=sprintf('interval%d',k); end
    ch='';
    if numel(opts.channels)>=k, ch=firstKnownChannel(app,opts.channels{k}); end
    app.iv(end+1)=oneInterval(opts.intervals(k,1),opts.intervals(k,2),nm,ch);
end
app.iv=sortIntervals(app.iv);
setApp(fig,app);
clampAll(fig);
end

function ch=firstKnownChannel(app,want)
%firstKnownChannel  The channel a pre-loaded window is analysed on, MATCHED AGAINST
%   THE CHANNELS THIS RECORDING ACTUALLY HAS.  Empty, or a list covering every
%   channel, means all of them; anything else takes the first name the recording
%   answers to, and nothing when it answers to none.  The match is what earns the
%   function: a set of windows carried over from another recording of the same
%   protocol may name a chamber this one does not have, and showing it would offer a
%   window on a channel that cannot be drawn.
ch='';
if isempty(want) || ~isWmyo(app.opts), return; end
want=reshape(cellstr(want),1,[]);
names=reshape({app.data.channels.name},1,[]);
if numel(want)>=numel(names) && all(ismember(lower(strtrim(names)),lower(strtrim(want))))
    return                                    % it covered the whole recording
end
k=find(ismember(lower(strtrim(names)),lower(strtrim(want))),1);
if ~isempty(k), ch=names{k}; end
end

function iv=oneInterval(t0,t1,name,channel)
iv=struct('t0',min(t0,t1),'t1',max(t0,t1),'name',char(name),'channel',char(channel));
end

function iv=sortIntervals(iv)
%sortIntervals  By start time - the order the table, the outputs and every index the
%   API takes are all in.
if isempty(iv), return; end
[~,o]=sort([iv.t0]); iv=iv(o);
end

function clampAll(fig)
%clampAll  A WINDOW CANNOT LEAVE THE ANALYSED SPAN.  The trace is exactly what was
%   analysed - after a time crop it is only the cropped stretch - so a window outside
%   it selects frames that do not exist.  A nudge or a typed time that would take a
%   window past either end is clipped to that end rather than refused, so the gesture
%   still does something and the result is always a real window.
app=getApp(fig);
sp=spanOf(app.data,app.opts);
if isempty(sp) || isempty(app.iv), return; end
lo=sp(1); hi=sp(2);
for k=1:1:numel(app.iv)
    w=min(app.iv(k).t1-app.iv(k).t0,hi-lo);   % longer than the recording IS the recording
    t0=min(max(app.iv(k).t0,lo),hi-w);
    app.iv(k).t0=t0; app.iv(k).t1=t0+w;
end
setApp(fig,app);
end

function onAddPressed(fig)
%onAddPressed  Add the yellow box as an interval, then move the box to the RIGHT of
%   it.  Leaving the box under the band it has just become is what made adding a
%   second interval start with dragging the box out from underneath the first.
app=getApp(fig);
if numel(app.iv)>=app.opts.maxIntervals, return; end
if isempty(app.selRoi) || ~isvalid(app.selRoi), return; end
p=app.selRoi.Position;
if p(3)<=0, return; end
nm=strtrim(app.c.name.Value);
if isempty(nm), nm=sprintf('interval%d',numel(app.iv)+1); end
app.iv(end+1)=oneInterval(p(1),p(1)+p(3),nm,defaultChannel(app));
app.iv=sortIntervals(app.iv);
app.c.name.Value='';
xl=xlim(app.c.axTrace); w=max(p(3),eps); x=p(1)+p(3);
if x+w>xl(2), x=max(xl(1),xl(2)-w); end           % tuck it against the right edge
np=app.selRoi.Position; np(1)=x; np(3)=w; app.selRoi.Position=np;
setApp(fig,app);
clampAll(fig);
modelChanged(fig);
end

function ch=defaultChannel(~)
%defaultChannel  A NEW WINDOW IS ANALYSED ON EVERY CHANNEL until somebody says
%   otherwise.  It used to inherit whatever the channel list happened to hold, which
%   meant the analysis depended on what the operator was looking at when they
%   pressed Add.
ch='';
end

function onDeletePressed(fig)
%onDeletePressed  Remove the interval selected in the table.  It is a button and not
%   a mode: the bands are not clickable, so there is nothing to arm and nothing to
%   forget to disarm.
app=getApp(fig);
if app.sel<1 || app.sel>numel(app.iv), return; end
apiRemove(fig,app.sel);
end

function modelChanged(fig)
%modelChanged  The one call every edit ends with: the bands and the table are both
%   views of the model, so they are rebuilt together and can never disagree.  The
%   SIGNAL layer is deliberately not touched - nothing about an interval changes
%   what the recording looks like.
drawBands(fig);
refreshTable(fig);
end

%% ===================== the table ==================================== %%
function refreshTable(fig)
%refreshTable  The table is a VIEW of the model, rebuilt whenever an interval
%   appears, moves or goes.
app=getApp(fig);
wire=isWmyo(app.opts);
D=cell(numel(app.iv),3+wire);
for k=1:1:numel(app.iv)
    D(k,1:3)={app.iv(k).t0,app.iv(k).t1,app.iv(k).name};
    if wire
        if isempty(app.iv(k).channel), D{k,4}=allChannelsWord();
        else,                          D{k,4}=app.iv(k).channel;
        end
    end
end
app.c.table.Data=D;
if isvalid(app.c.addBtn)
    app.c.addBtn.Enable=onOff(numel(app.iv)<app.opts.maxIntervals);
end
if isvalid(app.c.delBtn)
    app.c.delBtn.Enable=onOff(app.sel>=1 && app.sel<=numel(app.iv));
end
setApp(fig,app);
showWarning(fig);
end

function onTableEdit(fig,e)
%onTableEdit  THIS IS HOW AN INTERVAL IS EDITED.  A boundary, a name or - on a wire
%   myograph - the channel it is analysed on.  An impossible edit (an end before its
%   start) is simply not applied: the table redraws from the model, which still holds
%   the last good value.
if isempty(e.Indices), return; end
app=getApp(fig); r=e.Indices(1); col=e.Indices(2);
if r<1 || r>numel(app.iv), return; end
switch col
    case 1
        a=num(e.NewData,NaN);
        if isfinite(a) && app.iv(r).t1>a, app.iv(r).t0=a; end
    case 2
        b=num(e.NewData,NaN);
        if isfinite(b) && b>app.iv(r).t0, app.iv(r).t1=b; end
    case 3
        nm=strtrim(char(string(e.NewData)));
        if ~isempty(nm), app.iv(r).name=nm; end
    case 4
        app.iv(r).channel=channelFromCell(app,e.NewData);
end
sel=app.iv(r);
app.iv=sortIntervals(app.iv);
app.sel=find([app.iv.t0]==sel.t0 & strcmp({app.iv.name},sel.name),1);
if isempty(app.sel), app.sel=0; end
setApp(fig,app);
clampAll(fig);
modelChanged(fig);
end

function ch=channelFromCell(app,v)
%channelFromCell  The channel column's value as the model stores it: '' for every
%   channel, and otherwise the recording's own spelling of the name.
ch='';
v=strtrim(char(string(v)));
if isempty(v) || strcmpi(v,allChannelsWord()), return; end
names=reshape({app.data.channels.name},1,[]);
k=find(strcmpi(strtrim(names),v),1);
if ~isempty(k), ch=names{k}; end
end

function onRowSelected(fig,e)
%onRowSelected  Selecting a row selects that interval: it is the one the arrow keys
%   move, the one Delete removes, and the one the preview loads.
if isempty(e.Indices), return; end
apiSelect(fig,e.Indices(1));
end

%% ===================== keys and validation ========================== %%
function onKey(fig,e)
%onKey  Nudge or resize the selected interval.  Ctrl steps 5 s instead of 1 s, Shift
%   resizes instead of moving.  THE BOUNDARY IS NOT ROUNDED: a start set to 12.35 s
%   in the table stays 12.35 s after a nudge.
%
%   NOT WHILE THE NAME FIELD HAS THE KEYBOARD.  A uifigure delivers key presses to
%   the window whether or not a component is being typed into, so without this the
%   left arrow that moves the caret back through a half-typed name would also move
%   the selected interval a second earlier.
app=getApp(fig);
if typingAName(fig,app), return; end
k=app.sel;
if k<1 || k>numel(app.iv), return; end
step=1;
if any(strcmpi(e.Modifier,'control')), step=5; end
resize=any(strcmpi(e.Modifier,'shift'));
t0=app.iv(k).t0; w=app.iv(k).t1-app.iv(k).t0;
switch lower(e.Key)
    case 'leftarrow'
        if resize, w=max(step,w-step); else, t0=t0-step; end
    case 'rightarrow'
        if resize, w=w+step; else, t0=t0+step; end
    otherwise
        return
end
sel=app.iv(k); sel.t0=t0; sel.t1=t0+w;
app.iv(k)=sel;
app.iv=sortIntervals(app.iv);
app.sel=find([app.iv.t0]==sel.t0 & strcmp({app.iv.name},sel.name),1);
if isempty(app.sel), app.sel=0; end
setApp(fig,app);
clampAll(fig);
modelChanged(fig);
end

function showWarning(fig)
%showWarning  Overlaps and gaps are REPORTED, not refused: analysing two windows that
%   share a minute is legitimate, and so is ignoring the hour between them.  What is
%   not legitimate is not noticing, so the sentence sits under the table and updates
%   with every edit.
app=getApp(fig);
if ~isvalid(app.c.warn), return; end
app.c.warn.Text=warningText(fig);
end

function msg=warningText(fig)
app=getApp(fig);
msg='';
iv=app.iv;
if isempty(iv), return; end
notes={};
if any(cellfun(@(x) isempty(strtrim(char(x))),{iv.name}))
    notes{end+1}='an interval has no name';
end
% OVERLAP IS PER CHANNEL.  Two windows on two different chambers may share a minute
% and mean nothing by it; two on the same channel are the case worth mentioning.  A
% window analysed on every channel is in every comparison, because it really is.
if anyOverlap(iv)
    notes{end+1}='two intervals on the same channel overlap';
end
so=sortrows([[iv.t0]' [iv.t1]']);
gaps=so(2:end,1)-so(1:end-1,2);
if any(gaps>1)
    n=nnz(gaps>1);
    if n==1
        notes{end+1}=sprintf('%.0f s of the recording falls between two intervals',max(gaps));
    else
        notes{end+1}=sprintf('%d stretches of the recording fall between intervals, the longest %.0f s', ...
            n,max(gaps));
    end
end
if isempty(notes), return; end
msg=[upper(notes{1}(1)) notes{1}(2:end)];
if numel(notes)>1, msg=[msg '; ' strjoin(notes(2:end),'; ')]; end
msg=[msg '.'];
end

function tf=typingAName(fig,app)
%typingAName  Is the operator in the middle of typing the next interval's name?
%   CurrentObject is the last component the mouse went into, which is exactly what
%   "has the keyboard" means in a uifigure.
tf=false;
try
    co=fig.CurrentObject;
    tf=~isempty(co) && isfield(app.c,'name') && isvalid(app.c.name) && isequal(co,app.c.name);
catch
end
end

function tf=anyOverlap(iv)
tf=false;
for a=1:1:numel(iv)
    for b=a+1:1:numel(iv)
        if ~sameChannelScope(iv(a),iv(b)), continue; end
        if iv(a).t0 < iv(b).t1-eps && iv(b).t0 < iv(a).t1-eps, tf=true; return; end
    end
end
end

function tf=sameChannelScope(a,b)
tf=isempty(a.channel) || isempty(b.channel) || strcmpi(a.channel,b.channel);
end

%% ===================== the channels and the comments ================ %%
function onShownChannelsChanged(fig)
%onShownChannelsChanged  Ticking a channel changes WHAT IS DRAWN and nothing else.
%   The intervals are untouched: which channel a window is analysed on is the
%   table's own column.
app=getApp(fig);
app.shown=selectedChannelNames(app.c);
setApp(fig,app);
drawSignal(fig);
end

function onCommentSelected(fig,e)
%onCommentSelected  Clicking a comment SCROLLS THE TRACE TO IT, keeping the zoom
%   the operator has already set: the point of the list is to find the event, and
%   re-fitting the whole recording every time would undo the zoom that made the
%   event findable.
if isempty(e.Indices), return; end
apiClickComment(fig,e.Indices(1));
end

function apiClickComment(fig,k)
app=getApp(fig); d=app.data; ax=app.c.axTrace;
if k<1 || k>numel(d.comments), return; end
tc=double(d.comments(k).time);
sp=spanOf(d,app.opts);
if isempty(sp), return; end
xl=xlim(ax); w=max(xl(2)-xl(1),eps);
if w>=sp(2)-sp(1), w=max((sp(2)-sp(1))/10,eps); end   % zoomed all the way out: zoom in
lo=max(sp(1),tc-w/2); hi=min(sp(2),lo+w);
lo=max(sp(1),hi-w);
setXLim(fig,[lo hi]);
end

function setXLim(fig,xl)
%setXLim  Move the view.
app=getApp(fig);
xlim(app.c.axTrace,xl);
end

%% ===================== the video preview ============================ %%
function setupPreview(fig,t0,t1)
%setupPreview  Load the frames of one interval, ready to step or play through.
previewStop(fig);
app=getApp(fig);
if isempty(app.data.video) || ~isfield(app.c,'axPrev'), return; end
if isempty(app.preview.vr)
    try
        app.preview.vr=VideoReader(app.data.video);
    catch
        app.preview.vr=[];
    end
end
if isempty(app.preview.vr), setApp(fig,app); return; end
fps=app.preview.vr.FrameRate;
if ~isfinite(fps) || fps<=0, fps=max(app.data.fs,eps); end
app.preview.times=(t0:1/fps:t1)';
if isempty(app.preview.times), app.preview.times=t0; end
app.preview.curIdx=1;
setApp(fig,app);
n=max(2,numel(app.preview.times));
app.c.slider.Limits=[1 n]; app.c.slider.Value=1;
renderPreviewFrame(fig);
end

function onSlide(fig,val)
app=getApp(fig);
if isempty(app.preview.times), return; end
app.preview.curIdx=max(1,min(numel(app.preview.times),round(val)));
setApp(fig,app);
renderPreviewFrame(fig);
end

function renderPreviewFrame(fig)
%renderPreviewFrame  One frame of the recording with the detected walls over it.
%   THE WALLS GO RED when that frame is flagged invalid - one wall had dilated out
%   of the field of view, so its diameter is a lower bound rather than a measurement,
%   and the operator has to see that while choosing the window.
%
%   THE FRAME AND ITS WALLS COME FROM ONE PLACE, getMyographWallFrame, so they cannot
%   describe different frames.  This window is where a person decides where the
%   analysis windows go, and it decides it on whether the edges were found - so the
%   walls drawn have to be the walls of the frame shown, which means asking the
%   function that knows which frame of the ORIGINAL recording a stored row came from
%   rather than snapping to the nearest sample of the windows laid end to end.  The
%   editor's OWN reader is handed over so no reader is opened per frame.  A recording
%   with no measured walls (a crop chosen before the diameter step) falls back to
%   showing the frame alone.
if ~isvalid(fig), return; end
app=getApp(fig); pv=app.preview;
if ~isfield(app.c,'axPrev') || ~isvalid(app.c.axPrev), return; end
if isempty(pv.vr) || isempty(pv.times), return; end
ax=app.c.axPrev;
t=pv.times(min(pv.curIdx,numel(pv.times)));
d=app.data;

W=[];
if ~isempty(d.wallFrame)
    try
        W=d.wallFrame(t,pv.vr);
    catch
        W=[];
    end
end
if isempty(W)
    try
        pv.vr.CurrentTime=min(max(t,0),max(0,pv.vr.Duration-1/max(pv.vr.FrameRate,eps)));
        W=struct('frame',readFrame(pv.vr),'left',[],'right',[],'rows',[],'valid',true);
    catch
        return
    end
end
if isempty(W.frame), return; end

cla(ax); image(ax,W.frame); axis(ax,'image'); ax.XTick=[]; ax.YTick=[];
ok=W.valid;
if ~isempty(W.left)
    hold(ax,'on');
    col=[0.1 0.9 0.1]; if ~ok, col=[0.95 0.15 0.15]; end
    plot(ax,W.left, W.rows,'-','Color',col,'LineWidth',1.2);
    plot(ax,W.right,W.rows,'-','Color',col,'LineWidth',1.2);
    hold(ax,'off');
end
ttl=sprintf('%.1f s   (frame %d of %d)',t,pv.curIdx,numel(pv.times));
if ~ok, ttl=[ttl '   wall out of view']; end
title(ax,ttl);
end

function previewPlay(fig)
%previewPlay  Step through the interval on a timer.  ONE timer at a time, and it is
%   the object releaseResources always looks for.
app=getApp(fig);
if isempty(app.preview.vr) || isempty(app.preview.times), return; end
if ~isempty(app.preview.timer) && isvalid(app.preview.timer), return; end
tm=timer('ExecutionMode','fixedRate','Period',0.08,'BusyMode','drop', ...
    'Name','myographIntervalPreview','TimerFcn',@(~,~) previewTick(fig));
app.preview.timer=tm; setApp(fig,app);
start(tm);
end

function previewStop(fig)
app=getApp(fig);
if ~isstruct(app) || ~isfield(app,'preview'), return; end
tm=app.preview.timer;
if ~isempty(tm) && isvalid(tm)
    stop(tm); delete(tm);
end
app.preview.timer=[]; setApp(fig,app);
end

function previewTick(fig)
%previewTick  Advance by the chosen number of frames and wrap at the end.  A tick
%   that finds the window gone stops itself: a timer outliving its figure is exactly
%   how this editor would leak.
if ~isvalid(fig)
    try
        delete(timerfind('Name','myographIntervalPreview'));
    catch
    end
    return
end
app=getApp(fig);
if isempty(app.preview.times), previewStop(fig); return; end
step=2;
try
    step=max(1,round(app.c.step.Value));
catch
end
n=numel(app.preview.times);
idx=app.preview.curIdx+step; if idx>n, idx=1; end
app.preview.curIdx=idx; setApp(fig,app);
try
    app.c.slider.Value=min(idx,app.c.slider.Limits(2));
catch
end
renderPreviewFrame(fig);
end

%% ===================== finishing ==================================== %%
function requestDone(fig)
%requestDone  The one way out: set the sentinel the blocking waitfor watches.  The
%   window is NOT deleted here - the caller still has to read the intervals off it.
%
%   AND THE ONE PLACE THE OPERATOR IS ASKED BEFORE ANYTHING IS THROWN AWAY.  What
%   this window returns REPLACES the measurement: the stretches left outside every
%   window are not written back, and no result can bring them back - only measuring
%   the recording again can.  So Done asks, once, and only when there really is
%   something to lose.  There is no setting that turns it off.
if ~isvalid(fig), return; end
if ~confirmDiscard(fig), return; end
fig.UserData='done';
end

function ok=confirmDiscard(fig)
%confirmDiscard  Ask when these windows would discard measured data, and say how
%   much.  TWO CASES ARE SILENT AND BOTH ARE HONEST.  Nothing to lose is the
%   ordinary one - a window chosen before anything was measured, or one that still
%   covers everything - and there is nothing to warn about.  A window nobody can see
%   is the other: an invisible editor is being driven by a script rather than by a
%   person, and a modal dialog on it would wait for an answer that never comes.  The
%   wrapper's own refusal is what guards that path, and it guards it before the
%   editor ever opens.
ok=true;
app=getApp(fig);
[out,total]=discardedSamples(app);
if out<=0 || total<=0, return; end
if ~strcmpi(char(fig.Visible),'on'), return; end
word='frames'; if isWmyo(app.opts), word='samples'; end
msg=sprintf(['These windows leave out %d of the %d %s already measured, and only ' ...
    'the recording itself can bring them back.\n\nKeep the windows as they are?'], ...
    out,total,word);
try
    sel=uiconfirm(fig,msg,'Some of the measurement will be thrown away', ...
        'Options',{'Keep the windows','Go back'},'DefaultOption',2, ...
        'CancelOption',2,'Icon','warning');
catch
    return          % no dialog to be had: do not trap the operator in the window
end
ok=strcmp(sel,'Keep the windows');
end

function [out,total]=discardedSamples(app)
%discardedSamples  How much of what is already measured falls OUTSIDE the windows as
%   they now stand - which is exactly what saving them would throw away.  Counted on
%   the samples that EXIST rather than on the seconds a window covers: a stretch that
%   was never measured is not lost by being left out of a window.
%   NO WINDOWS AT ALL MEANS THE WHOLE RECORDING, the same rule the re-cut applies, so
%   it discards nothing rather than everything.
out=0; total=0;
if ~isstruct(app) || isempty(app.iv), return; end
if ~isfield(app.data,'measurement') || ~app.data.measurement, return; end
if isWmyo(app.opts)
    for i=1:1:numel(app.data.channels)
        t=app.data.channels(i).time;
        if isempty(t), continue; end
        total=total+numel(t);
        out=out+nnz(~coveredBy(t,app.iv,app.data.channels(i).name));
    end
    return
end
t=app.data.time;
if isempty(t), return; end
total=numel(t);
out=nnz(~coveredBy(t,app.iv,''));
end

function tf=coveredBy(t,iv,chanName)
%coveredBy  Which of these times fall inside some window.  A window assigned to one
%   channel covers only that channel's samples; one left on 'all channels' covers
%   every channel's, which is what it means.
tf=false(numel(t),1);
for k=1:1:numel(iv)
    if ~isempty(chanName) && ~isempty(iv(k).channel) && ...
            ~strcmpi(strtrim(iv(k).channel),strtrim(chanName))
        continue
    end
    tf=tf | (t>=iv(k).t0 & t<=iv(k).t1);
end
end

function [ivT,names]=collect(fig)
%collect  The intervals as they stand, sorted by start time.
ivT=zeros(0,2); names={};
if ~isvalid(fig), return; end
app=getApp(fig);
for i=1:1:numel(app.iv)
    ivT(end+1,:)=[app.iv(i).t0 app.iv(i).t1]; %#ok<AGROW>
    names{end+1}=app.iv(i).name; %#ok<AGROW>
end
end

function releaseResources(fig)
%releaseResources  The playback timer, the open recording and the window itself,
%   released however this editor ended.  A timer left running at a deleted figure
%   and a VideoReader left holding the file are the two ways this window leaks.
if ~isvalid(fig), return; end
try
    previewStop(fig);
catch
end
try
    app=getApp(fig);
    if isstruct(app) && isfield(app,'preview')
        app.preview.vr=[]; setApp(fig,app);
    end
catch
end
delete(fig);
end

function e=emptyExtra()
%emptyExtra  The per-interval side information.  The pressure myograph has none;
%   the wire myograph fills one element per interval with the channel it is analysed
%   on, which is why the output exists now rather than being invented then.
e=struct('channels',{});
end

function e=collectExtra(fig)
%collectExtra  The channel each interval is analysed on, IN THE SAME ORDER the
%   intervals come back in - sorted by start time - so extra(k) belongs to ivT(k,:)
%   and names{k} without anything having to pair them up afterwards.  It is a
%   cellstr because that is what the recut and the result tree store: EMPTY means
%   every channel, one name means that channel.
e=emptyExtra();
if ~isvalid(fig), return; end
app=getApp(fig);
if ~isWmyo(app.opts), return; end
for i=1:1:numel(app.iv)
    if isempty(app.iv(i).channel), e(i).channels={};
    else,                          e(i).channels={app.iv(i).channel};
    end
end
end

%% ===================== the programmatic API ========================= %%
function api=buildAPI(fig)
%buildAPI  The same callbacks the buttons use, as a struct of handles, so the editor
%   can be driven without gestures (tests, and any future scripted use).  Every index
%   is the row of the table - i.e. the intervals sorted by start time.
api=struct( ...
    'add',       @(t0,t1,nm) apiAdd(fig,t0,t1,nm), ...
    'remove',    @(k) apiRemove(fig,k), ...
    'rename',    @(k,nm) apiRename(fig,k,nm), ...
    'editCell',  @(k,col,v) apiEditCell(fig,k,col,v), ...
    'select',    @(k) apiSelect(fig,k), ...
    'selected',  @() apiSelected(fig), ...
    'press',     @(key,mods) apiPress(fig,key,mods), ...
    'intervals', @() apiIntervals(fig), ...
    'names',     @() apiNames(fig), ...
    'count',     @() numel(apiModel(fig)), ...
    'warning',   @() warningText(fig), ...
    'preview',   @(k) apiPreview(fig,k), ...
    'play',      @() previewPlay(fig), ...
    'stop',      @() previewStop(fig), ...
    'frame',     @() apiFrame(fig), ...
    'setShownChannels',@(nm) apiSetShownChannels(fig,nm), ...
    'shownChannels',@() apiShownChannels(fig), ...
    'setIntervalChannel',@(k,nm) apiSetIntervalChannel(fig,k,nm), ...
    'channels',  @() apiChannels(fig), ...
    'clickComment',@(k) apiClickComment(fig,k), ...
    'xlim',      @() apiXLim(fig), ...
    'setXlim',   @(xl) setXLim(fig,xl), ...
    'drawnPoints',@() apiDrawnPoints(fig), ...
    'done',      @() requestDone(fig), ...
    'getApp',    @() getApp(fig));
end

function iv=apiModel(fig)
app=getApp(fig); iv=app.iv;
end

function apiSetShownChannels(fig,nm)
%apiSetShownChannels  The channel list, set as a person sets it: pick the items,
%   which goes on to fire the same callback the mouse does.
app=getApp(fig);
if ~isfield(app.c,'chList') || ~isvalid(app.c.chList), return; end
if isempty(nm), nm={}; else, nm=reshape(cellstr(nm),1,[]); end
app.c.chList.Value=intersect(nm,app.c.chList.Items,'stable');
setApp(fig,app);
onShownChannelsChanged(fig);
end
function nm=apiShownChannels(fig)
app=getApp(fig); nm=app.shown;
end
function apiSetIntervalChannel(fig,k,nm)
%apiSetIntervalChannel  The table's channel cell, set as a person sets it.
if isempty(nm), nm=allChannelsWord(); end
if iscell(nm), if isempty(nm), nm=allChannelsWord(); else, nm=nm{1}; end, end
apiEditCell(fig,k,4,nm);
end
function ch=apiChannels(fig)
%apiChannels  The channel of every interval, in table order, as the recut stores it.
e=collectExtra(fig);
if isempty(e), ch={}; else, ch={e.channels}; end
end
function xl=apiXLim(fig)
app=getApp(fig); xl=xlim(app.c.axTrace);
end
function n=apiDrawnPoints(fig)
%apiDrawnPoints  How many points the signal layer actually put on screen.  It is the
%   number the decimation exists to hold down, so it is the number a test asserts on.
app=getApp(fig); n=0;
for i=1:1:numel(app.sigLayer)
    h=app.sigLayer(i);
    if isvalid(h) && isprop(h,'XData'), n=n+numel(h.XData); end
end
end

function apiAdd(fig,t0,t1,nm)
%apiAdd  The Add button: put the selection box where asked, name it, press Add.
app=getApp(fig);
if ~isempty(app.selRoi) && isvalid(app.selRoi)
    app.selRoi.Position=[min(t0,t1) app.yl(1) max(abs(t1-t0),eps) diff(app.yl)];
end
if nargin>=4 && ~isempty(nm), app.c.name.Value=char(nm); end
setApp(fig,app);
onAddPressed(fig);
end
function apiRemove(fig,k)
app=getApp(fig);
if k<1 || k>numel(app.iv), return; end
app.iv(k)=[];
if app.sel==k, app.sel=0; elseif app.sel>k, app.sel=app.sel-1; end
setApp(fig,app);
modelChanged(fig);
end
function apiRename(fig,k,nm)
apiEditCell(fig,k,3,nm);
end
function apiEditCell(fig,k,col,v)
e=struct('Indices',[k col],'NewData',v);
onTableEdit(fig,e);
end
function apiSelect(fig,k)
app=getApp(fig);
if k<1 || k>numel(app.iv), app.sel=0; setApp(fig,app); modelChanged(fig); return; end
app.sel=k;
setApp(fig,app);
modelChanged(fig);
if ~isWmyo(app.opts)
    setupPreview(fig,app.iv(k).t0,app.iv(k).t1);
end
end
function k=apiSelected(fig)
app=getApp(fig); k=app.sel;
end
function apiPress(fig,key,mods)
if nargin<3 || isempty(mods), mods={}; else, mods=reshape(cellstr(mods),1,[]); end
e=struct('Key',char(key),'Modifier',{mods});
onKey(fig,e);
end
function ivT=apiIntervals(fig), ivT=collect(fig); end
function nm=apiNames(fig), [~,nm]=collect(fig); end
function apiPreview(fig,k), apiSelect(fig,k); end
function i=apiFrame(fig)
app=getApp(fig); i=app.preview.curIdx;
end

%% ===================== small helpers ================================ %%
function app=getApp(fig), app=getappdata(fig,'app'); end
function setApp(fig,app), setappdata(fig,'app',app); end
function v=onOff(tf), if tf, v='on'; else, v='off'; end, end
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
function n=num(v,dflt)
if isnumeric(v) && ~isempty(v), n=double(v(1));
elseif ischar(v) || isstring(v), n=str2double(v);
else, n=dflt;
end
if isempty(n) || isnan(n), n=dflt; end
end
