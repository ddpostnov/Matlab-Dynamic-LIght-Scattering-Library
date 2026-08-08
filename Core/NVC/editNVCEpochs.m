%editNVCEpochs  Look at every stimulus repetition of one recording and keep or drop it
%
%   epochTrust = editNVCEpochs(ed,opts) opens ONE window on ONE recording, shows what
%   the analysis made of every stimulus repetition in it, lets the operator change any
%   of those decisions, and hands the set back when Done is pressed.  It is a CORE GUI
%   function - it takes data and returns a decision, opens no file and narrates nothing
%   - so the wrapper that calls it is the one that loops the recordings, exactly the way
%   setMyographIntervals calls editMyographIntervals.  There is deliberately NO file
%   dropdown: everything is decided for one recording, then the next window opens.
%
%   IT COMPUTES NOTHING.  Every number on the screen was measured by runNVC before this
%   window opened - the responding area, the two confidence summaries, the factors - and
%   this file only draws them and records what the operator did about them.  That is the
%   same split setRegions makes between s.nRegions and results.regionsMask, and it is why
%   the input is ONE struct: `ed` is also what the '_rep_nvccuts.jpg' page draws from, so
%   the page and this window cannot drift apart.
%
%   WHAT A REPETITION BEING KEPT MEANS.  The recording-wide trusted set is what the
%   metrics tables average over and what the representative repetition is built from, so
%   dropping one here removes it from every per-segment number the step writes.  It does
%   NOT change the per-repetition tree: every repetition is measured and stored either
%   way, which is what makes this window a decision rather than a filter.
%
%   WHAT THE WINDOW HOLDS
%     * the recording's mean trace, with every repetition's span shaded.  A KEPT
%       repetition is drawn in full colour and a dropped one is left faded, so the state
%       reads at a glance and it reads without relying on the difference between a green
%       tint and a red one;
%     * a TABLE, one row per repetition: how much of the imaged field responded, the two
%       confidence summaries of its median segment, WHICH CHECK THAT SEGMENT SCORED
%       LOWEST ON - in words - and a keep box.  The box is the edit;
%     * clicking the trace flips the repetition under the pointer, for the same gesture
%       the old picker had;
%     * Keep all / Drop all, and Restore, which puts the analysis's own proposal back;
%     * the dropped-frame series under the trace WHEN THERE IS ANY.  A flat zero line
%       across a good recording is clutter that teaches a reader to stop looking at it;
%     * DONE.  Closing the window means Done.
%
%   AND DONE ASKS WHEN NOTHING IS LEFT.  An empty set is not an error - the per-repetition
%   tree is the product and it is complete either way - but it makes every metrics-table
%   column NaN, and with the representative repetition switched on it stops the run with
%   nothing to average.  So Done asks once, and only then.  A window nobody can see is
%   silent: an invisible editor is being driven by a script, and a modal dialog on it
%   would wait for an answer that never comes.
%
%   WHICH CHECK WAS WEAKEST, AND WHICH QUESTION THAT ANSWERS.  The table names the factor
%   with the LOWEST MEDIAN across segments for that repetition, which is the median
%   segment's weakest check and is comparable with the confidence medians beside it.  It
%   is NOT the check that cancelled the most segments - those are different questions and
%   they have different answers - and '_rep_nvcconfidence.jpg' is where the second one is
%   answered, as a share per repetition.
%
%   HOW MANY CHECKS THERE ARE IS NOT A CONSTANT: nine by default and fourteen with the
%   across-repetition group switched on, and the timing family is named from
%   s.nvcAreaPcts.  So the factors are read from ed.factorNames and their sentences from
%   ed.factorWords, and no count and no name is written here.
%
%   Syntax:
%      epochTrust = editNVCEpochs(ed)
%      epochTrust = editNVCEpochs(ed,opts)
%
%   INPUTS
%     ed     what is drawn - the struct runNVC builds in its `editorData` subfunction.
%            Times are seconds on the RECORDING clock, so a span can always be located
%            back in the recording.
%       .signal        which trace the responding-area rule was applied to, for the title
%       .trace         [nT x 1] the recording's field-mean trace of that signal
%       .traceLabel    what it is a trace of, for the y-axis
%       .time          [nT x 1] the recording clock, seconds
%       .epochStart    [nEp x 1] where each repetition was cut, on that clock
%       .epochEnd      [nEp x 1] where each one ends
%       .timeLossSeries[nT x 1] per-sample dropped time, seconds
%       .epochTrust    [nEp x 1] logical - THE ANALYSIS'S PROPOSAL, and what comes back
%       .areaFrac      [nEp x 1] the share of the segmented field that responded
%       .conf .confMin [nSeg x nEp] the two confidence numbers of every pair
%       .trust         [nSeg x nEp] logical, which pairs cleared both thresholds
%       .factors       struct of [nSeg x nEp] arrays, one per check
%       .factorNames   {1 x nF} their order
%       .factorWords   {1 x nF} what a low score on each of them means, in words
%     opts   how it opens.
%       .title      the window title (default 'Stimulus repetitions')
%       .readyFcn   (optional) @(fig) called once, with the window built and drawn,
%                   BEFORE the blocking wait.  It is how the editor is driven without
%                   gestures (see the API below) - a hook that edits and presses Done
%       .visible    'on' (default) | 'off'
%
%   OUTPUTS
%     epochTrust  [nEp x 1] logical - the repetitions the operator kept.  It is what the
%                 caller assigns back over its own proposal, and everything downstream of
%                 the caller then reads the operator's set and not the rule's.
%
%   DRIVING IT WITHOUT GESTURES.  getappdata(fig,'nvcEpochsAPI') is a struct of handles
%   onto the same callbacks the widgets use - toggle / set / setAll / restore / select /
%   state / trust / summary / detail / done - so a test flips repetitions exactly as a
%   person would and reads back what comes out.  Same arrangement as guiExplore,
%   guiExport, guiResponse and editMyographIntervals.  It is reached through
%   opts.readyFcn, which is this window's own contract: nothing of it rides in a
%   wrapper's s.
%
%   IT RETURNS THE DECISION AND NOT THE WINDOW.  The figure is deleted on the way out
%   however this function ends, so a handle handed back would be a dead one; the caller
%   waits for Done and reads the set, which is the same contract editMyographIntervals
%   has.  A headless caller reaches the live figure through opts.readyFcn instead.
%
%   THE BLOCKING WAIT IS waitfor ON A PRIVATE SENTINEL, never uiwait: a modal dialog
%   inside this window runs its own uiwait/uiresume, and that inner uiresume would also
%   release an outer uiwait - the editor would return the moment Done asked a question.
%   editMyographIntervals and setRegions>editRegions document the identical trap.
%
%   THE TRACE IS DRAWN ONCE AND ONLY THE OVERLAY MOVES.  A recording is thousands of
%   samples, not the tens of millions a LabChart file carries, so there is no decimation
%   here and it was measured rather than assumed (see the session note).  What toggling a
%   repetition does touch is one line's data and one shaded span's colour, never the
%   trace - which is what makes twenty repetitions flip without a redraw.
%
%   EXAMPLE
%     ed = ...;                                   % runNVC>editorData builds it
%     keep = editNVCEpochs(ed,struct('title','Repetitions - Mouse1_t_BFI'));
%
% See also: runNVC, getNVCConfidence, getNVCEpochTrust, editMyographIntervals,
%           setMyographIntervals, guiResponse, setRegions
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 06-August-2026

%------------- BEGIN CODE --------------
function epochTrust = editNVCEpochs(ed,opts)

if nargin<2, opts=struct(); end
opts=withOptsDefaults(opts);
ed=withDataDefaults(ed);

fig=buildWindow(ed,opts);
% The window is released on the way out however this function ends - Done, a closed
% window, or an error inside a callback.
guard=onCleanup(@() releaseWindow(fig));

drawTrace(fig);                  % once; the state layers are drawn over it
refreshAll(fig);

% THE API GOES LIVE ONLY ONCE THERE IS SOMETHING TO DRIVE, and that is not fussiness.
% Building a uifigure runs the event queue, so anything already waiting on a timer can
% find this window - by name, before it has been drawn - the moment it exists.  Reached
% then, a toggle would index a layer that has not been made yet.  Registered here, the
% window is complete before it can be found.
setappdata(fig,'nvcEpochsAPI',buildAPI(fig));

if ~isempty(opts.readyFcn)
    opts.readyFcn(fig);          % headless: the caller drives the API and presses Done
end

waitfor(fig,'UserData','done');  % NOT uiwait - see the header

epochTrust=collect(fig);
end

%% ===================== defaults ===================================== %%
function d=withDataDefaults(ed)
%withDataDefaults  The drawn data, in the one shape the rest of the file assumes.
%   EVERY FIELD IS REQUIRED because every one of them is already computed by the caller:
%   this window measures nothing, so a missing field is a producer that stopped building
%   the contract rather than an option somebody did not pass.
if ~isstruct(ed)
    error('editNVCEpochs:badData','ed must be the struct runNVC builds for this window.');
end
need={'trace','time','epochStart','epochEnd','timeLossSeries','epochTrust','areaFrac', ...
      'conf','confMin','trust','factors','factorNames'};
miss=need(~isfield(ed,need));
if ~isempty(miss)
    error('editNVCEpochs:badData', ...
        'ed carries no %s.  It is the struct runNVC>editorData builds, whole.', ...
        strjoin(miss,' / '));
end
d=ed;
d.time      =double(d.time(:));
d.trace     =double(d.trace(:));
d.epochStart=double(d.epochStart(:));
d.epochEnd  =double(d.epochEnd(:));
d.areaFrac  =double(d.areaFrac(:));
d.epochTrust=logical(d.epochTrust(:));
d.timeLossSeries=double(d.timeLossSeries(:));
if numel(d.trace)~=numel(d.time)
    error('editNVCEpochs:badData','ed.time and ed.trace must be the same length.');
end
if numel(d.epochEnd)~=numel(d.epochStart) || numel(d.areaFrac)~=numel(d.epochStart) || ...
        numel(d.epochTrust)~=numel(d.epochStart)
    error('editNVCEpochs:badData', ...
        'ed.epochStart, .epochEnd, .areaFrac and .epochTrust must be one per repetition.');
end
if numel(d.timeLossSeries)~=numel(d.time), d.timeLossSeries=zeros(size(d.time)); end
d.factorNames=reshape(cellstr(d.factorNames),1,[]);
% THE SENTENCES ARE THE CORE'S and they arrive with the factors.  A caller that built ed
% before getNVCConfidence published them still gets a window, with the check named rather
% than explained - the one thing this file will not do is invent a second phrasebook.
if ~isfield(d,'factorWords') || numel(d.factorWords)~=numel(d.factorNames)
    d.factorWords=d.factorNames;
else
    d.factorWords=reshape(cellstr(d.factorWords),1,[]);
end
for f={'signal','traceLabel'}
    if ~isfield(d,f{1}) || isempty(d.(f{1})), d.(f{1})=''; else, d.(f{1})=char(d.(f{1})); end
end
end

function o=withOptsDefaults(opts)
%withOptsDefaults  How the window opens, with everything optional defaulted.
o=opts;
def=struct('title','Stimulus repetitions','readyFcn',[],'visible','on');
fn=fieldnames(def);
for i=1:1:numel(fn)
    if ~isfield(o,fn{i}) || isempty(o.(fn{i})), o.(fn{i})=def.(fn{i}); end
end
o.title=char(o.title);
o.visible=char(o.visible);
end

%% ===================== the window =================================== %%
function fig=buildWindow(ed,opts)
%buildWindow  The trace on the left, every decision on the right.
%
%   THE DROPPED-FRAME PANEL EXISTS ONLY WHEN FRAMES WERE DROPPED, and it is its own axes
%   rather than a second ruler on the trace's.  A right-hand y-axis is scaled from its own
%   data and only one side of a yyaxis pair answers `axis tight`, so the two curves would
%   be sharing a box and disagreeing about it; two axes in a grid row cannot.
fig=uifigure('Name',opts.title,'Visible',opts.visible, ...
    'Position',figurePosition(),'Color','w');
fig.UserData='editing';                          % the private wait sentinel
fig.CloseRequestFcn=@(~,~) requestDone(fig);     % closing the window == Done

gl=uigridlayout(fig,[1 2],'ColumnWidth',{'1.55x','1x'}, ...
    'Padding',[8 8 8 8],'ColumnSpacing',10,'BackgroundColor','w');

% ---- left: the trace, and the losses under it when there are any -------------
hasLoss=any(ed.timeLossSeries>0);
if hasLoss, rows={'3x','1x'}; else, rows={'1x'}; end
lp=uigridlayout(gl,[numel(rows) 1],'RowHeight',rows,'RowSpacing',6, ...
    'Padding',[0 0 0 0],'BackgroundColor','w');
c.axTrace=uiaxes(lp);
c.axTrace.ButtonDownFcn=@(~,~) onAxesClicked(fig);
c.axLoss=gobjects(0);
if hasLoss, c.axLoss=uiaxes(lp); end

% ---- right: what was decided, and what the operator decides instead ----------
rp=uigridlayout(gl,[4 1],'RowHeight',{'fit','1x','fit','fit'},'RowSpacing',7, ...
    'Padding',[0 0 0 0],'BackgroundColor','w');

uilabel(rp,'Text',hintText(),'WordWrap','on','FontSize',11);

% THE ROW NUMBER IS THE REPETITION NUMBER, so there is no column for it: the numbers on
% the shaded spans are the same numbers, and a column repeating the gutter beside it was
% one identity written twice.  THE WIDTHS ARE FIXED rather than 'auto' or 'fit' because
% the keep box is the only editable cell in the window and a sentence sized to its own
% content pushed it off the right-hand edge - an edit control nobody can reach is a worse
% failure than a truncated sentence, and the sentence is repeated in full underneath.
c.table=uitable(rp,'ColumnName', ...
    {'responding','confidence','lowest','weakest check','keep'}, ...
    'ColumnFormat',{'char','char','char','char','logical'}, ...
    'ColumnEditable',[false false false false true], ...
    'ColumnWidth',{78,78,55,240,42}, ...
    'Data',cell(0,5),'CellEditCallback',@(~,e) onTableEdit(fig,e), ...
    'CellSelectionCallback',@(~,e) onRowSelected(fig,e));

c.detail=uilabel(rp,'Text','','WordWrap','on','FontSize',11, ...
    'FontColor',[0.30 0.30 0.35]);

% DONE SITS APART from the three that change the set - it is the only one that ends the
% window, and a stretching gap is what says so without a separator line.
br=uigridlayout(rp,[1 5],'ColumnWidth',{'fit','fit','fit','1x','fit'}, ...
    'ColumnSpacing',5,'Padding',[0 0 0 0],'BackgroundColor','w');
c.allBtn=uibutton(br,'Text','Keep all','ButtonPushedFcn',@(~,~) apiSetAll(fig,true));
c.noneBtn=uibutton(br,'Text','Drop all','ButtonPushedFcn',@(~,~) apiSetAll(fig,false));
c.restoreBtn=uibutton(br,'Text','Restore','ButtonPushedFcn',@(~,~) apiRestore(fig), ...
    'Tooltip','put the analysis'' own decisions back');
c.doneBtn=uibutton(br,'Text','Done','FontWeight','bold', ...
    'BackgroundColor',[0.82 0.92 0.82],'ButtonPushedFcn',@(~,~) requestDone(fig));
c.doneBtn.Layout.Column=5;

% THE MODEL IS THE LOGICAL VECTOR and nothing else.  The shaded spans are drawn FROM it
% and the table is a view OF it, so there is never a graphics object to ask what the
% operator decided - the trap editMyographIntervals documents for its own bands.
app=struct('ed',ed,'opts',opts,'c',c,'keep',ed.epochTrust,'sel',0, ...
    'hSpan',gobjects(0),'hNum',gobjects(0),'hKept',gobjects(0),'yl',[0 1]);
setApp(fig,app);
end

function t=hintText()
%hintText  What this window is for and what to do in it, in the order it is done.
t=['Every stimulus repetition of this recording, and what the analysis made of it.  ' ...
   'Untick a repetition to leave it out of the per-segment results, or tick one back ' ...
   'in; clicking the trace does the same.  The repetitions kept here are the ones ' ...
   'the metrics are averaged over.  Press Done when ready.'];
end

function pos=figurePosition()
%figurePosition  A wide window, centred, that fits on a laptop screen.
sz=get(0,'ScreenSize');
w=min(1400,round(sz(3)*0.9)); h=min(760,round(sz(4)*0.80));
pos=[max(1,round((sz(3)-w)/2)) max(41,round((sz(4)-h)/2)) w h];
end

%% ===================== the drawing ================================== %%
function drawTrace(fig)
%drawTrace  The recording, once.  Everything that changes when a repetition is flipped is
%   a LAYER over this, so a toggle never re-renders the trace.
%
%   THE FADED CURVE IS THE WHOLE RECORDING and the strong one is only the samples inside
%   the repetitions being kept, broken over the gaps with NaN.  That is the "visibly
%   different state" this window exists for, and it survives being printed in grey: a
%   dropped repetition looks switched off rather than differently coloured.
app=getApp(fig);
ax=app.c.axTrace;
ed=app.ed;
cla(ax,'reset')
ax.ButtonDownFcn=@(~,~) onAxesClicked(fig);
hold(ax,'on')

lo=min(ed.trace); hi=max(ed.trace);
if ~(isfinite(lo)&&isfinite(hi)&&hi>lo), lo=0; hi=1; end
yl=[lo-0.04*(hi-lo) hi+0.10*(hi-lo)];

% the shaded spans first, so the curve is drawn over them
nEp=numel(ed.epochStart);
app.hSpan=gobjects(1,nEp);
app.hNum =gobjects(1,nEp);
for k=1:1:nEp
    app.hSpan(k)=xregion(ax,ed.epochStart(k),ed.epochEnd(k), ...
        'FaceColor',keptColour(true),'FaceAlpha',0.12,'HitTest','off');
    app.hNum(k)=text(ax,0.5*(ed.epochStart(k)+ed.epochEnd(k)),yl(2), ...
        sprintf('%d',k),'HorizontalAlignment','center','VerticalAlignment','top', ...
        'FontSize',9,'Color',[0.35 0.35 0.40],'HitTest','off','PickableParts','none');
end

% NOTHING DRAWN HERE TAKES THE CLICK.  The axes is what reports where the pointer was, so
% every curve, span and label is transparent to the mouse; a line that swallowed the
% press would make the gesture work everywhere except on the data.
plot(ax,ed.time,ed.trace,'-','Color',[0.72 0.74 0.78],'LineWidth',0.75, ...
    'HitTest','off','PickableParts','none');
app.hKept=plot(ax,NaN,NaN,'-','Color',[0.20 0.35 0.60],'LineWidth',1.2, ...
    'HitTest','off','PickableParts','none');

hold(ax,'off')
ylabel(ax,ed.traceLabel)
% the clock is labelled ONCE, under whichever panel is the bottom one
if isempty(app.c.axLoss), xlabel(ax,'Time, s'); end
% THE LIMITS COME FROM THE DATA and are not read back out of the axes: a uiaxes does not
% resolve `axis tight` before the next line runs, so a query here can still answer with
% the default on one end.
xlim(ax,[ed.time(1) ed.time(end)])
ylim(ax,yl)
app.yl=yl;
setApp(fig,app);

drawLosses(fig);
end

function drawLosses(fig)
%drawLosses  The seconds the camera did not deliver, under the trace, on the same x.
%   Drawn only when something was lost - a flat zero line across a good recording is a
%   panel that teaches a reader to stop looking at this one.
app=getApp(fig);
if isempty(app.c.axLoss) || ~isgraphics(app.c.axLoss), return; end
ax=app.c.axLoss;
cla(ax,'reset')
plot(ax,app.ed.time,app.ed.timeLossSeries,'-','Color',[0.85 0.20 0.15])
xlim(ax,[app.ed.time(1) app.ed.time(end)])
ylabel(ax,'Time lost, s')
xlabel(ax,'Time, s')
end

function refreshAll(fig)
%refreshAll  The one call every edit ends with: the spans, the strong curve, the table,
%   the title and the sentence under it are all views of the same logical vector, so they
%   are rebuilt together and cannot disagree.
refreshSpans(fig);
refreshTable(fig);
refreshText(fig);
end

function refreshSpans(fig)
%refreshSpans  The state layer: one colour per span and one line's data for the whole
%   kept set.  THIS IS WHAT MAKES A TOGGLE CHEAP - the trace underneath is untouched.
app=getApp(fig);
ed=app.ed;
x=[]; y=[];
for k=1:1:min(numel(app.keep),numel(app.hSpan))
    if ~isgraphics(app.hSpan(k)), continue; end
    app.hSpan(k).FaceColor=keptColour(app.keep(k));
    app.hSpan(k).FaceAlpha=0.12+0.10*(app.sel==k);
    if isgraphics(app.hNum(k))
        app.hNum(k).FontWeight=boldIf(app.sel==k);
    end
    if ~app.keep(k), continue; end
    in=ed.time>=ed.epochStart(k) & ed.time<=ed.epochEnd(k);
    % a NaN between the stretches, so one line object carries every kept repetition and
    % the pen is not dragged across the ones that were dropped
    x=[x; ed.time(in); NaN];      %#ok<AGROW>
    y=[y; ed.trace(in); NaN];     %#ok<AGROW>
end
if isempty(x), x=NaN; y=NaN; end
if isgraphics(app.hKept), set(app.hKept,'XData',x,'YData',y); end
setApp(fig,app);
end

function c=keptColour(tf)
%keptColour  The cuts page's own two colours, so the window and the page agree.
if tf, c=[0.15 0.60 0.25]; else, c=[0.85 0.20 0.15]; end
end

function w=boldIf(tf), if tf, w='bold'; else, w='normal'; end, end

function k=keepColumn(), k=5; end

%% ===================== the table and the words ====================== %%
function refreshTable(fig)
%refreshTable  The table is a VIEW of the model, rebuilt whenever a decision changes.
%   THE CELL CARRIES THE SENTENCE AND NOT THE CHECK'S NAME.  'fT90' on screen is a field
%   name leaking into a window a biologist reads; the sentence is what an operator can act
%   on, and it is the core's, not this file's.
app=getApp(fig);
n=numel(app.keep);
D=cell(n,5);
for k=1:1:n
    [~,words]=weakestOf(app.ed,k);
    D(k,:)={sprintf('%.0f %%',100*app.ed.areaFrac(k)), ...
        sprintf('%.2f',medianOf(app.ed.conf,k)), ...
        sprintf('%.2f',medianOf(app.ed.confMin,k)), ...
        words, app.keep(k)};
end
app.c.table.Data=D;
setApp(fig,app);
end

function refreshText(fig)
%refreshText  The title says how many are kept, the sentence says what is wrong with the
%   one that is selected.  Both are read off the model rather than off the table.
app=getApp(fig);
title(app.c.axTrace,sprintf('%d of %d repetitions kept%s', ...
    sum(app.keep),numel(app.keep),signalPhrase(app.ed)))
app.c.detail.Text=detailText(fig);
setApp(fig,app);
end

function p=signalPhrase(ed)
%signalPhrase  Which trace the responding area was measured on, when the caller said.
if isempty(ed.signal), p=''; return; end
p=sprintf('   (responding area of %s)',ed.signal);
end

function txt=detailText(fig)
%detailText  The selected repetition in one sentence: how much of the field responded,
%   how confident its median segment was, and WHY that segment scored as it did.  The
%   reason is the point - it is what turns "this repetition is poor" into something an
%   operator can act on.
app=getApp(fig);
k=app.sel;
if k<1 || k>numel(app.keep)
    txt='Select a repetition to see why it scored as it did.';
    return
end
[~,words]=weakestOf(app.ed,k);
txt=sprintf(['Repetition %d starts at %.0f s.  %.0f %% of the segmented field ' ...
    'responded, and its median segment scored %.2f overall.  Its weakest check, at ' ...
    '%.2f, is that %s.'], k,app.ed.epochStart(k),100*app.ed.areaFrac(k), ...
    medianOf(app.ed.conf,k),medianOf(app.ed.confMin,k),words);
end

function [name,words]=weakestOf(ed,k)
%weakestOf  The check the MEDIAN SEGMENT of one repetition scored lowest on, named and
%   explained.
%
%   THE MEDIAN OF EACH CHECK, then the smallest of those - not the check that was the
%   minimum for the most segments.  Those are different questions: a check can be the
%   weakest one for a third of the field without ever being the check that cancelled a
%   pair, which is exactly what the reference recording does with its slowest timing
%   level.  This column sits beside two medians and has to be comparable with them;
%   '_rep_nvcconfidence.jpg' answers the other question, as a share per repetition.
%
%   THE FACTOR SET IS READ, NEVER LISTED.  Nine by default, fourteen with the
%   across-repetition group on, and the timing family named from a setting.
name=''; words='every check passed';
best=Inf;
for j=1:1:numel(ed.factorNames)
    F=ed.factors.(ed.factorNames{j});
    if size(F,2)<k, continue; end
    v=median(double(F(:,k)),'omitnan');
    if isfinite(v) && v<best
        best=v; name=ed.factorNames{j}; words=ed.factorWords{j};
    end
end
if isempty(name), name='-'; words='there is nothing to judge this repetition on'; end
end

function v=medianOf(A,k)
%medianOf  One repetition's median across segments, as a plain double.
v=NaN;
if size(A,2)>=k, v=median(double(A(:,k)),'omitnan'); end
end

%% ===================== the edits ==================================== %%
function onTableEdit(fig,e)
%onTableEdit  THE KEEP BOX IS THE EDIT.  Every other column is a number this window was
%   handed and cannot change, so they are not editable and there is nothing else to
%   dispatch on.
if isempty(e.Indices) || e.Indices(2)~=keepColumn(), return; end
apiSet(fig,e.Indices(1),logical(e.NewData));
apiSelect(fig,e.Indices(1));
end

function onRowSelected(fig,e)
%onRowSelected  Selecting a row selects that repetition: it is the one the trace
%   highlights and the one the sentence under the table describes.
if isempty(e.Indices), return; end
apiSelect(fig,e.Indices(1));
end

function onAxesClicked(fig)
%onAxesClicked  Clicking the trace flips the repetition under the pointer, which is the
%   gesture the picker this window replaces had.  A click outside every span selects
%   nothing and changes nothing - a miss must not flip whichever repetition happened to
%   be nearest.
app=getApp(fig);
p=app.c.axTrace.CurrentPoint;
if isempty(p), return; end
x=p(1,1);
k=find(x>=app.ed.epochStart & x<=app.ed.epochEnd,1);
if isempty(k), return; end
apiToggle(fig,k);
apiSelect(fig,k);
end

%% ===================== finishing ==================================== %%
function requestDone(fig)
%requestDone  The one way out: set the sentinel the blocking waitfor watches.  The window
%   is NOT deleted here - the caller still has to read the decision off it.
if ~isvalid(fig), return; end
if ~confirmEmpty(fig), return; end
fig.UserData='done';
end

function ok=confirmEmpty(fig)
%confirmEmpty  Ask when nothing is left to average, and only then.  An empty set is legal
%   - every repetition is measured and stored either way - but every per-segment column
%   then comes out NaN, and with the representative repetition switched on the run stops
%   with nothing to average.  A window nobody can see is silent: an invisible editor is
%   being driven by a script, and a modal dialog on it would wait for an answer that never
%   comes.
ok=true;
app=getApp(fig);
if any(app.keep), return; end
if ~strcmpi(char(string(fig.Visible)),'on'), return; end
try
    sel=uiconfirm(fig,['No repetition is kept, so every per-segment result of this ' ...
        'recording will be empty.  Leave it that way?'], ...
        'Nothing is kept','Options',{'Leave it empty','Go back'}, ...
        'DefaultOption',2,'CancelOption',2,'Icon','warning');
catch
    return          % no dialog to be had: do not trap the operator in the window
end
ok=strcmp(sel,'Leave it empty');
end

function keep=collect(fig)
%collect  The decision as it stands, in the shape the caller assigns straight back.
keep=false(0,1);
if ~isvalid(fig), return; end
app=getApp(fig);
keep=logical(app.keep(:));
end

function releaseWindow(fig)
%releaseWindow  The window, released however this editor ended.  A parked handle onto a
%   figure that was never deleted is what keeps an onCleanup unfired somewhere else.
if ~isgraphics(fig), return; end
delete(fig);
end

%% ===================== the programmatic API ========================= %%
function api=buildAPI(fig)
%buildAPI  The same callbacks the widgets use, as a struct of handles, so the editor can
%   be driven without gestures (tests, and any future scripted use).  Every index is a
%   repetition number, which is also the table's row.
api=struct( ...
    'toggle',  @(k) apiToggle(fig,k), ...
    'set',     @(k,tf) apiSet(fig,k,tf), ...
    'setAll',  @(tf) apiSetAll(fig,tf), ...
    'restore', @() apiRestore(fig), ...
    'select',  @(k) apiSelect(fig,k), ...
    'selected',@() apiSelected(fig), ...
    'clickAt', @(x) apiClickAt(fig,x), ...
    'trust',   @() collect(fig), ...
    'state',   @() apiState(fig), ...
    'rows',    @() apiRows(fig), ...
    'summary', @() summaryText(fig), ...
    'detail',  @() detailText(fig), ...
    'drawnPoints',@() apiDrawnPoints(fig), ...
    'done',    @() requestDone(fig), ...
    'getApp',  @() getApp(fig));
end

function apiToggle(fig,k)
app=getApp(fig);
if k<1 || k>numel(app.keep), return; end
apiSet(fig,k,~app.keep(k));
end

function apiSet(fig,k,tf)
app=getApp(fig);
if k<1 || k>numel(app.keep), return; end
app.keep(k)=logical(tf);
setApp(fig,app);
refreshAll(fig);
end

function apiSetAll(fig,tf)
app=getApp(fig);
app.keep(:)=logical(tf);
setApp(fig,app);
refreshAll(fig);
end

function apiRestore(fig)
%apiRestore  The analysis's own proposal, back.  It is kept beside the operator's set for
%   exactly this - a window with no way back from a wrong click is one people stop using.
app=getApp(fig);
app.keep=app.ed.epochTrust;
setApp(fig,app);
refreshAll(fig);
end

function apiSelect(fig,k)
app=getApp(fig);
if k<1 || k>numel(app.keep), k=0; end
app.sel=k;
setApp(fig,app);
refreshSpans(fig);
refreshText(fig);
end

function k=apiSelected(fig)
app=getApp(fig); k=app.sel;
end

function apiClickAt(fig,x)
%apiClickAt  The trace click, at a time rather than at a pixel - the same dispatch the
%   mouse takes, without a pointer.
app=getApp(fig);
k=find(x>=app.ed.epochStart & x<=app.ed.epochEnd,1);
if isempty(k), return; end
apiToggle(fig,k);
apiSelect(fig,k);
end

function st=apiState(fig)
%apiState  WHAT THIS WINDOW IS SHOWING - the operator's set, the proposal it started
%   from, what is selected, and the words beside every repetition.  A nested read of the
%   live app and not a parked copy: a handle that captured the state by value would hand
%   every caller the set this window opened with.
app=getApp(fig);
n=numel(app.keep);
st=struct('trust',logical(app.keep(:)),'proposed',logical(app.ed.epochTrust(:)), ...
    'nEpochs',n,'selected',app.sel,'nKept',sum(app.keep), ...
    'summary',summaryText(fig),'detail',detailText(fig), ...
    'areaFrac',app.ed.areaFrac(:),'weakest',{cell(1,n)},'weakestWords',{cell(1,n)}, ...
    'conf',zeros(n,1),'confMin',zeros(n,1));
for k=1:1:n
    [st.weakest{k},st.weakestWords{k}]=weakestOf(app.ed,k);
    st.conf(k)   =medianOf(app.ed.conf,k);
    st.confMin(k)=medianOf(app.ed.confMin,k);
end
end

function D=apiRows(fig)
app=getApp(fig); D=app.c.table.Data;
end

function t=summaryText(fig)
app=getApp(fig);
t=sprintf('%d of %d repetitions kept',sum(app.keep),numel(app.keep));
end

function n=apiDrawnPoints(fig)
%apiDrawnPoints  How many points the trace layers actually put on screen.  There is no
%   decimation in this window, so this is the number that says whether one is needed.
app=getApp(fig); n=0;
kids=app.c.axTrace.Children;
for i=1:1:numel(kids)
    if isgraphics(kids(i)) && isprop(kids(i),'XData'), n=n+numel(kids(i).XData); end
end
end

%% ===================== small helpers ================================ %%
function app=getApp(fig), app=getappdata(fig,'app'); end
function setApp(fig,app), setappdata(fig,'app',app); end
