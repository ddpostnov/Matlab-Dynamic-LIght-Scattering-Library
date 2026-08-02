% guiExplore - Look at what a results file actually contains, and make a figure of it.
%
% WHAT THIS TOOL DOES
%   A single-window app for the *_r.mat RESULTS this library writes - the speckle
%   chain (runSegmentation -> runPulsatility / runVasomotion / runCTTH /
%   setVesselTypes) and the myograph one (runMyographDiameter ->
%   runMyographPropagation -> runMyographVasomotion) - on ONE recording or on many
%   at once, ending in a publication-ready figure at a chosen resolution.
%
%   IT READS THE FILES.  There is no list of plots in this file to fall out of date.
%   Every numeric leaf of the selected results is walked, classified by the axes it
%   lives on, and offered through a chain of menus that each narrow the next:
%
%       Kind -> Family -> Signal -> Variable -> Plot -> Where -> how to compare
%
%   Each step is built from the one above it out of what the files in hand actually
%   hold.  A step with nothing to offer is switched off rather than left showing the
%   previous file's answer; a step with exactly one answer fills itself and
%   disappears.  A quantity that is in the file can be drawn, one that is not is
%   never offered, and nobody has to edit this file when a results tree grows a new
%   quantity.
%
%   KIND IS WHAT YOU SEE ON THE AXES, not what the array is: one number each, curves,
%   a picture, or a picture with a dimension to scrub through.  One quantity can
%   appear under more than one - a per-segment number really is both a box plot
%   across groups and a map painted back into the image - and that is the question
%   the Kind menu is asking.
%
% THE SHAPE MODEL - WHY SOME PLOTS ARE OFFERED AND OTHERS ARE NOT
%   Everything measured belongs to some UNIT and varies along some AXES.  An axis is
%   a coordinate: it can go on the x of a plot.  A unit is a thing there are many of:
%   it has to be combined before anything can be compared.  Telling the two apart is
%   the whole of what decides which plots exist.
%
%     unit          one of them is             what a plot can do with it
%     -----------   ------------------------   ---------------------------------
%     segment       one vessel segment         COMBINE them - area-weighted, so a
%     vessel        one tracked vessel         large vessel counts for more than a
%                                              small one - or PAINT each one's value
%                                              back into the image through the
%                                              segment map
%     pixel         one pixel                  it is the picture
%     position      one line across a vessel   a COORDINATE: plot along the vessel,
%                                              or combine the lines
%     window        one analysed window        a COORDINATE: plot across the
%                                              protocol, or look at one window
%     recording     the recording itself       nothing to combine
%
%     axis          measured in                on a plot
%     -----------   ------------------------   ---------------------------------
%     time          seconds                    x
%     frequency     Hz                         x, logarithmic
%     percentile    0-100                      x, or the dimension to scrub
%     harmonic      cycle number               x, drawn as stems
%
%   SO THERE IS NO "VALUE AGAINST SEGMENT NUMBER" PLOT, and there never was a
%   meaningful one: a segment number is a name, not a position, and ordering by it
%   draws a shape that says nothing about the tissue.  Position along a vessel and
%   protocol order ARE positions, so "lag along the vessel" and "diameter across the
%   windows" are proper plots and are offered as such.
%
% ONE MODALITY AT A TIME
%   Speckle results and myograph results are never compared in one plot: they share
%   no unit, no axis and no y-scale, so a single figure holding both would be a
%   figure of two unrelated things.  The first file loaded decides which kind of
%   recording the window is about; anything else stays in the list, marked, excluded
%   and counted in the status line rather than quietly dropped.  "Clear" releases it.
%
%   ONE PIPELINE AT A TIME, for the same reason one step down.  A recording processed
%   both ways writes four results files ('_c_K', '_c_BFI', '_t_K', '_t_BFI') that mean
%   different things by the same variable names - the segment list is not even the
%   same length across them - so the Product menu shows one pipeline at a time.
%
%   MYOGRAPH RECORDINGS GO THROUGH THE SAME MACHINERY.  A myograph recording keeps
%   its results inside the WINDOWS it was analysed in, and one window's vasomotion
%   tree is the same shape as a speckle recording's, so the same cascade reaches
%   both.  What is different is the axis: the analysed WINDOW joins group, recording
%   index, animal and type as a fifth, independent way to organise a comparison, and
%   for a myograph protocol it is usually the one wanted on the x (baseline against
%   drug against washout).  The SIGNAL of a myograph comparison is the diameter
%   measure, or the recorded channel under the real name LabChart gave it.
%
%   IT IS ALSO THE ONLY PLACE A MYOGRAPH FIGURE COMES FROM.  No myograph step writes
%   a report page, by design, so the checks such a page would have carried live here:
%   the diameter along the vessel, which shows whether detection held; the detected
%   walls over a frame of the recording, which shows whether it found the right
%   edges; and the individual line traces, which show whether the lines moved
%   together - the assumption a propagation speed is fitted under.  Those three read
%   the recording beside the results rather than the results themselves, and are
%   offered as variables of their own under the Diameter family.
%
% WHEN TO USE IT
%   After processing recordings to *_r.mat results, to look at single-recording
%   detail or to compare experimental groups and sequences of recordings across many
%   files, and then export the figure for a paper or a talk.
%
%   It is STANDALONE.  It does not need - and never loads - the Processing Workbench.
%   The hand-off between them is the SESSION file, which carries the curated file
%   list together with each file's ANIMAL, recording TYPE and experimental GROUP; the
%   session's experimental group IS this tool's "group".  What is offered to plot
%   comes from the files themselves and not from the session's record of what ran:
%   contents beat records, and the contents are read anyway.  Its peer on the export
%   side is guiExport.
%
% HOW TO USE IT
%   1. Get files in, three ways:
%        - "Choose Folder / File" -> a FILE: everything is plotted for that one
%          recording; a FOLDER: searched recursively, then labelled by the three
%          regular-expression boxes:
%            (a) Include  : keep only files whose NAME matches this regexp
%                           (default '_r\.mat$' = results files).
%            (b) Group    : a regexp whose MATCH labels the experimental group of
%                           each file (e.g. 'a2|i2|k\d' for anaesthesia condition,
%                           or 'PSY\d+' for animal).  Files sharing a match are one
%                           group.
%            (c) Rec.index: a regexp whose MATCH is the recording index; files are
%                           stratified by it, on an axis of its own.  If the match
%                           contains a number, indices are ordered ASCENDING by it.
%          Click "Scan / apply" to (re)build the list.  Hover any of the three boxes
%          for worked pattern examples.
%        - "Load session..." -> a workbench session, whose recordings are resolved to
%          their results and arrive with every label already decided, so nothing is
%          re-derived from a file name.  The regexp boxes are switched off because
%          there is nothing left for them to guess.
%      Prefer to group by hand?  Select files in the list (Ctrl/Shift-click), type a
%      name in "New group" and press "Create group from selected files".  Hand-made
%      groups override BOTH the Group pattern and the session, and survive a re-scan.
%      If the folder holds more than one pipeline's results, pick which in "Product".
%   2. Choose WHAT to plot by working down section 2: Kind, Family, Signal, Variable,
%      Plot.  Signal and Variable take several at once as long as they are the same
%      shape - the three diameter measures together, or the two bands of one marker -
%      because two shapes on one axes would need two y-scales.  Picking something of
%      a different shape starts a fresh selection rather than refusing the click.
%   3. Choose WHERE in section 3: a vessel category, a named label, a region of the
%      map - offered from what the file actually holds, so a category it never had is
%      never listed.  "Units are" says how the segments or positions INSIDE one
%      recording are combined before recordings are compared; "X axis" and "Colour"
%      say how the recordings are then compared, over the FIVE independent axes -
%      group, recording index, animal, type, analysed window - which do not nest: an
%      animal may span groups and a group may span animals.  "Interval" narrows every
%      myograph plot to one analysed window, and is switched off when nothing loaded
%      has any.
%   4. If a slider appears under the picture, move it to look along the dimension the
%      picture is a slice of.  It is labelled in the recording's own coordinate
%      ("frequency = 0.12 Hz", never "frame 7"), where it sits is part of what an
%      export writes, and "Export video" writes the whole sweep to a file with the
%      colours held fixed so a recording that barely changes does not look alive.
%      There is no in-app playback.  A map AVERAGED over time and a single frame of
%      the same cube look identical, so the title always says which one it is.
%   5. Tweak titles / labels / legend, then "Plot".  Set DPI + format and "Export
%      image".
%
% PROGRAMMATIC USE
%   getappdata(fig,'exploreAPI') hands back the same logic as a struct of handles, so
%   the tool is fully headless-testable exactly like guiExport's API:
%
%     .loadPaths(paths, includePat, groupPat, recPat)   the folder / file route
%     .loadSession(sessionPath)                         the session route
%     .sessionPath()                                    where the file list came from
%     .index()                                          the union of what the selected
%                                                       files contain - the rows the
%                                                       cascade is built from
%     .setState('Key',value, ...)                       drive the cascade
%     .createGroup(name, fileIdx)                       a hand-made group
%     .render()                                         draw
%     .export(filename, dpi, format)                    write the figure
%     .exportVideo(filename, fps)                       write the whole sweep
%     .getApp()                                         the live application state
%
%   setState takes the cascade's own steps - 'Files' 'Product' 'Kind' 'Family'
%   'Signal' 'Variable' 'Plot' 'Frame' 'Where' 'Interval' 'Units' 'Baseline' 'XAxis'
%   'Colour' 'Stat' and the appearance keys - and applies them IN CASCADE ORDER
%   within one call, so one setState behaves like choosing them in the window.  An
%   entry that is not on offer raises an error rather than being ignored.
%
% NOTES
%   - The knowledge this window used to carry now lives in seven small modules beside
%     it, in GUIs/explore: exModality (which kind of recording), exAxes (the shared
%     coordinate vectors), exSchema (what each container holds), exScan (the walk
%     that produces the variable list), exPlotRules (which plots a quantity admits),
%     exFetch (the numbers behind one plot) and exMyograph (the three views that read
%     the recording rather than the results).  This file keeps the window, the
%     cascade, the renderers and the export.
%   - Only category-5 (lumen) segments are used for Artery/Vein/All-vessel selections
%     and, by default, for named vessel labels; parenchyma uses category 1.  See the
%     library README and getPixelCategories for the category definitions.
%   - Very large results files (per-pixel vasomotion, typically) are read leaf by leaf
%     rather than loaded, so time series, images and spectra can still be shown.  A
%     metrics table cannot be read that way, so for such a file every unit counts
%     equally and the status line says so.
%
% Syntax:
%    guiExplore                         % open the tool (folder / file route)
%    guiExplore(sessionPath)            % open it ON a workbench session
%    h = guiExplore('Visible','off')    % headless (programmatic drive / tests)
%
% Inputs:
%    sessionPath - optional, a workbench session file to open on.
%    'Visible'   - 'on' (default), or 'off' for a headless window.
%
% Outputs:
%    h - the uifigure handle, when one is requested.
%
% Example:
%    fig = guiExplore('Visible','off');
%    api = getappdata(fig,'exploreAPI');
%    api.loadPaths({'rec1_BFI_r.mat','rec2_BFI_r.mat'});
%    api.setState('Kind','Vector','Family','Vasomotion','Plot','curve.f');
%    api.export('spectrum.png', 300, 'PNG');
%
% See also: exScan, exPlotRules, exFetch, exAxes, exSchema, exModality, exMyograph,
%           wbSession, guiExport, guiWorkbench, exportToExcel, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------

function h = guiExplore(varargin)

% ---- inputs: an optional leading SESSION path, plus 'Visible' -------------
[sessionArg, vis] = parseArgs(varargin);

% ---- shared application state (seen by all nested functions) --------------
app = struct();
app.mode        = '';        % 'file' | 'folder' | 'session'
app.root        = '';        % chosen folder, or the folder of the chosen file
app.files       = emptyEntries();
app.cache       = containers.Map('KeyType','char','ValueType','any');
app.cacheOrder  = {};        % LRU order of cached paths
app.cacheLimit  = 6;         % max fully-loaded files kept in memory
app.sizeLimitGB = 3;         % above this a file is read field-by-field (HDF5)
app.groupOverride = containers.Map('KeyType','char','ValueType','char'); % hand-made groups (path->name)
app.manual      = struct('title',false,'xlab',false,'ylab',false); % user-edited?
app.sessionPath = '';        % the session this file list came from ('' = none)
% THE VARIABLE INDEX.  One exScan descriptor array per file, keyed by path and
% dropped with the file's entry in the LRU cache, plus the UNION over the files
% currently selected - which is what fills the menus, so a variable present in only
% some of them is still offered and the files that lack it simply contribute nothing.
app.index       = containers.Map('KeyType','char','ValueType','any');
app.vars        = emptyVars();
% ONE MODALITY AT A TIME (spec §6).  Locked by the first file a scan or a session
% produces; files of the other modality stay in the list, marked, and are excluded
% from every plot.  'Clear' or a new source releases it.
app.modality    = '';        % '' | 'speckle' | 'myograph'
app.excluded    = {};        % paths excluded by the modality lock
% THE BRANCH / PRODUCT FILTER (spec §5 step 1b).  One recording writes four results
% files - '_c_K', '_c_BFI', '_t_K', '_t_BFI' - and they disagree about what their own
% variable names mean, so they are filtered rather than pooled.
app.products    = {};        % the (branch, product) pairs present, in menu order
app.myoIntervals = {};       % the analysed windows, unioned over the file list
app.notes       = {};        % caveats the last fetch raised, for the status bar
app.hint        = '';        % a one-shot explanation of something the tool just did
% THE SCRUB POSITION (D3).  Which dimension the slider under the axes walks, its
% coordinate VALUES in the recording's own units, and where it currently sits.  It is
% part of what an export reproduces, so a written figure is the frame on screen.
app.scrub       = blankScrub();
app.want        = '';        % which of the candidate dims the user asked the slider to walk
app.wantKey     = '';        % ... and which variable and plot that answer was about
app.titleFrame  = [];        % during a video export: the frame the title should name

% ---- build the window (single instance, like guiExport) --------------------
delete(findall(groot,'Type','figure','Tag','guiExplore'));
fig = uifigure('Name','guiExplore - LSCI results explorer', ...
    'Color','w','Position',[80 60 1440 860],'Visible',vis,'Tag','guiExplore');
if strcmpi(vis,'on')
    try, fig.WindowState = 'maximized'; catch, end   %#ok<NOCOM> not every host allows it
end
app.fig = fig;

outer = uigridlayout(fig,[1 2],'ColumnWidth',{'2.4x','1x'}, ...
    'Padding',[6 6 6 6],'ColumnSpacing',8,'BackgroundColor','w');

% left: the single large plotting axes, with the scrub strip beneath it
axPanel = uipanel(outer,'BackgroundColor','w','BorderType','none');
axPanel.Layout.Column = 1;
axL = uigridlayout(axPanel,[2 1],'RowHeight',{'1x',0},'Padding',[2 2 2 2], ...
    'RowSpacing',4,'BackgroundColor','w');
ax = uiaxes(axL); ax.Color = 'w'; ax.Box = 'on';
title(ax,'Choose a folder or a file to begin'); app.ax = ax;
% THE SCRUB STRIP (D3).  A slider, and no timer and no play button: a timer against
% a shared axes races the render path for nothing a slider does not already give.  It
% is here rather than in the control column because it belongs to the PICTURE - and
% its row collapses to zero height whenever the chosen plot has nothing to walk.
app.axLayout = axL;
scrubBar = uigridlayout(axL,[1 3],'ColumnWidth',{220,130,'1x'},'Padding',[10 0 10 2], ...
    'ColumnSpacing',10,'BackgroundColor','w','Visible','off');
hScrubLbl = uilabel(scrubBar,'Text','','FontColor',[0.25 0.25 0.25],'WordWrap','on');
% A spectrum cube can be walked along its FREQUENCY or along its TIME, and those are
% two different questions - so when the leaf offers both, the choice is offered too.
hScrubDim = uidropdown(scrubBar,'Items',{},'Visible','off', ...
    'ValueChangedFcn',@(~,~)onScrubDim(), ...
    'Tooltip','Which dimension the slider walks.');
hScrub    = uislider(scrubBar,'Limits',[1 2],'Value',1,'MajorTicks',[], ...
    'MinorTicks',[],'ValueChangedFcn',@(src,~)onScrub(src));
hScrub.Tooltip = ['Move along the dimension this picture is a slice of.  The ' ...
    'position is part of what "Export image" writes, so an exported figure is the ' ...
    'frame on screen.'];
app.scrubBar = scrubBar;

% right: scrollable control column
ctrlPanel = uipanel(outer,'BackgroundColor','w','BorderType','none','Scrollable','on');
ctrlPanel.Layout.Column = 2;
C = uigridlayout(ctrlPanel,[1 1],'RowHeight',{'fit'},'Padding',[2 2 2 2], ...
    'BackgroundColor','w','Scrollable','on');   % the grid owns the scroll, not the panel
stack = uigridlayout(C,[8 1],'RowHeight',repmat({'fit'},1,8), ...
    'RowSpacing',8,'Padding',[0 0 0 0],'BackgroundColor','w');

c = struct();   % control handles
c.scrubLbl = hScrubLbl;   % the scrub strip, built with the axes it belongs to
c.scrubDim = hScrubDim;
c.scrub    = hScrub;

% --- section 1: data source ---
s1 = section(stack,'1 - Data source',14);
c.sourceBtn = uibutton(s1,'Text','Choose Folder / File','ButtonPushedFcn',@onChooseSource);
c.sourceBtn.Layout.Row = 1; c.sourceBtn.Layout.Column = 1;
c.sessionBtn = uibutton(s1,'Text','Load session...','ButtonPushedFcn',@onLoadSession, ...
    'BackgroundColor',[0.90 0.94 1.0], ...
    'Tooltip',['Read a Processing Workbench session: its recordings arrive as their ' ...
               'results files, with the animal, the recording type, the experimental ' ...
               'group and the recording index already decided.']);
c.sessionBtn.Layout.Row = 1; c.sessionBtn.Layout.Column = 2;
c.sourceLbl = uilabel(s1,'Text','(nothing selected)','FontColor',[0.35 0.35 0.35], ...
    'WordWrap','on'); c.sourceLbl.Layout.Row = 2; c.sourceLbl.Layout.Column = [1 2];
c.includeF = labelledField(s1,3,'Include',  '_r\.mat$');
c.groupF   = labelledField(s1,4,'Group',    '');
c.recF     = labelledField(s1,5,'Rec.index','');
c.includeF.Tooltip = tipInclude();          % hover helpers with regexp examples
c.groupF.Tooltip   = tipGroup();
c.recF.Tooltip     = tipRec();
c.scanBtn  = uibutton(s1,'Text','Scan / apply','ButtonPushedFcn',@onScan);
c.scanBtn.Layout.Row = 6; c.scanBtn.Layout.Column = 1;
c.clearBtn = uibutton(s1,'Text','Clear','ButtonPushedFcn',@onClear, ...
    'Tooltip',['Empty the file list and forget which kind of recording is loaded, ' ...
               'so a myograph set and a speckle set can be looked at one after the other.']);
c.clearBtn.Layout.Row = 6; c.clearBtn.Layout.Column = 2;
[c.productLbl, c.product] = labelledDropH(s1,7,'Product',{'(none)'},@(~,~)onProduct());
c.product.Tooltip = ['One recording is usually processed along more than one pipeline, ' ...
    'and each writes its own results file.  They measure different things under the ' ...
    'same names, so one is shown at a time rather than all of them together.'];
c.fileList = uilistbox(s1,'Items',{},'Multiselect','on', ...
    'Tooltip',['Ctrl/Shift-click to select several files, type a name below and press ' ...
               '"Create group" to group them by hand (overrides the Group pattern).'], ...
    'ValueChangedFcn',@(~,~)onFileSelection());
c.fileList.Layout.Row = [8 11]; c.fileList.Layout.Column = [1 2];
c.groupNameF = labelledField(s1,12,'New group','');
c.groupNameF.Tooltip = 'Name for a hand-made group (see the file list).';
c.createBtn = uibutton(s1,'Text','Create group from selected files', ...
    'ButtonPushedFcn',@onCreateGroup, ...
    'Tooltip','Assign the files selected above to the group named on the left.');
c.createBtn.Layout.Row = 13; c.createBtn.Layout.Column = [1 2];

% --- section 2: what to plot - the cascade (spec §5 steps 2-6) ---
% Each list is built from the one above it, out of what the selected files actually
% contain.  A step with nothing to offer is switched off rather than left showing the
% previous file's answer, and a step with exactly one answer fills itself and hides.
s2 = section(stack,'2 - What to plot',6,{'fit','fit',78,120,'fit','fit'});
[c.kindLbl,   c.kind]     = labelledDropH(s2,1,'Kind',    {},@(~,~)onKind());
[c.familyLbl, c.family]   = labelledDropH(s2,2,'Family',  {},@(~,~)onFamily());
[c.signalLbl, c.signal]   = labelledListH(s2,3,'Signal',  @(~,~)onSignal());
[c.variableLbl,c.variable]= labelledListH(s2,4,'Variable',@(~,~)onVariable());
[c.plotLbl,   c.plot]     = labelledDropH(s2,5,'Plot',    {},@(~,~)onPlot());
c.kind.Tooltip     = 'What you see on the axes: one number each, curves, a picture, or something to scrub through.';
c.signal.Tooltip   = ['Which signal the numbers belong to.  Several can be compared at ' ...
    'once as long as they are the same kind of thing; picking one that is not starts a fresh selection.'];
c.variable.Tooltip = ['Which quantity to draw.  Several can be compared at once as long ' ...
    'as they share a shape; picking one that does not starts a fresh selection.'];
c.count = uilabel(s2,'Text','','FontColor',[0.45 0.45 0.45],'FontSize',10,'WordWrap','on');
c.count.Layout.Row = 6; c.count.Layout.Column = [1 2];

% --- section 3: where, and how to compare (spec §5 steps 7-9) ---
s3 = section(stack,'3 - Where and how to compare',8,{110,'fit','fit','fit','fit','fit','fit','fit'});
c.selList = uilistbox(s3,'Items',{'All vessels (lumen)'},'Multiselect','on', ...
    'Tooltip',['Which part of the recording to read - offered from what this file ' ...
               'actually holds, so a category it never had is never listed.'], ...
    'ValueChangedFcn',@(~,~)onWhere());
c.selList.Layout.Row = 1; c.selList.Layout.Column = [1 2];
c.interval = labelledDrop(s3,2,'Interval',{'(all)'},@(~,~)requestRender());
c.interval.Tooltip = ['Which analysed window a myograph recording is plotted for.  ' ...
    'Switched off when nothing loaded has any.'];
c.interval.Enable = 'off';
c.units    = labelledDrop(s3,3,'Units are', ...
    {'Auto','Pooled (area-weighted)','One point each','One curve each'}, ...
    @(~,~)requestRender());
c.units.Tooltip = ['How the segments, vessels or positions INSIDE one recording are ' ...
    'combined before recordings are compared with each other.  Pooling weights by area, ' ...
    'so a large vessel counts for more than a small one.'];
c.baseline = uicheckbox(s3,'Text','As % of its own baseline','Value',false, ...
    'ValueChangedFcn',@(~,~)requestRender(), ...
    'Tooltip',['Divide every curve by its own mean, so vessels of different calibre ' ...
               'become comparable and a dilation is a number a reader recognises.']);
c.baseline.Layout.Row = 4; c.baseline.Layout.Column = [1 2];
c.xaxis    = labelledDrop(s3,5,'X axis',{'Auto'},@(~,~)requestRender());
c.xaxis.Tooltip = ['What separates the boxes, or the curves: the experimental group, the ' ...
    'recording index, the animal, the recording type, the analysed window, or what was ' ...
    'selected.  Only the ones that actually differ across what is loaded are offered.'];
c.colour   = labelledDrop(s3,6,'Colour',{'Auto'},@(~,~)requestRender());
c.colour.Tooltip = 'A second, independent way to split the same data - drawn as colour rather than as position.';
c.stat     = labelledDrop(s3,7,'Centre / error', ...
    {'Mean +/- SD','Mean +/- SEM','Median +/- IQR'},@(~,~)requestRender());

% --- section 4: appearance ---
s4 = section(stack,'4 - Appearance',8);
c.titleF = labelledField(s4,1,'Title','');
c.titleF.ValueChangedFcn = @(src,~)onManualEdit('title',src);
c.xlabF  = labelledField(s4,2,'X label','');
c.xlabF.ValueChangedFcn  = @(src,~)onManualEdit('xlab',src);
c.ylabF  = labelledField(s4,3,'Y label','');
c.ylabF.ValueChangedFcn  = @(src,~)onManualEdit('ylab',src);
c.legendChk = uicheckbox(s4,'Text','Show legend','Value',true, ...
    'ValueChangedFcn',@(~,~)requestRender());
c.legendChk.Layout.Row = 4; c.legendChk.Layout.Column = [1 2];
c.legendF = labelledField(s4,5,'Legend fmt','%s%g%r');
c.legendF.ValueChangedFcn = @(~,~)requestRender();
c.legHelp = uilabel(s4,'Text',['tokens: %s sel  %g group  %r index  %a animal  ' ...
    '%t type  %i interval  %f file  %v variable  %n signal'], ...
    'FontColor',[0.5 0.5 0.5],'FontSize',10,'WordWrap','on');
c.legHelp.Layout.Row = 6; c.legHelp.Layout.Column = [1 2];
c.cmap = labelledDrop(s4,7,'Image cmap',{'parula','turbo','hot','gray','jet'},@(~,~)requestRender());
lbFS = uilabel(s4,'Text','Font size'); lbFS.Layout.Row=8; lbFS.Layout.Column=1;
c.fontSize = uispinner(s4,'Value',12,'Limits',[6 40],'Step',1,'RoundFractionalValues','on', ...
    'ValueChangedFcn',@(~,~)requestRender(), ...
    'Tooltip','Base font size (points). Title and axis labels scale up from this; ticks, legend and colorbar match it.');
c.fontSize.Layout.Row=8; c.fontSize.Layout.Column=2;

% --- section 5: export ---
s5 = section(stack,'5 - Export',4);
c.dpi = labelledSpin(s5,1,'DPI',300);
c.fmt = labelledDrop(s5,2,'Format',{'PNG','TIFF','PDF (vector)','EPS (vector)','JPEG'},[]);
c.exportBtn = uibutton(s5,'Text','Export image','ButtonPushedFcn',@onExport);
c.exportBtn.Layout.Row = 3; c.exportBtn.Layout.Column = [1 2];
c.videoBtn = uibutton(s5,'Text','Export video','ButtonPushedFcn',@onExportVideo, ...
    'Enable','off', ...
    'Tooltip',['Write the whole sweep along the slider to a video file, one frame ' ...
               'per position, with the colours held fixed so a flat recording looks flat.']);
c.videoBtn.Layout.Row = 4; c.videoBtn.Layout.Column = [1 2];

% --- plot button + status ---
sPlot = uigridlayout(stack,[2 1],'RowHeight',{34,'fit'},'Padding',[0 0 0 0], ...
    'RowSpacing',4,'BackgroundColor','w');
c.plotBtn = uibutton(sPlot,'Text','Plot / refresh','FontWeight','bold', ...
    'BackgroundColor',[0.90 0.94 1.0],'ButtonPushedFcn',@(~,~)doRender());
c.status  = uilabel(sPlot,'Text','Ready.','FontColor',[0.30 0.30 0.30],'WordWrap','on');

app.c = c;
refreshCascade();                   % nothing loaded yet: every step switched off

% expose a small API so the tool can be driven programmatically (used for tests).
% 'DataType' is GONE rather than aliased: a silent alias would let a stale caller
% pass while plotting something else.
setappdata(fig,'exploreAPI',struct( ...
    'loadPaths',   @loadPathsProgrammatic, ...
    'loadSession', @(p) loadSessionFile(p), ...
    'sessionPath', @sessionPathLive, ...
    'index',       @indexLive, ...
    'setState',    @setStateProgrammatic, ...
    'createGroup', @createGroupProgrammatic, ...
    'render',      @doRender, ...
    'export',      @exportTo, ...
    'exportVideo', @exportVideoTo, ...
    'getApp',      @getAppLive));

% ---- a leading argument is a SESSION: the hand-off from the workbench ------
if ~isempty(sessionArg)
    try
        loadSessionFile(sessionArg);
    catch ME
        setStatus(['Session could not be read: ' ME.message]);
    end
end

if nargout>0, h = fig; end

%% ======================= CALLBACKS ====================================== %%
    function onChooseSource(~,~)
        choice = uiconfirm(fig,'Select a data source type:','Data source', ...
            'Options',{'Folder (many files)','Single file','Cancel'}, ...
            'DefaultOption',1,'CancelOption',3);
        switch choice
            case 'Folder (many files)'
                d = uigetdir(pwd,'Select the folder that contains your *_r.mat files');
                if isequal(d,0), return; end
                app.mode = 'folder'; app.root = d; app.sessionPath = '';
                c.sourceLbl.Text = ['Folder: ' d];
                setPatternFields('on');
                onScan();
            case 'Single file'
                [f,p] = uigetfile({'*_r.mat;*.mat','LSCI results (*_r.mat)'}, ...
                    'Select a results file');
                if isequal(f,0), return; end
                releaseModality();
                app.mode = 'file'; app.root = p; app.sessionPath = '';
                app.files = makeEntry(fullfile(p,f),'all','1');
                c.sourceLbl.Text = ['File: ' f];
                setPatternFields('off');
                applyFilters(); refreshFileList(); refreshIndex(); doRender();
        end
    end

    function onClear(~,~)
        %onClear  Empty the list and RELEASE THE MODALITY LOCK, which is the only way
        %   to look at a myograph set after a speckle one without restarting the tool.
        releaseModality();
        app.mode=''; app.root=''; app.sessionPath=''; app.files = emptyEntries();
        app.index = containers.Map('KeyType','char','ValueType','any');
        app.cache = containers.Map('KeyType','char','ValueType','any');
        app.cacheOrder = {};
        c.sourceLbl.Text = '(nothing selected)';
        setPatternFields('on');
        refreshFileList(); refreshIndex();
        cla(app.ax,'reset'); app.ax.Color='w'; title(app.ax,'Choose a folder or a file to begin');
        setStatus('Cleared.  Nothing is loaded, and either kind of recording can be opened next.');
    end

    function releaseModality()
        app.modality = ''; app.excluded = {}; app.products = {};
    end

    function onLoadSession(~,~)
        [f,p] = uigetfile({'*.mat','Workbench session (*.mat)'}, ...
            'Select a Processing Workbench session', app.root);
        if isequal(f,0), return; end
        try
            loadSessionFile(fullfile(p,f));
        catch ME
            uialert(fig,ME.message,'Session');
        end
    end

    function loadSessionFile(pth)
        %loadSessionFile  THE hand-off from the workbench (spec §5).  Read through
        %   wbSession('read',...) - never by reaching into the file layout - so the
        %   curated list arrives with its animal / type / index / experimental group
        %   already resolved, and with the record of what actually ran.  Each session
        %   recording is expanded to the *_r.mat RESULTS products it owns, which is
        %   what this tool plots; hand-made groups still win over the session's.
        pth = char(pth);
        S = wbSession('read', pth);
        entries = sessionEntries(S);
        releaseModality();
        app.mode = 'session'; app.sessionPath = pth;
        app.root = fileparts(pth);
        app.files = applyGroupOverrides(entries);
        c.sourceLbl.Text = sprintf('Session: %s   (%d recording(s) -> %d results file(s))', ...
            shortName(pth), numel(S.files), numel(app.files));
        setPatternFields('off');
        c.scanBtn.Enable = 'on';                 % re-reads the session (see onScan)
        c.scanBtn.Tooltip = 'Re-read the session file and rebuild the list.';
        applyFilters(); refreshFileList(); refreshIndex();
        if isempty(app.files)
            setStatus('The session has no *_r.mat results yet - process it first.');
            return
        end
        setStatus(sprintf('Session: %d results file(s), %d group(s), %d animal(s), %d type(s).%s', ...
            numel(app.files), nUnique('group'), nUnique('animal'), nUnique('type'), ...
            filterNote()));
        doRender();
    end

    function e = sessionEntries(S)
        %sessionEntries  Session -> file-list entries.  Every label is the SESSION's
        %   (expGroup IS this tool's group - spec §2/D1), never re-derived from the
        %   name, and every entry remembers which steps that recording completed so
        %   the plot menus can stop probing for trees that were never computed.
        e = emptyEntries();
        for i = 1:numel(S.files)
            f = S.files(i);
            rp = resolveResultFiles(f.path);
            for k = 1:numel(rp)
                x = makeEntryStruct(rp{k}, '', '');
                x.group  = labelOr(f.expGroup,'all');
                x.rec    = labelOr(f.index,'1');
                x.recnum = recToNum(x.rec);
                x.animal = charOf(f.animal);
                x.type   = charOf(f.type);
                e(end+1) = x; %#ok<AGROW>
            end
        end
        e = sortEntries(e);
    end

    function onScan(~,~)
        if strcmp(app.mode,'session')
            loadSessionFile(app.sessionPath);   % "apply" for a session = re-read it
            return
        end
        if ~strcmp(app.mode,'folder'), return; end
        setStatus('Scanning folder...');
        try
            entries = buildFileEntries(app.root, c.includeF.Value, c.groupF.Value, c.recF.Value);
        catch ME
            uialert(fig,ME.message,'Scan error'); setStatus('Scan failed.'); return;
        end
        if isempty(entries)
            app.files = entries; refreshFileList(); refreshIndex();
            setStatus('No files matched the Include pattern.'); return;
        end
        releaseModality();
        app.files = applyGroupOverrides(entries);   % keep any hand-made groups
        applyFilters(); refreshFileList(); refreshIndex();
        setStatus(sprintf('Found %d files, %d group(s), %d index level(s).%s', ...
            numel(app.files), numel(unique({app.files.group})), ...
            numel(unique({app.files.rec})), filterNote()));
        doRender();
    end

    function onProduct(~,~)
        %onProduct  The branch / product filter (spec §5 step 1b).  It changes WHICH
        %   FILES are in play, so the whole index is rebuilt behind it.
        applyFilters(); refreshFileList(); refreshIndex(); doRender();
    end

    function onFileSelection(~,~)
        %onFileSelection  A different set of files is a different UNION, so the menus
        %   are rebuilt before the plot is.
        refreshIndex(); requestRender();
    end

    function onCreateGroup(~,~)
        % Assign the files currently selected in the list to a hand-typed group.
        if ~any(strcmp(app.mode,{'folder','session'}))
            uialert(fig,'Manual groups apply when a folder or a session is loaded.','Create group'); return;
        end
        name = strtrim(c.groupNameF.Value);
        if isempty(name), uialert(fig,'Type a group name first.','Create group'); return; end
        idx = c.fileList.Value;
        if isempty(idx)
            uialert(fig,'Select one or more files in the list above first.','Create group'); return;
        end
        for i = idx(:)'
            app.files(i).group = name;
            app.groupOverride(app.files(i).path) = name;   % remember across re-scans
        end
        refreshFileList();                                  % show new labels, select all -> plot all
        setStatus(sprintf('Group "%s" now has %d file(s). Repeat to define more groups.', ...
            name, numel(idx)));
        requestRender();
    end

    function e = applyGroupOverrides(e)
        % Re-apply any hand-made group assignments to a freshly-scanned file list.
        for i = 1:numel(e)
            if app.groupOverride.isKey(e(i).path), e(i).group = app.groupOverride(e(i).path); end
        end
    end

    function onKind(~,~),   refreshFamilyItems();   onFamily();   end
    function onFamily(~,~), refreshSignalItems();   onSignal();   end
    function onSignal(~,~)
        enforceSignature(c.signal);      % D2: one shape per axes
        refreshVariableItems(); onVariable();
    end
    function onVariable(~,~)
        enforceSignature(c.variable);    % D2, again - the two lists share the rule
        refreshPlotItems(); onPlot();
    end
    function onPlot(~,~)
        refreshWhereItems();
        refreshUnitsDefault();
        refreshCompareItems();
        refreshScrub();
        autofillCosmetics();
        requestRender();
    end
    function onScrub(src)
        %onScrub  A new position along the scrubbed dimension.  The index is what the
        %   fetch is asked for, so the picture that comes back IS that frame - the
        %   renderer never picks one of its own.
        if isempty(app.scrub.name), return; end
        i = clampScrub(round(src.Value));
        src.Value = i;
        app.scrub.idx = i;
        c.scrubLbl.Text = scrubText(app.scrub, i);
        autofillCosmetics();          % the coordinate is part of the title
        requestRender();
    end
    function onWhere(~,~)
        refreshCompareItems();
        requestRender();
    end

    function onManualEdit(which,src)
        app.manual.(which) = ~isempty(strtrim(src.Value));
        requestRender();
    end

    function onExport(~,~)
        if isempty(app.files), uialert(fig,'Nothing to export yet.','Export'); return; end
        exts = {'*.png','*.tif','*.pdf','*.eps','*.jpg'};
        k = find(strcmp(c.fmt.Value,{'PNG','TIFF','PDF (vector)','EPS (vector)','JPEG'}),1);
        [f,p] = uiputfile(exts{k},'Export figure as');
        if isequal(f,0), return; end
        try
            exportTo(fullfile(p,f), c.dpi.Value, c.fmt.Value);
            setStatus(['Exported: ' fullfile(p,f)]);
        catch ME
            uialert(fig,ME.message,'Export error');
        end
    end

    function onExportVideo(~,~)
        if ~strcmp(currentPlot(),'video')
            uialert(fig,'Choose the video plot first, then export it.','Export video');
            return
        end
        [f,p] = uiputfile({'*.avi','AVI (Motion JPEG)';'*.mp4','MP4'}, ...
            'Export video as','sweep.avi');
        if isequal(f,0), return; end
        try
            exportVideoTo(fullfile(p,f));
        catch ME
            uialert(fig,ME.message,'Export error');
        end
    end

%% ==================== FILE LIST / MODEL ================================= %%
    function refreshFileList()
        items = cell(1,numel(app.files));
        for i=1:numel(app.files)
            % the bracket carries whichever axes are actually labelled: group and
            % index always, animal and type only when a session supplied them
            bits = {app.files(i).group, app.files(i).rec, ...
                    app.files(i).animal, app.files(i).type};
            bits = bits(~cellfun(@isempty,bits));
            items{i} = sprintf('%s   [%s]', app.files(i).name, strjoin(bits,' | '));
            % An off-modality file STAYS IN THE LIST and is marked: a missing file has
            % to be visible, so it is never silently dropped (spec §6).
            if isExcluded(app.files(i).path)
                items{i} = ['- not shown -   ' items{i}];
            end
        end
        c.fileList.Items = items;
        c.fileList.ItemsData = 1:numel(app.files);
        if isempty(items)
            % An empty list takes the EMPTY CELL and nothing else - a 1x0 double is
            % rejected outright, which is how clearing the tool used to throw.
            c.fileList.Value = {};
            return
        end
        keep = find(~arrayfun(@(f) isExcluded(f.path), app.files));
        if isempty(keep), c.fileList.Value = {}; return; end
        c.fileList.Value = reshape(keep,1,[]);     % everything usable, selected
    end

    function idx = selectedFileIdx()
        idx = c.fileList.Value;
        if isempty(idx), idx = 1:numel(app.files); end
        idx = idx(:)';
        idx = idx(~arrayfun(@(i) isExcluded(app.files(i).path), idx));
    end

    function tf = isExcluded(pth)
        tf = any(strcmp(app.excluded, char(pth)));
    end

%% ==================== THE MODALITY LOCK AND THE PRODUCT FILTER ========== %%
    function applyFilters()
        %applyFilters  The two things that decide which of the listed files are in
        %   play, in the order they have to happen: the PRODUCT filter narrows the
        %   list to one pipeline, and the MODALITY lock then takes the first survivor's
        %   kind of recording and excludes anything else.
        refreshProductItems();
        app.excluded = {};
        keep = productKeep();
        for i = 1:numel(app.files)
            if ~keep(i), app.excluded{end+1} = app.files(i).path; end
        end
        lockModality();
    end

    function keep = productKeep()
        %productKeep  Which files match the chosen (branch, product) pair.  With
        %   nothing chosen - or nothing parseable - every file stays.
        keep = true(1,numel(app.files));
        want = c.product.Value;
        if isempty(app.products) || isempty(want) || ~ischar(want), return; end
        for i = 1:numel(app.files)
            keep(i) = strcmp(productKeyOf(app.files(i).path), want);
        end
        if ~any(keep), keep = true(1,numel(app.files)); end
    end

    function refreshProductItems()
        %refreshProductItems  The (branch, product) pairs the listed files hold.  The
        %   parse is wbProducts' and wbFileModel's - never a prefix test on the name,
        %   for the reasons wbProducts' header sets out at length.
        n = numel(app.files);
        keys = cell(1,n); labels = cell(1,n); count = zeros(1,n); nk = 0;
        for i = 1:n
            k = productKeyOf(app.files(i).path);
            j = find(strcmp(keys(1:nk),k),1);
            if isempty(j)
                nk = nk+1; keys{nk} = k; labels{nk} = productLabelOf(k); count(nk) = 1;
            else
                count(j) = count(j)+1;
            end
        end
        keys = keys(1:nk); labels = labels(1:nk); count = count(1:nk);
        [keys,labels] = orderProducts(keys,labels,count);
        app.products = keys;
        keep = c.product.Value;
        if isempty(keys)
            c.product.Items = {'(none)'}; c.product.ItemsData = {''};
            c.product.Value = ''; c.product.Enable = 'off';
            setVisible(c.productLbl, c.product, false);
            return
        end
        c.product.Items = labels; c.product.ItemsData = keys;
        if ischar(keep) && any(strcmp(keys,keep)), c.product.Value = keep;
        else, c.product.Value = keys{1}; end
        c.product.Enable = 'on';
        setVisible(c.productLbl, c.product, numel(keys)>1);
    end

    function [keys,labels] = orderProducts(keys,labels,count)
        %orderProducts  The default is the MOST-PROCESSED product of the MOST-POPULATED
        %   branch - which on a folder holding a whole recording's output is the one a
        %   user came to look at, rather than the intermediate the pipeline passed
        %   through.  Ties on population fall back to the pipeline order.
        n = numel(keys); rank = zeros(n,3);
        for i = 1:n
            p = strsplit(keys{i},'|');
            br = p{1}; pr = ''; if numel(p)>1, pr = p{2}; end
            rank(i,:) = [-count(i), branchRank(br), -productRank(pr)];
        end
        [~,ord] = sortrows(rank,[1 2 3]);
        keys = keys(ord); labels = labels(ord);
    end

    function lockModality()
        %lockModality  ONE MODALITY AT A TIME (D1).  The first file still standing
        %   after the product filter decides; anything of the other kind is excluded
        %   and SAID SO, because a silently missing file is the failure this prevents.
        idx = find(~arrayfun(@(f) isExcluded(f.path), app.files));
        if isempty(idx), return; end
        if isempty(app.modality)
            app.modality = modalityOfFile(app.files(idx(1)).path);
        end
        for i = idx(:)'
            if ~strcmp(modalityOfFile(app.files(i).path), app.modality)
                app.excluded{end+1} = app.files(i).path;
            end
        end
    end

    function m = modalityOfFile(pth)
        R = tryLoad(pth);
        if isempty(R), m = app.modality; if isempty(m), m = 'speckle'; end, return; end
        m = exModality(R);
    end

    function s = filterNote()
        %filterNote  What the two filters took out, named and counted.
        s = '';
        nOff = 0;
        for i = 1:numel(app.files)
            if ~isExcluded(app.files(i).path), continue; end
            if strcmp(productKeyOf(app.files(i).path), c.product.Value), nOff = nOff+1; end
        end
        nProd = numel(app.excluded) - nOff;
        if nProd>0
            s = sprintf('  %d file(s) hidden: they belong to another pipeline.', nProd);
        end
        if nOff>0
            other = 'speckle'; if strcmp(app.modality,'speckle'), other = 'myograph'; end
            s = sprintf('%s  %d file(s) excluded: %s results cannot be compared with %s results.', ...
                s, nOff, other, app.modality);
        end
    end

%% ==================== THE VARIABLE INDEX ================================ %%
    function refreshIndex()
        %refreshIndex  Read what the selected files CONTAIN, and let that fill the
        %   menus.  This is the whole redesign in one function: exScan walks each
        %   file once, the descriptors are cached beside the file, and the UNION over
        %   the selected files becomes the cascade - so a variable present in only
        %   some of them is still offered and the others simply contribute nothing.
        idx = selectedFileIdx();
        rows = emptyVars();
        IV = {};
        for i = idx
            pth = app.files(i).path;
            [D,S] = indexOf(pth);
            R = tryLoad(pth);
            if strcmp(exModality(R),'myograph')
                nm = exMyograph('names',R); IV = [IV, nm(~ismember(nm,IV))]; %#ok<AGROW>
            end
            rows = addRows(rows, D, S);
        end
        app.vars = sortVars(rows);
        app.myoIntervals = IV;
        refreshIntervalItems();
        refreshCascade();
    end

    function [D,S] = indexOf(pth)
        %indexOf  One file's descriptors, scanned once and kept with the file.  The
        %   cache is dropped with the loaded struct, so a file that leaves memory
        %   leaves the index with it and cannot go stale against a re-read.
        if app.index.isKey(pth)
            v = app.index(pth); D = v.D; S = v.S; return
        end
        R = tryLoad(pth);
        D = []; S = emptyVars();
        if ~isempty(R)
            try
                D = exScan(R, pth);
            catch ME
                setStatus(sprintf('%s could not be read: %s', shortName(pth), ME.message));
                D = [];
            end
            S = sourceViews(R, pth);
        end
        app.index(pth) = struct('D',D,'S',S);
    end

    function rows = addRows(rows, D, S)
        %addRows  Fold one file's descriptors into the union.  A row is one
        %   (family, signal, variable) - the concrete leaf paths behind it are
        %   resolved again at draw time, per file and per analysed window.
        for k = 1:numel(D)
            d = D(k);
            if ~isempty(d.pairedWith), continue; end   % a confidence interval is an
                                                       % error bar, not a variable
            P = exPlotRules(d);
            if isempty(P), continue; end               % no legal plot: nothing to offer
            r = mkVarRow(d, P);
            rows = mergeRow(rows, r);
        end
        for k = 1:numel(S)
            rows = mergeRow(rows, S(k));
        end
    end

    function rows = mergeRow(rows, r)
        j = [];
        if ~isempty(rows)
            j = find(strcmp({rows.key}, r.key), 1);
        end
        if isempty(j)
            rows(end+1) = r;
        else
            rows(j).nFiles = rows(j).nFiles + 1;
            if ~isequal(rows(j).varLabel, r.varLabel)
                % Two signals spell the same quantity differently - outer, mid and
                % inner diameter all say '<measure> diameter' - so the shared name is
                % the one that does not claim to be one of them.
                rows(j).varLabel = genericLabel(rows(j));
            end
        end
    end

    function r = mkVarRow(d, P)
        r = emptyVars();
        r(1).family      = d.family;
        r(1).signal      = d.signal;
        r(1).signalLabel = d.signalLabel;
        r(1).varId       = varIdOf(d);
        r(1).varLabel    = d.label; if isempty(r(1).varLabel), r(1).varLabel = d.leaf; end
        r(1).leaf        = d.leaf;
        r(1).unit        = d.unit;
        r(1).dims        = d.dims;
        r(1).sig         = signatureOf(d.unit, d.dims);
        r(1).kinds       = exPlotRules('kinds', d);
        r(1).plots       = P;
        r(1).suspect     = d.suspect;
        r(1).srcView     = '';
        r(1).nFiles      = 1;
        r(1).key         = [d.family '||' d.signal '||' r(1).varId];
    end

    function id = varIdOf(d)
        %varIdOf  What makes two leaves THE SAME VARIABLE across signals and across
        %   analysed windows: everything except which signal it belongs to and which
        %   window it was measured in.  outer / mid / inner diameter share one id, so
        %   the Variable step holds one entry and the Signal step holds three.
        id = strjoin({d.family, d.branch, d.leaf, d.unit, strjoin(d.dims,',')}, '|');
    end

    function s = signatureOf(unit, dims)
        %signatureOf  THE D2 RULE, in one string.  Entries sharing it can be drawn on
        %   one axes with one y-scale; entries that do not, cannot.
        s = [char(unit) '|' strjoin(dims,',')];
    end

    function s = genericLabel(r)
        s = r.leaf;
        if isempty(s), s = r.family; end
    end

%% ==================== THE CASCADE ======================================= %%
    function refreshCascade()
        refreshKindItems();
        refreshFamilyItems();
        refreshSignalItems();
        refreshVariableItems();
        refreshPlotItems();
        refreshWhereItems();
        refreshUnitsDefault();
        refreshCompareItems();
        refreshScrub();
        autofillCosmetics();
    end

    function refreshKindItems()
        %refreshKindItems  WHAT YOU SEE ON THE AXES - one number each, curves, a
        %   picture, or something to scrub through.  Offered in that fixed order, and
        %   only where the loaded files hold something that can be drawn that way.
        order = {'Scalar','Vector','Image','Video'};
        items = order(cellfun(@(k) any(cellfun(@(K) any(strcmp(K,k)), {app.vars.kinds})), order));
        setStep(c.kindLbl, c.kind, items, items);
    end

    function refreshFamilyItems()
        rows = rowsForKind();
        items = uniqueStable({rows.family});
        setStep(c.familyLbl, c.family, items, items);
    end

    function refreshSignalItems()
        rows = rowsForFamily();
        [keys,labels] = signalItems(rows);
        setStepList(c.signalLbl, c.signal, labels, keys);
    end

    function refreshVariableItems()
        %refreshVariableItems  The quantities left once the kind, the family and the
        %   signal have been chosen.  A paired descriptor never reaches here - it was
        %   filtered at index time, because a confidence interval has no plot of its
        %   own and offering it would be a dead end.
        rows = rowsForSignal();
        [keys,labels] = variableItems(rows);
        setStepList(c.variableLbl, c.variable, labels, keys);
    end

    function refreshPlotItems()
        rows = rowsForVariable();
        P = {};
        if ~isempty(rows)
            P = rows(1).plots;
            for i = 2:numel(rows), P = P(ismember(P, rows(i).plots)); end
            P = P(strcmp(cellfun(@kindOfPlot, P,'uni',0), currentKind()));
        end
        setStep(c.plotLbl, c.plot, cellfun(@plotLabel,P,'uni',0), P);
    end

    function refreshWhereItems()
        %refreshWhereItems  THE SELECTION LIST, BUILT FROM THE CONTENTS.  exFetch
        %   knows what one file can be asked for; this unions that over the selected
        %   files so a label only one of them carries is still offered.
        items = {};
        rows = rowsForVariable();
        if ~isempty(rows)
            for i = selectedFileIdx()
                d = firstDescriptor(app.files(i).path, rows(1));
                if isempty(d), continue; end
                R = tryLoad(app.files(i).path);
                try
                    w = exFetch('where', R, d);
                catch
                    w = {};
                end
                items = [items, w(~ismember(w,items))]; %#ok<AGROW>
            end
        end
        if isempty(items), items = {'(everything)'}; end
        keep = c.selList.Value; if ischar(keep), keep = {keep}; end
        c.selList.Items = items;
        keep = keep(ismember(keep,items));
        if isempty(keep), keep = items(1); end
        c.selList.Value = keep;
    end

    function refreshUnitsDefault()
        %refreshUnitsDefault  A box wants the distribution, a bar wants the pooled
        %   value, a family of curves indexed by position wants one curve each.  The
        %   default therefore comes from the PLOT and never from a constant.
        rows = rowsForVariable();
        if isempty(rows) || isempty(currentPlot()), c.units.Enable='off'; return; end
        c.units.Enable = onOff(~strcmp(rows(1).unit,'whole'));
    end

% ---- the scrub strip (D3) --------------------------------------------------
    function refreshScrub()
        %refreshScrub  THE SLIDER IS SHOWN ONLY WHEN THERE IS SOMETHING TO WALK, and
        %   what it walks it takes from the rules and the recording's own axis
        %   registry - never from 1:size(...,3), which would label a frequency sweep
        %   in frame numbers and call it an axis.
        % A CHOICE BELONGS TO THE VARIABLE IT WAS MADE ABOUT.  Carrying it to the
        % next one would quietly override the documented default - time first, because
        % a wave passing over the tissue is what anyone wants to see first - on a leaf
        % the user has said nothing about.
        key = [strjoin(selectedKeys(c.variable),',') '||' currentPlot()];
        if ~strcmp(key, app.wantKey), app.want = ''; app.wantKey = key; end
        prev = app.scrub;
        app.scrub = scrubSpec(app.want);
        if strcmp(prev.name, app.scrub.name)
            app.scrub.idx = clampScrub(prev.idx);   % the same dim: keep the position
        end
        tf = ~isempty(app.scrub.name) && app.scrub.n > 1;
        app.axLayout.RowHeight = {'1x', scrubRowHeight(tf)};
        app.scrubBar.Visible = onOff(tf);
        c.videoBtn.Enable = onOff(strcmp(currentPlot(),'video') && tf);
        refreshScrubChoices();
        if ~tf, c.scrubLbl.Text = ''; return; end
        c.scrub.Limits = [1 app.scrub.n];
        c.scrub.Value  = app.scrub.idx;
        c.scrubLbl.Text = scrubText(app.scrub, app.scrub.idx);
    end

    function refreshScrubChoices()
        %refreshScrubChoices  A CUBE WITH TWO CANDIDATE DIMENSIONS IS TWO QUESTIONS.
        %   Sweeping a spectrum along frequency shows which rhythms live where;
        %   sweeping it along time shows a wave passing over the tissue.  Time is the
        %   default, and the other is offered rather than hidden.
        ch = app.scrub.choices;
        show = numel(ch)>1;
        c.scrubDim.Visible = onOff(show);
        if ~show, c.scrubDim.Items = {}; return; end
        c.scrubDim.Items = cellfun(@(L) ['along ' coordText(L,[])], ...
            app.scrub.choiceLabels, 'uni',0);
        c.scrubDim.ItemsData = ch;
        c.scrubDim.Value = app.scrub.name;
    end

    function onScrubDim()
        app.want = c.scrubDim.Value;
        refreshScrub();
        autofillCosmetics();
        requestRender();
    end

    function s = scrubSpec(want)
        %scrubSpec  Which dimension the current plot leaves for the slider, asked of
        %   exFetch rather than worked out here: the length of a dim is a fact about
        %   the leaf and the axis registry together, and a window that guessed it
        %   would be guessing about a lazy file it never read.
        s = blankScrub();
        P = currentPlot(); rows = rowsForVariable(); idx = selectedFileIdx();
        if isempty(P) || isempty(rows) || isempty(idx) || ~isempty(rows(1).srcView), return; end
        pth = app.files(idx(1)).path;
        d = firstDescriptor(pth, rows(1));
        if isempty(d), return; end
        R = tryLoad(pth);
        if isempty(R), return; end
        try
            q = exFetch('scrub', R, d, P, struct('path',pth,'scrub',want));
        catch
            return
        end
        s.name = q.name; s.label = q.label; s.values = q.values;
        s.n = q.n; s.choices = q.choices; s.choiceLabels = q.choiceLabels; s.idx = 1;
    end

    function i = clampScrub(i)
        i = round(double(i));
        if ~isfinite(i) || i<1, i = 1; end
        if app.scrub.n>=1 && i>app.scrub.n, i = app.scrub.n; end
    end

    function sl = sliceNow()
        %sliceNow  The slider's position, in the words exFetch's opts.slice takes.
        sl = struct();
        if isempty(app.scrub.name), return; end
        sl.(exFetch('slicekey', app.scrub.name)) = app.scrub.idx;
    end

    function refreshIntervalItems()
        %refreshIntervalItems  The interval filter, offered only when something
        %   loaded actually has windows - the axis is real or it is not there.
        items = [{'(all)'} app.myoIntervals];
        keep = c.interval.Value;
        c.interval.Items = items;
        if ismember(keep,items), c.interval.Value = keep; else, c.interval.Value = items{1}; end
        c.interval.Enable = onOff(~isempty(app.myoIntervals));
    end

    function refreshCompareItems()
        %refreshCompareItems  ONLY THE AXES THAT VARY.  The five stratification axes
        %   are independent - an animal may span groups and a group may span animals -
        %   so what is offered is simply whichever of them differs across what is
        %   loaded, plus the selection when more than one was picked.
        avail = {'Auto'};
        idx = selectedFileIdx();
        if numel(unique({app.files(idx).group}))>1,  avail{end+1}='Group';           end
        if numel(unique({app.files(idx).rec}))>1,    avail{end+1}='Recording index'; end
        if numel(unique({app.files(idx).animal}))>1, avail{end+1}='Animal';          end
        if numel(unique({app.files(idx).type}))>1,   avail{end+1}='Type';            end
        if numel(app.myoIntervals)>1 && strcmp(c.interval.Value,'(all)'), avail{end+1}='Interval'; end
        sels = c.selList.Value; if ischar(sels), sels={sels}; end
        if numel(sels)>1, avail{end+1}='Selection'; end
        if numel(selectedKeys(c.variable))>1, avail{end+1}='Variable'; end
        if numel(selectedKeys(c.signal))>1,   avail{end+1}='Signal';   end
        avail{end+1} = 'Pool all';
        setDropKeeping(c.xaxis,  avail);
        setDropKeeping(c.colour, [{'(none)'} avail(~strcmp(avail,'Pool all'))]);
    end

% ---- the rows the cascade is looking at, one step at a time ----------------
    function rows = rowsForKind()
        rows = app.vars;
        k = currentKind();
        if isempty(rows) || isempty(k), rows = emptyVars(); return; end
        rows = rows(cellfun(@(K) any(strcmp(K,k)), {rows.kinds}));
    end
    function rows = rowsForFamily()
        rows = rowsForKind();
        if isempty(rows) || isempty(c.family.Value), rows = emptyVars(); return; end
        rows = rows(strcmp({rows.family}, c.family.Value));
    end
    function rows = rowsForSignal()
        rows = rowsForFamily();
        sg = selectedKeys(c.signal);
        if isempty(rows) || isempty(sg), rows = emptyVars(); return; end
        rows = rows(ismember({rows.signal}, sg));
    end
    function rows = rowsForVariable()
        rows = rowsForSignal();
        vs = selectedKeys(c.variable);
        if isempty(rows) || isempty(vs), rows = emptyVars(); return; end
        rows = rows(ismember({rows.varId}, vs));
    end

% ---- the two multi-select lists, and the rule that keeps them comparable ----
    function enforceSignature(ctl)
        %enforceSignature  D2.  A selection spanning two shapes would put two y-scales
        %   on one axes, so the newcomer wins and the rest of the selection is CLEARED
        %   rather than the click being refused - a refusal with no visible reason is
        %   the worse failure.
        keys = selectedKeys(ctl);
        if numel(keys)<2, return; end
        sigs = cellfun(@(k) signatureOfKey(ctl,k), keys, 'uni',0);
        anchor = sigs{end};                       % the entry the user just added
        keep = keys(strcmp(sigs, anchor));
        if numel(keep)<numel(keys)
            ctl.Value = keep;
            % A hint rather than a status line: the render that follows would
            % otherwise overwrite the one visible explanation of what just happened.
            app.hint = ['Those cannot share one axes, so the earlier choice was ' ...
                        'dropped.  Pick entries of the same shape to compare them.'];
        end
    end

    function s = signatureOfKey(ctl, key)
        if isequal(ctl, c.signal), rows = rowsForFamily(); f='signal'; else, rows = rowsForSignal(); f='varId'; end
        j = find(strcmp({rows.(f)}, key),1);
        if isempty(j), s = ''; else, s = rows(j).sig; end
    end

    function [keys,labels] = signalItems(rows)
        keys = uniqueStable({rows.signal}); labels = cell(1,numel(keys));
        for i = 1:numel(keys)
            j = find(strcmp({rows.signal}, keys{i}),1);
            labels{i} = signalName(rows(j));
        end
    end

    function [keys,labels] = variableItems(rows)
        keys = uniqueStable({rows.varId}); labels = cell(1,numel(keys));
        for i = 1:numel(keys)
            j = find(strcmp({rows.varId}, keys{i}),1);
            labels{i} = rows(j).varLabel;
            if rows(j).suspect, labels{i} = [labels{i} '  (check this one)']; end
        end
    end

    function nm = signalName(r)
        %signalName  The signal as the user should see it: a wire channel's REAL name
        %   when the tree kept it, the measure spelled out when it is a diameter, and
        %   the field name only when nothing better exists.
        nm = r.signalLabel;
        if ~isempty(nm), return; end
        switch r.signal
            case 'sData',       nm = 'segments';
            case 'gsData',      nm = 'segments, full time resolution';
            case 'dvsData',     nm = 'tracked vessels';
            case 'dvsDiameter', nm = 'tracked vessel diameter';
            case 'ppx',         nm = 'every pixel';
            case 'sMetrics',    nm = 'segments';
            case 'dvsMetrics',  nm = 'tracked vessels';
            case '',            nm = '(this recording)';
            otherwise,          nm = exMyograph('measure',r.signal);
        end
    end

%% ==================== RENDER DISPATCH =================================== %%
    function requestRender()
        % Cheap guard so we do not attempt to plot before a source is chosen.
        if isempty(app.files), return; end
        doRender();
    end

    function doRender()
        if isempty(app.files), return; end
        cla(app.ax,'reset'); app.ax.Color='w';
        delete(findobj(app.fig,'Type','ColorBar'));   % drop a colorbar left by a previous image
        hold(app.ax,'on');
        try
            renderInto(app.ax);
        catch ME
            hold(app.ax,'off');
            title(app.ax,'Could not plot this combination');
            setStatus(['Plot error: ' ME.message]);
            app.hint = '';
            return;
        end
        hold(app.ax,'off');
        app.hint = '';   % a hint lives for exactly the render that follows it
    end

    function renderInto(ax)
        %renderInto  THE one dispatch from a PLOT to its renderer, shared by the live
        %   axes and the export axes so the exported figure can never be a different
        %   plot from the one on screen.  It dispatches on the plot the cascade
        %   chose - not on a data type - and every renderer is handed numbers that
        %   exFetch has already shaped, so adding a plot does not touch the cascade.
        app.notes = {};
        P = currentPlot();
        if isempty(P), title(ax,'Choose what to plot'); return; end
        switch P
            case {'curve.time','curve.f','curve.position'}, renderCurve(ax);
            case {'box','bar'},       renderMetric(ax);
            case 'curve.harmonic',    renderHarmonics(ax);
            case 'curve.pct',         renderAmpPct(ax, onePayload());
            case 'curves.family',     renderPctSpectra(ax, onePayload());
            case {'heat.ft','heat.fpct'}, renderHeatMap(ax, onePayload());
            case 'heat.positionTime', renderDiameterMap(ax, onePayload());
            case {'image','image.fromSegments','image.timeAverage','image.frame'}
                                      renderImage(ax, onePayload());
            case 'video',             renderVideo(ax, onePayload());
            case 'lines.overlay',     renderPerLineDiameter(ax, onePayload());
            case 'image.walls',       renderWalls(ax, onePayload());
            otherwise
                error('guiExplore:noRenderer','%s has no renderer.', char(P));
        end
    end

%% ==================== OBSERVATION BUILDERS ============================== %%
    function [obs,meta] = gatherCurves()
        %gatherCurves  Every curve the current cascade asks for, flattened and tagged.
        %   ONE observation per (file, signal, variable, analysed window, selection,
        %   unit-curve), which is exactly the product the "Compare by" machinery then
        %   folds back up into series.
        P = currentPlot();
        sels = selectionList();
        rows = rowsForVariable();
        obs = emptyObs();
        meta = struct('xlabel','','ylabel','','xlog',false);
        for i = selectedFileIdx()
            pth = app.files(i).path;
            R = tryLoad(pth);
            if isempty(R), continue; end
            for r = 1:numel(rows)
                ds = descriptorsFor(pth, rows(r));
                for k = 1:numel(ds)
                    iv = intervalNameOf(R, ds(k));
                    if ~intervalWanted(iv), continue; end
                    for s = 1:numel(sels)
                        p = fetchPayload(R, pth, ds(k), rows(r), P, sels{s});
                        noteFrom(p);
                        if ~p.ok || isempty(p.Y), continue; end
                        meta.xlabel = p.xlab; meta.ylabel = p.ylab; meta.xlog = p.xlog;
                        [Y,meta.ylabel] = asBaseline(p.Y, meta.ylabel);
                        for col = 1:size(Y,2)
                            obs(end+1) = mkObs(p.x, Y(:,col), sels{s}, i, iv, rows(r)); %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    function [vals,tags,meta] = gatherScalars()
        %gatherScalars  The same walk, for a plot whose payload is one number per
        %   unit rather than a curve.
        P = currentPlot();
        sels = selectionList();
        rows = rowsForVariable();
        vals = []; tags = emptyTags();
        meta = struct('xlabel','','ylabel','','xlog',false);
        for i = selectedFileIdx()
            pth = app.files(i).path;
            R = tryLoad(pth);
            if isempty(R), continue; end
            for r = 1:numel(rows)
                ds = descriptorsFor(pth, rows(r));
                for k = 1:numel(ds)
                    iv = intervalNameOf(R, ds(k));
                    if ~intervalWanted(iv), continue; end
                    for s = 1:numel(sels)
                        p = fetchPayload(R, pth, ds(k), rows(r), P, sels{s});
                        noteFrom(p);
                        if ~p.ok, continue; end
                        meta.ylabel = p.ylab;
                        v = double(p.Y(:)); v = v(isfinite(v));
                        for j = 1:numel(v)
                            vals(end+1,1) = v(j); %#ok<AGROW>
                            tags(end+1)   = mkTag(sels{s}, i, iv, rows(r)); %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    function P = onePayload()
        %onePayload  THE single-recording views: one file, one signal, one variable,
        %   one window, one selection.  Several selected is not an error - the first
        %   is shown and the status bar says so, which is what the tool has always
        %   done for a picture.
        P = [];
        rows = rowsForVariable();
        sels = selectionList();
        idx = selectedFileIdx();
        if isempty(rows) || isempty(idx), return; end
        nDrawn = 0;
        for i = idx
            pth = app.files(i).path;
            R = tryLoad(pth);
            if isempty(R), continue; end
            ds = descriptorsFor(pth, rows(1));
            for k = 1:numel(ds)
                iv = intervalNameOf(R, ds(k));
                if ~intervalWanted(iv), continue; end
                nDrawn = nDrawn + 1;
                if nDrawn>1, continue; end
                P = fetchPayload(R, pth, ds(k), rows(1), currentPlot(), sels{1});
                P.title = seriesTitle(app.files(i), iv, rows(1), sels{1});
                noteFrom(P);
            end
        end
        more = numel(idx)*max(numel(rows),1)*numel(sels);
        if nDrawn>1 || more>1
            setStatus('This view shows one recording and one selection; showing the first.');
        end
    end

% ---- resolving a menu row back to the leaves behind it ---------------------
    function ds = descriptorsFor(pth, row)
        %descriptorsFor  The concrete leaves one menu row stands for, in one file.  A
        %   myograph recording answers with one per analysed WINDOW, which is what
        %   puts the window on an axis; everything else answers with exactly one.
        [D,S] = indexOf(pth);
        if ~isempty(row.srcView)
            ds = S; if isempty(ds), return; end
            ds = ds(strcmp({ds.key}, row.key));
            return
        end
        ds = D;
        if isempty(ds), return; end
        keep = false(1,numel(ds));
        for i = 1:numel(ds)
            keep(i) = strcmp([ds(i).family '||' ds(i).signal '||' varIdOf(ds(i))], row.key);
        end
        ds = ds(keep);
    end

    function d = firstDescriptor(pth, row)
        ds = descriptorsFor(pth, row);
        if isempty(ds), d = []; else, d = ds(1); end
    end

    function P = fetchPayload(R, pth, d, row, plotId, sel)
        %fetchPayload  THE seam.  Everything a renderer draws comes through here, so a
        %   renderer never reaches into a results tree and session 4 can add one
        %   without touching the cascade.
        if ~isempty(row.srcView)
            P = fetchSourceView(R, pth, d, plotId);
            return
        end
        opts = struct('select',sel, 'units',unitsChoice(), 'path',pth, ...
                      'rows',[], 'interval',[], 'slice',sliceNow(), 'scrub',app.want);
        try
            P = exFetch(R, d, plotId, opts);
        catch ME
            P = struct('ok',false,'note',ME.message,'Y',[],'x',[],'img',[],'cdata',[], ...
                'yvals',[],'frames',[],'fvals',[],'names',{{}},'w',[], ...
                'xlab','','ylab','','clab','','flab','','xlog',false,'xscale','linear', ...
                'nUnits',0,'plot',plotId,'kind','');
        end
        if ~isfield(P,'title'), P.title = ''; end
    end

    function nm = intervalNameOf(R, d)
        %intervalNameOf  Which analysed window a leaf belongs to, '' when the
        %   recording has none.  A wire myograph hangs its windows off the CHANNEL, so
        %   the flat index has to be recomposed - two chambers both calling their
        %   first window 'interval1' would otherwise be averaged together.
        nm = '';
        if isfield(d,'srcView') && ~isempty(d.srcView)
            nm = exMyograph('name', R, d.srcInterval); return
        end
        k = exMyograph('index', R, d.path);
        if k>0, nm = exMyograph('name', R, k); end
    end

    function tf = intervalWanted(nm)
        want = c.interval.Value;
        tf = isempty(want) || strcmp(want,'(all)') || strcmp(nm, want);
    end

    function s = selectionList()
        s = c.selList.Value;
        if ischar(s), s = {s}; end
        if isempty(s), s = {''}; end
    end

    function u = unitsChoice()
        switch c.units.Value
            case 'Pooled (area-weighted)', u = 'pooled';
            case 'One point each',         u = 'points';
            case 'One curve each',         u = 'curves';
            otherwise,                     u = '';      % let the plot decide
        end
    end

    function [Y,ylab] = asBaseline(Y, ylab)
        %asBaseline  Each curve against ITS OWN mean, which is how a myograph result
        %   is usually read: vessels of different calibre become comparable and a
        %   dilation is a number a reader recognises.
        if ~c.baseline.Value, return; end
        for j = 1:size(Y,2)
            b = mean(Y(:,j),'omitnan');
            if isfinite(b) && b~=0, Y(:,j) = 100*Y(:,j)/b; else, Y(:,j) = NaN; end
        end
        ylab = [ylab ', % of baseline'];
    end

    function noteFrom(p)
        %noteFrom  A CAVEAT ON A DRAWABLE PAYLOAD IS STILL A CAVEAT.  A file too big
        %   to load has no readable metrics table, so its units come back unselected
        %   and equally weighted - and swallowing that is how a user comes to believe
        %   a big file was area-weighted.
        if isempty(p) || ~isstruct(p) || ~isfield(p,'note') || isempty(p.note), return; end
        n = strtrim(p.note);
        if ~any(strcmp(app.notes,n)), app.notes{end+1} = n; end
    end

    function o = mkObs(x,y,sel,i,interval,row)
        o = struct('x',x(:),'y',y(:),'sel',prettySel(sel), ...
            'group',app.files(i).group,'rec',app.files(i).rec, ...
            'animal',app.files(i).animal,'type',app.files(i).type, ...
            'interval',interval,'file',app.files(i).name, ...
            'var',row.varLabel,'signal',signalName(row));
    end

    function t = mkTag(sel,i,interval,row)
        t = struct('sel',prettySel(sel),'group',app.files(i).group, ...
            'rec',app.files(i).rec,'animal',app.files(i).animal, ...
            'type',app.files(i).type,'interval',interval,'file',app.files(i).name, ...
            'var',row.varLabel,'signal',signalName(row));
    end

%% ==================== RENDERERS ========================================= %%
% Every renderer below is handed numbers that exFetch has already shaped.  None of
% them reaches into a results tree, and none of them knows which container the
% numbers came from - which is what lets a new plot be added without touching the
% cascade, and what lets the same renderer serve a speckle recording and a myograph
% one.

    function renderCurve(ax)
        [obs,meta] = gatherCurves();
        if isempty(obs), title(ax,'No data for this selection'); setStatus(noteText()); return; end
        [key,leg] = seriesKeys(obs);
        uk = unique(key,'stable');
        cols = seriesColours(leg(firstIndex(key,uk)));
        [statFun,loFun,hiFun,statName] = statFuns();
        xg = commonGrid(obs);
        hLeg = gobjects(1,numel(uk)); names = cell(1,numel(uk));
        for m = 1:numel(uk)
            sel = strcmp(key,uk{m});
            Y = resampleObs(obs(sel), xg);
            mu = statFun(Y,2); lo = loFun(Y,2); hi = hiFun(Y,2);
            band = [lo;flipud(hi)]; xx=[xg(:);flipud(xg(:))];
            good = isfinite(xx)&isfinite(band);
            fill(ax, xx(good), band(good), cols(m,:), 'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
            if size(Y,2)>1              % show the individual constituents faintly
                plot(ax, xg, Y, '-','Color',[cols(m,:) 0.25],'LineWidth',0.5, ...
                    'HandleVisibility','off');
            end
            hLeg(m) = plot(ax, xg, mu, '-','Color',cols(m,:),'LineWidth',1.8);
            names{m} = legendName(leg{firstIndex(key,uk(m))});
        end
        if numel(uk)<2, hLeg=gobjects(0); names={}; end   % nothing to distinguish -> no legend
        finishAxes(ax, meta.xlabel, meta.ylabel, hLeg, names);
        if meta.xlog, set(ax,'XScale','log'); end
        axis(ax,'tight'); yl=ylim(ax); ylim(ax, yl+[-1 1]*0.03*range(yl)+[ -eps eps ]);
        setStatus(sprintf('%d series, %d observations (%s).%s%s', ...
            numel(uk), numel(obs), statName, seriesWarning(numel(uk)), noteText()));
    end

    function renderHarmonics(ax)
        %renderHarmonics  A HARMONIC NUMBER IS ORDINAL, NOT METRIC.  There is no
        %   half-harmonic, so nothing is interpolated between the stems and the ticks
        %   are the integers themselves; a line drawn through them would suggest a
        %   continuum that the cardiac series does not have.  Across recordings they
        %   are grouped exactly as the curves are, centre and spread from the same
        %   Centre / error control, with the series dodged so two do not overprint.
        [obs,meta] = gatherCurves();
        if isempty(obs), title(ax,'No data for this selection'); setStatus(noteText()); return; end
        [key,leg] = seriesKeys(obs);
        uk = unique(key,'stable');
        cols = seriesColours(leg(firstIndex(key,uk)));
        [statFun,loFun,hiFun,statName] = statFuns();
        xg = commonGrid(obs);
        n = numel(uk);
        step = 0; if n>1, step = 0.6/n; end
        hLeg = gobjects(1,n); names = cell(1,n);
        for m = 1:n
            sel = strcmp(key,uk{m});
            Y  = resampleObs(obs(sel), xg);
            mu = statFun(Y,2); lo = loFun(Y,2); hi = hiFun(Y,2);
            xx = xg(:) + (m-(n+1)/2)*step;
            hLeg(m) = stem(ax, xx, mu, 'filled','Color',cols(m,:), ...
                'MarkerFaceColor',cols(m,:),'MarkerSize',5,'LineWidth',1.4);
            if any(isfinite(hi-lo))
                errorbar(ax, xx, mu, mu-lo, hi-mu, 'LineStyle','none', ...
                    'Color',cols(m,:),'CapSize',6,'HandleVisibility','off');
            end
            names{m} = legendName(leg{firstIndex(key,uk(m))});
        end
        if n<2, hLeg = gobjects(0); names = {}; end
        finishAxes(ax, meta.xlabel, meta.ylabel, hLeg, names);
        set(ax,'XScale','linear');
        ticks = unique(round(xg(:)));
        if ~isempty(ticks) && numel(ticks)<=30, set(ax,'XTick',ticks); end
        if numel(xg)>0, xlim(ax, [min(xg)-0.6, max(xg)+0.6]); end
        setStatus(sprintf('%d series, %d observations (%s).%s%s', ...
            n, numel(obs), statName, seriesWarning(n), noteText()));
    end

    function renderMetric(ax)
        [vals,tags,meta] = gatherScalars();
        if isempty(vals), title(ax,'No data for this selection'); legend(ax,'off'); ...
                setStatus(noteText()); return; end
        dims = activeDims(tags);
        xName = dims{1}; colDims = dims(2:end);
        xc  = arrayfun(@(t) t.(xName), tags, 'uni',0);
        xcats = orderedCats(xc);
        colName = meta.ylabel;
        if ~isempty(colDims)
            cc = arrayfun(@(t) strjoin(cellfun(@(d) t.(d), colDims,'uni',0),' | '), tags, 'uni',0);
            ccats = unique(cc,'stable'); ccats = ccats(orderNumeric(ccats));
            cols  = seriesColours(ccats);
            ax.ColorOrder = cols;
            boxchart(ax, categorical(xc,xcats,'Ordinal',true), vals, ...
                'GroupByColor', categorical(cc,ccats), 'MarkerStyle','none');
            overlayPoints(ax, xc, cc, vals, xcats, ccats, cols, true);
            if c.legendChk.Value
                lg = legend(ax, cellfun(@legendDyn,ccats,'uni',0), ...
                    'Interpreter','none','Location','best'); lg.Box='off'; lg.FontSize=c.fontSize.Value;
            else
                legend(ax,'off');
            end
        else
            cols = seriesColours(xcats);
            for k = 1:numel(xcats)
                selk = strcmp(xc, xcats{k}); if ~any(selk), continue; end
                b = boxchart(ax, categorical(xc(selk),xcats,'Ordinal',true), vals(selk), ...
                    'MarkerStyle','none');
                b.BoxFaceColor = cols(k,:); b.WhiskerLineColor = cols(k,:); b.BoxFaceAlpha = 0.45;
            end
            overlayPoints(ax, xc, xc, vals, xcats, xcats, cols, false);
            legend(ax,'off');
        end
        finishAxes(ax, xLabelFor(xName), colName, gobjects(0), {});
        yl = ylim(ax); if diff(yl)>0, ylim(ax, yl+[-0.06 0.06]*diff(yl)); end  % a little head/foot room
        setStatus(sprintf('%d x-category(ies), %d points.%s%s', numel(xcats), numel(vals), ...
            seriesWarning(numel(xcats)*max(numel(colDims),1)), noteText()));
    end

    function renderHeatMap(ax, P)
        %renderHeatMap  A HEAT MAP WITH FREQUENCY UP THE SIDE - time along the bottom
        %   for a spectrogram, percentile for the flare/silence view.  The two are the
        %   same picture of the same array and share this one setup, so the frequency
        %   axis of one cannot drift away from the frequency axis of the other.
        %   ROWS ARE THE Y AXIS AND COLUMNS THE X - the orientation exFetch already
        %   returns - so nothing is transposed on the way in.
        if ~drawable(ax,P,'No heat map for this selection'), return; end
        imagesc(ax, P.x, P.yvals, P.cdata); set(ax,'YDir','normal');
        if all(P.yvals(:)>0)
            try, set(ax,'YScale','log'); catch, end
        end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String = P.clab;
        finishAxes(ax, P.xlab, P.ylab, gobjects(0), {}); legend(ax,'off');
        grid(ax,'off'); axis(ax,'tight');
        setStatus([P.title noteText()]);
    end

    function renderPctSpectra(ax, P)
        %renderPctSpectra  A FAMILY OF CURVES on one axes - one per percentile bin of
        %   the envelope amplitude, or one per position along a vessel.  Which of the
        %   two it is was decided by the rules; here they are simply columns of Y with
        %   names beside them.
        if ~drawable(ax,P,'No family of curves for this selection'), return; end
        n = size(P.Y,2);
        hLeg = gobjects(1,n); names = cell(1,n);
        cols = distinctPalette(n);
        for p = 1:n
            hLeg(p) = plot(ax, P.x, P.Y(:,p), '-', 'Color', cols(p,:), 'LineWidth',1.2);
            names{p} = seriesLabel(P.names, p);
        end
        if n>12, hLeg = hLeg([]); names = {}; end   % a legend of 480 lines is not a legend
        finishAxes(ax, P.xlab, P.ylab, hLeg, names);
        if P.xlog, set(ax,'XScale','log'); end
        axis(ax,'tight');
        setStatus(sprintf('%s%d curve(s).%s', prefixOf(P.title), n, noteText()));
    end

    function renderAmpPct(ax, P)
        %renderAmpPct  A curve against the percentile LEVELS - which are read from the
        %   recording's own settings file rather than assumed, so a run with a
        %   non-default set of percentiles is drawn on its own x and not on 0:10:100.
        if ~drawable(ax,P,'No percentile curve for this selection'), return; end
        n = size(P.Y,2);
        hLeg = gobjects(1,n); names = cell(1,n);
        cols = distinctPalette(n);
        for p = 1:n
            hLeg(p) = plot(ax, P.x, P.Y(:,p), '-o', 'Color', cols(p,:), 'LineWidth',1.4, ...
                'MarkerSize',4);
            names{p} = seriesLabel(P.names, p);
        end
        if n<2 || n>12, hLeg = gobjects(0); names = {}; end
        finishAxes(ax, P.xlab, P.ylab, hLeg, names);
        set(ax,'XScale','linear');
        setStatus([P.title noteText()]);
    end

    function renderImage(ax, P)
        %renderImage  A PICTURE, however it was arrived at - a per-pixel map, a
        %   per-segment number painted back through sMap, one frame of a cube, or a
        %   cube averaged over time.  All four are one [Y x X] array by the time they
        %   reach here, which is the point: the difference between them is a question
        %   about the DATA and belongs in exFetch, and the only thing that must not be
        %   lost is which one it was - so the title says so and never lets a
        %   mean-over-time map pass for a single frame.
        %   Whatever the selection did not pick came back NaN, which is drawn
        %   transparent rather than as a hole of zeros.
        if ~drawable(ax,P,'No map for this selection'), return; end
        img = double(P.img);
        m = isfinite(img);
        imagesc(ax, img, 'AlphaData', m); axis(ax,'image');
        if any(m(:))
            lims = prctile(img(m),[2 99]);
            if lims(2)>lims(1), clim(ax,lims); end
        end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String = P.clab;
        finishAxes(ax, '', '', gobjects(0), {}); legend(ax,'off');
        set(ax,'XTick',[],'YTick',[]); grid(ax,'off');
        setStatus([P.title imageNote(P) noteText()]);
    end

    function s = imageNote(P)
        %imageNote  What this particular picture is, said in the status bar as well as
        %   in the title - the two views that look alike are the ones worth naming.
        switch char(P.plot)
            case 'image.fromSegments'
                s = sprintf('%d segment(s) painted back into the map.  ', P.nUnits);
            case 'image.timeAverage', s = 'averaged over time.  ';
            case 'image.frame',       s = [scrubText(app.scrub, app.scrub.idx, false) '.  '];
            otherwise,                s = '';
        end
    end

    function renderVideo(ax, P)
        %renderVideo  ONE FRAME OF THE SWEEP, with the slider under the axes saying
        %   where it is and "Export video" writing the whole of it.  THE COLOUR LIMITS
        %   ARE FIXED OVER THE WHOLE CUBE and not over the frame: rescaling each frame
        %   to its own range makes a recording that barely changes look alive, which is
        %   the one thing a movie of a vascular field must not do.
        if ~drawable(ax,P,'No video for this selection'), return; end
        nF = size(P.frames,3);
        k  = min(max(app.scrub.idx,1),nF);
        img = double(P.frames(:,:,k));
        m = isfinite(img);
        imagesc(ax, img, 'AlphaData', m); axis(ax,'image');
        lims = cubeLimits(P.frames);
        if ~isempty(lims), clim(ax,lims); end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String = P.clab;
        finishAxes(ax, '', '', gobjects(0), {}); legend(ax,'off');
        set(ax,'XTick',[],'YTick',[]); grid(ax,'off');
        setStatus(sprintf('%s%s, frame %d of %d.  "Export video" writes the whole sweep.%s', ...
            prefixOf(P.title), frameText(P,k), k, nF, noteText()));
    end

    function renderDiameterMap(ax, P)
        %renderDiameterMap  POSITION AGAINST TIME - the kymograph.  For a myograph
        %   recording this is the view that says whether detection held over the whole
        %   vessel for the whole window: a line that wandered shows as a stripe, and a
        %   wall that was lost shows as a band of missing values.  No myograph step
        %   writes a report page, so this is where that check lives.
        if ~drawable(ax,P,'No position-time map for this selection'), return; end
        imagesc(ax, P.x, P.yvals, P.cdata); set(ax,'YDir','normal'); axis(ax,'tight');
        m = isfinite(P.cdata);
        if any(m(:))
            lims = prctile(P.cdata(m),[2 98]);
            if lims(2)>lims(1), clim(ax,lims); end
        end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String = P.clab;
        finishAxes(ax, P.xlab, P.ylab, gobjects(0),{});
        legend(ax,'off'); grid(ax,'off');
        setStatus(sprintf('%s%d position(s).%s', prefixOf(P.title), size(P.cdata,1), noteText()));
    end

    function renderPerLineDiameter(ax, P)
        %renderPerLineDiameter  EVERY LINE'S OWN DIAMETER against time, rather than
        %   their average or their map.  The map says whether detection HELD; this says
        %   what the individual traces DID - whether the lines move together, which is
        %   the assumption a propagation speed is fitted under, and whether one stray
        %   line is carrying an average the rest do not support.  Lines are coloured by
        %   position so they can be read against the map, and the line average is drawn
        %   over them in black.
        if ~drawable(ax,P,'The per-line diameter is not beside this result - keep the _MYO_d.mat with it'), return; end
        M = P.Y;
        cm = colormap(ax, c.cmap.Value);
        ci = round(linspace(1, size(cm,1), max(size(M,2),1)));
        for j = 1:size(M,2)
            plot(ax, P.x, M(:,j), '-', 'Color',[cm(ci(j),:) 0.45],'LineWidth',0.5, ...
                'HandleVisibility','off');
        end
        hMean = plot(ax, P.x, mean(M,2,'omitnan'), '-','Color',[0 0 0],'LineWidth',1.8);
        if numel(P.yvals)>1
            clim(ax,[P.yvals(1) P.yvals(end)]); cb = colorbar(ax); cb.Label.String = 'position';
        end
        finishAxes(ax, P.xlab, P.ylab, hMean, {'line average'});
        axis(ax,'tight');
        setStatus(sprintf('%s%d line(s).%s', prefixOf(P.title), size(M,2), noteText()));
    end

    function renderWalls(ax, P)
        %renderWalls  ONE FRAME OF THE RECORDING WITH THE DETECTED WALLS ON IT - the
        %   "did it find the right edges" check, and the other half of what a report
        %   page would have carried.  The walls are drawn RED when that frame was
        %   flagged invalid, which is a wall that had left the field of view: its
        %   diameter is a lower bound rather than a measurement, and that has to be
        %   visible rather than inferred.
        if isempty(P) || ~isstruct(P) || ~isfield(P,'walls') || isempty(P.walls.frame)
            title(ax,'The recording is not beside this result - keep the video with it');
            setStatus(noteText()); return
        end
        W = P.walls;
        image(ax, W.frame); axis(ax,'image');
        col = [0.1 0.9 0.1]; if ~W.valid, col = [0.95 0.15 0.15]; end
        plot(ax, W.left,  W.rows, '-', 'Color', col, 'LineWidth', 1.2);
        plot(ax, W.right, W.rows, '-', 'Color', col, 'LineWidth', 1.2);
        finishAxes(ax,'','',gobjects(0),{}); legend(ax,'off');
        set(ax,'XTick',[],'YTick',[]); grid(ax,'off');
        note = ''; if ~W.valid, note = '   wall out of view'; end
        setStatus(sprintf('%s%.1f s%s%s', prefixOf(P.title), W.time, note, noteText()));
    end

    function tf = drawable(ax, P, msg)
        tf = ~isempty(P) && isstruct(P) && isfield(P,'ok') && P.ok;
        if tf, return; end
        why = msg;
        if ~isempty(P) && isstruct(P) && isfield(P,'note') && ~isempty(P.note), why = P.note; end
        title(ax, why, 'Interpreter','none');
        setStatus([why noteText()]);
    end

%% ============ THE VIEWS THAT READ THE RECORDING'S SOURCE ================ %%
% Three views of a myograph recording are NOT in the results tree and so cannot be
% in the variable index: the per-line diameters, their position-time map and the
% detected walls all live once, in the SOURCE beside the results (a window carries
% only the line-averaged trace).  They are folded into the same cascade by being
% offered as variables of their own, and answered here in the payload shape every
% other renderer already eats.

    function S = sourceViews(R, pth)
        %sourceViews  What the recording beside this result can also show.  Offered
        %   only when the source is actually there, so a results file moved away from
        %   its triplet simply does not offer them rather than offering and failing.
        S = emptyVars();
        if isempty(R) || ~strcmp(exModality(R),'myograph'), return; end
        if ~isfile(strrep(char(pth),'_r.mat','_d.mat')), return; end
        % THE SIGNAL IS THE TREE'S OWN FIELD, not the measure's display name.  These
        % three views sit under the same Signal step as the diameter the tree stores,
        % and a second key spelling the same measure differently would put every
        % measure in that list TWICE, under one name each time.
        sig  = exMyograph('signals', R);
        meas = {sig(strcmp({sig.kind},'measure')).field};
        if isempty(meas), return; end
        nIv = exMyograph('count', R);
        spec = {'perLine','Every line''s own diameter','Vector',{'lines.overlay'}; ...
                'map',    'Diameter along the vessel', 'Image', {'heat.positionTime'}; ...
                'walls',  'Detected walls',            'Image', {'image.walls'}};
        for m = 1:numel(meas)
            for v = 1:size(spec,1)
                for k = 1:nIv
                    r = emptyVars();
                    r(1).family   = 'Diameter';
                    r(1).signal   = meas{m};
                    r(1).signalLabel = '';
                    r(1).suspect  = false;
                    r(1).varId    = ['src|' spec{v,1}];
                    r(1).varLabel = spec{v,2};
                    r(1).leaf     = spec{v,1};
                    r(1).unit     = 'line';
                    r(1).dims     = {'time'};
                    r(1).sig      = ['src|' spec{v,1}];
                    r(1).kinds    = spec(v,3);
                    r(1).plots    = spec{v,4};
                    r(1).srcView  = spec{v,1};
                    r(1).srcInterval = k;
                    r(1).nFiles   = 1;
                    r(1).key      = ['Diameter||' meas{m} '||src|' spec{v,1}];
                    S(end+1) = r; %#ok<AGROW>
                end
            end
        end
    end

    function P = fetchSourceView(R, pth, d, plotId)
        %fetchSourceView  exFetch's counterpart for the three source-backed views: the
        %   same payload struct, filled from the recording rather than from the tree.
        P = struct('ok',false,'note','','plot',plotId,'kind','','x',[],'Y',[],'w',[], ...
            'names',{{}},'img',[],'cdata',[],'yvals',[],'frames',[],'fvals',[], ...
            'xlab','','ylab','','clab','','flab','','xlog',false,'xscale','linear', ...
            'nUnits',0,'title','','walls',[]);
        k = d.srcInterval;
        if strcmp(d.srcView,'walls')
            W = exMyograph('walls', R, pth, k, d.signal);
            P.walls = W;
            P.ok = ~isempty(W.frame);
            if ~P.ok, P.note = 'The recording is not beside this result - keep the video with it.'; end
            return
        end
        [M,t,y,unit] = exMyograph('diametermap', R, pth, k, d.signal);
        if isempty(M)
            P.note = 'The per-line diameter is not beside this result - keep the _MYO_d.mat with it.';
            return
        end
        P.xlab = 'time (s)'; P.nUnits = size(M,2);
        P.x = double(t(:)); P.yvals = double(y(:));
        switch d.srcView
            case 'perLine'
                [P.Y, P.ylab] = asBaseline(double(M), ['diameter (' unit ')']);
                P.names = arrayfun(@(i) sprintf('position %d',i), y(:).', 'UniformOutput', false);
            case 'map'
                P.cdata = double(M).';                 % rows = position, columns = time
                P.ylab  = 'position along the vessel';
                P.clab  = ['diameter (' unit ')'];
        end
        P.ok = true;
    end

    function overlayPoints(ax, xc, cg, vals, xcats, cgcats, cols, dodge)
        % Individual observations as jittered points, colour-matched to the boxes.
        nC = numel(cgcats); bw = 0.75/max(nC,1);
        for j = 1:nC
            sel = strcmp(cg, cgcats{j}); if ~any(sel), continue; end
            xi  = double(categorical(xc(sel), xcats, 'Ordinal',true));
            if dodge, off = (j-(nC+1)/2)*bw; else, off = 0; end
            jit = (rand(numel(xi),1)-0.5)*bw*0.5;
            scatter(ax, xi(:)+off+jit, vals(sel), 16, cols(j,:), 'filled', ...
                'MarkerFaceAlpha',0.55,'MarkerEdgeColor','none','HandleVisibility','off');
        end
    end

%% ==================== COSMETICS ========================================= %%
    function finishAxes(ax, xlab, ylab, hLeg, names)
        set(ax,'Color','w'); grid(ax,'on'); ax.GridAlpha=0.12; box(ax,'on');
        ax.TickLabelInterpreter='none';   % keep '_' in labels literal (no TeX subscripts)
        ax.FontSize = c.fontSize.Value;   % base size: ticks match, title/labels scale up (x1.1)
        title(ax, currentTitle(),'Interpreter','none','FontWeight','bold');
        if ~isempty(strtrim(c.xlabF.Value)), xlab = c.xlabF.Value; end
        if ~isempty(strtrim(c.ylabF.Value)), ylab = c.ylabF.Value; end
        xlabel(ax, xlab,'Interpreter','none'); ylabel(ax, ylab,'Interpreter','none');
        cb = findobj(ancestor(ax,'figure'),'Type','ColorBar');   % keep colorbar text in step
        for k=1:numel(cb), cb(k).FontSize=c.fontSize.Value; cb(k).Label.FontSize=c.fontSize.Value; end
        % Legend is only managed here when the caller supplies entries (time series /
        % spectra). Box/metric renderers set their own legend; images set none.
        if ~isempty(names)
            if c.legendChk.Value
                lg = legend(ax, hLeg, names,'Interpreter','none','Location','best'); ...
                    lg.Box='off'; lg.FontSize=c.fontSize.Value;
            else
                legend(ax,'off');
            end
        end
    end

    function dims = activeDims(tg)
        % Ordered list of the tag dimensions that separate the data, given the X axis
        % and Colour choices. The first is the natural x-axis (box plots) and the
        % combination defines a series (curves). Selection, variable and signal are
        % always kept separate when several are shown, so two different quantities are
        % never pooled into one box.
        % GROUP, INDEX, ANIMAL, TYPE and INTERVAL are five INDEPENDENT axes (spec §2)
        % - they never nest, so Auto simply takes whichever of them actually varies.
        % Animal and type are empty unless a session filled them, and interval unless
        % the recording is a myograph one, which is why a folder scan of speckle
        % results behaves exactly as it always did.
        prim = {};
        x = dimKeyOf(c.xaxis.Value);
        if strcmp(c.xaxis.Value,'Pool all')
            prim = {};
        elseif ~isempty(x)
            prim = {x};
        else                                           % Auto: use whatever varies
            if nDim(tg,'rec')>1,      prim{end+1}='rec';      end
            if nDim(tg,'group')>1,    prim{end+1}='group';    end
            if nDim(tg,'animal')>1,   prim{end+1}='animal';   end
            if nDim(tg,'type')>1,     prim{end+1}='type';     end
            if nDim(tg,'interval')>1, prim{end+1}='interval'; end
        end
        col = dimKeyOf(c.colour.Value);
        if ~isempty(col) && ~any(strcmp(prim,col)), prim{end+1}=col; end
        dims=prim;
        for f = {'sel','var','signal'}
            if nDim(tg,f{1})>1 && ~any(strcmp(dims,f{1})), dims{end+1}=f{1}; end %#ok<AGROW>
        end
        if isempty(dims), dims={'sel'}; end            % single group -> one x category
    end

    function k = dimKeyOf(name)
        %dimKeyOf  The X axis / Colour menu entry, as the tag field it stratifies on.
        switch char(name)
            case 'Group',           k='group';
            case 'Recording index', k='rec';
            case 'Animal',          k='animal';
            case 'Type',            k='type';
            case 'Interval',        k='interval';
            case 'Selection',       k='sel';
            case 'Variable',        k='var';
            case 'Signal',          k='signal';
            otherwise,              k='';             % 'Auto', '(none)', 'Pool all'
        end
    end

    function s = seriesWarning(n)
        %seriesWarning  Twelve series is about as many as a reader can follow.  Past
        %   that the status bar says so and NAMES the axis that is multiplying, which
        %   is the one to switch off.
        s = '';
        if n<=12, return; end
        nSel = numel(selectionList());
        nVar = numel(selectedKeys(c.variable));
        nSig = numel(selectedKeys(c.signal));
        nFil = numel(selectedFileIdx());
        [~,w] = max([nSel nVar nSig nFil]);
        axName = {'the selection','the variable','the signal','the file'};
        s = sprintf('  %d series is a lot to read - %s is what multiplies them.', ...
            n, axName{w});
    end

    function s = noteText()
        %noteText  The caveats the last fetch raised, appended to whatever the
        %   renderer had to say.  A payload that is drawable can still carry one.
        bits = app.notes;
        if ~isempty(app.hint), bits = [{app.hint} bits]; end
        if isempty(bits), s = ''; return; end
        s = ['  ' strjoin(bits,'  ')];
    end

    function nm = legendDyn(val)
        % Apply the legend pattern to a single dynamic category (box colour groups).
        L = emptyLegendTag(); L.dyn = val;
        nm = legendName(L);
    end

    function [key,leg] = seriesKeys(obs)
        % A series (one line+band) = unique combination of the ACTIVE dimensions.
        dims=activeDims(obs);
        key=cell(1,numel(obs)); leg=cell(1,numel(obs));
        for i=1:numel(obs)
            key{i}=strjoin(cellfun(@(d) obs(i).(d), dims,'uni',0),'|');
            L=emptyLegendTag(); L.file=obs(i).file;
            for d=1:numel(dims), L.(dims{d})=obs(i).(dims{d}); end
            leg{i}=L;
        end
    end

    function t = currentTitle()
        if app.manual.title && ~isempty(strtrim(c.titleF.Value))
            t = c.titleF.Value; return;
        end
        t = autoTitle(); c.titleF.Value = t;   % keep the box in sync when auto
    end

    function autofillCosmetics()
        app.manual.title=false; app.manual.xlab=false; app.manual.ylab=false;
        c.titleF.Value = autoTitle(); c.xlabF.Value=''; c.ylabF.Value='';
    end

    function t = autoTitle()
        v = currentVarLabel();
        sels = selectionList();
        selTxt = strjoin(cellfun(@prettySel,sels,'uni',0),', ');
        if isempty(v), t = ''; else, t = sprintf('%s - %s', v, selTxt); end
        t = [t plotQualifier()];
        if numel(app.files)>1, t = [t sprintf('  (n=%d files)',numel(selectedFileIdx()))]; end
    end

    function s = plotQualifier()
        %plotQualifier  WHAT HAPPENED TO THE DIMENSION THAT IS NOT ON SCREEN.  A map
        %   averaged over time and a single frame of the same cube are the same
        %   picture to look at and different results to read, so the title says which
        %   one this is rather than leaving the reader to remember.
        s = '';
        if isempty(currentVarLabel()), return; end
        switch currentPlot()
            case 'image.timeAverage',  s = ', mean over time';
            case 'image.fromSegments', s = ', painted from the segments';
            case {'image.frame','video'}
                if ~isempty(app.scrub.name)
                    k = app.scrub.idx;
                    if ~isempty(app.titleFrame), k = app.titleFrame; end
                    s = [', ' scrubText(app.scrub, k, false)];
                end
        end
    end

    function s = currentVarLabel()
        %currentVarLabel  What the plot is OF, in the words the menu used - several
        %   names when several variables are being compared.
        rows = rowsForVariable();
        if isempty(rows), s = ''; return; end
        s = strjoin(uniqueStable({rows.varLabel}), ', ');
    end

    function s = seriesTitle(entry, iv, row, sel)
        %seriesTitle  What a single-recording view is showing, for the status bar.
        bits = {shortName(entry.path), iv, signalName(row), prettySel(sel)};
        bits = bits(~cellfun(@isempty,bits));
        s = [strjoin(bits,' - ') ':  '];
    end

    function s = prefixOf(t)
        s = ''; if ~isempty(t), s = t; end
    end

    function nm = legendName(tagLike)
        pat = c.legendF.Value; if isempty(pat), pat='%s%g%r'; end
        if isfield(tagLike,'dyn') && ~isempty(tagLike.dyn)
            % dynamic single-dimension legend (box colour groups)
            nm = pat;
            for tk = {'%s','%g','%r','%a','%t','%i','%f','%v'}, nm = strrep(nm,tk{1},''); end
            nm = strtrim([strtrim(nm) ' ' tagLike.dyn]);
            if isempty(strtrim(nm)), nm=tagLike.dyn; end
            return;
        end
        nm = pat;
        nm = strrep(nm,'%s', pad0(tagLike.sel));
        nm = strrep(nm,'%g', pad0(tagLike.group));
        nm = strrep(nm,'%r', pad0(tagLike.rec));
        nm = strrep(nm,'%a', pad0(tagLike.animal));
        nm = strrep(nm,'%t', pad0(tagLike.type));
        nm = strrep(nm,'%i', pad0(tagLike.interval));
        nm = strrep(nm,'%f', pad0(tagLike.file));
        nm = strrep(nm,'%n', pad0(tagLike.signal));
        if isfield(tagLike,'var') && ~isempty(tagLike.var)
            nm = strrep(nm,'%v', pad0(tagLike.var));
        else
            nm = strrep(nm,'%v', currentVarLabel());
        end
        nm = strtrim(regexprep(nm,'\s+',' '));
        if isempty(nm), nm = tagLike.sel; end
    end

%% ==================== LOADING (cache + lazy HDF5) ======================= %%
    function R = tryLoad(path)
        try, R = loadResults(path); catch, R = []; end
    end

    function R = loadResults(path)
        if app.cache.isKey(path)
            app.cacheOrder = [path, setdiff(app.cacheOrder,path,'stable')];
            R = app.cache(path); return;
        end
        if isLazyFile(path)
            R = struct('x_lazy_',true,'x_path_',path);  % handle; fields read on demand
        else
            S = load(path,'results');
            if ~isfield(S,'results'), error('File has no "results" variable: %s',path); end
            R = S.results;
        end
        app.cache(path) = R;
        app.cacheOrder = [path, app.cacheOrder];
        while numel(app.cacheOrder) > app.cacheLimit
            drop = app.cacheOrder{end}; app.cacheOrder(end)=[];
            if app.cache.isKey(drop), app.cache.remove(drop); end
            % THE INDEX GOES WITH THE FILE.  A descriptor array outliving the struct
            % it was scanned from is how a menu comes to offer a variable that the
            % re-read file no longer has.
            if app.index.isKey(drop), app.index.remove(drop); end
        end
    end

    function tf = isLazyFile(path)
        d = dir(path); tf = ~isempty(d) && d(1).bytes > app.sizeLimitGB*1e9;
    end

%% ==================== PROGRAMMATIC API (for testing) =================== %%
    function loadPathsProgrammatic(paths, includePat, groupPat, recPat)
        if nargin<2||isempty(includePat), includePat='_r\.mat$'; end
        if nargin<3, groupPat=''; end
        if nargin<4, recPat=''; end
        app.mode='folder'; app.sessionPath='';
        if ~isempty(paths), app.root = fileparts(paths{1}); end
        e = emptyEntries();
        for i=1:numel(paths)
            [~,nm,ex]=fileparts(paths{i});
            if ~isempty(includePat) && isempty(regexp([nm ex],includePat,'once')), continue; end
            e(end+1) = makeEntryStruct(paths{i}, groupPat, recPat); %#ok<AGROW>
        end
        releaseModality();
        app.files = applyGroupOverrides(sortEntries(e));
        c.includeF.Value=includePat; c.groupF.Value=groupPat; c.recF.Value=recPat;
        setPatternFields('on');
        applyFilters(); refreshFileList(); refreshIndex();
        setStatus(sprintf('%d file(s) loaded.%s', numel(app.files), filterNote()));
    end

    function createGroupProgrammatic(name, idx)
        c.fileList.Value = idx; c.groupNameF.Value = name; onCreateGroup();
    end
    function a = getAppLive(),      a = app;              end
    function p = sessionPathLive(), p = app.sessionPath;  end
    function v = indexLive(),       v = app.vars;         end

    function setStateProgrammatic(varargin)
        %setStateProgrammatic  The cascade, driven from code.  The keys are the
        %   cascade's own steps and they are applied IN CASCADE ORDER within one call,
        %   so setState('Kind','Vector','Family','Vasomotion',...) behaves exactly like
        %   choosing them one after the other in the window.
        order = {'Files','Product','Kind','Family','Signal','Variable','Plot','Frame', ...
                 'Where','Interval','Units','Baseline','XAxis','Colour','Stat', ...
                 'Title','Legend','LegendPattern','XLabel','YLabel','Cmap','FontSize'};
        args = reshape(varargin,2,[])';
        [~,ord] = ismember(args(:,1), order);
        ord(ord==0) = numel(order)+1;
        [~,ix] = sort(ord);
        args = args(ix,:);
        for k = 1:size(args,1)
            key = args{k,1}; val = args{k,2};
            switch key
                case 'Files',         c.fileList.Value=val; onFileSelection();
                case 'Product',       c.product.Value=val;  onProduct();
                case 'Kind',          setPick(c.kind,val);     onKind();
                case 'Family',        setPick(c.family,val);   onFamily();
                case 'Signal',        setPick(c.signal,val);   onSignal();
                case 'Variable',      setPick(c.variable,val); onVariable();
                case 'Plot',          setPick(c.plot,val);     onPlot();
                case 'Frame',         setFrame(val);
                case 'Where',         c.selList.Items=union(c.selList.Items,cellstr(val),'stable');
                                      c.selList.Value=val;     onWhere();
                case 'Interval'
                    if ~ismember(val,c.interval.Items), c.interval.Items=[c.interval.Items {char(val)}]; end
                    c.interval.Value=val;
                case 'Units',         c.units.Value=val;
                case 'Baseline',      c.baseline.Value=logical(val);
                case 'XAxis'
                    if ~ismember(val,c.xaxis.Items), c.xaxis.Items=[c.xaxis.Items {char(val)}]; end
                    c.xaxis.Value=val;
                case 'Colour'
                    if ~ismember(val,c.colour.Items), c.colour.Items=[c.colour.Items {char(val)}]; end
                    c.colour.Value=val;
                case 'Stat',          c.stat.Value=val;
                case 'Title',         c.titleF.Value=val; app.manual.title=~isempty(val);
                case 'Legend',        c.legendChk.Value=logical(val);
                case 'LegendPattern', c.legendF.Value=val;
                case 'XLabel',        c.xlabF.Value=val; app.manual.xlab=~isempty(val);
                case 'YLabel',        c.ylabF.Value=val; app.manual.ylab=~isempty(val);
                case 'Cmap',          c.cmap.Value=val;
                case 'FontSize',      c.fontSize.Value=val;
                otherwise
                    error('guiExplore:setState','Unknown state key ''%s''.',char(key));
            end
        end
    end

    function setPick(ctl, val)
        %setPick  Choose an entry of a cascade step by its KEY or by the text shown -
        %   whichever the caller has to hand - so a test can name a plot id and a user
        %   can read a phrase.
        want = val; if ischar(want)||isstring(want), want = {char(want)}; end
        data = ctl.ItemsData; if isempty(data), data = ctl.Items; end
        pick = {};
        for i = 1:numel(want)
            j = find(strcmp(data, want{i}),1);
            if isempty(j), j = find(strcmp(ctl.Items, want{i}),1); end
            if isempty(j)
                error('guiExplore:setState','%s is not on offer here.', char(want{i}));
            end
            pick{end+1} = data{j}; %#ok<AGROW>
        end
        if isa(ctl,'matlab.ui.control.DropDown'), ctl.Value = pick{1};
        else,                                     ctl.Value = pick;
        end
    end

    function setFrame(v)
        %setFrame  The slider, driven from code.  It refuses on a plot that has
        %   nothing to scrub rather than silently accepting a position that means
        %   nothing - a test that thinks it moved a slider must be told it did not.
        if isempty(app.scrub.name)
            error('guiExplore:setState','This plot has nothing to scrub along.');
        end
        c.scrub.Value = clampScrub(v);
        onScrub(c.scrub);
    end

    function exportTo(filename, dpi, fmt)
        if nargin<2||isempty(dpi), dpi=c.dpi.Value; end
        if nargin<3||isempty(fmt), fmt=c.fmt.Value; end
        tmp = figure('Visible','off','Color','w','Units','pixels', ...
            'Position',[100 100 1100 800]);
        axe = axes(tmp,'Color','w'); hold(axe,'on'); %#ok<LAXES>
        renderInto(axe);
        hold(axe,'off');
        switch fmt
            case 'PNG',          exportgraphics(axe, filename, 'Resolution', dpi, 'BackgroundColor','white');
            case 'TIFF',         exportgraphics(axe, filename, 'Resolution', dpi, 'BackgroundColor','white');
            case 'JPEG',         exportgraphics(axe, filename, 'Resolution', dpi, 'BackgroundColor','white');
            case 'PDF (vector)', exportgraphics(axe, filename, 'ContentType','vector','BackgroundColor','white');
            case 'EPS (vector)', exportgraphics(axe, filename, 'ContentType','vector','BackgroundColor','white');
        end
        delete(tmp);
    end

    function exportVideoTo(filename, fps)
        %exportVideoTo  THE WHOLE SWEEP, one frame per position along the scrubbed
        %   dimension, written through VideoWriter.  The colour limits are computed
        %   ONCE over the entire cube - a per-frame rescale would make a flat
        %   recording look alive - and each frame carries its own coordinate in the
        %   title, so a still lifted out of the video still says where it came from.
        %
        %   VIDEOWRITER APPENDS ITS PROFILE'S EXTENSION.  Handed 'x.avi.partial' it
        %   writes 'x.avi.partial.avi', which has cost this library time before, so it
        %   is written to a temporary name that ALREADY ends in that extension and the
        %   result is moved to the name the caller asked for.
        if nargin<2||isempty(fps), fps = 8; end
        P = onePayload();
        if isempty(P) || ~isstruct(P) || ~isfield(P,'frames') || isempty(P.frames)
            error('guiExplore:noVideo','There is nothing to write: choose a video plot first.');
        end
        V = P.frames; nF = size(V,3);
        lims = cubeLimits(V);
        prof = videoProfile(filename);
        tmpName = [tempname videoExtension(prof)];
        w = VideoWriter(tmpName, prof);
        w.FrameRate = fps;
        open(w);
        % A modest, EVEN frame size: MPEG-4 refuses odd dimensions, and the DPI control
        % is honoured within the range a video can carry rather than literally - 300 dpi
        % on a figure this size is a 2800-pixel frame nobody asked for.
        tmpFig = figure('Visible','off','Color','w','Units','pixels', ...
            'Position',[100 100 900 700]);
        axe = axes(tmpFig,'Color','w');
        res = min(max(round(c.dpi.Value),72),150);
        savedTitle = c.titleF.Value;
        try
            for k = 1:nF
                cla(axe,'reset');
                img = double(V(:,:,k));
                imagesc(axe, img, 'AlphaData', isfinite(img)); axis(axe,'image');
                if ~isempty(lims), clim(axe,lims); end
                colormap(axe, c.cmap.Value);
                cb = colorbar(axe); cb.Label.String = P.clab;
                set(axe,'XTick',[],'YTick',[],'FontSize',c.fontSize.Value);
                app.titleFrame = k;          % the title names THIS frame, not the slider's
                title(axe, currentTitle(), 'Interpreter','none','FontWeight','bold');
                frame = print(tmpFig,'-RGBImage',sprintf('-r%d',res));
                writeVideo(w, evenSized(frame));
                if mod(k,5)==0 || k==nF
                    setStatus(sprintf('Writing video: frame %d of %d...', k, nF));
                end
            end
        catch ME
            app.titleFrame = []; c.titleF.Value = savedTitle;
            try, close(w); catch, end
            delete(tmpFig);
            if isfile(tmpName), delete(tmpName); end
            rethrow(ME);
        end
        app.titleFrame = []; c.titleF.Value = savedTitle;
        close(w); delete(tmpFig);
        movefile(tmpName, filename, 'f');
        setStatus(sprintf('Exported %d frame(s): %s', nF, filename));
    end

%% ==================== small nested helpers ============================= %%
    function setStatus(msg), c.status.Text = msg; drawnow limitrate; end

    function k = currentKind()
        k = c.kind.Value; if ~ischar(k), k = ''; end
    end

    function p = currentPlot()
        p = c.plot.Value;
        if ~ischar(p) || isempty(c.plot.Items), p = ''; end
    end

    function setPatternFields(state)
        %setPatternFields  The three regexp boxes guess labels from file NAMES, so
        %   they are switched off whenever the labels came from somewhere better - a
        %   single picked file, or a session that already resolved them.
        c.includeF.Enable=state; c.groupF.Enable=state; c.recF.Enable=state;
        c.scanBtn.Enable=state;
    end

    function n = nUnique(field)
        if isempty(app.files), n = 0; return; end
        v = {app.files.(field)};
        n = numel(unique(v(~cellfun(@isempty,v))));
    end

    function n = nDim(tg,field)
        n = numel(unique({tg.(field)}));
    end

    function [statFun,loFun,hiFun,name] = statFuns()
        % Centre + lower/upper band generators chosen from the "Centre/error" box.
        switch c.stat.Value
            case 'Mean +/- SEM'
                name='Mean +/- SEM';
                statFun=@(Y,d) mean(Y,d,'omitnan');
                sem=@(Y,d) std(Y,0,d,'omitnan')./sqrt(max(sum(isfinite(Y),d),1));
                loFun=@(Y,d) statFun(Y,d)-sem(Y,d);
                hiFun=@(Y,d) statFun(Y,d)+sem(Y,d);
            case 'Median +/- IQR'
                name='Median +/- IQR';
                statFun=@(Y,d) median(Y,d,'omitnan');
                loFun=@(Y,d) prctileDim(Y,25,d);
                hiFun=@(Y,d) prctileDim(Y,75,d);
            otherwise
                name='Mean +/- SD';
                statFun=@(Y,d) mean(Y,d,'omitnan');
                sd=@(Y,d) std(Y,0,d,'omitnan');
                loFun=@(Y,d) statFun(Y,d)-sd(Y,d);
                hiFun=@(Y,d) statFun(Y,d)+sd(Y,d);
        end
    end







    function e = makeEntry(path,grp,rec)
        e = makeEntryStruct(path,'','');
        e.group = grp; e.rec = rec; e.recnum = recToNum(rec);
    end
end

%% ========================= LOCAL HELPERS ================================ %%
function s = section(parent, titleText, nBodyRows, rowH)
% A titled sub-panel holding a small 2-column grid of controls.  A row height may be
% given explicitly, which the cascade needs: a LIST in a 'fit' row collapses to one
% line, and a one-line list is not a list.
if nargin<4 || isempty(rowH), rowH = repmat({'fit'},1,nBodyRows); end
p = uipanel(parent,'Title',titleText,'FontWeight','bold','BackgroundColor','w', ...
    'ForegroundColor',[0.15 0.25 0.45]);
s = uigridlayout(p,[nBodyRows 2],'RowHeight',rowH, ...
    'ColumnWidth',{'fit','1x'},'RowSpacing',5,'ColumnSpacing',6, ...
    'Padding',[6 6 6 6],'BackgroundColor','w');
end

function f = labelledField(g,row,name,val)
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
f=uieditfield(g,'text','Value',val); f.Layout.Row=row; f.Layout.Column=2;
end

function [lb,d] = labelledDropH(g,row,name,items,cb)
%labelledDropH  A dropdown that hands back its LABEL as well, so a cascade step with
%   nothing to choose can hide both halves of itself rather than showing a label over
%   an empty box.
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
d=uidropdown(g,'Items',items); d.Layout.Row=row; d.Layout.Column=2;
if ~isempty(cb), d.ValueChangedFcn=cb; end
end

function [lb,f] = labelledListH(g,row,name,cb)
%labelledListH  The same, for the two MULTI-SELECT steps.
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
lb.VerticalAlignment='top';
f=uilistbox(g,'Items',{},'Multiselect','on'); f.Layout.Row=row; f.Layout.Column=2;
if ~isempty(cb), f.ValueChangedFcn=cb; end
end

function t = tipInclude()
% Hover helper for the Include field.
t = sprintf(['Regular expression matched against each file NAME (not the folder).\n' ...
    'Only files whose name matches are kept.  Examples:\n' ...
    '   _r\\.mat$          all results files (default)\n' ...
    '   _c_BFI_r\\.mat$    cardiac BFI results only\n' ...
    '   PSY0[678]         only animals PSY06 / 07 / 08\n' ...
    '   (i2|k1)           only i2 or k1 recordings\n' ...
    'Tip: \\.  = a literal dot,   $ = end of name,   (a|b) = a or b.']);
end
function t = tipGroup()
% Hover helper for the Group field.
t = sprintf(['Regexp whose MATCH labels each file''s experimental group.\n' ...
    'Files sharing a match form one group (mean +/- SD, points = files).\n' ...
    'Anchor it so you don''t accidentally match the date digits!  Examples:\n' ...
    '   PSY\\d+            group by animal id (PSY06, PSY07, ...)\n' ...
    '   WT|KO             genotype taken from the name\n' ...
    '   Ctrl|Stroke       condition taken from the name\n' ...
    'Leave empty for a single group.  You can also select files in the\n' ...
    'list and press "Create group" to set groups by hand.']);
end
function t = tipRec()
% Hover helper for the Rec.index field.
t = sprintf(['Regexp whose MATCH is the recording index - a second, independent\n' ...
    'grouping axis (stratifier).  If the match contains a number, indices are\n' ...
    'ordered ASCENDING by it.  Examples:\n' ...
    '   [aik]\\d           a condition+repeat token, e.g. i2, k1, a2\n' ...
    '   run(\\d+)          the number after "run"\n' ...
    '   \\d+(?=_c_)        the digits just before _c_\n' ...
    'Leave empty to treat every file as the same index.']);
end
function d = labelledDrop(g,row,name,items,cb)
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
d=uidropdown(g,'Items',items); d.Layout.Row=row; d.Layout.Column=2;
if ~isempty(cb), d.ValueChangedFcn=cb; end
end
function sp = labelledSpin(g,row,name,val)
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
sp=uispinner(g,'Value',val,'Limits',[50 1200],'Step',50); sp.Layout.Row=row; sp.Layout.Column=2;
end

function entries = buildFileEntries(root, includePat, groupPat, recPat)
% Recursively list *.mat under root, keep those whose name matches includePat,
% then tag each with an experimental group and a recording index via regexp.
if isempty(includePat), includePat='_r\.mat$'; end
d = dir(fullfile(root,'**','*.mat')); d = d(~[d.isdir]);
entries = emptyEntries();
for i=1:numel(d)
    if isempty(regexp(d(i).name, includePat, 'once')), continue; end
    entries(end+1) = makeEntryStruct(fullfile(d(i).folder,d(i).name), groupPat, recPat); %#ok<AGROW>
end
entries = sortEntries(entries);
end

function e = emptyEntries()
%emptyEntries  The 0x0 shape of the file list.  Declared ONCE so every input route
%   - folder scan, single file, programmatic paths, session - produces the same
%   fields and the struct arrays concatenate.  animal and type are empty unless a
%   session filled them.
e = struct('path',{},'name',{},'group',{},'rec',{},'recnum',{}, ...
    'animal',{},'type',{});
end

function e = makeEntryStruct(path, groupPat, recPat)
[~,nm,ex]=fileparts(path); name=[nm ex];
grp = regexpMatch(name, groupPat, 'all');
rec = regexpMatch(name, recPat,  '1');
e = struct('path',path,'name',name,'group',grp,'rec',rec,'recnum',recToNum(rec), ...
    'animal','','type','');
end

function v = emptyVars()
%emptyVars  The 0x0 shape of a CASCADE ROW - one (family, signal, variable) the
%   loaded files offer.  Declared once, because the rows built from the index and the
%   rows built from the recording's source have to concatenate.
v = struct('key',{},'family',{},'signal',{},'signalLabel',{},'varId',{}, ...
    'varLabel',{},'leaf',{},'unit',{},'dims',{},'sig',{},'kinds',{},'plots',{}, ...
    'suspect',{},'srcView',{},'srcInterval',{},'nFiles',{});
end

function v = sortVars(v)
%sortVars  Menu order: the families in the order a user thinks about them, and a
%   SUSPECT leaf last - offered, but not offered first.
if isempty(v), return; end
fam = {'Flow','Pulsatility','Vasomotion','Diameter','Propagation','Maps','Metrics','Other'};
[~,fi] = ismember({v.family}, fam); fi(fi==0) = numel(fam)+1;
% one element per row even when a row left .suspect unset - a shorter column here
% would not sort wrongly, it would fail to concatenate at all
sus = arrayfun(@(r) double(~isempty(r.suspect) && r.suspect), v);
[~,ord] = sortrows([sus(:), fi(:), (1:numel(v))'],[1 2 3]);
v = v(ord);
end

function o = emptyObs()
o = struct('x',{},'y',{},'sel',{},'group',{},'rec',{},'animal',{},'type',{}, ...
    'interval',{},'file',{},'var',{},'signal',{});
end

function t = emptyTags()
t = struct('sel',{},'group',{},'rec',{},'animal',{},'type',{},'interval',{}, ...
    'file',{},'var',{},'signal',{});
end

function u = uniqueStable(cs)
if isempty(cs), u = {}; return; end
u = unique(cs,'stable');
end

function k = selectedKeys(ctl)
%selectedKeys  What a multi-select step has picked, as a cellstr of KEYS.
k = ctl.Value;
if isempty(k), k = {}; return; end
if ischar(k), k = {k}; end
end

function setStep(lb, ctl, labels, keys)
%setStep  Fill one cascade step, and say what its state means.  A step with NOTHING
%   to offer is switched OFF rather than left showing the previous file's answer -
%   that staleness is the bug this redesign exists to fix.  A step with exactly ONE
%   answer fills itself and HIDES: a myograph with one measure should not be asked
%   which measure.
keep = ctl.Value;
if isempty(labels)
    ctl.Items = {}; ctl.ItemsData = {}; ctl.Enable = 'off';
    setVisible(lb, ctl, true);
    return
end
ctl.Items = labels; ctl.ItemsData = keys; ctl.Enable = 'on';
if ischar(keep) && any(strcmp(keys,keep)), ctl.Value = keep; else, ctl.Value = keys{1}; end
setVisible(lb, ctl, numel(labels)>1);
end

function setStepList(lb, ctl, labels, keys)
%setStepList  The same, for a multi-select step.  A selection that survives the
%   rebuild is kept, so narrowing the step above does not throw away a choice the
%   user already made.
keep = selectedKeys(ctl);
if isempty(labels)
    ctl.Items = {}; ctl.ItemsData = {}; ctl.Enable = 'off';
    setVisible(lb, ctl, true);
    return
end
ctl.Items = labels; ctl.ItemsData = keys; ctl.Enable = 'on';
keep = keep(ismember(keep, keys));
if isempty(keep), keep = keys(1); end
ctl.Value = keep;
setVisible(lb, ctl, numel(labels)>1);
end

function setVisible(lb, ctl, tf)
%setVisible  Show or hide one cascade step - and COLLAPSE ITS ROW when it hides.  The
%   two list steps sit in rows of a fixed height, so hiding the control alone would
%   leave a hole the height of a list where the step used to be.
s = 'off'; if tf, s = 'on'; end
lb.Visible = s; ctl.Visible = s;
g = ctl.Parent;
if ~isa(g,'matlab.ui.container.GridLayout'), return; end
r = ctl.Layout.Row;
if ~isscalar(r) || r<1 || r>numel(g.RowHeight), return; end
if isempty(lb.UserData), lb.UserData = g.RowHeight{r}; end   % the height it wants
if tf, g.RowHeight{r} = lb.UserData; else, g.RowHeight{r} = 0; end
end

function setDropKeeping(ctl, items)
%setDropKeeping  Rebuild a plain dropdown's list, keeping the current choice when it
%   survives.
keep = ctl.Value;
ctl.Items = items;
if any(strcmp(items,keep)), ctl.Value = keep; else, ctl.Value = items{1}; end
end

function rp = resolveResultFiles(pth)
%resolveResultFiles  The *_r.mat RESULTS a session recording owns.  A session lists
%   the RECORDINGS the workbench curated (often raw .rls), and this tool plots the
%   RESULTS - so one recording expands to every result product it carries ('_t_BFI_r',
%   '_c_BFI_r', ...).  Identity comes from wbFileModel, so a 'Roi2_' crop of the same
%   stem stays a different recording.  Mirrors guiExport>resolveProducts on the '_r'
%   side of the triplet.
rp = {};
pth = char(pth);
if isempty(pth), return; end
if endsWith(pth,'_r.mat')
    if isfile(pth), rp = {pth}; end
    return
end
if endsWith(pth,'_d.mat') && isfile(strrep(pth,'_d.mat','_r.mat'))
    rp = {strrep(pth,'_d.mat','_r.mat')};
    return
end
m = wbFileModel(pth);
if isempty(m.folder) || ~isfolder(m.folder), return; end
d = dir(fullfile(m.folder,[m.roiPrefix m.stem '*_r.mat']));
for k = 1:numel(d)
    p = fullfile(d(k).folder,d(k).name);
    cm = wbFileModel(p);
    if ~strcmp(cm.identity, m.identity), continue; end
    rp{end+1} = p; %#ok<AGROW>
end
end

%% ==================== THE PLOTS, AS A READER SEES THEM ================== %%
function s = plotLabel(plotId)
%plotLabel  The plot ids of exPlotRules, spelled for a biologist.  Two of them are
%   this tool's own - the per-line overlay and the detected walls read the
%   recording's SOURCE rather than its results, so they are not in the rule table.
switch char(plotId)
    case 'box',                s = 'Box plot';
    case 'bar',                s = 'Bar';
    case 'curve.time',         s = 'Against time';
    case 'curve.f',            s = 'Against frequency';
    case 'curve.pct',          s = 'Against percentile';
    case 'curve.harmonic',     s = 'Against harmonic';
    case 'curve.position',     s = 'Along the vessel';
    case 'curves.family',      s = 'A family of curves';
    case 'image',              s = 'Map';
    case 'image.fromSegments', s = 'Map, painted from the segments';
    case 'image.timeAverage',  s = 'Map, averaged over time';
    case 'image.frame',        s = 'Map, one frame';
    case 'heat.ft',            s = 'Time against frequency';
    case 'heat.fpct',          s = 'Percentile against frequency';
    case 'heat.positionTime',  s = 'Position against time';
    case 'video',              s = 'Video';
    case 'lines.overlay',      s = 'Every line';
    case 'image.walls',        s = 'Detected walls';
    otherwise,                 s = char(plotId);
end
end

function k = kindOfPlot(plotId)
%kindOfPlot  Which kind a plot appears under.  Delegated to exPlotRules for everything
%   it knows about, so the two answers cannot drift apart.
switch char(plotId)
    case 'lines.overlay', k = 'Vector';
    case 'image.walls',   k = 'Image';
    otherwise,            k = exPlotRules('kindof', plotId);
end
end

function s = seriesLabel(names, i)
if iscell(names) && numel(names)>=i && ~isempty(names{i}), s = names{i};
else, s = sprintf('%d',i);
end
end

%% ============ THE SCRUBBED DIMENSION, AS A READER SEES IT =============== %%
% A slice of a cube is only meaningful if the reader can say WHERE it was taken, so
% every caption below is written in the recording's own coordinate - never in the
% index that happens to address it.

function s = scrubText(spec, i, withCount)
%scrubText  The slider's caption.  'frequency = 0.12 Hz' is a place in the
%   recording; 'index 7' is a place in an array, and the two are not the same claim.
%   The count belongs beside the slider and not in a title, which is why it is
%   optional rather than baked in.
if nargin<3, withCount = true; end
s = '';
if isempty(spec.name), return; end
v = [];
if ~isempty(spec.values) && i>=1 && i<=numel(spec.values), v = spec.values(i); end
s = coordText(spec.label, v);
if withCount && spec.n>1, s = sprintf('%s   (%d of %d)', s, i, spec.n); end
end

function s = frameText(P, k)
%frameText  Where one frame of a video sits on the scrubbed axis.
s = '';
if ~isfield(P,'fvals') || isempty(P.fvals) || k<1 || k>numel(P.fvals), return; end
s = coordText(P.flab, P.fvals(k));
end

function s = coordText(label, v)
%coordText  'frequency = 0.12 Hz', out of the axis label 'frequency (Hz)' and a
%   number - the unit is where the label already keeps it, in brackets at the end.
lab = strtrim(char(label)); unit = '';
t = regexp(lab,'^(.*?)\s*\(([^()]*)\)\s*$','tokens','once');
if ~isempty(t), lab = strtrim(t{1}); unit = strtrim(t{2}); end
if isempty(lab), lab = 'frame'; end
if isempty(v), s = lab; return; end
s = sprintf('%s = %s', lab, num2str(double(v),'%.4g'));
if ~isempty(unit), s = [s ' ' unit]; end
end

function lims = cubeLimits(V)
%cubeLimits  ONE COLOUR SCALE FOR EVERY FRAME, taken over the whole cube.  Rescaling
%   each frame to its own range is what makes a recording that barely changes look
%   alive, and that is a lie a movie tells very convincingly.
lims = [];
v = double(V(:)); v = v(isfinite(v));
if isempty(v), return; end
lims = prctile(v,[2 99]);
if ~(lims(2)>lims(1)), lims = [min(v) max(v)]; end
if ~(lims(2)>lims(1)), lims = []; end
end

function s = blankScrub()
% Nothing to walk.  Built in one place so the empty answer and exFetch's full one
% cannot grow different fields - a struct missing one is an error at render time.
s = struct('name','','label','','values',[],'n',0,'idx',1, ...
           'choices',{{}},'choiceLabels',{{}});
end

function h = scrubRowHeight(tf)
% A hidden strip must not leave a gap under the axes, so the row itself collapses.
if tf, h = 'fit'; else, h = 1; end
end

function p = videoProfile(filename)
%videoProfile  Which VideoWriter profile the requested NAME asks for.  Anything the
%   profiles do not recognise is written as Motion JPEG AVI, which every machine can
%   read - and the file still ends up under the name that was asked for, see below.
[~,~,e] = fileparts(char(filename));
switch lower(e)
    case '.mp4', p = 'MPEG-4';
    case '.mj2', p = 'Motion JPEG 2000';
    otherwise,   p = 'Motion JPEG AVI';
end
end

function e = videoExtension(profile)
%videoExtension  THE EXTENSION VIDEOWRITER WILL APPEND WHATEVER IT IS GIVEN.  Writing
%   to a temporary name that already carries it is what stops 'sweep.partial' from
%   becoming 'sweep.partial.avi' on the way to its real name.
switch char(profile)
    case 'MPEG-4',             e = '.mp4';
    case 'Motion JPEG 2000',   e = '.mj2';
    otherwise,                 e = '.avi';
end
end

function F = evenSized(F)
% MPEG-4 refuses an odd frame dimension, so the last row or column is dropped rather
% than the export failing on the frame size of a figure nobody chose.
if mod(size(F,1),2), F = F(1:end-1,:,:); end
if mod(size(F,2),2), F = F(:,1:end-1,:); end
end

%% ==================== THE BRANCH / PRODUCT FILTER ====================== %%
function k = productKeyOf(pth)
%productKeyOf  Which pipeline and which product one results file is.  The parse is
%   wbProducts' and wbFileModel's, NEVER a prefix test on the name: a derived product
%   SUBSTITUTES the product token rather than extending it ('_t_K' becomes '_t_BFI'),
%   so there is no prefix for such a test to find.  See wbProducts' header.
m = wbFileModel(pth);
st = m.stage;
try
    s2 = wbProducts('stageOf', pth);
    if ~isempty(s2) && ~strcmp(s2, wbProducts('unknown')), st = s2; end
catch
end
k = [wbFileModel('branch', st, m.product) '|' m.product];
end

function s = productLabelOf(key)
p = strsplit(char(key),'|');
br = p{1}; pr = ''; if numel(p)>1, pr = p{2}; end
if isempty(br) && isempty(pr), s = '(unlabelled)'; return; end
if isempty(br), s = pr; return; end
if isempty(pr), s = br; return; end
s = [br ' / ' pr];
end

function r = branchRank(br)
%branchRank  The pipeline order, for breaking a tie on how many files each holds.
switch char(br)
    case 'contrast', r = 1;
    case 'cardiac',  r = 2;
    case 'epoch',    r = 3;
    case 'bolus',    r = 4;
    case 'myograph', r = 5;
    otherwise,       r = 6;
end
end

function r = productRank(pr)
%productRank  How far down the pipeline a product is, so the DEFAULT is the finished
%   article rather than the intermediate it was computed from.
switch char(pr)
    case 'MYO', r = 4;
    case 'BFI', r = 3;
    case 'K',   r = 2;
    case {'I','g'}, r = 1;
    otherwise,  r = 0;
end
end

function L = emptyLegendTag()
%emptyLegendTag  The blank legend record: one field per token legendName expands.
L = struct('sel','','group','','rec','','animal','','type','','interval','','file','', ...
    'var','','signal','');
end

function s = labelOr(v, dflt)
s = charOf(v);
if isempty(s), s = dflt; end
end

function s = charOf(v)
if ischar(v), s = v; elseif isstring(v), s = char(v);
elseif isnumeric(v) || islogical(v), s = num2str(v); else, s = ''; end
end

function s = shortName(p)
[~,n,e] = fileparts(char(p)); s = [n e];
end

function m = regexpMatch(str, pat, dflt)
if isempty(pat), m=dflt; return; end
tok = regexp(str, pat, 'match', 'once');
if isempty(tok), m=dflt; else, m=tok; end
end

function n = recToNum(rec)
t = regexp(rec,'\d+','match','once');
if isempty(t), n=Inf; else, n=str2double(t); end
end

function entries = sortEntries(entries)
if isempty(entries), return; end
[~,ord] = sortrows([ [entries.recnum]', (1:numel(entries))' ], [1 2]);
% keep alpha order within equal numbers, and push non-numeric (Inf) to the end
recnum = [entries.recnum]';
[~,ord2] = sort(recnum,'ascend');
entries = entries(ord2);
end


function s = onOff(tf), if tf, s = 'on'; else, s = 'off'; end, end

%% ======================================================================= %%








function i = firstIndex(key,uk)
i=zeros(1,numel(uk));
for m=1:numel(uk), i(m)=find(strcmp(key,uk{m}),1); end
end

function xl = xLabelFor(name)
switch name
    case 'rec',      xl='recording index';
    case 'group',    xl='group';
    case 'animal',   xl='animal';
    case 'type',     xl='recording type';
    case 'interval', xl='interval';
    otherwise,       xl='selection';
end
end

function cats = orderedCats(vals)
u = unique(vals,'stable');
cats = u(orderNumeric(u));
end
function ord = orderNumeric(u)
% Order categories by an embedded number ascending; categories without a number
% keep their first-appearance order (so a manual selection list stays as picked).
num = nan(1,numel(u));
for i=1:numel(u), t=regexp(u{i},'\d+','match','once'); if ~isempty(t), num(i)=str2double(t); end, end
if all(isnan(num))
    ord = 1:numel(u);                       % no numbers anywhere -> stable order
else
    [~,ord]=sortrows([num(:), (1:numel(u))'],[1 2]);
    nanmask=isnan(num(ord)); ord=[ord(~nanmask); ord(nanmask)];   % numbered first, then rest
end
ord=ord(:)';
end

function cols = seriesColours(legOrNames)
% Distinct, colour-blind-aware colours. Semantic tints (artery=red, vein=blue,
% parenchyma=green, ...) are used ONLY when they stay unique; if two series would
% map to the same tint (e.g. arteries across several groups) the whole set falls
% back to a distinct qualitative palette so no two series ever share a colour.
n = numel(legOrNames);
sem = nan(n,3);
for i=1:n, sem(i,:) = semanticColour(seriesSelText(legOrNames{i})); end
hasSem = ~any(isnan(sem),2);
S = sem(hasSem,:);
clash = size(S,1)>1 && size(unique(round(S*1000),'rows'),1) < size(S,1);
if clash
    cols = distinctPalette(n);                 % semantic tints would repeat -> palette
    return
end
cols = sem;                                    % keep the (unique) semantic tints
pal  = qualPalette(); used = S; pk = 1;
for j = find(~hasSem)'                          % fill the rest from the palette, no repeats
    while pk<=size(pal,1) && ~isempty(used) && any(all(abs(used-pal(pk,:))<0.03,2)), pk=pk+1; end
    if pk>size(pal,1), cols(j,:)=hsvColour(j); else, cols(j,:)=pal(pk,:); used=[used;pal(pk,:)]; pk=pk+1; end %#ok<AGROW>
end
end
function cols = distinctPalette(n)
pal = qualPalette(); cols = zeros(n,3);
for i=1:n
    if i<=size(pal,1), cols(i,:)=pal(i,:); else, cols(i,:)=hsvColour(i); end
end
end
function s = seriesSelText(L)
if ischar(L)||isstring(L), s=char(L); return; end
if isfield(L,'sel')&&~isempty(L.sel), s=L.sel; elseif isfield(L,'dyn'), s=L.dyn; else, s=''; end
end
function c = semanticColour(s)
s=lower(s);
if     contains(s,'arter'),      c=[0.80 0.16 0.16];
elseif contains(s,'vein'),       c=[0.16 0.34 0.63];
elseif contains(s,'parenchyma'), c=[0.20 0.55 0.30];
elseif contains(s,'uncert'),     c=[0.50 0.50 0.50];
elseif contains(s,'all vessel')||contains(s,'all dvs'), c=[0.45 0.20 0.60];
else,  c=[NaN NaN NaN];
end
end
function pal = qualPalette()
pal = [0.12 0.47 0.71; 0.85 0.37 0.01; 0.17 0.63 0.17; 0.84 0.15 0.16; ...
       0.58 0.40 0.74; 0.55 0.34 0.29; 0.89 0.47 0.76; 0.45 0.45 0.45; ...
       0.74 0.74 0.13; 0.09 0.75 0.81; 0.20 0.30 0.60; 0.95 0.65 0.10];
end
function c = hsvColour(i,n)
c = hsv2rgb([mod((i-1)*0.61803,1) 0.65 0.85]);
if nargin<2, n=1; end
end

function p = prctileDim(Y,q,dim)
% Percentile along one dimension, NaN-tolerant, returned with that dim collapsed.
if dim==2, p = nan(size(Y,1),1); for i=1:size(Y,1), p(i)=prctile(Y(i,:),q); end
else,      p = nan(1,size(Y,2)); for i=1:size(Y,2), p(i)=prctile(Y(:,i),q); end
end
end

function xg = commonGrid(obs)
% Shared x-grid for pooling curves. Identical grids are used as-is; otherwise a
% linear grid over the overlapping range with the median sample count.
X = {obs.x};
if all(cellfun(@(x)isequal(x,X{1}),X)), xg=X{1}(:); return; end
lo = max(cellfun(@min,X)); hi = min(cellfun(@max,X));
if ~(hi>lo), lo=min(cellfun(@min,X)); hi=max(cellfun(@max,X)); end
n = round(median(cellfun(@numel,X)));
xg = linspace(lo,hi,max(n,2))';
end
function Y = resampleObs(obs, xg)
Y = nan(numel(xg),numel(obs));
for i=1:numel(obs)
    x=obs(i).x(:); y=obs(i).y(:);
    if isequal(x,xg), Y(:,i)=y;
    else, Y(:,i)=interp1(x,y,xg,'linear',NaN);
    end
end
end

function s = prettySel(sel)
if startsWith(sel,'Label:'), s=strtrim(sel(7:end));
else, s=regexprep(sel,'\s*\((lumen|DVS)\)','');
end
end
function s = pad0(x), if isempty(x), s=''; else, s=char(string(x)); end, end

function [sess, vis] = parseArgs(args)
% guiExplore(sessionPath) / guiExplore('Visible',v) / both.  A LEADING char that is
% not an option name is the session path - the shape guiWorkbench hands over and the
% shape guiExport takes, so the two tools are called the same way.
sess = ''; vis = 'on'; first = 1;
if ~isempty(args) && (ischar(args{1}) || isstring(args{1})) && ~isOptionName(args{1})
    sess = char(args{1}); first = 2;
end
for i = first:2:numel(args)-1
    if (ischar(args{i}) || isstring(args{i})) && strcmpi(args{i},'Visible')
        vis = char(string(args{i+1}));
    end
end
end
function tf = isOptionName(x)
tf = any(strcmpi(char(string(x)), {'Visible'}));
end
