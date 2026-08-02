% guiExplore - Interactive GUI to explore & export publication figures from results.
%
% WHAT THIS TOOL DOES
%   A single-window app for browsing the *_r.mat RESULTS produced by this library
%   (runSegmentation -> runPulsatility / runVasomotion / runCTTH / setVesselTypes,
%   and the myograph chain runMyographDiameter -> runMyographPropagation ->
%   runMyographVasomotion) and turning them into clean, publication-ready plots. It
%   works on ONE file or on MANY files at once, and it always matches the plot type
%   to the kind of data:
%
%     * time series (sData / dvsData / dvsDiameter)      -> mean +/- error band vs time
%     * scalar metrics (sMetrics / dvsMetrics columns)   -> grouped box plots + points
%     * vasomotion spectra (results.vasomotion.*)        -> amplitude vs frequency
%     * vasomotion time-frequency (...spectrum.amp)      -> 2-D spectrogram
%     * vasomotion percentile spectra (fVectors.VB.ampMeanPct) -> spectra split by
%                                                          envelope-amplitude percentile bin
%     * vasomotion amplitude percentiles (scalars.VB/CB.ampPct) -> band amplitude vs
%                                                          percentile level (VB & CB)
%     * vasomotion markers (scalars.VB/CB.<marker>)      -> grouped box plots + points
%     * diameter traces (myograph)                       -> mean +/- error band vs time,
%                                                          in pixels or % of baseline
%     * the individual line diameters (myograph)         -> every line's own trace
%     * diameter statistics / propagation (myograph)     -> grouped box plots + points
%     * the per-line propagation lag (myograph)          -> lag vs position along vessel
%     * the diameter along the vessel (myograph)         -> position vs time map
%     * the detected walls (myograph)                    -> one video frame + overlay
%     * image maps (imgBFI / pulsatility.ppx maps / mapType) -> axis-image display
%
%   THE VASOMOTION MARKER FAMILY SERVES BOTH KINDS OF RECORDING.  scalars.<band>.*
%   holds the ~29 band markers (ampMean, fCentMean, durFlareMean, ...) in the same
%   place in a speckle tree and in a myograph one, so one box-plot family compares
%   them across files for either - by group, animal, type or interval like any other
%   scalar.  Percentile leaves are not offered here: they are a vector per unit, and
%   the amplitude-percentile family draws them properly.
%
%   MYOGRAPH RECORDINGS ARE PLOTTED BY THE SAME MACHINERY.  A myograph recording
%   holds its results inside the WINDOWS it was analysed in, and
%   results.intervals(k).vasomotion.<signal> is the SAME tree shape as
%   results.vasomotion.<signal> - so the spectra, the percentile spectra, the
%   spectrogram and the band amplitude percentiles are the plots that already
%   existed, pointed at an interval's sub-tree.  What is new is the axis: the
%   INTERVAL joins group / recording index / animal / type as a fifth, independent
%   way to organise a comparison, and it is usually the one a myograph protocol
%   wants on the x-axis (baseline vs drug vs washout).  The signal a myograph
%   comparison is about is the analysed diameter MEASURE, or the recorded CHANNEL -
%   which is what the selection list offers for those files, in place of the vessel
%   categories they do not have.
%
%   IT IS ALSO THE ONLY PLACE A MYOGRAPH FIGURE COMES FROM.  No myograph step
%   writes a report page, by design, so the checks a page would have carried live
%   here: the diameter map, which shows whether detection held along the whole
%   vessel; the detected walls over a frame of the recording, which shows whether it
%   found the right edges; and the individual line traces, which show whether the
%   lines moved together - the assumption a propagation speed is fitted under.  All
%   three are about ONE recording, so they sit with the single-recording views rather
%   than with the comparisons.  The per-line propagation LAG is the fourth: it is the
%   evidence behind the speed, and it compares across files like any other curve.
%
% WHEN TO USE IT
%   After you have processed recordings to *_BFI_r.mat (or *_I_r.mat / *_MYO_r.mat)
%   results and want to (a) look at single-recording detail, or (b) compare
%   experimental groups and/or a sequence of recordings across many files, then
%   export the figure at a chosen resolution for a paper or a talk.
%
%   It is STANDALONE. It does not need - and never loads - the Processing
%   Workbench: the hand-off between them is the SESSION file (wbSession), read
%   through wbSession('read',...), which carries the curated file list together
%   with each file's ANIMAL, recording TYPE and experimental GROUP, and a record of
%   which steps actually ran. Its peer on the export side is guiExport.
%
% HOW TO USE IT
%   1. Get files in, three ways:
%        - "Choose Folder / File" -> a FILE: everything is plotted for that single
%          recording; a FOLDER: searched recursively, then labelled by the three
%          regular-expression boxes:
%            (a) Include  : keep only files whose NAME matches this regexp
%                           (default '_r\.mat$' = results files).
%            (b) Group    : a regexp whose MATCH labels the experimental group of
%                           each file (e.g. 'a2|i2|k\d' for anaesthesia condition,
%                           or 'PSY\d+' for animal). Files sharing a match are one
%                           group and are drawn as mean +/- std with each file a point.
%            (c) Rec.index: a regexp whose MATCH is the recording index; files are
%                           stratified by it (like groups, but a separate axis).
%                           If the match contains a number, indices are ordered
%                           ASCENDING by that number. Group and Rec.index are
%                           independent - either can be the x-axis or the colour.
%          Click "Scan / apply" to (re)build the file list. Hover any of the three
%          regexp boxes for a popup with worked pattern examples.
%        - "Load session..." -> a workbench session. Its curated recordings are
%          resolved to their *_r.mat RESULTS, and every file arrives with its
%          ANIMAL, recording TYPE, experimental GROUP and recording INDEX already
%          decided, so nothing is re-derived from the file name: the session's
%          expGroup IS this tool's "group". The regexp boxes are switched off
%          because there is nothing left for them to guess, and "Scan / apply"
%          re-reads the session. The session also says which steps actually ran, so
%          a result tree that was never computed is not offered as a plot.
%      Prefer to group by hand? Select files in the list (Ctrl/Shift-click), type a
%      name in "New group" and press "Create group from selected files". Hand-made
%      groups override BOTH the Group pattern and the session, and survive a re-scan.
%   2. Choose WHAT to plot with the dropdowns (data type -> variable), pick one or
%      more SELECTIONS (arteries / veins / parenchyma / all / a named label; for a
%      myograph recording, the analysed diameters or channels), tweak titles /
%      labels / legend, then "Plot". Labels that span several segments are
%      area-weighted into one trace/value automatically. "Organise by" stratifies on
%      any of the FIVE independent axes - group, recording index, animal, type,
%      interval - which do not nest: an animal may span groups and a group may span
%      animals. "Interval" narrows every myograph plot to one analysed window, and
%      is switched off when nothing loaded has any.
%   3. Set DPI + format and click "Export image".
%
% PROGRAMMATIC USE
%   getappdata(fig,'exploreAPI') exposes the same logic as a struct of handles -
%   .loadPaths .loadSession .setState .createGroup .render .export .getApp - so the
%   tool is fully headless-testable (claude-tests/Workbench/testExploreTool.m),
%   exactly like guiExport's API.
%
% NOTES
%   - Only category-5 (lumen) segments are used for Artery/Vein/All-vessel selections
%     and, by default, for named vessel labels; parenchyma uses category 1. See the
%     library README and getPixelCategories for the category definitions.
%   - Very large results files (e.g. temporal vasomotion files with per-pixel maps)
%     are read field-by-field (HDF5) so time-series / images / spectra can still be
%     shown without loading the whole struct; per-segment SELECTION metadata (type /
%     label) needs a full load, so for such files only the "All ..." selection is
%     offered. A confirmation is shown before a multi-GB full load.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar tool)
%   "Write a self-contained MATLAB uifigure app that loads one or many *_r.mat files,
%    each holding a struct 'results' with a label map results.sMap, a per-segment
%    table results.sMetrics (columns include idx, category [5=lumen,1=parenchyma],
%    area, type ['Artery'/'Vein'/'Uncertain'], label [free text], BFI, psPI, ...),
%    a time vector results.time, per-segment time series results.sData [nTime x nSeg]
%    aligned to the table rows, optional dynamic-vessel fields dvsData/dvsDiameter/
%    dvsMetrics, optional vasomotion results.vasomotion.<sig>.{fVectors.ampMean/ampStd,
%    fVectors.VB.ampMeanPct, scalars.VB/CB.ampPct, spectrum.amp/phase} + root axes
%    f/timeDWT/pctCenters, and maps imgBFI/mapType plus the relocated
%    per-pixel pulsatility maps results.pulsatility.ppx.scalars.*. Provide a folder/file
%    picker with three regexp boxes (include-filter, experimental-group, recording-
%    index sorted ascending by any number), dropdowns to choose data type + variable
%    + one-or-more selections (arteries/veins/parenchyma/all/each label, labels area-
%    weighted across their segments), auto-picked plot types (time-series mean+/-band,
%    grouped boxplots with individual points, amplitude-vs-frequency spectra, 2-D
%    spectrograms, axis-image maps), colour-blind-safe non-clashing colours, an
%    editable title, a single editable legend pattern with tokens, tight/padded axes
%    per plot type, a white background, and high-resolution export via exportgraphics."
%
% Syntax:
%    guiExplore                         % open the tool (folder / file route)
%    guiExplore(sessionPath)            % open it ON a workbench session
%    h = guiExplore('Visible','off')    % headless (programmatic drive / tests)
%
% See also: wbSession, guiExport, guiWorkbench, exportToExcel, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

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
app.labels      = {};        % union of vessel/ROI labels across the file list
app.groupOverride = containers.Map('KeyType','char','ValueType','char'); % hand-made groups (path->name)
app.manual      = struct('title',false,'xlab',false,'ylab',false); % user-edited?
app.sessionPath = '';        % the session this file list came from ('' = none)
% What the loaded myograph recordings offer, unioned over the file list: the names
% of their analysed WINDOWS and of the SIGNALS analysed inside them (the diameter
% measures of a pressure myograph, the channels of a wire one).  Empty when nothing
% loaded is a myograph recording, which is what switches those menus off.
app.myo         = struct('intervals',{{}},'signals',{{}},'only',false);

CAT = struct('PARENCHYMA',1,'UNSEG',2,'WALL',3,'INNER',4,'LUMEN',5);

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

% left: the single large plotting axes
axPanel = uipanel(outer,'BackgroundColor','w','BorderType','none');
axPanel.Layout.Column = 1;
axL = uigridlayout(axPanel,[1 1],'Padding',[2 2 2 2],'BackgroundColor','w');
ax = uiaxes(axL); ax.Color = 'w'; ax.Box = 'on';
title(ax,'Choose a folder or a file to begin'); app.ax = ax;

% right: scrollable control column
ctrlPanel = uipanel(outer,'BackgroundColor','w','BorderType','none','Scrollable','on');
ctrlPanel.Layout.Column = 2;
C = uigridlayout(ctrlPanel,[1 1],'RowHeight',{'fit'},'Padding',[2 2 2 2], ...
    'BackgroundColor','w','Scrollable','on');   % the grid owns the scroll, not the panel
stack = uigridlayout(C,[8 1],'RowHeight',repmat({'fit'},1,8), ...
    'RowSpacing',8,'Padding',[0 0 0 0],'BackgroundColor','w');

c = struct();   % control handles

% --- section 1: data source ---
s1 = section(stack,'1 - Data source',12);
c.sourceBtn = uibutton(s1,'Text','Choose Folder / File','ButtonPushedFcn',@onChooseSource);
c.sourceBtn.Layout.Row = 1; c.sourceBtn.Layout.Column = 1;
c.sessionBtn = uibutton(s1,'Text','Load session...','ButtonPushedFcn',@onLoadSession, ...
    'BackgroundColor',[0.90 0.94 1.0], ...
    'Tooltip',['read a Processing Workbench session: its recordings arrive as their ' ...
               '*_r.mat results, with animal / recording type / experimental group / ' ...
               'index already resolved and only the plots it actually computed offered']);
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
c.scanBtn.Layout.Row = 6; c.scanBtn.Layout.Column = [1 2];
c.fileList = uilistbox(s1,'Items',{},'Multiselect','on', ...
    'Tooltip',['Ctrl/Shift-click to select several files, type a name below and press ' ...
               '"Create group" to group them by hand (overrides the Group pattern).'], ...
    'ValueChangedFcn',@(~,~)requestRender());
c.fileList.Layout.Row = [7 10]; c.fileList.Layout.Column = [1 2];
c.groupNameF = labelledField(s1,11,'New group','');
c.groupNameF.Tooltip = 'Name for a hand-made group (see the file list).';
c.createBtn = uibutton(s1,'Text','Create group from selected files', ...
    'ButtonPushedFcn',@onCreateGroup, ...
    'Tooltip','Assign the files selected above to the group named on the left.');
c.createBtn.Layout.Row = 12; c.createBtn.Layout.Column = [1 2];

% --- section 2: what to plot ---
s2 = section(stack,'2 - What to plot',6);
c.dataType = labelledDrop(s2,1,'Data type', allDataTypes(), @onDataType);
c.variable = labelledDrop(s2,2,'Variable',{'(none)'},@(~,~)onVariable());
c.interval = labelledDrop(s2,3,'Interval',{'(all)'},@(~,~)requestRender());
c.interval.Tooltip = ['Which analysed window a myograph recording is plotted for.  ' ...
    'Switched off when nothing loaded has any.'];
c.interval.Enable = 'off';
c.organize = labelledDrop(s2,4,'Organise by', ...
    {'Auto','Group','Recording index','Animal','Type','Group x Index','Pool all'}, ...
    @(~,~)requestRender());
c.organize.Tooltip = ['The stratification axes are INDEPENDENT: an animal may span ' ...
    'experimental groups and a group may span animals.  Animal and Type are filled by ' ...
    'a loaded session, Interval by a myograph recording; Auto uses whichever of them ' ...
    'actually varies.'];
c.points   = labelledDrop(s2,5,'Points are', ...
    {'Auto (files if >1)','Files','Segments'},@(~,~)requestRender());
c.stat     = labelledDrop(s2,6,'Centre / error', ...
    {'Mean +/- SD','Mean +/- SEM','Median +/- IQR'},@(~,~)requestRender());

% --- section 3: selections ---
s3 = section(stack,'3 - Selection (multi-select)',4);
c.selList = uilistbox(s3,'Items',{'All vessels (lumen)'},'Multiselect','on', ...
    'ValueChangedFcn',@(~,~)requestRender());
c.selList.Layout.Row = [1 4]; c.selList.Layout.Column = [1 2];

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
    '%t type  %i interval  %f file  %v var'], ...
    'FontColor',[0.5 0.5 0.5],'FontSize',10,'WordWrap','on');
c.legHelp.Layout.Row = 6; c.legHelp.Layout.Column = [1 2];
c.cmap = labelledDrop(s4,7,'Image cmap',{'parula','turbo','hot','gray','jet'},@(~,~)requestRender());
lbFS = uilabel(s4,'Text','Font size'); lbFS.Layout.Row=8; lbFS.Layout.Column=1;
c.fontSize = uispinner(s4,'Value',12,'Limits',[6 40],'Step',1,'RoundFractionalValues','on', ...
    'ValueChangedFcn',@(~,~)requestRender(), ...
    'Tooltip','Base font size (points). Title and axis labels scale up from this; ticks, legend and colorbar match it.');
c.fontSize.Layout.Row=8; c.fontSize.Layout.Column=2;

% --- section 5: export ---
s5 = section(stack,'5 - Export',3);
c.dpi = labelledSpin(s5,1,'DPI',300);
c.fmt = labelledDrop(s5,2,'Format',{'PNG','TIFF','PDF (vector)','EPS (vector)','JPEG'},[]);
c.exportBtn = uibutton(s5,'Text','Export image','ButtonPushedFcn',@onExport);
c.exportBtn.Layout.Row = 3; c.exportBtn.Layout.Column = [1 2];

% --- plot button + status ---
sPlot = uigridlayout(stack,[2 1],'RowHeight',{34,'fit'},'Padding',[0 0 0 0], ...
    'RowSpacing',4,'BackgroundColor','w');
c.plotBtn = uibutton(sPlot,'Text','Plot / refresh','FontWeight','bold', ...
    'BackgroundColor',[0.90 0.94 1.0],'ButtonPushedFcn',@(~,~)doRender());
c.status  = uilabel(sPlot,'Text','Ready.','FontColor',[0.30 0.30 0.30],'WordWrap','on');

app.c = c;
onDataType();                       % initialise variable/selection item lists

% expose a small API so the tool can be driven programmatically (used for tests)
setappdata(fig,'exploreAPI',struct( ...
    'loadPaths',   @loadPathsProgrammatic, ...
    'loadSession', @(p) loadSessionFile(p), ...
    'sessionPath', @sessionPathLive, ...
    'dataTypes',   @dataTypesLive, ...
    'setState',    @setStateProgrammatic, ...
    'createGroup', @createGroupProgrammatic, ...
    'render',      @doRender, ...
    'export',      @exportTo, ...
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
                app.mode = 'file'; app.root = p; app.sessionPath = '';
                app.files = makeEntry(fullfile(p,f),'all','1');
                c.sourceLbl.Text = ['File: ' f];
                setPatternFields('off');
                refreshFileList(); refreshLabelsAndItems(); doRender();
        end
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
        app.mode = 'session'; app.sessionPath = pth;
        app.root = fileparts(pth);
        app.files = applyGroupOverrides(entries);
        c.sourceLbl.Text = sprintf('Session: %s   (%d recording(s) -> %d results file(s))', ...
            shortName(pth), numel(S.files), numel(app.files));
        setPatternFields('off');
        c.scanBtn.Enable = 'on';                 % re-reads the session (see onScan)
        c.scanBtn.Tooltip = 'Re-read the session file and rebuild the list.';
        refreshFileList(); refreshLabelsAndItems();
        if isempty(app.files)
            setStatus('The session has no *_r.mat results yet - process it first.');
            return
        end
        setStatus(sprintf('Session: %d results file(s), %d group(s), %d animal(s), %d type(s).', ...
            numel(app.files), nUnique('group'), nUnique('animal'), nUnique('type')));
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
            done = completedSteps(S, f.path);
            rp = resolveResultFiles(f.path);
            for k = 1:numel(rp)
                x = makeEntryStruct(rp{k}, '', '');
                x.group  = labelOr(f.expGroup,'all');
                x.rec    = labelOr(f.index,'1');
                x.recnum = recToNum(x.rec);
                x.animal = charOf(f.animal);
                x.type   = charOf(f.type);
                x.steps  = done;
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
            app.files = entries; refreshFileList();
            setStatus('No files matched the Include pattern.'); return;
        end
        app.files = applyGroupOverrides(entries);   % keep any hand-made groups
        refreshFileList(); refreshLabelsAndItems();
        setStatus(sprintf('Found %d files, %d group(s), %d index level(s).', ...
            numel(app.files), numel(unique({app.files.group})), numel(unique({app.files.rec}))));
        doRender();
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

    function onDataType(~,~)
        refreshVariableItems();
        onVariable();
    end

    function onVariable(~,~)
        refreshSelectionItems();
        autofillCosmetics();
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
        end
        c.fileList.Items = items;
        c.fileList.ItemsData = 1:numel(app.files);
        c.fileList.Value = 1:numel(app.files);     % all selected by default
    end

    function idx = selectedFileIdx()
        idx = c.fileList.Value;
        if isempty(idx), idx = 1:numel(app.files); end
        idx = idx(:)';
    end

    function refreshLabelsAndItems()
        % Collect the union of vessel/ROI labels across the (selected) files so the
        % selection list can offer "Label: X" entries. Reads only metrics tables.
        % The same pass collects what the myograph recordings offer - their analysed
        % windows and the signals inside them - because both answers come out of the
        % files that are already being opened here.
        L = {}; IV = {}; SG = {}; nSeg = 0; hasDia = false; hasProp = false;
        idx = 1:numel(app.files);
        for i = idx
            R = tryLoad(app.files(i).path);
            if isempty(R), continue; end
            for tn = {'sMetrics','dvsMetrics'}
                if isfield(R,tn{1}) && istable(R.(tn{1})) && ...
                        ismember('label',R.(tn{1}).Properties.VariableNames)
                    lv = strtrim(string(R.(tn{1}).label));
                    L = [L, cellstr(lv(strlength(lv)>0))']; %#ok<AGROW>
                end
            end
            if ~isMyoResults(R), nSeg = nSeg + 1; continue; end
            nm = myoIntervalNames(R);  IV = [IV, nm(~ismember(nm,IV))]; %#ok<AGROW>
            sg = {myoSignalList(R).name};
            SG = [SG, sg(~ismember(sg,SG))]; %#ok<AGROW>
            hasDia  = hasDia  || myoHasBranch(R,'diameter');
            hasProp = hasProp || myoHasBranch(R,'propagation');
        end
        app.labels = unique(L);
        % '.only' is what lets a myograph-only list stop offering vessel categories
        % nothing in it has, while a mixed list still offers both.  The two branch
        % flags do the same job one level down: a WIRE myograph has no diameter and
        % no propagation, so those plot families are not offered for it.
        app.myo = struct('intervals',{IV},'signals',{SG}, ...
            'only', ~isempty(IV) && nSeg==0, 'diameter',hasDia, 'propagation',hasProp);
        refreshIntervalItems();
        refreshDataTypeItems();
        refreshVariableItems();
        refreshSelectionItems();
        autofillCosmetics();
    end

    function refreshIntervalItems()
        %refreshIntervalItems  The interval filter, offered only when something
        %   loaded actually has windows - the axis is real or it is not there.
        items = [{'(all)'} app.myo.intervals];
        keep = c.interval.Value;
        c.interval.Items = items;
        if ismember(keep,items), c.interval.Value = keep; else, c.interval.Value = items{1}; end
        c.interval.Enable = onOff(~isempty(app.myo.intervals));
        % the Organise-by list grows an Interval entry for the same reason
        base = {'Auto','Group','Recording index','Animal','Type','Group x Index','Pool all'};
        if ~isempty(app.myo.intervals)
            base = [base(1:5) {'Interval'} base(6:end)];
        end
        keepO = c.organize.Value;
        c.organize.Items = base;
        if ismember(keepO,base), c.organize.Value = keepO; else, c.organize.Value = base{1}; end
    end

    function refreshDataTypeItems()
        %refreshDataTypeItems  Offer only the plot families the loaded files can
        %   actually produce.  With a session loaded that answer is KNOWN - the
        %   completion record says which steps ran - so a result tree that was never
        %   computed is dropped from the menu instead of being probed for and then
        %   found empty.  Without a session nothing is known and the full list stays,
        %   which is the behaviour the folder / file routes have always had.
        %
        %   THE MYOGRAPH FAMILIES ARE GATED TWICE, because a session is not the only
        %   way in: they also need a loaded file that actually carries analysed
        %   windows, so a folder of speckle results never offers them.
        items = allDataTypes();
        spec  = dataTypeSteps();
        for k = 1:size(spec,1)
            if ~sessionOffers(spec{k,2})
                items = items(~strcmp(items, spec{k,1}));
            end
        end
        if isempty(app.myo.intervals)
            items = items(~ismember(items, myoDataTypes()));
        else
            if ~app.myo.diameter
                items = items(~ismember(items, myoDiameterTypes()));
            end
            if ~app.myo.propagation
                items = items(~ismember(items, myoPropagationTypes()));
            end
        end
        if isempty(items), items = allDataTypes(); end
        keep = c.dataType.Value;
        c.dataType.Items = items;
        if ismember(keep,items), c.dataType.Value = keep; else, c.dataType.Value = items{1}; end
    end

    function tf = sessionOffers(stepIds)
        %sessionOffers  Did ANY loaded file complete ANY of these steps?  Tri-state:
        %   a file with no completion record contributes no knowledge, and when
        %   nothing is known at all the answer is yes (probe as before).  A LIST of
        %   step ids rather than one, because the same plot is produced by different
        %   steps in different modalities - the vasomotion of a speckle recording and
        %   of a myograph one draw the identical tree.
        stepIds = cellstr(stepIds);
        known = false;
        for i = 1:numel(app.files)
            if isempty(app.files(i).steps), continue; end
            known = true;
            if any(ismember(stepIds, app.files(i).steps)), tf = true; return; end
        end
        tf = ~known;
    end

%% ==================== DYNAMIC DROPDOWN ITEMS ============================ %%
    function refreshVariableItems()
        %refreshVariableItems  WHICH QUANTITY of the chosen family to draw.  For a
        %   myograph family the signal is picked in the selection list (it is what a
        %   comparison is about), so what is left for this menu is the quantity -
        %   which statistic, which propagation number - and the families that have
        %   only one offer exactly one entry rather than a misleading choice.
        dt = c.dataType.Value;
        idx = selectedFileIdx();
        R = firstAvailable(idx);
        items = {};
        switch dt
            case 'Time series'
                if hasField(R,'sData'),       items{end+1} = 'Segment BFI (sData)'; end
                if hasField(R,'dvsData'),     items{end+1} = 'DVS flow (dvsData)'; end
                if hasField(R,'dvsDiameter'), items{end+1} = 'DVS diameter (dvsDiameter)'; end
                if isempty(items), items = {'Segment BFI (sData)'}; end
            case 'Scalar metric'
                items = [ prefixCols(R,'sMetrics','[seg] '), prefixCols(R,'dvsMetrics','[dvs] ') ];
                if isempty(items), items = {'[seg] BFI'}; end
            case 'Vasomotion spectrum'
                items = vsmSignalItems(R);
            case 'Vasomotion time-freq'
                items = vsmSignalItems(R,'spectrum.amp');
            case 'Vasomotion percentile spectra'
                items = vsmSignalItems(R,'fVectors.VB.ampMeanPct');
            case 'Vasomotion amplitude percentiles'
                items = vsmSignalItems(R,'scalars.VB.ampPct');
            case 'Vasomotion marker'
                items = vsmMarkerItems(R, firstMyo(idx));
            case 'Diameter stats'
                items = myoStatItems(firstMyo(idx));
            case 'Propagation'
                items = myoPropItems(firstMyo(idx));
            case 'Diameter trace'
                items = {'diameter','% of baseline'};
            case {'Per-line diameter','Diameter map','Detected walls'}
                items = {'diameter'};
            case 'Propagation lag'
                items = {'lag'};
            case 'Image map'
                items = imageItems(R);
                if isempty(items), items = {'imgBFI'}; end
        end
        keep = c.variable.Value;
        c.variable.Items = items;
        if ismember(keep,items), c.variable.Value = keep; else, c.variable.Value = items{1}; end
    end

    function refreshSelectionItems()
        %refreshSelectionItems  WHAT IS COMPARED.  For a segmented recording that is
        %   a vessel category or a named label; a myograph recording has neither, so
        %   for those files it is the analysed diameter MEASURE or the recorded
        %   CHANNEL - the same idea one level down.  A mixed list offers both, and a
        %   file simply contributes nothing to a selection that is not about it.
        dom = currentDomain();
        idx = selectedFileIdx();
        lazyAny = any(arrayfun(@(i)isLazyFile(app.files(i).path), idx));
        if strcmp(c.dataType.Value,'Image map')
            items = {'(whole image)'};
        elseif isMyoDataType(c.dataType.Value) || app.myo.only
            items = app.myo.signals;
            if isempty(items), items = {'(nothing analysed)'}; end
        elseif strcmp(dom,'dvs')
            items = {'All DVS','Arteries (DVS)','Veins (DVS)','Uncertain (DVS)'};
        else
            items = {'All vessels (lumen)','Arteries (lumen)','Veins (lumen)', ...
                'Uncertain (lumen)','Parenchyma'};
            if ~isempty(app.myo.signals), items = [items app.myo.signals]; end
        end
        if ~lazyAny && ~strcmp(c.dataType.Value,'Image map')
            for i=1:numel(app.labels), items{end+1} = ['Label: ' app.labels{i}]; end %#ok<AGROW>
        end
        keep = c.selList.Value;
        c.selList.Items = items;
        keep = keep(ismember(keep,items));
        if isempty(keep), keep = items(1); end
        c.selList.Value = keep;
    end

    function R = firstMyo(idx)
        %firstMyo  The first loaded file that is a myograph recording, for filling a
        %   menu from what one of them actually carries.
        R = [];
        for i = idx
            Ri = tryLoad(app.files(i).path);
            if isMyoResults(Ri), R = Ri; return; end
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
            return;
        end
        hold(app.ax,'off');
    end

    function renderInto(ax)
        %renderInto  THE one dispatch from a plot family to its renderer, shared by
        %   the live axes and the export axes so the exported figure can never be a
        %   different plot from the one on screen.
        switch c.dataType.Value
            case 'Time series',           renderCurve('ts', ax);
            case 'Scalar metric',         renderMetric(ax);
            case 'Vasomotion spectrum',   renderCurve('spct', ax);
            case 'Vasomotion time-freq',  renderSpectrogram(ax);
            case 'Vasomotion percentile spectra',    renderPctSpectra(ax);
            case 'Vasomotion amplitude percentiles', renderAmpPct(ax);
            case 'Vasomotion marker',     renderMetric(ax);
            case 'Diameter trace',        renderCurve('diam', ax);
            case 'Per-line diameter',     renderPerLineDiameter(ax);
            case 'Diameter stats',        renderMetric(ax);
            case 'Propagation',           renderMetric(ax);
            case 'Propagation lag',       renderCurve('lag', ax);
            case 'Diameter map',          renderDiameterMap(ax);
            case 'Detected walls',        renderWalls(ax);
            case 'Image map',             renderImage(ax);
        end
    end

%% ==================== OBSERVATION BUILDERS ============================== %%
    function [obs,meta] = gatherCurveObservations(kind)
        % Build a flat list of "curve" observations (time series, spectrum or
        % diameter trace), each tagged with selection/group/rec/interval/file,
        % honouring the Points-are setting.  A MYOGRAPH file contributes one
        % observation per analysed WINDOW, which is what puts the interval on an
        % axis; a segmented file contributes one per file or per segment, as before,
        % with an empty interval tag that no axis ever separates on.
        idx  = selectedFileIdx();
        sels = c.selList.Value; if ischar(sels), sels = {sels}; end
        pooledPerFile = pointsAreFiles(numel(idx));
        obs = struct('x',{},'y',{},'sel',{},'group',{},'rec',{}, ...
            'animal',{},'type',{},'interval',{},'file',{});
        meta = struct('xlabel','','ylabel','','xlog',false);
        for s = 1:numel(sels)
            for i = idx
                R = tryLoad(app.files(i).path);
                if isempty(R), continue; end
                if isMyoResults(R)
                    [obs,meta] = gatherMyoCurves(obs,meta,R,i,kind,sels{s});
                    continue
                end
                if any(strcmp(kind,{'diam','lag'})), continue; end   % a segmented file has neither
                [x,Y,w,ylab,xlab,xlog] = curveMatrix(R, app.files(i).path, kind, sels{s});
                if isempty(Y), continue; end
                meta.xlabel=xlab; meta.ylabel=ylab; meta.xlog=xlog;
                isLabelSel = startsWith(sels{s},'Label:');
                if pooledPerFile || isLabelSel
                    y = wmean(Y,w,2);                    % one representative curve/file
                    obs(end+1) = mkObs(x,y,sels{s},i,''); %#ok<AGROW>
                else
                    for k = 1:size(Y,2)                 % each segment is an observation
                        obs(end+1) = mkObs(x,Y(:,k),sels{s},i,''); %#ok<AGROW>
                    end
                end
            end
        end
    end

    function [obs,meta] = gatherMyoCurves(obs,meta,R,i,kind,sel)
        %gatherMyoCurves  One curve per analysed WINDOW of one myograph recording.
        %   The units inside a window (the lines of a per-line vasomotion) are
        %   averaged into the window's curve: they are locations along one vessel,
        %   not independent recordings, so a window is one observation.
        for k = myoIntervalIdx(R, c.interval.Value)
            switch kind
                case 'diam'
                    [x,y] = myoDiameterTrace(R,k,sel);
                    if isempty(y), continue; end
                    meta.xlabel='time (s)'; meta.ylabel='diameter (px)'; meta.xlog=false;
                    if strcmp(c.variable.Value,'% of baseline')
                        % Each window against ITS OWN mean, which is how a myograph
                        % result is usually read: vessels of different calibre become
                        % comparable, and a dilation is a number a reader recognises.
                        base = mean(y,'omitnan');
                        if ~isfinite(base) || base==0, continue; end
                        y = 100*y/base; meta.ylabel='diameter (% of baseline)';
                    end
                case 'lag'
                    [x,y] = myoPropLag(R,k,sel);
                    if isempty(y), continue; end
                    meta.xlabel='position along the vessel (line)';
                    meta.ylabel='lag (s)'; meta.xlog=false;
                case 'spct'
                    V = myoVsmTree(R,k,sel);
                    if isempty(V) || ~isfield(V,'fVectors') || ~isfield(V.fVectors,'ampMean'), continue; end
                    x = double(V.f(:));
                    y = mean(double(V.fVectors.ampMean),1,'omitnan')';
                    meta.xlabel='frequency (Hz)'; meta.ylabel='amplitude (a.u.)'; meta.xlog=true;
                otherwise
                    continue
            end
            obs(end+1) = mkObs(x,y,sel,i,myoNameOf(R,k)); %#ok<AGROW>
        end
    end

    function o = mkObs(x,y,sel,i,interval)
        o = struct('x',x(:),'y',y(:),'sel',prettySel(sel), ...
            'group',app.files(i).group,'rec',app.files(i).rec, ...
            'animal',app.files(i).animal,'type',app.files(i).type, ...
            'interval',interval,'file',app.files(i).name);
    end

    function [vals,tags] = gatherScalarObservations()
        idx  = selectedFileIdx();
        sels = c.selList.Value; if ischar(sels), sels = {sels}; end
        [colDom,colName] = parseMetricVar(c.variable.Value);
        pooledPerFile = pointsAreFiles(numel(idx));
        vals = []; tags = struct('sel',{},'group',{},'rec',{}, ...
            'animal',{},'type',{},'interval',{},'file',{});
        myoKind = myoScalarKind(c.dataType.Value);
        for s = 1:numel(sels)
            for i = idx
                R = tryLoad(app.files(i).path);
                if isempty(R), continue; end
                if ~isempty(myoKind)
                    % a myograph scalar: one value per analysed window
                    if ~isMyoResults(R), continue; end
                    for k = myoIntervalIdx(R, c.interval.Value)
                        v = myoScalar(R,k,sels{s},myoKind,c.variable.Value);
                        if ~isfinite(v), continue; end
                        vals(end+1,1) = v; %#ok<AGROW>
                        tags(end+1)   = mkTag(sels{s},i,myoNameOf(R,k)); %#ok<AGROW>
                    end
                    continue
                end
                if strcmp(c.dataType.Value,'Vasomotion marker')
                    [vals,tags] = addMarkerValues(vals,tags,R,i,sels{s},pooledPerFile);
                    continue
                end
                if ~hasField(R,colDom), continue; end
                T = R.(colDom);
                if ~istable(T) || ~ismember(colName,T.Properties.VariableNames), continue; end
                rows = selectRows(T, sels{s}, domainOf(colDom));
                v = double(T.(colName)); v = v(rows);
                w = weightCol(T); w = w(rows);
                good = isfinite(v); v=v(good); w=w(good);
                if isempty(v), continue; end
                isLabelSel = startsWith(sels{s},'Label:');
                if pooledPerFile || isLabelSel
                    vals(end+1,1) = wmean(v,w,1); %#ok<AGROW>
                    tags(end+1)   = mkTag(sels{s},i,''); %#ok<AGROW>
                else
                    for k=1:numel(v)
                        vals(end+1,1) = v(k); %#ok<AGROW>
                        tags(end+1)   = mkTag(sels{s},i,''); %#ok<AGROW>
                    end
                end
            end
        end
    end

    function [vals,tags] = addMarkerValues(vals,tags,R,i,sel,pooledPerFile)
        %addMarkerValues  The band scalars of ONE recording, for one selection.  This
        %   is the one scalar family that serves BOTH kinds of recording, because
        %   scalars.<band>.<marker> is the identical leaf in both trees.  A segmented
        %   recording contributes the segments the selection picks, area-weighted when
        %   they are pooled, exactly as the Scalar metric family does; a myograph one
        %   contributes one value per analysed WINDOW, its units averaged, because they
        %   are locations along a single vessel rather than independent measurements -
        %   the same rule every other myograph family already follows.
        [sig,band,marker] = parseVsmMarkerVar(c.variable.Value);
        if isempty(marker), return; end
        if isMyoResults(R)
            for k = myoIntervalIdx(R, c.interval.Value)
                u = vsmMarkerValues(myoVsmTree(R,k,sel), band, marker);
                u = u(isfinite(u));
                if isempty(u), continue; end
                vals(end+1,1) = mean(u); %#ok<AGROW>
                tags(end+1)   = mkTag(sel,i,myoNameOf(R,k)); %#ok<AGROW>
            end
            return
        end
        tName = domainTable(sig);
        if ~hasField(R,tName), return; end
        v = vsmMarkerValues(getNested(R,['vasomotion.' sig]), band, marker);
        T = R.(tName);
        if isempty(v) || ~istable(T) || numel(v)~=height(T), return; end
        rows = selectRows(T, sel, domainOf(tName));
        w = weightCol(T); w = w(rows); v = v(rows);
        good = isfinite(v); v = v(good); w = w(good);
        if isempty(v), return; end
        if pooledPerFile || startsWith(sel,'Label:')
            vals(end+1,1) = wmean(v,w,1);
            tags(end+1)   = mkTag(sel,i,'');
        else
            for k = 1:numel(v)
                vals(end+1,1) = v(k); %#ok<AGROW>
                tags(end+1)   = mkTag(sel,i,''); %#ok<AGROW>
            end
        end
    end

    function t = mkTag(sel,i,interval)
        t = struct('sel',prettySel(sel),'group',app.files(i).group, ...
            'rec',app.files(i).rec,'animal',app.files(i).animal, ...
            'type',app.files(i).type,'interval',interval,'file',app.files(i).name);
    end

%% ==================== RENDERERS ========================================= %%
    function renderCurve(kind, ax)
        [obs,meta] = gatherCurveObservations(kind);
        if isempty(obs), title(ax,'No data for this selection'); return; end
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
            if pointsAreSegmentsNow()   % show the individual constituents faintly
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
        setStatus(sprintf('%d series, %d observations (%s).',numel(uk),numel(obs),statName));
    end

    function renderMetric(ax)
        [vals,tags] = gatherScalarObservations();
        if isempty(vals), title(ax,'No data for this selection'); legend(ax,'off'); return; end
        dims = activeDims(tags);
        xName = dims{1}; colDims = dims(2:end);
        xc  = arrayfun(@(t) t.(xName), tags, 'uni',0);
        xcats = orderedCats(xc);
        [~,colName] = parseMetricVar(c.variable.Value);
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
        setStatus(sprintf('%d x-category(ies), %d points.', numel(xcats), numel(vals)));
    end

    function X = vsmContext()
        %vsmContext  THE vasomotion tree the single-recording views draw, and which
        %   of its units go into the average.  ONE resolver for both kinds of
        %   recording, which is what lets those views be written once: a segmented
        %   recording answers with results.vasomotion.<signal> and the segments the
        %   selection picks, area-weighted as they have always been; a myograph
        %   recording answers with the chosen WINDOW's sub-tree - the SAME tree shape
        %   - and every unit in it, unweighted, because its units are locations along
        %   one vessel rather than segments of different sizes.
        X = struct('R',[],'V',[],'rows',[],'w',[],'label','');
        idx = selectedFileIdx(); sels = c.selList.Value; if ischar(sels), sels={sels}; end
        if isempty(idx) || isempty(sels), return; end
        X.R = tryLoad(app.files(idx(1)).path);
        if isempty(X.R), return; end
        if isMyoResults(X.R)
            k = myoIntervalIdx(X.R, c.interval.Value);
            if isempty(k), return; end
            X.V = myoVsmTree(X.R,k(1),sels{1});
            if isempty(X.V), return; end
            n = vsmUnits(X.V);
            X.rows = (1:n)'; X.w = ones(n,1); X.label = myoNameOf(X.R,k(1));
        else
            sig = vsmSignalOf(c.variable.Value);
            X.V = getNested(X.R, ['vasomotion.' sig]);
            if isempty(X.V), return; end
            T = X.R.(domainTable(sig));
            X.rows = selectRows(T, sels{1}, domainOf(domainTable(sig)));
            X.w = weightCol(T); X.w = X.w(X.rows);
        end
        if numel(idx)>1 || numel(sels)>1
            setStatus('This view shows one recording and one selection; showing the first.');
        end
    end

    function v = vsmAxis(X,name)
        %vsmAxis  A shared axis of the tree (f / timeDWT / pctCenters): the sub-tree's
        %   own copy when it has one, else the one at the results.vasomotion root.
        if isfield(X.V,name), v = double(X.V.(name)); else, v = double(getNested(X.R,['vasomotion.' name])); end
        v = v(:);
    end

    function renderSpectrogram(ax)
        X = vsmContext();
        if isempty(X.V) || ~isfield(X.V,'spectrum') || ~isfield(X.V.spectrum,'amp') ...
                || isempty(X.V.spectrum.amp)
            title(ax,'No stored time-frequency spectrum for this selection'); return;
        end
        S = double(X.V.spectrum.amp(X.rows,:,:));          % [nUnit x nF x nD] wavelet amplitude
        M = squeeze(wmean(permute(S,[2 3 1]), X.w, 3));    % [nF x nD] weighted mean
        f = vsmAxis(X,'f'); t = vsmAxis(X,'timeDWT');
        imagesc(ax, t, f, M); set(ax,'YDir','normal');
        try, set(ax,'YScale','log'); catch, end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String='amplitude (a.u.)';
        finishAxes(ax,'time (s)','frequency (Hz)',gobjects(0),{}); legend(ax,'off');
        grid(ax,'off'); axis(ax,'tight');
    end

    function renderPctSpectra(ax)
        % Percentile-resolved spectra: amplitude vs frequency, one line per VB-envelope
        % amplitude percentile bin (area-weighted mean over the selected segments), the
        % legend giving each line's bin centre. Mirrors Myograph plotPctSpectra; single
        % file / first selection, like the time-freq view. Bins are keyed to root
        % pctCenters (bin CENTRES) - distinct from the ampPct LEVELS below.
        X = vsmContext();
        if isempty(X.V) || ~isfield(X.V,'fVectors') || ~isfield(X.V.fVectors,'VB') || ...
                ~isfield(X.V.fVectors.VB,'ampMeanPct') || isempty(X.V.fVectors.VB.ampMeanPct)
            title(ax,'No percentile-resolved spectra for this selection'); return;
        end
        S = double(X.V.fVectors.VB.ampMeanPct(X.rows,:,:)); % [nUnit x nF x nB] spectrum per amplitude bin
        M = squeeze(wmean(permute(S,[2 3 1]), X.w, 3));     % [nF x nB] weighted mean over the units
        if isvector(M), M=M(:); end
        f = vsmAxis(X,'f'); pc = vsmAxis(X,'pctCenters');
        nB = min(numel(pc), size(M,2)); hLeg = gobjects(1,nB); names = cell(1,nB);
        for p = 1:nB                                       % one curve per envelope-amplitude bin (colours cycle the default order)
            hLeg(p) = plot(ax, f, M(:,p), 'DisplayName', sprintf('%gth pct', pc(p)));
            names{p} = sprintf('%gth pct', pc(p));
        end
        finishAxes(ax,'frequency (Hz)','amplitude (a.u.)',hLeg,names);
        set(ax,'XScale','log'); axis(ax,'tight');
    end

    function renderAmpPct(ax)
        % Scalar amplitude percentiles: band-envelope amplitude vs percentile LEVEL, VB
        % (and CB if stored), area-weighted mean over the selected segments. Mirrors
        % Myograph plotAmpPct; single file / first selection. The x-axis is the
        % percentile LEVELS s.pcts (0..100 %) - NOT the fVector bin CENTRES (pctCenters)
        % that key the percentile SPECTRA above.
        X = vsmContext();
        if isempty(X.V) || ~isfield(X.V,'scalars') || ~isfield(X.V.scalars,'VB') || ...
                ~isfield(X.V.scalars.VB,'ampPct') || isempty(X.V.scalars.VB.ampPct)
            title(ax,'No band amplitude percentiles for this selection'); return;
        end
        aVB = wmean(double(X.V.scalars.VB.ampPct(X.rows,:)), X.w, 1); aVB=aVB(:);
        nP = numel(aVB);
        % Percentile LEVELS: a *_r.mat results file carries no settings struct, so reuse
        % a stored s.pcts if a file provides one, else the s.pcts default (linspace(0,100,
        % nP) = 0:10:100 at nP=11); a level/count mismatch falls back to a plain index
        % (mirrors Myograph plotAmpPct).
        pcts = vsmAxis(X,'pcts');
        if isempty(pcts), pcts = linspace(0,100,nP)'; end
        if numel(pcts)~=nP, pcts = (1:nP)'; end
        hLeg = plot(ax, pcts, aVB, '-o', 'DisplayName','VB'); names = {'VB'};
        if isfield(X.V.scalars,'CB') && isfield(X.V.scalars.CB,'ampPct') && ~isempty(X.V.scalars.CB.ampPct)
            aCB = wmean(double(X.V.scalars.CB.ampPct(X.rows,:)), X.w, 1); aCB=aCB(:);
            if numel(aCB)==nP
                hLeg(2) = plot(ax, pcts, aCB, '-s', 'DisplayName','CB'); names{2}='CB';
            end
        end
        finishAxes(ax,'amplitude percentile (%)','amplitude (a.u.)',hLeg,names);
        set(ax,'XScale','linear');
    end

    function renderImage(ax)
        idx = selectedFileIdx(); i = idx(1);
        R = tryLoad(app.files(i).path);
        name = c.variable.Value;
        img = double(fetchField(R, app.files(i).path, name));
        if isempty(img), title(ax,['Field not found: ' name]); return; end
        if strcmp(name,'sMap')||strcmp(name,'pMap')||strcmp(name,'dvsMap')
            img(img==0)=NaN;
        end
        m = isfinite(img);
        imagesc(ax, img, 'AlphaData', m); axis(ax,'image');
        if any(m(:))
            lims = prctile(img(m),[2 99]);
            if lims(2)>lims(1), clim(ax,lims); end
        end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String = name;
        finishAxes(ax, '', '', gobjects(0), {}); legend(ax,'off');
        set(ax,'XTick',[],'YTick',[]); grid(ax,'off');
        if numel(idx)>1, setStatus('Image shows one file; showing the first selected.'); end
    end

    function renderDiameterMap(ax)
        %renderDiameterMap  THE DIAMETER ALONG THE VESSEL, position against time.
        %   This is the view that says whether detection held over the whole vessel
        %   for the whole window: a line that wandered shows as a stripe, and a wall
        %   that was lost shows as a band of missing values.  No myograph step writes
        %   a report page, so this is where that check lives.
        %
        %   It reads the per-line arrays, which live ONCE in the recording's SOURCE
        %   (the interval carries only the line-averaged trace), and it decimates in
        %   time on the way in - a two-hour recording is not drawn at full rate onto
        %   an axes a thousand pixels wide.
        [R,path] = firstMyoFile();
        if isempty(R), title(ax,'No myograph recording selected'); return; end
        k = myoIntervalIdx(R, c.interval.Value);
        if isempty(k), title(ax,'No analysed window to show'); return; end
        sels = c.selList.Value; if ischar(sels), sels={sels}; end
        [M,t,y,unit] = myoDiameterMap(R, path, k(1), sels{1});
        if isempty(M)
            title(ax,'The per-line diameter is not beside this result - keep the _MYO_d.mat with it');
            return
        end
        imagesc(ax, t, y, M'); set(ax,'YDir','normal'); axis(ax,'tight');
        m = isfinite(M);
        if any(m(:))
            lims = prctile(M(m),[2 98]);
            if lims(2)>lims(1), clim(ax,lims); end
        end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String = ['diameter (' unit ')'];
        finishAxes(ax,'time (s)','position along the vessel (line)',gobjects(0),{});
        legend(ax,'off'); grid(ax,'off');
        setStatus(sprintf('%s - %s, %d lines.', shortName(path), myoNameOf(R,k(1)), numel(y)));
    end

    function renderPerLineDiameter(ax)
        %renderPerLineDiameter  EVERY LINE'S OWN DIAMETER against time, rather than
        %   their average or their map.  The map says whether detection HELD; this says
        %   what the individual traces DID - whether the lines move together, which is
        %   the assumption a propagation speed is fitted under, and whether one stray
        %   line is carrying an average the rest do not support.  Lines are coloured by
        %   position so they can be read against the map, and the window's line average
        %   is drawn over them in black.
        [R,path] = firstMyoFile();
        if isempty(R), title(ax,'No myograph recording selected'); return; end
        k = myoIntervalIdx(R, c.interval.Value);
        if isempty(k), title(ax,'No analysed window to show'); return; end
        sels = c.selList.Value; if ischar(sels), sels={sels}; end
        [M,t,y,unit] = myoDiameterMap(R, path, k(1), sels{1});
        if isempty(M)
            title(ax,'The per-line diameter is not beside this result - keep the _MYO_d.mat with it');
            return
        end
        cm = colormap(ax, c.cmap.Value);
        ci = round(linspace(1, size(cm,1), max(size(M,2),1)));
        for j = 1:size(M,2)
            plot(ax, t, M(:,j), '-', 'Color',[cm(ci(j),:) 0.45],'LineWidth',0.5, ...
                'HandleVisibility','off');
        end
        hMean = plot(ax, t, mean(M,2,'omitnan'), '-','Color',[0 0 0],'LineWidth',1.8);
        if numel(y)>1
            clim(ax,[y(1) y(end)]); cb = colorbar(ax); cb.Label.String = 'line';
        end
        finishAxes(ax,'time (s)',['diameter (' unit ')'], hMean, {'line average'});
        axis(ax,'tight');
        setStatus(sprintf('%s - %s, %d lines.', shortName(path), myoNameOf(R,k(1)), size(M,2)));
    end

    function renderWalls(ax)
        %renderWalls  ONE FRAME OF THE RECORDING WITH THE DETECTED WALLS ON IT - the
        %   "did it find the right edges" check, and the other half of what a report
        %   page would have carried.  The walls are drawn RED when that frame was
        %   flagged invalid, which is a wall that had left the field of view: its
        %   diameter is a lower bound rather than a measurement, and that has to be
        %   visible rather than inferred.  The frame shown is the middle of the chosen
        %   window, which is the one most likely to be representative of it.
        [R,path] = firstMyoFile();
        if isempty(R), title(ax,'No myograph recording selected'); return; end
        k = myoIntervalIdx(R, c.interval.Value);
        if isempty(k), title(ax,'No analysed window to show'); return; end
        sels = c.selList.Value; if ischar(sels), sels={sels}; end
        W = myoWallFrame(R, path, k(1), sels{1});
        if isempty(W.frame)
            title(ax,'The recording is not beside this result - keep the video with it');
            return
        end
        image(ax, W.frame); axis(ax,'image');
        col = [0.1 0.9 0.1]; if ~W.valid, col = [0.95 0.15 0.15]; end
        plot(ax, W.left,  W.rows, '-', 'Color', col, 'LineWidth', 1.2);
        plot(ax, W.right, W.rows, '-', 'Color', col, 'LineWidth', 1.2);
        finishAxes(ax,'','',gobjects(0),{}); legend(ax,'off');
        set(ax,'XTick',[],'YTick',[]); grid(ax,'off');
        note = ''; if ~W.valid, note = '   wall out of view'; end
        setStatus(sprintf('%s - %s, %.1f s%s', shortName(path), myoNameOf(R,k(1)), W.time, note));
    end

    function [R,path] = firstMyoFile()
        %firstMyoFile  The myograph recording a single-recording view draws: the first
        %   selected one, with the path it came from, because these two views read the
        %   SOURCE beside it as well as the results.
        R = []; path = '';
        idx = selectedFileIdx();
        for i = idx
            Ri = tryLoad(app.files(i).path);
            if isMyoResults(Ri), R = Ri; path = app.files(i).path; break; end
        end
        if ~isempty(R) && numel(idx)>1
            setStatus('This view shows one recording; showing the first myograph one.');
        end
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
        % Ordered list of the tag dimensions that separate the data, given the
        % "Organise by" choice. The first is the natural x-axis (box plots) and the
        % combination defines a series (curves). Selection is always kept separate
        % when several are shown so different vessel types are never pooled together.
        % GROUP, INDEX, ANIMAL, TYPE and INTERVAL are five INDEPENDENT axes (spec §2)
        % - they never nest, so Auto simply takes whichever of them actually varies.
        % Animal and type are empty unless a session filled them, and interval unless
        % the recording is a myograph one, which is why a folder scan of speckle
        % results behaves exactly as it always did.
        switch c.organize.Value
            case 'Group',           prim={'group'};
            case 'Recording index', prim={'rec'};
            case 'Animal',          prim={'animal'};
            case 'Type',            prim={'type'};
            case 'Interval',        prim={'interval'};
            case 'Group x Index',   prim={'rec','group'};
            case 'Pool all',        prim={};
            otherwise                                  % Auto: use whatever varies
                prim={};
                if nDim(tg,'rec')>1,      prim{end+1}='rec';      end
                if nDim(tg,'group')>1,    prim{end+1}='group';    end
                if nDim(tg,'animal')>1,   prim{end+1}='animal';   end
                if nDim(tg,'type')>1,     prim{end+1}='type';     end
                if nDim(tg,'interval')>1, prim{end+1}='interval'; end
        end
        dims=prim;
        if nDim(tg,'sel')>1, dims{end+1}='sel'; end
        if isempty(dims), dims={'sel'}; end            % single group -> one x category
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
        dt = c.dataType.Value; v = c.variable.Value;
        sels = c.selList.Value; if ischar(sels), sels={sels}; end
        selTxt = strjoin(cellfun(@prettySel,sels,'uni',0),', ');
        t = sprintf('%s - %s', v, selTxt);
        if numel(app.files)>1, t = [t sprintf('  (n=%d files)',numel(selectedFileIdx()))]; end
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
        nm = strrep(nm,'%v', c.variable.Value);
        nm = strtrim(regexprep(nm,'\s+',' '));
        if isempty(nm), nm = tagLike.sel; end
    end

%% ==================== DATA EXTRACTION ================================== %%
    function [x,Y,w,ylab,xlab,xlog] = curveMatrix(R, path, kind, sel)
        x=[]; Y=[]; w=[]; ylab=''; xlab=''; xlog=false;
        if strcmp(kind,'ts')
            [dataField,ylab] = tsFieldOf(c.variable.Value);
            D = fetchField(R,path,dataField);
            t = fetchField(R,path,'time');
            if isempty(D)||isempty(t), return; end
            T = metricsTableFor(R,dataField);
            rows = selectRowsMaybeLazy(R,T,sel,dataField, size(D,2));
            Y = double(D(:,rows)); x = double(t(:)); xlab='time (s)';
            w = rowWeights(T,rows,numel(rows));
            keep = any(isfinite(Y),1); Y=Y(:,keep); w=w(keep);
        else % spectrum
            sig = vsmSignalOf(c.variable.Value);
            V = getNested(R, ['vasomotion.' sig]);
            if isempty(V)||~isfield(V,'fVectors'), return; end
            spName = spectrumField();
            Sp = double(getNested(V,spName));          % [nSeg x nF] mean amplitude spectrum
            if ndims(Sp)==3, Sp = Sp(:,:,1); end
            T = metricsTableFor(R, domainTable(sig));
            rows = selectRowsMaybeLazy(R,T,sel,domainTable(sig), size(Sp,1));
            Y = Sp(rows,:)';                           % [nF x nSel]
            % shared frequency axis lives at the results.vasomotion root; gsData carries its own
            if isfield(V,'f'), x = double(V.f(:)); else, x = double(getNested(R,'vasomotion.f')); x=x(:); end
            w = rowWeights(T,rows,numel(rows));
            ylab='amplitude (a.u.)'; xlab='frequency (Hz)'; xlog=true;
            keep = any(isfinite(Y),1); Y=Y(:,keep); w=w(keep);
        end
    end

    function rows = selectRowsMaybeLazy(R,T,sel,~,nCols)
        if istable(T)
            rows = find(selectRows(T, sel, tableDomain(T)));
        else
            rows = (1:nCols)';                     % lazy file: "All" over every column
        end
    end

    function w = rowWeights(T,rows,n)
        if istable(T), w = weightCol(T); w=w(rows); else, w=ones(n,1); end
        w = w(:);
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
        app.files = applyGroupOverrides(sortEntries(e));
        c.includeF.Value=includePat; c.groupF.Value=groupPat; c.recF.Value=recPat;
        setPatternFields('on');
        refreshFileList(); refreshLabelsAndItems();
    end

    function createGroupProgrammatic(name, idx)
        c.fileList.Value = idx; c.groupNameF.Value = name; onCreateGroup();
    end
    function a = getAppLive(),      a = app;              end
    function p = sessionPathLive(), p = app.sessionPath;  end
    function d = dataTypesLive(),   d = c.dataType.Items; end

    function setStateProgrammatic(varargin)
        for k=1:2:numel(varargin)
            key=varargin{k}; val=varargin{k+1};
            switch key
                case 'DataType',      c.dataType.Value=val; onDataType();
                case 'Variable',      c.variable.Value=val; onVariable();
                case 'Selection',     c.selList.Items=union(c.selList.Items,val,'stable'); c.selList.Value=val;
                case 'Interval'
                    if ~ismember(val,c.interval.Items), c.interval.Items=[c.interval.Items {char(val)}]; end
                    c.interval.Value=val;
                case 'Organise'
                    if ~ismember(val,c.organize.Items), c.organize.Items=[c.organize.Items {char(val)}]; end
                    c.organize.Value=val;
                case 'Points',        c.points.Value=val;
                case 'Stat',          c.stat.Value=val;
                case 'Title',         c.titleF.Value=val; app.manual.title=~isempty(val);
                case 'Legend',        c.legendChk.Value=logical(val);
                case 'LegendPattern', c.legendF.Value=val;
                case 'XLabel',        c.xlabF.Value=val; app.manual.xlab=~isempty(val);
                case 'YLabel',        c.ylabF.Value=val; app.manual.ylab=~isempty(val);
                case 'Cmap',          c.cmap.Value=val;
                case 'FontSize',      c.fontSize.Value=val;
                case 'Files',         c.fileList.Value=val;
            end
        end
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

%% ==================== small nested helpers ============================= %%
    function setStatus(msg), c.status.Text = msg; drawnow limitrate; end

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

    function n = spectrumField()   % which vasomotion spectrum to draw (kept simple: time-mean)
        n = 'fVectors.ampMean';
    end

    function v = fetchField(R, path, name)
        % Resolve results.<dotted.name> for a full struct OR a lazy HDF5 handle.
        if isstruct(R) && isfield(R,'x_lazy_')
            v = [];
            try, v = h5read(path, ['/results/' strrep(name,'.','/')]); catch, end
        else
            v = getNested(R, name);
        end
    end

    function tf = pointsAreFiles(nFiles)
        switch c.points.Value
            case 'Files',    tf=true;
            case 'Segments', tf=false;
            otherwise,       tf = nFiles>1;
        end
    end
    function tf = pointsAreSegmentsNow()
        tf = ~pointsAreFiles(numel(selectedFileIdx()));
    end

    function dom = currentDomain()
        switch c.dataType.Value
            case 'Time series',          dom = tsDomain(c.variable.Value);
            case 'Scalar metric',        [d,~]=parseMetricVar(c.variable.Value); dom=domainOf(d);
            case 'Vasomotion spectrum',  dom = domainOf(domainTable(vsmSignalOf(c.variable.Value)));
            case 'Vasomotion time-freq', dom = domainOf(domainTable(vsmSignalOf(c.variable.Value)));
            case 'Vasomotion percentile spectra',    dom = domainOf(domainTable(vsmSignalOf(c.variable.Value)));
            case 'Vasomotion amplitude percentiles', dom = domainOf(domainTable(vsmSignalOf(c.variable.Value)));
            case 'Vasomotion marker'
                sig = parseVsmMarkerVar(c.variable.Value);
                dom = domainOf(domainTable(sig));
            otherwise,                   dom = 'seg';
        end
    end

    function R = firstAvailable(idx)
        R = [];
        for i = idx
            R = tryLoad(app.files(i).path);
            if ~isempty(R), return; end
        end
    end

    function items = vsmSignalItems(R,need)
        % List the vasomotion signals present. `need` (optional) is a field path into
        % the <sig> sub-tree - a branch ('spectrum.amp') or a deeper path ('fVectors.VB.ampMeanPct',
        % 'scalars.VB.ampPct') - tested non-empty via getNested, so a signal is offered
        % only when that product was actually stored.
        % A MYOGRAPH recording keeps one tree per analysed SIGNAL inside each window,
        % and the signal is chosen in the selection list (it is what a comparison is
        % about), so there is nothing left for this menu to choose and it says so
        % with one entry rather than offering a choice that does nothing.
        if isMyoResults(R), items = {'Vasomotion'}; return; end
        items = {};
        sigs = {'sData','Segments (vasomotion.sData)'; 'dvsData','DVS flow (vasomotion.dvsData)'; ...
                'dvsDiameter','DVS diameter (vasomotion.dvsDiameter)'};
        for k=1:size(sigs,1)
            V = getNested(R, ['vasomotion.' sigs{k,1}]);
            if ~isempty(V) && (nargin<2 || ~isempty(getNested(V,need))), items{end+1}=sigs{k,2}; end %#ok<AGROW>
        end
        if isempty(items), items = {'Segments (vasomotion.sData)'}; end
    end

    function e = makeEntry(path,grp,rec)
        e = makeEntryStruct(path,'','');
        e.group = grp; e.rec = rec; e.recnum = recToNum(rec);
    end
end

%% ========================= LOCAL HELPERS ================================ %%
function s = section(parent, titleText, nBodyRows)
% A titled sub-panel holding a small 2-column grid of controls.
p = uipanel(parent,'Title',titleText,'FontWeight','bold','BackgroundColor','w', ...
    'ForegroundColor',[0.15 0.25 0.45]);
s = uigridlayout(p,[nBodyRows 2],'RowHeight',repmat({'fit'},1,nBodyRows), ...
    'ColumnWidth',{'fit','1x'},'RowSpacing',5,'ColumnSpacing',6, ...
    'Padding',[6 6 6 6],'BackgroundColor','w');
end

function f = labelledField(g,row,name,val)
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
f=uieditfield(g,'text','Value',val); f.Layout.Row=row; f.Layout.Column=2;
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
%   fields and the struct arrays concatenate.  animal / type / steps are empty
%   unless a session filled them.
e = struct('path',{},'name',{},'group',{},'rec',{},'recnum',{}, ...
    'animal',{},'type',{},'steps',{});
end

function e = makeEntryStruct(path, groupPat, recPat)
[~,nm,ex]=fileparts(path); name=[nm ex];
grp = regexpMatch(name, groupPat, 'all');
rec = regexpMatch(name, recPat,  '1');
e = struct('path',path,'name',name,'group',grp,'rec',rec,'recnum',recToNum(rec), ...
    'animal','','type','','steps',{{}});
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

function done = completedSteps(S, path)
%completedSteps  The step ids a session records as FINISHED for one recording.
%   The record is keyed 'path||stepId' (wbSession); anything that errored, was
%   skipped or never ran is simply not in the answer.  An empty answer means "this
%   session knows nothing about that file", which the caller reads as "probe".
done = {};
if ~isfield(S,'completed') || ~isa(S.completed,'containers.Map'), return; end
pre = [char(path) '||'];
k = keys(S.completed);
for i = 1:numel(k)
    if ~startsWith(k{i}, pre), continue; end
    v = S.completed(k{i});
    if isstruct(v) && isfield(v,'state') && ~strcmpi(char(v.state),'done'), continue; end
    done{end+1} = k{i}(numel(pre)+1:end); %#ok<AGROW>
end
done = unique(done);
end

function items = allDataTypes()
%allDataTypes  Every plot family this tool can draw, in menu order.
items = {'Time series','Scalar metric','Vasomotion spectrum','Vasomotion time-freq', ...
         'Vasomotion percentile spectra','Vasomotion amplitude percentiles', ...
         'Vasomotion marker','Diameter trace','Per-line diameter','Diameter stats', ...
         'Propagation','Propagation lag','Diameter map','Detected walls','Image map'};
end

function items = myoDataTypes()
%myoDataTypes  The families only a MYOGRAPH recording can produce, so they can be
%   dropped from the menu when nothing loaded is one.  The vasomotion families are
%   deliberately NOT here: they draw the same tree from either kind of recording -
%   including the marker box plot, which compares the band scalars of a speckle
%   recording exactly as readily as those of a myograph one.
items = [myoDiameterTypes() myoPropagationTypes()];
end

function items = myoDiameterTypes()
%myoDiameterTypes  The families that need a measured DIAMETER, so they are dropped
%   for a wire myograph, which records channels and never measures one.
items = {'Diameter trace','Per-line diameter','Diameter stats','Diameter map', ...
         'Detected walls'};
end

function items = myoPropagationTypes()
%myoPropagationTypes  The families that need a PROPAGATION estimate: the scalars it
%   concluded, and the per-line lags those were fitted to.
items = {'Propagation','Propagation lag'};
end

function tf = isMyoDataType(name)
tf = any(strcmp(name, myoDataTypes()));
end

function spec = dataTypeSteps()
%dataTypeSteps  Which PIPELINE STEP each plot family needs (wbStepRegistry ids).
%   Only families whose data is a whole result TREE are listed: a session that never
%   ran that step cannot offer them.  Families read straight off the standard
%   RESULTS members (time series, scalar metrics, image maps) are always offered.
%
%   EACH ENTRY IS A LIST, because the same plot comes out of different steps in
%   different modalities: the vasomotion of a speckle recording and of a myograph
%   one are the identical tree, written by 'vasomotion' and by 'myoVasomotion'.
spec = {'Vasomotion spectrum',              {'vasomotion','myoVasomotion'}; ...
        'Vasomotion time-freq',             {'vasomotion','myoVasomotion'}; ...
        'Vasomotion percentile spectra',    {'vasomotion','myoVasomotion'}; ...
        'Vasomotion amplitude percentiles', {'vasomotion','myoVasomotion'}; ...
        'Vasomotion marker',                {'vasomotion','myoVasomotion'}; ...
        'Diameter trace',                   {'myoDiameter'}; ...
        'Per-line diameter',                {'myoDiameter'}; ...
        'Diameter stats',                   {'myoDiameter'}; ...
        'Diameter map',                     {'myoDiameter'}; ...
        'Detected walls',                   {'myoDiameter'}; ...
        'Propagation',                      {'myoPropagation'}; ...
        'Propagation lag',                  {'myoPropagation'}};
end

function L = emptyLegendTag()
%emptyLegendTag  The blank legend record: one field per token legendName expands.
L = struct('sel','','group','','rec','','animal','','type','','interval','','file','');
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

function tf = hasField(R,name)
tf = ~isempty(R) && isstruct(R) && isfield(R,name);
end

function items = prefixCols(R, tname, prefix)
items = {};
if hasField(R,tname) && istable(R.(tname))
    T = R.(tname); vn = T.Properties.VariableNames;
    for i=1:numel(vn)
        if isnumeric(T.(vn{i})) && isvector(T.(vn{i}))
            items{end+1} = [prefix vn{i}]; %#ok<AGROW>
        end
    end
end
end

function items = imageItems(R)
items = {};
if isempty(R)||~isstruct(R), return; end
cand = {'imgBFI','imgK','imgI','mapType','sMap','pMap','dvsMap','cMask','mask'};
for i=1:numel(cand)
    if isfield(R,cand{i}) && ismatrix(R.(cand{i})) && ~isvector(R.(cand{i}))
        items{end+1}=cand{i}; %#ok<AGROW>
    end
end
if isfield(R,'extendedMetrics') && isstruct(R.extendedMetrics)
    fn=fieldnames(R.extendedMetrics);
    for i=1:numel(fn)
        v=R.extendedMetrics.(fn{i});
        if ismatrix(v)&&~isvector(v), items{end+1}=['extendedMetrics.' fn{i}]; end %#ok<AGROW>
    end
end
% relocated per-pixel pulsatility maps (results.pulsatility.ppx.scalars.*): the
% [Y x X] marker maps that used to be top-level imgPI/imgRI; hAmp/hPhase cubes are
% 3-D and skipped.  fetchField resolves the dotted name for eager and lazy files.
if isfield(R,'pulsatility') && isstruct(R.pulsatility) && isfield(R.pulsatility,'ppx') ...
        && isstruct(R.pulsatility.ppx) && isfield(R.pulsatility.ppx,'scalars') && isstruct(R.pulsatility.ppx.scalars)
    fn=fieldnames(R.pulsatility.ppx.scalars);
    for i=1:numel(fn)
        v=R.pulsatility.ppx.scalars.(fn{i});
        if ismatrix(v)&&~isvector(v), items{end+1}=['pulsatility.ppx.scalars.' fn{i}]; end %#ok<AGROW>
    end
end
end

%% ==================== THE MYOGRAPH ADAPTER ============================== %%
% Everything this tool needs to know about a myograph recording, in one place.  The
% PLOTS are the ones that already existed - these functions only answer "which tree,
% which trace, which unit", so that a myograph file can be handed to them.

function iv = myoFlat(R)
%myoFlat  THE WINDOWS OF A MYOGRAPH RECORDING, in the one shape this adapter reads.
%   A pressure myograph keeps them flat in results.intervals; a wire myograph splits
%   them by CHANNEL, because one LabChart file is several chambers.  myographIntervals
%   knows about both and hands back one list with .channelName written onto every
%   element, so every k below is an index into THAT list and the rest of this file
%   never has to ask which myograph it is looking at.
iv = [];
if isempty(R) || ~isstruct(R) || isfield(R,'x_lazy_'), return; end
if ~isfield(R,'intervals') && ~isfield(R,'channel'), return; end
iv = myographIntervals(R);
end

function iv = myoIntervalAt(R,k)
%myoIntervalAt  One window of the flat list, or [] when the index is not one.
iv = [];
IV = myoFlat(R);
if k>=1 && k<=numel(IV), iv = IV(k); end
end

function tf = isMyoResults(R)
%isMyoResults  Is this a myograph recording?  Decided by the DATA: a myograph result
%   holds its analysis inside the WINDOWS it was measured in, and nothing else does.
tf = ~isempty(myoFlat(R));
end

function nm = myoIntervalNames(R)
%myoIntervalNames  What the analysed windows are called, in recording order.  A wire
%   myograph names them PER CHANNEL: two chambers whose baseline window is called
%   'baseline' are two different windows, and a filter that could not tell them apart
%   would silently average the chambers together.
nm = {};
iv = myoFlat(R);
for k = 1:numel(iv), nm{end+1} = myoNameOf(R,k); end %#ok<AGROW>
end

function nm = myoNameOf(R,k)
%myoNameOf  One window's name, or what it is when nobody named it, qualified by the
%   channel it belongs to when there is one.
iv = myoFlat(R);
nm = '';
if k<1 || k>numel(iv), return; end
if isfield(iv(k),'name'), nm = char(string(iv(k).name)); end
if isempty(nm), nm = sprintf('interval %d',k); end
ch = '';
if isfield(iv(k),'channelName'), ch = char(string(iv(k).channelName)); end
if ~isempty(ch), nm = [nm ' - ' ch]; end
end

function idx = myoIntervalIdx(R,want)
%myoIntervalIdx  Which windows the Interval filter leaves: all of them, or the one
%   named.  A file that does not have the named window contributes nothing, which is
%   how one filter serves a list of recordings whose windows differ.
idx = [];
if ~isMyoResults(R), return; end
nm = myoIntervalNames(R);
if isempty(want) || strcmp(want,'(all)'), idx = 1:numel(nm); return; end
idx = find(strcmp(nm, want));
end

function S = myoSignalList(R)
%myoSignalList  WHAT WAS ANALYSED in this recording, as the user should see it and
%   as the tree stores it.  A pressure myograph analyses diameter MEASURES; a wire
%   myograph analyses CHANNELS, whose real names ('Force 1 (mN)') cannot be struct
%   field names and are therefore kept beside their trees.  Unioned over the windows
%   in the order they first appear.
S = struct('name',{},'field',{},'kind',{});
IV = myoFlat(R);
for k = 1:numel(IV)
    iv = IV(k);
    for b = {'diameter','propagation','vasomotion'}
        if ~isfield(iv,b{1}), continue; end
        B = iv.(b{1});
        if ~isstruct(B) || ~isscalar(B), continue; end
        if strcmp(b{1},'diameter')
            if ~isfield(B,'stats') || ~isstruct(B.stats), continue; end
            flds = fieldnames(B.stats)';
        else
            flds = fieldnames(B)';
        end
        for j = 1:numel(flds)
            nm = measureLabel(flds{j}); kind = 'measure';
            if strcmp(b{1},'vasomotion') && isstruct(B.(flds{j})) && ...
                    isfield(B.(flds{j}),'channelName') && ~isempty(B.(flds{j}).channelName)
                nm = char(string(B.(flds{j}).channelName)); kind = 'channel';
            end
            if any(strcmp(nm,{S.name})), continue; end
            S(end+1) = struct('name',nm,'field',flds{j},'kind',kind); %#ok<AGROW>
        end
    end
end
end

function tf = myoHasBranch(R,branch)
%myoHasBranch  Did ANY analysed window of this recording get this branch written?
tf = false;
iv = myoFlat(R);
if isempty(iv) || ~isfield(iv,branch), return; end
for k = 1:numel(iv)
    if ~isempty(iv(k).(branch)), tf = true; return; end
end
end

function nm = measureLabel(field)
%measureLabel  The diameter measures, named for a reader rather than for the code.
switch field
    case 'outer', nm = 'outer diameter';
    case 'mid',   nm = 'wall-centre diameter';
    case 'inner', nm = 'luminal diameter';
    otherwise,    nm = field;
end
end

function f = myoSignalField(R,name)
%myoSignalField  The tree field behind a selection, '' when this recording has none.
f = '';
S = myoSignalList(R);
k = find(strcmp({S.name}, name),1);
if ~isempty(k), f = S(k).field; end
end

function meas = myoMeasureNames(R)
%myoMeasureNames  The diameter measures this recording carries, in the order they
%   were measured - which is the order of the third dimension of source.data.
meas = {};
iv = myoFlat(R);
for k = 1:numel(iv)
    d = iv(k).diameter;
    if isstruct(d) && isscalar(d) && isfield(d,'stats') && isstruct(d.stats)
        meas = fieldnames(d.stats)'; return
    end
end
end

function m = myoMeasureIndex(R,name)
%myoMeasureIndex  Which plane of the per-line arrays a selection is, read off the
%   tree's own order rather than assumed.  Falls back to the wall centre.
m = min(2,max(1,numel(myoMeasureNames(R))));
k = find(strcmp(myoMeasureNames(R), myoSignalField(R,name)),1);
if ~isempty(k), m = k; end
end

function V = myoVsmTree(R,k,name)
%myoVsmTree  ONE window's vasomotion sub-tree for one signal.  This is the whole
%   lever: it is the same <VSM> tree shape as results.vasomotion.<signal>, so every
%   vasomotion plot works on it unchanged.
V = [];
iv = myoFlat(R);
if k<1 || k>numel(iv), return; end
T = iv(k).vasomotion;
if ~isstruct(T) || ~isscalar(T), return; end
f = myoSignalField(R,name);
if ~isempty(f) && isfield(T,f) && isstruct(T.(f)), V = T.(f); end
end

function n = vsmUnits(V)
%vsmUnits  How many results one signal's tree holds: one for a line-averaged trace
%   or a channel, one per image row when the vasomotion ran per line.
n = 0;
if ~isstruct(V) || ~isfield(V,'scalars'), return; end
if isfield(V.scalars,'VB') && isfield(V.scalars.VB,'ampMean')
    n = numel(V.scalars.VB.ampMean);
elseif isfield(V.scalars,'CB') && isfield(V.scalars.CB,'ampMean')
    n = numel(V.scalars.CB.ampMean);
end
end

function [x,y] = myoDiameterTrace(R,k,name)
%myoDiameterTrace  One window's line-averaged diameter, for one measure.
x = []; y = [];
iv = myoFlat(R);
if k<1 || k>numel(iv), return; end
d = iv(k).diameter;
if ~isstruct(d) || ~isscalar(d) || ~isfield(d,'time'), return; end
f = myoSignalField(R,name);
if isempty(f) || ~isfield(d,f), return; end
x = double(d.time(:)); y = double(d.(f)); y = y(:);
if numel(y)~=numel(x), x = []; y = []; end
end

function [x,y] = myoPropLag(R,k,name)
%myoPropLag  One window's per-line lag - THE EVIDENCE the propagation speed was
%   fitted to.  A straight run of lags against position is a wave crossing the field
%   of view; a scatter is a speed that should not be believed whatever its R2 says.
%   Stored as [row lag_s] and holding ONLY the lines the fit accepted, so the x axis
%   is the true position along the vessel and rejected lines are simply absent.
x = []; y = [];
iv = myoFlat(R);
if k<1 || k>numel(iv), return; end
p = iv(k).propagation;
if ~isstruct(p) || ~isscalar(p), return; end
f = myoSignalField(R,name);
if isempty(f) || ~isfield(p,f) || ~isstruct(p.(f)), return; end
L = p.(f).lagByRow;
if isempty(L) || ~isnumeric(L) || size(L,2)~=2, return; end
x = double(L(:,1)); y = double(L(:,2));
end

function [sig,band,marker] = parseVsmMarkerVar(varItem)
%parseVsmMarkerVar  Split a marker menu entry into the signal it belongs to, its band
%   and its name.  A MYOGRAPH entry carries no signal prefix - its signal is chosen in
%   the selection list, as it is for every myograph family - so the prefix is simply
%   absent and the returned signal is never used.
sig = 'sData'; s = char(varItem);
pre = {'[dvs diam] ','dvsDiameter'; '[dvs] ','dvsData'; '[seg] ','sData'};
for k = 1:size(pre,1)
    if startsWith(s,pre{k,1}), sig = pre{k,2}; s = s(numel(pre{k,1})+1:end); break; end
end
t = strsplit(strtrim(s),' ');
band = t{1}; marker = '';
if numel(t)>1, marker = strjoin(t(2:end),' '); end
end

function v = vsmMarkerValues(V,band,marker)
%vsmMarkerValues  One band scalar out of one <VSM> tree, as a column of units -
%   segments for a speckle recording, lines (or the single averaged trace) for a
%   myograph one.  Empty when that analysis level was never requested, which is how a
%   lean tree contributes nothing instead of erroring.
v = [];
if ~isstruct(V) || ~isscalar(V) || ~isfield(V,'scalars') || ~isstruct(V.scalars), return; end
S = V.scalars;
if ~isfield(S,band) || ~isstruct(S.(band)) || ~isfield(S.(band),marker), return; end
x = S.(band).(marker);
if ~isnumeric(x) && ~islogical(x), return; end
v = double(x(:));
end

function items = vsmMarkerItems(R,Rmyo)
%vsmMarkerItems  Every band SCALAR the loaded trees actually carry, band-qualified and
%   - for a speckle recording - prefixed with the signal it belongs to, the way the
%   Scalar metric menu already names its two tables.  A myograph tree needs no prefix.
%   Percentile leaves are left out on purpose: they are a vector per unit rather than a
%   number, and the amplitude-percentile family already draws them properly.
items = {};
sigs = {'sData','[seg] '; 'dvsData','[dvs] '; 'dvsDiameter','[dvs diam] '};
for k = 1:size(sigs,1)
    items = [items, vsmMarkerNames(getNested(R,['vasomotion.' sigs{k,1}]), sigs{k,2})]; %#ok<AGROW>
end
IVm = myoFlat(Rmyo);
if ~isempty(IVm)
    for k = 1:numel(IVm)
        T = IVm(k).vasomotion;
        if ~isstruct(T) || ~isscalar(T), continue; end
        f = fieldnames(T)';
        for j = 1:numel(f)
            items = [items, vsmMarkerNames(T.(f{j}), '')]; %#ok<AGROW>
        end
    end
end
items = unique(items,'stable');
if isempty(items), items = {'VB ampMean'}; end
end

function items = vsmMarkerNames(V,prefix)
%vsmMarkerNames  The scalar leaves of one <VSM> tree, as menu entries.  One column is
%   what makes a leaf a scalar per unit; a percentile leaf is [unit x percentile] and
%   is skipped by that test rather than by name.
items = {};
if ~isstruct(V) || ~isscalar(V) || ~isfield(V,'scalars') || ~isstruct(V.scalars), return; end
bands = fieldnames(V.scalars)';
for b = 1:numel(bands)
    B = V.scalars.(bands{b});
    if ~isstruct(B), continue; end
    f = fieldnames(B)';
    for j = 1:numel(f)
        x = B.(f{j});
        if isempty(x) || ~(isnumeric(x)||islogical(x)) || size(x,2)~=1, continue; end
        items{end+1} = [prefix bands{b} ' ' f{j}]; %#ok<AGROW>
    end
end
end

function items = myoStatItems(R)
%myoStatItems  The diameter statistics this recording carries.
items = {'mean'};
iv = myoFlat(R);
for k = 1:numel(iv)
    d = iv(k).diameter;
    if ~isstruct(d) || ~isscalar(d) || ~isfield(d,'stats'), continue; end
    fn = fieldnames(d.stats);
    if isempty(fn), continue; end
    items = fieldnames(d.stats.(fn{1}))'; return
end
end

function items = myoPropItems(R)
%myoPropItems  The propagation numbers worth comparing across recordings - the
%   answer and the two numbers that say whether to believe it.
want = {'speed','R2','confidence','pValue','domFreq','nRows'};
items = want;
iv = myoFlat(R);
for k = 1:numel(iv)
    p = iv(k).propagation;
    if ~isstruct(p) || ~isscalar(p), continue; end
    fn = fieldnames(p);
    if isempty(fn) || ~isstruct(p.(fn{1})), continue; end
    items = want(ismember(want, fieldnames(p.(fn{1}))'));
    if isempty(items), items = want; end
    return
end
end

function kind = myoScalarKind(dataType)
%myoScalarKind  Which myograph branch a scalar plot family reads ('' = not one).
switch dataType
    case 'Diameter stats', kind = 'stat';
    case 'Propagation',    kind = 'prop';
    otherwise,             kind = '';
end
end

function v = myoScalar(R,k,name,kind,varName)
%myoScalar  One number out of one window: a diameter statistic or a propagation
%   result, for the selected measure.
v = NaN;
f = myoSignalField(R,name);
if isempty(f), return; end
iv = myoFlat(R);
if k<1 || k>numel(iv), return; end
switch kind
    case 'stat'
        d = iv(k).diameter;
        if ~isstruct(d) || ~isscalar(d) || ~isfield(d,'stats'), return; end
        if ~isfield(d.stats,f) || ~isfield(d.stats.(f),varName), return; end
        v = scalarOf(d.stats.(f).(varName));
    case 'prop'
        p = iv(k).propagation;
        if ~isstruct(p) || ~isscalar(p) || ~isfield(p,f), return; end
        if ~isfield(p.(f),varName), return; end
        v = scalarOf(p.(f).(varName));
end
end

function v = scalarOf(x)
v = NaN;
if ~isempty(x) && (isnumeric(x)||islogical(x)), x = double(x); v = x(1); end
end

% ---- the two views that read the recording's SOURCE, not just its results ----
% The per-line arrays and the wall positions live ONCE, in the SOURCE beside the
% results (the window carries only the line-averaged trace), so these two read it -
% and read it in BLOCKS when the file is large, because a two-hour recording's
% arrays are not something to pull into memory to draw one picture.

function S = myoSourceHandle(rPath)
%myoSourceHandle  The SOURCE beside a myograph RESULTS file: the whole struct when
%   it is small, an HDF5 handle when it is not - the same rule the results
%   themselves already follow.
S = [];
d = strrep(char(rPath),'_r.mat','_d.mat');
if ~isfile(d), return; end
info = dir(d);
if ~isempty(info) && info(1).bytes > 1.5e9
    S = struct('x_lazy_',true,'x_path_',d);
else
    L = load(d,'source');
    if isfield(L,'source'), S = L.source; end
end
end

function v = myoSourceField(S,name)
%myoSourceField  One small SOURCE field (a time base, a row range, a flag).
v = [];
if isempty(S), return; end
if isfield(S,'x_lazy_')
    try, v = h5read(S.x_path_, ['/source/' name]); catch, end
else
    if isfield(S,name), v = S.(name); end
end
end

function sz = myoSourceSize(S,name)
%myoSourceSize  The shape of one big SOURCE array, without reading it.
sz = [];
if isempty(S), return; end
if isfield(S,'x_lazy_')
    try
        I = h5info(S.x_path_, ['/source/' name]); sz = I.Dataspace.Size;
    catch
    end
else
    if isfield(S,name), sz = size(S.(name)); end
end
end

function A = myoSlab(S,name,start,count,stride)
%myoSlab  One block of a big SOURCE array: `count` elements per dimension, spaced
%   `stride` apart, starting at `start`.  Read straight out of the file when the
%   SOURCE is a handle, indexed out of the struct when it is not.
A = [];
if isempty(S), return; end
if isfield(S,'x_lazy_')
    try, A = h5read(S.x_path_, ['/source/' name], start, count, stride); catch, end
    return
end
if ~isfield(S,name), return; end
subs = cell(1,numel(start));
for i = 1:numel(start)
    subs{i} = start(i) : stride(i) : start(i)+(count(i)-1)*stride(i);
end
A = S.(name)(subs{:});
end

function p = myoRecordingPath(S,rPath)
%myoRecordingPath  THE recording this product was made from: the path the entry step
%   wrote down, else the video sitting beside the product - the same order the
%   diameter step uses, so this finds one exactly when the pipeline would.
p = '';
v = myoSourceField(S,'fName');
if ~isempty(v)
    if isnumeric(v), v = char(v(:)'); end
    p = char(string(v));
end
if ~isempty(p) && isfile(p), return; end
[fPath,stem] = fileparts(regexprep(char(rPath),'_MYO_r\.mat$',''));
for ext = {'.avi','.mp4','.mov','.mkv'}
    cand = fullfile(fPath,[stem ext{1}]);
    if isfile(cand), p = cand; return; end
end
p = '';
end

function rows = myoMeasuredRows(S,nY)
%myoMeasuredRows  The image rows the vessel was actually measured on.  Rows outside
%   them carry an interpolated fill and are not a measurement, so nothing draws them.
rr = double(myoSourceField(S,'rowRange'));
rows = 1:nY;
if numel(rr)==2
    r0 = max(1,round(rr(1))); r1 = min(nY,round(rr(2)));
    if r1>=r0, rows = r0:r1; end
end
end

function [M,t,y,unit] = myoDiameterMap(R,rPath,k,name)
%myoDiameterMap  The per-line diameter of one window as [time x line], decimated in
%   time so a long recording draws at the size of the axes rather than of the file.
M = []; t = []; y = []; unit = 'px';
S = myoSourceHandle(rPath);
if isempty(S), return; end
iv = myoIntervalAt(R,k);
if isempty(iv) || ~isfield(iv,'frames') || numel(iv.frames)~=2, return; end
i0 = double(iv.frames(1)); i1 = double(iv.frames(2));
sz = myoSourceSize(S,'data');
if numel(sz)<3 || i1>sz(1) || i1<i0, return; end
rows = myoMeasuredRows(S,sz(2));
m = min(myoMeasureIndex(R,name), sz(3));
nT = i1-i0+1;
dec = max(1, ceil(nT/2000));
n = floor((nT-1)/dec)+1;
M = double(myoSlab(S,'data',[i0 rows(1) m],[n numel(rows) 1],[dec 1 1]));
if isempty(M), return; end
tAll = double(myoSourceField(S,'time'));
if numel(tAll)>=i1, t = tAll(i0:dec:i0+(n-1)*dec); else, t = (0:n-1)'; end
t = t(:)'; y = rows;
if ~isempty(myoSourceField(S,'pixelSize')), unit = 'px'; end
end

function W = myoWallFrame(R,rPath,k,name)
%myoWallFrame  One frame of the recording, with the detected walls of the selected
%   measure over it.  The MIDDLE frame of the window, as the one most likely to be
%   representative of it, and the walls come back flagged when that frame was
%   invalid - a wall had left the field of view, so its diameter is a lower bound.
W = struct('frame',[],'left',[],'right',[],'rows',[],'valid',true,'time',NaN);
S = myoSourceHandle(rPath);
if isempty(S), return; end
iv = myoIntervalAt(R,k);
if isempty(iv) || ~isfield(iv,'frames') || numel(iv.frames)~=2, return; end
sz = myoSourceSize(S,'wallL');
if numel(sz)<3, return; end
fi = min(max(1,round(mean(double(iv.frames)))), sz(1));
rows = myoMeasuredRows(S,sz(2));
m = min(myoMeasureIndex(R,name), sz(3));
W.left  = double(myoSlab(S,'wallL',[fi rows(1) m],[1 numel(rows) 1],[1 1 1])); W.left = W.left(:);
W.right = double(myoSlab(S,'wallR',[fi rows(1) m],[1 numel(rows) 1],[1 1 1])); W.right = W.right(:);
W.rows  = rows(:);
if numel(W.left)~=numel(W.rows) || isempty(W.left), W = resetWalls(W); return; end
tAll = double(myoSourceField(S,'time'));
if numel(tAll)>=fi, W.time = tAll(fi); end
vv = myoSourceField(S,'valid');
if ~isempty(vv) && numel(vv)>=fi, W.valid = logical(vv(fi)); end
video = myoRecordingPath(S,rPath);
if isempty(video), return; end
try
    vr = VideoReader(video);
    vr.CurrentTime = min(max(W.time,0), max(0, vr.Duration - 1/max(vr.FrameRate,eps)));
    W.frame = readFrame(vr);
catch
    W.frame = [];
end
end

function W = resetWalls(W)
W.left = []; W.right = []; W.rows = []; W.frame = [];
end

function s = onOff(tf), if tf, s = 'on'; else, s = 'off'; end, end

%% ======================================================================= %%
function v = getNested(R, name)
% Fetch results.<a>.<b>... (dotted) from a full struct; [] if missing/lazy.
v = [];
if isempty(R)||~isstruct(R)||isfield(R,'x_lazy_'), return; end
parts = strsplit(name,'.'); v=R;
for i=1:numel(parts)
    if isstruct(v)&&isfield(v,parts{i}), v=v.(parts{i}); else, v=[]; return; end
end
end

function [field,ylab] = tsFieldOf(varItem)
if     contains(varItem,'sData'),       field='sData';       ylab='BFI (a.u.)';
elseif contains(varItem,'dvsData'),     field='dvsData';     ylab='BFI (a.u.)';
elseif contains(varItem,'dvsDiameter'), field='dvsDiameter'; ylab='diameter (px)';
else,  field='sData'; ylab='BFI (a.u.)';
end
end
function dom = tsDomain(varItem)
if contains(varItem,'dvs'), dom='dvs'; else, dom='seg'; end
end

function sig = vsmSignalOf(varItem)
if     contains(varItem,'dvsDiameter'), sig='dvsDiameter';
elseif contains(varItem,'dvsData'),     sig='dvsData';
else,  sig='sData';
end
end
function t = domainTable(sig)
if strcmp(sig,'sData'), t='sMetrics'; else, t='dvsMetrics'; end
end
function dom = domainOf(tableOrField)
if any(strcmp(tableOrField,{'dvsMetrics','dvsData','dvsDiameter'})), dom='dvs'; else, dom='seg'; end
end

function [dom,name] = parseMetricVar(varItem)
if startsWith(varItem,'[dvs] '), dom='dvsMetrics'; name=varItem(7:end);
elseif startsWith(varItem,'[seg] '), dom='sMetrics'; name=varItem(7:end);
else, dom='sMetrics'; name=varItem;
end
end

function T = metricsTableFor(R,field)
switch field
    case {'sData','sMetrics'},                 T=tbl(R,'sMetrics');
    case {'dvsData','dvsDiameter','dvsMetrics'}, T=tbl(R,'dvsMetrics');
    otherwise, T=[];
end
end
function T = tbl(R,name)
if hasField(R,name)&&istable(R.(name)), T=R.(name); else, T=[]; end
end
function dom = tableDomain(T)
if ismember('category',T.Properties.VariableNames), dom='seg'; else, dom='dvs'; end
end

function rows = selectRows(T, sel, dom)
% Logical row mask into a metrics table for a named selection.
n = height(T); rows = false(n,1);
vn = T.Properties.VariableNames;
hasCat = ismember('category',vn); hasType = ismember('type',vn); hasLbl = ismember('label',vn);
if hasCat, isLum = (T.category==5); else, isLum = true(n,1); end   % DVS rows are all lumen
switch true
    case startsWith(sel,'Label:')
        if hasLbl
            want = strtrim(sel(7:end));
            m = strtrim(string(T.label))==want;
            if strcmp(dom,'seg') && hasCat
                r = m & (T.category==5);
                if ~any(r), r = m & (T.category==1); end  % parenchyma-only label
                rows = r;
            else
                rows = m;
            end
        end
    case contains(sel,'Arter'), if hasType, rows = isLum & (string(T.type)=="Artery"); end
    case contains(sel,'Vein'),  if hasType, rows = isLum & (string(T.type)=="Vein"); end
    case contains(sel,'Uncert'),if hasType, rows = isLum & (string(T.type)=="Uncertain"); end
    case contains(sel,'Parenchyma'), if hasCat, rows = (T.category==1); end
    case contains(sel,'All DVS'), rows = true(n,1);
    case contains(sel,'All vessels'), rows = isLum;
    otherwise, rows = isLum;
end
% drop rows with a non-finite idx (empty table slots)
if ismember('idx',vn)
    rows = rows & isfinite(double(T.idx));
end
end

function w = weightCol(T)
vn = T.Properties.VariableNames;
if ismember('area',vn), w=double(T.area);
elseif ismember('length',vn), w=double(T.length);
else, w=ones(height(T),1);
end
w(~isfinite(w)|w<=0)=eps;
end

function m = wmean(Y,w,dim)
% Weighted mean along dim, ignoring NaNs (equal weights if w empty). The weight
% vector is reshaped onto `dim` so it broadcasts against 2-D or 3-D Y.
if nargin<3, dim=2; end
if isempty(w), m=mean(Y,dim,'omitnan'); return; end
shp=ones(1,max(dim,2)); shp(dim)=numel(w); w=reshape(w(:),shp);
W = w.*isfinite(Y); Yz=Y; Yz(~isfinite(Y))=0;
m = sum(Yz.*W,dim)./max(sum(W,dim),eps);
end

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
