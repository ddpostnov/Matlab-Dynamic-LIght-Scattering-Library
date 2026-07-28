% guiExplore - Interactive GUI to explore & export publication figures from LSCI results.
%
% WHAT THIS TOOL DOES
%   A single-window app for browsing the *_r.mat RESULTS produced by this library
%   (runSegmentation -> runPulsatility / runVasomotion / runCTTH / setVesselTypes)
%   and turning them into clean, publication-ready plots. It works on ONE file or
%   on MANY files at once, and it always matches the plot type to the kind of data:
%
%     * time series (sData / dvsData / dvsDiameter)      -> mean +/- error band vs time
%     * scalar metrics (sMetrics / dvsMetrics columns)   -> grouped box plots + points
%     * vasomotion spectra (results.vasomotion.*)        -> amplitude vs frequency
%     * vasomotion time-frequency (...spectrum.amp)      -> 2-D spectrogram
%     * vasomotion percentile spectra (fVectors.VB.ampMeanPct) -> spectra split by
%                                                          envelope-amplitude percentile bin
%     * vasomotion amplitude percentiles (scalars.VB/CB.ampPct) -> band amplitude vs
%                                                          percentile level (VB & CB)
%     * image maps (imgBFI / pulsatility.ppx maps / mapType) -> axis-image display
%
% WHEN TO USE IT
%   After you have processed recordings to *_BFI_r.mat (or *_I_r.mat) results and
%   want to (a) look at single-recording detail, or (b) compare experimental groups
%   and/or a sequence of recordings across many files, then export the figure at a
%   chosen resolution for a paper or a talk.
%
% HOW TO USE IT
%   1. Normally you do NOT call this directly: open  guiWorkbench  and switch to its
%      Explore tab, which hosts this app and seeds it with the workbench's loaded
%      recordings and groups. To browse results without the workbench, run
%      guiExplore  standalone (with the library on the path) and continue from 2.
%   2. Click "Choose Folder / File".
%        - Pick a FILE  -> everything is plotted for that single recording.
%        - Pick a FOLDER -> the folder is searched recursively; then use the three
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
%      Click "Scan / apply" to (re)build the file list. Hover any of the three
%      regexp boxes for a popup with worked pattern examples.
%      Prefer to group by hand? Select files in the list (Ctrl/Shift-click), type a
%      name in "New group" and press "Create group from selected files". Hand-made
%      groups override the Group pattern and survive a re-scan.
%   3. Choose WHAT to plot with the dropdowns (data type -> variable), pick one or
%      more SELECTIONS (arteries / veins / parenchyma / all / a named label), tweak
%      titles / labels / legend, then "Plot". Labels that span several segments are
%      area-weighted into one trace/value automatically.
%   4. Set DPI + format and click "Export image".
%
% EMBEDDING
%   guiExplore('Parent',container) builds the whole app INTO an existing uitab or
%   uipanel instead of its own window (used by the Processing Workbench's Explore
%   tab, seeded with the workbench's file set / groups via the programmatic API);
%   with no 'Parent' it opens standalone exactly as described above.
%
%   This file therefore lives in GUIs/workbench/ next to the other workbench
%   components: it is the implementation of guiWorkbench's Explore tab rather than
%   a separate entry point. It stays fully usable standalone (genpath puts the whole
%   tree on the path), but guiWorkbench is the documented way in.
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
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

function h = guiExplore(varargin)

% ---- optional host container: with 'Parent',<uitab/uipanel> the whole app is
%      built INTO that container (e.g. the Processing Workbench's Explore tab)
%      instead of a standalone window; dialogs and figure-scoped lookups target
%      the host figure.  With no 'Parent' the behaviour is unchanged.
parent = parseParent(varargin);

% ---- shared application state (seen by all nested functions) --------------
app = struct();
app.mode        = '';        % 'file' | 'folder'
app.root        = '';        % chosen folder, or the folder of the chosen file
app.files       = struct('path',{},'name',{},'group',{},'rec',{},'recnum',{});
app.cache       = containers.Map('KeyType','char','ValueType','any');
app.cacheOrder  = {};        % LRU order of cached paths
app.cacheLimit  = 6;         % max fully-loaded files kept in memory
app.sizeLimitGB = 3;         % above this a file is read field-by-field (HDF5)
app.labels      = {};        % union of vessel/ROI labels across the file list
app.groupOverride = containers.Map('KeyType','char','ValueType','char'); % hand-made groups (path->name)
app.manual      = struct('title',false,'xlab',false,'ylab',false); % user-edited?

CAT = struct('PARENCHYMA',1,'UNSEG',2,'WALL',3,'INNER',4,'LUMEN',5);

% ---- build the window (standalone) or host inside the given container ------
if isempty(parent)
    fig = uifigure('Name','guiExplore - LSCI results explorer', ...
        'Color','w','Position',[80 60 1440 860]);
    try, fig.WindowState = 'maximized'; catch, end
    root = fig;                              % standalone: build on the figure
else
    root = parent;                           % hosted: build on the given container
    fig  = ancestor(parent,'figure');        % dialogs / findobj target the host figure
end
app.fig = fig;

outer = uigridlayout(root,[1 2],'ColumnWidth',{'2.4x','1x'}, ...
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
    'BackgroundColor','w');
stack = uigridlayout(C,[8 1],'RowHeight',repmat({'fit'},1,8), ...
    'RowSpacing',8,'Padding',[0 0 0 0],'BackgroundColor','w');

c = struct();   % control handles

% --- section 1: data source ---
s1 = section(stack,'1 - Data source',12);
c.sourceBtn = uibutton(s1,'Text','Choose Folder / File','ButtonPushedFcn',@onChooseSource);
c.sourceBtn.Layout.Row = 1; c.sourceBtn.Layout.Column = [1 2];
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
s2 = section(stack,'2 - What to plot',5);
c.dataType = labelledDrop(s2,1,'Data type', ...
    {'Time series','Scalar metric','Vasomotion spectrum','Vasomotion time-freq', ...
     'Vasomotion percentile spectra','Vasomotion amplitude percentiles','Image map'}, ...
    @onDataType);
c.variable = labelledDrop(s2,2,'Variable',{'(none)'},@(~,~)onVariable());
c.organize = labelledDrop(s2,3,'Organise by', ...
    {'Auto','Group','Recording index','Group x Index','Pool all'},@(~,~)requestRender());
c.points   = labelledDrop(s2,4,'Points are', ...
    {'Auto (files if >1)','Files','Segments'},@(~,~)requestRender());
c.stat     = labelledDrop(s2,5,'Centre / error', ...
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
c.legHelp = uilabel(s4,'Text','tokens: %s sel  %g group  %r index  %f file  %v var', ...
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
    'setState',    @setStateProgrammatic, ...
    'createGroup', @createGroupProgrammatic, ...
    'render',      @doRender, ...
    'export',      @exportTo, ...
    'getApp',      @getAppLive));

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
                app.mode = 'folder'; app.root = d;
                c.sourceLbl.Text = ['Folder: ' d];
                c.includeF.Enable='on'; c.groupF.Enable='on'; c.recF.Enable='on'; c.scanBtn.Enable='on';
                onScan();
            case 'Single file'
                [f,p] = uigetfile({'*_r.mat;*.mat','LSCI results (*_r.mat)'}, ...
                    'Select a results file');
                if isequal(f,0), return; end
                app.mode = 'file'; app.root = p;
                app.files = makeEntry(fullfile(p,f),'all','1');
                c.sourceLbl.Text = ['File: ' f];
                c.includeF.Enable='off'; c.groupF.Enable='off'; c.recF.Enable='off'; c.scanBtn.Enable='off';
                refreshFileList(); refreshLabelsAndItems(); doRender();
        end
    end

    function onScan(~,~)
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
        if ~strcmp(app.mode,'folder')
            uialert(fig,'Manual groups apply when a folder is loaded.','Create group'); return;
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
            items{i} = sprintf('%s   [%s | %s]', app.files(i).name, ...
                app.files(i).group, app.files(i).rec);
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
        L = {};
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
        end
        app.labels = unique(L);
        refreshVariableItems();
        refreshSelectionItems();
        autofillCosmetics();
    end

%% ==================== DYNAMIC DROPDOWN ITEMS ============================ %%
    function refreshVariableItems()
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
            case 'Image map'
                items = imageItems(R);
                if isempty(items), items = {'imgBFI'}; end
        end
        keep = c.variable.Value;
        c.variable.Items = items;
        if ismember(keep,items), c.variable.Value = keep; else, c.variable.Value = items{1}; end
    end

    function refreshSelectionItems()
        dom = currentDomain();
        idx = selectedFileIdx();
        lazyAny = any(arrayfun(@(i)isLazyFile(app.files(i).path), idx));
        if strcmp(c.dataType.Value,'Image map')
            items = {'(whole image)'};
        elseif strcmp(dom,'dvs')
            items = {'All DVS','Arteries (DVS)','Veins (DVS)','Uncertain (DVS)'};
        else
            items = {'All vessels (lumen)','Arteries (lumen)','Veins (lumen)', ...
                'Uncertain (lumen)','Parenchyma'};
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
            switch c.dataType.Value
                case 'Time series',           renderCurve('ts', app.ax);
                case 'Scalar metric',         renderMetric(app.ax);
                case 'Vasomotion spectrum',   renderCurve('spct', app.ax);
                case 'Vasomotion time-freq',  renderSpectrogram(app.ax);
                case 'Vasomotion percentile spectra',    renderPctSpectra(app.ax);
                case 'Vasomotion amplitude percentiles', renderAmpPct(app.ax);
                case 'Image map',             renderImage(app.ax);
            end
        catch ME
            hold(app.ax,'off');
            title(app.ax,'Could not plot this combination');
            setStatus(['Plot error: ' ME.message]);
            return;
        end
        hold(app.ax,'off');
    end

%% ==================== OBSERVATION BUILDERS ============================== %%
    function [obs,meta] = gatherCurveObservations(kind)
        % Build a flat list of "curve" observations (time series or spectrum), each
        % tagged with selection/group/rec/file, honouring the Points-are setting.
        idx  = selectedFileIdx();
        sels = c.selList.Value; if ischar(sels), sels = {sels}; end
        pooledPerFile = pointsAreFiles(numel(idx));
        obs = struct('x',{},'y',{},'sel',{},'group',{},'rec',{},'file',{});
        meta = struct('xlabel','','ylabel','','xlog',false);
        for s = 1:numel(sels)
            for i = idx
                R = tryLoad(app.files(i).path);
                if isempty(R), continue; end
                [x,Y,w,ylab,xlab,xlog] = curveMatrix(R, app.files(i).path, kind, sels{s});
                if isempty(Y), continue; end
                meta.xlabel=xlab; meta.ylabel=ylab; meta.xlog=xlog;
                isLabelSel = startsWith(sels{s},'Label:');
                if pooledPerFile || isLabelSel
                    y = wmean(Y,w,2);                    % one representative curve/file
                    obs(end+1) = mkObs(x,y,sels{s},i); %#ok<AGROW>
                else
                    for k = 1:size(Y,2)                 % each segment is an observation
                        obs(end+1) = mkObs(x,Y(:,k),sels{s},i); %#ok<AGROW>
                    end
                end
            end
        end
    end

    function o = mkObs(x,y,sel,i)
        o = struct('x',x(:),'y',y(:),'sel',prettySel(sel), ...
            'group',app.files(i).group,'rec',app.files(i).rec,'file',app.files(i).name);
    end

    function [vals,tags] = gatherScalarObservations()
        idx  = selectedFileIdx();
        sels = c.selList.Value; if ischar(sels), sels = {sels}; end
        [colDom,colName] = parseMetricVar(c.variable.Value);
        pooledPerFile = pointsAreFiles(numel(idx));
        vals = []; tags = struct('sel',{},'group',{},'rec',{},'file',{});
        for s = 1:numel(sels)
            for i = idx
                R = tryLoad(app.files(i).path);
                if isempty(R) || ~hasField(R,colDom), continue; end
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
                    tags(end+1)   = mkTag(sels{s},i); %#ok<AGROW>
                else
                    for k=1:numel(v)
                        vals(end+1,1) = v(k); %#ok<AGROW>
                        tags(end+1)   = mkTag(sels{s},i); %#ok<AGROW>
                    end
                end
            end
        end
    end

    function t = mkTag(sel,i)
        t = struct('sel',prettySel(sel),'group',app.files(i).group, ...
            'rec',app.files(i).rec,'file',app.files(i).name);
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
        nSel=numel(unique({tags.sel})); nGrp=numel(unique({tags.group})); nRec=numel(unique({tags.rec}));
        dims = activeDims(nSel,nGrp,nRec);
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

    function renderSpectrogram(ax)
        idx = selectedFileIdx(); sels = c.selList.Value; if ischar(sels), sels={sels}; end
        i = idx(1); sel = sels{1};
        R = tryLoad(app.files(i).path);
        sig = vsmSignalOf(c.variable.Value);
        V = getNested(R, ['vasomotion.' sig]);
        if isempty(V) || ~isfield(V,'spectrum') || ~isfield(V.spectrum,'amp')
            title(ax,'No stored 2-D spectrum (run runVasomotion with ''spectrum'' in s.segVsmReturn)'); return;
        end
        rows = selectRows(R.(domainTable(sig)), sel, domainOf(domainTable(sig)));
        w = weightCol(R.(domainTable(sig))); w=w(rows);
        S = double(V.spectrum.amp(rows,:,:));              % [nSeg x nF x nD] wavelet amplitude (already |CWT|)
        M = squeeze(wmean(permute(S,[2 3 1]), w, 3));      % [nF x nD] weighted mean
        % shared axes live at the results.vasomotion root; gsData carries its own
        if isfield(V,'f'), f=double(V.f(:)); else, f=double(getNested(R,'vasomotion.f')); f=f(:); end
        if isfield(V,'timeDWT'), t=double(V.timeDWT(:)); else, t=double(getNested(R,'vasomotion.timeDWT')); t=t(:); end
        imagesc(ax, t, f, M); set(ax,'YDir','normal');
        try, set(ax,'YScale','log'); catch, end
        colormap(ax, c.cmap.Value); cb=colorbar(ax); cb.Label.String='amplitude (a.u.)';
        if numel(idx)>1 || numel(sels)>1
            setStatus('Time-freq shows one file/selection; showing the first.');
        end
        finishAxes(ax,'time (s)','frequency (Hz)',gobjects(0),{}); legend(ax,'off');
        grid(ax,'off'); axis(ax,'tight');
    end

    function renderPctSpectra(ax)
        % Percentile-resolved spectra: amplitude vs frequency, one line per VB-envelope
        % amplitude percentile bin (area-weighted mean over the selected segments), the
        % legend giving each line's bin centre. Mirrors Myograph plotPctSpectra; single
        % file / first selection, like the time-freq view. Bins are keyed to root
        % pctCenters (bin CENTRES) - distinct from the ampPct LEVELS below.
        idx = selectedFileIdx(); sels = c.selList.Value; if ischar(sels), sels={sels}; end
        i = idx(1); sel = sels{1};
        R = tryLoad(app.files(i).path);
        sig = vsmSignalOf(c.variable.Value);
        V = getNested(R, ['vasomotion.' sig]);
        if isempty(V) || ~isfield(V,'fVectors') || ~isfield(V.fVectors,'VB') || ...
                ~isfield(V.fVectors.VB,'ampMeanPct') || isempty(V.fVectors.VB.ampMeanPct)
            title(ax,'No percentile spectra (fVectors.VB.ampMeanPct) - run runVasomotion with ''moments'' in s.segVsmReturn'); return;
        end
        rows = selectRows(R.(domainTable(sig)), sel, domainOf(domainTable(sig)));
        w = weightCol(R.(domainTable(sig))); w=w(rows);
        S = double(V.fVectors.VB.ampMeanPct(rows,:,:));    % [nSel x nF x nB] mean spectrum per VB amplitude bin
        M = squeeze(wmean(permute(S,[2 3 1]), w, 3));      % [nF x nB] area-weighted mean over the selected segments
        if isvector(M), M=M(:); end
        % shared axes live at the results.vasomotion root; gsData carries its own
        if isfield(V,'f'), f=double(V.f(:)); else, f=double(getNested(R,'vasomotion.f')); f=f(:); end
        if isfield(V,'pctCenters'), pc=double(V.pctCenters(:)); else, pc=double(getNested(R,'vasomotion.pctCenters')); pc=pc(:); end
        nB = min(numel(pc), size(M,2)); hLeg = gobjects(1,nB); names = cell(1,nB);
        for p = 1:nB                                       % one curve per envelope-amplitude bin (colours cycle the default order)
            hLeg(p) = plot(ax, f, M(:,p), 'DisplayName', sprintf('%gth pct', pc(p)));
            names{p} = sprintf('%gth pct', pc(p));
        end
        if numel(idx)>1 || numel(sels)>1
            setStatus('Percentile spectra show one file/selection; showing the first.');
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
        idx = selectedFileIdx(); sels = c.selList.Value; if ischar(sels), sels={sels}; end
        i = idx(1); sel = sels{1};
        R = tryLoad(app.files(i).path);
        sig = vsmSignalOf(c.variable.Value);
        V = getNested(R, ['vasomotion.' sig]);
        if isempty(V) || ~isfield(V,'scalars') || ~isfield(V.scalars,'VB') || ...
                ~isfield(V.scalars.VB,'ampPct') || isempty(V.scalars.VB.ampPct)
            title(ax,'No amplitude percentiles (scalars.VB.ampPct) - run runVasomotion with ''bands'' in s.segVsmReturn'); return;
        end
        rows = selectRows(R.(domainTable(sig)), sel, domainOf(domainTable(sig)));
        w = weightCol(R.(domainTable(sig))); w=w(rows);
        aVB = wmean(double(V.scalars.VB.ampPct(rows,:)), w, 1); aVB=aVB(:);   % [nP x 1] area-weighted mean over the selected segments
        nP = numel(aVB);
        % Percentile LEVELS: a *_r.mat results file carries no settings struct, so reuse
        % a stored s.pcts if a file provides one, else the s.pcts default (linspace(0,100,
        % nP) = 0:10:100 at nP=11); a level/count mismatch falls back to a plain index
        % (mirrors Myograph plotAmpPct).
        pcts = double(getNested(R,'vasomotion.pcts'));
        if isempty(pcts), pcts = linspace(0,100,nP)'; end
        pcts = pcts(:);
        if numel(pcts)~=nP, pcts = (1:nP)'; end
        hLeg = plot(ax, pcts, aVB, '-o', 'DisplayName','VB'); names = {'VB'};
        if isfield(V.scalars,'CB') && isfield(V.scalars.CB,'ampPct') && ~isempty(V.scalars.CB.ampPct)
            aCB = wmean(double(V.scalars.CB.ampPct(rows,:)), w, 1); aCB=aCB(:);
            if numel(aCB)==nP
                hLeg(2) = plot(ax, pcts, aCB, '-s', 'DisplayName','CB'); names{2}='CB';
            end
        end
        if numel(idx)>1 || numel(sels)>1
            setStatus('Amplitude percentiles show one file/selection; showing the first.');
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

    function dims = activeDims(nSel,nGrp,nRec)
        % Ordered list of the tag dimensions that separate the data, given the
        % "Organise by" choice. The first is the natural x-axis (box plots) and the
        % combination defines a series (curves). Selection is always kept separate
        % when several are shown so different vessel types are never pooled together.
        switch c.organize.Value
            case 'Group',           prim={'group'};
            case 'Recording index', prim={'rec'};
            case 'Group x Index',   prim={'rec','group'};
            case 'Pool all',        prim={};
            otherwise                                  % Auto: use whatever varies
                prim={};
                if nRec>1, prim{end+1}='rec';   end
                if nGrp>1, prim{end+1}='group'; end
        end
        dims=prim;
        if nSel>1, dims{end+1}='sel'; end
        if isempty(dims), dims={'sel'}; end            % single group -> one x category
    end

    function nm = legendDyn(val)
        % Apply the legend pattern to a single dynamic category (box colour groups).
        nm = legendName(struct('sel','','group','','rec','','file','','dyn',val));
    end

    function [key,leg] = seriesKeys(obs)
        % A series (one line+band) = unique combination of the ACTIVE dimensions.
        nSel=numel(unique({obs.sel})); nGrp=numel(unique({obs.group})); nRec=numel(unique({obs.rec}));
        dims=activeDims(nSel,nGrp,nRec);
        key=cell(1,numel(obs)); leg=cell(1,numel(obs));
        for i=1:numel(obs)
            key{i}=strjoin(cellfun(@(d) obs(i).(d), dims,'uni',0),'|');
            L=struct('sel','','group','','rec','','file',obs(i).file);
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
            nm = strtrim(strrep(strrep(strrep(strrep(strrep(pat, ...
                '%s',''),'%g',''),'%r',''),'%f',''),'%v',''));
            nm = strtrim([nm ' ' tagLike.dyn]); if isempty(strtrim(nm)), nm=tagLike.dyn; end
            return;
        end
        nm = pat;
        nm = strrep(nm,'%s', pad0(tagLike.sel));
        nm = strrep(nm,'%g', pad0(tagLike.group));
        nm = strrep(nm,'%r', pad0(tagLike.rec));
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
        app.mode='folder';
        if ~isempty(paths), app.root = fileparts(paths{1}); end
        e = struct('path',{},'name',{},'group',{},'rec',{},'recnum',{});
        for i=1:numel(paths)
            [~,nm,ex]=fileparts(paths{i});
            if ~isempty(includePat) && isempty(regexp([nm ex],includePat,'once')), continue; end
            e(end+1) = makeEntryStruct(paths{i}, groupPat, recPat); %#ok<AGROW>
        end
        app.files = applyGroupOverrides(sortEntries(e));
        c.includeF.Value=includePat; c.groupF.Value=groupPat; c.recF.Value=recPat;
        refreshFileList(); refreshLabelsAndItems();
    end

    function createGroupProgrammatic(name, idx)
        c.fileList.Value = idx; c.groupNameF.Value = name; onCreateGroup();
    end
    function a = getAppLive(), a = app; end

    function setStateProgrammatic(varargin)
        for k=1:2:numel(varargin)
            key=varargin{k}; val=varargin{k+1};
            switch key
                case 'DataType',      c.dataType.Value=val; onDataType();
                case 'Variable',      c.variable.Value=val; onVariable();
                case 'Selection',     c.selList.Items=union(c.selList.Items,val,'stable'); c.selList.Value=val;
                case 'Organise',      c.organize.Value=val;
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
        switch c.dataType.Value
            case 'Time series',          renderCurve('ts', axe);
            case 'Scalar metric',        renderMetric(axe);
            case 'Vasomotion spectrum',  renderCurve('spct', axe);
            case 'Vasomotion time-freq', renderSpectrogram(axe);
            case 'Vasomotion percentile spectra',    renderPctSpectra(axe);
            case 'Vasomotion amplitude percentiles', renderAmpPct(axe);
            case 'Image map',            renderImage(axe);
        end
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
        [~,nm,ex]=fileparts(path);
        e = struct('path',path,'name',[nm ex],'group',grp,'rec',rec,'recnum',recToNum(rec));
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
entries = struct('path',{},'name',{},'group',{},'rec',{},'recnum',{});
for i=1:numel(d)
    if isempty(regexp(d(i).name, includePat, 'once')), continue; end
    entries(end+1) = makeEntryStruct(fullfile(d(i).folder,d(i).name), groupPat, recPat); %#ok<AGROW>
end
entries = sortEntries(entries);
end

function e = makeEntryStruct(path, groupPat, recPat)
[~,nm,ex]=fileparts(path); name=[nm ex];
grp = regexpMatch(name, groupPat, 'all');
rec = regexpMatch(name, recPat,  '1');
e = struct('path',path,'name',name,'group',grp,'rec',rec,'recnum',recToNum(rec));
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
switch name, case 'rec', xl='recording index'; case 'group', xl='group'; otherwise, xl='selection'; end
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

function p = parseParent(args)
% Pull an optional 'Parent',<container> pair out of varargin ([] if absent).
p = [];
for i = 1:2:numel(args)-1
    if (ischar(args{i}) || isstring(args{i})) && strcmpi(args{i},'Parent')
        p = args{i+1};
    end
end
end
