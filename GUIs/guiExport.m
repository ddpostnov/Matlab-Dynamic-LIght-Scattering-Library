%guiExport - Standalone export tool: pick results, choose the averaging, write workbooks.
%
% WHAT THIS TOOL DOES
%   A single-window app in front of  exportToExcel  (Utilities/).  It answers the
%   three questions a spreadsheet export actually poses - WHICH recordings, WHICH
%   parameters, and averaged HOW - and then calls the one routine that owns the
%   science.  It computes nothing itself: every number in a preview or a workbook
%   comes back out of exportToExcel.
%
%   It is STANDALONE.  It does not need - and never loads - the Processing
%   Workbench: the hand-off between them is the SESSION file (wbSession), read
%   through wbSession('read',...), which carries the curated file list together
%   with each file's ANIMAL, recording TYPE and experimental GROUP.  Those three
%   labels are what turn a pile of workbooks into a table you can do statistics on.
%
% WHEN TO USE IT
%   After processing, whenever you want results in Excel: one workbook per
%   recording as the launchers have always written them, or ONE merged workbook
%   whose rows carry file / animal / type / group and are ready for a statistics
%   package.
%
% HOW TO USE IT
%   1. Get files in, three ways (mix them freely):
%        - "Add files..."     pick *_r.mat products by hand;
%        - "Scan folder..."   a recursive glob (getFileNamesList), e.g. '*_BFI_r.mat';
%        - "Load session..."  a workbench session - the file list arrives with its
%                             animal / type / group labels already resolved.
%      Every input is resolved to the exportable products of that recording: its
%      *_r.mat RESULTS files.
%   2. Curate the list.  Untick 'use' to leave a file out; edit animal / type /
%      group in place for files that came without labels; filter the table by any
%      of the three labels to export one experimental group at a time.  The filter
%      values are built FROM THE DATA - the vocabularies are yours, not the code's.
%   3. Choose parameters.  The picker lists the union of the sMetrics columns the
%      selected files actually contain, each annotated with how many of them have
%      it ("psPI - 42/60 files"), so availability is visible before anything is
%      written rather than discovered afterwards.  A column a file lacks is simply
%      skipped for that file.
%   4. Choose the averaging: per segment, over LABELS (the classic area-weighted
%      ROI sheets), or over VESSEL TYPE (artery / vein / ..., from setVesselTypes).
%      Area weighting is forced for labels - that is what those sheets have always
%      been - and free for vessel type.
%      (VESSEL type is not the RECORDING type of the Files list: one describes a
%      segment, the other describes the recording.)
%   5. Choose the output: one workbook per file, or merged into a single workbook
%      at a destination you pick.  "Preview" writes the first selected file to a
%      scratch workbook with exactly the chosen options and shows you the first
%      rows of the resulting sheet BEFORE anything real is written.
%   6. "Export" writes them, one file at a time, with a progress line.
%
% PROGRAMMATIC USE
%   getappdata(fig,'exportAPI') exposes the same logic as a struct of handles -
%   .loadFiles .loadSession .scanFolder .setUse .setFilter .setOpts .parameters
%   .files .opts .previewTable .run - so the tool is fully headless-testable
%   (claude-tests/Workbench/testExportTool.m), exactly like guiExplore's API.
%
% NO LAUNCHER IS NEEDED TO OPEN IT
%   Every time this window opens it puts the library on the MATLAB path itself, from
%   its own location (setLibraryPath), so typing guiExport in a fresh MATLAB is the
%   whole of the setup - and a stale copy of the library left on the path by an
%   earlier session is corrected rather than quietly used.
%
% Syntax:
%    guiExport                        % open the tool
%    guiExport(sessionPath)           % open it ON a workbench session
%    h = guiExport('Visible','off')   % headless (programmatic drive / tests)
%
% See also: exportToExcel, wbSession, guiExplore, guiWorkbench, getFileNamesList,
%           setLibraryPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

%------------- BEGIN CODE --------------
function h = guiExport(varargin)

% ---- THE LIBRARY PATH, SET FROM THIS FILE'S OWN LOCATION -------------------
% A front window is usually the FIRST thing typed into a fresh MATLAB, with no
% launcher run and often nothing on the path but the folder this file is in, so it
% puts the library there itself instead of assuming somebody already did.
% Utilities goes on first because that is where setLibraryPath lives; setLibraryPath
% then does the rest, and keeps the .claude tooling copies OFF the path - they hold
% whole checkouts of this library at older commits and SHADOW it silently.
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(libraryFolder,'Utilities'));
setLibraryPath(libraryFolder);

[sessionArg, vis] = parseArgs(varargin);

% ---- shared application state (seen by every nested function) --------------
app = struct();
app.items       = emptyItems();     % ONE ENTRY PER EXPORTABLE PRODUCT
app.view        = [];               % item indices currently shown (after filters)
app.probe       = containers.Map('KeyType','char','ValueType','any'); % rpath -> probe
app.colOn       = containers.Map('KeyType','char','ValueType','any'); % column -> ticked
app.colChk      = struct();         % column -> checkbox handle
app.sessionPath = '';
app.lastRoot    = pwd;
app.previewRows = 12;

AXES = {'animal','type','group'};   % the three label axes (open vocabularies)

delete(findall(groot,'Type','figure','Tag','guiExport'));
fig = uifigure('Name','guiExport - export processed results to Excel', ...
    'Color','w','Position',[70 60 1380 840],'Visible',vis,'Tag','guiExport');
app.fig = fig;

outer = uigridlayout(fig,[2 2],'RowHeight',{'fit','1x'},'ColumnWidth',{'2.1x','1x'}, ...
    'Padding',[6 6 6 6],'RowSpacing',8,'ColumnSpacing',8,'BackgroundColor','w');

c = struct();

% ===================== 1 - data source (spans both columns) ================
s1p = uipanel(outer,'Title','1 - Data source','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
s1p.Layout.Row = 1; s1p.Layout.Column = [1 2];
s1 = uigridlayout(s1p,[2 7],'RowHeight',{'fit','fit'}, ...
    'ColumnWidth',{'fit','fit','fit','fit','1x','fit','fit'}, ...
    'Padding',[6 6 6 6],'RowSpacing',5,'ColumnSpacing',6,'BackgroundColor','w');

c.addBtn = uibutton(s1,'Text','Add files...','ButtonPushedFcn',@onAddFiles, ...
    'Tooltip','pick processed *_r.mat products by hand');
c.scanBtn = uibutton(s1,'Text','Scan folder...','ButtonPushedFcn',@onScanFolder, ...
    'Tooltip','search a folder recursively for the glob on the right');
c.sessBtn = uibutton(s1,'Text','Load session...','ButtonPushedFcn',@onLoadSession, ...
    'BackgroundColor',[0.90 0.94 1.0], ...
    'Tooltip',['read a Processing Workbench session: the file list arrives with ' ...
               'its animal / recording type / experimental group already resolved']);
uilabel(s1,'Text','Glob');
c.globF = uieditfield(s1,'text','Value','*_BFI_r.mat', ...
    'Tooltip',['what "Scan folder..." looks for, searched recursively.  ' ...
               '*_BFI_r.mat for processed speckle recordings, *_MYO_r.mat for ' ...
               'myograph ones']);
c.clearBtn = uibutton(s1,'Text','Clear list','ButtonPushedFcn',@(~,~)clearList(), ...
    'Tooltip','empty the file list');
c.reprobeBtn = uibutton(s1,'Text','Re-read results','ButtonPushedFcn',@(~,~)reprobeAll(), ...
    'Tooltip','re-open every RESULTS file to refresh "has" and the parameter list');

c.srcLbl = uilabel(s1,'Text','(no files)','FontColor',[0.35 0.35 0.35],'WordWrap','on');
c.srcLbl.Layout.Row = 2; c.srcLbl.Layout.Column = [1 7];

% ===================== 2 - the file list (left column) ====================
leftP = uipanel(outer,'Title','2 - Files (tick "use", edit a label, filter a column)', ...
    'FontWeight','bold','BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
leftP.Layout.Row = 2; leftP.Layout.Column = 1;
L = uigridlayout(leftP,[4 1],'RowHeight',{'fit','1.6x','fit','1x'}, ...
    'Padding',[6 6 6 6],'RowSpacing',6,'BackgroundColor','w');

fr = uigridlayout(L,[1 9],'RowHeight',{'fit'}, ...
    'ColumnWidth',{'fit','1x','fit','1x','fit','1x','fit','fit','fit'}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',5,'BackgroundColor','w');
uilabel(fr,'Text','Animal');
c.fAnimal = uidropdown(fr,'Items',{'(all)'},'ValueChangedFcn',@(~,~)refreshTable());
uilabel(fr,'Text','Type');
c.fType   = uidropdown(fr,'Items',{'(all)'},'ValueChangedFcn',@(~,~)refreshTable());
uilabel(fr,'Text','Group');
c.fGroup  = uidropdown(fr,'Items',{'(all)'},'ValueChangedFcn',@(~,~)refreshTable());
uibutton(fr,'Text','Use sel.','ButtonPushedFcn',@(~,~)useSelected(true), ...
    'Tooltip','tick "use" on the rows selected in the table');
uibutton(fr,'Text','Drop sel.','ButtonPushedFcn',@(~,~)useSelected(false), ...
    'Tooltip','untick "use" on the rows selected in the table');
uibutton(fr,'Text','Remove','ButtonPushedFcn',@(~,~)removeSelected(), ...
    'Tooltip','remove the selected rows from the list');

c.fileTbl = uitable(L,'Data',emptyTableData(), ...
    'ColumnName',{'use','file','animal','type','group','has'}, ...
    'ColumnEditable',[true false true true true false], ...
    'ColumnWidth',{42,'auto',80,70,80,'auto'}, ...
    'Multiselect','on','SelectionType','row', ...
    'CellEditCallback',@onCellEdit);

c.previewLbl = uilabel(L,'Text','Preview: press "Preview" to see what would be written.', ...
    'FontColor',[0.30 0.30 0.30],'WordWrap','on');
c.previewTbl = uitable(L,'Data',cell(0,1),'ColumnName',{''});

% ===================== right column: parameters / averaging / output =======
rightP = uipanel(outer,'BackgroundColor','w','BorderType','none','Scrollable','on');
rightP.Layout.Row = 2; rightP.Layout.Column = 2;
R = uigridlayout(rightP,[4 1],'RowHeight',{'1x','fit','fit','fit'}, ...
    'Padding',[0 0 0 0],'RowSpacing',8,'BackgroundColor','w','Scrollable','on');

% --- 3 - parameters ---
p3 = uipanel(R,'Title','3 - Parameters (sMetrics columns)','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45],'Scrollable','on');
g3 = uigridlayout(p3,[2 1],'RowHeight',{'fit','fit'},'Padding',[6 6 6 6], ...
    'RowSpacing',4,'BackgroundColor','w','Scrollable','on');   % the grid owns the scroll
b3 = uigridlayout(g3,[1 3],'RowHeight',{'fit'},'ColumnWidth',{'fit','fit','1x'}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',5,'BackgroundColor','w');
uibutton(b3,'Text','All','ButtonPushedFcn',@(~,~)setAllColumns(true));
uibutton(b3,'Text','None','ButtonPushedFcn',@(~,~)setAllColumns(false));
c.colLbl = uilabel(b3,'Text','(no files)','FontColor',[0.5 0.5 0.5],'FontSize',10,'WordWrap','on');
c.colBox = uigridlayout(g3,[1 1],'RowHeight',{'fit'},'Padding',[0 0 0 0], ...
    'RowSpacing',1,'BackgroundColor','w');

% --- 4 - averaging ---
p4 = uipanel(R,'Title','4 - Averaging','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
g4 = uigridlayout(p4,[3 2],'RowHeight',{'fit','fit','fit'},'ColumnWidth',{'fit','1x'}, ...
    'Padding',[6 6 6 6],'RowSpacing',5,'ColumnSpacing',6,'BackgroundColor','w');
uilabel(g4,'Text','Rows are');
c.groupBy = uidropdown(g4,'Items',{'per segment (none)','average over labels', ...
    'average over vessel type'},'ValueChangedFcn',@(~,~)onGroupBy(), ...
    'Tooltip',['per segment = one row per ROI (the classic sheets).  ' ...
               'labels = the area-weighted ROI averages.  ' ...
               'vessel type = artery / vein / ... from setVesselTypes.']);
c.weight = uicheckbox(g4,'Text','weight by segment area','Value',true, ...
    'Tooltip','weights are sMetrics.area - the weights the ROI sheets have always used');
c.weight.Layout.Row = 2; c.weight.Layout.Column = [1 2];
c.avgLbl = uilabel(g4,'Text','','FontColor',[0.5 0.5 0.5],'FontSize',10,'WordWrap','on');
c.avgLbl.Layout.Row = 3; c.avgLbl.Layout.Column = [1 2];

% --- 5 - output ---
p5 = uipanel(R,'Title','5 - Output','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
g5 = uigridlayout(p5,[4 3],'RowHeight',{'fit','fit','fit','fit'}, ...
    'ColumnWidth',{'fit','1x','fit'},'Padding',[6 6 6 6],'RowSpacing',5, ...
    'ColumnSpacing',6,'BackgroundColor','w');
uilabel(g5,'Text','Workbooks');
c.merge = uidropdown(g5,'Items',{'one per file','merged into one'}, ...
    'ValueChangedFcn',@(~,~)onMerge(), ...
    'Tooltip',['merged writes ONE workbook whose rows carry file / animal / type / ' ...
               'group - directly usable for statistics']);
c.merge.Layout.Column = [2 3];
uilabel(g5,'Text','Format');
c.fmt = uidropdown(g5,'Items',{'xlsx','xls'},'Value','xlsx');
c.fmt.Layout.Column = [2 3];
uilabel(g5,'Text','Destination');
c.outF = uieditfield(g5,'text','Value','','Enable','off', ...
    'Tooltip','the merged workbook (one-per-file mode writes beside each input)');
c.outBtn = uibutton(g5,'Text','...','Enable','off','ButtonPushedFcn',@onPickOut);
c.prevBtn = uibutton(g5,'Text','Preview','ButtonPushedFcn',@(~,~)previewTable());
c.prevBtn.Layout.Row = 4; c.prevBtn.Layout.Column = 1;
c.runBtn = uibutton(g5,'Text','Export','FontWeight','bold', ...
    'BackgroundColor',[0.82 0.92 0.82],'ButtonPushedFcn',@(~,~)doRun());
c.runBtn.Layout.Row = 4; c.runBtn.Layout.Column = [2 3];

% --- log ---
p6 = uipanel(R,'Title','Log','FontWeight','bold','BackgroundColor','w', ...
    'ForegroundColor',[0.15 0.25 0.45]);
g6 = uigridlayout(p6,[2 1],'RowHeight',{'fit',110},'Padding',[6 6 6 6], ...
    'RowSpacing',4,'BackgroundColor','w');
c.status = uilabel(g6,'Text','Ready.','FontColor',[0.30 0.30 0.30],'WordWrap','on');
c.log    = uitextarea(g6,'Value',{''},'Editable','off');

app.c = c;
onGroupBy(); onMerge(); refreshTable();

% ---- programmatic API (drives the same logic as the UI) -------------------
setappdata(fig,'exportAPI',struct( ...
    'loadFiles',   @(p) loadFiles(p), ...
    'loadSession', @(p) loadSessionFile(p), ...
    'scanFolder',  @(root,glob) scanFolder(root,glob), ...
    'files',       @apiFiles, ...
    'exportSet',   @apiExportSet, ...
    'setUse',      @(p,tf) setUse(p,tf), ...
    'setFilter',   @(axis,v) setFilter(axis,v), ...
    'filterValues',@(axis) axisValues(axis), ...
    'parameters',  @() parameterList(), ...
    'setColumns',  @(cols) setColumns(cols), ...
    'setOpts',     @(varargin) setOpts(varargin{:}), ...
    'opts',        @() buildOpts(), ...
    'previewTable',@(varargin) previewTable(varargin{:}), ...
    'run',         @() doRun(), ...
    'status',      @() c.status.Text, ...
    'log',         @() c.log.Value, ...
    'getApp',      @apiApp));

% ---- a leading argument is a SESSION: the hand-off from the workbench ------
if ~isempty(sessionArg)
    try
        loadSessionFile(sessionArg);
    catch ME
        setStatus(['Session could not be read: ' ME.message]);
        logLine(['SESSION ERROR ' ME.message]);
    end
end

if nargout>0, h = fig; end

%% ========================== API ACCESSORS ============================== %%
% NESTED, not anonymous: an anonymous handle would capture a SNAPSHOT of `app`
% taken when the window was built (an empty list), so every state read has to go
% through a nested function that shares the live workspace - guiExplore's getApp
% idiom.
    function f = apiFiles(),     f = app.items;                  end
    function f = apiExportSet(), f = app.items(exportIdx());     end
    function a = apiApp(),       a = app;                        end

%% ========================== INPUT ====================================== %%
    function onAddFiles(~,~)
        [f,p] = uigetfile({'*_r.mat','Processed results files (*_r.mat)'; ...
                           '*.mat','MATLAB files (*.mat)'}, ...
            'Select processed recordings to export','MultiSelect','on',app.lastRoot);
        if isequal(f,0), return; end
        if ischar(f), f = {f}; end
        app.lastRoot = p;
        loadFiles(cellfun(@(x) fullfile(p,x), f, 'UniformOutput', false));
    end

    function onScanFolder(~,~)
        d = uigetdir(app.lastRoot,'Select the folder that holds your processed results');
        if isequal(d,0), return; end
        app.lastRoot = d;
        scanFolder(d, c.globF.Value);
    end

    function onLoadSession(~,~)
        [f,p] = uigetfile({'*.mat','Workbench session (*.mat)'}, ...
            'Select a Processing Workbench session',app.lastRoot);
        if isequal(f,0), return; end
        app.lastRoot = p;
        try
            loadSessionFile(fullfile(p,f));
        catch ME
            uialert(fig,ME.message,'Session');
        end
    end

    function scanFolder(root, glob)
        %scanFolder  Recursive glob through the library's own file-list helper.
        if nargin<2 || isempty(glob), glob = c.globF.Value; end
        c.globF.Value = glob;
        setStatus(['Scanning ' root ' ...']);
        paths = {};
        try
            paths = getFileNamesList(root, glob);        % recursive, flat list
        catch ME
            setStatus(['Scan failed: ' ME.message]);
        end
        if isempty(paths), setStatus('Nothing matched the glob.'); return; end
        loadFiles(paths(:)');
    end

    function loadFiles(paths)
        %loadFiles  Resolve each path to the exportable products of its recording
        %   and append the new ones (labels left empty - they are the session's).
        if isempty(paths), return; end
        if ischar(paths) || isstring(paths), paths = cellstr(paths); end
        addItems(paths, repmat(emptyLabels(),1,numel(paths)));
    end

    function loadSessionFile(pth)
        %loadSessionFile  THE hand-off from the workbench (spec §5).  Read through
        %   wbSession('read',...) - never by reaching into the file layout - so the
        %   curated list arrives with its animal / type / experimental group.
        S = wbSession('read', char(pth));
        app.sessionPath = char(pth);
        if ~isempty(S.files)
            paths = {S.files.path};
            lab   = repmat(emptyLabels(),1,numel(paths));
            for is = 1:numel(S.files)
                lab(is).animal = charOf(S.files(is).animal);
                lab(is).type   = charOf(S.files(is).type);
                lab(is).group  = charOf(S.files(is).expGroup);
            end
        else                                            % pre-schema-5 sidecar
            paths = S.paths(:)';
            lab   = repmat(emptyLabels(),1,numel(paths));
        end
        addItems(paths, lab);
        logLine(sprintf('session %s: %d file(s) -> %d exportable product(s)', ...
            shortName(pth), numel(paths), numel(app.items)));
    end

    function addItems(paths, lab)
        %addItems  Expand every input path to its exportable *_r.mat products, drop
        %   duplicates, probe each RESULTS file once, then rebuild the whole UI.
        n0 = numel(app.items);
        have = containers.Map('KeyType','char','ValueType','logical');
        for ia = 1:numel(app.items), have(app.items(ia).path) = true; end
        for ia = 1:numel(paths)
            prods = resolveProducts(paths{ia});
            for k = 1:numel(prods)
                if isKey(have, prods{k}), continue; end
                have(prods{k}) = true;
                app.items(end+1) = makeItem(prods{k}, lab(min(ia,numel(lab))));
            end
        end
        probeAll();
        refreshTable();
        n = numel(app.items);
        setStatus(sprintf('%d exportable recording(s) (%d new).', n, n-n0));
        if n==n0
            logLine('no new exportable products (a recording needs an _r.mat results file)');
        end
    end

    function it = makeItem(rPath, lab)
        %makeItem  One row of the list.  'path' is what exportToExcel is handed and
        %   'rpath' is what the probe reads - the SAME file since the library refers
        %   to a product by its RESULTS member; the field is kept because the probe
        %   cache and the session sidecar are both keyed by it.
        [~,nm,ex] = fileparts(rPath);
        it = struct('path',rPath,'rpath',rPath, ...
            'name',[nm ex],'animal',lab.animal,'type',lab.type,'group',lab.group, ...
            'use',true,'cols',{{}},'has','','note','');
    end

    function probeAll()
        %probeAll  Open each RESULTS file once (cached) to learn which sMetrics
        %   columns and which result trees it carries - the parameter picker's
        %   availability counts and the table's "has" column come from here.
        for ib = 1:numel(app.items)
            setStatus(sprintf('Reading results %d/%d ...', ib, numel(app.items)));
            drawnow limitrate
            P = probeFor(app.items(ib).rpath);
            app.items(ib).cols = P.cols;
            app.items(ib).has  = strjoin(P.trees, ', ');
            app.items(ib).note = P.note;
        end
    end

    function P = probeFor(rpath)
        if isKey(app.probe, rpath), P = app.probe(rpath); return; end
        P = probeResults(rpath);
        app.probe(rpath) = P;
    end

    function reprobeAll()
        app.probe = containers.Map('KeyType','char','ValueType','any');
        probeAll(); refreshTable();
        setStatus(sprintf('Re-read %d results file(s).', numel(app.items)));
    end

    function clearList()
        app.items = emptyItems(); app.sessionPath = '';
        refreshTable(); setStatus('List cleared.');
    end

%% ====================== THE FILE TABLE ================================= %%
    function refreshTable()
        %refreshTable  Rebuild the filter vocabularies FROM THE DATA, apply them,
        %   redraw the table, then rebuild the parameter picker for what is left.
        for a = 1:numel(AXES), refreshFilter(AXES{a}); end
        app.view = viewIdx();
        c.fileTbl.Data = tableData(app.view);
        refreshColumns();
        updateCounts();
    end

    function refreshFilter(axisName)
        d = filterCtrl(axisName);
        vals = axisValues(axisName);
        items = [{'(all)'} vals];
        keep = d.Value;
        d.Items = items;
        if any(strcmp(keep,items)), d.Value = keep; else, d.Value = '(all)'; end
    end

    function v = axisValues(axisName)
        %axisValues  The values this axis actually takes in the loaded set.  Built
        %   from the data every time: the vocabularies are the user's, open and
        %   never validated against a list (spec §2 / D1).
        v = {};
        if isempty(app.items), return; end
        raw = cellfun(@(x) labelOrNone(x), {app.items.(axisName)}, 'UniformOutput', false);
        v = unique(raw,'stable'); v = sort(v);
    end

    function idx = viewIdx()
        idx = [];
        for iv = 1:numel(app.items)
            if passesFilters(app.items(iv)), idx(end+1) = iv; end %#ok<AGROW>
        end
    end

    function tf = passesFilters(it)
        tf = true;
        for a = 1:numel(AXES)
            want = filterCtrl(AXES{a}).Value;
            if strcmp(want,'(all)'), continue; end
            if ~strcmp(labelOrNone(it.(AXES{a})), want), tf = false; return; end
        end
    end

    function d = filterCtrl(axisName)
        switch axisName
            case 'animal', d = c.fAnimal;
            case 'type',   d = c.fType;
            otherwise,     d = c.fGroup;
        end
    end

    function D = tableData(idx)
        if isempty(idx), D = emptyTableData(); return; end
        it = app.items(idx);
        D = table(logical([it.use]'), {it.name}', {it.animal}', {it.type}', ...
            {it.group}', {it.has}', 'VariableNames', ...
            {'use','file','animal','type','group','has'});
    end

    function onCellEdit(~,ev)
        r = ev.Indices(1); col = ev.Indices(2);
        if r>numel(app.view), return; end
        i = app.view(r);
        switch col
            case 1, app.items(i).use = logical(ev.NewData);
            case 3, app.items(i).animal = charOf(ev.NewData);
            case 4, app.items(i).type   = charOf(ev.NewData);
            case 5, app.items(i).group  = charOf(ev.NewData);
        end
        refreshTable();
    end

    function useSelected(tf)
        r = c.fileTbl.Selection;
        if isempty(r), setStatus('Select rows in the table first.'); return; end
        r = unique(r(:,1));
        for k = 1:numel(r)
            if r(k)<=numel(app.view), app.items(app.view(r(k))).use = tf; end
        end
        refreshTable();
    end

    function removeSelected()
        r = c.fileTbl.Selection;
        if isempty(r), setStatus('Select rows in the table first.'); return; end
        kill = app.view(unique(r(:,1)));
        app.items(kill) = [];
        refreshTable(); setStatus(sprintf('Removed %d row(s).', numel(kill)));
    end

    function setUse(paths, tf)
        if ischar(paths) || isstring(paths), paths = cellstr(paths); end
        for iu = 1:numel(app.items)
            if any(strcmp(app.items(iu).path, paths)), app.items(iu).use = logical(tf); end
        end
        refreshTable();
    end

    function setFilter(axisName, v)
        d = filterCtrl(axisName);
        if ~any(strcmp(v,d.Items)), d.Items = [d.Items {char(v)}]; end
        d.Value = char(v);
        refreshTable();
    end

    function idx = exportIdx()
        %exportIdx  What "Export" will run over: visible under the filters AND used.
        idx = app.view(arrayfun(@(i) app.items(i).use, app.view));
    end

%% ====================== THE PARAMETER PICKER =========================== %%
    function P = parameterList()
        %parameterList  Union of the sMetrics columns of the export set, each with
        %   how many of those files actually carry it - availability made visible
        %   before anything is written.
        P = struct('name',{},'count',{},'total',{},'on',{});
        idx = exportIdx();
        if isempty(idx), return; end
        names = {};
        for ip = idx, names = union(names, app.items(ip).cols, 'stable'); end
        names = sort(names);
        for k = 1:numel(names)
            n = 0;
            for ip = idx, n = n + any(strcmp(names{k}, app.items(ip).cols)); end
            P(end+1) = struct('name',names{k},'count',n,'total',numel(idx), ...
                'on',columnOn(names{k})); %#ok<AGROW>
        end
    end

    function tf = columnOn(name)
        if isKey(app.colOn,name), tf = app.colOn(name); else, tf = true; end
    end

    function refreshColumns()
        %refreshColumns  Redraw the checkbox list for the current export set.
        P = parameterList();
        delete(c.colBox.Children);
        c.colBox.RowHeight = repmat({'fit'},1,max(numel(P),1));
        app.colChk = struct();
        for k = 1:numel(P)
            cb = uicheckbox(c.colBox,'Text',sprintf('%s  -  %d/%d files', ...
                P(k).name, P(k).count, P(k).total), 'Value',P(k).on, ...
                'Tooltip','untick to leave this sMetrics column out of the workbook', ...
                'ValueChangedFcn',@(src,~)onColumnTick(P(k).name,src.Value));
            app.colChk.(matlab.lang.makeValidName(P(k).name)) = cb;
        end
        if isempty(P)
            c.colLbl.Text = 'no sMetrics columns in the current selection';
        else
            c.colLbl.Text = sprintf('%d column(s) offered', numel(P));
        end
    end

    function onColumnTick(name, tf)
        app.colOn(name) = logical(tf);
        updateCounts();
    end

    function setAllColumns(tf)
        P = parameterList();
        for k = 1:numel(P), app.colOn(P(k).name) = logical(tf); end
        refreshColumns(); updateCounts();
    end

    function setColumns(cols)
        %setColumns  Tick exactly these columns (programmatic; [] = all).
        P = parameterList();
        if isempty(cols)
            for k = 1:numel(P), app.colOn(P(k).name) = true; end
        else
            cols = cellstr(string(cols));
            for k = 1:numel(P), app.colOn(P(k).name) = any(strcmp(P(k).name,cols)); end
        end
        refreshColumns();
    end

    function cols = tickedColumns()
        %tickedColumns  The .columns value for exportToExcel: {} when everything on
        %   offer is ticked, so a column we never saw (an unprobed large file) is
        %   not silently dropped.
        P = parameterList();
        cols = {};
        if isempty(P), return; end
        on = [P.on];
        if all(on), return; end
        cols = {P(on).name};
    end

%% ====================== AVERAGING / OUTPUT ============================= %%
    function onGroupBy()
        g = groupByValue();
        switch g
            case 'label'
                c.weight.Value = true; c.weight.Enable = 'off';
                c.avgLbl.Text = ['Area-weighted averages over results.sMetrics.label - ' ...
                    'the sMetricsROI / sDataROI sheets, unchanged.'];
            case 'type'
                c.weight.Enable = 'on';
                c.avgLbl.Text = ['Averages over the VESSEL type (results.sMetrics.type, ' ...
                    'from setVesselTypes) - not the recording type of the file list.'];
            otherwise
                c.weight.Enable = 'off';
                c.avgLbl.Text = 'One row per segment: the full nine-sheet workbook.';
        end
        updateCounts();
    end

    function g = groupByValue()
        switch c.groupBy.Value
            case 'average over labels',      g = 'label';
            case 'average over vessel type', g = 'type';
            otherwise,                       g = '';
        end
    end

    function onMerge()
        on = mergeOn();
        c.outF.Enable = onOff(on); c.outBtn.Enable = onOff(on);
        if on && isempty(c.outF.Value), c.outF.Value = defaultOut(); end
        updateCounts();
    end

    function tf = mergeOn(), tf = strcmp(c.merge.Value,'merged into one'); end

    function o = defaultOut()
        idx = exportIdx(); o = '';
        if isempty(idx), return; end
        o = fullfile(fileparts(app.items(idx(1)).path), ['mergedExport.' c.fmt.Value]);
    end

    function onPickOut(~,~)
        start = c.outF.Value; if isempty(start), start = defaultOut(); end
        [f,p] = uiputfile({['*.' c.fmt.Value],'Excel workbook'},'Merged workbook',start);
        if isequal(f,0), return; end
        c.outF.Value = fullfile(p,f);
    end

    function O = buildOpts()
        %buildOpts  THE only thing this GUI decides: the opts struct.  Every number
        %   in the workbook is computed by exportToExcel from it.
        O = struct();
        O.format       = ['.' c.fmt.Value];
        O.groupBy      = groupByValue();
        O.weightByArea = c.weight.Value;
        cols = tickedColumns();
        if ~isempty(cols), O.columns = cols; end
        if mergeOn()
            O.merge   = true;
            O.outFile = c.outF.Value;
            if isempty(O.outFile), O.outFile = defaultOut(); end
            O.labels  = labelStruct(exportIdx());
        end
    end

    function L = labelStruct(idx)
        L = struct('animal',{},'type',{},'group',{});
        for il = idx
            L(end+1) = struct('animal',app.items(il).animal, ...
                'type',app.items(il).type,'group',app.items(il).group); %#ok<AGROW>
        end
    end

    function setOpts(varargin)
        %setOpts  Programmatic mirror of the averaging / output controls.
        for k = 1:2:numel(varargin)
            v = varargin{k+1};
            switch lower(varargin{k})
                case 'groupby'
                    switch char(string(v))
                        case 'label', c.groupBy.Value = 'average over labels';
                        case 'type',  c.groupBy.Value = 'average over vessel type';
                        otherwise,    c.groupBy.Value = 'per segment (none)';
                    end
                    onGroupBy();
                case 'weightbyarea', c.weight.Value = logical(v);
                case 'columns',      setColumns(v);
                case 'merge'
                    if logical(v), c.merge.Value = 'merged into one';
                    else,          c.merge.Value = 'one per file'; end
                    onMerge();
                case 'outfile',      c.outF.Value = char(string(v));
                case 'format'
                    e = char(string(v)); if ~isempty(e) && e(1)=='.', e = e(2:end); end
                    c.fmt.Value = e;
                case 'previewrows',  app.previewRows = double(v);
            end
        end
        updateCounts();
    end

%% ====================== PREVIEW / RUN ================================== %%
    function T = previewTable(nRows)
        %previewTable  Build the preview AND show it, so a programmatic call and the
        %   button see exactly the same thing.
        if nargin<1 || isempty(nRows), nRows = app.previewRows; end
        T = buildPreview(nRows);
        renderPreview(T);
    end

    function T = buildPreview(nRows)
        %buildPreview  What WOULD be written, for the first file of the export set -
        %   produced by exportToExcel itself into a scratch workbook, so the preview
        %   can never disagree with the result.  The merged row prefix is dropped
        %   again when the user is exporting one workbook per file.
        if nargin<1 || isempty(nRows), nRows = app.previewRows; end
        T = table();
        idx = exportIdx();
        if isempty(idx), setStatus('Nothing selected to preview.'); return; end
        O = buildOpts();
        work = [tempname '_guiExportPreview']; mkdir(work);
        cleaner = onCleanup(@() rmdirSafe(work));
        O.merge   = true;                                  % the only way to redirect
        O.outFile = fullfile(work,['preview.' c.fmt.Value]);
        O.labels  = labelStruct(idx(1));
        try
            exportToExcel({app.items(idx(1)).path}, O);
        catch ME
            setStatus(['Preview failed: ' ME.message]); logLine(['PREVIEW ERROR ' ME.message]);
            return
        end
        if ~isfile(O.outFile), setStatus('Preview produced no sheet.'); return; end
        sh = sheetnames(O.outFile);
        T = readtable(O.outFile,'Sheet',sh(1),'VariableNamingRule','preserve');
        if ~mergeOn()
            % undo the merged framing so the preview shows the per-file schema:
            % the label prefix goes, and the vessel type gets its own name back
            v = T.Properties.VariableNames;
            T = T(:, v(~ismember(v,{'file','animal','type','group'})));
            v = T.Properties.VariableNames;
            k = strcmp(v,'vesselType');
            if any(k) && ~any(strcmp(v,'type'))
                v{k} = 'type'; T.Properties.VariableNames = v;
            end
        end
        if height(T) > nRows, T = T(1:nRows,:); end
        setStatus(sprintf('Preview: sheet "%s" of %s.', sh(1), app.items(idx(1)).name));
    end

    function renderPreview(T)
        if isempty(T)
            c.previewTbl.Data = cell(0,1); c.previewTbl.ColumnName = {''};
            c.previewLbl.Text = 'Preview: nothing to show.';
            return
        end
        c.previewTbl.Data = T;
        % explicit: the table's own names, because ColumnName was set once at build
        c.previewTbl.ColumnName = T.Properties.VariableNames;
        c.previewLbl.Text = sprintf(['Preview of the FIRST selected file - %d row(s) of ' ...
            'the first sheet that would be written.'], height(T));
    end

    function paths = doRun()
        %doRun  Call exportToExcel over the export set: once per file (one workbook
        %   each, with a progress line) or once for the whole set (merged).
        paths = {};
        idx = exportIdx();
        if isempty(idx), setStatus('Nothing selected to export.'); return; end
        O = buildOpts();
        files = {app.items(idx).path};
        logLine(sprintf('=== EXPORT: %d file(s), groupBy=''%s'', %s ===', numel(files), ...
            O.groupBy, tern(mergeOn(),'merged','one workbook per file')));
        if mergeOn()
            setStatus(sprintf('Writing the merged workbook (%d file(s)) ...', numel(files)));
            drawnow limitrate
            try
                exportToExcel(files, O);
                paths = {O.outFile};
                logLine(['  wrote ' O.outFile]);
                setStatus(['Wrote ' O.outFile]);
            catch ME
                logLine(['  ERROR ' ME.message]);
                setStatus(['Export failed: ' ME.message]);
            end
            return
        end
        nok = 0; nfail = 0;
        for k = 1:numel(files)
            setStatus(sprintf('Exporting %d/%d: %s', k, numel(files), shortName(files{k})));
            drawnow limitrate
            try
                exportToExcel(files(k), O);
                out = strrep(files{k},'_r.mat',['.' c.fmt.Value]);
                paths{end+1} = out; %#ok<AGROW>
                nok = nok + 1;
                logLine(['  wrote ' out]);
            catch ME
                nfail = nfail + 1;
                logLine(sprintf('  ERROR %s - %s', shortName(files{k}), ME.message));
            end
        end
        logLine(sprintf('=== EXPORT complete: %d ok, %d failed ===', nok, nfail));
        setStatus(sprintf('Exported %d workbook(s), %d failed.', nok, nfail));
    end

%% ====================== small nested helpers =========================== %%
    function updateCounts()
        n = numel(exportIdx());
        if isempty(app.items)
            c.srcLbl.Text = '(no files)';
        else
            c.srcLbl.Text = sprintf(['%d exportable product(s) loaded, %d selected for ' ...
                'export%s.'], numel(app.items), n, sessionNote());
        end
        if mergeOn() && isempty(c.outF.Value), c.outF.Value = defaultOut(); end
    end

    function s = sessionNote()
        if isempty(app.sessionPath), s = '';
        else, s = ['   [session: ' shortName(app.sessionPath) ']']; end
    end

    function setStatus(msg), c.status.Text = msg; drawnow limitrate; end

    function logLine(msg)
        v = c.log.Value;
        if isscalar(v) && isempty(v{1}), v = {}; end
        v{end+1} = msg;
        c.log.Value = v;
    end
end

%% ========================= LOCAL HELPERS ================================ %%
function [sess, vis] = parseArgs(args)
%parseArgs  guiExport(sessionPath) / guiExport('Visible',v) / both.  A LEADING char
%   that is not an option name is the session path - the shape guiWorkbench hands
%   over, and the same shape guiExplore takes.
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

function it = emptyItems()
%emptyItems  The 0x0 shape of the file list (one entry per exportable product).
it = struct('path',{},'rpath',{},'name',{},'animal',{},'type',{},'group',{}, ...
    'use',{},'cols',{},'has',{},'note',{});
end

function L = emptyLabels()
L = struct('animal','','type','','group','');
end

function D = emptyTableData()
D = table(false(0,1),cell(0,1),cell(0,1),cell(0,1),cell(0,1),cell(0,1), ...
    'VariableNames',{'use','file','animal','type','group','has'});
end

function prods = resolveProducts(pth)
%resolveProducts  The exportable products of ONE input path.  exportToExcel takes
%   the RESULTS member of a triplet ('*_r.mat'), so a product is exportable exactly
%   when it has one.  A raw recording (or any other member of the same recording)
%   expands to EVERY product of that recording - '_t_BFI', '_c_BFI', ... - which is
%   how a session of .rls files becomes a list of workbooks.  Identity comes from
%   wbFileModel, so a 'Roi2_' crop stays distinct.
prods = {};
pth = char(pth);
if isempty(pth), return; end
if endsWith(pth,'_r.mat') && isfile(pth)
    prods = {pth};
    return
end
m = wbFileModel(pth);
if isempty(m.folder) || ~isfolder(m.folder), return; end
d = dir(fullfile(m.folder,[m.roiPrefix m.stem '*_r.mat']));
for k = 1:numel(d)
    p = fullfile(d(k).folder,d(k).name);
    cm = wbFileModel(p);
    if ~strcmp(cm.identity, m.identity), continue; end
    prods{end+1} = p; %#ok<AGROW>
end
end

function P = probeResults(rpath)
%probeResults  What one RESULTS file offers: its sMetrics column names and which
%   result trees it carries.  A very large file (a per-pixel vasomotion result) is
%   NOT loaded just to fill a table cell - it is reported as unprobed and simply
%   contributes no columns to the availability counts.
P = struct('ok',false,'cols',{{}},'trees',{{}},'note','');
if isempty(rpath) || ~isfile(rpath), P.note = 'no _r.mat results file'; return; end
d = dir(rpath);
if ~isempty(d) && d(1).bytes > 1.5e9
    P.trees = {'(large - not read)'}; P.note = 'file too large to probe'; return
end
try
    S = load(rpath,'results');
catch ME
    P.note = ME.message; return
end
if ~isfield(S,'results'), P.note = 'no results struct'; return; end
R = S.results;
t = {};
if isfield(R,'sMetrics') && istable(R.sMetrics)
    P.cols = R.sMetrics.Properties.VariableNames;
    t{end+1} = 'sMetrics';
    if isfield(R,'sData'), t{end+1} = 'sData'; end
    if ismember('label',P.cols), t{end+1} = 'labels'; end
    if ismember('type', P.cols), t{end+1} = 'vessel types'; end
end
if isfield(R,'dvsMetrics'),  t{end+1} = 'dvs';         end
if isfield(R,'pulsatility'), t{end+1} = 'pulsatility'; end
if isfield(R,'vasomotion'),  t{end+1} = 'vasomotion';  end
if isfield(R,'CTTH'),        t{end+1} = 'CTTH';        end
% A myograph recording keeps its results inside its analysed WINDOWS rather than in
% a per-segment table, so what it offers is read from there.  myographIntervals is
% what knows the two shapes those windows are stored in - flat for a pressure
% myograph, split by channel for a wire one.  This is a probe, not a branch:
% exportToExcel decides which sheets a file writes, and this only says what the
% "has" column and the availability counts should show for it.
iv = myographIntervals(R);
if ~isempty(iv)
    t{end+1} = sprintf('%d intervals', numel(iv));
    for b = {'diameter','propagation','vasomotion'}
        if anyInterval(iv, b{1}), t{end+1} = b{1}; end %#ok<AGROW>
    end
end
P.trees = unique(t,'stable');
P.ok = true;
end

function tf = anyInterval(iv, branch)
%anyInterval  Did ANY analysed window get this branch written into it?
tf = false;
if ~isfield(iv,branch), return; end
for i = 1:numel(iv)
    if ~isempty(iv(i).(branch)), tf = true; return; end
end
end

function s = labelOrNone(v)
%labelOrNone  A label axis value for display / filtering; '(none)' when unset.
s = charOf(v);
if isempty(s), s = '(none)'; end
end

function s = charOf(v)
if ischar(v), s = v; elseif isstring(v), s = char(v);
elseif isnumeric(v) || islogical(v), s = num2str(v); else, s = ''; end
end

function s = shortName(p)
[~,n,e] = fileparts(char(p)); s = [n e];
end

function s = onOff(tf), if tf, s = 'on'; else, s = 'off'; end, end
function s = tern(tf,a,b), if tf, s = a; else, s = b; end, end

function rmdirSafe(d)
try
    if exist(d,'dir'), rmdir(d,'s'); end
catch
end
end
