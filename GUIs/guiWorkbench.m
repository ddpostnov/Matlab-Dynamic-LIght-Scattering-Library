%guiWorkbench - Processing Workbench: a file x step matrix for the LSCI pipeline.
%
%   A programmatic uifigure that turns the launcher workflow into a spreadsheet:
%   ROWS are recordings (grouped by animal, reference first), COLUMNS are the
%   pipeline steps from wbStepRegistry, and each CELL shows - and lets you queue -
%   the state of one step for one recording (ready / done / stale / unavailable).
%   The window is a thin controller: every decision is delegated to the headless
%   Phase-2 brain (wbStepRegistry, wbDiscoverFiles, wbFileModel, wbStateEngine,
%   wbSettingsModel, wbInvalidate), to wbExecutor / wbModalGuard / wbArtifacts,
%   and to wbSession.  Run executes every checked cell in dependency order
%   through the real wrappers (wbExecutor), streaming per-cell progress and a
%   command-window mirror into the log, with cooperative Stop and continue-on-
%   error; Preview order lists the same plan without calling anything.  Finished
%   cells turn into clickable 'done' buttons that open their report images.
%
%   Tabs: 1 Files (the three loaders + a grouped, sortable file list) *
%         2 Process (the matrix + per-step settings + reports + a log pane) *
%         3 Export  (a sheet/format selection UI over exportToExcel) *
%         4 Explore (GUIs/workbench/guiExplore hosted in-tab, seeded from the
%                    loaded files/groups) *
%
% Syntax:
%    guiWorkbench                       % open the workbench
%    h = guiWorkbench('Visible','off')  % headless (for programmatic drive/tests)
%
%   getappdata(h,'workbenchAPI') exposes a struct of function handles that drive
%   the same internal logic as the UI (used by testWorkbenchShell).
%
% See also: wbStepRegistry, wbDiscoverFiles, wbStateEngine, wbSettingsModel,
%           wbInvalidate, wbExecutor, wbModalGuard, wbArtifacts, wbSession,
%           guiExplore, guiMyograph
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function h = guiWorkbench(varargin)

vis = 'on';
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i},'Visible'), vis = varargin{i+1}; end
end

app = struct();
app.reg          = wbStepRegistry('LSCI');       % v1 columns (the LSCI steps)
app.modality     = 'LSCI';
app.sm           = wbSettingsModel('new');        % settings bag + overrides
app.disc         = emptyDisc();
app.rows         = emptyRows();                    % flat, group-major, reference-first
app.groupNames   = {};
app.modelArr     = wbFileModel('x.rls'); app.modelArr(1) = [];   % 1x0 model array
app.base         = containers.Map('KeyType','char','ValueType','any');   % identity -> struct(stepId->state)
app.checked      = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> true
app.stale        = containers.Map('KeyType','char','ValueType','any');   % session edit overlay
app.runState     = containers.Map('KeyType','char','ValueType','any');   % transient run overlay: 'running'|'done'|'error'
app.cellMsg      = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> error/tooltip text
app.rendered     = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> last drawn state
app.cellComp     = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> handle
app.gridRowOf    = containers.Map('KeyType','char','ValueType','double');% identity -> matrix grid row
app.selStep      = app.reg(1).id;                  % step shown in the settings panel
app.presetDir    = presetDirDefault();
app.presetRef    = '';
app.maxWidgets   = 800;                            % >this many rows -> uitable fallback
app.running      = false;
app.cancel       = false;
app.exportSrcs   = struct('path',{},'rpath',{},'group',{},'label',{});  % exportable BFI sources
app.exploreAPI   = [];                             % handle bundle of the hosted guiExplore
app.exploreHosted= false;                          % true = in-tab, false = child window / none
app.exploreSeedKey = '';                           % file-set fingerprint of the last seed

% single instance: replace any stale window (by Tag AND Name so older,
% Tag-less windows are cleared too)
delete(findall(groot,'Type','figure','Tag','guiWorkbench'));
delete(findall(groot,'Type','figure','Name','Processing Workbench'));
fig = uifigure('Name','Processing Workbench','Position',[50 50 1500 880], ...
    'Visible',vis,'Tag','guiWorkbench');
app.fig = fig;
g  = uigridlayout(fig,[1 1],'Padding',[4 4 4 4]);
tg = uitabgroup(g);
app.tabs.files   = uitab(tg,'Title','1 - Files');
app.tabs.process = uitab(tg,'Title','2 - Process');
app.tabs.export  = uitab(tg,'Title','3 - Export');
app.tabs.explore = uitab(tg,'Title','4 - Explore');
app.tg = tg;
setappdata(fig,'app',app);

buildFilesTab(fig);
buildProcessTab(fig);
buildExportTab(fig);
buildExploreTab(fig);
tg.SelectionChangedFcn = @(~,~) onTabSelected(fig);   % lazily seed Export/Explore on show

fig.CloseRequestFcn = @(~,~) requestExit(fig);

% ---- programmatic API (drives the same internal logic as the UI) ----
api = struct( ...
    'load',        @(varargin) apiLoad(fig,varargin{:}), ...
    'columns',     @() {getApp(fig).reg.id}, ...
    'rowCount',    @() numel(getApp(fig).rows), ...
    'groups',      @() getApp(fig).groupNames, ...
    'cellState',   @(id,stepId) resolveCellState(getApp(fig),id,stepId), ...
    'check',       @(id,stepId,tf) apiCheck(fig,id,stepId,tf), ...
    'checkColumn', @(stepId,tf) checkColumn(fig,stepId,tf), ...
    'checkGroup',  @(g,tf) checkGroup(fig,g,tf), ...
    'checkModality',@(m,tf) checkModality(fig,m,tf), ...
    'checkAll',    @(tf) checkAll(fig,tf), ...
    'checkedList', @() checkedList(fig), ...
    'editSetting', @(stepId,field,value) onSettingEdit(fig,stepId,field,value), ...
    'getSetting',  @(stepId,field) getSetting(fig,stepId,field), ...
    'selectStep',  @(stepId) selectStep(fig,stepId), ...
    'savePreset',  @(p) savePreset(fig,p), ...
    'loadPreset',  @(p) loadPreset(fig,p), ...
    'saveSession', @(p) saveSessionTo(fig,p), ...
    'loadSession', @(p) loadSessionFrom(fig,p), ...
    'runOrder',    @() buildRunOrder(fig), ...
    'dryRun',      @() dryRun(fig), ...
    'run',         @() runChecked(fig), ...
    'log',         @() getLog(fig), ...
    'exportRefresh',@() refreshExportFiles(fig), ...
    'exportSources',@() exportSourcePaths(fig), ...
    'exportAvailableSheets',@(rpaths) availableExportSheets(rpaths), ...
    'runExport',   @(files,opts) runExport(fig,files,opts), ...
    'seedExplore', @() seedExplore(fig,true), ...
    'exploreApi',  @() getApp(fig).exploreAPI, ...
    'identities',  @() {getApp(fig).rows.identity}, ...
    'getApp',      @() getApp(fig), ...
    'exit',        @() onClose(fig));
setappdata(fig,'workbenchAPI',api);

renderMatrix(fig);
selectStep(fig,app.selStep);
if nargout>0, h = fig; end
end

%% ===================== app-state helpers ============================ %%
function app = getApp(fig), app = getappdata(fig,'app'); end
function setApp(fig,app), setappdata(fig,'app',app); end
function d = emptyDisc()
d = struct('fNames',{{}},'models',{{}},'flat',[],'groups',[],'referenceMode',false);
end
function r = emptyRows()
r = struct('model',{},'identity',{},'group',{},'groupIdx',{},'rowInGroup',{},'isRef',{},'label',{});
end
function d = presetDirDefault()
d = fullfile(prefdir,'guiWorkbenchPresets');
end
function s = stepById(reg,id), s = reg(strcmp({reg.id},id)); end
function k = cellKey(identity,stepId), k = [identity '||' stepId]; end

%% ===================== FILES tab ==================================== %%
function buildFilesTab(fig)
app = getApp(fig); t = app.tabs.files;
gl = uigridlayout(t,[1 1],'Padding',[6 6 6 6]);
lp = uigridlayout(gl,[7 1],'RowHeight',{'fit','fit','fit','fit','1x','fit','fit'},'RowSpacing',6);

uilabel(lp,'Text','Load recordings (all three loaders group by animal, reference first):','FontWeight','bold');

% -- structured loader row (getFileNamesList) --
sp = uigridlayout(lp,[1 6],'ColumnWidth',{'fit','1.6x','fit','1x','fit','fit'},'Padding',[0 0 0 0]);
uilabel(sp,'Text','root');
c.root = uieditfield(sp,'text','Value','','Tooltip','folder searched recursively');
uilabel(sp,'Text','glob / animal / ref');
c.glob = uieditfield(sp,'text','Value','*_t_K_d.mat','Tooltip','dir() glob, e.g. *_t_K_d.mat');
c.animal = uieditfield(sp,'text','Value','[A-Z]+\d+','Tooltip','animal-id regexp (empty = one group)');
c.ref = uieditfield(sp,'text','Value','','Tooltip','reference regexp -> pinned to column 1 (empty = none)');

% -- loader buttons --
bp = uigridlayout(lp,[1 4],'Padding',[0 0 0 0]);
uibutton(bp,'Text','Load structured','ButtonPushedFcn',@(~,~)uiLoadStructured(fig), ...
    'Tooltip','recurse root for the glob, group by animal, pin the reference (getFileNamesList)');
uibutton(bp,'Text','Load folder...','ButtonPushedFcn',@(~,~)uiLoadFolder(fig), ...
    'Tooltip','pick a folder, recurse for the glob, group by the animal regexp');
uibutton(bp,'Text','Add files...','ButtonPushedFcn',@(~,~)uiLoadManual(fig), ...
    'Tooltip','pick files by hand (multiselect), group by the animal regexp');
uibutton(bp,'Text','Clear','ButtonPushedFcn',@(~,~)uiClear(fig));

uilabel(lp,'Text','Loaded files (grouped, reference first; column headers are sortable - a view only, execution walks group -> reference -> row):','WordWrap','on');
c.fileTbl = uitable(lp,'ColumnName',{'group','file','ref','modality','stage'}, ...
    'ColumnSortable',[true true true true true],'RowName',{});

% -- session row --
ssp = uigridlayout(lp,[1 4],'Padding',[0 0 0 0]);
uibutton(ssp,'Text','Save session...','ButtonPushedFcn',@(~,~)uiSaveSession(fig));
uibutton(ssp,'Text','Load session...','ButtonPushedFcn',@(~,~)uiLoadSession(fig));
c.status = uilabel(ssp,'Text','No files loaded.');
uibutton(ssp,'Text','Exit workbench','ButtonPushedFcn',@(~,~)requestExit(fig),'BackgroundColor',[1 0.82 0.82]);

app.c.files = c; setApp(fig,app);
end

%% ===================== PROCESS tab ================================== %%
function buildProcessTab(fig)
app = getApp(fig); t = app.tabs.process;
gl = uigridlayout(t,[1 2],'ColumnWidth',{'2.4x','1x'},'Padding',[6 6 6 6],'ColumnSpacing',8);

% -- LEFT: toolbar + matrix + log --
left = uigridlayout(gl,[3 1],'RowHeight',{'fit','1x','fit'},'RowSpacing',6);
tb = uigridlayout(left,[1 9], ...
    'ColumnWidth',{'fit','fit','fit','fit','fit','1x','fit','fit','fit'},'Padding',[0 0 0 0]);
uibutton(tb,'Text','Check all','ButtonPushedFcn',@(~,~)checkAll(fig,true));
uibutton(tb,'Text','Clear checks','ButtonPushedFcn',@(~,~)checkAll(fig,false));
c.presetDrop = uidropdown(tb,'Items',{'(launcher defaults)'},'Value','(launcher defaults)', ...
    'ValueChangedFcn',@(s,~)onPresetPick(fig,s.Value),'Tooltip','seed the settings bag from a saved preset');
uibutton(tb,'Text','Save preset...','ButtonPushedFcn',@(~,~)uiSavePreset(fig));
uibutton(tb,'Text','Load preset...','ButtonPushedFcn',@(~,~)uiLoadPreset(fig));
c.progLabel = uilabel(tb,'Text','Ready.','FontAngle','italic');
c.previewBtn = uibutton(tb,'Text','Preview order','ButtonPushedFcn',@(~,~)dryRun(fig), ...
    'Tooltip','list, in execution order, every checked cell that would run (no wrapper is called)');
c.runBtn = uibutton(tb,'Text','Run','ButtonPushedFcn',@(~,~)runChecked(fig), ...
    'BackgroundColor',[0.82 0.92 0.82],'Tooltip','run every checked cell in dependency order');
c.stopBtn = uibutton(tb,'Text','Stop','Enable','off','ButtonPushedFcn',@(~,~)cancelRun(fig), ...
    'BackgroundColor',[1 0.86 0.86],'Tooltip','stop after the current cell finishes');

c.matrixPanel = uipanel(left,'BorderType','none');   % the matrix uigridlayout owns the scroll
c.log = uitextarea(left,'Value',{'Ready.'},'Editable','off');

% -- RIGHT: settings panel for the selected step + a reports panel --
right = uigridlayout(gl,[4 1],'RowHeight',{'fit','fit','1x','1x'},'RowSpacing',6);
srow = uigridlayout(right,[1 2],'ColumnWidth',{'fit','1x'},'Padding',[0 0 0 0]);
uilabel(srow,'Text','Step','FontWeight','bold');
c.stepDrop = uidropdown(srow,'Items',{app.reg.label},'ItemsData',{app.reg.id}, ...
    'ValueChangedFcn',@(s,~)selectStep(fig,s.Value),'Tooltip','the step whose settings are shown below');
c.stepInfo = uilabel(right,'Text','','WordWrap','on','FontAngle','italic');
c.paramPanel = uipanel(right,'Scrollable','on','Title','Settings','FontWeight','bold');
c.reportPanel = uipanel(right,'Scrollable','on','Title','Reports (click a done cell)','FontWeight','bold');

app.c.process = c; setApp(fig,app);
end

%% ===================== EXPORT tab ================================== %%
function buildExportTab(fig)
%buildExportTab  A sheet/format selection UI in front of exportToExcel: pick which
%   of the loaded recordings to export, which sheets to write, and the format.
app = getApp(fig); t = app.tabs.export;
gl = uigridlayout(t,[1 2],'ColumnWidth',{'1.4x','1x'},'Padding',[8 8 8 8],'ColumnSpacing',10);

% -- LEFT: exportable recordings + run --
left = uigridlayout(gl,[4 1],'RowHeight',{'fit','1x','fit','fit'},'RowSpacing',6);
uilabel(left,'Text','Export processed results to Excel (wraps exportToExcel).','FontWeight','bold');
c.fileList = uilistbox(left,'Items',{},'Multiselect','on', ...
    'Tooltip',['The *_BFI_d.mat recordings discovered for the loaded files.  ' ...
               'Selected ones are exported (one workbook each).'], ...
    'ValueChangedFcn',@(~,~)refreshExportSheets(fig));
rr = uigridlayout(left,[1 3],'ColumnWidth',{'fit','fit','1x'},'Padding',[0 0 0 0]);
uibutton(rr,'Text','Refresh files','ButtonPushedFcn',@(~,~)refreshExportFiles(fig), ...
    'Tooltip','re-scan the loaded recordings for exportable *_BFI results');
uibutton(rr,'Text','Export selected','BackgroundColor',[0.82 0.92 0.82], ...
    'ButtonPushedFcn',@(~,~)uiRunExport(fig),'Tooltip','write one workbook per selected recording');
c.exportStatus = uilabel(rr,'Text','No files.','FontAngle','italic');

% -- RIGHT: format + which sheets --
right = uigridlayout(gl,[3 1],'RowHeight',{'fit','1x','fit'},'RowSpacing',6);
fr = uigridlayout(right,[1 2],'ColumnWidth',{'fit','1x'},'Padding',[0 0 0 0]);
uilabel(fr,'Text','Format');
c.fmt = uidropdown(fr,'Items',{'xlsx','xls'},'Value','xlsx', ...
    'Tooltip','output workbook format (both hold multiple sheets)');
sp = uipanel(right,'Title','Sheets to write','FontWeight','bold','Scrollable','on');
names = exportSheetNames();
sg = uigridlayout(sp,[numel(names) 1],'RowHeight',repmat({'fit'},1,numel(names)), ...
    'RowSpacing',2,'Padding',[6 6 6 6],'Scrollable','on');    % the grid owns the scroll, not the panel
for i = 1:numel(names)
    c.sheetChk.(names{i}) = uicheckbox(sg,'Text',names{i},'Value',true, ...
        'Tooltip','tick to include this sheet (greyed = not produced by any selected file)');
end
brow = uigridlayout(right,[1 2],'ColumnWidth',{'1x','1x'},'Padding',[0 0 0 0]);
uibutton(brow,'Text','All','ButtonPushedFcn',@(~,~)setAllSheets(fig,true));
uibutton(brow,'Text','None','ButtonPushedFcn',@(~,~)setAllSheets(fig,false));

app.c.export = c; setApp(fig,app);
refreshExportFiles(fig);
end

function names = exportSheetNames()
%exportSheetNames  The canonical sheets exportToExcel can produce, in write order.
names = {'sMetrics','sData','sMetricsROI','sDataROI','dvsMetrics','dvsData', ...
         'dvsDiameter','pulsatility','dvsPulsatility'};
end

function refreshExportFiles(fig)
%refreshExportFiles  Rebuild the exportable-recording list from the loaded rows.
app = getApp(fig); c = app.c.export;
srcs = resolveExportSources(app);
app.exportSrcs = srcs; setApp(fig,app);
items = cell(1,numel(srcs));
for i = 1:numel(srcs), items{i} = sprintf('%s   [%s]', srcs(i).label, srcs(i).group); end
c.fileList.Items = items;
if isempty(srcs)
    c.fileList.ItemsData = {};
    c.fileList.Value = {};                            % empty Items -> Value must be {}
    c.exportStatus.Text = 'No exportable *_BFI results for the loaded files.';
else
    c.fileList.ItemsData = 1:numel(srcs);
    c.fileList.Value = 1:numel(srcs);                 % all selected by default
    c.exportStatus.Text = sprintf('%d exportable recording(s).', numel(srcs));
end
refreshExportSheets(fig);
end

function refreshExportSheets(fig)
%refreshExportSheets  Grey / untick sheets no selected file can produce.
app = getApp(fig); c = app.c.export;
sel = selectedExportSrcs(app,c);
avail = availableExportSheets({sel.rpath});
names = exportSheetNames();
for i = 1:numel(names)
    nm = names{i};
    if ~isfield(c.sheetChk,nm) || ~isgraphics(c.sheetChk.(nm)), continue; end
    chk = c.sheetChk.(nm);
    if any(strcmp(nm,avail))
        if strcmp(chk.Enable,'off'), chk.Value = true; end     % re-enabled -> re-tick
        chk.Enable = 'on';
    else
        chk.Value = false; chk.Enable = 'off';                 % can't write it
    end
end
end

function srcs = resolveExportSources(app)
%resolveExportSources  Exportable *_BFI_d.mat sources for the loaded recordings.
%   For each loaded row, glob the recording's base name for *_BFI_d.mat, keep the
%   identity matches whose _r.mat sibling (exportToExcel's real input) exists.
srcs = struct('path',{},'rpath',{},'group',{},'label',{});
seen = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(app.rows)
    r = app.rows(i); m = r.model;
    if isempty(m.folder) || ~isfolder(m.folder), continue; end
    base = [m.roiPrefix m.stem];
    d = dir(fullfile(m.folder,[base '*_BFI_d.mat']));
    for k = 1:numel(d)
        p = fullfile(d(k).folder, d(k).name);
        if isKey(seen,p), continue; end
        cm = wbFileModel(p);
        if ~strcmp(cm.identity, m.identity), continue; end
        rp = strrep(p,'_d.mat','_r.mat');
        if ~isfile(rp), continue; end                          % export needs the RESULTS file
        seen(p) = true;
        [~,nm,ex] = fileparts(p);
        srcs(end+1) = struct('path',p,'rpath',rp,'group',r.group,'label',[nm ex]); %#ok<AGROW>
    end
end
end

function sel = selectedExportSrcs(app,c)
%selectedExportSrcs  The exportable sources at the currently-selected list rows.
srcs = app.exportSrcs;
if isempty(srcs), sel = srcs; return; end
idx = c.fileList.Value;
if isempty(idx), idx = 1:numel(srcs); end
sel = srcs(idx);
end

function names = availableExportSheets(rpaths)
%availableExportSheets  Union (canonical order) of sheets producible from the given
%   RESULTS files - mirrors exportToExcel's own presence tests, so the tick list
%   matches what a full export would actually write.
all = exportSheetNames();
found = {};
for i = 1:numel(rpaths), found = union(found, availOne(rpaths{i}), 'stable'); end
names = all(ismember(all,found));
end

function a = availOne(rpath)
%availOne  Sheets a single RESULTS file can produce.  A very large file (a per-pixel
%   vasomotion result) is not fully loaded just to grey a checkbox - assume the base
%   sheets and let exportToExcel skip any that turn out absent.
a = {};
if isempty(rpath) || ~isfile(rpath), return; end
d = dir(rpath);
if ~isempty(d) && d(1).bytes > 1.5e9, a = {'sMetrics','sData'}; return; end
try
    S = load(rpath,'results');
catch
    return
end
if ~isfield(S,'results'), return; end
R = S.results;
if isfield(R,'sMetrics') && istable(R.sMetrics)
    a{end+1} = 'sMetrics';
    if isfield(R,'sData') && isfield(R,'time'), a{end+1} = 'sData'; end
    if ismember('label', R.sMetrics.Properties.VariableNames)
        a{end+1} = 'sMetricsROI'; a{end+1} = 'sDataROI';
    end
end
if isfield(R,'dvsMetrics') && istable(R.dvsMetrics)
    a{end+1} = 'dvsMetrics';
    if isfield(R,'dvsData'),     a{end+1} = 'dvsData';     end
    if isfield(R,'dvsDiameter'), a{end+1} = 'dvsDiameter'; end
end
if isfield(R,'pulsatility') && isstruct(R.pulsatility)
    P = R.pulsatility;
    if isfield(P,'sData') && isstruct(P.sData) && isfield(P.sData,'scalars')
        a{end+1} = 'pulsatility';
    end
    if isfield(R,'dvsMetrics') && isfield(P,'dvsData') && isstruct(P.dvsData) ...
            && isfield(P.dvsData,'scalars')
        a{end+1} = 'dvsPulsatility';
    end
end
end

function s = tickedSheets(c)
%tickedSheets  The sheet names whose checkbox is on (and enabled).
names = exportSheetNames(); s = {};
for i = 1:numel(names)
    nm = names{i};
    if isfield(c.sheetChk,nm) && isgraphics(c.sheetChk.(nm)) ...
            && strcmp(c.sheetChk.(nm).Enable,'on') && c.sheetChk.(nm).Value
        s{end+1} = nm; %#ok<AGROW>
    end
end
end

function setAllSheets(fig,tf)
%setAllSheets  Tick / untick every currently-enabled sheet checkbox.
app = getApp(fig); c = app.c.export; names = exportSheetNames();
for i = 1:numel(names)
    nm = names{i};
    if isfield(c.sheetChk,nm) && isgraphics(c.sheetChk.(nm)) && strcmp(c.sheetChk.(nm).Enable,'on')
        c.sheetChk.(nm).Value = tf;
    end
end
end

function uiRunExport(fig)
%uiRunExport  Gather the selection from the tab and run exportToExcel over it.
app = getApp(fig); c = app.c.export;
sel = selectedExportSrcs(app,c);
if isempty(sel),                setExportStatus(fig,'Nothing selected to export.'); return; end
sheets = tickedSheets(c);
if isempty(sheets),             setExportStatus(fig,'No sheets ticked.'); return; end
opts = struct('sheets',{sheets}, 'format',['.' c.fmt.Value]);
runExport(fig, {sel.path}, opts);
end

function paths = runExport(fig, files, opts)
%runExport  Call exportToExcel per file (per-file logging + continue-on-error),
%   mirror the outcome into the Process-tab log, and refresh the matrix.
if nargin<3 || isempty(opts), opts = struct(); end
paths = {}; nok = 0; nfail = 0;
wbLog(fig, sprintf('=== EXPORT: %d file(s) ===', numel(files)));
for i = 1:numel(files)
    f = files{i};
    try
        exportToExcel({f}, opts);
        outp = strrep(f,'_d.mat', optExt(opts));
        paths{end+1} = outp; %#ok<AGROW>
        nok = nok + 1;
        wbLog(fig, ['  wrote ' outp]);
    catch ME
        nfail = nfail + 1;
        wbLog(fig, sprintf('  ERROR %s - %s', shortName(f), ME.message));
    end
end
wbLog(fig, sprintf('=== EXPORT complete: %d ok, %d failed ===', nok, nfail));
setExportStatus(fig, sprintf('Exported %d file(s), %d failed. See the Process-tab log.', nok, nfail));
recomputeBase(fig); refreshCells(fig,{});          % a new .xlsx flips the export cell to done
end

function e = optExt(opts)
%optExt  '.xlsx' (default) or '.xls' from an export opts struct.
e = '.xlsx';
if isstruct(opts) && isfield(opts,'format') && ~isempty(opts.format)
    e = lower(char(string(opts.format))); if e(1)~='.', e = ['.' e]; end
end
end

function p = exportSourcePaths(fig)
%exportSourcePaths  The exportable *_BFI_d.mat paths currently listed (for tests).
app = getApp(fig);
if isfield(app,'exportSrcs') && ~isempty(app.exportSrcs), p = {app.exportSrcs.path}; else, p = {}; end
end

function setExportStatus(fig,msg)
app = getApp(fig);
if isfield(app.c,'export') && isfield(app.c.export,'exportStatus') && isgraphics(app.c.export.exportStatus)
    app.c.export.exportStatus.Text = msg;
end
drawnow limitrate;
end

%% ===================== EXPLORE tab ================================= %%
function buildExploreTab(fig)
%buildExploreTab  Host guiExplore inside the Explore tab (a light toolbar + the
%   embedded app).  If hosting fails on some MATLAB build, fall back to a button
%   that launches guiExplore in its own window (still seeded from the workbench).
app = getApp(fig); t = app.tabs.explore;
gl = uigridlayout(t,[2 1],'RowHeight',{'fit','1x'},'Padding',[6 6 6 6],'RowSpacing',6);
tb = uigridlayout(gl,[1 2],'ColumnWidth',{'fit','1x'},'Padding',[0 0 0 0]);
uibutton(tb,'Text','Load workbench files & groups','ButtonPushedFcn',@(~,~)seedExplore(fig,true), ...
    'Tooltip','seed the explorer with the workbench''s loaded recordings (_r.mat) and groups');
c.exploreStatus = uilabel(tb,'Text','(switch to this tab to seed from the loaded files)','FontAngle','italic');
host = uipanel(gl,'BorderType','none');
app.c.explore = c; app.exploreSeedKey = '';
try
    guiExplore('Parent',host);                     % embed the whole explorer in the tab
    app.exploreAPI    = getappdata(fig,'exploreAPI');
    app.exploreHosted = true;
catch ME
    app.exploreAPI    = [];
    app.exploreHosted = false;
    delete(host.Children);
    fb = uigridlayout(host,[2 1],'RowHeight',{'fit','fit'},'Padding',[16 16 16 16],'RowSpacing',8);
    uilabel(fb,'Text',['Could not host guiExplore in the tab (' ME.message ').'],'WordWrap','on');
    uibutton(fb,'Text','Open guiExplore in its own window','ButtonPushedFcn',@(~,~)launchExploreChild(fig));
end
setApp(fig,app);
end

function launchExploreChild(fig)
%launchExploreChild  Fallback: open guiExplore standalone and seed it.
h = guiExplore;
app = getApp(fig);
app.exploreAPI = getappdata(h,'exploreAPI');
app.exploreHosted = false; app.exploreSeedKey = '';
setApp(fig,app);
seedExplore(fig,true);
end

function seedExplore(fig, force)
%seedExplore  Push the workbench's exportable RESULTS files + groups into the
%   hosted guiExplore (only when the file set actually changed, unless forced).
if nargin<2, force = false; end
app = getApp(fig);
if ~isfield(app,'exploreAPI') || isempty(app.exploreAPI)
    setExploreStatus(fig,'Explore is not hosted - use the button to open it in its own window.');
    return
end
srcs = resolveExportSources(app);
rpaths = {srcs.rpath};
key = strjoin(sort(rpaths),'|');
if ~force && strcmp(key, app.exploreSeedKey), return; end       % nothing changed
app.exploreSeedKey = key; setApp(fig,app);
if isempty(rpaths)
    setExploreStatus(fig,'No processed *_r.mat results among the loaded files yet.');
    return
end
app.exploreAPI.loadPaths(rpaths, '_r\.mat$', '', '');           % seed the file list
a = app.exploreAPI.getApp(); order = {a.files.path};            % explore's file order
groups = unique({srcs.group},'stable');
for gi = 1:numel(groups)
    mem = {srcs(strcmp({srcs.group},groups{gi})).rpath};
    idx = find(ismember(order, mem));
    if ~isempty(idx), app.exploreAPI.createGroup(groups{gi}, idx); end
end
setExploreStatus(fig, sprintf('Seeded %d file(s) in %d group(s) from the workbench.', ...
    numel(rpaths), numel(groups)));
end

function setExploreStatus(fig,msg)
app = getApp(fig);
if isfield(app.c,'explore') && isfield(app.c.explore,'exploreStatus') && isgraphics(app.c.explore.exploreStatus)
    app.c.explore.exploreStatus.Text = msg;
end
end

%% ===================== tab-change seeding ========================== %%
function onTabSelected(fig)
%onTabSelected  Lazily refresh the Export list / seed Explore when its tab shows.
app = getApp(fig);
if ~isfield(app,'tg') || ~isgraphics(app.tg), return; end
tab = app.tg.SelectedTab;
if     isequal(tab, app.tabs.export),  refreshExportFiles(fig);
elseif isequal(tab, app.tabs.explore), seedExplore(fig,false);
end
end

%% ===================== loaders ===================================== %%
function uiLoadStructured(fig)
c = getApp(fig).c.files;
root = strtrim(c.root.Value);
if isempty(root) || ~isfolder(root), alert(fig,'Set a valid root folder first.'); return; end
disc = wbDiscoverFiles('structured', root, strtrim(c.glob.Value), ...
    strtrim(c.animal.Value), strtrim(c.ref.Value));
applyDiscovery(fig,disc);
end
function uiLoadFolder(fig)
c = getApp(fig).c.files;
root = uigetdir(defaultDir(c.root.Value),'Pick a folder to scan recursively');
if isequal(root,0), return; end
c.root.Value = root;
disc = wbDiscoverFiles('folder', root, strtrim(c.glob.Value), strtrim(c.animal.Value));
applyDiscovery(fig,disc);
end
function uiLoadManual(fig)
c = getApp(fig).c.files;
[f,p] = uigetfile({'*.mat;*.rls;*.cxd;*.avi','Recordings & products';'*.*','All files'}, ...
    'Pick files (multiselect)','MultiSelect','on',defaultDir(c.root.Value));
if isequal(f,0), return; end
if ischar(f), f = {f}; end
paths = cellfun(@(x)fullfile(p,x),f,'UniformOutput',false);
disc = wbDiscoverFiles('manual', paths, strtrim(c.animal.Value));
applyDiscovery(fig,disc);
end
function uiClear(fig)
applyDiscovery(fig,emptyDisc());
end
function apiLoad(fig,mode,varargin)
%apiLoad  Programmatic loader: run wbDiscoverFiles(mode,...) and adopt the result.
disc = wbDiscoverFiles(mode, varargin{:});
applyDiscovery(fig,disc);
end
function d = defaultDir(v)
if ~isempty(v) && isfolder(v), d = v; else, d = pwd; end
end

function applyDiscovery(fig,disc)
%applyDiscovery  Adopt a discovery result: flatten to rows, set modality, render.
app = getApp(fig);
app.disc = disc;
[app.rows, app.groupNames, app.modelArr] = flattenDisc(disc);
if ~isempty(app.modelArr)
    mods = {app.modelArr.modality};
    app.modality = modeStr(mods);
else
    app.modality = 'LSCI';
end
app.reg = wbStepRegistry(app.modality);
% a fresh load starts with no session overlay
app.checked = containers.Map('KeyType','char','ValueType','any');
app.stale   = containers.Map('KeyType','char','ValueType','any');
setApp(fig,app);
recomputeBase(fig);
refreshFileTable(fig);
refreshStepSelector(fig);
renderMatrix(fig); selectStep(fig,getApp(fig).selStep);
refreshStatus(fig);
end

function refreshStepSelector(fig)
%refreshStepSelector  Keep the settings-panel step dropdown valid for the columns.
app = getApp(fig); c = app.c.process;
if isempty(app.reg)
    c.stepDrop.Items = {'(no steps)'}; c.stepDrop.ItemsData = {''};
    app.selStep = ''; setApp(fig,app); return
end
c.stepDrop.Items = {app.reg.label}; c.stepDrop.ItemsData = {app.reg.id};
if ~any(strcmp(app.selStep,{app.reg.id})), app.selStep = app.reg(1).id; end
c.stepDrop.Value = app.selStep; setApp(fig,app);
end

function [rows, groupNames, modelArr] = flattenDisc(disc)
%flattenDisc  Discovery grid -> flat rows in group-major, reference-first order.
%   ONE ROW PER RECORDING IDENTITY.  A recording can own several branch products
%   with the SAME identity (e.g. a _t_K and a _c_K from one .rls), so a glob that
%   matches more than one branch (*_K_d.mat) would otherwise emit duplicate rows.
%   The first model seen anchors the row; every step still resolves its own branch
%   file on disk (wbExecutor), and the disk state unions all of the recording's
%   _s files (wbStateEngine), so the anchor branch does not change what runs.  In
%   the normal single-entry-glob workflow no identity ever repeats, so this is a
%   no-op there.
rows = emptyRows(); groupNames = {}; modelArr = wbFileModel('x.rls'); modelArr(1) = [];
if isempty(disc.models), return; end
seen = containers.Map('KeyType','char','ValueType','logical');
[nr,nc] = size(disc.models);
for r = 1:nr
    if ~isempty(disc.groups) && numel(disc.groups)>=r, gname = disc.groups(r).name; else, gname = sprintf('group%d',r); end
    groupNames{end+1} = gname; %#ok<AGROW>
    for cc = 1:nc
        m = disc.models{r,cc};
        if isempty(m), continue; end
        if isKey(seen, m.identity), continue; end     % dedup: one row per recording
        seen(m.identity) = true;
        row = struct('model',m,'identity',m.identity,'group',gname,'groupIdx',r, ...
            'rowInGroup',cc,'isRef',logical(disc.referenceMode && cc==1), ...
            'label',fileLabel(m));
        rows(end+1) = row; %#ok<AGROW>
        modelArr(end+1) = m; %#ok<AGROW>
    end
end
end
function s = fileLabel(m), s = m.name; if isempty(s), s = [m.roiPrefix m.stem]; end, end
function s = modeStr(c)
if isempty(c), s = 'LSCI'; return; end
u = unique(c); n = cellfun(@(x)sum(strcmp(x,c)),u); [~,i] = max(n); s = u{i};
end

function refreshFileTable(fig)
app = getApp(fig); c = app.c.files;
n = numel(app.rows);
D = cell(n,5);
for i = 1:n
    r = app.rows(i); m = r.model;
    D(i,:) = {r.group, fileLabel(m), ternary(r.isRef,'ref',''), m.modality, stageLabel(m)};
end
c.fileTbl.Data = D;
end
function s = stageLabel(m)
if m.isRaw, s = 'raw'; elseif isempty(m.flags), s = m.product; else, s = [strjoin(m.flags,'_') '_' m.product]; end
end
function refreshStatus(fig)
app = getApp(fig);
if isempty(app.rows)
    txt = 'No files loaded.';
else
    txt = sprintf('%d files in %d group(s); modality %s.', ...
        numel(app.rows), numel(app.groupNames), app.modality);
end
app.c.files.status.Text = txt;
end

%% ===================== state derivation ============================ %%
function recomputeBase(fig)
%recomputeBase  Disk baseline (+ current-settings fingerprint) for every row.
app = getApp(fig);
app.base = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(app.rows)
    m = app.rows(i).model;
    cs = curSettingsFor(app,m);
    st = wbStateEngine(m, app.reg, cs);
    bs = struct();
    for k = 1:numel(st), bs.(st(k).id) = st(k).state; end
    app.base(m.identity) = bs;
end
setApp(fig,app);
end
function cs = curSettingsFor(app,model)
cs = struct();
for k = 1:numel(app.reg)
    cs.(app.reg(k).id) = wbSettingsModel('resolve', app.sm, app.reg(k), model);
end
end

function s = resolveCellState(app,identity,stepId)
%resolveCellState  Disk baseline with the live session overlay (checked/stale).
key = cellKey(identity,stepId);
% transient run overlay (running/done/error) is authoritative during & after a run
if isfield(app,'runState') && isKey(app.runState,key)
    s = app.runState(key); return
end
if ~isKey(app.base,identity), s = 'unavailable'; return; end
bs = app.base(identity);
if ~isfield(bs,stepId), s = 'unavailable'; return; end
s = bs.(stepId);
if strcmp(s,'done') && isKey(app.stale,key)
    s = 'stale';                                   % an in-session edit pushed it stale
end
if strcmp(s,'unavailable')
    s = projectCheckable(app,identity,stepId);     % promotable if prereqs are queued
end
if any(strcmp(s,{'ready','stale'})) && isKey(app.checked,key) && app.checked(key)
    s = 'checked';                                 % queued for the next run
end
end

function s = projectCheckable(app,identity,stepId)
%projectCheckable  Look-ahead gating: an 'unavailable' cell becomes 'ready' when
%   it is blocked ONLY by prerequisites that are themselves already done on disk or
%   QUEUED (checked) for this recording - so a whole chain (contrast -> setRegions
%   -> segmentation -> ...) can be queued in one pass, provided it is checked in
%   dependency order (the executor then runs it in that order).  Entry steps (no
%   requires) and modality-inapplicable steps are never promoted: an entry step
%   needs its real input, and an inapplicable step's prereqs are inapplicable too,
%   so this returns 'unavailable' for them.
s = 'unavailable';
step = stepById(app.reg,stepId);
if isempty(step) || isempty(step.requires), return; end
for i = 1:numel(step.requires)
    reqState = resolveCellState(app,identity,step.requires{i});   % recurse over the DAG
    if ~any(strcmp(reqState,{'done','stale','checked'})), return; end
end
s = 'ready';
end

%% ===================== the matrix =================================== %%
function renderMatrix(fig)
app = getApp(fig);
c = app.c.process;
delete(c.matrixPanel.Children);
app.cellComp  = containers.Map('KeyType','char','ValueType','any');
app.gridRowOf = containers.Map('KeyType','char','ValueType','double');
app.rendered  = containers.Map('KeyType','char','ValueType','any');
setApp(fig,app);

if isempty(app.rows)
    gl = uigridlayout(c.matrixPanel,[1 1],'Padding',[20 20 20 20]);
    uilabel(gl,'Text','Load recordings on the Files tab to populate the matrix.', ...
        'HorizontalAlignment','center','FontAngle','italic');
    return
end
if isempty(app.reg)
    gl = uigridlayout(c.matrixPanel,[1 1],'Padding',[20 20 20 20]);
    uilabel(gl,'Text',sprintf('No pipeline steps are defined for modality "%s" yet.',app.modality), ...
        'HorizontalAlignment','center','FontAngle','italic');
    return
end
if numel(app.rows) > app.maxWidgets
    renderMatrixTable(fig);
    return
end
renderMatrixWidgets(fig);
end

function renderMatrixWidgets(fig)
app = getApp(fig); c = app.c.process; reg = app.reg;
nSteps = numel(reg);
nCols  = 1 + nSteps;
nGroups = numel(app.groupNames);
nBody  = numel(app.rows);
nRows  = 1 + nGroups + nBody;                      % header + group headers + files

colW = [{190}, repmat({96},1,nSteps)];
rowH = [{58}, num2cell(repmat(26,1,nRows-1))];
grid = uigridlayout(c.matrixPanel,[nRows nCols],'ColumnWidth',colW,'RowHeight',rowH, ...
    'RowSpacing',2,'ColumnSpacing',2,'Padding',[2 2 2 2],'Scrollable','on');

% ---- header row ----
hc = uilabel(grid,'Text','Recording  \  Step','FontWeight','bold'); hc.Layout.Row = 1; hc.Layout.Column = 1;
for s = 1:nSteps
    hp = uigridlayout(grid,[2 1],'RowHeight',{'1x','fit'},'Padding',[1 1 1 1],'RowSpacing',1);
    hp.Layout.Row = 1; hp.Layout.Column = s+1;
    cb = uicheckbox(hp,'Text',reg(s).label,'WordWrap','on','FontSize',10, ...
        'Tooltip',headerTip(reg(s)),'ValueChangedFcn',@(o,~)checkColumn(fig,reg(s).id,o.Value));
    gr = uigridlayout(hp,[1 1],'Padding',[0 0 0 0]);
    uibutton(gr,'Text','settings','FontSize',9,'ButtonPushedFcn',@(~,~)selectStep(fig,reg(s).id));
    app.c.process.colCheck.(reg(s).id) = cb;
end
setApp(fig,app);

% ---- body: group header + one row per file ----
gridRow = 1;
lastGroup = -1;
app = getApp(fig);
for i = 1:numel(app.rows)
    r = app.rows(i);
    if r.groupIdx ~= lastGroup
        gridRow = gridRow + 1;
        gh = uigridlayout(grid,[1 3],'ColumnWidth',{'1x','fit','fit'},'Padding',[0 0 0 0]);
        gh.Layout.Row = gridRow; gh.Layout.Column = [1 nCols];
        uilabel(gh,'Text',['  ' r.group],'FontWeight','bold','BackgroundColor',[0.92 0.92 0.96]);
        uibutton(gh,'Text','check group','FontSize',9,'ButtonPushedFcn',@(~,~)checkGroup(fig,r.group,true));
        uibutton(gh,'Text','clear','FontSize',9,'ButtonPushedFcn',@(~,~)checkGroup(fig,r.group,false));
        lastGroup = r.groupIdx;
    end
    gridRow = gridRow + 1;
    app.gridRowOf(r.identity) = gridRow;
    lc = uigridlayout(grid,[1 2],'ColumnWidth',{'1x','fit'},'Padding',[0 0 0 0]);
    lc.Layout.Row = gridRow; lc.Layout.Column = 1;
    uilabel(lc,'Text',[r.label ternary(r.isRef,'  (ref)','')],'Tooltip',r.model.path,'FontSize',11);
    uibutton(lc,'Text','check','FontSize',9,'ButtonPushedFcn',@(~,~)checkRow(fig,r.identity,true));
    for s = 1:nSteps
        makeCell(fig,grid,gridRow,s+1,r.identity,reg(s));
    end
end
setApp(fig,app);
syncHeaderChecks(fig);
end

function makeCell(fig,grid,gridRow,gridCol,identity,step)
%makeCell  Create (or replace) the one component for a (file,step) cell.
app = getApp(fig);
key = cellKey(identity,step.id);
state = resolveCellState(app,identity,step.id);
if isKey(app.cellComp,key)
    old = app.cellComp(key);
    if isgraphics(old), delete(old); end
end
comp = cellComponent(fig,grid,identity,step,state);
comp.Layout.Row = gridRow; comp.Layout.Column = gridCol;
app.cellComp(key) = comp;
app.rendered(key) = state;
setApp(fig,app);
end

function comp = cellComponent(fig,parent,identity,step,state)
%cellComponent  A compact stateful widget per cell (checkbox / glyph label).
switch state
    case {'ready','stale','checked'}
        isChk = strcmp(state,'checked');
        txt = ''; fc = [0 0 0]; tip = 'ready - tick to queue';
        if strcmp(state,'stale'), txt = '!'; fc = [0.75 0.45 0]; tip = 'stale - tick to re-run with the current settings'; end
        if isChk, tip = 'queued for the next run'; end
        comp = uicheckbox(parent,'Text',txt,'FontColor',fc,'Value',isChk,'Tooltip',tip, ...
            'ValueChangedFcn',@(o,~)onCellCheck(fig,identity,step.id,o.Value));
    case 'done'
        comp = uibutton(parent,'Text','done','FontColor',[0 0.5 0], ...
            'Tooltip','done - click to open its report image(s)', ...
            'ButtonPushedFcn',@(~,~)showReports(fig,identity,step.id));
    case 'running'
        comp = uilabel(parent,'Text','...','FontColor',[0.85 0.5 0],'HorizontalAlignment','center','Tooltip','running');
    case 'error'
        msg = 'error'; app = getApp(fig); key = cellKey(identity,step.id);
        if isKey(app.cellMsg,key) && ~isempty(app.cellMsg(key)), msg = ['error: ' app.cellMsg(key)]; end
        comp = uilabel(parent,'Text','X','FontColor',[0.8 0 0],'HorizontalAlignment','center','Tooltip',msg);
    otherwise   % unavailable - a disabled (greyed) checkbox reads as "not yet checkable"
        reason = unavailableReason(getApp(fig),identity,step.id);
        comp = uicheckbox(parent,'Text','','Value',false,'Enable','off','Tooltip',reason);
end
end

function s = unavailableReason(app,identity,stepId)
%unavailableReason  A short tooltip for a greyed cell (modality / prerequisites).
s = 'not available yet';
step = stepById(app.reg,stepId);
if isKey(app.base,identity) && ~any(strcmp(app.modality,step.modalities))
    s = ['not for ' app.modality];
elseif ~isempty(step.requires)
    s = ['needs: ' strjoin(step.requires,', ')];
end
end

function renderMatrixTable(fig)
%renderMatrixTable  Escape hatch for very large sessions: a colored summary table.
app = getApp(fig); c = app.c.process; reg = app.reg;
gl = uigridlayout(c.matrixPanel,[2 1],'RowHeight',{'fit','1x'},'Padding',[4 4 4 4]);
uilabel(gl,'Text',sprintf(['%d files exceed the widget limit (%d): showing a read-only ' ...
    'state table. Use the toolbar / column / group checks to queue.'],numel(app.rows),app.maxWidgets), ...
    'WordWrap','on','FontAngle','italic');
n = numel(app.rows);
D = cell(n, numel(reg)+1);
for i = 1:n
    r = app.rows(i);
    D{i,1} = [r.group ' / ' r.label];
    for s = 1:numel(reg)
        D{i,s+1} = stateGlyph(resolveCellState(app,r.identity,reg(s).id));
    end
end
tbl = uitable(gl,'Data',D,'ColumnName',[{'recording'},{reg.id}],'RowName',{});
% color the done/stale cells
styGreen = uistyle('FontColor',[0 0.5 0]);
styAmber = uistyle('FontColor',[0.75 0.45 0]);
for i = 1:n
    for s = 1:numel(reg)
        g = D{i,s+1};
        if strcmp(g,'done'), addStyle(tbl,styGreen,'cell',[i s+1]);
        elseif any(strcmp(g,{'stale','[x]'})), addStyle(tbl,styAmber,'cell',[i s+1]); end
    end
end
end
function g = stateGlyph(state)
switch state
    case 'checked', g = '[x]';
    case 'ready',   g = '[ ]';
    case 'done',    g = 'done';
    case 'stale',   g = 'stale';
    case 'running', g = '...';
    case 'error',   g = 'X';
    otherwise,      g = '-';
end
end

function refreshCells(fig,identities)
%refreshCells  Redraw only the cells whose resolved state changed (widget mode).
app = getApp(fig);
if isempty(app.rows) || numel(app.rows) > app.maxWidgets
    renderMatrix(fig); return              % table mode: cheapest to rebuild
end
grid = firstMatrixGrid(app.c.process.matrixPanel);
if isempty(grid), renderMatrix(fig); return; end
if nargin<2 || isempty(identities), identities = {app.rows.identity}; end
reg = app.reg;
for i = 1:numel(identities)
    id = identities{i};
    if ~isKey(app.gridRowOf,id), continue; end
    gridRow = app.gridRowOf(id);
    for s = 1:numel(reg)
        key = cellKey(id,reg(s).id);
        newState = resolveCellState(app,id,reg(s).id);
        old = ''; if isKey(app.rendered,key), old = app.rendered(key); end
        if ~strcmp(newState,old)
            makeCell(fig,grid,gridRow,s+1,id,reg(s));
            app = getApp(fig);             % makeCell mutated appdata
        end
    end
end
syncHeaderChecks(fig);
end
function grid = firstMatrixGrid(panel)
grid = [];
ch = panel.Children;
for i = 1:numel(ch)
    if isa(ch(i),'matlab.ui.container.GridLayout'), grid = ch(i); return; end
end
end

function syncHeaderChecks(fig)
%syncHeaderChecks  Tick a column header when every eligible cell in it is checked.
app = getApp(fig);
if ~isfield(app.c.process,'colCheck'), return; end
for s = 1:numel(app.reg)
    id = app.reg(s).id;
    if ~isfield(app.c.process.colCheck,id), continue; end
    cb = app.c.process.colCheck.(id);
    if ~isgraphics(cb), continue; end
    [elig,chk] = columnTally(app,id);
    cb.Value = elig>0 && chk==elig;
end
end
function [elig,chk] = columnTally(app,stepId)
elig = 0; chk = 0;
for i = 1:numel(app.rows)
    st = resolveCellState(app,app.rows(i).identity,stepId);
    if any(strcmp(st,{'ready','stale','checked'})), elig = elig+1; end
    if strcmp(st,'checked'), chk = chk+1; end
end
end

%% ===================== check / bulk-select ========================= %%
function onCellCheck(fig,identity,stepId,tf)
setChecked(fig,identity,stepId,tf);
refreshCells(fig,{identity});
end
function apiCheck(fig,identity,stepId,tf)
setChecked(fig,identity,stepId,tf);
refreshCells(fig,{identity});
end
function setChecked(fig,identity,stepId,tf)
app = getApp(fig);
key = cellKey(identity,stepId);
st = resolveCellState(app,identity,stepId);
if tf
    if any(strcmp(st,{'ready','stale','checked'})), app.checked(key) = true; end
else
    if isKey(app.checked,key), remove(app.checked,key); end
end
setApp(fig,app);
end
function checkColumn(fig,stepId,tf)
app = getApp(fig);
for i = 1:numel(app.rows), setCheckedQuiet(app,app.rows(i).identity,stepId,tf); end
setApp(fig,app); refreshCells(fig,{});
end
function checkGroup(fig,groupName,tf)
app = getApp(fig);
for i = 1:numel(app.rows)
    if strcmp(app.rows(i).group,groupName)
        for s = 1:numel(app.reg), setCheckedQuiet(app,app.rows(i).identity,app.reg(s).id,tf); end
    end
end
setApp(fig,app); refreshCells(fig,{});
end
function checkRow(fig,identity,tf)
app = getApp(fig);
for s = 1:numel(app.reg), setCheckedQuiet(app,identity,app.reg(s).id,tf); end
setApp(fig,app); refreshCells(fig,{identity});
end
function checkModality(fig,modality,tf)
app = getApp(fig);
for i = 1:numel(app.rows)
    if strcmp(app.rows(i).model.modality,modality)
        for s = 1:numel(app.reg), setCheckedQuiet(app,app.rows(i).identity,app.reg(s).id,tf); end
    end
end
setApp(fig,app); refreshCells(fig,{});
end
function checkAll(fig,tf)
app = getApp(fig);
for i = 1:numel(app.rows)
    for s = 1:numel(app.reg), setCheckedQuiet(app,app.rows(i).identity,app.reg(s).id,tf); end
end
setApp(fig,app); refreshCells(fig,{});
end
function setCheckedQuiet(app,identity,stepId,tf)
%setCheckedQuiet  Mutate the checked map on an app struct without saving/rendering.
key = cellKey(identity,stepId);
st = resolveCellState(app,identity,stepId);
if tf
    if any(strcmp(st,{'ready','stale','checked'})), app.checked(key) = true; end
else
    if isKey(app.checked,key), remove(app.checked,key); end
end
end
function L = checkedList(fig)
app = getApp(fig);
k = keys(app.checked);
L = cell(numel(k),2);
for i = 1:numel(k)
    parts = strsplit(k{i},'||');
    L(i,:) = {parts{1}, parts{2}};
end
end

%% ===================== settings panel ============================== %%
function selectStep(fig,stepId)
app = getApp(fig);
if ~any(strcmp(stepId,{app.reg.id})), return; end
app.selStep = stepId; setApp(fig,app);
c = app.c.process;
if isprop(c.stepDrop,'Value'), c.stepDrop.Value = stepId; end
step = stepById(app.reg,stepId);
c.stepInfo.Text = stepInfoText(step);
buildSettingsPanel(fig,c.paramPanel,step);
end
function s = stepInfoText(step)
bits = {};
if ~isempty(step.gatingField), bits{end+1} = ['gates on settings.' step.gatingField];
else, bits{end+1} = 'done by output (no settings field)'; end
if strcmp(step.arity,'perGroup'), bits{end+1} = 'per-group (reference in column 1)'; end
if ~isequal(step.interactive,false), bits{end+1} = 'interactive'; end
if step.needsRaw, bits{end+1} = 'also needs the raw recording'; end
if ~isempty(step.requires), bits{end+1} = ['requires ' strjoin(step.requires,', ')]; end
s = strjoin(bits,'  |  ');
end

function buildSettingsPanel(fig,panel,step)
%buildSettingsPanel  Generalised buildParamEditor: render step.settingGroups.
delete(panel.Children);
app = getApp(fig);
groups = step.settingGroups;
if isempty(groups)
    gl = uigridlayout(panel,[1 1],'Padding',[8 8 8 8]);
    uilabel(gl,'Text','This step has no tunable parameters (interactive or automatic).', ...
        'WordWrap','on','FontAngle','italic');
    return
end
s = wbSettingsModel('resolve', app.sm, step);      % step-level resolution (no file)
stack = uigridlayout(panel,[size(groups,1) 1],'RowHeight',repmat({'fit'},1,size(groups,1)), ...
    'RowSpacing',8,'Padding',[6 6 6 6],'Scrollable','on');   % the grid owns the scroll, not the panel
for gi = 1:size(groups,1)
    p = uipanel(stack,'Title',groups{gi,1},'FontWeight','bold');
    flds = groups{gi,2};
    gg = uigridlayout(p,[numel(flds) 2],'ColumnWidth',{'fit','1x'}, ...
        'RowHeight',repmat({'fit'},1,numel(flds)),'RowSpacing',3);
    for fi = 1:numel(flds)
        f = flds{fi};
        lbl = uilabel(gg,'Text',f);
        if isfield(step.tips,f), lbl.Tooltip = step.tips.(f); end
        val = getfieldOr(s,f,[]);
        makeParamControl(fig,gg,step,f,val);
    end
end
end

function makeParamControl(fig,parent,step,field,val)
%makeParamControl  One control sized to the field's type (enum/logical/cell/num).
if isfield(step.enums,field)
    items = step.enums.(field);
    v = char(string(val)); if ~any(strcmp(v,items)), v = items{1}; end
    uidropdown(parent,'Items',items,'Value',v, ...
        'ValueChangedFcn',@(o,~)onSettingEdit(fig,step.id,field,o.Value));
elseif islogical(defaultType(step,field))
    uicheckbox(parent,'Text','','Value',logical(firstTrue(val)), ...
        'ValueChangedFcn',@(o,~)onSettingEdit(fig,step.id,field,o.Value));
elseif iscell(val)
    uieditfield(parent,'text','Value',cellToStr(val), ...
        'ValueChangedFcn',@(o,~)onSettingEdit(fig,step.id,field,o.Value), ...
        'Tooltip','comma-separated list');
else
    uieditfield(parent,'text','Value',val2str(val), ...
        'ValueChangedFcn',@(o,~)onSettingEdit(fig,step.id,field,o.Value));
end
end
function t = defaultType(step,field)
%defaultType  A representative value of the field from the step's default preset.
t = 0;
if isfield(step.presets,'default') && isfield(step.presets.default,field)
    t = step.presets.default.(field);
end
end
function v = firstTrue(x), if islogical(x)&&~isempty(x), v = x(1); elseif isnumeric(x)&&~isempty(x), v = x(1)~=0; else, v = false; end, end

function onSettingEdit(fig,stepId,field,rawval)
%onSettingEdit  Route an edit through wbSettingsModel + wbInvalidate + re-render.
app = getApp(fig);
step = stepById(app.reg,stepId);
value = coerceValue(step,field,rawval);
if ismember(field,step.sharedKeys)
    app.sm = wbSettingsModel('setShared', app.sm, field, value);
    seed = field;                                  % shared: invalidate every reader
else
    app.sm = wbSettingsModel('setStep', app.sm, stepId, field, value);
    seed = stepId;                                 % per-step: invalidate this step only
end
setApp(fig,app);
recomputeBase(fig);                                % new fingerprint (edited step may go stale)
applyInvalidation(fig,seed);                       % forward cascade over done cells
refreshCells(fig,{});
% refresh the open panel so the control shows the coerced value
selectStep(fig,app.selStep);
wbLog(fig,sprintf('edit %s.%s -> %s',stepId,field,val2str(value)));
end

function applyInvalidation(fig,seed)
%applyInvalidation  Mark every DONE cell in the forward set of the edit stale.
app = getApp(fig);
if isempty(app.modelArr), return; end
cells = wbInvalidate(app.reg, seed, app.modelArr);
for i = 1:size(cells,1)
    id = cells{i,1}; stepId = cells{i,2};
    if ~isKey(app.base,id), continue; end
    bs = app.base(id);
    if isfield(bs,stepId) && strcmp(bs.(stepId),'done')
        app.stale(cellKey(id,stepId)) = true;
    end
end
setApp(fig,app);
end

function v = getSetting(fig,stepId,field)
app = getApp(fig);
s = wbSettingsModel('resolve', app.sm, stepById(app.reg,stepId));
v = getfieldOr(s,field,[]);
end

function value = coerceValue(step,field,rawval)
%coerceValue  Turn a control's raw value into the typed settings value.
if isfield(step.enums,field)
    value = char(rawval); return
end
if islogical(defaultType(step,field))
    value = logical(firstTrue(rawval)); return
end
if iscell(defaultType(step,field)) || (ischar(rawval) && contains(rawval,','))
    value = strToCell(rawval); return
end
if ischar(rawval) || isstring(rawval)
    value = str2val(char(rawval));
else
    value = rawval;
end
end

%% ===================== presets ===================================== %%
function uiSavePreset(fig)
app = getApp(fig);
ensureDir(app.presetDir);
[f,p] = uiputfile('*.mat','Save preset',fullfile(app.presetDir,'preset.mat'));
if isequal(f,0), return; end
savePreset(fig,fullfile(p,f));
end
function uiLoadPreset(fig)
app = getApp(fig);
[f,p] = uigetfile('*.mat','Load preset',defaultDir(app.presetDir));
if isequal(f,0), return; end
loadPreset(fig,fullfile(p,f));
end
function savePreset(fig,pth)
app = getApp(fig);
wbSettingsModel('save', app.sm, pth);
app.presetRef = pth; setApp(fig,app);
refreshPresetDrop(fig);
wbLog(fig,['saved preset ' pth]);
end
function loadPreset(fig,pth)
app = getApp(fig);
app.sm = wbSettingsModel('load', pth);
app.presetRef = pth; setApp(fig,app);
recomputeBase(fig);
% an edit changes fingerprints everywhere; clearest is to reset the stale overlay
app = getApp(fig); app.stale = containers.Map('KeyType','char','ValueType','any'); setApp(fig,app);
refreshCells(fig,{});
selectStep(fig,getApp(fig).selStep);
refreshPresetDrop(fig);
wbLog(fig,['loaded preset ' pth]);
end
function refreshPresetDrop(fig)
app = getApp(fig); c = app.c.process;
items = {'(launcher defaults)'};
if isfolder(app.presetDir)
    d = dir(fullfile(app.presetDir,'*.mat'));
    items = [items, {d.name}];
end
c.presetDrop.Items = items;
if ~isempty(app.presetRef)
    [~,nm,ex] = fileparts(app.presetRef);
    if any(strcmp([nm ex],items)), c.presetDrop.Value = [nm ex]; end
end
end
function onPresetPick(fig,name)
if strcmp(name,'(launcher defaults)')
    app = getApp(fig); app.sm = wbSettingsModel('new'); setApp(fig,app);
    recomputeBase(fig);
    app = getApp(fig); app.stale = containers.Map('KeyType','char','ValueType','any'); setApp(fig,app);
    refreshCells(fig,{}); selectStep(fig,getApp(fig).selStep);
    wbLog(fig,'reset to launcher defaults');
    return
end
loadPreset(fig,fullfile(getApp(fig).presetDir,name));
end
function ensureDir(d), if ~isfolder(d), mkdir(d); end, end

%% ===================== session ===================================== %%
function uiSaveSession(fig)
[f,p] = uiputfile('*.mat','Save session','workbench_session.mat');
if isequal(f,0), return; end
saveSessionTo(fig,fullfile(p,f));
end
function uiLoadSession(fig)
[f,p] = uigetfile('*.mat','Load session');
if isequal(f,0), return; end
loadSessionFrom(fig,fullfile(p,f));
end
function saveSessionTo(fig,pth)
app = getApp(fig);
session = struct();
session.fNames        = app.disc.fNames;
session.referenceMode = app.disc.referenceMode;
session.groupNames    = app.groupNames;
session.modality      = app.modality;
session.rowOrder      = [];
session.bag           = app.sm.bag;
session.stepOverrides = app.sm.stepOverrides;
session.fileOverrides = app.sm.fileOverrides;
session.checked       = app.checked;
session.staleOverlay  = app.stale;
session.presetRef     = app.presetRef;
wbSession('save', pth, session);
wbLog(fig,['saved session ' pth]);
end
function loadSessionFrom(fig,pth)
session = wbSession('load', pth);
app = getApp(fig);
% rebuild the discovery grid + rows from the stored paths (no disk scan needed)
disc = discFromGrid(session.fNames, session.referenceMode, session.groupNames);
app.disc = disc;
[app.rows, app.groupNames, app.modelArr] = flattenDisc(disc);
app.modality = session.modality;
app.reg = wbStepRegistry(app.modality);
% restore the settings model
app.sm = wbSettingsModel('new', session.bag);
app.sm.stepOverrides = session.stepOverrides;
foK = keys(session.fileOverrides); foV = values(session.fileOverrides);
for i = 1:numel(foK), app.sm.fileOverrides(foK{i}) = foV{i}; end
% restore the session overlay
app.checked = session.checked;
app.stale   = session.staleOverlay;
app.presetRef = session.presetRef;
setApp(fig,app);
recomputeBase(fig);                                % re-derive the disk baseline
refreshFileTable(fig);
refreshStepSelector(fig);
renderMatrix(fig); selectStep(fig,getApp(fig).selStep);
refreshStatus(fig); refreshPresetDrop(fig);
wbLog(fig,['loaded session ' pth]);
end

function disc = discFromGrid(fNames,referenceMode,groupNames)
%discFromGrid  Rebuild a wbDiscoverFiles-shaped struct from a stored grid.
[nr,nc] = size(fNames);
models = cell(nr,nc);
flat = wbFileModel('x.rls'); flat(1) = [];
groups = struct('name',{},'rowIndex',{});
for r = 1:nr
    if numel(groupNames)>=r, gname = groupNames{r}; else, gname = sprintf('group%d',r); end
    groups(r) = struct('name',gname,'rowIndex',r);
    for cc = 1:nc
        if isempty(fNames{r,cc}), continue; end
        m = wbFileModel(fNames{r,cc});
        m.group = gname;
        m.isReference = referenceMode && cc==1;
        models{r,cc} = m;
        flat(end+1) = m; %#ok<AGROW>
    end
end
disc = struct('fNames',{fNames},'models',{models},'flat',flat, ...
    'groups',groups,'referenceMode',referenceMode);
end

%% ===================== dry-run Run ================================= %%
function entries = buildRunOrder(fig)
%buildRunOrder  Ordered checked cells: group -> reference-first -> row -> step.
app = getApp(fig); reg = app.reg;
entries = struct('groupIdx',{},'rowInGroup',{},'stepIdx',{},'stepId',{}, ...
    'stepLabel',{},'identity',{},'label',{},'group',{},'isRef',{},'arity',{});
seenPG = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(app.rows)
    r = app.rows(i);
    for s = 1:numel(reg)
        step = reg(s);
        state = resolveCellState(app,r.identity,step.id);
        if ~strcmp(state,'checked'), continue; end
        if strcmp(step.arity,'perGroup')
            k = [num2str(r.groupIdx) '||' step.id];
            if isKey(seenPG,k), continue; end       % once per group, at the reference position
            seenPG(k) = true;
            e = mkEntry(r.groupIdx,1,s,step,groupRef(app,r.groupIdx),r.group,true);
        else
            e = mkEntry(r.groupIdx,r.rowInGroup,s,step,r.identity,r.group,r.isRef);
        end
        entries(end+1) = e; %#ok<AGROW>
    end
end
if isempty(entries), return; end
K = [[entries.groupIdx]', [entries.rowInGroup]', [entries.stepIdx]'];
[~,ord] = sortrows(K);
entries = entries(ord);
end
function e = mkEntry(gi,rig,si,step,identity,group,isRef)
lbl = identity; [~,nm] = fileparts(identity); if ~isempty(nm), lbl = nm; end
e = struct('groupIdx',gi,'rowInGroup',rig,'stepIdx',si,'stepId',step.id, ...
    'stepLabel',step.label,'identity',identity,'label',lbl,'group',group, ...
    'isRef',isRef,'arity',step.arity);
end
function id = groupRef(app,groupIdx)
id = '';
for i = 1:numel(app.rows)
    if app.rows(i).groupIdx==groupIdx, id = app.rows(i).identity; return; end
end
end

function lines = dryRun(fig)
entries = buildRunOrder(fig);
lines = {sprintf('=== DRY RUN: %d step(s) would execute (no wrapper is called) ===',numel(entries))};
for i = 1:numel(entries)
    e = entries(i);
    if strcmp(e.arity,'perGroup'), who = ['GROUP ' e.group]; else, who = e.label; end
    lines{end+1} = sprintf('  %2d. [%s] %-26s :: %s', i, e.group, who, e.stepLabel); %#ok<AGROW>
end
if numel(entries)==0
    lines{end+1} = '  (nothing checked - tick cells in the matrix first)';
end
setLog(fig,lines);
end

%% ===================== real Run (executor) ========================= %%
function runChecked(fig)
%runChecked  Execute every checked cell in dependency order via wbExecutor.
app = getApp(fig);
if isfield(app,'running') && app.running, return; end       % already running
entries = buildRunOrder(fig);
if isempty(entries)
    setLog(fig,{'Nothing checked - tick cells in the matrix first.'});
    return
end
% a fresh run starts with a clean transient overlay
app.runState = containers.Map('KeyType','char','ValueType','any');
app.cellMsg  = containers.Map('KeyType','char','ValueType','any');
app.running  = true;  app.cancel = false;
setApp(fig,app);
setRunningUI(fig,true);
cleaner = onCleanup(@() finishRun(fig)); % restore UI + fold results even on error
ctx = buildExecContext(fig);
wbExecutor(entries, ctx);
end

function ctx = buildExecContext(fig)
%buildExecContext  The callback bundle wbExecutor drives the figure through.
ctx = struct();
ctx.reg           = getApp(fig).reg;
ctx.modelOf       = @(id) modelByIdentity(fig,id);
ctx.groupModels   = @(gi) groupModelArr(fig,gi);
ctx.resolve       = @(step,mdl) wbSettingsModel('resolve', getApp(fig).sm, step, mdl);
ctx.contrastStage = @(mdl) contrastStageOf(fig,mdl);
ctx.setState      = @(id,sid,state,msg) execSetState(fig,id,sid,state,msg);
ctx.progress      = @(id,sid,f,l) execProgress(fig,id,sid,f,l);
ctx.log           = @(msg) wbLog(fig,msg);
ctx.isCancelled   = @() execIsCancelled(fig);
ctx.modalGuard    = @(fcn) wbModalGuard(fig,fcn);
ctx.afterDone     = @(id,sid,mdl) execAfterDone(fig,id,sid,mdl);
end

function execSetState(fig,identity,stepId,state,msg)
%execSetState  Set a cell's transient run overlay ('' clears it) and redraw it.
if ~isvalid(fig), return; end
if nargin<5, msg = ''; end
app = getApp(fig);
key = cellKey(identity,stepId);
if isempty(state)
    if isKey(app.runState,key), remove(app.runState,key); end
    if isKey(app.cellMsg,key),  remove(app.cellMsg,key);  end
else
    app.runState(key) = state;
    if ~isempty(msg), app.cellMsg(key) = msg; end
end
setApp(fig,app);
refreshCells(fig,{identity});
drawnow limitrate;
end

function execProgress(fig,identity,stepId,frac,label)
%execProgress  Update the running cell's percent + the global progress label.
if ~isvalid(fig), return; end
if nargin<5 || ~ischar(label), label = ''; end
app = getApp(fig);
pct = max(0,min(100,round(100*frac)));
key = cellKey(identity,stepId);
if isKey(app.cellComp,key)
    comp = app.cellComp(key);
    if isgraphics(comp) && isprop(comp,'Text'), comp.Text = sprintf('%d%%',pct); end
end
if isfield(app.c.process,'progLabel') && isgraphics(app.c.process.progLabel)
    app.c.process.progLabel.Text = sprintf('%s  %s (%d%%)', shortId(identity), label, pct);
end
drawnow limitrate;
end

function tf = execIsCancelled(fig)
%execIsCancelled  The cooperative cancel flag (read safely from a hook).
tf = false;
try
    if isvalid(fig)
        a = getApp(fig);
        tf = isfield(a,'cancel') && a.cancel;
    end
catch
end
end

function execAfterDone(fig,identity,stepId,model)
%execAfterDone  A cell finished: clear its own staleness, push this recording's
%   downstream done cells to stale, drop it from the queue, and show its reports.
if ~isvalid(fig), return; end
app = getApp(fig);
selfKey = cellKey(identity,stepId);
if isKey(app.stale,selfKey), remove(app.stale,selfKey); end
cells = wbInvalidate(app.reg, stepId, model);
for i = 1:size(cells,1)
    sp = cells{i,2};
    if strcmp(sp,stepId), continue; end
    if isKey(app.base,identity)
        bs = app.base(identity);
        if isfield(bs,sp) && strcmp(bs.(sp),'done'), app.stale(cellKey(identity,sp)) = true; end
    end
end
if isKey(app.checked,selfKey), remove(app.checked,selfKey); end
setApp(fig,app);
showReports(fig,identity,stepId);
end

function finishRun(fig)
%finishRun  Reset the run UI, fold successful cells into the disk baseline, and
%   keep only error overlays (errors are not an on-disk state).
if ~isvalid(fig), return; end
app = getApp(fig);
wasCancel = isfield(app,'cancel') && app.cancel;
app.running = false; app.cancel = false;
setApp(fig,app);
recomputeBase(fig);
app = getApp(fig);
ks = keys(app.runState);
for i = 1:numel(ks)
    if ~strcmp(app.runState(ks{i}),'error'), remove(app.runState,ks{i}); end
end
ck = keys(app.checked);
for i = 1:numel(ck)
    parts = strsplit(ck{i},'||');
    if numel(parts)<2, continue; end
    if strcmp(resolveCellState(app,parts{1},parts{2}),'done'), remove(app.checked,ck{i}); end
end
setApp(fig,app);
setRunningUI(fig,false);
refreshCells(fig,{});
c = getApp(fig).c.process;
if isfield(c,'progLabel') && isgraphics(c.progLabel)
    c.progLabel.Text = ternary(wasCancel,'Stopped.','Done.');
end
drawnow limitrate;
end

function cancelRun(fig)
%cancelRun  Request a cooperative stop after the current cell.
app = getApp(fig);
if isfield(app,'running') && app.running
    app.cancel = true; setApp(fig,app);
    if isfield(app.c.process,'progLabel') && isgraphics(app.c.process.progLabel)
        app.c.process.progLabel.Text = 'Stopping after the current cell...';
    end
    wbLog(fig,'Stop requested - halting after the current cell.');
end
end

function setRunningUI(fig,running)
%setRunningUI  Toggle the Run/Preview/Stop buttons for a run in progress.
if ~isvalid(fig), return; end
c = getApp(fig).c.process;
onoff = ternary(running,'off','on');
if isfield(c,'runBtn')     && isgraphics(c.runBtn),     c.runBtn.Enable     = onoff; end
if isfield(c,'previewBtn') && isgraphics(c.previewBtn), c.previewBtn.Enable = onoff; end
if isfield(c,'stopBtn')    && isgraphics(c.stopBtn),    c.stopBtn.Enable    = ternary(running,'on','off'); end
drawnow limitrate;
end

function m = modelByIdentity(fig,identity)
%modelByIdentity  The wbFileModel of a loaded row, by recording identity.
app = getApp(fig); m = [];
for i = 1:numel(app.rows)
    if strcmp(app.rows(i).identity,identity), m = app.rows(i).model; return; end
end
end

function ms = groupModelArr(fig,groupIdx)
%groupModelArr  A group's models ordered reference-first (rowInGroup order).
app = getApp(fig);
sel = app.rows([app.rows.groupIdx]==groupIdx);
if isempty(sel), ms = []; return; end
[~,ord] = sort([sel.rowInGroup]);
sel = sel(ord);
ms = [sel.model];
end

function st = contrastStageOf(fig,model)
%contrastStageOf  't' or 's' - which contrast flag this project's files carry.
st = 't';
app = getApp(fig);
cstep = stepById(app.reg,'contrast');
if isempty(cstep), return; end
s = wbSettingsModel('resolve', app.sm, cstep, model);
if isfield(s,'contrastType') && strcmpi(char(string(s.contrastType)),'spatial'), st = 's'; end
end

%% ===================== reports (artifacts) ========================= %%
function showReports(fig,identity,stepId)
%showReports  Populate the Reports panel with a done cell's artifact images.
if ~isvalid(fig), return; end
app = getApp(fig);
if ~isfield(app.c.process,'reportPanel') || ~isgraphics(app.c.process.reportPanel), return; end
panel = app.c.process.reportPanel;
delete(panel.Children);
model = modelByIdentity(fig,identity);
step  = stepById(app.reg,stepId);
if isempty(model) || isempty(step)
    uilabel(uigridlayout(panel,[1 1],'Padding',[8 8 8 8]),'Text','No reports.','FontAngle','italic');
    return
end
files = wbArtifacts(model,step);
if isempty(files)
    gl = uigridlayout(panel,[1 1],'Padding',[8 8 8 8]);
    uilabel(gl,'Text',sprintf('%s: no report image on disk yet.',step.label),'WordWrap','on','FontAngle','italic');
    return
end
n = numel(files);
gl = uigridlayout(panel,[n+1 1],'RowHeight',[{'fit'}, repmat({150},1,n)], ...
    'RowSpacing',6,'Padding',[6 6 6 6],'Scrollable','on');
uilabel(gl,'Text',sprintf('%s - %s (click to open full size)',shortId(identity),step.label), ...
    'FontWeight','bold','WordWrap','on');
for i = 1:n, addArtifactThumb(gl,files{i}); end
end

function addArtifactThumb(parent,pth)
%addArtifactThumb  One clickable image thumbnail (falls back to a link button).
try
    uiimage(parent,'ImageSource',pth,'ScaleMethod','fit','Tooltip',pth, ...
        'ImageClickedFcn',@(~,~)openArtifactViewer(pth));
catch
    uibutton(parent,'Text',shortName(pth),'ButtonPushedFcn',@(~,~)openArtifactViewer(pth));
end
end

function openArtifactViewer(pth)
%openArtifactViewer  Open one report image full-size, preferring the DESKTOP's
%   own image viewer (winopen on Windows, open/xdg-open elsewhere).  That viewer
%   is a separate process, so the report stays resizable, zoomable and scrollable
%   while the workbench is busy inside a wrapper - a MATLAB window would only
%   repaint when the main thread yields, and uiimage offers no zoom at all.
%   Falls back to an in-MATLAB figure (axes toolbar = zoom/pan) if the hand-off
%   fails, e.g. no file association.
if isempty(pth) || ~isfile(pth), return; end
if ~openInDesktopViewer(pth)
    showArtifactInMatlab(pth);
end
end

function ok = openInDesktopViewer(pth)
%openInDesktopViewer  Hand the file to the OS-associated application.  Returns
%   as soon as the app is launched (never blocks the run loop); ok=false lets the
%   caller fall back.  The trailing & keeps the unix shells from blocking.
try
    if ispc
        winopen(pth);                                    % ShellExecute; non-blocking
        ok = true;
    elseif ismac
        ok = system(['open "' pth '" &']) == 0;
    else
        ok = system(['xdg-open "' pth '" >/dev/null 2>&1 &']) == 0;
    end
catch
    ok = false;
end
end

function showArtifactInMatlab(pth)
%showArtifactInMatlab  Fallback viewer: a classic figure sized to the image, so
%   the axes toolbar's zoom/pan are available (uiimage has neither).
[~,nm,ex] = fileparts(pth);
try
    img = imread(pth);
catch
    v = uifigure('Name',['Report - ' nm ex],'Position',[140 140 760 640]);
    uilabel(uigridlayout(v,[1 1]),'Text',['Cannot display ' pth], ...
        'HorizontalAlignment','center');
    return
end
scr = get(0,'ScreenSize');
w = min(size(img,2), 0.9*scr(3)); h = min(size(img,1), 0.9*scr(4));
f = figure('Name',['Report - ' nm ex],'NumberTitle','off','Color','w', ...
    'Position',[140 140 max(320,w) max(240,h)]);
ax = axes('Parent',f,'Position',[0 0 1 1]);
image(ax,img); axis(ax,'image','off');
end

function s = shortId(identity)
[~,s] = fileparts(char(identity));
if isempty(s), s = char(identity); end
end
function s = shortName(p)
[~,n,e] = fileparts(char(p));
s = [n e];
end

%% ===================== log ========================================= %%
function wbLog(fig,msg)
app = getApp(fig); c = app.c.process;
v = c.log.Value; if ischar(v), v = {v}; end
v{end+1} = msg;
if numel(v)>400, v = v(end-400:end); end
c.log.Value = v; drawnow limitrate;
end
function setLog(fig,lines)
app = getApp(fig);
if ischar(lines), lines = {lines}; end
if numel(lines)>400, lines = lines(end-400:end); end
app.c.process.log.Value = lines; drawnow limitrate;
end
function v = getLog(fig)
v = getApp(fig).c.process.log.Value; if ischar(v), v = {v}; end
end

%% ===================== exit ======================================== %%
function requestExit(fig)
if ~isvalid(fig), return; end
app = getApp(fig);
if isfield(app,'running') && app.running
    app.cancel = true; setApp(fig,app);
    wbLog(fig,'Stopping... click Exit again to close.');
else
    onClose(fig);
end
end
function onClose(fig)
if isvalid(fig), delete(fig); end
end

%% ===================== small value helpers ========================= %%
function alert(fig,msg)
try
    uialert(fig,msg,'Processing Workbench');
catch
    warning(msg);
end
end
function o = ternary(c,a,b), if c, o = a; else, o = b; end, end
function v = getfieldOr(s,f,d), if isstruct(s)&&isfield(s,f), v = s.(f); else, v = d; end, end
function t = headerTip(step)
t = step.label;
if ~isempty(step.gatingField), t = [t ' (gates on settings.' step.gatingField ')']; end
end
function s = val2str(v)
if ischar(v), s = v;
elseif isstring(v), s = char(v);
elseif islogical(v), s = mat2str(v);
elseif iscell(v), s = cellToStr(v);
elseif isempty(v), s = '';
elseif isscalar(v), s = num2str(v);
else, s = mat2str(v);
end
end
function v = str2val(str)
str = strtrim(str);
if isempty(str), v = []; return; end
n = str2num(str); %#ok<ST2NM>
if ~isempty(n), v = n; else, v = str; end
end
function s = cellToStr(c)
parts = cell(1,numel(c));
for i = 1:numel(c), parts{i} = char(string(c{i})); end
s = strjoin(parts,', ');
end
function c = strToCell(str)
str = char(str);
parts = strsplit(str,',');
c = cell(1,0);
for i = 1:numel(parts)
    p = strtrim(parts{i});
    if ~isempty(p), c{end+1} = p; end %#ok<AGROW>
end
end
