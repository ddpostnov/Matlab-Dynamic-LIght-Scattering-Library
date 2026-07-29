%guiWorkbench - Processing Workbench: configure per TYPE, run, watch.
%
%   A programmatic uifigure that turns the launcher workflow into three questions:
%   WHICH files (Files tab), WHAT each recording TYPE runs (Constructor tab), and
%   then GO (Process tab).  The window is a thin controller: every decision is
%   delegated to the headless brain (wbStepRegistry, wbDiscoverFiles, wbFileModel,
%   wbStateEngine, wbPrereqs, wbTypeModel, wbTypeSelection, wbSettingsModel,
%   wbRefBranch, wbInvalidate), to wbExecutor / wbModalGuard / wbArtifacts, and to
%   wbSession.
%
%   RUN IS AN EXPANSION, NOT A SELECTION.  Pressing Run turns the per-type
%   configuration into an ordered list of (file, step) work items - each type's
%   (type, branch) rows x each of its files, plus the two per-animal steps once
%   per animal - and hands it to wbExecutor, which resolves each step's actual
%   branch products from disk (branchScope).  Nothing is ticked per file, so a
%   protocol of 200 recordings is configured in the same handful of clicks as one
%   of 6.  The list runs STEP-MAJOR, exactly as a launcher does: every recording
%   is contrasted, then every one segmented, and so on, so the steps that read
%   across a file set (registration, vessel typing, split regions) always find the
%   whole set at the same level.  Progress streams back through the hook seam into a READ-ONLY
%   state table with cooperative Stop and continue-on-error; Preview order lists
%   the same plan without calling anything.  Report images the run writes collect
%   in a link list that opens them in the desktop viewer.
%
%   THE SESSION IS THE DELIVERABLE.  wbSession is written on every state change
%   and at the end of every run, so a crash, a Stop or a MATLAB restart always
%   leaves a resumable project behind - and it is the hand-off to guiExport and
%   guiExplore, which read it without this window being open.
%
%   THREE LABEL AXES.  Every file carries an ANIMAL (the subject - the scope of
%   registration, vessel typing and the reference recording), a TYPE (its
%   experimental role) and an experimental GROUP (a comparison label Export and
%   Explore use and processing ignores), plus a recording INDEX.  They are
%   independent - a group may span animals and an animal may span groups - and
%   none has a fixed vocabulary: each is a regexp match over the file name, hand-
%   overridable in the Files table (wbTypeModel).  Each animal owns at most ONE
%   reference RECORDING (stored as an identity, so each step still resolves the
%   branch it needs), pinned by the reference regexp or by hand.
%
%   CONFIGURATION IS A PROPERTY OF THE TYPE - AND OF THE PRODUCT.  The steps that
%   read the raw recording (contrast, internal cycle) each write a NEW, independent
%   triplet, and everything later APPENDS to one of them, so one recording can
%   carry TWO independent result sets.  The Constructor therefore asks two
%   questions: which RAW steps a type runs (top, once per type - each tick creates
%   that type's product row), and which DERIVED steps each (type, product) row runs
%   (bottom, one row per pipeline, cells greyed with their reason where a branch
%   does not offer the step, and shown as inherited where a step is drawn once and
%   copied to the other branches).  Both halves of a row label are data: the type
%   token comes from the Files tab and the product flag from the producing step plus
%   that type's own settings ('_t_K', '_s_K', '_c_K').  Settings stay keyed by
%   (step, TYPE) - two rows of one type share them - so a divergent animal is
%   simply a second type.  Rows are built FROM THE DATA: two types or eleven work
%   with no code change, and nothing anywhere branches on a type's name.  The two
%   per-animal steps (registration, vessel typing) span the animal instead and get
%   their own box; which BRANCH FILE of each animal's reference recording every
%   reference-taking step resolves to is reported in the selection summary.
%
%   Tabs: 1 Files (scan a root recursively for an extension - or add files by
%                    hand - then CURATE: sort, delete, label every file (animal /
%                    type / index / group / modality, in place or in bulk over a
%                    selection), and pin ONE reference recording per animal.  The
%                    other tabs stay locked until every field is assigned and the
%                    file names are unique - see fileProblems) *
%         2 Constructor (type x raw step on top with the per-animal steps beside
%                    it, one (type, product) row per pipeline below, and the
%                    per-(step,type) settings on the right; the summary carries the
%                    resolved reference branches and the warnings) *
%         3 Process (READ-ONLY: one state table, rows = files and columns =
%                    steps, showing '.'/queued/running NN%/done/error/skipped; the
%                    selected row's file below it; the run log; and the result
%                    links on the right.  All selection lives in the Constructor -
%                    spec D6/D7) *
%         4 Export  (a sheet/format selection UI over exportToExcel) *
%         5 Explore (GUIs/workbench/guiExplore hosted in-tab, seeded from the
%                    loaded files/animals) *
%
% Syntax:
%    guiWorkbench                       % open the workbench
%    h = guiWorkbench('Visible','off')  % headless (for programmatic drive/tests)
%
%   getappdata(h,'workbenchAPI') exposes a struct of function handles that drive
%   the same internal logic as the UI (used by testWorkbenchShell).
%
% See also: wbStepRegistry, wbDiscoverFiles, wbTypeModel, wbTypeSelection,
%           wbFileModel, wbStateEngine, wbSettingsModel, wbInvalidate,
%           wbRefBranch, wbExecutor, wbModalGuard, wbArtifacts, wbSession,
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
app.rows         = emptyRows();                    % flat, animal-major, reference-first
app.animalNames   = {};
% ---- the curated working set (Files tab) ----
app.files        = emptyFiles();                   % ONE ENTRY PER FILE (rows dedup by identity)
app.root         = '';                             % scan root
app.glob         = '*.rls';                        % scan glob: the raw recordings
app.patterns     = wbTypeModel('emptyPatterns');   % animal/type/index/expGroup regexps
app.patterns.ref = '';                             % + the reference regexp (not a label axis)
app.overrides    = wbTypeModel('emptyOverrides');  % hand labels, path -> value, per axis
app.modalityOvr  = containers.Map('KeyType','char','ValueType','char'); % path -> modality
app.labelsAuto   = [];                             % regexp-only labels (to spot a real override)
app.labels       = [];                             % labels actually in force
app.autoRef      = containers.Map('KeyType','char','ValueType','char'); % animal -> identity (regexp)
app.animalRefMan = containers.Map('KeyType','char','ValueType','char'); % animal -> identity (hand)
app.animalRef    = containers.Map('KeyType','char','ValueType','char'); % animal -> identity (effective)
app.modelArr     = wbFileModel('x.rls'); app.modelArr(1) = [];   % 1x0 model array
app.base         = containers.Map('KeyType','char','ValueType','any');   % identity -> struct(stepId->state)
app.fileState    = containers.Map('KeyType','char','ValueType','any');   % PATH     -> struct(stepId->state)  (D6)
app.branchState  = containers.Map('KeyType','char','ValueType','any');   % 'identity||branch' -> struct(...)   (D6)
app.checked      = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> true
app.stale        = containers.Map('KeyType','char','ValueType','any');   % session edit overlay
app.runState     = containers.Map('KeyType','char','ValueType','any');   % transient run overlay: 'running'|'done'|'error'
app.cellMsg      = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> error/tooltip text
app.cellPct      = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> 0..100 while running
app.completed    = containers.Map('KeyType','char','ValueType','any');   % 'path||stepId' -> completion record (session)
app.pRows        = emptyProgressRows();            % the monitor's rows: ONE PER FILE
app.pRowOf       = containers.Map('KeyType','char','ValueType','double');% path -> progress-table row
app.pSel         = 0;                              % the monitor row whose detail is shown
app.results      = emptyResults();                 % report images this session produced
app.sessionPath  = '';                             % where the durable session is autosaved
app.selStep      = app.reg(1).id;                  % step last selected programmatically
% ---- the Constructor: configuration per recording TYPE (never per file) ----
app.typeSel      = wbTypeSelection('new');         % 'type||branch||stepId' -> true
app.animalSel    = containers.Map('KeyType','char','ValueType','any'); % stepId -> true
app.selType      = '';                             % type shown in the Constructor panel
app.selCStep     = '';                             % step shown in the Constructor panel
app.presetDir    = presetDirDefault();
app.presetRef    = '';
app.running      = false;
app.cancel       = false;
app.exportSrcs   = struct('path',{},'rpath',{},'animal',{},'label',{}); % exportable BFI sources
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
app.tabs.files       = uitab(tg,'Title','1 - Files');
app.tabs.constructor = uitab(tg,'Title','2 - Constructor');
app.tabs.process     = uitab(tg,'Title','3 - Process');
app.tabs.export      = uitab(tg,'Title','4 - Export');
app.tabs.explore     = uitab(tg,'Title','5 - Explore');
app.tg = tg;
setappdata(fig,'app',app);

buildFilesTab(fig);
buildConstructorTab(fig);
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
    'animals',     @() getApp(fig).animalNames, ...
    ... % ---- Files tab: scan, curate, label, reference ----
    'setSource',   @(root,glob) setSource(fig,root,glob), ...
    'setPattern',  @(axis,rx) setPattern(fig,axis,rx), ...
    'patterns',    @() getApp(fig).patterns, ...
    'scan',        @() doScan(fig), ...
    'addPaths',    @(p) addPaths(fig,p), ...
    'deletePaths', @(p) deletePaths(fig,p), ...
    'files',       @() getApp(fig).files, ...
    'filePaths',   @() {getApp(fig).files.path}, ...
    'fileTable',   @() fileTableData(getApp(fig)), ...
    'setLabel',    @(p,axis,v) setLabel(fig,p,axis,v), ...
    'labelValues', @(axis) labelValues(getApp(fig),axis), ...
    'setModality', @(p,v) setModality(fig,p,v), ...
    'quickAssign', @(p,field,v) quickAssign(fig,p,field,v), ...
    'selectRows',  @(p) selectRows(fig,p), ...
    'selectedPaths',@() selectedPaths(getApp(fig),getApp(fig).c.files), ...
    'fileProblems',@() fileProblems(getApp(fig)), ...
    'filesValid',  @() filesValid(getApp(fig)), ...
    'filesBanner', @() getApp(fig).c.files.problems.Text, ...
    'setReference',@(p) setReference(fig,p,true), ...
    ... % ---- Constructor tab: raw producers, (type,branch) rows, animal steps ----
    'constructorTypes',@() constructorTypes(getApp(fig)), ...
    'typeStepIds', @() wbTypeSelection('typeSteps',   getApp(fig).reg), ...
    'rawStepIds',  @() wbTypeSelection('rawSteps',    getApp(fig).reg), ...
    'derivedStepIds',@() wbTypeSelection('derivedSteps',getApp(fig).reg), ...
    'animalStepIds',@() wbTypeSelection('animalSteps',getApp(fig).reg), ...
    'branchIds',   @() wbTypeSelection('branches',    getApp(fig).reg), ...
    'typeFileCount',@(ty) typeFileCount(getApp(fig),ty), ...
    'tickRaw',     @(ty,stepId,tf) tickRaw(fig,ty,stepId,tf), ...
    'rawChecked',  @(ty,stepId) rawChecked(getApp(fig),ty,stepId), ...
    'tickRow',     @(ty,br,stepId,tf) tickRow(fig,ty,br,stepId,tf), ...
    'rowChecked',  @(ty,br,stepId) wbTypeSelection('isOn',getApp(fig).typeSel,ty,br,stepId), ...
    'rowShown',    @(ty,br,stepId) rowShown(getApp(fig),ty,br,stepId), ...
    'rowInherited',@(ty,br,stepId) rowInherited(getApp(fig),ty,br,stepId), ...
    'rowSelection',@(ty,br) wbTypeSelection('steps',getApp(fig).typeSel,getApp(fig).reg,ty,br), ...
    'rowInheritedSteps',@(ty,br) wbTypeSelection('inherited',getApp(fig).typeSel,getApp(fig).reg,ty,br), ...
    'rowOffers',   @(br,stepId) wbTypeSelection('offers',getApp(fig).reg,br,stepId), ...
    'rowWhyNot',   @(br,stepId) wbTypeSelection('why',getApp(fig).reg,br,stepId), ...
    'constructorRows',@() constructorRows(getApp(fig)), ...
    'rowFlag',     @(ty,br) rowFlagFor(getApp(fig),ty,br), ...
    'typeRows',    @(ty) wbTypeSelection('rows',getApp(fig).typeSel,getApp(fig).reg,ty), ...
    'derivedEnabled',@() ~isempty(constructorRows(getApp(fig))), ...
    'tickAnimalStep',@(stepId,tf) tickAnimalStep(fig,stepId,tf), ...
    'animalStepChecked',@(stepId) isKey(getApp(fig).animalSel,stepId), ...
    'animalPlan',  @() animalStepPlan(getApp(fig)), ...
    'animalPlanLines',@() animalPlanLines(getApp(fig)), ...
    'constructorWarnings',@() constructorWarnings(getApp(fig)), ...
    'summaryLines',@() summaryLines(getApp(fig)), ...
    'editTypeSetting',@(ty,stepId,field,value) onTypeSettingEdit(fig,ty,stepId,field,value), ...
    'getTypeSetting',@(ty,stepId,field) getTypeSetting(fig,ty,stepId,field), ...
    'copyTypeConfig',@(src,dst) copyTypeConfig(fig,src,dst), ...
    'presetNames', @() wbTypePresets('names'), ...
    'applyPreset', @(name,ty) applyPreset(fig,name,ty), ...
    'resetTypeConfig',@(ty) resetTypeConfig(fig,ty), ...
    'selectType',  @(ty) selectConstructorType(fig,ty), ...
    'clearReference',@(p) setReference(fig,p,false), ...
    'animalRef',   @(a) animalRefOf(getApp(fig),a), ...
    'animalRefs',  @() getApp(fig).animalRef, ...
    'filesStatus', @() getApp(fig).c.files.status.Text, ...
    'tabsLocked',  @() ~filesValid(getApp(fig)), ...
    'nextEnabled', @() strcmp(getApp(fig).c.files.nextBtn.Enable,'on'), ...
    'goNext',      @() guardTabSwitchTo(fig,getApp(fig).tabs.constructor), ...
    'currentTab',  @() getApp(fig).tg.SelectedTab.Title, ...
    'cellState',   @(id,stepId) resolveCellState(getApp(fig),id,stepId), ...
    'check',       @(id,stepId,tf) apiCheck(fig,id,stepId,tf), ...
    'checkColumn', @(stepId,tf) checkColumn(fig,stepId,tf), ...
    'checkAnimal', @(a,tf) checkAnimal(fig,a,tf), ...
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
    ... % ---- Process tab: the read-only monitor, its results, the session ----
    'progressRows',@() {getApp(fig).pRows.path}, ...
    'progressData',@() getApp(fig).c.process.stateTable.Data, ...
    'progressCell',@(p,stepId) progressCellOf(fig,p,stepId), ...
    'fileState',   @(p,stepId) fileStateOf(getApp(fig),p,stepId), ...
    'plannedSteps',@(p) plannedStepsOf(fig,p), ...
    'resultLinks', @() resultLinks(fig), ...
    'sessionPath', @() getApp(fig).sessionPath, ...
    'autosave',    @() autosaveSession(fig), ...
    'completed',   @() getApp(fig).completed, ...
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

renderProgress(fig);
renderConstructor(fig);
if nargout>0, h = fig; end
end

%% ===================== app-state helpers ============================ %%
function app = getApp(fig), app = getappdata(fig,'app'); end
function setApp(fig,app), setappdata(fig,'app',app); end
function d = emptyDisc()
d = struct('fNames',{{}},'models',{{}},'flat',[],'animals',[],'referenceMode',false, ...
    'patterns',wbTypeModel('emptyPatterns'));
end
function r = emptyRows()
r = struct('model',{},'identity',{},'animal',{},'animalIdx',{},'rowInAnimal',{},'isRef',{},'label',{});
end
function f = emptyFiles()
%emptyFiles  The curated working set: one entry per FILE (see buildFileEntries).
f = struct('model',{},'path',{},'name',{},'animal',{},'type',{},'index',{}, ...
    'expGroup',{},'modality',{},'isRef',{});
end
function r = emptyProgressRows()
%emptyProgressRows  The Process tab's rows - one per FILE, run-ordered.
r = struct('path',{},'model',{},'identity',{},'label',{},'animal',{},'type',{}, ...
    'expGroup',{},'branch',{},'isRef',{},'animalIdx',{},'rowInAnimal',{});
end
function r = emptyResults()
%emptyResults  The result-link list (spec D5): report images, newest first.
r = struct('path',{},'label',{},'when',{});
end
function d = presetDirDefault()
d = fullfile(prefdir,'guiWorkbenchPresets');
end
function s = stepById(reg,id), s = reg(strcmp({reg.id},id)); end
function k = cellKey(identity,stepId), k = [identity '||' stepId]; end

%% ===================== FILES tab ==================================== %%
function buildFilesTab(fig)
%buildFilesTab  The curation surface: scan a tree (or pick by hand), then sort,
%   delete and label every file on the animal / type / index / group axes and pin
%   ONE reference recording per animal.
app = getApp(fig); t = app.tabs.files;
gl = uigridlayout(t,[1 1],'Padding',[6 6 6 6]);
%   EIGHT rows, one per child: a child beyond the declared row count lands in an
%   auto-added row of height '1x', which is what blew the session buttons up.
lp = uigridlayout(gl,[8 1], ...
    'RowHeight',{'fit','fit','fit','1x','fit','fit','fit','fit'},'RowSpacing',6);

uilabel(lp,'Text',['Scan a root folder for an extension (or add files by hand), then curate: ' ...
    'sort, delete, label and pin one reference recording per ANIMAL.'],'FontWeight','bold','WordWrap','on');

% -- source row: root + glob + the loaders --
sp = uigridlayout(lp,[1 9], ...
    'ColumnWidth',{'fit','2x','fit','fit','0.9x','fit','fit','fit','fit'},'Padding',[0 0 0 0]);
uilabel(sp,'Text','root');
c.root = uieditfield(sp,'text','Value','','Tooltip','folder searched RECURSIVELY by Scan', ...
    'ValueChangedFcn',@(s,~)onSourceEdit(fig));
uibutton(sp,'Text','Browse...','ButtonPushedFcn',@(~,~)uiBrowseRoot(fig));
uilabel(sp,'Text','files');
c.glob = uieditfield(sp,'text','Value',app.glob,'ValueChangedFcn',@(s,~)onSourceEdit(fig), ...
    'Tooltip',tipGlob());
uibutton(sp,'Text','Scan','BackgroundColor',[0.82 0.92 0.82],'ButtonPushedFcn',@(~,~)uiScan(fig), ...
    'Tooltip','recurse the root for the glob and REPLACE the working set');
uibutton(sp,'Text','Add files...','ButtonPushedFcn',@(~,~)uiAddFiles(fig), ...
    'Tooltip','pick files by hand (multiselect) and ADD them to the working set');
uibutton(sp,'Text','Add folder...','ButtonPushedFcn',@(~,~)uiAddFolder(fig), ...
    'Tooltip','recurse another folder for the glob and ADD what it finds');
uibutton(sp,'Text','Clear','ButtonPushedFcn',@(~,~)uiClear(fig),'Tooltip','empty the working set (no file is deleted)');

% -- regexp row: one box per label axis + the reference rule --
rp = uigridlayout(lp,[2 5],'RowHeight',{'fit','fit'},'Padding',[0 0 0 0],'RowSpacing',1);
axes5 = {'animal','Animal','[A-Z]+\d+',tipAnimal(); 'type','Type','',tipType(); ...
         'index','Rec. index','',tipIndex(); 'expGroup','Group (experimental)','',tipExpGroup(); ...
         'ref','Reference','',tipRef()};
for k = 1:size(axes5,1)
    lb = uilabel(rp,'Text',axes5{k,2},'FontSize',10); lb.Layout.Row = 1; lb.Layout.Column = k;
    f  = patField(rp,fig,axes5{k,1},axes5{k,3},axes5{k,4});
    f.Layout.Row = 2; f.Layout.Column = k;
    c.pat.(axes5{k,1}) = f;
end

% -- the curation table --
%   MODALITY is read-only here on purpose: a uitable's ColumnFormat dropdown is
%   one list for the whole COLUMN, and what a file may legally be depends on ITS
%   extension (wbFileModel owns that rule).  Rather than offering illegal values
%   and refusing them afterwards, the modality is set through Quick assign below,
%   whose dropdown is narrowed to what every selected row can actually be.
c.fileTbl = uitable(lp, ...
    'ColumnName',fileTableColumns(), ...
    'ColumnFormat',repmat({'char'},1,numel(fileTableColumns())), ...
    'ColumnEditable',fileTableEditable(), ...
    'ColumnWidth',fileTableWidths(), ...
    'ColumnSortable',true(1,numel(fileTableColumns())),'RowName',{}, ...
    'SelectionType','row','Multiselect','on', ...
    'CellEditCallback',@(s,e)onFileTableEdit(fig,s,e));
c.fileTbl.ColumnFormat{1} = 'logical';              % the reference tick

% -- quick assign: give EVERY selected row the same value in one go -----------
cp = uigridlayout(lp,[1 7], ...
    'ColumnWidth',{'fit','fit','fit','1.2x','fit','fit','2x'},'Padding',[0 0 0 0]);
uilabel(cp,'Text','Quick assign: set','FontWeight','bold');
c.assignAxis = uidropdown(cp,'Items',assignableColumns(),'Value','type', ...
    'Tooltip','which field the selected rows get', ...
    'ValueChangedFcn',@(~,~)refreshAssignItems(fig));
uilabel(cp,'Text','=');
c.assignVal = uidropdown(cp,'Items',{''},'Value','','Editable','on', ...
    'Tooltip','pick a value already in use, or TYPE A NEW ONE - the vocabulary is yours');
uibutton(cp,'Text','Apply to selected rows','BackgroundColor',[0.85 0.9 1], ...
    'ButtonPushedFcn',@(~,~)uiAssignLabel(fig), ...
    'Tooltip','set this field on EVERY selected row at once (select rows by clicking / shift-clicking the table)');
uibutton(cp,'Text','Delete selected','BackgroundColor',[1 0.88 0.82], ...
    'ButtonPushedFcn',@(~,~)uiDeleteSelected(fig), ...
    'Tooltip','remove the selected rows from the WORKING SET only - nothing is ever deleted from disk');
c.curStatus = uilabel(cp,'Text','','FontAngle','italic','FontSize',10);

c.problems = uilabel(lp,'Text','','WordWrap','on','FontWeight','bold');
c.status   = uilabel(lp,'Text','No files loaded.','WordWrap','on');

% -- session row (RowHeight fit: buttons keep their natural height) --
ssp = uigridlayout(lp,[1 5],'RowHeight',{'fit'}, ...
    'ColumnWidth',{'fit','fit','1x','fit','fit'},'Padding',[0 0 0 0]);
uibutton(ssp,'Text','Save session...','ButtonPushedFcn',@(~,~)uiSaveSession(fig));
uibutton(ssp,'Text','Load session...','ButtonPushedFcn',@(~,~)uiLoadSession(fig));
uilabel(ssp,'Text','the session carries the file list, its labels, the references and the settings', ...
    'FontAngle','italic','FontSize',10);
c.nextBtn = uibutton(ssp,'Text','Next: Constructor','BackgroundColor',[0.82 0.92 0.82], ...
    'ButtonPushedFcn',@(~,~)goToConstructor(fig), ...
    'Tooltip','continue to the Constructor (enabled once every file is labelled)');
uibutton(ssp,'Text','Exit workbench','ButtonPushedFcn',@(~,~)requestExit(fig),'BackgroundColor',[1 0.82 0.82]);

app.c.files = c; setApp(fig,app);
end

function f = patField(parent,fig,axis,dflt,tip)
%patField  One regexp box, wired to the pattern struct.
f = uieditfield(parent,'text','Value',dflt,'Tooltip',tip, ...
    'ValueChangedFcn',@(s,~)onPatternEdit(fig,axis,s.Value));
end

function cols = fileTableColumns()
%fileTableColumns  The curation table's HEADERS.  'file' (the bare name, without
%   its extension) plus 'file type' (the extension) is the row KEY, so file names
%   must be unique across the scanned tree - the workbench refuses to go on when
%   they are not (see fileProblems), which is also the library's own naming rule.
cols = {'reference','file','animal','type','index','group','file type','modality','stage'};
end
function keys = fileTableKeys()
%fileTableKeys  Field-safe names for the same columns (headers carry a space).
keys = {'reference','file','animal','type','index','group','filetype','modality','stage'};
end
function e = fileTableEditable()
%fileTableEditable  What may be typed in place.  The name, its extension, the
%   parsed stage and the modality are properties OF THE FILE, not labels.
e = [true false true true true true false false false];
end
function w = fileTableWidths()
%fileTableWidths  'reference' is just wide enough for its own header; the rest
%   size to their content.
w = [{78}, repmat({'auto'},1,numel(fileTableColumns())-1)];
end
function cols = assignableColumns()
%assignableColumns  The fields Quick assign (and a direct cell edit) can set.
cols = {'animal','type','index','group','modality'};
end

%% ---- Files-tab tooltips (the guiExplore idiom: examples, not prose) ---- %%
function t = tipGlob()
t = sprintf(['dir() glob matched against the file NAME, searched recursively.\n' ...
    '   *.rls           raw LSCI recordings\n' ...
    '   *_t_K_d.mat     temporal-contrast SOURCE files\n' ...
    '   *_BFI_d.mat     every BFI branch\n' ...
    'This is a GLOB (shell wildcards), not a regexp - the boxes below are regexps.']);
end
function t = tipAnimal()
t = sprintf(['Regexp whose MATCH is the ANIMAL id - the subject.  One animal owns one\n' ...
    'reference recording and is the scope of registration / vessel typing.\n' ...
    '   [A-Z]+\\d+        PSY01, REG12, AB3 ...\n' ...
    '   m\\d+             m07, m11 ...\n' ...
    '   (mouse|rat)\\d+   a species-prefixed id\n' ...
    'No match -> the "(unassigned)" bucket.  Empty = every file in one row.']);
end
function t = tipType()
t = sprintf(['Regexp whose MATCH is the recording TYPE - its experimental role.  The type\n' ...
    'owns the processing configuration, so files of one type are processed alike.\n' ...
    '   BV|BN|BP         one lab''s slow / NVC / pulsatile recordings\n' ...
    '   ctrl|stim|washout   another lab''s protocol\n' ...
    '   \\d(?=BP)         whatever your names actually carry\n' ...
    'There is NO built-in list of types: the values are whatever this matches.\n' ...
    'No match -> the "(untyped)" bucket, which is a type like any other.']);
end
function t = tipIndex()
t = sprintf(['Regexp whose MATCH is the recording index inside an animal (a stratifier).\n' ...
    '   [aik]\\d           a condition+repeat token, e.g. i2, k1, a2\n' ...
    '   run(\\d+)          the number after "run"\n' ...
    '   \\d+(?=_c_)        the digits just before _c_\n' ...
    'Leave empty to treat every file as the same index.']);
end
function t = tipExpGroup()
t = sprintf(['Regexp whose MATCH is the EXPERIMENTAL group - a comparison label used by\n' ...
    'Export and Explore.  Processing ignores it entirely.\n' ...
    '   KO|WT             genotype taken from the name\n' ...
    '   Ctrl|Stroke       condition taken from the name\n' ...
    '   pre|post          a timepoint\n' ...
    'It is INDEPENDENT of the animal: a group may span animals and one animal\n' ...
    'may span groups.  No match -> the "(ungrouped)" bucket.']);
end
function t = tipRef()
t = sprintf(['Regexp picking the REFERENCE recording of each animal, valid for ALL of\n' ...
    'that animal''s files whatever their type or group.  Passed to getFileNamesList,\n' ...
    'which forces the match into column 1.\n' ...
    '   1BP_c_BFI_d\\.mat  the animal''s first pulsatile recording\n' ...
    '   _ref_             a name you mark by hand\n' ...
    'Applied by Scan.  You can always override it by ticking "ref" in the table -\n' ...
    'a hand-pinned reference wins, and an animal may legally have none.']);
end

%% ===================== CONSTRUCTOR tab ============================== %%
function buildConstructorTab(fig)
%buildConstructorTab  Configure the pipeline per recording TYPE - and, below, per
%   (type, BRANCH) row, because ONE RAW RECORDING CAN DRIVE TWO INDEPENDENT RESULT
%   SETS.  The steps that read the recording itself (contrast, internal cycle) each
%   write a new triplet, and everything later APPENDS to one of those products:
%     TOP    what each type computes FROM THE RECORDING (the raw producers), with
%            the two per-animal steps beside it;
%     BOTTOM what is then done to each product - one row per (type, product), so
%            "segment the cardiac branch but not the contrast one" is sayable;
%     RIGHT  the settings of one (step, type) - configuration is per TYPE, so both
%            rows of a type share it.
%   Every row is built FROM THE DATA (two types or eleven, no code change, nothing
%   branching on a token's value), and the per-animal reference plan and the
%   warnings live in the selection summary - the one place status is read.
app = getApp(fig); t = app.tabs.constructor;
gl = uigridlayout(t,[1 2],'ColumnWidth',{'2.2x','1x'},'Padding',[6 6 6 6],'ColumnSpacing',8);

% -- LEFT: header + raw panel + derived rows + summary --
% the raw panel is sized to its content by renderRawPanel (a 'fit' row cannot
% measure the scrollable grids inside it), the rest shares what is left
left = uigridlayout(gl,[4 1],'RowHeight',{'fit',170,'1.3x','1x'},'RowSpacing',6);
c.left = left;

hdr = uigridlayout(left,[2 1],'RowHeight',{'fit','fit'},'Padding',[0 0 0 0],'RowSpacing',3);
uilabel(hdr,'Text',['Tick the steps each recording TYPE runs.  A type owns its ' ...
    'processing configuration and every file of that type is processed alike - if one ' ...
    'animal needs different parameters, make it a second type on the Files tab.'], ...
    'WordWrap','on','FontWeight','bold');
tb = uigridlayout(hdr,[1 7], ...
    'ColumnWidth',{'fit','1x','fit','1x','fit','fit','1.2x'},'Padding',[0 0 0 0]);
uilabel(tb,'Text','Copy configuration from');
c.copySrc = uidropdown(tb,'Items',{''},'Value','', ...
    'Tooltip',['a standard protocol, or another type you have already configured.  ' ...
               'A protocol only POPULATES the boxes - everything stays editable.']);
uilabel(tb,'Text','to');
c.copyDst = uidropdown(tb,'Items',{''},'Value','','Tooltip','the type to overwrite');
uibutton(tb,'Text','Copy','BackgroundColor',[0.85 0.9 1],'ButtonPushedFcn',@(~,~)uiCopyType(fig), ...
    'Tooltip','give the second type the first one''s selected steps AND settings');
uibutton(tb,'Text','Reset selected type','ButtonPushedFcn',@(~,~)uiResetType(fig), ...
    'Tooltip','drop the type shown on the right back to the launcher defaults');
c.status = uilabel(tb,'Text','','FontAngle','italic','FontSize',10);

c.rawPanel = uipanel(left,'Title','1 - Raw processing: what each TYPE computes from the recording', ...
    'FontWeight','bold');
c.derivedPanel = uipanel(left,'Title','2 - Derived processing: one row per TYPE and PRODUCT', ...
    'FontWeight','bold');
c.summaryPanel = uipanel(left,'Title','Selection summary','FontWeight','bold');

% -- RIGHT: settings for one (step, type) --
right = uigridlayout(gl,[4 1],'RowHeight',{'fit','fit','fit','1x'},'RowSpacing',6);
sel = uigridlayout(right,[2 2],'ColumnWidth',{'fit','1x'},'RowHeight',{'fit','fit'}, ...
    'Padding',[0 0 0 0],'RowSpacing',3);
uilabel(sel,'Text','Type','FontWeight','bold');
c.typeDrop = uidropdown(sel,'Items',{'(no types)'},'Value','(no types)', ...
    'Tooltip','the recording type whose settings are edited below', ...
    'ValueChangedFcn',@(s,~)selectConstructorType(fig,s.Value));
uilabel(sel,'Text','Step','FontWeight','bold');
c.stepDrop = uidropdown(sel,'Items',{''},'Value','', ...
    'Tooltip','the step whose settings are shown below', ...
    'ValueChangedFcn',@(s,~)selectConstructorStep(fig,s.Value));
c.stepInfo = uilabel(right,'Text','','WordWrap','on','FontAngle','italic');
c.scopeInfo = uilabel(right,'Text','','WordWrap','on','FontSize',10);
c.paramPanel = uipanel(right,'BorderType','none');   % the section stack owns the scroll

app.c.constructor = c; setApp(fig,app);
end

%% ---- Constructor: raw producers, then the (type,branch) rows ------- %%
function renderConstructor(fig)
%renderConstructor  Rebuild every Constructor surface from the current data.
app = getApp(fig);
if ~isfield(app,'c') || ~isfield(app.c,'constructor'), return; end
renderRawPanel(fig);
renderDerivedPanel(fig);
refreshConstructorSelectors(fig);
refreshSummary(fig);
end

function types = constructorTypes(app)
%constructorTypes  The matrix rows: the type values ACTUALLY IN THE DATA, in order
%   of first appearance.  Never a list - '(untyped)' is a row like any other, and a
%   type that no file carries any more simply has no row (its configuration is kept
%   so re-typing a file brings the row back unchanged).
types = {};
if isempty(app.files) || isempty(app.labels), return; end
types = wbTypeModel('values', app.labels, 'type');
end
function n = typeFileCount(app,type)
n = 0;
if isempty(app.files), return; end
n = nnz(strcmp({app.files.type}, char(type)));
end

function renderRawPanel(fig)
%renderRawPanel  The raw producers (type x step, small) with the per-animal steps
%   beside them.  Ticking a producer is what brings its branch ROW into existence
%   below - the two panels are one decision seen twice.
app = getApp(fig); c = app.c.constructor;
delete(c.rawPanel.Children);
gl = uigridlayout(c.rawPanel,[1 2],'ColumnWidth',{'1.7x','1x'}, ...
    'Padding',[4 4 4 4],'ColumnSpacing',8);
renderRawMatrix(fig,gl);
renderAnimalBox(fig,gl);
% no dead space under two types, no clipping under eleven: the panel is exactly as
% tall as the taller of its two halves, and only then does the matrix scroll
if isfield(c,'left') && isgraphics(c.left)
    nT = numel(constructorTypes(app));
    nA = numel(wbTypeSelection('animalSteps', app.reg));
    h  = max(66 + 26*nA, 96 + 26*nT);
    rh = c.left.RowHeight; rh{2} = min(320,max(120,h)); c.left.RowHeight = rh;
end
end

function renderRawMatrix(fig,parent)
%renderRawMatrix  One checkbox per (type, raw producer).  The producers are the
%   steps whose input is the recording itself - derived from the registry, never
%   listed - and each writes its own independent product.
app = getApp(fig);
types = constructorTypes(app);
cols  = wbTypeSelection('rawSteps', app.reg);
if isempty(types) || isempty(cols)
    g = uigridlayout(parent,[1 1],'Padding',[16 16 16 16]);
    uilabel(g,'Text','Load and label recordings on the Files tab to populate the type matrix.', ...
        'HorizontalAlignment','center','FontAngle','italic');
    return
end

grid = uigridlayout(parent,[1+numel(types) 1+numel(cols)], ...
    'ColumnWidth',[{170}, repmat({120},1,numel(cols))], ...
    'RowHeight',[{50}, repmat({26},1,numel(types))], ...
    'RowSpacing',2,'ColumnSpacing',2,'Padding',[2 2 2 2],'Scrollable','on');
h = uilabel(grid,'Text','Type  \  Raw step','FontWeight','bold');
h.Layout.Row = 1; h.Layout.Column = 1;
for s = 1:numel(cols)
    step = stepById(app.reg,cols{s});
    hp = uigridlayout(grid,[2 1],'RowHeight',{'1x','fit'},'Padding',[1 1 1 1],'RowSpacing',1);
    hp.Layout.Row = 1; hp.Layout.Column = s+1;
    uilabel(hp,'Text',step.label,'WordWrap','on','FontSize',10,'Tooltip',rawStepTip(step));
    gr = uigridlayout(hp,[1 1],'Padding',[0 0 0 0]);
    uibutton(gr,'Text','settings','FontSize',9,'ButtonPushedFcn',@(~,~)selectConstructorStep(fig,step.id));
end
for r = 1:numel(types)
    ty = types{r};
    lb = uilabel(grid,'Text',sprintf('%s - %d files', ty, typeFileCount(app,ty)), ...
        'FontWeight','bold','Tooltip','every file of this type is processed alike');
    lb.Layout.Row = r+1; lb.Layout.Column = 1;
    for s = 1:numel(cols)
        step = stepById(app.reg,cols{s});
        cb = uicheckbox(grid,'Text',rowFlagFor(app,ty,step.branch), ...
            'Value',wbTypeSelection('isOn',app.typeSel,ty,step.branch,cols{s}), ...
            'FontSize',10,'Tooltip',rawCellTip(app,ty,step), ...
            'ValueChangedFcn',@(o,~)tickRaw(fig,ty,cols{s},o.Value));
        cb.Layout.Row = r+1; cb.Layout.Column = s+1;
    end
end
end

function t = rawStepTip(step)
t = sprintf(['%s reads the raw recording and writes a NEW, independent product.\n' ...
    'Ticking it gives this type a %s row below; everything done to that product is ' ...
    'ticked there.'], step.label, step.branch);
end
function t = rawCellTip(app,type,step)
t = sprintf('%s produces %s%s for every %s recording.', step.label, ...
    rowFlagFor(app,type,step.branch), ' ', type);
end

function lines = tickRaw(fig,type,stepId,tf)
%tickRaw  Switch one raw producer on/off for a TYPE - i.e. create or remove that
%   type's branch row.  The row's own configuration is KEPT when it goes away, so
%   re-ticking the producer brings it back exactly as it was, and the OTHER row is
%   never touched.
app = getApp(fig);
type = char(type); stepId = char(stepId);
step = stepById(app.reg,stepId);
if isempty(step), lines = {}; return; end
br = step.branch;
app.typeSel = wbTypeSelection('set', app.typeSel, app.reg, type, br, stepId, logical(tf));
setApp(fig,app);
lines = {sprintf('constructor: [%s] %s %s -> %s row %s', type, ternary(logical(tf),'+','-'), ...
    stepId, rowFlagFor(getApp(fig),type,br), ternary(logical(tf),'added','removed'))};
if ~tf
    kept = wbTypeSelection('steps', app.typeSel, app.reg, type, br);
    if isempty(kept)
        lines{end+1} = sprintf('  its configuration is kept - re-tick %s to get the row back', stepId);
    end
end
for i = 1:numel(lines), wbLog(fig,lines{i}); end
renderConstructor(fig);
end

function lines = tickRow(fig,type,branch,stepId,tf)
%tickRow  Tick / untick one cell of a (type,branch) row and run the cascade.
%   TICK pulls the prerequisites in, so a chain is constructible in one click;
%   UNTICK pushes the dependants out, since they could no longer run.  Both stay
%   INSIDE the row - the other branch of the same recording is a separate pipeline -
%   and both are logged, because a box the user did not click must never move
%   silently.
app = getApp(fig);
type = char(type); branch = char(branch); stepId = char(stepId);
flag = rowFlagFor(app,type,branch);
[app.typeSel, changed, animalIds] = wbTypeSelection('set', app.typeSel, app.reg, ...
    type, branch, stepId, logical(tf));
lines = {};
if tf
    lines{end+1} = sprintf('constructor: [%s %s] + %s', type, flag, stepId);
    if ~isempty(changed)
        lines{end+1} = sprintf('  auto-ticked prerequisite(s) of %s: %s', stepId, strjoin(changed,', '));
    end
    for i = 1:numel(animalIds)                     % a prerequisite that is an ANIMAL step
        app.animalSel(animalIds{i}) = true;
    end
    if ~isempty(animalIds)
        lines{end+1} = sprintf('  also ticked animal step(s) %s (%s requires them)', ...
            strjoin(animalIds,', '), stepId);
    end
else
    lines{end+1} = sprintf('constructor: [%s %s] - %s', type, flag, stepId);
    if ~isempty(changed)
        lines{end+1} = sprintf('  unticked dependant(s) of %s: %s', stepId, strjoin(changed,', '));
    end
    animalIds = animalIds(cellfun(@(id) isKey(app.animalSel,id), animalIds));
    if ~isempty(animalIds)
        % animal steps span every type and both rows, so ONE row dropping a
        % prerequisite must not switch them off for the rest - say so instead.
        lines{end+1} = sprintf('  note: animal step(s) %s stay selected, but %s %s no longer provides %s', ...
            strjoin(animalIds,', '), type, flag, stepId);
    end
end
setApp(fig,app);
for i = 1:numel(lines), wbLog(fig,lines{i}); end
renderConstructor(fig);
end

%% ---- Constructor: the derived (type, branch) rows ------------------ %%
function rows = constructorRows(app)
%constructorRows  The bottom panel's rows: one per (type, BRANCH) the user has
%   asked for.  BOTH halves of the label are data - the type is whatever token the
%   Files tab produced, and the product flag is resolved from the producing step
%   plus THIS TYPE's own settings, so a spatial protocol reads '_s_K' where a
%   temporal one reads '_t_K', in the same session.
rows = struct('type',{},'branch',{},'flag',{},'label',{},'files',{});
types = constructorTypes(app);
brs   = wbTypeSelection('branches', app.reg);
for i = 1:numel(types)
    for b = 1:numel(brs)
        if ~wbTypeSelection('rowOn', app.typeSel, app.reg, types{i}, brs{b}), continue; end
        rows(end+1) = mkRow(app,types{i},brs{b}); %#ok<AGROW>
    end
end
end
function r = mkRow(app,type,branch)
f = rowFlagFor(app,type,branch);
r = struct('type',type,'branch',branch,'flag',f, ...
           'label',sprintf('%s (%s)',type,f),'files',typeFileCount(app,type));
end

function f = rowFlagFor(app,type,branch)
%rowFlagFor  The PRODUCT FLAG a (type,branch) row's files carry - '_t_K', '_s_K',
%   '_c_K'.  Never a literal: the producing step's own outSuffix supplies the
%   product and the default stage, and a producer whose stage is a SETTING (the
%   contrast step's contrastType) is asked for this type's answer.  Row identity
%   is (type, BRANCH), so switching a project from temporal to spatial contrast
%   changes the flag shown and the files resolved, not the configuration.
f = '';
pid = wbTypeSelection('producer', app.reg, branch);
if isempty(pid), return; end
step = stepById(app.reg,pid);
if isempty(step.outSuffix), return; end
m  = wbFileModel(['x' step.outSuffix{1} '.mat']);   % '_t_K_d' -> stage t, product K
st = m.stage;
alt = settingStage(app,type,step);
if ~isempty(alt), st = alt; end
f = ['_' st '_' m.product];
end

function st = settingStage(app,type,step)
%settingStage  The stage flag a producer writes when its own SETTINGS choose it.
%   Only the contrast step does today (temporal -> _t, spatial -> _s); the rule
%   lives here once and both the row labels and the file resolution read it.
st = '';
s = wbSettingsModel('resolve', app.sm, step, char(type));
if isfield(s,'contrastType'), st = stageOfContrastType(s.contrastType); end
end
function st = stageOfContrastType(v)
st = 't';
if strcmpi(char(string(v)),'spatial'), st = 's'; end
end

function tf = rawChecked(app,type,stepId)
%rawChecked  Is this raw producer on for this type - i.e. does its row exist?
tf = wbTypeSelection('isOn', app.typeSel, type, branchOfStep(app.reg,stepId), stepId);
end
function tf = rowShown(app,type,branch,stepId)
%rowShown  What the cell DISPLAYS: its own tick, or the one it inherits.
tf = wbTypeSelection('effective', app.typeSel, app.reg, type, branch, stepId);
end
function tf = rowInherited(app,type,branch,stepId)
%rowInherited  Is this cell showing a tick that belongs to the copy-source row?
[~,tf] = wbTypeSelection('effective', app.typeSel, app.reg, type, branch, stepId);
end

function renderDerivedPanel(fig)
%renderDerivedPanel  One row per (type, product), columns = the derived steps.  A
%   cell is greyed with its reason when the row's branch does not offer that step
%   (vasomotion is contrast-only, pulsatility cardiac-only), and shown as an
%   INHERITED tick when the step is drawn once on the contrast side and copied to
%   the other branches (setRegions, segmentation - branchScope 'copy').
app = getApp(fig); c = app.c.constructor;
delete(c.derivedPanel.Children);
cols = wbTypeSelection('derivedSteps', app.reg);
rows = constructorRows(app);
types = constructorTypes(app);
if isempty(cols) || (isempty(rows) && isempty(types))
    g = uigridlayout(c.derivedPanel,[1 1],'Padding',[16 16 16 16]);
    uilabel(g,'Text','Load and label recordings on the Files tab first.', ...
        'HorizontalAlignment','center','FontAngle','italic');
    return
end

gl = uigridlayout(c.derivedPanel,[2 1],'RowHeight',{'fit','1x'},'RowSpacing',4,'Padding',[4 4 4 4]);
uilabel(gl,'Text',derivedHint(rows),'WordWrap','on','FontAngle','italic','FontSize',10);

% with no producer ticked there is nothing to append to: show the shape, disabled
preview = isempty(rows);
if preview
    for i = 1:numel(types), rows(end+1) = mkRow(app,types{i},''); end %#ok<AGROW>
end

grid = uigridlayout(gl,[1+numel(rows) 1+numel(cols)], ...
    'ColumnWidth',[{190}, repmat({88},1,numel(cols))], ...
    'RowHeight',[{56}, repmat({28},1,numel(rows))], ...
    'RowSpacing',2,'ColumnSpacing',2,'Padding',[2 2 2 2],'Scrollable','on');
h = uilabel(grid,'Text','Type (product)  \  Step','FontWeight','bold');
h.Layout.Row = 1; h.Layout.Column = 1;
for s = 1:numel(cols)
    step = stepById(app.reg,cols{s});
    hp = uigridlayout(grid,[2 1],'RowHeight',{'1x','fit'},'Padding',[1 1 1 1],'RowSpacing',1);
    hp.Layout.Row = 1; hp.Layout.Column = s+1;
    uilabel(hp,'Text',step.label,'WordWrap','on','FontSize',10,'Tooltip',headerTip(step));
    gr = uigridlayout(hp,[1 1],'Padding',[0 0 0 0]);
    uibutton(gr,'Text','settings','FontSize',9,'ButtonPushedFcn',@(~,~)selectConstructorStep(fig,step.id));
end
for r = 1:numel(rows)
    row = rows(r);
    txt = row.label;
    if preview, txt = sprintf('%s - tick a raw step above', row.type); end
    lb = uilabel(grid,'Text',sprintf('%s - %d files', txt, row.files),'FontWeight','bold', ...
        'Tooltip','one pipeline: this type''s recordings, this product');
    lb.Layout.Row = r+1; lb.Layout.Column = 1;
    if preview, lb.FontColor = [0.45 0.45 0.45]; end
    for s = 1:numel(cols)
        cb = makeRowCell(fig,grid,app,row,cols{s},preview);
        cb.Layout.Row = r+1; cb.Layout.Column = s+1;
    end
end
end

function cb = makeRowCell(fig,grid,app,row,stepId,preview)
%makeRowCell  One derived cell: live, greyed-with-a-reason, or inherited.
why = wbTypeSelection('why', app.reg, row.branch, stepId);
if preview || ~isempty(why)
    if isempty(why), why = 'tick a raw step above to give this type a product to work on'; end
    cb = uicheckbox(grid,'Text','','Value',false,'Enable','off','Tooltip',why);
    return
end
[tf,inh] = wbTypeSelection('effective', app.typeSel, app.reg, row.type, row.branch, stepId);
if inh
    src = wbTypeSelection('anchorBranch', app.reg);
    cb = uicheckbox(grid,'Text','','Value',tf,'Enable','off', ...
        'Tooltip',sprintf('%s  Tick it on the %s row.', ...
        sharedReason(stepById(app.reg,stepId)), rowFlagFor(app,row.type,src)));
    return
end
cb = uicheckbox(grid,'Text','','Value',tf,'Tooltip',cellTip(app.reg,stepId), ...
    'ValueChangedFcn',@(o,~)tickRow(fig,row.type,row.branch,stepId,o.Value));
end

function t = sharedReason(step)
%sharedReason  Why a step is ticked once for the whole recording instead of per
%   row - straight from its branchScope, which is also what the executor obeys.
switch step.branchScope
    case 'copy'
        t = sprintf(['%s is done ONCE, on the contrast product, and the result is ' ...
            'inherited by the other products of the same recording.'], step.label);
    otherwise
        t = sprintf(['%s runs over EVERY product of a recording in one call, so it ' ...
            'cannot be on for one product and off for another.'], step.label);
end
end

function t = derivedHint(rows)
if isempty(rows)
    t = ['Nothing to configure yet: tick a raw step above.  Each raw step gives its ' ...
         'type a product, and each product is its own pipeline row here.'];
else
    t = ['One row = one pipeline (a type''s recordings, one product).  A type running ' ...
         'both raw steps gets two independent rows.'];
end
end

function t = cellTip(reg,stepId)
step = stepById(reg,stepId);
t = step.label;
req = wbPrereqs('describe',step);
if ~isempty(req)
    t = [t sprintf('\nneeds %s - ticking this ticks them too, in this row', req)];
end
end

%% ---- Constructor: the per-animal steps ----------------------------- %%
function renderAnimalBox(fig,parent)
%renderAnimalBox  The steps that span the ANIMAL (registration, vessel typing):
%   ONE box, one line each - a checkbox and a 'settings' button, exactly like the
%   step columns beside it.  Their parameters open in the settings panel on the
%   right rather than inline, and WHICH BRANCH FILE of each animal's reference they
%   resolve to is reported in the selection summary, the one place status is read.
app = getApp(fig);
ids = wbTypeSelection('animalSteps', app.reg);
p  = uipanel(parent,'Title','Applied per animal','FontWeight','bold','FontSize',10);
n  = max(1,numel(ids));
gl = uigridlayout(p,[n+1 2],'ColumnWidth',{'1x','fit'}, ...
    'RowHeight',[repmat({'fit'},1,n), {'1x'}], ...
    'RowSpacing',6,'ColumnSpacing',6,'Padding',[8 8 8 8]);
for i = 1:numel(ids)
    step = stepById(app.reg,ids{i});
    uicheckbox(gl,'Text',step.label,'Value',isKey(app.animalSel,step.id), ...
        'Tooltip',animalStepTip(step), ...
        'ValueChangedFcn',@(o,~)tickAnimalStep(fig,step.id,o.Value));
    uibutton(gl,'Text','settings','FontSize',9, ...
        'Tooltip','open this step''s parameters in the settings panel', ...
        'ButtonPushedFcn',@(~,~)selectConstructorStep(fig,step.id));
end
uilabel(gl,'Text',['One run per animal, over all of its files whatever their type ' ...
    'or product.  These parameters are global.'],'WordWrap','on','FontSize',9, ...
    'FontColor',[0.4 0.4 0.4]);
end

function t = animalStepTip(step)
t = sprintf(['%s runs ONCE per animal over all of its files, whatever their type.\n' ...
    'It reads the animal''s reference recording''s %s branch.'], step.label, refBranchWord(step));
end
function w = refBranchWord(step)
switch step.refBranch
    case 'contrast', w = 'contrast (_t/_s)';
    case 'cardiac',  w = 'cardiac (_c)';
    otherwise,       w = 'any';
end
end
function s = warnLine(w)
if isempty(w), s = 'Every selected step can resolve its reference.'; else, s = ['Warning: ' strjoin(w,'   |   ')]; end
end

function tickAnimalStep(fig,stepId,tf)
%tickAnimalStep  Switch one per-animal step on/off.  It is NOT a row of the type
%   matrix: it spans the animal, so it is one flag for the whole session.
%   Ticking it pulls its prerequisites in FOR EVERY TYPE - the step runs over all
%   of an animal's files whatever their type, so each of them must produce what it
%   reads (vessel typing needs a BFI everywhere, registration a contrast).
app = getApp(fig);
stepId = char(stepId);
if ~any(strcmp(stepId, wbTypeSelection('animalSteps',app.reg))), return; end
if tf, app.animalSel(stepId) = true;
elseif isKey(app.animalSel,stepId), remove(app.animalSel,stepId);
end
setApp(fig,app);
wbLog(fig,sprintf('constructor: animal step %s %s', stepId, ternary(logical(tf),'on','off')));
if tf, autoTickPrereqs(fig,stepId); end
renderConstructor(fig);
end

function autoTickPrereqs(fig,stepId)
%autoTickPrereqs  Give every type the steps this animal step consumes.  The step
%   runs over ALL of an animal's files whatever their type or product, so each type
%   must produce what it reads - and a type with no row yet gets the raw producer
%   that creates one (wbPrereqs names the default producer of a requiresAny list).
app = getApp(fig);
step = stepById(app.reg,stepId);
if isempty(step) || isempty(wbPrereqs('all',step)), return; end
types = constructorTypes(app);
ticked = {};
for t = 1:numel(types)
    ty  = types{t};
    brs = wbTypeSelection('rows', app.typeSel, app.reg, ty);
    if isempty(brs)                                  % no product yet: make one
        pid = defaultProducerFor(app.reg, step.id);
        if ~isempty(pid)
            app.typeSel = wbTypeSelection('tick', app.typeSel, app.reg, ty, ...
                branchOfStep(app.reg,pid), pid);
            ticked = [ticked, {pid}]; %#ok<AGROW>
        end
        brs = wbTypeSelection('rows', app.typeSel, app.reg, ty);
    end
    for b = 1:numel(brs)
        have = wbTypeSelection('steps', app.typeSel, app.reg, ty, brs{b});
        need = wbPrereqs('missing', step, have);     % nothing when the row already feeds it
        for r = 1:numel(need)
            if ~wbTypeSelection('offers', app.reg, brs{b}, need{r}), continue; end
            app.typeSel = wbTypeSelection('tick', app.typeSel, app.reg, ty, brs{b}, need{r});
            ticked = [ticked, need(r)]; %#ok<AGROW>
        end
    end
end
setApp(fig,app);
if ~isempty(ticked)
    wbLog(fig,sprintf('  ticked %s for every type - %s runs over all of an animal''s files', ...
        strjoin(unique(ticked,'stable'),', '), stepId));
end
end

function id = defaultProducerFor(reg,stepId)
%defaultProducerFor  The raw producer a step's prerequisite chain ends at, so a
%   type with no row yet can be given one.  wbPrereqs already names the DEFAULT
%   producer of a requiresAny list (its first entry); this just follows the chain
%   until it reaches a step that reads the recording itself.
id = ''; seen = {}; frontier = {char(stepId)};
while ~isempty(frontier)
    cur = frontier{1}; frontier(1) = [];
    if any(strcmp(cur,seen)), continue; end
    seen{end+1} = cur; %#ok<AGROW>
    if ~isempty(branchOfStep(reg,cur)), id = cur; return; end
    frontier = [frontier, wbPrereqs('missing', stepById(reg,cur), {})]; %#ok<AGROW>
end
end

function b = branchOfStep(reg,stepId)
%branchOfStep  The row a raw producer creates ('' when the step is not one).
b = '';
step = stepById(reg,stepId);
if isempty(step), return; end
if any(strcmp(stepId, wbTypeSelection('rawSteps',reg))), b = step.branch; end
end

function plan = animalStepPlan(app)
%animalStepPlan  Per animal x reference-taking step: the pinned reference recording
%   and the file that step resolves to on it (wbRefBranch).  This is the surface
%   that makes D5 legible - the user pins a RECORDING, each step takes the branch
%   of it that it needs, and a missing branch falls back rather than failing.
%   THE STEPS ARE THOSE DECLARING A refBranch, not "the per-animal ones": vessel
%   typing paints the reference's _c, registration templates on its _t, and the
%   vascular tree - per FILE, but derived on the cardiac product - prefers _c too.
plan = struct('animal',{},'refIdentity',{},'refLabel',{},'stepId',{},'stepLabel',{}, ...
              'selected',{},'path',{},'name',{},'status',{},'note',{});
if isempty(app.files) || isempty(app.labels), return; end
animals = wbTypeModel('values', app.labels, 'animal');
ids     = refStepsOf(app.reg);
for a = 1:numel(animals)
    an    = animals{a};
    refId = animalRefOf(app,an);
    rm    = modelForIdentity(app,refId);
    cst   = contrastStageForModel(app,rm);
    ty    = typeOfIdentity(app,refId);
    for i = 1:numel(ids)
        step = stepById(app.reg,ids{i});
        r = wbRefBranch(step, rm, cst);
        plan(end+1) = struct('animal',an,'refIdentity',refId,'refLabel',shortId(refId), ...
            'stepId',step.id,'stepLabel',step.label,'selected',stepIsSelected(app,step,ty), ...
            'path',r.path,'name',r.name,'status',r.status,'note',r.note); %#ok<AGROW>
    end
end
end

function ids = refStepsOf(reg)
%refStepsOf  The steps that resolve a branch OF THE REFERENCE recording, in
%   registry order - derived from refBranch, never a name list.
ids = reshape({reg(~cellfun(@isempty,{reg.refBranch})).id},1,[]);
end

function tf = stepIsSelected(app,step,type)
%stepIsSelected  Whether a step will run on THIS reference recording: the animal
%   panel decides for a per-animal step (it spans every type), while a per-file one
%   is on only if a row of the reference's OWN type ticked it.
if strcmp(step.arity,'perAnimal'), tf = isKey(app.animalSel,step.id); return; end
tf = false;
if isempty(type), return; end
brs = wbTypeSelection('rows', app.typeSel, app.reg, type);
for b = 1:numel(brs)
    if wbTypeSelection('isOn', app.typeSel, type, brs{b}, step.id), tf = true; return; end
end
end

function m = modelForIdentity(app,identity)
%modelForIdentity  The wbFileModel of any loaded FILE of that recording ([] if none).
m = [];
if isempty(identity), return; end
for i = 1:numel(app.files)
    if strcmp(app.files(i).model.identity, identity), m = app.files(i).model; return; end
end
m = wbFileModel(identity);          % not loaded, but the identity still locates it
end

function lines = animalPlanLines(app)
%animalPlanLines  One line per animal: its reference and every animal step's file.
lines = {};
plan = animalStepPlan(app);
if isempty(plan), return; end
animals = unique({plan.animal},'stable');
for a = 1:numel(animals)
    p = plan(strcmp({plan.animal},animals{a}));
    if isempty(p(1).refIdentity)
        lines{end+1} = sprintf('%s: NO reference recording pinned', animals{a}); %#ok<AGROW>
        continue
    end
    bits = cell(1,numel(p));
    for i = 1:numel(p)
        switch p(i).status
            case 'ok',       bits{i} = sprintf('%s %s', p(i).stepId, p(i).name);
            case 'fallback', bits{i} = sprintf('%s %s (fallback)', p(i).stepId, p(i).name);
            otherwise,       bits{i} = sprintf('%s (no file)', p(i).stepId);
        end
    end
    lines{end+1} = sprintf('%s: ref = %s -> %s', animals{a}, p(1).refLabel, strjoin(bits,' | ')); %#ok<AGROW>
end
end

function w = constructorWarnings(app)
%constructorWarnings  What a SELECTED reference-taking step cannot resolve.  Warn,
%   never block: an animal may legally have no reference until the step is run.
%   A FILE THAT IS ABOUT TO BE PRODUCED IS NOT A PROBLEM.  Vessel typing reads a
%   _BFI file; complaining that it is missing while BFI is ticked for that type is
%   noise, since the run order puts BFI first.  So a missing/fallback branch is
%   only reported when the step that would produce it is NOT selected.
w = {};
plan = animalStepPlan(app);
for i = 1:numel(plan)
    p = plan(i);
    if ~p.selected, continue; end
    if strcmp(p.status,'noref')
        w{end+1} = sprintf('%s: %s is selected but the animal has no reference recording', ...
            p.animal, p.stepId); %#ok<AGROW>
    elseif any(strcmp(p.status,{'fallback','missing'})) && ~willBeProduced(app,p)
        w{end+1} = sprintf('%s: %s - %s', p.animal, p.stepId, p.note); %#ok<AGROW>
    end
end
end

function tf = willBeProduced(app,p)
%willBeProduced  Whether this run will create the file the step wants, so that
%   its absence today is not worth a warning.  The reference recording's own TYPE
%   decides, and the two failure modes ask different questions:
%     'missing'  - no readable file at all: every prerequisite of the step must be
%                  ticked for that type (vessel typing warning about an absent
%                  _BFI is noise when BFI is queued to run first);
%     'fallback' - a file exists but on the WRONG branch: the entry step that
%                  creates the wanted branch must be ticked as well, because no
%                  amount of downstream work will produce a cardiac product from a
%                  contrast one.
tf = false;
step = stepById(app.reg,p.stepId);
if isempty(step) || isempty(wbPrereqs('all',step)), return; end
ty = typeOfIdentity(app,p.refIdentity);
if isempty(ty), return; end
sel = animalStepsOn(app);
brs = wbTypeSelection('rows', app.typeSel, app.reg, ty);
for b = 1:numel(brs)
    sel = [sel, wbTypeSelection('steps', app.typeSel, app.reg, ty, brs{b}), ...
                wbTypeSelection('inherited', app.typeSel, app.reg, ty, brs{b})]; %#ok<AGROW>
end
tf  = wbPrereqs('met', step, unique(sel,'stable'));
if tf && strcmp(p.status,'fallback')
    % the wanted branch does not exist yet: only its raw producer can create it
    tf = any(strcmp(wbTypeSelection('producer', app.reg, step.refBranch), sel));
end
end

function ty = typeOfIdentity(app,identity)
ty = '';
for i = 1:numel(app.files)
    if strcmp(app.files(i).model.identity, identity), ty = app.files(i).type; return; end
end
end
function ids = animalStepsOn(app)
ids = wbTypeSelection('animalSteps', app.reg);
ids = ids(cellfun(@(id) isKey(app.animalSel,id), ids));
end

%% ---- Constructor: settings per (step, type) ------------------------ %%
function refreshConstructorSelectors(fig)
%refreshConstructorSelectors  Keep the type / step dropdowns valid for the data.
app = getApp(fig); c = app.c.constructor;
types = constructorTypes(app);
srcItems = [presetItems(), types];      % the protocols first, then your own types
if isempty(types)
    c.typeDrop.Items = {'(no types)'}; c.typeDrop.Value = '(no types)';
    c.copySrc.Items  = srcItems; c.copySrc.Value = srcItems{1};
    c.copyDst.Items  = {''}; c.copyDst.Value = '';
    app.selType = ''; setApp(fig,app);
else
    c.typeDrop.Items = types;
    if ~any(strcmp(app.selType,types)), app.selType = types{1}; end
    c.typeDrop.Value = app.selType;
    c.copySrc.Items = srcItems; c.copyDst.Items = types;
    if ~any(strcmp(c.copySrc.Value,srcItems)), c.copySrc.Value = srcItems{1}; end
    c.copyDst.Value = app.selType;
    setApp(fig,app);
end
% every step is editable here, INCLUDING the per-animal ones - their box carries a
% 'settings' button like any other column and its parameters open in this panel
cols = {app.reg.id};
app = getApp(fig);
if isempty(cols)
    c.stepDrop.Items = {''}; c.stepDrop.ItemsData = {''}; app.selCStep = '';
else
    lbls = cell(1,numel(cols));
    for i = 1:numel(cols)
        st = stepById(app.reg,cols{i});
        lbls{i} = st.label;
        if strcmp(st.arity,'perAnimal'), lbls{i} = [st.label '  (per animal)']; end
    end
    c.stepDrop.Items = lbls; c.stepDrop.ItemsData = cols;
    if ~any(strcmp(app.selCStep,cols)), app.selCStep = cols{1}; end
    c.stepDrop.Value = app.selCStep;
end
setApp(fig,app);
refreshConstructorSettings(fig);
end

function selectConstructorType(fig,type)
app = getApp(fig);
types = constructorTypes(app);
if ~any(strcmp(type,types)), return; end
app.selType = char(type); setApp(fig,app);
c = app.c.constructor;
if isgraphics(c.typeDrop), c.typeDrop.Value = app.selType; end
refreshConstructorSettings(fig);
end
function selectConstructorStep(fig,stepId)
app = getApp(fig);
if ~any(strcmp(stepId, {app.reg.id})), return; end
app.selCStep = char(stepId); setApp(fig,app);
c = app.c.constructor;
if isgraphics(c.stepDrop), c.stepDrop.Value = app.selCStep; end
refreshConstructorSettings(fig);
end

function refreshConstructorSettings(fig)
%refreshConstructorSettings  Render the settings panel for the selected step.
%   SCOPE depends on the step: a per-file step is configured per TYPE, while a
%   per-animal step spans every type by definition and so edits the global layer.
app = getApp(fig); c = app.c.constructor;
if ~isgraphics(c.paramPanel), return; end
step = [];
if ~isempty(app.selCStep), step = stepById(app.reg,app.selCStep); end
isAnimal = ~isempty(step) && strcmp(step.arity,'perAnimal');
if isempty(step) || (isempty(app.selType) && ~isAnimal)
    delete(c.paramPanel.Children);
    g = uigridlayout(c.paramPanel,[1 1],'Padding',[8 8 8 8]);
    uilabel(g,'Text','Label the files with a TYPE on the Files tab first.','WordWrap','on','FontAngle','italic');
    c.stepInfo.Text = ''; c.scopeInfo.Text = '';
    return
end
c.stepInfo.Text = stepInfoText(step);
if isAnimal
    c.scopeInfo.Text = sprintf(['%s runs once per ANIMAL, over all of its files whatever ' ...
        'their type or product, so these values are GLOBAL - there is no per-type ' ...
        'version of them.'], step.label);
    buildSettingsPanel(fig,c.paramPanel,step,'');
else
    c.scopeInfo.Text = sprintf(['These values apply to EVERY file of type "%s" (%d file(s)) - both of ' ...
        'its product rows if it has two - and to no other type.  Editing one marks this ' ...
        'type''s downstream steps stale.'], app.selType, typeFileCount(app,app.selType));
    buildSettingsPanel(fig,c.paramPanel,step,app.selType);
end
end

function v = getTypeSetting(fig,type,stepId,field)
app = getApp(fig);
s = wbSettingsModel('resolve', app.sm, stepById(app.reg,stepId), char(type));
v = getfieldOr(s,field,[]);
end

function onTypeSettingEdit(fig,type,stepId,field,rawval)
%onTypeSettingEdit  An edit scoped to ONE type: it changes that type's resolved
%   settings, its files' fingerprints, and nothing of any other type.
app = getApp(fig);
type = char(type);
step = stepById(app.reg,stepId);
if isempty(step), return; end
value = coerceValue(step,field,rawval);
if ismember(field,step.sharedKeys)
    app.sm = wbSettingsModel('setTypeShared', app.sm, type, field, value);
    seed = field;                      % shared: every step of THIS type that reads it
else
    app.sm = wbSettingsModel('setTypeStep', app.sm, type, stepId, field, value);
    seed = stepId;                     % STEP-ID keyed, never field-keyed: 'deleteOriginal'
end                                    % on BFI must not disturb splitRegions
setApp(fig,app);
recomputeBase(fig);                    % new fingerprint for this type's files
applyInvalidation(fig,seed,type);      % forward cascade, restricted to this type
refreshCells(fig,{});
% a setting can change the PRODUCT this type writes (contrastType -> _t_K | _s_K),
% and that flag is on every row label and checkbox - so redraw them, do not leave
% the panels showing the flag the type had before the edit
renderRawPanel(fig);
renderDerivedPanel(fig);
refreshConstructorSettings(fig);
refreshSummary(fig);
wbLog(fig,sprintf('edit [%s] %s.%s -> %s',type,stepId,field,val2str(value)));
end

function items = presetItems()
%presetItems  The standard protocols, tagged so they cannot collide with a type
%   token (a user is free to call a type 'NVC').
names = wbTypePresets('names');
items = cellfun(@(n) ['[protocol] ' n], names, 'UniformOutput', false);
end
function name = presetNameOf(item)
%presetNameOf  '' when the picked item is a type rather than a protocol.
name = '';
item = char(item);
tag = '[protocol] ';
if startsWith(item, tag) && wbTypePresets('exists', item(numel(tag)+1:end))
    name = item(numel(tag)+1:end);
end
end

function copyTypeConfig(fig,src,dst)
%copyTypeConfig  "Build a BP out of a BV", or out of a standard protocol: same
%   steps, same settings, then diverge.  A protocol source populates the boxes and
%   leaves the settings at the launcher defaults - it defines a pipeline, not
%   parameters.
app = getApp(fig);
src = char(src); dst = char(dst);
if isempty(dst), setConstructorStatus(fig,'Pick a type to configure.'); return; end
preset = presetNameOf(src);
if ~isempty(preset)
    applyPreset(fig,preset,dst); return
end
if isempty(src) || strcmp(src,dst)
    setConstructorStatus(fig,'Pick two different types.'); return
end
app.typeSel = wbTypeSelection('copy', app.typeSel, app.reg, src, dst);
app.sm      = wbSettingsModel('copyType', app.sm, src, dst);
setApp(fig,app);
recomputeBase(fig);
wbLog(fig,sprintf('constructor: copied configuration %s -> %s', src, dst));
setConstructorStatus(fig,sprintf('%s now matches %s.', dst, src));
renderConstructor(fig); refreshCells(fig,{});
end

function applyPreset(fig,presetName,type)
%applyPreset  Populate one type's checkboxes from a standard protocol.
%   The RAW steps go first, because they are what create the rows; every other step
%   then lands on each row that offers it - a protocol listing both entry steps
%   therefore configures both pipelines, and one listing a single metric leaves the
%   other branch alone.  The per-animal steps go to the animal panel, never into a
%   row, and each tick still runs the ordinary prerequisite cascade - so a protocol
%   can never produce a selection the Constructor would reject.
app = getApp(fig);
p = wbTypePresets('get', presetName);
if isempty(p), return; end
app.typeSel = wbTypeSelection('clear', app.typeSel, app.reg, type);
raw = wbTypeSelection('rawSteps', app.reg);
for i = 1:numel(p.steps)
    if ~any(strcmp(p.steps{i},raw)), continue; end
    b = branchOfStep(app.reg,p.steps{i});
    app.typeSel = wbTypeSelection('tick', app.typeSel, app.reg, type, b, p.steps{i});
end
brs = wbTypeSelection('rows', app.typeSel, app.reg, type);
for i = 1:numel(p.steps)
    if any(strcmp(p.steps{i},raw)), continue; end
    for b = 1:numel(brs)
        if ~wbTypeSelection('offers', app.reg, brs{b}, p.steps{i}), continue; end
        app.typeSel = wbTypeSelection('tick', app.typeSel, app.reg, type, brs{b}, p.steps{i});
    end
end
for i = 1:numel(p.animalSteps), app.animalSel(p.animalSteps{i}) = true; end
setApp(fig,app);
recomputeBase(fig);
wbLog(fig,sprintf('constructor: %s configured from the "%s" protocol (%s)', ...
    type, p.name, p.note));
setConstructorStatus(fig,sprintf('%s configured from "%s".', type, p.name));
renderConstructor(fig); refreshCells(fig,{});
end

function resetTypeConfig(fig,type)
%resetTypeConfig  Drop one type back to the launcher defaults (steps and settings).
app = getApp(fig);
type = char(type);
if isempty(type), setConstructorStatus(fig,'No type selected.'); return; end
app.typeSel = wbTypeSelection('clear', app.typeSel, app.reg, type);
app.sm      = wbSettingsModel('resetType', app.sm, type);
setApp(fig,app);
recomputeBase(fig);
wbLog(fig,sprintf('constructor: reset %s to the launcher defaults', type));
setConstructorStatus(fig,sprintf('%s reset to the launcher defaults.', type));
renderConstructor(fig); refreshCells(fig,{});
end

function uiCopyType(fig)
c = getApp(fig).c.constructor;
copyTypeConfig(fig, c.copySrc.Value, c.copyDst.Value);
end
function uiResetType(fig)
resetTypeConfig(fig, getApp(fig).selType);
end
function setConstructorStatus(fig,msg)
c = getApp(fig).c.constructor;
if isfield(c,'status') && isgraphics(c.status), c.status.Text = msg; end
end

%% ---- Constructor: the selection summary ---------------------------- %%
function lines = summaryLines(app)
%summaryLines  The sanity check before switching to Processing: what each ROW will
%   actually run, and how much of it there is.  One line per (type, product), so a
%   type driving both branches is shown as the two pipelines it really is.  An
%   INHERITED 'copy' step is named on the row that inherits it but is NOT counted
%   there - it runs once, on the row it is drawn on.
lines = {};
types = constructorTypes(app);
for i = 1:numel(types)
    ty  = types{i};
    n   = typeFileCount(app,ty);
    brs = wbTypeSelection('rows', app.typeSel, app.reg, ty);
    if isempty(brs)
        lines{end+1} = sprintf('%s: (no raw step selected) - %d file(s)', ty, n); %#ok<AGROW>
        continue
    end
    for b = 1:numel(brs)
        ids = wbTypeSelection('steps', app.typeSel, app.reg, ty, brs{b});
        inh = wbTypeSelection('inherited', app.typeSel, app.reg, ty, brs{b});
        lbl = sprintf('%s (%s)', ty, rowFlagFor(app,ty,brs{b}));
        txt = strjoin(ids,', ');
        if ~isempty(inh)
            txt = sprintf('%s [+ inherited: %s]', txt, strjoin(inh,', '));
        end
        lines{end+1} = sprintf('%s: %s (%d %s x %d files = %d cell-runs)', ...
            lbl, txt, numel(ids), plural(numel(ids),'step'), n, numel(ids)*n); %#ok<AGROW>
    end
end
aids = {};
allAnimal = wbTypeSelection('animalSteps', app.reg);
for i = 1:numel(allAnimal)
    if isKey(app.animalSel,allAnimal{i}), aids{end+1} = allAnimal{i}; end %#ok<AGROW>
end
nA = 0;
if ~isempty(app.labels), nA = numel(wbTypeModel('values', app.labels, 'animal')); end
if isempty(aids)
    lines{end+1} = sprintf('Animal steps: none selected - %d animal(s)', nA);
else
    lines{end+1} = sprintf('Animal steps: %s (%d %s x %d %s = %d runs)', ...
        strjoin(aids,', '), numel(aids), plural(numel(aids),'step'), ...
        nA, plural(nA,'animal'), numel(aids)*nA);
end
end
function w = plural(n,word)
w = word; if n~=1, w = [word 's']; end
end

function refreshSummary(fig)
%refreshSummary  The one place status is read: what every row runs, what the animal
%   steps will run over, WHICH FILE of each animal's reference each
%   reference-taking step resolves to, and what cannot be resolved yet.
app = getApp(fig); c = app.c.constructor;
if ~isfield(c,'summaryPanel') || ~isgraphics(c.summaryPanel), return; end
delete(c.summaryPanel.Children);
g = uigridlayout(c.summaryPanel,[2 1],'RowHeight',{'1x','fit'}, ...
    'RowSpacing',3,'Padding',[4 4 4 4]);
lines = summaryLines(app);
ref   = animalPlanLines(app);
if isempty(ref), ref = {'(no animals loaded)'}; end
lines = [lines, {''}, {'References:'}, ref];
uitextarea(g,'Value',lines,'Editable','off','FontName','monospaced', ...
    'Tooltip',['what each row runs, and the reference RECORDING of each animal with ' ...
               'the branch file every reference-taking step resolves to']);
w = constructorWarnings(app);
lb = uilabel(g,'Text',warnLine(w),'WordWrap','on');
if isempty(w), lb.FontColor = [0 0.45 0]; else, lb.FontColor = [0.75 0.35 0]; end
end

%% ===================== PROCESS tab ================================== %%
function buildProcessTab(fig)
%buildProcessTab  The MONITOR (spec D6): press Run, watch, open results.  Nothing
%   here is selectable - what runs was decided per TYPE on the Constructor tab -
%   so the file x step picture is ONE read-only uitable rather than a live widget
%   per cell.  200 files x 15 steps is 3000 widgets and was exactly where the old
%   matrix fell over; a table costs one component whatever the file count.
app = getApp(fig); t = app.tabs.process;
gl = uigridlayout(t,[1 2],'ColumnWidth',{'2.6x','1x'},'Padding',[6 6 6 6],'ColumnSpacing',8);

% -- LEFT: toolbar + the state table + the selected row + the log --
left = uigridlayout(gl,[4 1],'RowHeight',{'fit','1.8x','fit','1x'},'RowSpacing',6);
tb = uigridlayout(left,[1 7], ...
    'ColumnWidth',{'fit','fit','fit','1x','fit','fit','fit'},'Padding',[0 0 0 0]);
c.previewBtn = uibutton(tb,'Text','Preview order','ButtonPushedFcn',@(~,~)dryRun(fig), ...
    'Tooltip','list, in execution order, every step the configuration would run (no wrapper is called)');
c.runBtn = uibutton(tb,'Text','Run','ButtonPushedFcn',@(~,~)runChecked(fig), ...
    'BackgroundColor',[0.82 0.92 0.82],'Tooltip','run the configured steps, in dependency order');
c.stopBtn = uibutton(tb,'Text','Stop','Enable','off','ButtonPushedFcn',@(~,~)cancelRun(fig), ...
    'BackgroundColor',[1 0.86 0.86],'Tooltip','stop after the current step finishes');
c.progLabel = uilabel(tb,'Text','Ready.','FontAngle','italic');
c.presetDrop = uidropdown(tb,'Items',{'(launcher defaults)'},'Value','(launcher defaults)', ...
    'ValueChangedFcn',@(s,~)onPresetPick(fig,s.Value),'Tooltip','seed the settings bag from a saved preset');
uibutton(tb,'Text','Save preset...','ButtonPushedFcn',@(~,~)uiSavePreset(fig));
uibutton(tb,'Text','Load preset...','ButtonPushedFcn',@(~,~)uiLoadPreset(fig));

c.stateTable = uitable(left,'Data',{},'ColumnName',{'File'},'RowName',{}, ...
    'ColumnEditable',false,'CellSelectionCallback',@(~,ev)onProgressSelect(fig,ev), ...
    'Tooltip',progressLegend());
c.detail = uilabel(left,'Text','Select a row to see its file.','WordWrap','on','FontAngle','italic');
c.log = uitextarea(left,'Value',{'Ready.'},'Editable','off');

% -- RIGHT: the result links (spec D5) --
right = uigridlayout(gl,[3 1],'RowHeight',{'fit','1x','fit'},'RowSpacing',6);
uilabel(right,'Text','Results (newest first)','FontWeight','bold');
c.resultList = uilistbox(right,'Items',{'(no report images yet)'},'ItemsData',{''}, ...
    'Tooltip',['every report image the run produced, newest first - double-click to open it ' ...
               'in the desktop viewer, which stays zoomable while processing continues'], ...
    'DoubleClickedFcn',@(s,~)openArtifactViewer(s.Value));
c.openBtn = uibutton(right,'Text','Open selected', ...
    'ButtonPushedFcn',@(~,~)openArtifactViewer(getApp(fig).c.process.resultList.Value));

app.c.process = c; setApp(fig,app);
end

function t = progressLegend()
%progressLegend  What the state words in the table mean, in one tooltip.
t = ['read-only state of every (file, step): ' char(183) ' = the configuration does not ' ...
     'run it | queued = it will | running NN% = it is | done = on disk | error = it threw ' ...
     '| skipped = configured but its input cannot be produced'];
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
for i = 1:numel(srcs), items{i} = sprintf('%s   [%s]', srcs(i).label, srcs(i).animal); end
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
srcs = struct('path',{},'rpath',{},'animal',{},'label',{});
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
        srcs(end+1) = struct('path',p,'rpath',rp,'animal',r.animal,'label',[nm ex]); %#ok<AGROW>
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
uibutton(tb,'Text','Load workbench files & animals','ButtonPushedFcn',@(~,~)seedExplore(fig,true), ...
    'Tooltip','seed the explorer with the workbench''s loaded recordings (_r.mat), one explorer group per animal');
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
%seedExplore  Push the workbench's exportable RESULTS files + animals into the
%   hosted guiExplore (only when the file set actually changed, unless forced).
%   NOTE the explorer's own axis is called 'group' and is the EXPERIMENTAL group;
%   this seed fills it with the workbench's ANIMAL, which is what it carried before
%   the two axes were told apart (Phase 5 will feed it the real expGroup).
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
animals = unique({srcs.animal},'stable');
for ai = 1:numel(animals)
    mem = {srcs(strcmp({srcs.animal},animals{ai})).rpath};
    idx = find(ismember(order, mem));
    if ~isempty(idx), app.exploreAPI.createGroup(animals{ai}, idx); end
end
setExploreStatus(fig, sprintf('Seeded %d file(s) in %d animal(s) from the workbench.', ...
    numel(rpaths), numel(animals)));
end

function setExploreStatus(fig,msg)
app = getApp(fig);
if isfield(app.c,'explore') && isfield(app.c.explore,'exploreStatus') && isgraphics(app.c.explore.exploreStatus)
    app.c.explore.exploreStatus.Text = msg;
end
end

%% ===================== tab-change seeding ========================== %%
function onTabSelected(fig)
%onTabSelected  Gate the move off the Files tab, then lazily refresh the Export
%   list / seed Explore when its tab shows.
if ~guardTabSwitch(fig), return; end
app = getApp(fig);
if ~isfield(app,'tg') || ~isgraphics(app.tg), return; end
tab = app.tg.SelectedTab;
if     isequal(tab, app.tabs.export),  refreshExportFiles(fig);
elseif isequal(tab, app.tabs.explore), seedExplore(fig,false);
end
end

%% ===================== loaders ===================================== %%
function onSourceEdit(fig)
%onSourceEdit  Remember the root/glob boxes (they are session state, not just UI).
app = getApp(fig); c = app.c.files;
app.root = strtrim(c.root.Value);
app.glob = strtrim(c.glob.Value);
setApp(fig,app);
end
function onPatternEdit(fig,axis,value)
%onPatternEdit  A label regexp changed -> re-derive every label (overrides win).
setPattern(fig,axis,value);
end
function setPattern(fig,axis,value)
app = getApp(fig);
if ~isfield(app.patterns,axis), return; end
app.patterns.(axis) = strtrim(char(value));
setApp(fig,app);
syncFilesControls(fig);
if strcmp(axis,'ref'), return; end                  % the reference rule applies on Scan
rebuildWorkingSet(fig, {app.files.path}, true);
end
function setSource(fig,root,glob)
app = getApp(fig);
if nargin>=2 && ~isempty(root), app.root = char(root); end
if nargin>=3 && ~isempty(glob), app.glob = char(glob); end
setApp(fig,app); syncFilesControls(fig);
end

function uiBrowseRoot(fig)
app = getApp(fig);
d = uigetdir(defaultDir(app.root),'Pick the root folder to scan recursively');
if isequal(d,0), return; end
setSource(fig,d,'');
end
function uiScan(fig)
app = getApp(fig);
if isempty(app.root) || ~isfolder(app.root)
    alert(fig,'Set a valid root folder first (Browse... or type one).'); return
end
doScan(fig);
end
function n = doScan(fig)
%doScan  Recurse the root for the glob and REPLACE the working set.
%   With a Reference regexp this goes through getFileNamesList's reference mode,
%   which forces the matching file into column 1 - that file's recording IDENTITY
%   becomes the animal's default reference (a hand-pinned one still wins).
app = getApp(fig);
p = app.patterns;
if isempty(p.ref)
    disc = wbDiscoverFiles('folder', app.root, app.glob, p.animal, p.type, p.expGroup);
else
    disc = wbDiscoverFiles('structured', app.root, app.glob, p.animal, p.ref, p.type, p.expGroup);
end
app.autoRef = autoRefsFrom(disc);
setApp(fig,app);
n = rebuildWorkingSet(fig, gridPaths(disc), false);
end
function uiAddFiles(fig)
app = getApp(fig);
[f,p] = uigetfile({'*.mat;*.rls;*.cxd;*.avi','Recordings & products';'*.*','All files'}, ...
    'Pick files (multiselect)','MultiSelect','on',defaultDir(app.root));
if isequal(f,0), return; end
if ischar(f), f = {f}; end
addPaths(fig, cellfun(@(x)fullfile(p,x),f,'UniformOutput',false));
end
function uiAddFolder(fig)
app = getApp(fig);
d = uigetdir(defaultDir(app.root),'Pick a folder to scan recursively and ADD');
if isequal(d,0), return; end
disc = wbDiscoverFiles('folder', d, app.glob, app.patterns.animal, ...
    app.patterns.type, app.patterns.expGroup);
addPaths(fig, gridPaths(disc));
end
function n = addPaths(fig,paths)
%addPaths  Union new paths into the working set (manual is the escape hatch, so
%   it must COEXIST with a scan rather than replace it).
app = getApp(fig);
if ischar(paths), paths = {paths}; end
n = rebuildWorkingSet(fig, [reshape({app.files.path},1,[]), reshape(paths,1,[])], true);
end
function n = deletePaths(fig,paths)
%deletePaths  Drop rows from the WORKING SET.  Nothing is removed from disk.
app = getApp(fig);
if ischar(paths), paths = {paths}; end
keep = setdiff(reshape({app.files.path},1,[]), reshape(paths,1,[]), 'stable');
n = rebuildWorkingSet(fig, keep, true);
end
function uiClear(fig)
app = getApp(fig);
app.autoRef      = containers.Map('KeyType','char','ValueType','char');
app.animalRefMan = containers.Map('KeyType','char','ValueType','char');
setApp(fig,app);
rebuildWorkingSet(fig, {}, false);
end
function apiLoad(fig,mode,varargin)
%apiLoad  Programmatic loader: run wbDiscoverFiles(mode,...) and adopt the result
%   as the working set (the labels then come from the current regexp boxes).
disc = wbDiscoverFiles(mode, varargin{:});
app  = getApp(fig);
if isfield(disc,'patterns')
    f = intersect(fieldnames(app.patterns), fieldnames(disc.patterns));
    for i = 1:numel(f), app.patterns.(f{i}) = disc.patterns.(f{i}); end
end
app.autoRef      = autoRefsFrom(disc);
app.animalRefMan = containers.Map('KeyType','char','ValueType','char');
setApp(fig,app); syncFilesControls(fig);
rebuildWorkingSet(fig, gridPaths(disc), false);
end
function d = defaultDir(v)
if ~isempty(v) && isfolder(v), d = v; else, d = pwd; end
end

function paths = gridPaths(disc)
%gridPaths  The discovery grid's files as a flat list, row-major (animal order).
paths = {};
if ~isfield(disc,'fNames') || isempty(disc.fNames), return; end
g = disc.fNames.';                                   % transpose -> row-major read
paths = reshape(g(~cellfun(@isempty,g)),1,[]);
end

function m = autoRefsFrom(disc)
%autoRefsFrom  Animal -> reference recording IDENTITY, taken from column 1 of a
%   reference-mode grid (that IS getFileNamesList's answer to the ref regexp).
%   getFileNamesList puts SOMETHING in column 1 for every animal, matched or not,
%   so the match is re-checked here: an animal whose files match nothing simply
%   has no reference, which is legal.
m = containers.Map('KeyType','char','ValueType','char');
if ~isfield(disc,'referenceMode') || ~disc.referenceMode, return; end
rx = '';
if isfield(disc,'patterns') && isfield(disc.patterns,'ref'), rx = disc.patterns.ref; end
for r = 1:size(disc.models,1)
    mdl = disc.models{r,1};
    if isempty(mdl), continue; end
    if ~isempty(rx) && isempty(regexp(mdl.name, rx, 'once')), continue; end
    m(mdl.animal) = mdl.identity;                   % identity: never a branch flag
end
end

%% ===================== the curation loop ============================ %%
function n = rebuildWorkingSet(fig, paths, keepOverlay)
%rebuildWorkingSet  THE Files-tab loop.  Curated paths -> labels (regexp, then
%   hand overrides) -> per-animal references -> the animal grid -> rows + table.
%   Every curation action (scan, add, delete, label edit, reference tick, pattern
%   edit) funnels through here, so there is one definition of the working set.
app = getApp(fig);
paths = uniqueStable(cleanPathList(paths));

% ---- labels: the regexp answer, then the hand overrides on top -------------
app.labelsAuto = wbTypeModel('derive', paths, app.patterns);
app.labels     = wbTypeModel('applyOverrides', app.labelsAuto, app.overrides);
% modality is PARSED from the extension/product rather than matched, so it rides
% alongside the regexp axes with its own override map (same path->value contract)
app.labels.modality = modalityLabels(paths, app.modalityOvr);

% ---- effective per-animal reference: hand-pinned wins over the regexp ------
animalsNow = wbTypeModel('values', app.labels, 'animal');
app.animalRef = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(animalsNow)
    a = animalsNow{i};
    if isKey(app.animalRefMan,a),  app.animalRef(a) = app.animalRefMan(a);
    elseif isKey(app.autoRef,a),   app.animalRef(a) = app.autoRef(a);
    end
end
app.animalRef = pruneRefs(app.animalRef, paths);     % a ref whose file is gone is no ref

% ---- the animal grid + the deduped matrix rows -----------------------------
disc = wbDiscoverFiles('curated', paths, app.labels, app.animalRef);
disc.patterns = app.patterns;
app.files = buildFileEntries(app, paths, disc);
setApp(fig,app);
adoptDiscovery(fig, disc, keepOverlay);
n = numel(paths);
end

function labels = modalityLabels(paths, ovr)
%modalityLabels  The parsed modality of each file, with any hand override on top.
%   An override is only honoured if the EXTENSION allows it (wbFileModel owns
%   that rule), so a stale session can never resurrect an impossible pairing.
labels = cell(1,numel(paths));
for i = 1:numel(paths)
    m = wbFileModel(paths{i});
    v = m.modality;
    if isa(ovr,'containers.Map') && isKey(ovr,paths{i})
        cand = ovr(paths{i});
        if any(strcmp(cand, wbFileModel('modalities', m.ext))), v = cand; end
    end
    labels{i} = v;
end
end

function files = buildFileEntries(app, paths, disc)
%buildFileEntries  One entry per FILE: its model and its five curated fields.
%   (app.rows stays one entry per RECORDING - a recording owns several branch
%   products and the matrix must not show it twice.)
files = emptyFiles();
byPath = containers.Map('KeyType','char','ValueType','any');
for r = 1:size(disc.models,1)
    for c = 1:size(disc.models,2)
        if isempty(disc.models{r,c}), continue; end
        byPath(disc.models{r,c}.path) = disc.models{r,c};
    end
end
for i = 1:numel(paths)
    p = paths{i};
    if isKey(byPath,p), m = byPath(p); else, m = wbFileModel(p); end
    isRef = isKey(app.animalRef, m.animal) && strcmp(app.animalRef(m.animal), m.identity);
    files(end+1) = struct('model',m,'path',p,'name',m.name, ...
        'animal',m.animal,'type',m.type,'index',m.index,'expGroup',m.expGroup, ...
        'modality',m.modality,'isRef',logical(isRef)); %#ok<AGROW>
end
end

function m = pruneRefs(m, paths)
%pruneRefs  Drop any reference whose recording no longer has a file in the set.
if m.Count==0, return; end
live = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(paths)
    mm = wbFileModel(paths{i}); live(mm.identity) = true;
end
k = keys(m);
for i = 1:numel(k)
    if ~isKey(live, m(k{i})), remove(m,k{i}); end
end
end

function adoptDiscovery(fig,disc,keepOverlay)
%adoptDiscovery  Adopt a discovery grid: flatten to rows, set modality, render.
app = getApp(fig);
app.disc = disc;
[app.rows, app.animalNames, app.modelArr] = flattenDisc(disc);
if ~isempty(app.modelArr)
    mods = {app.modelArr.modality};
    app.modality = modeStr(mods);
else
    app.modality = 'LSCI';
end
app.reg = wbStepRegistry(app.modality);
if nargin<3 || ~keepOverlay
    % a fresh load starts with no session overlay
    app.checked = containers.Map('KeyType','char','ValueType','any');
    app.stale   = containers.Map('KeyType','char','ValueType','any');
end
setApp(fig,app);
recomputeBase(fig);
refreshFileTable(fig);
renderProgress(fig);
logPerFileFlips(fig);         % D6: say so when a cell stopped reading as done
renderConstructor(fig);       % types are data: a label edit adds/removes a row live
refreshStatus(fig);
end

function n = logPerFileFlips(fig)
%logPerFileFlips  Report how many cells the PER-FILE gating re-evaluated (D6/D8).
%   Existing projects WILL show cells fall back from done: unioning a recording's
%   settings across its pipelines made a half-processed recording read as finished
%   (only the cardiac product ever got a BFI, yet BFI showed done for the contrast
%   file too).  That is the bug surfacing, not a regression - so it is counted and
%   said out loud on load rather than left to be noticed.  The old answer is
%   recovered for free by asking about the same FILE with its pipeline flag
%   cleared, which is exactly the whole-recording union wbStateEngine used to do.
app = getApp(fig);
n = 0;
seen = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(app.files)
    f = app.files(i);
    if isempty(f.model.branch) || isKey(seen,f.path), continue; end
    seen(f.path) = true;
    was = statesOf(app, wholeRecordingModel(f.model), f.type);
    now_ = app.fileState(f.path);
    for k = 1:numel(app.reg)
        id = app.reg(k).id;
        if isfield(was,id) && strcmp(was.(id),'done') && ...
           isfield(now_,id) && ~strcmp(now_.(id),'done')
            n = n + 1;
        end
    end
end
if n>0
    wbLog(fig,sprintf(['%d cell(s) re-evaluated per-file: their step ran on ANOTHER ' ...
        'product of the same recording, not on this file.'], n));
end
end

function m = wholeRecordingModel(model)
%wholeRecordingModel  The same FILE asked about as the whole recording: identity,
%   extension and modality intact, pipeline flag cleared.  Clearing the flag is
%   what re-opens the union across every product (wbStateEngine>samePipeline);
%   re-parsing the bare identity would not do - it has no extension, so it would
%   come back with no modality and gate every step off.
m = model;
m.stage  = '';
m.branch = '';
end
%% ---- curation actions (table edits, assignment, references) -------- %%
function setLabel(fig,path,axis,value)
%setLabel  Hand-assign one label.  A value that equals what the regexp already
%   says clears the override instead of freezing it, so re-tuning the regexp
%   still works; anything else is remembered path->value and survives a re-scan.
app = getApp(fig);
axis = axisFieldOf(axis);
if isempty(axis) || ~isfield(app.overrides,axis), return; end
value = strtrim(char(value));
k = find(strcmp(app.labelsAuto.path, path), 1);
if isempty(value) || (~isempty(k) && strcmp(value, app.labelsAuto.(axis){k}))
    if isKey(app.overrides.(axis), path), remove(app.overrides.(axis), path); end
else
    app.overrides.(axis)(path) = value;
end
setApp(fig,app);
rebuildWorkingSet(fig, {app.files.path}, true);
end
function a = axisFieldOf(axis)
%axisFieldOf  Accept the table's column name ('group') for the axis ('expGroup').
switch char(axis)
    case {'group','expGroup'}, a = 'expGroup';
    case {'animal','type','index'}, a = char(axis);
    otherwise, a = '';
end
end
function ok = setModality(fig,path,value)
%setModality  Hand-correct a file's modality.  The EXTENSION decides what is
%   possible (a .rls cannot be a myograph video), so an impossible value is
%   refused rather than stored - wbFileModel owns that rule.
ok = false;
app = getApp(fig);
i = find(strcmp({app.files.path}, path), 1);
if isempty(i), return; end
value = strtrim(char(value));
allowed = wbFileModel('modalities', app.files(i).model.ext);
if ~any(strcmp(value,allowed)), return; end
if strcmp(value, wbFileModel(path).modality)
    if isKey(app.modalityOvr,path), remove(app.modalityOvr,path); end   % back to parsed
else
    app.modalityOvr(path) = value;
end
setApp(fig,app);
rebuildWorkingSet(fig, {app.files.path}, true);
ok = true;
end
function setReference(fig,path,tf)
%setReference  Tick/untick 'ref'.  RADIO SEMANTICS SCOPED TO THE ANIMAL: pinning
%   a recording replaces that animal's previous reference and leaves every other
%   animal alone.  What is stored is the recording IDENTITY, never a file path -
%   each step resolves the branch (_t / _c / ...) it needs at run time.
app = getApp(fig);
i = find(strcmp({app.files.path}, path), 1);
if isempty(i), return; end
m = app.files(i).model;
if tf
    app.animalRefMan(m.animal) = m.identity;
else
    if isKey(app.animalRefMan,m.animal), remove(app.animalRefMan,m.animal); end
    if isKey(app.autoRef,m.animal),      remove(app.autoRef,m.animal); end
end
setApp(fig,app);
rebuildWorkingSet(fig, {app.files.path}, true);
end
function id = animalRefOf(app,animal)
id = '';
if isKey(app.animalRef,animal), id = app.animalRef(animal); end
end
function v = labelValues(app,axis)
axis = axisFieldOf(axis);
if isempty(axis) || isempty(app.labels), v = {}; return; end
v = wbTypeModel('values', app.labels, axis);
end

function onFileTableEdit(fig,src,~)
%onFileTableEdit  Reconcile the whole table after any cell edit.
%   Reading the FULL Data, keyed by the non-editable 'file' column, instead of
%   the edited index keeps this correct whatever the table's sort order is - a
%   sorted uitable renumbers what the user sees, not its Data.  The key works
%   because file names are unique across a scanned tree; a name that is NOT
%   unique is reported by fileProblems and its rows are left alone here.
D = src.Data;
if isempty(D), return; end
col = columnIndex();
newRef = ''; unRef = {}; refused = {};
for i = 1:size(D,1)
    p = pathOfName(getApp(fig), charOf(D{i,col.file}), charOf(D{i,col.filetype}));
    if isempty(p), continue; end                            % unknown or ambiguous
    for ax = {'animal','type','index','group'}
        f = fileEntry(getApp(fig), p);                      % labels may have moved
        v = charOf(D{i,col.(ax{1})});
        if ~strcmp(v, f.(axisFieldOf(ax{1}))), setLabel(fig,p,ax{1},v); end
    end
    f = fileEntry(getApp(fig), p);
    if  logical(D{i,col.reference}) && ~f.isRef
        if isempty(animalRefOf(getApp(fig), f.animal))
            newRef = p;
        else
            refused{end+1} = sprintf('%s already has a reference', f.animal); %#ok<AGROW>
        end
    end
    if ~logical(D{i,col.reference}) &&  f.isRef, unRef{end+1} = p; end %#ok<AGROW>
end
% the reference is a RADIO scoped to the animal: pin one where there is none, or
% unpin the current one first (the other rows of that animal are greyed meanwhile)
if ~isempty(newRef)
    setReference(fig,newRef,true);
else
    for i = 1:numel(unRef), setReference(fig,unRef{i},false); end
end
if ~isempty(refused)
    setCurStatus(fig, ['refused: ' strjoin(unique(refused),'; ') ...
        ' - untick its current reference first.']);
end
refreshFileTable(fig); refreshStatus(fig);
end
function f = fileEntry(app,path)
%fileEntry  The working-set entry of one path ([] when it is not in the set).
f = [];
k = find(strcmp({app.files.path}, path), 1);
if ~isempty(k), f = app.files(k); end
end
function p = pathOfName(app,name,ext)
%pathOfName  The path behind a table row, keyed by its name + extension columns -
%   '' when the pair is not in the set or is AMBIGUOUS (two files of that name,
%   which fileProblems flags as a hard block).
p = '';
if nargin<3, ext = ''; end
k = find(strcmp({app.files.name}, [name ext]));
if isscalar(k), p = app.files(k).path; end
end
function setCurStatus(fig,msg)
c = getApp(fig).c.files;
if isfield(c,'curStatus') && isgraphics(c.curStatus), c.curStatus.Text = msg; end
end
function c = columnIndex()
%columnIndex  Field-safe name -> column number of the curation table.
names = fileTableKeys();
c = struct();
for i = 1:numel(names), c.(names{i}) = i; end
end
function s = charOf(v)
if isempty(v), s = ''; elseif ischar(v), s = v; else, s = char(string(v)); end
end

function uiAssignLabel(fig)
%uiAssignLabel  Give EVERY selected row the same value in one action.
app = getApp(fig); c = app.c.files;
sel = selectedPaths(app,c);
if isempty(sel)
    setCurStatus(fig,'Select one or more rows in the table first (click, then shift/ctrl-click).');
    return
end
v = strtrim(c.assignVal.Value);
if isempty(v), setCurStatus(fig,'Type or pick a value first.'); return; end
field = c.assignAxis.Value;
[n,refused] = quickAssign(fig, sel, field, v);
msg = sprintf('%d of %d selected row(s) -> %s "%s".', n, numel(sel), field, v);
if refused > 0
    msg = [msg sprintf('  %d refused (the file extension does not allow it).', refused)];
end
setCurStatus(fig,msg);
end
function [n,refused] = quickAssign(fig, paths, field, value)
%quickAssign  Set one field on many files at once - the bulk-curation workhorse
%   behind "Apply to selected rows" (and the programmatic entry point for it).
n = 0; refused = 0;
if ischar(paths), paths = {paths}; end
for i = 1:numel(paths)
    if strcmp(field,'modality')
        if setModality(fig,paths{i},value), n = n + 1; else, refused = refused + 1; end
    else
        setLabel(fig,paths{i},field,value);
        n = n + 1;
    end
end
refreshFileTable(fig); refreshStatus(fig);
end
function uiDeleteSelected(fig)
app = getApp(fig); c = app.c.files;
sel = selectedPaths(app,c);
if isempty(sel), c.curStatus.Text = 'Select one or more rows first.'; return; end
deletePaths(fig,sel);
c.curStatus.Text = sprintf('%d row(s) removed from the working set (files untouched on disk).', numel(sel));
end
function selectRows(fig,paths)
%selectRows  Select the table rows holding these paths (the programmatic twin of
%   clicking / shift-clicking, so bulk curation is drivable and testable).
app = getApp(fig); t = app.c.files.fileTbl;
if ~isgraphics(t), return; end
if ischar(paths), paths = {paths}; end
D = t.DisplayData; if isempty(D), D = t.Data; end
col = columnIndex(); rows = [];
for i = 1:size(D,1)
    p = pathOfName(app, charOf(D{i,col.file}), charOf(D{i,col.filetype}));
    if ~isempty(p) && any(strcmp(p,paths)), rows(end+1) = i; end %#ok<AGROW>
end
t.Selection = reshape(rows,1,[]);   % row-selection tables demand a 1-by-N vector
refreshAssignItems(fig);
end

function sel = selectedPaths(app,c)
%selectedPaths  The selected rows as PATHS, read through the DISPLAYED data so a
%   sorted table selects what the user actually clicked.
sel = {};
t = c.fileTbl;
if ~isgraphics(t), return; end
idx = [];
if isprop(t,'DisplaySelection'), idx = t.DisplaySelection; end
if isempty(idx) && isprop(t,'Selection'), idx = t.Selection; end
if isempty(idx), return; end
if strcmp(t.SelectionType,'row')
    idx = unique(idx(:));            % row selection: a 1-by-N list of ROW indices
else
    idx = unique(idx(:,1));          % cell selection: N-by-2 [row col] pairs
end
D = t.DisplayData;
if isempty(D), D = t.Data; end
col = columnIndex();
for i = 1:numel(idx)
    if idx(i) > size(D,1), continue; end
    p = pathOfName(app, charOf(D{idx(i),col.file}), charOf(D{idx(i),col.filetype}));
    if ~isempty(p), sel{end+1} = p; end %#ok<AGROW>
end
end

%% ---- the gate: nothing proceeds on a half-curated file set --------- %%
function [problems, warnings] = fileProblems(app)
%fileProblems  What still stands between this file set and processing.
%   BLOCKING problems: file names that are not unique across the scanned tree
%   (the workbench keys rows, and the pipeline keys products, by name), and any
%   file whose animal / type / index / group is still sitting in its no-match
%   bucket.  Every field must be assigned before the next stage - use the regexp
%   boxes, or select the rows and Quick-assign them.
%   WARNINGS do not block: an animal with no reference recording is legal until
%   a per-animal step (registration, vessel typing) is actually selected.
problems = {}; warnings = {};
if isempty(app.files), problems = {'No files loaded - scan a folder or add files.'}; return; end

names = {app.files.name};
[u,~,ic] = unique(names);
dup = u(accumarray(ic,1) > 1);
if ~isempty(dup)
    problems{end+1} = sprintf(['%d file name(s) appear more than once (%s). ' ...
        'The workbench and the pipeline identify recordings BY NAME, so a scanned ' ...
        'tree must not repeat one - rename them or scan a narrower root.'], ...
        numel(dup), strjoin(shortList(dup),', '));
end

%   'index' is the one axis whose no-match bucket ('1') is a real answer - with no
%   index regexp every file simply IS index 1 - so only an empty one is missing.
ax = {'animal','type','index','expGroup'};
nm = {'animal','type','index','group'};
for i = 1:numel(ax)
    bucket = wbTypeModel('default', ax{i});
    if strcmp(ax{i},'index'), bucket = char(0);  end        % never matches: '1' is fine
    bad = strcmp({app.files.(ax{i})}, bucket) | cellfun(@isempty,{app.files.(ax{i})});
    if any(bad)
        problems{end+1} = sprintf('%d file(s) have no %s assigned (%s).', ...
            nnz(bad), nm{i}, strjoin(shortList({app.files(bad).name}),', ')); %#ok<AGROW>
    end
end

noRef = {};
animalsNow = wbTypeModel('values', app.labels, 'animal');
for i = 1:numel(animalsNow)
    if ~isKey(app.animalRef, animalsNow{i}), noRef{end+1} = animalsNow{i}; end %#ok<AGROW>
end
if ~isempty(noRef)
    warnings{end+1} = sprintf(['%d animal(s) have no reference recording (%s) - ' ...
        'legal, but registration and vessel typing need one.'], ...
        numel(noRef), strjoin(shortList(noRef),', '));
end
end
function s = shortList(c)
%shortList  At most four names, so a status line stays a status line.
s = reshape(c,1,[]);
if numel(s) > 4, s = [s(1:4), {sprintf('+%d more',numel(s)-4)}]; end
end
function tf = filesValid(app)
tf = isempty(fileProblems(app));
end
function refreshProblems(fig)
%refreshProblems  The banner above the status line: the gate, in words.
app = getApp(fig); c = app.c.files;
if ~isfield(c,'problems') || ~isgraphics(c.problems), return; end
[p,w] = fileProblems(app);
syncTabLock(fig);                       % no files = nothing to move on to, either
if isempty(app.files)
    c.problems.Text = ''; return
end
if isempty(p)
    c.problems.Text = ['Ready for processing.' warnText(w)];
    c.problems.FontColor = [0 0.45 0];
else
    c.problems.Text = ['Not ready - fix before processing:  ' strjoin(p,'   |   ') warnText(w)];
    c.problems.FontColor = [0.75 0.2 0];
end
syncTabLock(fig);
end

function syncTabLock(fig)
%syncTabLock  Grey the other tabs (and disable Next) while the file set is
%   incomplete.  A uitab has no Enable property, so the lock is its title colour
%   plus guardTabSwitch, which bounces a click back to the Files tab.
app = getApp(fig);
ok = filesValid(app);
locked = [0.62 0.62 0.62]; open = [0 0 0];
names = fieldnames(app.tabs);
for i = 1:numel(names)
    tb = app.tabs.(names{i});
    if strcmp(names{i},'files') || ~isgraphics(tb) || ~isprop(tb,'ForegroundColor'), continue; end
    tb.ForegroundColor = ternary(ok, open, locked);
end
if isfield(app.c,'files') && isfield(app.c.files,'nextBtn') && isgraphics(app.c.files.nextBtn)
    app.c.files.nextBtn.Enable = ternary(ok,'on','off');
end
end

function goToConstructor(fig)
%goToConstructor  The Files tab's "Next": move on once the set is complete.
if ~guardTabSwitchTo(fig, getApp(fig).tabs.constructor), return; end
end
function s = warnText(w)
if isempty(w), s = ''; else, s = ['   (' strjoin(w,'; ') ')']; end
end
function tf = guardTabSwitch(fig)
%guardTabSwitch  Refuse to leave the Files tab while the file set is incomplete.
%   The author's rule: no half-curated set reaches the next stage.
tf = true;
app = getApp(fig);
if ~isfield(app,'tg') || ~isgraphics(app.tg), return; end
if isequal(app.tg.SelectedTab, app.tabs.files), return; end
if filesValid(app), return; end
app.tg.SelectedTab = app.tabs.files;
tf = false;
alertIncomplete(fig,app);
end

function tf = guardTabSwitchTo(fig,tab)
%guardTabSwitchTo  Move to a tab, or refuse and say what is still missing.
tf = false;
app = getApp(fig);
if ~isfield(app,'tg') || ~isgraphics(app.tg), return; end
if ~filesValid(app), alertIncomplete(fig,app); return; end
app.tg.SelectedTab = tab;
tf = true;
end

function alertIncomplete(fig,app)
alert(fig, sprintf(['Finish the Files tab first:\n\n%s\n\nEvery file needs an animal, ' ...
    'type, index and group.  Select rows and use Quick assign to set them in bulk.'], ...
    strjoin(fileProblems(app), newline)));
end

function p = cleanPathList(p)
if isempty(p), p = {}; return; end
if ischar(p), p = {p}; end
p = reshape(p,1,[]);
p = p(~cellfun(@isempty,p));
p = cellfun(@char,p,'UniformOutput',false);
end
function u = uniqueStable(c)
if isempty(c), u = {}; return; end
[~,ia] = unique(c,'stable'); u = c(sort(ia));
end


function [rows, animalNames, modelArr] = flattenDisc(disc)
%flattenDisc  Discovery grid -> flat rows in animal-major, reference-first order.
%   ONE ROW PER RECORDING IDENTITY.  A recording can own several branch products
%   with the SAME identity (e.g. a _t_K and a _c_K from one .rls), so a glob that
%   matches more than one branch (*_K_d.mat) would otherwise emit duplicate rows.
%   The first model seen anchors the row; every step still resolves its own branch
%   file on disk (wbExecutor), and the disk state unions all of the recording's
%   _s files (wbStateEngine), so the anchor branch does not change what runs.  In
%   the normal single-entry-glob workflow no identity ever repeats, so this is a
%   no-op there.
rows = emptyRows(); animalNames = {}; modelArr = wbFileModel('x.rls'); modelArr(1) = [];
if isempty(disc.models), return; end
seen = containers.Map('KeyType','char','ValueType','logical');
[nr,nc] = size(disc.models);
for r = 1:nr
    if ~isempty(disc.animals) && numel(disc.animals)>=r, aname = disc.animals(r).name; else, aname = sprintf('animal%d',r); end
    animalNames{end+1} = aname; %#ok<AGROW>
    for cc = 1:nc
        m = disc.models{r,cc};
        if isempty(m), continue; end
        if isKey(seen, m.identity), continue; end     % dedup: one row per recording
        seen(m.identity) = true;
        row = struct('model',m,'identity',m.identity,'animal',aname,'animalIdx',r, ...
            'rowInAnimal',cc,'isRef',logical(disc.referenceMode && cc==1), ...
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

function D = fileTableData(app)
%fileTableData  The curation table's Data: ONE ROW PER FILE, in working-set order.
n = numel(app.files);
D = cell(n,numel(fileTableColumns()));
for i = 1:n
    f = app.files(i);
    [~,bare,ext] = fileparts(f.path);
    D(i,:) = {f.isRef, bare, f.animal, f.type, f.index, f.expGroup, ...
              ext, f.modality, stageLabel(f.model)};
end
end
function refreshFileTable(fig)
app = getApp(fig); c = app.c.files;
c.fileTbl.Data = fileTableData(app);
styleReferenceColumn(fig);
refreshAssignItems(fig);
refreshProblems(fig);
end

function styleReferenceColumn(fig)
%styleReferenceColumn  Grey the reference ticks of an animal that already has one.
%   The reference is a RADIO scoped to the animal: once it is pinned, the other
%   rows of that animal are not available until it is unpinned, and the table says
%   so (a uitable cannot disable a single cell, so this is colour + a refusal in
%   the edit callback - onFileTableEdit).
app = getApp(fig); t = app.c.files.fileTbl;
if ~isgraphics(t), return; end
removeStyle(t);
if isempty(app.files), return; end
col = columnIndex();
rows = [];
for i = 1:numel(app.files)
    f = app.files(i);
    if ~f.isRef && ~isempty(animalRefOf(app,f.animal)), rows(end+1) = i; end %#ok<AGROW>
end
if isempty(rows), return; end
sty = uistyle('FontColor',[0.65 0.65 0.65],'BackgroundColor',[0.94 0.94 0.94]);
addStyle(t, sty, 'cell', [rows(:), repmat(col.reference,numel(rows),1)]);
end
function refreshAssignItems(fig)
%refreshAssignItems  Seed the assign combo from the values ACTUALLY IN USE on the
%   chosen field - the vocabulary is discovered, never a built-in list, and the
%   box stays editable so a brand-new value can simply be typed.  Modality is the
%   one exception: its vocabulary IS fixed (wbFileModel owns it) and is offered
%   narrowed to what every selected row's extension can be.
app = getApp(fig); c = app.c.files;
if ~isfield(c,'assignVal') || ~isgraphics(c.assignVal), return; end
if strcmp(c.assignAxis.Value,'modality')
    v = allowedModalitiesFor(app, selectedPaths(app,c));
    c.assignVal.Editable = 'off';
else
    v = labelValues(app, c.assignAxis.Value);
    c.assignVal.Editable = 'on';
end
if isempty(v), v = {''}; end
cur = c.assignVal.Value;
c.assignVal.Items = v;
if any(strcmp(cur,v)), c.assignVal.Value = cur; else, c.assignVal.Value = v{1}; end
end
function v = allowedModalitiesFor(app, paths)
%allowedModalitiesFor  What EVERY one of these files could be (the intersection).
v = wbFileModel('modalities');
for i = 1:numel(paths)
    f = fileEntry(app, paths{i});
    if isempty(f), continue; end
    v = intersect(v, wbFileModel('modalities', f.model.ext), 'stable');
end
end
function s = stageLabel(m)
if m.isRaw, s = 'raw'; elseif isempty(m.flags), s = m.product; else, s = [strjoin(m.flags,'_') '_' m.product]; end
end
function syncFilesControls(fig)
%syncFilesControls  Push root/glob/patterns back into their boxes (session load).
app = getApp(fig);
if ~isfield(app,'c') || ~isfield(app.c,'files'), return; end
c = app.c.files;
if isgraphics(c.root), c.root.Value = app.root; end
if isgraphics(c.glob), c.glob.Value = app.glob; end
f = fieldnames(c.pat);
for i = 1:numel(f)
    if isfield(app.patterns,f{i}) && isgraphics(c.pat.(f{i}))
        c.pat.(f{i}).Value = app.patterns.(f{i});
    end
end
end
function refreshStatus(fig)
%refreshStatus  N files / A animals / T types / G groups / U untyped / R refless.
app = getApp(fig);
if isempty(app.files)
    txt = 'No files loaded.';
else
    animalsNow = wbTypeModel('values', app.labels, 'animal');
    nNoRef = 0;
    for i = 1:numel(animalsNow)
        if ~isKey(app.animalRef, animalsNow{i}), nNoRef = nNoRef + 1; end
    end
    nUntyped = sum(strcmp({app.files.type}, wbTypeModel('default','type')));
    txt = sprintf(['%d files - %d animal(s) - %d type(s) - %d group(s) - ' ...
        '%d untyped - %d animal(s) without a reference.  Modality %s.'], ...
        numel(app.files), numel(animalsNow), ...
        numel(wbTypeModel('values', app.labels, 'type')), ...
        numel(wbTypeModel('values', app.labels, 'expGroup')), ...
        nUntyped, nNoRef, app.modality);
end
app.c.files.status.Text = txt;
end

%% ===================== state derivation ============================ %%
function recomputeBase(fig)
%recomputeBase  The disk picture (+ the current-settings fingerprint), at the
%   three scopes the workbench asks about:
%     .base        per RECORDING - the legacy row view the look-ahead uses;
%     .fileState   per FILE       - what the monitor shows (D6: each file is
%                                   gated on its OWN pipeline, so a recording
%                                   whose cardiac product alone got a BFI no
%                                   longer reads as fully done);
%     .branchState per (recording, PIPELINE) - what the run expansion asks, since
%                                   a step is only skippable when it is done on
%                                   every product it would touch (D8).
%   All three come from the same wbStateEngine call shape; only the model handed
%   to it differs, which is exactly where the per-file rule lives.
app = getApp(fig);
app.base        = containers.Map('KeyType','char','ValueType','any');
app.fileState   = containers.Map('KeyType','char','ValueType','any');
app.branchState = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(app.rows)
    m = app.rows(i).model;
    app.base(m.identity) = statesOf(app, m, modelType(m));
end
for i = 1:numel(app.files)
    f = app.files(i);
    app.fileState(f.path) = statesOf(app, f.model, f.type);
end
brs = wbTypeSelection('branches', app.reg);
for i = 1:numel(app.rows)
    id = app.rows(i).identity;
    ty = typeOfIdentity(app, id);
    for b = 1:numel(brs)
        bm = branchModelOf(id, brs{b});
        if isempty(bm), continue; end                % that pipeline has no file yet
        app.branchState([id '||' brs{b}]) = statesOf(app, bm, ty);
    end
end
setApp(fig,app);
end

function bs = statesOf(app, model, type)
%statesOf  wbStateEngine's answer for one model, folded into a stepId->state struct.
if nargin<3, type = modelType(model); end
cs = curSettingsFor(app, model, type);
st = wbStateEngine(model, app.reg, cs);
bs = struct();
for k = 1:numel(st), bs.(st(k).id) = st(k).state; end
end

function m = branchModelOf(identity, branch)
%branchModelOf  ANY data file of one recording sitting on one pipeline ([] if the
%   pipeline has produced nothing yet).  Which file does not matter: wbStateEngine
%   unions the settings ALONG the pipeline, so '_t_K_d' and '_t_BFI_d' give the
%   same answer - and asking for a fixed one would break the moment a step deleted
%   its original (runBFI's deleteOriginal).
m = [];
d = dir([identity '*_d.mat']);
for i = 1:numel(d)
    cm = wbFileModel(fullfile(d(i).folder, d(i).name));
    if strcmp(cm.identity, identity) && strcmp(cm.branch, branch), m = cm; return; end
end
end

function cs = curSettingsFor(app,model,type)
%curSettingsFor  The settings a recording would run with - resolved for ITS TYPE,
%   so the staleness fingerprint follows the type's configuration.  A bare parsed
%   model carries no type, so the caller passes the one the Files tab assigned.
if nargin<3 || isempty(type), ty = modelType(model); else, ty = char(type); end
cs = struct();
for k = 1:numel(app.reg)
    cs.(app.reg(k).id) = wbSettingsModel('resolve', app.sm, app.reg(k), model, ty);
end
end
function ty = modelType(model)
ty = '';
if isstruct(model) && isfield(model,'type'), ty = char(model.type); end
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
req  = wbPrereqs('all', step);
if isempty(step) || isempty(req), return; end
have = {};
for i = 1:numel(req)
    reqState = resolveCellState(app,identity,req{i});             % recurse over the DAG
    if any(strcmp(reqState,{'done','stale','checked'})), have{end+1} = req{i}; end %#ok<AGROW>
end
if ~wbPrereqs('met', step, have), return; end
s = 'ready';
end

%% ===================== the progress table (read-only) =============== %%
function rows = progressRows(app)
%progressRows  The monitor's rows: ONE PER FILE (spec D4/D6), ordered the way the
%   run goes - animal, reference first, then the file's position in its animal -
%   so watching the table reads top to bottom.  A recording with two product files
%   contributes two rows, which is the whole point: each is gated on its OWN
%   pipeline (D8) and each carries only the steps its branch runs.
rows = emptyProgressRows();
for i = 1:numel(app.files)
    f = app.files(i);
    k = find(strcmp({app.rows.identity}, f.model.identity),1);
    if isempty(k), ai = inf; ria = inf; isRef = f.isRef; else
        ai = app.rows(k).animalIdx; ria = app.rows(k).rowInAnimal; isRef = app.rows(k).isRef;
    end
    rows(end+1) = struct('path',f.path,'model',f.model,'identity',f.model.identity, ...
        'label',f.name,'animal',f.animal,'type',f.type,'expGroup',f.expGroup, ...
        'branch',f.model.branch,'isRef',logical(isRef),'animalIdx',ai,'rowInAnimal',ria); %#ok<AGROW>
end
if isempty(rows), return; end
[~,ord] = sortrows([[rows.animalIdx]', [rows.rowInAnimal]', (1:numel(rows))']);
rows = rows(ord);
end

function renderProgress(fig)
%renderProgress  (Re)build the state table from scratch: rows, columns, contents.
app = getApp(fig); c = app.c.process;
if ~isfield(c,'stateTable') || ~isgraphics(c.stateTable), return; end
app.pRows  = progressRows(app);
app.pRowOf = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(app.pRows), app.pRowOf(app.pRows(i).path) = i; end
setApp(fig,app);

c.stateTable.ColumnName = [{'File'}, {app.reg.label}];
c.stateTable.ColumnWidth = [{230}, repmat({94},1,numel(app.reg))];
refreshCells(fig,{});
refreshProgressDetail(fig);
end

function refreshCells(fig,identities)
%refreshCells  Repaint the state table.  It writes table DATA, never components:
%   the whole point of the read-only monitor is that a 200-file project costs one
%   uitable instead of 3000 live widgets.  'identities' narrows the repaint to the
%   rows of those recordings ({} = all), which is what keeps a per-cell progress
%   tick cheap during a run.
app = getApp(fig); c = app.c.process;
if ~isfield(c,'stateTable') || ~isgraphics(c.stateTable), return; end
if isempty(app.pRows)
    c.stateTable.Data = {};
    return
end
D = c.stateTable.Data;
nCol = numel(app.reg)+1;
if ~iscell(D) || size(D,1)~=numel(app.pRows) || size(D,2)~=nCol
    D = cell(numel(app.pRows), nCol);
    identities = {};                                   % a fresh grid: fill it all
end
if nargin<2 || isempty(identities)
    which = 1:numel(app.pRows);
else
    which = find(ismember({app.pRows.identity}, identities));
end
for i = which
    r = app.pRows(i);
    D{i,1} = [r.label ternary(r.isRef,'  (ref)','')];
    planned = plannedStepsFor(app, r);
    for s = 1:numel(app.reg)
        D{i,s+1} = progressCellText(app, r, app.reg(s), planned);
    end
end
c.stateTable.Data = D;
end

function ids = plannedStepsFor(app, r)
%plannedStepsFor  Which steps THIS FILE's configuration runs on it, in registry
%   order.  A file sitting on one of its type's product rows gets that row's
%   selection - what it runs plus what it inherits from the anchor row - so the
%   cardiac file of a two-pipeline recording never shows the contrast row's steps.
%   A file with no stage flag (the raw recording, the usual case) stands for the
%   whole recording and gets the union.  The per-animal steps are added for the
%   animal, not the row: they span every type by definition.
ids = animalStepsOn(app);
if isempty(r.type), ids = orderIds(app.reg, ids); return; end
brs = wbTypeSelection('rows', app.typeSel, app.reg, r.type);
if ~isempty(r.branch) && any(strcmp(r.branch, brs)), brs = {r.branch}; end
for b = 1:numel(brs)
    ids = [ids, wbTypeSelection('steps',     app.typeSel, app.reg, r.type, brs{b}), ...
                wbTypeSelection('inherited', app.typeSel, app.reg, r.type, brs{b})]; %#ok<AGROW>
end
ids = orderIds(app.reg, ids);
end

function out = orderIds(reg, ids)
%orderIds  Unique step ids in registry (dependency) order.
out = {reg(ismember({reg.id}, ids)).id};
end

function txt = progressCellText(app, r, step, planned)
%progressCellText  ONE cell of the monitor (spec D4).  The words are the states a
%   watcher cares about, in this order of authority: what the run is doing now,
%   then what the disk says, then what the configuration intends.
key = cellKey(r.identity, step.id);
if isKey(app.runState,key) && stepTouchesFile(step, r.branch)
    switch app.runState(key)
        case 'running'
            pct = 0;
            if isKey(app.cellPct,key), pct = app.cellPct(key); end
            txt = sprintf('running %d%%', round(pct)); return
        case 'done',  txt = 'done';  return
        case 'error', txt = 'error'; return
    end
end
st = fileStateOf(app, r.path, step.id);
if strcmp(st,'done'),                        txt = 'done'; return; end
if ~any(strcmp(step.id, planned)),           txt = char(183); return; end   % not configured
if any(strcmp(st,{'ready','stale'})),        txt = 'queued'; return; end
% configured, but blocked: queued only if everything it needs is itself queued
have = {};
req = wbPrereqs('all', step);
for i = 1:numel(req)
    rs = fileStateOf(app, r.path, req{i});
    if strcmp(rs,'done') || any(strcmp(req{i}, planned)), have{end+1} = req{i}; end %#ok<AGROW>
end
if wbPrereqs('met', step, have), txt = 'queued'; else, txt = 'skipped'; end
end

function tf = stepTouchesFile(step, fileBranch)
%stepTouchesFile  Whether a step's work lands on this file.  A file with no stage
%   flag is the whole recording; a branch-agnostic step covers every product of it
%   ('all' in one call, 'copy' by propagation); anything else is branch-local.
tf = isempty(fileBranch) || isempty(step.branch) || strcmp(step.branch, fileBranch);
end

function s = fileStateOf(app, path, stepId)
%fileStateOf  The PER-FILE disk state of one (file,step) - D6's whole point.
s = 'unavailable';
if ~isKey(app.fileState,path), return; end
bs = app.fileState(path);
if isfield(bs,stepId), s = bs.(stepId); end
if strcmp(s,'done') && isKey(app.stale, cellKey(fileIdentityOf(app,path), stepId))
    s = 'stale';                                   % an in-session edit pushed it stale
end
end
function id = fileIdentityOf(app, path)
id = path;
if isKey(app.pRowOf,path), id = app.pRows(app.pRowOf(path)).identity; end
end

function onProgressSelect(fig, ev)
%onProgressSelect  Remember which row the user clicked and describe it below.
if isempty(ev) || ~isfield(struct(ev),'Indices') || isempty(ev.Indices), return; end
app = getApp(fig);
app.pSel = ev.Indices(1,1);
setApp(fig,app);
refreshProgressDetail(fig);
end

function refreshProgressDetail(fig)
%refreshProgressDetail  The selected file's identity card: where it is, how it is
%   labelled, and what went wrong if something did.
app = getApp(fig); c = app.c.process;
if ~isfield(c,'detail') || ~isgraphics(c.detail), return; end
i = 0;
if isfield(app,'pSel'), i = app.pSel; end
if i<1 || i>numel(app.pRows)
    c.detail.Text = 'Select a row to see its file.';
    return
end
r = app.pRows(i);
bits = {sprintf('%s   animal %s | type %s | group %s', r.path, ...
    dashIfEmpty(r.animal), dashIfEmpty(r.type), dashIfEmpty(r.expGroup))};
blocked = firstSkipped(app, r);
if ~isempty(blocked)
    bits{end+1} = sprintf('%s skipped - %s', blocked, unavailableReason(app,r.identity,blocked));
end
err = lastErrorFor(app, r.identity);
if ~isempty(err), bits{end+1} = ['last error: ' err]; end
c.detail.Text = strjoin(bits, '   |   ');
end

function id = firstSkipped(app, r)
%firstSkipped  The first configured step this file cannot reach ('' if none) - the
%   one worth explaining, since everything after it is skipped for the same reason.
id = '';
planned = plannedStepsFor(app, r);
for k = 1:numel(app.reg)
    if strcmp(progressCellText(app, r, app.reg(k), planned),'skipped')
        id = app.reg(k).id; return
    end
end
end
function s = dashIfEmpty(v), s = char(v); if isempty(s), s = '-'; end, end

function msg = lastErrorFor(app, identity)
%lastErrorFor  The most recent error text recorded against this recording.
msg = '';
k = keys(app.cellMsg);
for i = 1:numel(k)
    p = strsplit(k{i},'||');
    if numel(p)>=2 && strcmp(p{1},identity) && ~isempty(app.cellMsg(k{i}))
        msg = sprintf('%s: %s', p{2}, app.cellMsg(k{i}));
    end
end
end

function s = unavailableReason(app,identity,stepId)
%unavailableReason  A short line for a step that cannot run (modality / prereqs).
s = 'not available yet';
step = stepById(app.reg,stepId);
if isKey(app.base,identity) && ~any(strcmp(app.modality,step.modalities))
    s = ['not for ' app.modality];
elseif ~isempty(wbPrereqs('all',step))
    s = ['needs: ' wbPrereqs('describe',step)];
end
end

%% ===================== check / bulk-select ========================= %%
% The per-cell queue predates the type model and no longer drives the run - the
% Constructor's (type,branch) selection does (buildRunOrder).  What survives here
% is the programmatic surface the API still exposes and the look-ahead it shares
% with the monitor; Phase 6 retires the map itself along with the old matrix.
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
function checkAnimal(fig,animalName,tf)
app = getApp(fig);
for i = 1:numel(app.rows)
    if strcmp(app.rows(i).animal,animalName)
        for s = 1:numel(app.reg), setCheckedQuiet(app,app.rows(i).identity,app.reg(s).id,tf); end
    end
end
setApp(fig,app); refreshCells(fig,{});
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
%selectStep  Which step the PROGRAMMATIC settings calls act on.  It no longer
%   renders anything: settings are a per-TYPE question and live on the Constructor
%   tab (spec D4), while the Process tab is a read-only monitor (D6).
app = getApp(fig);
if ~any(strcmp(stepId,{app.reg.id})), return; end
app.selStep = stepId; setApp(fig,app);
end
function s = stepInfoText(step)
bits = {};
if ~isempty(step.gatingField), bits{end+1} = ['gates on settings.' step.gatingField];
else, bits{end+1} = 'done by output (no settings field)'; end
if strcmp(step.arity,'perAnimal'), bits{end+1} = 'per-animal (reference in column 1)'; end
if ~isequal(step.interactive,false), bits{end+1} = 'interactive'; end
if step.needsRaw, bits{end+1} = 'also needs the raw recording'; end
if ~isempty(wbPrereqs('all',step)), bits{end+1} = ['requires ' wbPrereqs('describe',step)]; end
s = strjoin(bits,'  |  ');
end

function buildSettingsPanel(fig,panel,step,type)
%buildSettingsPanel  Generalised buildParamEditor: render step.settingGroups.
%   SCOPE.  An empty type edits the GLOBAL layer (the Process tab, and the
%   per-animal steps, which span types by definition); a type edits that type's
%   own layer and no other's - wbSettingsModel resolves and stores both.
if nargin<4, type = ''; end
delete(panel.Children);
groups = step.settingGroups;
if isempty(groups)
    gl = uigridlayout(panel,[1 1],'Padding',[8 8 8 8]);
    uilabel(gl,'Text','This step has no tunable parameters (interactive or automatic).', ...
        'WordWrap','on','FontAngle','italic');
    return
end
renderSections(fig,panel,step,type,true);
end

function renderSections(fig,parent,step,type,scroll)
%renderSections  BASIC first, ADVANCED below - the only structure of a settings
%   panel.  Both halves behave identically; the split is purely about what a user
%   has to read before they can run something.
%   TWO BOXES, AND NO BOXES INSIDE THEM (author, 2026-07-28).  A settings group is
%   a bold heading with its fields under it, separated by space rather than by yet
%   another border: nested frames made the panel unreadable at a glance.
app = getApp(fig);
s = wbSettingsModel('resolve', app.sm, step, type);   % step(+type)-level, no file
[basic, advanced] = splitBasicAdvanced(step);
sections = {'Basic', basic; 'Advanced', advanced};
sections = sections(~cellfun(@isempty, sections(:,2)), :);
stack = uigridlayout(parent,[size(sections,1) 1],'RowHeight',repmat({'fit'},1,size(sections,1)), ...
    'RowSpacing',10,'Padding',[4 4 4 4],'Scrollable',ternary(scroll,'on','off'));
for si = 1:size(sections,1)
    sec = uipanel(stack,'Title',[sections{si,1} ' settings'],'FontWeight','bold', ...
        'ForegroundColor',ternary(si==1,[0 0 0],[0.42 0.42 0.42]), ...
        'Tooltip',sectionTip(sections{si,1}));
    g = sections{si,2};
    inner = uigridlayout(sec,[2*size(g,1) 1],'RowHeight',repmat({'fit'},1,2*size(g,1)), ...
        'RowSpacing',4,'Padding',[8 6 8 8]);
    for gi = 1:size(g,1)                            % heading, then its fields
        uilabel(inner,'Text',g{gi,1},'FontWeight','bold','FontColor',[0.25 0.25 0.25]);
        fieldGrid(fig,inner,step,g{gi,2},s,type);
    end
end
end

function [basic, advanced] = splitBasicAdvanced(step)
%splitBasicAdvanced  Deal a step's setting GROUPS into the two display sections.
%   BASIC holds the fields a protocol is actually written around (the registry's
%   basicFields); ADVANCED holds the rest, because most users never touch them.  A
%   step that names no basic fields shows everything as Basic rather than hiding
%   all of it under Advanced.
groups = step.settingGroups;
bf = {};
if isfield(step,'basicFields'), bf = step.basicFields; end
if isempty(bf), basic = groups; advanced = {}; return; end
basic = {}; advanced = {};
for gi = 1:size(groups,1)
    flds = groups{gi,2};
    isB  = ismember(flds, bf);
    if any(isB),  basic(end+1,:)    = {groups{gi,1}, flds(isB)};  end %#ok<AGROW>
    if any(~isB), advanced(end+1,:) = {groups{gi,1}, flds(~isB)}; end %#ok<AGROW>
end
end

function t = sectionTip(name)
if strcmp(name,'Basic')
    t = 'the settings a protocol is normally written around';
else
    t = 'rarely changed - the launcher defaults suit most recordings';
end
end

function fieldGrid(fig,parent,step,flds,s,type)
%fieldGrid  label + control per field, laid out two columns and indented under the
%   group heading (the indent replaces the frame the group used to be drawn in).
gg = uigridlayout(parent,[numel(flds) 2],'ColumnWidth',{'fit','1x'}, ...
    'RowHeight',repmat({'fit'},1,numel(flds)),'RowSpacing',3,'Padding',[12 0 0 2]);
for fi = 1:numel(flds)
    f = flds{fi};
    lbl = uilabel(gg,'Text',f);
    if isfield(step.tips,f), lbl.Tooltip = step.tips.(f); end
    makeParamControl(fig,gg,step,f,getfieldOr(s,f,[]),type);
end
end

function makeParamControl(fig,parent,step,field,val,type)
%makeParamControl  One control sized to the field's type (enum/logical/cell/num).
if nargin<6, type = ''; end
cb = @(o) onSettingEditScoped(fig,type,step.id,field,o.Value);
if isfield(step.enums,field)
    items = step.enums.(field);
    v = char(string(val)); if ~any(strcmp(v,items)), v = items{1}; end
    uidropdown(parent,'Items',items,'Value',v,'ValueChangedFcn',@(o,~)cb(o));
elseif islogical(defaultType(step,field))
    uicheckbox(parent,'Text','','Value',logical(firstTrue(val)),'ValueChangedFcn',@(o,~)cb(o));
elseif iscell(val)
    uieditfield(parent,'text','Value',cellToStr(val),'ValueChangedFcn',@(o,~)cb(o), ...
        'Tooltip','comma-separated list');
else
    uieditfield(parent,'text','Value',val2str(val),'ValueChangedFcn',@(o,~)cb(o));
end
end

function onSettingEditScoped(fig,type,stepId,field,rawval)
%onSettingEditScoped  Route an edit to the global layer or to one type's layer.
if isempty(type)
    onSettingEdit(fig,stepId,field,rawval);
else
    onTypeSettingEdit(fig,type,stepId,field,rawval);
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
wbLog(fig,sprintf('edit %s.%s -> %s',stepId,field,val2str(value)));
end

function applyInvalidation(fig,seed,type)
%applyInvalidation  Mark every DONE cell in the forward set of the edit stale.
%   With a type, only that type's recordings are touched - which is the whole
%   point of keying the settings by type: a BP edit must leave BV's done cells
%   alone.  The seed is a STEP ID (or a shared-key name), never a bare field, so
%   'deleteOriginal' edited on BFI cannot invalidate splitRegions.
app = getApp(fig);
if isempty(app.modelArr), return; end
models = app.modelArr;
if nargin>=3 && ~isempty(type)
    models = models(strcmp({models.type}, char(type)));
    if isempty(models), return; end
end
cells = wbInvalidate(app.reg, seed, models);
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
renderConstructor(fig);            % a preset carries the per-type layers too
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
    refreshCells(fig,{});
    renderConstructor(fig);        % every per-type layer went with the bag
    wbLog(fig,'reset to launcher defaults (global AND per-type settings)');
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
%saveSessionTo  Write the durable session (spec §5).  Everything the workbench
%   knows goes in, including the RESOLVED per-file labels and the completion
%   record - guiExport and guiExplore read this file and never the workbench.
app = getApp(fig);
saveSessionInto(app, wbSession('empty'), pth);
app.sessionPath = pth; setApp(fig,app);   % later state changes autosave to here
wbLog(fig,['saved session ' pth]);
end

function files = sessionFileList(app)
%sessionFileList  The curated set with its labels RESOLVED, in table order.
%   THE EXPERIMENTAL GROUP IS CARRIED HERE AND NOWHERE ELSE IN THE RUN: processing
%   ignores it by design (spec §2), and Export/Explore are its only consumers - so
%   it has to survive into the session or it is lost.
files = struct('path',{},'name',{},'animal',{},'type',{},'index',{}, ...
    'expGroup',{},'modality',{},'identity',{},'branch',{},'isRef',{},'use',{});
for i = 1:numel(app.files)
    f = app.files(i);
    files(end+1) = struct('path',f.path,'name',f.name,'animal',f.animal, ...
        'type',f.type,'index',f.index,'expGroup',f.expGroup,'modality',f.modality, ...
        'identity',f.model.identity,'branch',f.model.branch, ...
        'isRef',logical(f.isRef),'use',true); %#ok<AGROW>
end
end

function recordCompletion(fig,identity,stepId,state,msg)
%recordCompletion  One line of the session's per-FILE, per-STEP completion record.
%   Keyed by PATH because that is what the read-only tools hold; a step that
%   covered several branch products of a recording writes one line per file it
%   touched, which is the only reading that survives D8.
if ~isvalid(fig), return; end
app = getApp(fig);
step = stepById(app.reg,stepId);
if isempty(step), return; end
stamp = datestr(now,'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
for i = 1:numel(app.pRows)
    r = app.pRows(i);
    if ~strcmp(r.identity,identity) || ~stepTouchesFile(step, r.branch), continue; end
    app.completed([r.path '||' stepId]) = struct('state',char(state), ...
        'when',stamp,'fingerprint',settingsFingerprint(app,step,r), ...
        'message',char(msg));
end
setApp(fig,app);
end

function h = settingsFingerprint(app,step,r)
%settingsFingerprint  A short, stable hash of the settings a step ran with, so a
%   resumed session can tell "already done" from "done with OTHER parameters".
%   Field order is normalised and the transport hooks are dropped, so the same
%   configuration always hashes the same way whatever order it was built in.
try
    s = wbSettingsModel('resolve', app.sm, step, r.model, r.type);
    f = fieldnames(s);
    bits = {};
    for i = 1:numel(f)
        v = s.(f{i});
        if isa(v,'function_handle'), continue; end     % transport hooks are not settings
        bits{end+1} = [f{i} '=' val2str(v)]; %#ok<AGROW>
    end
    h = shortHash(strjoin(sort(bits),';'));
catch
    h = '';                                            % unresolvable = no fingerprint
end
end

function h = shortHash(txt)
%shortHash  A plain 32-bit polynomial hash of a char row - no toolbox, no Java.
d = double(char(txt));
acc = 0;
for i = 1:numel(d), acc = mod(acc*31 + d(i), 4294967296); end
h = sprintf('%08X', acc);
end

function autosaveSession(fig)
%autosaveSession  Keep the durable session current WITHOUT being asked (spec §5).
%   A run that is killed, cancelled or crashes must still leave something
%   resumable behind, so the session is rewritten on every state change rather
%   than only when the user remembers to save.  It goes wherever the session was
%   last saved or loaded from, or beside the scanned root when there is no such
%   file yet; with neither, there is nowhere to put it and nothing is written.
if ~isvalid(fig), return; end
app = getApp(fig);
pth = app.sessionPath;
if isempty(pth)
    if isempty(app.root) || ~isfolder(app.root), return; end
    pth = fullfile(app.root,'workbench_session.mat');
end
try
    session = wbSession('empty');
    saveSessionInto(app, session, pth);
    app.sessionPath = pth; setApp(fig,app);
catch ME
    wbLog(fig,['session autosave failed: ' ME.message]);
end
end

function saveSessionInto(app, session, pth)
%saveSessionInto  Fill a blank session from the app struct and write it.  ONE
%   definition of what a session contains, used by both the explicit Save and the
%   autosave, so the two can never drift apart.
session.root          = app.root;
session.glob          = app.glob;
session.patterns      = app.patterns;
session.paths         = {app.files.path};
session.files         = sessionFileList(app);
session.completed     = app.completed;
session.overrides     = app.overrides;
session.modalityOvr   = app.modalityOvr;
session.animalRef     = app.animalRef;
session.animalRefMan  = app.animalRefMan;
session.fNames        = app.disc.fNames;
session.referenceMode = app.disc.referenceMode;
session.animalNames   = app.animalNames;
session.modality      = app.modality;
session.rowOrder      = [];
session.bag           = app.sm.bag;
session.stepOverrides = app.sm.stepOverrides;
session.fileOverrides = app.sm.fileOverrides;
session.typeBag       = app.sm.typeBag;
session.typeOverrides = app.sm.typeOverrides;
session.typeSel       = app.typeSel;
session.animalSel     = app.animalSel;
session.checked       = app.checked;
session.staleOverlay  = app.stale;
session.presetRef     = app.presetRef;
wbSession('save', pth, session);
end
function loadSessionFrom(fig,pth)
%loadSessionFrom  Restore a session WITHOUT a re-scan: the curated path list, its
%   labels (regexps + hand overrides + modality corrections), the per-animal
%   references and the settings/overlay all come back from the sidecar.
session = wbSession('load', pth);
app = getApp(fig);
% ---- curation state ---------------------------------------------------------
app.root         = session.root;
if ~isempty(session.glob), app.glob = session.glob; end
app.patterns     = session.patterns;
app.overrides    = session.overrides;
app.modalityOvr  = session.modalityOvr;
app.animalRefMan = session.animalRefMan;
app.autoRef      = session.animalRef;              % the effective refs, minus the
k = keys(app.animalRefMan);                        % hand ones, are the auto ones
for i = 1:numel(k)
    if isKey(app.autoRef,k{i}), remove(app.autoRef,k{i}); end
end
% ---- settings model + overlay ----------------------------------------------
app.modality = session.modality;
app.sm = wbSettingsModel('new', session.bag);
app.sm.stepOverrides = session.stepOverrides;
foK = keys(session.fileOverrides); foV = values(session.fileOverrides);
for i = 1:numel(foK), app.sm.fileOverrides(foK{i}) = foV{i}; end
% ---- the Constructor's output (absent from a pre-schema-3 sidecar) ---------
tbK = keys(session.typeBag); tbV = values(session.typeBag);
for i = 1:numel(tbK), app.sm.typeBag(tbK{i}) = tbV{i}; end
toK = keys(session.typeOverrides); toV = values(session.typeOverrides);
for i = 1:numel(toK), app.sm.typeOverrides(toK{i}) = toV{i}; end
% the row set is rebuilt THROUGH wbTypeSelection, so a sidecar written before the
% (type,branch) rows existed is upgraded to them instead of being read as empty
app.typeSel   = wbTypeSelection('fromCells', keys(session.typeSel), app.reg);
app.animalSel = session.animalSel;
app.checked   = session.checked;
app.stale     = session.staleOverlay;
app.completed = session.completed;                 % what already ran, per file/step
app.presetRef = session.presetRef;
app.files     = emptyFiles();
app.sessionPath = pth;                             % keep autosaving where it came from
setApp(fig,app);
syncFilesControls(fig);

% the stored path list is authoritative; older sidecars only carry the grid
paths = session.paths;
if isempty(paths), paths = gridPaths(struct('fNames',{session.fNames})); end
rebuildWorkingSet(fig, paths, true);
refreshPresetDrop(fig);
wbLog(fig,sprintf('loaded session %s (schema %g, %d file(s), %d completed cell(s))', ...
    pth, session.schema, numel(session.paths), session.completed.Count));
end

%% ===================== the run expansion =========================== %%
function entries = buildRunOrder(fig)
%buildRunOrder  THE expansion: a per-TYPE configuration becomes an ordered list of
%   (file, step) work items.  It asks three questions and derives nothing itself.
%
%   1. WHICH ROWS DOES THIS FILE'S TYPE HAVE?  wbTypeSelection('rows') - a
%      (type,branch) row exists exactly while its raw producer is ticked.  Nothing
%      here enumerates branches, step ids or type names.
%   2. WHAT DOES EACH ROW RUN?  wbTypeSelection('steps') - the raw producer first,
%      then the derived steps ticked on that row, already in registry (dependency)
%      order.  An INHERITED recording-level step ('copy' / 'all') is deliberately
%      absent from every row but its anchor, which is what queues setRegions,
%      segmentation and BFI ONCE for a two-pipeline recording instead of once per
%      pipeline - re-expanding them would split every branch twice and re-open the
%      ROI editor.
%   3. WHICH BRANCH PRODUCTS DOES A STEP ACTUALLY TOUCH?  Not asked here at all -
%      that is branchScope, resolved from disk inside wbExecutor>buildFNames
%      (spec §6).  This function names files and steps; the executor names files'
%      products.
%
%   The two per-ANIMAL steps are not part of any type: they span the animal's
%   files whatever their type OR experimental group, so they are emitted once per
%   animal, and the reference FILE each of them resolves to comes from
%   animalStepPlan (wbRefBranch), handed to the executor through ctx.refFile.  The
%   experimental group plays no part in execution at all - it rides in the session
%   for Export and Explore.
%
%   ORDER IS STEP-MAJOR: registry step -> animal -> reference-first -> file.  This
%   is the LAUNCHER'S order (author, 2026-07-29), and the reason is not cosmetic: a
%   launcher cell hands ONE step the whole fNames list, so every recording reaches
%   a given level before anything moves past it.  Steps that read across a file set
%   depend on that - registration templates the animal's files onto one reference,
%   vessel typing paints one reference and inherits to the rest, split regions
%   crops every product - and running a file all the way down before starting the
%   next one would reach them while the others were still raw.  Within a step the
%   animal's own order is kept (reference first), so an interactive step still
%   walks an animal together.  A step already done on every product it would touch
%   is left out, so a re-run resumes rather than repeats (spec §5).
app = getApp(fig); reg = app.reg;
entries = emptyEntries();
seenPA  = containers.Map('KeyType','char','ValueType','logical');
seen    = containers.Map('KeyType','char','ValueType','logical');
aIds    = animalStepsOn(app);
for i = 1:numel(app.rows)
    r  = app.rows(i);
    ty = typeOfIdentity(app, r.identity);

    % ---- the per-animal steps: once per animal, at its reference -------------
    for a = 1:numel(aIds)
        step = stepById(reg, aIds{a});
        if isempty(step), continue; end
        k = [num2str(r.animalIdx) '||' step.id];
        if isKey(seenPA,k), continue; end
        seenPA(k) = true;
        entries(end+1) = mkEntry(r.animalIdx, 1, stepIndexOf(reg,step.id), step, ...
            animalRefIdentity(app,r.animalIdx), r.animal, true, ty); %#ok<AGROW>
    end

    % ---- the type's rows, and what each of them runs on this file ------------
    if isempty(ty), continue; end
    brs = wbTypeSelection('rows', app.typeSel, reg, ty);
    for b = 1:numel(brs)
        ids = wbTypeSelection('steps', app.typeSel, reg, ty, brs{b});
        for j = 1:numel(ids)
            step = stepById(reg, ids{j});
            if isempty(step) || strcmp(step.arity,'perAnimal'), continue; end
            k = cellKey(r.identity, ids{j});
            if isKey(seen,k), continue; end             % one work item per (file,step)
            if stepAlreadyDone(app, r.identity, ty, brs{b}, step), continue; end
            seen(k) = true;
            entries(end+1) = mkEntry(r.animalIdx, r.rowInAnimal, stepIndexOf(reg,step.id), ...
                step, r.identity, r.animal, r.isRef, ty); %#ok<AGROW>
        end
    end
end
if isempty(entries), return; end
K = [[entries.stepIdx]', [entries.animalIdx]', [entries.rowInAnimal]'];
[~,ord] = sortrows(K);                       % step-major: the launcher's order
entries = entries(ord);
end

function e = emptyEntries()
e = struct('animalIdx',{},'rowInAnimal',{},'stepIdx',{},'stepId',{}, ...
    'stepLabel',{},'identity',{},'label',{},'animal',{},'isRef',{},'arity',{},'type',{});
end
function e = mkEntry(ai,ria,si,step,identity,animal,isRef,type)
lbl = identity; [~,nm] = fileparts(identity); if ~isempty(nm), lbl = nm; end
e = struct('animalIdx',ai,'rowInAnimal',ria,'stepIdx',si,'stepId',step.id, ...
    'stepLabel',step.label,'identity',identity,'label',lbl,'animal',animal, ...
    'isRef',isRef,'arity',step.arity,'type',char(type));
end
function k = stepIndexOf(reg,stepId)
k = find(strcmp({reg.id},stepId),1);
if isempty(k), k = numel(reg)+1; end
end
function id = animalRefIdentity(app,animalIdx)
id = '';
for i = 1:numel(app.rows)
    if app.rows(i).animalIdx==animalIdx, id = app.rows(i).identity; return; end
end
end

function tf = stepAlreadyDone(app, identity, type, branch, step)
%stepAlreadyDone  May the expansion leave this step out because the work is there?
%   Only when it is done on EVERY product it would touch.  A row-local ('one')
%   step touches its own row's pipeline, so its own pipeline decides; a
%   recording-level one covers the lot in a single call, so every active row must
%   already carry it - which is the D8 case again: a recording whose cardiac
%   product alone got a BFI must still get one on the contrast side.
brs = {branch};
if ~strcmp(step.branchScope,'one')
    brs = wbTypeSelection('rows', app.typeSel, app.reg, type);
end
tf = false;
for b = 1:numel(brs)
    key = [identity '||' brs{b}];
    if ~isKey(app.branchState,key), return; end          % nothing produced there yet
    bs = app.branchState(key);
    if ~isfield(bs,step.id) || ~strcmp(bs.(step.id),'done'), return; end
    if isKey(app.stale, cellKey(identity,step.id)), return; end   % an edit re-opened it
end
tf = true;
end

function lines = dryRun(fig)
entries = buildRunOrder(fig);
lines = {sprintf('=== DRY RUN: %d step(s) would execute (no wrapper is called) ===',numel(entries))};
for i = 1:numel(entries)
    e = entries(i);
    if strcmp(e.arity,'perAnimal'), who = ['ANIMAL ' e.animal]; else, who = e.label; end
    lines{end+1} = sprintf('  %2d. [%s] %-26s :: %s', i, e.animal, who, e.stepLabel); %#ok<AGROW>
end
if numel(entries)==0
    lines{end+1} = '  (nothing to do - configure the types on the Constructor tab)';
end
setLog(fig,lines);
end

%% ===================== real Run (executor) ========================= %%
function runChecked(fig)
%runChecked  Run what the per-type configuration says, through wbExecutor.
app = getApp(fig);
if isfield(app,'running') && app.running, return; end       % already running
entries = buildRunOrder(fig);
if isempty(entries)
    setLog(fig,{'Nothing to do - configure the types on the Constructor tab.'});
    return
end
if ~confirmRun(fig, entries), return; end
% a fresh run starts with a clean transient overlay
app.runState = containers.Map('KeyType','char','ValueType','any');
app.cellMsg  = containers.Map('KeyType','char','ValueType','any');
app.cellPct  = containers.Map('KeyType','char','ValueType','any');
app.running  = true;  app.cancel = false;
setApp(fig,app);
setRunningUI(fig,true);
autosaveSession(fig);                    % a run that never returns still leaves a session
cleaner = onCleanup(@() finishRun(fig)); % restore UI + fold results even on error
ctx = buildExecContext(fig);
wbExecutor(entries, ctx);
end

function ok = confirmRun(fig, entries)
%confirmRun  The pre-run summary.  A type-level tick can expand into hundreds of
%   file-steps, and the number is the only honest warning about how long this will
%   take - so it is shown before anything runs, not discovered afterwards.
ok = true;
n  = numel(entries);
tys = unique({entries.type},'stable');
tys = tys(~cellfun(@isempty,tys));
sts = unique({entries.stepId},'stable');
msg = sprintf('About to run %s across %s = %d file-step(s).  Continue?', ...
    plural(numel(sts),'step'), plural(numel(tys),'type'), n);
wbLog(fig,msg);
if ~isvalid(fig) || strcmp(fig.Visible,'off'), return; end   % headless: never block
sel = uiconfirm(fig, msg, 'Run', 'Options',{'Run','Cancel'}, ...
    'DefaultOption',1, 'CancelOption',2, 'Icon','question');
ok = strcmp(sel,'Run');
if ~ok, wbLog(fig,'Run cancelled at the summary.'); end
end

function ctx = buildExecContext(fig)
%buildExecContext  The callback bundle wbExecutor drives the figure through.
ctx = struct();
ctx.reg           = getApp(fig).reg;
ctx.modelOf       = @(id) modelByIdentity(fig,id);
ctx.animalModels   = @(gi) animalModelArr(fig,gi);
% the run must use the SAME layers the staleness fingerprint compares against
% (curSettingsFor), or a cell would read stale the moment it finished
ctx.resolve       = @(step,mdl) wbSettingsModel('resolve', getApp(fig).sm, step, mdl, modelType(mdl));
ctx.contrastStage = @(mdl) contrastStageOf(fig,mdl);
% which FILE of an animal's pinned reference recording a reference-taking step
% templates on: the Constructor already resolves it (animalStepPlan > wbRefBranch)
% and shows it in the summary, so the executor consumes that answer rather than
% re-deriving one that could disagree with what the user was shown
ctx.refFile       = @(ai,sid) animalRefFile(fig,ai,sid);
ctx.setState      = @(id,sid,state,msg) execSetState(fig,id,sid,state,msg);
ctx.progress      = @(id,sid,f,l) execProgress(fig,id,sid,f,l);
ctx.log           = @(msg) wbLog(fig,msg);
ctx.isCancelled   = @() execIsCancelled(fig);
ctx.modalGuard    = @(fcn) wbModalGuard(fig,fcn);
ctx.afterDone     = @(id,sid,mdl) execAfterDone(fig,id,sid,mdl);
end

function p = animalRefFile(fig,animalIdx,stepId)
%animalRefFile  The resolved reference FILE for one (animal, step), '' if none.
p = '';
app = getApp(fig);
an = '';
for i = 1:numel(app.rows)
    if app.rows(i).animalIdx==animalIdx, an = app.rows(i).animal; break; end
end
if isempty(an), return; end
plan = animalStepPlan(app);
k = find(strcmp({plan.animal},an) & strcmp({plan.stepId},stepId),1);
if ~isempty(k), p = plan(k).path; end
end

function execSetState(fig,identity,stepId,state,msg)
%execSetState  Set a cell's transient run overlay ('' clears it), repaint the
%   table DATA (never a per-cell component - that is what the read-only monitor
%   buys) and fold the change into the durable session, so a crash or a Stop still
%   leaves a resumable record behind.
if ~isvalid(fig), return; end
if nargin<5, msg = ''; end
app = getApp(fig);
key = cellKey(identity,stepId);
if isempty(state)
    if isKey(app.runState,key), remove(app.runState,key); end
    if isKey(app.cellMsg,key),  remove(app.cellMsg,key);  end
    if isKey(app.cellPct,key),  remove(app.cellPct,key);  end
else
    app.runState(key) = state;
    if ~isempty(msg), app.cellMsg(key) = msg; end
end
setApp(fig,app);
if any(strcmp(state,{'done','error'})), recordCompletion(fig,identity,stepId,state,msg); end
refreshCells(fig,{identity});
refreshProgressDetail(fig);
drawnow limitrate;
end

function execProgress(fig,identity,stepId,frac,label)
%execProgress  Update the running cell's percent + the global progress label.
%   The drawnow yields are what keep the window alive through a long wrapper.
if ~isvalid(fig), return; end
if nargin<5 || ~ischar(label), label = ''; end
app = getApp(fig);
pct = max(0,min(100,round(100*frac)));
key = cellKey(identity,stepId);
app.cellPct(key) = pct;
setApp(fig,app);
refreshCells(fig,{identity});
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
%   downstream done cells to stale, drop it from the queue, and add whatever
%   report images it wrote to the result list.
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
addResults(fig,identity,stepId,model);
autosaveSession(fig);
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
renderProgress(fig);
refreshResultList(fig);
autosaveSession(fig);            % a Stop or a crash still leaves a resumable session
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

function ms = animalModelArr(fig,animalIdx)
%animalModelArr  An animal's models ordered reference-first (rowInAnimal order).
app = getApp(fig);
sel = app.rows([app.rows.animalIdx]==animalIdx);
if isempty(sel), ms = []; return; end
[~,ord] = sort([sel.rowInAnimal]);
sel = sel(ord);
ms = [sel.model];
end

function st = contrastStageOf(fig,model)
%contrastStageOf  't' or 's' - which contrast flag this project's files carry.
st = contrastStageForModel(getApp(fig),model);
end

function st = contrastStageForModel(app,model)
%contrastStageForModel  The contrast flag of one recording, per ITS TYPE's setting.
%   The same rule the Constructor's row labels use (rowFlagFor > settingStage), so
%   a spatial protocol reads '_s' everywhere - rows, files and reference alike.
st = 't';
pid = wbTypeSelection('producer', app.reg, 'contrast');
if isempty(pid), return; end
s = wbSettingsModel('resolve', app.sm, stepById(app.reg,pid), model, modelType(model));
if isfield(s,'contrastType'), st = stageOfContrastType(s.contrastType); end
end

%% ===================== results (artifact links) ==================== %%
function addResults(fig,identity,stepId,model)
%addResults  Fold a finished step's report images into the result list (spec D5).
%   A LIST, NOT THUMBNAILS.  The old panel drew a live uiimage per report of the
%   cell you happened to click; across a 200-file run that is both unaffordable
%   and useless - you want the last thing that finished, and you want it in a real
%   viewer.  So the run accumulates links, newest first, and a double-click hands
%   the file to the desktop (openArtifactViewer), which keeps it zoomable while
%   the next step is already running.
if ~isvalid(fig), return; end
app = getApp(fig);
step = stepById(app.reg,stepId);
if nargin<4 || isempty(model), model = modelByIdentity(fig,identity); end
if isempty(model) || isempty(step), return; end
files = wbArtifacts(model,step);
have  = {};
if ~isempty(app.results), have = {app.results.path}; end
added = false;
stamp = datetime('now');
for i = 1:numel(files)
    if any(strcmp(files{i},have)), continue; end
    app.results(end+1) = struct('path',files{i}, ...
        'label',sprintf('%s - %s - %s', shortId(identity), step.label, shortName(files{i})), ...
        'when',stamp);
    added = true;
end
if ~added, return; end
setApp(fig,app);
refreshResultList(fig);
end

function refreshResultList(fig)
%refreshResultList  Repaint the link list, newest first.
if ~isvalid(fig), return; end
app = getApp(fig); c = app.c.process;
if ~isfield(c,'resultList') || ~isgraphics(c.resultList), return; end
if isempty(app.results)
    c.resultList.Items = {'(no report images yet)'};
    c.resultList.ItemsData = {''};
    return
end
[~,ord] = sort([app.results.when],'descend');
r = app.results(ord);
c.resultList.Items     = {r.label};
c.resultList.ItemsData = {r.path};
c.resultList.Value     = r(1).path;
end

function txt = progressCellOf(fig,path,stepId)
%progressCellOf  The monitor's word for one (file,step) - the headless view of
%   the table, so a test can assert what a watcher would read.
txt = '';
app = getApp(fig);
if ~isKey(app.pRowOf,path), return; end
r = app.pRows(app.pRowOf(path));
step = stepById(app.reg,stepId);
if isempty(step), return; end
txt = progressCellText(app, r, step, plannedStepsFor(app,r));
end

function ids = plannedStepsOf(fig,path)
%plannedStepsOf  Which steps the configuration runs on one file (headless view).
ids = {};
app = getApp(fig);
if ~isKey(app.pRowOf,path), return; end
ids = plannedStepsFor(app, app.pRows(app.pRowOf(path)));
end

function L = resultLinks(fig)
%resultLinks  The result list as {label, path} - the programmatic view of D5.
app = getApp(fig);
L = cell(0,2);
if isempty(app.results), return; end
[~,ord] = sort([app.results.when],'descend');
r = app.results(ord);
L = [{r.label}', {r.path}'];
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
