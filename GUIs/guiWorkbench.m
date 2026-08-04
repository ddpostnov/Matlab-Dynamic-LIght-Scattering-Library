%guiWorkbench - Processing Workbench: configure per TYPE, run, watch.
%
%   A programmatic uifigure that turns the launcher workflow into three questions:
%   WHICH files (Files tab), WHAT each recording TYPE runs (Constructor tab), and
%   then GO (Process tab).  The window is a thin controller: every decision is
%   delegated to the headless brain (wbStepRegistry, wbDiscoverFiles, wbFileModel,
%   wbStateEngine, wbPrereqs, wbTypeModel, wbTypeSelection, wbSettingsModel,
%   wbRefBranch, wbRunRange, wbInvalidate), to wbExecutor / wbModalGuard /
%   wbUiGuard / wbArtifacts / makeReportPdf / wbPool, and to wbSession.
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
%   whole set at the same level.  A step already done on every product it would
%   touch is left out, so a re-run RESUMES rather than repeats - and that now
%   includes the two per-animal steps, which are done once every configured product
%   of every recording of the animal carries them (animalStepAlreadyDone).  Naming
%   a From column is how you deliberately run one again.  Progress streams back
%   through the hook seam into READ-ONLY state tables with cooperative Stop and
%   continue-on-error; Preview order lists the same plan without calling anything.
%   Report images the run writes collect in a link list that opens them in the
%   desktop viewer.
%
%   A REPORT IS A BY-PRODUCT (spec D9/D10).  With 'Create PDF reports' on, the
%   images a COLUMN produced across the whole run are appended into ONE PDF - a
%   page per image, each fitted to its own image (makeReportPdf) - the moment that
%   column's last entry lands, so a 60-file step is one document to page through
%   instead of 60 links; a Stop or an error still leaves the columns that finished.
%   The list tells images and PDFs apart with a kind filter, and hands either to
%   the desktop exactly the same way.  None of it is part of the run: it cannot
%   change what a wrapper writes, cannot mark a step in error, and is a window
%   preference the session does not carry.
%   THE DOCUMENT IS MADE BETWEEN THE STEPS, NEVER INSIDE ONE.  A wrapper writes
%   its pages and stops there, whoever ran it; a launcher assembles them in a cell
%   of its own after the step, and here the assembly is per COLUMN because a column
%   spans several batched calls.  Both hosts call the same Core/Reporting builder.
%   A second switch,
%   'Delete the images once the PDF is written', is off by default: every entry in
%   this list is a path resolved from disk (wbArtifacts), so deleting the images
%   empties the list of them, and the PDF becomes the only copy of a page that
%   exists.  When it is on, the deletion and the removal of the matching links
%   happen in the same step (flushPdfColumn > dropSourceImages).
%
%   FROM / TO IS A RANGE OVER THAT LIST, NOT A SECOND SELECTION (spec D1/D2).  The
%   two dropdowns beside Run slice the expansion by registry column, so a protocol
%   can be taken to segmentation, tuned, and carried on without unticking anything.
%   'Last valid' RESUMES - it starts at the first configured column that is not
%   finished and leaves out whatever is already on disk, which is the behaviour that
%   was there before - while naming a column FORCES: the whole range runs again,
%   done or not, so segmentation can be re-run with new settings.  Both dropdowns
%   are built from the expansion's own output (wbRunRange), so nothing here lists a
%   step or branches on one.
%
%   THE SESSION IS THE DELIVERABLE - AND THE ONLY COUPLING.  wbSession is written
%   on every state change - a scan, a label, a Constructor tick, a settings edit,
%   each finished cell, the end of a run, and the exit itself - so a crash, a Stop
%   or a MATLAB restart always leaves a resumable project behind.  It always goes to
%   THE LAST SESSION, one fixed file in the gitignored 'workbench-sessions' folder
%   beside the code (sessionDir / lastSessionPath), and additionally to wherever a
%   session was last saved or loaded from, so choosing a location never stops the
%   automatic one.  It is also the whole
%   hand-off to the two READ-ONLY tools: neither is a tab any more, and this window
%   does not host, drive or hold a handle to either.  'Export...' and 'Explore...'
%   at the bottom of the Files tab make the session on disk current and then open
%   guiExport / guiExplore ON IT (handOffSession); both open perfectly well on their
%   own with no session at all, and closing one leaves this window untouched.
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
%   token comes from the Files tab and the STAGE flag from the producing step plus
%   that type's own settings ('_t', '_s', '_c', '_e' - the stage alone, since the
%   product letter changes down the pipeline while the row's identity, the branch,
%   does not; spec D8 - falling back to the product token for a product that carries
%   no stage at all).  It is a label and never resolves a file.  Settings stay keyed by
%   (step, TYPE) - two rows of one type share them - so a divergent animal is
%   simply a second type.  Rows are built FROM THE DATA: two types or eleven work
%   with no code change, and nothing anywhere branches on a type's name.  The two
%   per-animal steps (registration, vessel typing) span the animal instead and get
%   their own box; which BRANCH FILE of each animal's reference recording every
%   reference-taking step resolves to is reported in the selection summary.
%
%   ONE SESSION MAY HOLD SEVERAL MODALITIES, AND THE QUESTION IS ASKED PER SCOPE.
%   The step registry is the UNION over the modalities the working set actually
%   holds, so a folder of speckle recordings and videos keeps both sets of columns
%   (it used to be filtered by the statistical MODE of the two, which silently threw
%   the minority away).  Below that, each scope asks for itself: the CONSTRUCTOR
%   configures a type with its own modality's registry (regFor), which is what makes
%   a row resolve to exactly one producer and a step serving two modalities pull in
%   the right one; a cell whose column belongs to another modality is dead and says
%   so.  Per FILE, wbStateEngine gates on the file's own modality - the one the user
%   set on the Files tab.  Nothing outside wbFileModel and the registry's
%   'modalities' field ever names a modality.
%
%   THREE TABS, THREE PROGRAMS.  This window PROCESSES; it neither exports nor
%   plots.  Export and Explore are separate windows (guiExport, guiExplore), so
%   browsing yesterday's results never means loading the heavy processing GUI.
%
%   Tabs: 1 Files (scan a root recursively for an extension - or add files by
%                    hand - then CURATE: sort, delete, label every file (animal /
%                    type / index / group / modality, in place or in bulk over a
%                    selection), and pin ONE reference recording per animal.  The
%                    other tabs stay locked until every field is assigned and the
%                    file names are unique - see fileProblems.  The session lives
%                    here, and so do the two buttons that hand it to the read-only
%                    tools) *
%         2 Constructor (type x raw step on top with the per-animal steps beside
%                    it, one (type, product) row per pipeline below, and the
%                    per-(step,type) settings on the right; the summary carries the
%                    resolved reference branches and the warnings) *
%         3 Process (READ-ONLY, and shaped like the Constructor: the From/To range
%                    beside Run (and Wipe all), then TWO state tables - recordings
%                    x raw steps on top BESIDE the run log, (recording, product)
%                    rows x derived steps below, the latter built from the
%                    CONFIGURATION so a product row appears with its producer rather
%                    than with its file (spec D7) - each cell showing
%                    '.'/queued/running NN%/done/error/skipped; the selected row's
%                    file below them; and on the right the Reports panel - a kind
%                    filter, the 'Create PDF reports' switch and the link list.  All
%                    selection lives in the Constructor - spec D6/D7) *
%
%   THE MONITOR'S COLUMNS SIZE THEMSELVES (spec D5).  Both state tables hand their
%   widths to the layout engine rather than naming pixels: every real column is
%   'fit', so a header like 'Dynamic segmentation' and a row label like a file name
%   plus its product flag are never cropped, and one EMPTY trailing column carries a
%   '1x' that swallows whatever is left over - so a table with room to spare covers
%   its block instead of stopping halfway across it, and one whose content is wider
%   scrolls sideways with nothing cropped.  The spacer is empty because a '1x' has
%   no minimum and clips whatever column it is put on.  The engine re-runs all of it
%   on every resize, so the window has no resize code of its own.
%
%   NEXT AND EXIT NEVER MOVE.  'Next' (Files -> Constructor -> Process) and 'Exit
%   workbench' sit in the same bottom-right corner of every tab - the last tab has
%   no Next - so the settings and Reports panels each stop one row short of the
%   bottom to leave that strip free (tabFooter).  Exit is one decision, not three:
%   it requests the same cooperative Stop the Stop button does, writes the session,
%   and closes the moment the executor has let go.
%
%   WIPE ALL throws away the processing and keeps the recordings: every _d/_r/_s.mat
%   of the files listed on the Files tab, counted and confirmed before anything goes
%   (wipeTargets fences on the whole identity, so 'Foo' never carries off 'Foo2').
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
%           wbRefBranch, wbRunRange, wbExecutor, wbModalGuard, wbUiGuard,
%           wbArtifacts, makeReportPdf, wbPool, wbSession, guiExport, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function h = guiWorkbench(varargin)

vis = 'on';
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i},'Visible'), vis = varargin{i+1}; end
end

app = struct();
app.reg          = wbStepRegistry('LSCI');       % the UNION over the modalities loaded
app.modality     = 'LSCI';                       % the working set's most common one
app.typeMod      = containers.Map('KeyType','char','ValueType','char'); % type -> modality
app.regCache     = containers.Map('KeyType','char','ValueType','any');  % modality -> registry
app.sm           = wbSettingsModel('new');        % settings bag + overrides
app.disc         = emptyDisc();
app.rows         = emptyRows();                    % flat, animal-major, reference-first
app.animalNames   = {};
% ---- the curated working set (Files tab) ----
app.files        = emptyFiles();                   % ONE ENTRY PER FILE (rows dedup by identity)
app.root         = '';                             % scan root
app.resultsRoot  = '';                             % where the results go ('' = beside them)
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
app.stale        = containers.Map('KeyType','char','ValueType','any');   % session INVALIDATION overlay
app.runState     = containers.Map('KeyType','char','ValueType','any');   % transient run overlay: 'running'|'done'|'error'
app.cellMsg      = containers.Map('KeyType','char','ValueType','any');   % 'identity||stepId' -> error/tooltip text
app.completed    = containers.Map('KeyType','char','ValueType','any');   % 'path||stepId' -> completion record (session)
app.pRows        = emptyProgressRows();            % TOP table: ONE PER FILE of the working set
app.pRowOf       = containers.Map('KeyType','char','ValueType','double');% path -> raw-table row
app.dRows        = emptyProgressRows();            % BOTTOM table: one per (recording, PRODUCT)
app.dRowOf       = containers.Map('KeyType','char','ValueType','double');% 'identity||branch' -> row
app.pSel         = struct('kind','raw','idx',0);   % the monitor row whose detail is shown
app.rangeFrom    = wbRunRange('lastValid');        % Process tab: From (sentinel = the frontier)
app.rangeTo      = '';                             % Process tab: To ('' = the last column)
app.reprocess    = false;                          % Process tab: re-run finished work (D6)
app.results      = emptyResults();                 % report artifacts this session produced
app.resultFilter = 'all';                          % result list: 'all' | 'images' | 'pdfs'
app.pdfReports   = false;                          % Create PDF reports (UI ONLY - spec D9)
app.deleteAfterPdf = false;                        % ... and then remove the images (Q3: default OFF)
app.colTotal     = containers.Map('KeyType','char','ValueType','double'); % stepId -> entries THIS run
app.colDone      = containers.Map('KeyType','char','ValueType','double'); % stepId -> entries finished
app.colArt       = containers.Map('KeyType','char','ValueType','any');    % stepId -> its artifact paths
app.reportPdfs   = {};                             % the per-column PDFs assembled so far
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
app.exitAfterStop = false;                         % Exit pressed during a run

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
app.tg = tg;
setappdata(fig,'app',app);

buildFilesTab(fig);
buildConstructorTab(fig);
buildProcessTab(fig);
% NOT wrapped in wbUiGuard (spec D4): this fires AFTER MATLAB has already moved to
% the other tab, so it is a notification rather than a command.  Dropping it would
% leave the tab switched and the monitor showing the previous tab's picture, which
% is worse than either doing it or not doing it.
tg.SelectionChangedFcn = @(~,~) onTabSelected(fig);   % keep the monitor in step

% NOT wrapped in wbUiGuard (spec D4), for the same reason the Exit button is not:
% closing the window must stay possible while a run holds the latch, and
% requestExit already turns a mid-run close into a cooperative stop.
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
    'renameRecording',@(p,stem) renameRecording(fig,p,stem), ...
    'renamePlan',  @(p,stem) wbRename('plan', ...
                        wbFileModel(p,getApp(fig).root,getApp(fig).resultsRoot),stem), ...
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
    'rowSelection',@(ty,br) wbTypeSelection('steps',getApp(fig).typeSel,regFor(getApp(fig),ty),ty,br), ...
    'rowInheritedSteps',@(ty,br) wbTypeSelection('inherited',getApp(fig).typeSel,regFor(getApp(fig),ty),ty,br), ...
    'rowOffers',   @(br,stepId) wbTypeSelection('offers',getApp(fig).reg,br,stepId), ...
    'rowWhyNot',   @(br,stepId) wbTypeSelection('why',getApp(fig).reg,br,stepId), ...
    ... % ---- what ONE TYPE is configured with: its modality and its own steps ----
    'typeModality',@(ty) modalityOfType(getApp(fig),ty), ...
    'typeStepsFor',@(ty) {regFor(getApp(fig),ty).id}, ...
    'typeBranches',@(ty) wbTypeSelection('branches',regFor(getApp(fig),ty)), ...
    'typeOffers',  @(ty,br,stepId) stepInReg(regFor(getApp(fig),ty),stepId) && ...
                                   wbTypeSelection('offers',regFor(getApp(fig),ty),br,stepId), ...
    'modalities',  @() uniqueStable({getApp(fig).modelArr.modality}), ...
    'constructorRows',@() constructorRows(getApp(fig)), ...
    'rowFlag',     @(ty,br) rowFlagFor(getApp(fig),ty,br), ...
    'typeRows',    @(ty) wbTypeSelection('rows',getApp(fig).typeSel,regFor(getApp(fig),ty),ty), ...
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
    ... % ---- Process tab: the From/To run range (spec D1-D5) ----
    'runColumns',  @() runRangeColumns(fig), ...
    'frontier',    @() runRangeFrontier(fig), ...
    'setRange',    @(f,t) setRunRange(fig,f,t), ...
    'setReprocess',@(tf) setReprocess(fig,tf), ...
    'reprocess',   @() getApp(fig).reprocess, ...
    'range',       @() runRangeValue(fig), ...
    'runOrderIn',  @(f,t) runOrderIn(fig,f,t), ...
    'lastValid',   @() wbRunRange('lastValid'), ...
    ... % ---- Process tab: the read-only monitor, its results, the session ----
    'progressRows',@() {getApp(fig).pRows.path}, ...
    'progressData',@() getApp(fig).c.process.rawTable.Data, ...
    'rawProgressData',@() getApp(fig).c.process.rawTable.Data, ...
    'derivedProgressData',@() getApp(fig).c.process.derivedTable.Data, ...
    'derivedProgressRows',@() getApp(fig).dRows, ...
    'monitorWidths',@(kind) monitorWidthPolicy(fig,kind), ...
    'progressCell',@(p,stepId) progressCellOf(fig,p,stepId), ...
    'derivedCell', @(id,br,stepId) derivedCellOf(fig,id,br,stepId), ...
    'fileState',   @(p,stepId) fileStateOf(getApp(fig),p,stepId), ...
    'plannedSteps',@(p) plannedStepsOf(fig,p), ...
    'resultLinks', @() resultLinks(fig), ...
    ... % ---- Process tab: the reports panel (spec D9/D10) ----
    'setPdfReports',@(tf) setPdfReports(fig,tf), ...
    'pdfReports',  @() getApp(fig).pdfReports, ...
    'setDeleteAfterPdf',@(tf) setDeleteAfterPdf(fig,tf), ...
    'deleteAfterPdf',@() getApp(fig).deleteAfterPdf, ...
    'resultKinds', @() resultKinds(fig), ...
    'setResultFilter',@(kind) setResultFilter(fig,kind), ...
    'resultFilter',@() getApp(fig).resultFilter, ...
    'reportPdfs',  @() reportPdfList(fig), ...
    'sessionPath', @() getApp(fig).sessionPath, ...
    'autosave',    @() autosaveSession(fig), ...
    'lastSessionPath',@() lastSessionPath(), ...
    'wipeTargets', @() wipeTargets(getApp(fig)), ...
    'wipeAll',     @() wipeAll(fig), ...
    'completed',   @() getApp(fig).completed, ...
    ... % ---- the hand-off to the two standalone tools (spec §5) ----
    'handOff',     @(tool) handOffSession(fig,tool), ...
    'identities',  @() {getApp(fig).rows.identity}, ...
    'getApp',      @() getApp(fig), ...
    'exit',        @() requestExit(fig));   % the BUTTON's contract: stop, save, close
setappdata(fig,'workbenchAPI',api);

renderProgress(fig);
renderConstructor(fig);
if nargout>0, h = fig; end
end

%% ===================== app-state helpers ============================ %%
function app = getApp(fig), app = getappdata(fig,'app'); end
function setApp(fig,app), setappdata(fig,'app',app); end
function h = ui(fig, fcn)
%ui  ONE USER ACTION AT A TIME (spec D4).  Every command callback in this window
%   goes through wbUiGuard, so a click made while another one is still running -
%   which, during a run, means every click the user made at a frozen window and
%   the next drawnow then replayed - is dropped rather than queued.  The inner
%   handle keeps its own signature; this one forwards whatever MATLAB passes.
%   Stop, Exit and the close button are deliberately NOT wrapped: see their sites.
h = wbUiGuard('wrap', fig, fcn);
end
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
%emptyProgressRows  A row of EITHER monitor table, run-ordered.  'kind' says which:
%   'raw'     - one per FILE of the working set (what happens to the recording);
%   'derived' - one per (recording, PRODUCT) row the Constructor built, so a
%               two-pipeline recording appears twice, once per branch (spec D6).
r = struct('kind',{},'path',{},'model',{},'identity',{},'label',{},'rowLabel',{}, ...
    'animal',{},'type',{},'expGroup',{},'branch',{},'flag',{},'isRef',{}, ...
    'animalIdx',{},'rowInAnimal',{},'branchIdx',{});
end
function r = emptyResults()
%emptyResults  The result-link list (spec D5/D10): report artifacts, newest first.
%   An entry carries the STEP that produced it and its KIND, so the panel's filter
%   is a property of the data rather than a string test in the render function -
%   and so a column's own images can be collected without asking the file name.
r = struct('path',{},'label',{},'when',{},'stepId',{},'kind',{});
end
function k = artifactKind(pth)
%artifactKind  'pdf' or 'image' - the only distinction the result list draws (D10).
[~,~,e] = fileparts(char(pth));
k = ternary(strcmpi(e,'.pdf'),'pdf','image');
end
function d = presetDirDefault()
d = fullfile(prefdir,'guiWorkbenchPresets');
end
function d = sessionDir()
%sessionDir  THE workbench sessions folder, created on demand.  Sessions belong
%   with the library rather than with the data - a protocol is re-scanned from
%   several roots and a session must not go missing with one of them - so they live
%   in 'workbench-sessions' beside the code and the folder is gitignored.  A
%   read-only install falls back to the MATLAB preferences directory, which is the
%   only place a session is guaranteed to be writable.
d = fullfile(fileparts(fileparts(mfilename('fullpath'))),'workbench-sessions');
if isfolder(d), return; end
[ok,~] = mkdir(d);
if ~ok
    d = fullfile(prefdir,'guiWorkbenchSessions');
    if ~isfolder(d), mkdir(d); end
end
end
function p = lastSessionPath()
%lastSessionPath  The LAST SESSION: one fixed file, rewritten on every state change
%   and on exit, so re-opening the workbench after a crash, a Stop or a plain close
%   always has somewhere to resume from without anyone having pressed Save.
p = fullfile(sessionDir(),'lastSession.mat');
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

uilabel(lp,'Text',['Collect the recordings you want to process - scan a folder, or add files ' ...
    'by hand.  Then sort them, delete what you do not need, label every one, and mark ' ...
    'one reference recording for each animal.  Typing over a file name renames that ' ...
    'recording on your disk.'],'FontWeight','bold','WordWrap','on');

% -- source row: root + results root + glob + the loaders --
%   TWO FOLDERS, and the second one is optional: leave 'results go in' matching
%   'look in' and everything is written beside the recordings, which is what this
%   library has always done.  Point it somewhere else and the recordings are never
%   written to - the subfolders inside it mirror the ones inside 'look in'.
sp = uigridlayout(lp,[1 11], ...
    'ColumnWidth',{'fit','2x','fit','fit','2x','fit','fit','0.9x','fit','fit','fit'}, ...
    'Padding',[0 0 0 0]);
uilabel(sp,'Text','look in');
c.root = uieditfield(sp,'text','Value','','Tooltip', ...
    'The folder to search.  Every subfolder inside it is searched too.', ...
    'ValueChangedFcn',ui(fig,@(s,~)onSourceEdit(fig)));
uibutton(sp,'Text','Browse...','ButtonPushedFcn',ui(fig,@(~,~)uiBrowseRoot(fig)));
uilabel(sp,'Text','results go in');
c.resultsRoot = uieditfield(sp,'text','Value','','Tooltip', ...
    ['Where the processed files and the reports are written.  It follows the folder ' ...
     'you search until you change it; set it elsewhere to leave your recordings ' ...
     'untouched, and the subfolders inside it will match.'], ...
    'ValueChangedFcn',ui(fig,@(s,~)onSourceEdit(fig)));
uibutton(sp,'Text','Browse...','ButtonPushedFcn',ui(fig,@(~,~)uiBrowseResults(fig)));
uilabel(sp,'Text','for files named');
c.glob = uieditfield(sp,'text','Value',app.glob,'ValueChangedFcn',ui(fig,@(s,~)onSourceEdit(fig)), ...
    'Tooltip',tipGlob());
uibutton(sp,'Text','Scan','BackgroundColor',[0.82 0.92 0.82],'ButtonPushedFcn',ui(fig,@(~,~)uiScan(fig)), ...
    'Tooltip','Search the folder and start a new list from what it finds.');
uibutton(sp,'Text','Add files...','ButtonPushedFcn',ui(fig,@(~,~)uiAddFiles(fig)), ...
    'Tooltip','Pick files yourself and add them to the list.  You can select several at once.');
uibutton(sp,'Text','Add folder...','ButtonPushedFcn',ui(fig,@(~,~)uiAddFolder(fig)), ...
    'Tooltip','Search another folder and add what it finds to the list.');
uibutton(sp,'Text','Clear','ButtonPushedFcn',ui(fig,@(~,~)uiClear(fig)), ...
    'Tooltip','Empty the list.  No file is deleted from your disk.');

% -- regexp row: one box per label axis + the reference rule --
rp = uigridlayout(lp,[2 5],'RowHeight',{'fit','fit'},'Padding',[0 0 0 0],'RowSpacing',1);
axes5 = {'animal','Animal','[A-Z]+\d+',tipAnimal(); 'type','Recording type','',tipType(); ...
         'index','Recording number','',tipIndex(); 'expGroup','Experimental group','',tipExpGroup(); ...
         'ref','Reference recording','',tipRef()};
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
    'CellEditCallback',ui(fig,@(s,e)onFileTableEdit(fig,s,e)));
c.fileTbl.ColumnFormat{1} = 'logical';              % the reference tick

% -- quick assign: give EVERY selected row the same value in one go -----------
cp = uigridlayout(lp,[1 7], ...
    'ColumnWidth',{'fit','fit','fit','1.2x','fit','fit','2x'},'Padding',[0 0 0 0]);
uilabel(cp,'Text','Quick assign: set','FontWeight','bold');
c.assignAxis = uidropdown(cp,'Items',assignableColumns(),'Value','type', ...
    'Tooltip','Which label the selected rows should get.', ...
    'ValueChangedFcn',ui(fig,@(~,~)refreshAssignItems(fig)));
uilabel(cp,'Text','=');
c.assignVal = uidropdown(cp,'Items',{''},'Value','','Editable','on', ...
    'Tooltip',['Pick a value you have used before, or type a new one.  The names are ' ...
               'entirely yours - nothing here has a fixed list.']);
uibutton(cp,'Text','Apply to selected rows','BackgroundColor',[0.85 0.9 1], ...
    'ButtonPushedFcn',ui(fig,@(~,~)uiAssignLabel(fig)), ...
    'Tooltip',['Give every selected row this value at once.  Select rows by clicking ' ...
               'them in the table, shift-click for a range.']);
uibutton(cp,'Text','Delete selected','BackgroundColor',[1 0.88 0.82], ...
    'ButtonPushedFcn',ui(fig,@(~,~)uiDeleteSelected(fig)), ...
    'Tooltip','Take the selected rows out of this list.  Nothing is deleted from your disk.');
c.curStatus = uilabel(cp,'Text','','FontAngle','italic','FontSize',10);

c.problems = uilabel(lp,'Text','','WordWrap','on','FontWeight','bold');
c.status   = uilabel(lp,'Text','No files loaded.','WordWrap','on');

% -- session row (RowHeight fit: buttons keep their natural height) -----------
%   The two hand-off buttons live HERE rather than beside Run: they are session
%   CONSUMERS, this is where the session lives, and the Process toolbar is a run
%   control - nothing that opens another window belongs in it.
ssp = uigridlayout(lp,[1 7],'RowHeight',{'fit'}, ...
    'ColumnWidth',{'fit','fit','fit','fit','1x','fit','fit'},'Padding',[0 0 0 0]);
uibutton(ssp,'Text','Save session...','ButtonPushedFcn',ui(fig,@(~,~)uiSaveSession(fig)));
uibutton(ssp,'Text','Load session...','ButtonPushedFcn',ui(fig,@(~,~)uiLoadSession(fig)));
uibutton(ssp,'Text','Export...','BackgroundColor',[0.90 0.94 1.0], ...
    'ButtonPushedFcn',ui(fig,@(~,~)handOffSession(fig,'guiExport')), ...
    'Tooltip',['Save the session and open the export tool on it, to write your results ' ...
               'to Excel.  It is a separate window and this one keeps running.']);
uibutton(ssp,'Text','Explore...','BackgroundColor',[0.90 0.94 1.0], ...
    'ButtonPushedFcn',ui(fig,@(~,~)handOffSession(fig,'guiExplore')), ...
    'Tooltip',['Save the session and open the results explorer on it, to plot and compare ' ...
               'your results.  It is a separate window and this one keeps running.']);
uilabel(ssp,'Text','the session remembers your file list, its labels, the references and the settings', ...
    'FontAngle','italic','FontSize',10);
c.nextBtn = nextButton(ssp,fig,'Next: Constructor',app.tabs.constructor, ...
    'Go on to the Constructor.  It becomes available once every file is fully labelled.');
c.exitBtn = exitButton(ssp,fig);

app.c.files = c; setApp(fig,app);
end

%% ---- the strip every tab ends with (spec: Next / Exit never move) ---- %%
function b = nextButton(parent,fig,text,tab,tip)
%nextButton  'Next' - the same button, in the same corner, on every tab but the
%   last.  It is the ONLY forward move the window offers, and it goes through the
%   same gate a tab click does, so an incomplete file set is refused identically.
b = uibutton(parent,'Text',text,'BackgroundColor',[0.82 0.92 0.82], ...
    'ButtonPushedFcn',ui(fig,@(~,~)guardTabSwitchTo(fig,tab)),'Tooltip',tip);
end
function b = exitButton(parent,fig)
%exitButton  'Exit workbench' - PERSISTENT: the same button in the same corner of
%   all three tabs, so leaving never means hunting for the tab that owns it.
%
%   NOT wrapped in ui() (spec D4).  Exit is the one command that must work while a
%   run holds the latch: requestExit already handles the mid-run case by asking for
%   the same cooperative stop the Stop button does, and dropping it would leave the
%   user with a window they cannot leave for the length of the batch.
b = uibutton(parent,'Text','Exit workbench','BackgroundColor',[1 0.82 0.82], ...
    'ButtonPushedFcn',@(~,~)requestExit(fig), ...
    'Tooltip',['Stop anything still running, save your session and close the window.  ' ...
               'The session is always saved, so you can pick up where you left off.']);
end
function c = tabFooter(parent,fig,nextText,nextTab,nextTip)
%tabFooter  The bottom strip of the Constructor and Process tabs: a stretch, then
%   Next (omitted on the last tab), then Exit - laid out so both land exactly where
%   the Files tab's own pair does.  The panel above simply ends one row early to
%   leave the strip free, which is the whole cost of the two buttons persisting.
n = 2 + double(~isempty(nextText));
fp = uigridlayout(parent,[1 n],'RowHeight',{'fit'}, ...
    'ColumnWidth',[{'1x'}, repmat({'fit'},1,n-1)],'Padding',[0 0 0 0],'ColumnSpacing',4);
uilabel(fp,'Text','');                                   % the stretch
c = struct();
if ~isempty(nextText), c.nextBtn = nextButton(fp,fig,nextText,nextTab,nextTip); end
c.exitBtn = exitButton(fp,fig);
end

function f = patField(parent,fig,axis,dflt,tip)
%patField  One regexp box, wired to the pattern struct.
f = uieditfield(parent,'text','Value',dflt,'Tooltip',tip, ...
    'ValueChangedFcn',ui(fig,@(s,~)onPatternEdit(fig,axis,s.Value)));
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
%fileTableEditable  What may be typed in place.  The extension, the parsed stage and
%   the modality are properties OF THE FILE, not labels, and are read-only.
%   THE NAME IS EDITABLE AND MEANS SOMETHING ELSE ENTIRELY: typing over it renames
%   the RECORDING on disk - every file named after it moves together - which is why
%   it is the one edit that asks first (renameRecording).
e = [true true true true true true false false false];
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

%% ---- Files-tab tooltips.  Written for the person at the microscope: what the
%%      box is for, then real examples.  Examples teach a pattern far faster than
%%      prose about patterns, so every one of these ends in a worked list. ---- %%
function t = tipGlob()
t = sprintf(['Which file names to look for.  * stands for any text.\n' ...
    '   *.rls           raw speckle recordings\n' ...
    '   *.avi           myograph videos\n' ...
    '   *.adicht        LabChart recordings\n' ...
    '   *_t_K_r.mat     contrast files\n' ...
    '   *_BFI_r.mat     blood flow files\n' ...
    'The five boxes below work differently - they are search patterns, not names.']);
end
function f = recordingFilter()
%recordingFilter  The Add-files dialog's own list: processed products plus every
%   raw container the name parser knows, so a modality added to wbFileModel shows
%   up here without a second list to remember.
exts = wbFileModel('extensions');
f = strjoin([{'*.mat'}, cellfun(@(e) ['*' e], exts, 'UniformOutput', false)], ';');
end
function t = tipAnimal()
t = sprintf(['Finds the ANIMAL - the subject - in each file name.  Each animal has one\n' ...
    'reference recording, and registration and vessel typing work animal by animal.\n' ...
    '   [A-Z]+\\d+        matches PSY01, REG12, AB3 ...\n' ...
    '   m\\d+             matches m07, m11 ...\n' ...
    '   (mouse|rat)\\d+   matches mouse4, rat12 ...\n' ...
    'A file this does not find anything in is listed as "(unassigned)" and you can\n' ...
    'type its animal in the table.  Leave the box empty to treat every file as one\n' ...
    'animal.']);
end
function t = tipType()
t = sprintf(['Finds the RECORDING TYPE - what the recording is for in your experiment.\n' ...
    'You set up the processing once per type, and every recording of that type is\n' ...
    'then processed the same way.\n' ...
    '   BV|BN|BP         one lab''s slow, stimulation and pulsatile recordings\n' ...
    '   ctrl|stim|washout   another lab''s protocol\n' ...
    '   \\d(?=BP)         whatever your own file names carry\n' ...
    'There is no built-in list of types - they are whatever you call them.  A file\n' ...
    'this does not find anything in is listed as "(untyped)", which is a type like\n' ...
    'any other, and you can type its real one in the table.']);
end
function t = tipIndex()
t = sprintf(['Finds the RECORDING NUMBER within an animal, so repeats can be told apart.\n' ...
    '   [aik]\\d           matches i2, k1, a2 ...\n' ...
    '   run(\\d+)          matches the number after "run"\n' ...
    '   \\d+(?=_c_)        matches the digits just before _c_\n' ...
    'Leave the box empty to treat every recording as the same number.']);
end
function t = tipExpGroup()
t = sprintf(['Finds the EXPERIMENTAL GROUP - the label you will compare by.  Export and\n' ...
    'Explore use it; processing ignores it, so it never changes a result.\n' ...
    '   KO|WT             genotype taken from the name\n' ...
    '   Ctrl|Stroke       condition taken from the name\n' ...
    '   pre|post          a timepoint\n' ...
    'It is independent of the animal: a group can span several animals, and one\n' ...
    'animal can appear in several groups.  A file this does not find anything in is\n' ...
    'listed as "(ungrouped)".']);
end
function t = tipRef()
t = sprintf(['Finds each animal''s REFERENCE RECORDING - the one every other recording of\n' ...
    'that animal is aligned to, whatever its type or group.\n' ...
    '   1BP_c_BFI_d\\.mat  the animal''s first pulsatile recording\n' ...
    '   _ref_             a name you mark yourself\n' ...
    'It takes effect as soon as you type it, on the files already listed - there is\n' ...
    'no need to scan again.  You can always change it afterwards by ticking the\n' ...
    'reference box in the table - your own choice wins.  An animal may have no\n' ...
    'reference at all; steps that need one will say so.']);
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
uilabel(hdr,'Text',['Tick the processing steps each RECORDING TYPE should run.  You set a ' ...
    'type up once and every recording of that type is processed the same way.  If one ' ...
    'animal needs different parameters, give it its own type on the Files tab.'], ...
    'WordWrap','on','FontWeight','bold');
tb = uigridlayout(hdr,[1 7], ...
    'ColumnWidth',{'fit','1x','fit','1x','fit','fit','1.2x'},'Padding',[0 0 0 0]);
uilabel(tb,'Text','Copy configuration from');
c.copySrc = uidropdown(tb,'Items',{''},'Value','', ...
    'Tooltip',['A ready-made protocol, or another type you have already set up.  ' ...
               'A protocol only fills the boxes in - you can still change everything.']);
uilabel(tb,'Text','to');
c.copyDst = uidropdown(tb,'Items',{''},'Value','','Tooltip','The type that will be overwritten.');
uibutton(tb,'Text','Copy','BackgroundColor',[0.85 0.9 1],'ButtonPushedFcn',ui(fig,@(~,~)uiCopyType(fig)), ...
    'Tooltip','Give the second type the same steps and the same parameters as the first.');
uibutton(tb,'Text','Reset selected type','ButtonPushedFcn',ui(fig,@(~,~)uiResetType(fig)), ...
    'Tooltip','Put the type shown on the right back to the standard parameters.');
c.status = uilabel(tb,'Text','','FontAngle','italic','FontSize',10);

c.rawPanel = uipanel(left,'Title','1 - Processing the recording itself', ...
    'FontWeight','bold');
c.derivedPanel = uipanel(left,'Title','2 - Processing each result: one row per recording type and result', ...
    'FontWeight','bold');
c.summaryPanel = uipanel(left,'Title','What will run','FontWeight','bold');

% -- RIGHT: settings for one (step, type), then the persistent Next/Exit strip --
%   the settings panel stops one row short of the bottom on purpose: that row is
%   where Next and Exit sit on every tab (tabFooter)
right = uigridlayout(gl,[5 1],'RowHeight',{'fit','fit','fit','1x','fit'},'RowSpacing',6);
sel = uigridlayout(right,[2 2],'ColumnWidth',{'fit','1x'},'RowHeight',{'fit','fit'}, ...
    'Padding',[0 0 0 0],'RowSpacing',3);
uilabel(sel,'Text','Type','FontWeight','bold');
c.typeDrop = uidropdown(sel,'Items',{'(no types)'},'Value','(no types)', ...
    'Tooltip','The recording type whose parameters you are editing below.', ...
    'ValueChangedFcn',ui(fig,@(s,~)selectConstructorType(fig,s.Value)));
uilabel(sel,'Text','Step','FontWeight','bold');
c.stepDrop = uidropdown(sel,'Items',{''},'Value','', ...
    'Tooltip','The processing step whose parameters are shown below.', ...
    'ValueChangedFcn',ui(fig,@(s,~)selectConstructorStep(fig,s.Value)));
c.stepInfo = uilabel(right,'Text','','WordWrap','on','FontAngle','italic');
c.scopeInfo = uilabel(right,'Text','','WordWrap','on','FontSize',10);
c.paramPanel = uipanel(right,'BorderType','none');   % the section stack owns the scroll
f = tabFooter(right,fig,'Next: Process',app.tabs.process, ...
    'Go on to the Process tab and run what you have set up here.');
c.nextBtn = f.nextBtn; c.exitBtn = f.exitBtn;

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
%
%   THE COLUMNS ARE THE SESSION'S, THE CELLS ARE THE TYPE'S.  A mixed working set
%   shows every modality's entry steps, so the matrix stays rectangular and the
%   user can see what the other half of the folder is being offered; a cell whose
%   step does not run on that type's recordings is dead and says why.
app = getApp(fig);
types = constructorTypes(app);
cols  = wbTypeSelection('rawSteps', app.reg);
if isempty(types) || isempty(cols)
    g = uigridlayout(parent,[1 1],'Padding',[16 16 16 16]);
    uilabel(g,'Text','Load your recordings and label them on the Files tab, then come back here.', ...
        'HorizontalAlignment','center','FontAngle','italic');
    return
end

grid = uigridlayout(parent,[1+numel(types) 1+numel(cols)], ...
    'ColumnWidth',[{170}, repmat({120},1,numel(cols))], ...
    'RowHeight',[{50}, repmat({26},1,numel(types))], ...
    'RowSpacing',2,'ColumnSpacing',2,'Padding',[2 2 2 2],'Scrollable','on');
h = uilabel(grid,'Text','Recording type','FontWeight','bold');
h.Layout.Row = 1; h.Layout.Column = 1;
for s = 1:numel(cols)
    step = stepById(app.reg,cols{s});
    hp = uigridlayout(grid,[2 1],'RowHeight',{'1x','fit'},'Padding',[1 1 1 1],'RowSpacing',1);
    hp.Layout.Row = 1; hp.Layout.Column = s+1;
    uilabel(hp,'Text',step.label,'WordWrap','on','FontSize',10,'Tooltip',rawStepTip(step));
    gr = uigridlayout(hp,[1 1],'Padding',[0 0 0 0]);
    uibutton(gr,'Text','settings','FontSize',9,'ButtonPushedFcn',ui(fig,@(~,~)selectConstructorStep(fig,step.id)));
end
for r = 1:numel(types)
    ty = types{r};
    lb = uilabel(grid,'Text',sprintf('%s - %d files', ty, typeFileCount(app,ty)), ...
        'FontWeight','bold','Tooltip','Every recording of this type is processed the same way.');
    lb.Layout.Row = r+1; lb.Layout.Column = 1;
    for s = 1:numel(cols)
        step = stepById(app.reg,cols{s});
        if ~stepInReg(regFor(app,ty),cols{s})
            cb = uicheckbox(grid,'Text','','Value',false,'Enable','off', ...
                'FontSize',10,'Tooltip',wrongModalityTip(app,ty,step));
        else
            cb = uicheckbox(grid,'Text',rowFlagFor(app,ty,step.branch), ...
                'Value',wbTypeSelection('isOn',app.typeSel,ty,step.branch,cols{s}), ...
                'FontSize',10,'Tooltip',rawCellTip(app,ty,step), ...
                'ValueChangedFcn',ui(fig,@(o,~)tickRaw(fig,ty,cols{s},o.Value)));
        end
        cb.Layout.Row = r+1; cb.Layout.Column = s+1;
    end
end
end

function t = rawStepTip(step)
t = sprintf(['%s reads the raw recording and writes a new, separate result.\n' ...
    'Tick it and this type gets its own %s row below, where you choose everything ' ...
    'that is then done to that result.'], step.label, step.branch);
end
function t = rawCellTip(app,type,step)
t = sprintf('%s produces %s for every %s recording.', step.label, ...
    rowFlagFor(app,type,step.branch), type);
end

function lines = tickRaw(fig,type,stepId,tf)
%tickRaw  Switch one raw producer on/off for a TYPE - i.e. create or remove that
%   type's branch row.  The row's own configuration is KEPT when it goes away, so
%   re-ticking the producer brings it back exactly as it was, and the OTHER row is
%   never touched.
app = getApp(fig);
type = char(type); stepId = char(stepId);
reg  = regFor(app,type);                  % this type's own steps, not the session's
step = stepById(reg,stepId);
if isempty(step), lines = {}; return; end
br = step.branch;
app.typeSel = wbTypeSelection('set', app.typeSel, reg, type, br, stepId, logical(tf));
setApp(fig,app);
lines = {sprintf('constructor: [%s] %s %s -> %s row %s', type, ternary(logical(tf),'+','-'), ...
    stepId, rowFlagFor(getApp(fig),type,br), ternary(logical(tf),'added','removed'))};
if ~tf
    kept = wbTypeSelection('steps', app.typeSel, reg, type, br);
    if isempty(kept)
        lines{end+1} = sprintf('  its configuration is kept - re-tick %s to get the row back', stepId);
    end
end
for i = 1:numel(lines), wbLog(fig,lines{i}); end
renderConstructor(fig);
recomputeBase(fig);       % the configuration is an INPUT to the states (see statesOf)
renderProgress(fig);      % D7: the product row appears with its producer, not with its file
autosaveSession(fig);
end

function lines = tickRow(fig,type,branch,stepId,tf)
%tickRow  Tick / untick one cell of a (type,branch) row and run the cascade.
%   TICK pulls the prerequisites in, so a chain is constructible in one click;
%   UNTICK pushes the dependants out, since they could no longer run.  Both stay
%   INSIDE the row - the other branch of the same recording is a separate pipeline -
%   and both are logged, because a box the user did not click must never move
%   silently.  A TICK CAN ALSO MOVE A BOX THE OTHER WAY: a step the registry says
%   cannot run beside this one is unticked, so what came back is sorted by where
%   each box ended up rather than assumed to be all prerequisites.
app = getApp(fig);
type = char(type); branch = char(branch); stepId = char(stepId);
reg  = regFor(app,type);                  % this type's own steps, not the session's
flag = rowFlagFor(app,type,branch);
[app.typeSel, changed, animalIds] = wbTypeSelection('set', app.typeSel, reg, ...
    type, branch, stepId, logical(tf));
lines = {};
if tf
    lines{end+1} = sprintf('constructor: [%s %s] + %s', type, flag, stepId);
    % 'effective', not 'isOn': a recording-level prerequisite is stored on the
    % anchor row, so asking this row's own key about it would read as unticked
    on  = changed(cellfun(@(id) wbTypeSelection('effective',app.typeSel,reg,type,branch,id), changed));
    off = changed(~ismember(changed,on));
    if ~isempty(off)
        % The one line here written in LABELS rather than ids: a box the user did
        % not click has just been cleared, and the sentence explaining why is read
        % by the person who ticked the other one.
        lines{end+1} = sprintf('  unticked %s - it cannot be selected together with %s', ...
            stepWords(reg,off), stepWords(reg,{stepId}));
    end
    changed = on;
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
recomputeBase(fig);       % the configuration is an INPUT to the states (see statesOf)
renderProgress(fig);      % what a row is queued for changed; so did the run columns
autosaveSession(fig);
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
for i = 1:numel(types)
    reg = regFor(app,types{i});                 % a type's rows are ITS modality's
    brs = wbTypeSelection('branches', reg);
    for b = 1:numel(brs)
        if ~wbTypeSelection('rowOn', app.typeSel, reg, types{i}, brs{b}), continue; end
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
%rowFlagFor  The FLAG a (type,branch) row's files carry - '_t', '_s', '_c', '_e',
%   '_MYO' (spec D8).  Never a literal: it is read off the PRODUCING STEP's own
%   outSuffix, and a producer whose stage is a SETTING (the contrast step's
%   contrastType) is asked for this type's answer.  Row identity is (type, BRANCH),
%   so switching a project from temporal to spatial contrast changes the flag shown,
%   not the configuration.
%
%   THE STAGE ALONE, NOT STAGE+PRODUCT.  The product letter changes down the
%   pipeline ('_t_K' becomes '_t_BFI' at the BFI step) while the row's identity is
%   the BRANCH, which does not - so showing '_t_K' on a row that will end up holding
%   a '_t_BFI' was over-specified.
%
%   UNLESS THERE IS NO STAGE.  A product that carries only a product token and no
%   stage flag ('_MYO': one PAIR per recording, appended to in place) would come
%   out of the stage rule as a bare '_', so the product token stands for it.  Same
%   rule, one fall-back, and it is what keeps this from being contrast-specific.
%
%   This is a DISPLAY label and nothing else: no file is ever resolved through it
%   (that goes through branchScope / wbFileModel / contrastStageForModel), so
%   neither the shortening nor the fall-back can move a file.
f = '';
reg = regFor(app,type);
pid = wbTypeSelection('producer', reg, branch);
if isempty(pid), return; end
step = stepById(reg,pid);
if isempty(step.outSuffix), return; end
m  = wbFileModel(['x' step.outSuffix{1} '.mat']);   % '_t_K_d' -> stage t, '_MYO_r' -> product MYO
st = m.stage;
alt = settingStage(app,type,step);
if ~isempty(alt), st = alt; end
if isempty(st), st = m.product; end                 % no stage to show: name the product
if isempty(st), return; end
f = ['_' st];
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
tf = wbTypeSelection('isOn', app.typeSel, type, branchOfStep(regFor(app,type),stepId), stepId);
end
function tf = rowShown(app,type,branch,stepId)
%rowShown  What the cell DISPLAYS: its own tick, or the one it inherits.
tf = wbTypeSelection('effective', app.typeSel, regFor(app,type), type, branch, stepId);
end
function tf = rowInherited(app,type,branch,stepId)
%rowInherited  Is this cell showing a tick that belongs to the copy-source row?
[~,tf] = wbTypeSelection('effective', app.typeSel, regFor(app,type), type, branch, stepId);
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
    uilabel(g,'Text','Load your recordings and label them on the Files tab first.', ...
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
h = uilabel(grid,'Text','Recording type and result','FontWeight','bold');
h.Layout.Row = 1; h.Layout.Column = 1;
for s = 1:numel(cols)
    step = stepById(app.reg,cols{s});
    hp = uigridlayout(grid,[2 1],'RowHeight',{'1x','fit'},'Padding',[1 1 1 1],'RowSpacing',1);
    hp.Layout.Row = 1; hp.Layout.Column = s+1;
    uilabel(hp,'Text',step.label,'WordWrap','on','FontSize',10,'Tooltip',headerTip(step));
    gr = uigridlayout(hp,[1 1],'Padding',[0 0 0 0]);
    uibutton(gr,'Text','settings','FontSize',9,'ButtonPushedFcn',ui(fig,@(~,~)selectConstructorStep(fig,step.id)));
end
for r = 1:numel(rows)
    row = rows(r);
    txt = row.label;
    if preview, txt = sprintf('%s - tick a step above first', row.type); end
    lb = uilabel(grid,'Text',sprintf('%s - %d files', txt, row.files),'FontWeight','bold', ...
        'Tooltip','One pipeline - this type''s recordings, and this one result of them.');
    lb.Layout.Row = r+1; lb.Layout.Column = 1;
    if preview, lb.FontColor = [0.45 0.45 0.45]; end
    for s = 1:numel(cols)
        cb = makeRowCell(fig,grid,app,row,cols{s},preview);
        cb.Layout.Row = r+1; cb.Layout.Column = s+1;
    end
end
end

function cb = makeRowCell(fig,grid,app,row,stepId,preview)
%makeRowCell  One derived cell: live, greyed-with-a-reason, or inherited.  Every
%   question is asked of THIS ROW'S TYPE registry, so a column belonging to another
%   modality is dead here and says so, and the row rule ('why') never has to answer
%   for a step this type does not have.
reg = regFor(app,row.type);
if ~stepInReg(reg,stepId)
    cb = uicheckbox(grid,'Text','','Value',false,'Enable','off', ...
        'Tooltip',wrongModalityTip(app,row.type,stepById(app.reg,stepId)));
    return
end
why = wbTypeSelection('why', reg, row.branch, stepId);
if preview || ~isempty(why)
    if isempty(why), why = 'Tick a step in the panel above first, to give this type a result to work on.'; end
    cb = uicheckbox(grid,'Text','','Value',false,'Enable','off','Tooltip',why);
    return
end
[tf,inh] = wbTypeSelection('effective', app.typeSel, reg, row.type, row.branch, stepId);
if inh
    src = wbTypeSelection('anchorBranch', reg);
    cb = uicheckbox(grid,'Text','','Value',tf,'Enable','off', ...
        'Tooltip',sprintf('%s  Tick it on the %s row instead.', ...
        sharedReason(stepById(reg,stepId)), rowFlagFor(app,row.type,src)));
    return
end
cb = uicheckbox(grid,'Text','','Value',tf,'Tooltip',cellTip(reg,stepId), ...
    'ValueChangedFcn',ui(fig,@(o,~)tickRow(fig,row.type,row.branch,stepId,o.Value)));
end

function t = sharedReason(step)
%sharedReason  Why a step is ticked once for the whole recording instead of per
%   row - straight from its branchScope, which is also what the executor obeys.
switch step.branchScope
    case 'copy'
        t = sprintf(['%s is done once, on the contrast result, and carried over to the ' ...
            'other results of the same recording.'], step.label);
    otherwise
        t = sprintf(['%s covers every result of a recording in one go, so it cannot be ' ...
            'on for one result and off for another.'], step.label);
end
end

function t = derivedHint(rows)
if isempty(rows)
    t = ['Nothing to set up yet - tick a step in the panel above.  Each one gives its ' ...
         'type a result, and each result gets its own row here.'];
else
    t = ['One row is one pipeline: a type''s recordings and one result of them.  A type ' ...
         'running both steps above gets two rows, and they are independent.'];
end
end

function t = cellTip(reg,stepId)
step = stepById(reg,stepId);
t = step.label;
req = wbPrereqs('describe',step);
if ~isempty(req)
    t = [t sprintf('\nIt needs %s first.  Ticking this ticks those too, on this row.', req)];
end
cw = step.conflictsWith;
if ~isempty(cw)
    t = [t sprintf('\nIt is an alternative to %s - ticking this unticks that.', ...
        stepWords(reg,cw))];
end
if isfield(step,'note') && ~isempty(step.note)
    % The one thing about this step a person has to know before ticking it, in the
    % registry's own words - there is nowhere else the Constructor could read it.
    t = [t sprintf('\n%s', step.note)];
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
p  = uipanel(parent,'Title','Done once per animal','FontWeight','bold','FontSize',10);
n  = max(1,numel(ids));
gl = uigridlayout(p,[n+1 2],'ColumnWidth',{'1x','fit'}, ...
    'RowHeight',[repmat({'fit'},1,n), {'1x'}], ...
    'RowSpacing',6,'ColumnSpacing',6,'Padding',[8 8 8 8]);
for i = 1:numel(ids)
    step = stepById(app.reg,ids{i});
    uicheckbox(gl,'Text',step.label,'Value',isKey(app.animalSel,step.id), ...
        'Tooltip',animalStepTip(step), ...
        'ValueChangedFcn',ui(fig,@(o,~)tickAnimalStep(fig,step.id,o.Value)));
    uibutton(gl,'Text','settings','FontSize',9, ...
        'Tooltip','Open this step''s parameters in the panel on the right.', ...
        'ButtonPushedFcn',ui(fig,@(~,~)selectConstructorStep(fig,step.id)));
end
uilabel(gl,'Text',['These run once per animal, over all of its recordings whatever ' ...
    'their type.  Their parameters are the same for every animal.'],'WordWrap','on', ...
    'FontSize',9,'FontColor',[0.4 0.4 0.4]);
end

function t = animalStepTip(step)
t = sprintf(['%s runs once per animal, over all of its recordings whatever their type.\n' ...
    'It reads the %s of that animal''s reference recording.'], step.label, refBranchWord(step));
end
function w = refBranchWord(step)
switch step.refBranch
    case 'contrast', w = 'contrast result, _t or _s';
    case 'cardiac',  w = 'cardiac-cycle result, _c';
    otherwise,       w = 'first result available';
end
end
function s = warnLine(w)
if isempty(w), s = 'Every step you have ticked can find the reference recording it needs.';
else,          s = ['Warning: ' strjoin(w,'   |   ')]; end
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
recomputeBase(fig);       % the configuration is an INPUT to the states (see statesOf)
renderProgress(fig);      % an animal step is a monitor column and a run column too
autosaveSession(fig);
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
    reg = regFor(app,ty);                            % each type answers for itself
    if ~stepInReg(reg, step.id), continue; end       % it does not run on this type at all
    brs = wbTypeSelection('rows', app.typeSel, reg, ty);
    if isempty(brs)                                  % no product yet: make one
        pid = defaultProducerFor(reg, step.id);
        if ~isempty(pid)
            app.typeSel = wbTypeSelection('tick', app.typeSel, reg, ty, ...
                branchOfStep(reg,pid), pid);
            ticked = [ticked, {pid}]; %#ok<AGROW>
        end
        brs = wbTypeSelection('rows', app.typeSel, reg, ty);
    end
    for b = 1:numel(brs)
        have = wbTypeSelection('steps', app.typeSel, reg, ty, brs{b});
        need = wbPrereqs('missing', stepById(reg,step.id), have);  % nothing when the row feeds it
        for r = 1:numel(need)
            if ~wbTypeSelection('offers', reg, brs{b}, need{r}), continue; end
            app.typeSel = wbTypeSelection('tick', app.typeSel, reg, ty, brs{b}, need{r});
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
brs = wbTypeSelection('rows', app.typeSel, regFor(app,type), type);
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
% not loaded, but the identity still locates it - and it is given the two project
% folders like any other model, or the reference lookup below it would go hunting
% for products in the raw tree
m = wbFileModel(identity, app.root, app.resultsRoot);
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
        lines{end+1} = sprintf('%s: no reference recording marked', animals{a}); %#ok<AGROW>
        continue
    end
    bits = cell(1,numel(p));
    for i = 1:numel(p)
        switch p(i).status
            case 'ok',       bits{i} = sprintf('%s %s', p(i).stepLabel, p(i).name);
            case 'fallback', bits{i} = sprintf('%s %s - second choice', p(i).stepLabel, p(i).name);
            otherwise,       bits{i} = sprintf('%s - no file yet', p(i).stepLabel);
        end
    end
    lines{end+1} = sprintf('%s: reference %s -> %s', animals{a}, p(1).refLabel, strjoin(bits,' | ')); %#ok<AGROW>
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
        w{end+1} = sprintf('%s: %s is ticked but this animal has no reference recording', ...
            p.animal, p.stepLabel); %#ok<AGROW>
    elseif any(strcmp(p.status,{'fallback','missing'})) && ~willBeProduced(app,p)
        w{end+1} = sprintf('%s: %s - %s', p.animal, p.stepLabel, p.note); %#ok<AGROW>
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
sel = plannedStepIds(app, ty);
tf  = wbPrereqs('met', step, sel);
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

function ids = plannedStepIds(app, type, branches)
%plannedStepIds  WHAT THIS RUN WILL PRODUCE, for one recording type: the per-animal
%   steps that are on, plus - for each (type,branch) row - the steps ticked on that
%   row and the ones it inherits from its anchor, in registry order.
%
%   ONE DEFINITION, because three different questions are asked of it and they may
%   never disagree: the Constructor's "is that missing file about to be made?"
%   (constructorWarnings > willBeProduced), the monitor's "what is queued for this
%   row?" (plannedStepsFor), and the state engine's look-ahead - a step whose
%   producer runs earlier in the same sequence has its input by the time it is
%   reached, and must not be reported as having none (statesOf).
%
%   'branches' narrows it to particular rows and may legitimately be EMPTY (a
%   monitor row whose branch its type does not run); omit the argument entirely for
%   the whole type.
ids = animalStepsOn(app);
if isempty(type), ids = orderIds(app.reg, ids); return; end
reg = regFor(app,type);
if nargin<3
    branches = wbTypeSelection('rows', app.typeSel, reg, char(type));
end
for b = 1:numel(branches)
    ids = [ids, wbTypeSelection('steps',     app.typeSel, reg, type, branches{b}), ...
                wbTypeSelection('inherited', app.typeSel, reg, type, branches{b})]; %#ok<AGROW>
end
ids = orderIds(app.reg, ids);            % ordered by the SESSION's registry: one order
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
    c.scopeInfo.Text = sprintf(['%s runs once per animal, over all of its recordings ' ...
        'whatever their type, so these values are the same everywhere - there is no ' ...
        'separate set per type.'], step.label);
    buildSettingsPanel(fig,c.paramPanel,step,'');
else
    c.scopeInfo.Text = sprintf(['These values apply to every recording of type "%s" - %d ' ...
        'of them, and to both its result rows if it has two - and to no other type.  ' ...
        'Changing one marks the later steps of this type as out of date.'], ...
        app.selType, typeFileCount(app,app.selType));
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
    seed = stepId;                     % STEP-ID keyed, never field-keyed: 'method'
end                                    % on BFI must not disturb the internal cycle
setApp(fig,app);
recomputeBase(fig);                    % new fingerprint for this type's files
applyInvalidation(fig,seed,type);      % forward cascade, restricted to this type
renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
% a setting can change the PRODUCT this type writes (contrastType -> _t_K | _s_K),
% and that flag is on every row label and checkbox - so redraw them, do not leave
% the panels showing the flag the type had before the edit
renderRawPanel(fig);
renderDerivedPanel(fig);
refreshConstructorSettings(fig);
refreshSummary(fig);
wbLog(fig,sprintf('edit [%s] %s.%s -> %s',type,stepId,field,val2str(value)));
autosaveSession(fig);
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
app.typeSel = wbTypeSelection('copy', app.typeSel, regFor(app,dst), src, dst);
app.sm      = wbSettingsModel('copyType', app.sm, src, dst);
setApp(fig,app);
recomputeBase(fig);
wbLog(fig,sprintf('constructor: copied configuration %s -> %s', src, dst));
setConstructorStatus(fig,sprintf('%s now matches %s.', dst, src));
renderConstructor(fig); renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
autosaveSession(fig);
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
reg = regFor(app,type);            % a protocol lands on the steps this type has
app.typeSel = wbTypeSelection('clear', app.typeSel, reg, type);
raw = wbTypeSelection('rawSteps', reg);
for i = 1:numel(p.steps)
    if ~any(strcmp(p.steps{i},raw)), continue; end
    b = branchOfStep(reg,p.steps{i});
    app.typeSel = wbTypeSelection('tick', app.typeSel, reg, type, b, p.steps{i});
end
brs = wbTypeSelection('rows', app.typeSel, reg, type);
for i = 1:numel(p.steps)
    if any(strcmp(p.steps{i},raw)), continue; end
    for b = 1:numel(brs)
        if ~wbTypeSelection('offers', reg, brs{b}, p.steps{i}), continue; end
        app.typeSel = wbTypeSelection('tick', app.typeSel, reg, type, brs{b}, p.steps{i});
    end
end
for i = 1:numel(p.animalSteps), app.animalSel(p.animalSteps{i}) = true; end
setApp(fig,app);
recomputeBase(fig);
wbLog(fig,sprintf('constructor: %s configured from the "%s" protocol (%s)', ...
    type, p.name, p.note));
setConstructorStatus(fig,sprintf('%s configured from "%s".', type, p.name));
renderConstructor(fig); renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
autosaveSession(fig);
end

function resetTypeConfig(fig,type)
%resetTypeConfig  Drop one type back to the launcher defaults (steps and settings).
app = getApp(fig);
type = char(type);
if isempty(type), setConstructorStatus(fig,'No type selected.'); return; end
app.typeSel = wbTypeSelection('clear', app.typeSel, regFor(app,type), type);
app.sm      = wbSettingsModel('resetType', app.sm, type);
setApp(fig,app);
recomputeBase(fig);
wbLog(fig,sprintf('constructor: reset %s to the launcher defaults', type));
setConstructorStatus(fig,sprintf('%s reset to the launcher defaults.', type));
renderConstructor(fig); renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
autosaveSession(fig);
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
    reg = regFor(app,ty);
    if typeIsMixed(app,ty)
        % legal but worth saying: the type is configured for the kind of recording
        % MOST of its files are, and the odd one out will show its steps as skipped
        lines{end+1} = sprintf(['%s: its recordings are not all the same kind - it is set up ' ...
            'for %s, and the others will be skipped'], ty, modalityOfType(app,ty)); %#ok<AGROW>
    end
    brs = wbTypeSelection('rows', app.typeSel, reg, ty);
    if isempty(brs)
        lines{end+1} = sprintf('%s: nothing selected yet - %d recording(s)', ty, n); %#ok<AGROW>
        continue
    end
    for b = 1:numel(brs)
        ids = wbTypeSelection('steps', app.typeSel, reg, ty, brs{b});
        inh = wbTypeSelection('inherited', app.typeSel, reg, ty, brs{b});
        lbl = sprintf('%s (%s)', ty, rowFlagFor(app,ty,brs{b}));
        txt = stepWords(reg,ids);
        if ~isempty(inh)
            txt = sprintf('%s + carried over: %s', txt, stepWords(reg,inh));
        end
        lines{end+1} = sprintf('%s: %s - %d %s on %d recordings = %d runs', ...
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
    lines{end+1} = sprintf('Once per animal: nothing selected - %d %s', nA, plural(nA,'animal'));
else
    lines{end+1} = sprintf('Once per animal: %s - %d %s on %d %s = %d runs', ...
        stepWords(app.reg,aids), numel(aids), plural(numel(aids),'step'), ...
        nA, plural(nA,'animal'), numel(aids)*nA);
end
end
function w = plural(n,word)
w = word; if n~=1, w = [word 's']; end
end
function s = stepWords(reg,ids)
%stepWords  A step list as the READER knows it: the registry's own labels, never
%   the internal ids.  'setRegions, dynamicSegmentation' is a programmer's list;
%   'Regions, Dynamic segmentation' is what the columns and the checkboxes say, so
%   the summary and the warnings say it too.
ids = reshape(ids,1,[]);
out = cell(1,numel(ids));
for i = 1:numel(ids)
    st = stepById(reg,ids{i});
    if isempty(st), out{i} = ids{i}; else, out{i} = st.label; end
end
s = strjoin(out,', ');
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
if isempty(ref), ref = {'no animals loaded yet'}; end
lines = [lines, {''}, {'Reference recordings:'}, ref];
uitextarea(g,'Value',lines,'Editable','off','FontName','monospaced', ...
    'Tooltip',['What each row will run, and each animal''s reference recording with ' ...
               'the exact file every step that needs one will read.']);
w = constructorWarnings(app);
lb = uilabel(g,'Text',warnLine(w),'WordWrap','on');
if isempty(w), lb.FontColor = [0 0.45 0]; else, lb.FontColor = [0.75 0.35 0]; end
end

%% ===================== PROCESS tab ================================== %%
function buildProcessTab(fig)
%buildProcessTab  The MONITOR (spec D6): choose how much to run, press Run, watch,
%   open results.  Nothing here is selectable - what runs was decided per TYPE on
%   the Constructor tab, and the From/To pair only says how much of THAT to execute
%   (spec D1) - so the picture is read-only uitables rather than a live widget per
%   cell.  200 files x 15 steps is 3000 widgets and was exactly where the old matrix
%   fell over; a table costs one component whatever the file count.
%
%   TWO TABLES, SHAPED LIKE THE CONSTRUCTOR (spec D6/D7).  The Constructor asks two
%   questions - which RAW steps a type runs, then which DERIVED steps each
%   (type, product) row runs - so the monitor answers in the same two halves:
%   above, one row per recording under the raw producers; below, one row per
%   (recording, product) under the derived steps.  With a type driving two pipelines
%   the single flat table could not say which column landed on which product, which
%   is exactly the ambiguity Phase 2b removed upstream.
%
%   THE LOG SITS BESIDE THE RAW TABLE, NOT UNDER EVERYTHING.  The raw table carries
%   only the handful of RAW producer columns, so it never fills the width it was
%   given, while the derived table carries every later column and wants all the
%   height it can get.  So the run log - the preview order, the progress and the
%   errors - moves into the space the raw table was wasting, and the products table
%   takes the full height that frees up.
app = getApp(fig); t = app.tabs.process;
gl = uigridlayout(t,[1 2],'ColumnWidth',{'2.6x','1x'},'Padding',[6 6 6 6],'ColumnSpacing',8);

% -- LEFT: toolbar + (recordings | log) + the products table + the selected row --
left = uigridlayout(gl,[4 1], ...
    'RowHeight',{'fit','1x','1.5x','fit'},'RowSpacing',4);
tb = uigridlayout(left,[1 13], ...
    'ColumnWidth',{'fit',130,'fit',130,'fit','fit','fit','fit','1x','fit','fit','fit','fit'}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',4);
uilabel(tb,'Text','From','FontWeight','bold','Tooltip',tipRangeFrom());
c.fromDrop = uidropdown(tb,'Items',{noRangeItem()},'ItemsData',{''}, ...
    'ValueChangedFcn',ui(fig,@(s,~)onRangeEdit(fig,'from',s.Value)),'Tooltip',tipRangeFrom());
uilabel(tb,'Text','To','FontWeight','bold','Tooltip',tipRangeTo());
c.toDrop = uidropdown(tb,'Items',{noRangeItem()},'ItemsData',{''}, ...
    'ValueChangedFcn',ui(fig,@(s,~)onRangeEdit(fig,'to',s.Value)),'Tooltip',tipRangeTo());
c.previewBtn = uibutton(tb,'Text','Preview order','ButtonPushedFcn',ui(fig,@(~,~)dryRun(fig)), ...
    'Tooltip','List everything that would run, in the order it would run.  Nothing is processed.');
c.runBtn = uibutton(tb,'Text','Run','ButtonPushedFcn',ui(fig,@(~,~)runChecked(fig)), ...
    'BackgroundColor',[0.82 0.92 0.82], ...
    'Tooltip','Run everything from the From step to the To step, in the right order.');
% Stop is NOT wrapped in ui() (spec D4): it is the one control that has to work
% BECAUSE the latch is held, and the run holds it from the first click to the last
% line of finishRun.  Dropping it would mean a run could never be stopped.
c.stopBtn = uibutton(tb,'Text','Stop','Enable','off','ButtonPushedFcn',@(~,~)cancelRun(fig), ...
    'BackgroundColor',[1 0.86 0.86],'Tooltip','Stop once the step now running has finished.');
c.reproCheck = uicheckbox(tb,'Text','Re-process finished files','Value',false, ...
    'ValueChangedFcn',ui(fig,@(s,~)onReprocessEdit(fig,s.Value)),'Tooltip',tipReprocess());
c.progLabel = uilabel(tb,'Text','Ready.','FontAngle','italic');
c.presetDrop = uidropdown(tb,'Items',{'(launcher defaults)'},'Value','(launcher defaults)', ...
    'ValueChangedFcn',ui(fig,@(s,~)onPresetPick(fig,s.Value)), ...
    'Tooltip','Start every parameter from a set you saved earlier.');
uibutton(tb,'Text','Save preset...','ButtonPushedFcn',ui(fig,@(~,~)uiSavePreset(fig)));
uibutton(tb,'Text','Load preset...','ButtonPushedFcn',ui(fig,@(~,~)uiLoadPreset(fig)));
c.wipeBtn = uibutton(tb,'Text','Wipe all...','BackgroundColor',[1 0.82 0.82], ...
    'ButtonPushedFcn',ui(fig,@(~,~)uiWipeAll(fig)),'Tooltip',tipWipeAll());

% -- the recordings table and the run log share one row (see the header) --
%    HALF EACH, and the row is a child of 'left' just as the results table is, so
%    the two of them plus the gap between come to exactly the results table's width.
topRow = uigridlayout(left,[1 2],'ColumnWidth',{'1x','1x'}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',8);
rawBox = uigridlayout(topRow,[2 1],'RowHeight',{'fit','1x'},'Padding',[0 0 0 0],'RowSpacing',2);
uilabel(rawBox,'Text','Recordings - what is done to the recording itself, and it goes first', ...
    'FontWeight','bold','Tooltip',progressLegend());
c.rawTable = uitable(rawBox,'Data',{},'ColumnName',{'Recording'},'RowName',{}, ...
    'ColumnEditable',false,'CellSelectionCallback',ui(fig,@(~,ev)onProgressSelect(fig,'raw',ev)), ...
    'Tooltip',progressLegend());
logBox = uigridlayout(topRow,[2 1],'RowHeight',{'fit','1x'},'Padding',[0 0 0 0],'RowSpacing',2);
uilabel(logBox,'Text','Run log - the preview order, the progress and the errors', ...
    'FontWeight','bold');
c.log = uitextarea(logBox,'Value',{'Ready.'},'Editable','off');

prodBox = uigridlayout(left,[2 1],'RowHeight',{'fit','1x'},'Padding',[0 0 0 0],'RowSpacing',2);
uilabel(prodBox,'Text','Results - one row for each result your setup creates', ...
    'FontWeight','bold','Tooltip',derivedLegend());
c.derivedTable = uitable(prodBox,'Data',{},'ColumnName',{'Recording and result'},'RowName',{}, ...
    'ColumnEditable',false,'CellSelectionCallback',ui(fig,@(~,ev)onProgressSelect(fig,'derived',ev)), ...
    'Tooltip',derivedLegend());
c.detail = uilabel(left,'Text','Select a row to see its file.','WordWrap','on','FontAngle','italic');

% -- RIGHT: the report links (spec D5), what to show, the PDF switches (D9), Exit --
right = uigridlayout(gl,[6 1],'RowHeight',{'fit','fit','fit','1x','fit','fit'},'RowSpacing',6);
head = uigridlayout(right,[1 2],'ColumnWidth',{'1x','fit'}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',4);
uilabel(head,'Text','Reports (newest first)','FontWeight','bold');
c.kindDrop = uidropdown(head,'Items',resultKindItems(),'ItemsData',resultKindIds(), ...
    'Value','all','ValueChangedFcn',ui(fig,@(s,~)setResultFilter(fig,s.Value)), ...
    'Tooltip',tipResultKind());
c.pdfCheck = uicheckbox(right,'Text','Create PDF reports','Value',false, ...
    'ValueChangedFcn',ui(fig,@(s,~)setPdfReports(fig,s.Value)),'Tooltip',tipPdfReports());
c.delCheck = uicheckbox(right,'Text','Delete the images once the PDF is written', ...
    'Value',false,'Enable','off', ...
    'ValueChangedFcn',ui(fig,@(s,~)setDeleteAfterPdf(fig,s.Value)),'Tooltip',tipDeleteAfterPdf());
c.resultList = uilistbox(right,'Items',{emptyResultItem()},'ItemsData',{''}, ...
    'Tooltip',['Every report the run has produced, newest first.  Double-click one to ' ...
               'open it in your image viewer, where you can zoom while processing ' ...
               'carries on.'], ...
    'DoubleClickedFcn',ui(fig,@(s,~)openArtifactViewer(fig,s.Value)));
c.openBtn = uibutton(right,'Text','Open selected', ...
    'ButtonPushedFcn',ui(fig,@(~,~)openArtifactViewer(fig,getApp(fig).c.process.resultList.Value)));
f = tabFooter(right,fig,'',[],'');        % the last tab: Exit only
c.exitBtn = f.exitBtn;

app.c.process = c; setApp(fig,app);
end

function t = progressLegend()
%progressLegend  What the state words in the table mean, in one tooltip.  It is a
%   PICTURE, not a control: what runs is decided on the Constructor tab.
t = ['This table shows you what is happening; you cannot change anything in it.  ' ...
     'Each cell reads ' char(183) ' when your setup does not run that step, queued when ' ...
     'it is waiting its turn - which includes waiting for a file an earlier step of ' ...
     'the same run will make - running while it works, done when the result is ' ...
     'saved, error when it failed, and skipped when it was asked for but nothing in ' ...
     'your setup will make what it needs.'];
end

function t = derivedLegend()
%derivedLegend  The bottom table's extra promise: its rows come from the
%   CONFIGURATION, not from disk (spec D7), so a product appears here the moment its
%   producer is ticked - which is the only way to see what is queued FOR it.
t = [progressLegend() '  Each row is one recording and one of its results, taken ' ...
     'from what you set up on the Constructor tab - so a row appears as soon as you ' ...
     'tick the step that will create it, before the file exists.'];
end

function s = noRangeItem()
%noRangeItem  What the two dropdowns read with nothing configured.
s = '(nothing set up yet)';
end
function t = tipRangeFrom()
t = ['Where the run STARTS.  "Last valid" carries on where you left off: it begins at ' ...
     'the first step that is not finished.  Choosing a step by name starts there ' ...
     'instead, however much is done before it.  Whether work that is already ' ...
     'finished gets done again is the tick beside Run, not this.'];
end
function t = tipRangeTo()
t = ['Where the run STOPS.  That step is included.  Only steps at or after From are ' ...
     'offered, and this moves up with From if you push From past it.'];
end
function t = tipReprocess()
t = ['Off, a run leaves alone whatever is already finished and only does what is ' ...
     'still missing.  On, every step between From and To runs again on every ' ...
     'recording, finished or not, and overwrites what is there.  That is how you ' ...
     'redo segmentation with new parameters, and how you re-open the vessel-type ' ...
     'painter on an animal you have already done.'];
end

function s = emptyResultItem()
%emptyResultItem  What the result list reads with nothing (yet) to show.
s = '(no reports yet)';
end
function items = resultKindItems()
items = {'All','Images','PDFs'};
end
function ids = resultKindIds()
%resultKindIds  The filter vocabulary.  It filters what is SHOWN and never what is
%   stored, so switching back to All brings every entry back (spec D10).
ids = {'all','images','pdfs'};
end
function t = tipResultKind()
t = ['What to show: everything, only the report images, or only the PDFs.  This ' ...
     'hides entries, it never throws them away - switch back to All and they are ' ...
     'all still there.'];
end
function t = tipWipeAll()
t = ['Delete every processed file belonging to the ' ...
     'recordings listed on the Files tab, so you can process the whole set again ' ...
     'from scratch.  Your raw recordings are never touched, nothing outside the ' ...
     'listed recordings is touched, and you are shown how many files it is and ' ...
     'asked before anything is removed.  It cannot be undone.'];
end
function t = tipPdfReports()
t = ['With this on, all the report images one step produces are collected into a ' ...
     'single PDF, one page per image, as soon as that step finishes its last ' ...
     'recording.  You page through one document instead of opening 60 JPGs.  The ' ...
     'PDFs go into a "workbenchReports" folder beside the folder you scanned, or ' ...
     'beside the images themselves if you added files by hand, and they appear in ' ...
     'this list.  Stopping a run still leaves you the PDFs it already finished.'];
end
function t = tipDeleteAfterPdf()
t = ['With this on, the report images that went into a PDF are deleted once that ' ...
     'PDF has been written, so a folder keeps one document per step instead of ' ...
     'hundreds of JPGs.  Off by default, and worth thinking about before you turn ' ...
     'it on: their links leave this list with them, the PDF becomes the only copy ' ...
     'you have, and anything that goes looking for the images later - this list ' ...
     'when you come back to the recordings, or your own scripts - will not find ' ...
     'them.  Only images that really made it into the document are removed, and ' ...
     'nothing is removed at all while PDFs are off.'];
end

%% ===================== hand-off to the standalone tools ============ %%
function handOffSession(fig, tool)
%handOffSession  THE hand-off contract (spec §5), shared by both Files-tab buttons.
%   Export and Explore are no longer tabs: they are separate programs, and the ONLY
%   thing that travels between them and this window is the session FILE.  So the
%   sequence is always the same - make the session on disk current, then open the
%   tool on that path:
%
%       autosaveSession(fig)  ->  guiExport(sessionPath) / guiExplore(sessionPath)
%
%   Nothing is captured on the way back.  The return value is deliberately dropped,
%   so this window holds no handle to the tool, never waits for it, and is entirely
%   unaffected when the user closes it - and the tool, for its part, has no idea the
%   workbench exists.  When there is no session to write yet (nothing loaded at all)
%   the tool is simply opened empty: both work perfectly well from their own file /
%   folder pickers, which is the point of un-hosting them.
autosaveSession(fig);
app = getApp(fig);
pth = app.sessionPath;
if isempty(pth) && ~isempty(app.files)
    % autosave needs a scanned root to invent a location; files added BY HAND have
    % none, so the hand-off writes one beside the first file rather than handing
    % over nothing.  This is the only reason the button ever writes on its own.
    saveSessionTo(fig, fullfile(fileparts(app.files(1).path), 'workbench_session.mat'));
    pth = getApp(fig).sessionPath;
end
try
    if isempty(pth) || ~isfile(pth)
        wbLog(fig, [tool ': no session on disk yet - opening it empty.']);
        feval(tool);
    else
        wbLog(fig, [tool ' <- session ' pth]);
        feval(tool, pth);
    end
catch ME
    wbLog(fig, [tool ' could not be opened: ' ME.message]);
end
end

%% ===================== tab-change refresh ========================== %%
function onTabSelected(fig)
%onTabSelected  Gate the move off the Files tab, then repaint the monitor when the
%   Process tab shows: its rows and its From/To range are derived from the
%   Constructor's configuration, so they have to be rebuilt after every edit.
if ~guardTabSwitch(fig), return; end
app = getApp(fig);
if ~isfield(app,'tg') || ~isgraphics(app.tg), return; end
tab = app.tg.SelectedTab;
if isequal(tab, app.tabs.process), renderProgress(fig); end   % rows + range follow the config
end

%% ===================== loaders ===================================== %%
function onSourceEdit(fig)
%onSourceEdit  Remember the root/results/glob boxes (session state, not just UI).
app = getApp(fig); c = app.c.files;
app.root        = strtrim(c.root.Value);
app.resultsRoot = strtrim(c.resultsRoot.Value);
app.glob        = strtrim(c.glob.Value);
setApp(fig,app);
end
function onPatternEdit(fig,axis,value)
%onPatternEdit  A label regexp changed -> re-derive every label (overrides win).
setPattern(fig,axis,value);
end
function setPattern(fig,axis,value)
%setPattern  Store one regexp and re-derive what it decides.  THE REFERENCE RULE IS
%   LIVE, exactly like the four label axes: it used to be applied only by Scan, so
%   typing it afterwards changed nothing at all and the reference ticks never moved.
%   It is now answered from the working set that is already loaded (autoRefsFor), so
%   the rule can be tuned and re-tuned without touching the disk.  A pattern that is
%   not a valid regexp must not throw out of an edit box: it is reported in the
%   curation status line and the previous references stand.
app = getApp(fig);
if ~isfield(app.patterns,axis), return; end
app.patterns.(axis) = strtrim(char(value));
if strcmp(axis,'ref')
    try
        app.autoRef = autoRefsFor(app, {app.files.path});
    catch ME
        setApp(fig,app); syncFilesControls(fig);
        setCurStatus(fig,['That reference pattern cannot be used yet: ' oneLine(ME.message)]);
        return
    end
    setCurStatus(fig,'');
end
setApp(fig,app);
syncFilesControls(fig);
rebuildWorkingSet(fig, {app.files.path}, true);
end
function s = oneLine(msg)
%oneLine  A multi-line MATLAB error message folded onto one status line.
s = strtrim(regexprep(char(msg),'\s+',' '));
end
function setSource(fig,root,glob,resultsRoot)
%setSource  Set the two folders and the glob.  THE RESULTS FOLDER FOLLOWS THE ROOT
%   UNTIL IT IS TOLD OTHERWISE: while the two are equal there is no mapping at all
%   (getResultsPath rule 1), which is the default and what every existing project
%   gets, so a user who never thinks about it never sees a difference.  Once they
%   point it somewhere else the two differ and it stops following - retyping the
%   root would otherwise silently throw their answer away.
app = getApp(fig);
followed = ~resultsApart(app);
if nargin>=2 && ~isempty(root)
    app.root = char(root);
    if followed, app.resultsRoot = app.root; end
end
if nargin>=3 && ~isempty(glob), app.glob = char(glob); end
if nargin>=4 && ~isempty(resultsRoot), app.resultsRoot = char(resultsRoot); end
setApp(fig,app); syncFilesControls(fig);
end

function uiBrowseRoot(fig)
app = getApp(fig);
d = uigetdir(defaultDir(app.root),'Pick the root folder to scan recursively');
if isequal(d,0), return; end
setSource(fig,d,'');
end
function uiBrowseResults(fig)
app = getApp(fig);
d = uigetdir(defaultDir(resultsDirOf(app)),'Pick the folder the results should go in');
if isequal(d,0), return; end
setSource(fig,'','',d);
end
function d = resultsDirOf(app)
%resultsDirOf  Where the results browser opens: the results folder if one is set,
%   otherwise the root, so the pick starts beside the recordings rather than at
%   whatever MATLAB's working folder happens to be.
if ~isempty(app.resultsRoot), d = app.resultsRoot; else, d = app.root; end
end
function tf = resultsApart(app)
%resultsApart  Do the results live somewhere other than the recordings?  This is
%   getResultsPath's rule 1 asked from the window's side, written down ONCE so the
%   two places that care - the results folder following the root, and the rename
%   confirmation saying what will not move - cannot disagree about whether anything
%   is being mapped at all.
tf = ~isempty(app.resultsRoot) && ~strcmp(app.resultsRoot, app.root);
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
%   which forces the matching file into column 1 of the grid.  Which recording that
%   makes the animal's default reference is decided by autoRefsFor, the same helper
%   the edit box uses - the rule has ONE definition, and a scan is not a special
%   case of it (a hand-pinned reference still wins over both).
app = getApp(fig);
p = app.patterns;
if isempty(p.ref)
    disc = wbDiscoverFiles('folder', app.root, app.glob, p.animal, p.type, p.expGroup, app.resultsRoot);
else
    disc = wbDiscoverFiles('structured', app.root, app.glob, p.animal, p.ref, p.type, p.expGroup, app.resultsRoot);
end
paths = gridPaths(disc);
app.autoRef = autoRefsFor(app, paths);
setApp(fig,app);
n = rebuildWorkingSet(fig, paths, false);
end
function uiAddFiles(fig)
app = getApp(fig);
[f,p] = uigetfile({recordingFilter(),'Recordings & products';'*.*','All files'}, ...
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
    app.patterns.type, app.patterns.expGroup, app.resultsRoot);
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
paths = gridPaths(disc);
app.autoRef      = autoRefsFor(app, paths);
app.animalRefMan = containers.Map('KeyType','char','ValueType','char');
setApp(fig,app); syncFilesControls(fig);
rebuildWorkingSet(fig, paths, false);
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

function m = autoRefsFor(app, paths)
%autoRefsFor  THE reference rule: animal -> reference recording IDENTITY, answered
%   from a path list and the CURRENT Reference regexp.  No disk scan, and no
%   dependence on which loader produced the paths - that dependence is exactly what
%   made the rule inert unless a Scan had just run.
%
%   The regexp is matched against the bare file NAME WITH ITS EXTENSION, which is
%   what getFileNamesList matches in its own reference mode, so a pattern written
%   for a scan keeps meaning the same thing when it is retyped afterwards.  Where
%   several files of an animal match, the FIRST in working-set order wins, and what
%   is stored is that file's recording IDENTITY - never a branch file, since each
%   step resolves the branch it needs at run time.  An empty pattern is no rule at
%   all: every automatic reference goes, and the hand-pinned ones stand.
%
%   MATLAB's regexp does not reject a malformed pattern - it simply matches nothing
%   - so a half-typed rule costs the automatic references and no more.  It can still
%   throw (wbTypeModel converts a pattern its own axes choke on into an error), and
%   the caller is expected to catch that rather than let it out of an edit box.
m = containers.Map('KeyType','char','ValueType','char');
rx = strtrim(char(app.patterns.ref));
paths = cleanPathList(paths);
if isempty(rx) || isempty(paths), return; end

labels = wbTypeModel('applyOverrides', ...
    wbTypeModel('derive', paths, app.patterns), app.overrides);
for i = 1:numel(paths)
    mdl = wbFileModel(paths{i});
    if isempty(regexp(mdl.name, rx, 'once')), continue; end
    a = labels.animal{i};
    if isKey(m,a), continue; end                    % first match of the animal wins
    m(a) = mdl.identity;                            % identity: never a branch flag
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
disc = wbDiscoverFiles('curated', paths, app.labels, app.animalRef, app.root, app.resultsRoot);
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
    if isKey(byPath,p), m = byPath(p); else, m = wbFileModel(p, app.root, app.resultsRoot); end
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
%adoptDiscovery  Adopt a discovery grid: flatten to rows, set the modalities, render.
%
%   A SESSION IS NOT ONE MODALITY.  It used to be: the working set's modalities
%   were reduced to their statistical MODE and the whole registry filtered by it,
%   so a folder holding twelve speckle recordings and four videos silently lost
%   every video step - and one holding more videos than recordings lost the lot.
%   app.reg is now the UNION over the modalities actually present, in registry
%   order, and the two questions that are NOT session-wide are asked one scope
%   down: per FILE by wbStateEngine (a step is applicable when the file's own
%   modality exposes it) and per TYPE by regFor, which is what keeps a row's
%   producer, its branch and its requiresAny resolving to exactly one answer.
%   app.modality survives as the one-word summary in the Files-tab status line.
app = getApp(fig);
app.disc = disc;
[app.rows, app.animalNames, app.modelArr] = flattenDisc(disc);
mods = {};
if ~isempty(app.modelArr), mods = uniqueStable({app.modelArr.modality}); end
mods = mods(~cellfun(@isempty,mods));
if isempty(mods), mods = {'LSCI'}; end
app.modality = modeStr({app.modelArr.modality});
app.reg      = wbStepRegistry(mods);
app.regCache = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(mods), app.regCache(mods{i}) = wbStepRegistry(mods{i}); end
app.typeMod  = typeModalities(app);
if nargin<3 || ~keepOverlay
    % a fresh load starts with no invalidation overlay
    app.stale = containers.Map('KeyType','char','ValueType','any');
end
setApp(fig,app);
recomputeBase(fig);
refreshFileTable(fig);
renderProgress(fig);
logPerFileFlips(fig);         % D6: say so when a cell stopped reading as done
renderConstructor(fig);       % types are data: a label edit adds/removes a row live
refreshStatus(fig);
autosaveSession(fig);         % the state changed -> the last session follows it
end

function m = typeModalities(app)
%typeModalities  Each recording TYPE's modality: the modality of its files.
%   A type is a group of recordings processed the same way, so in practice its
%   files are one modality.  When they are not, the MODE stands for the type (and
%   the selection summary says so) - the Constructor has to configure the type
%   somehow, and the majority is the only defensible reading.  Nothing is stored
%   for a type with no files: regFor then falls back to the union.
m = containers.Map('KeyType','char','ValueType','char');
if isempty(app.files), return; end
tys = uniqueStable({app.files.type});
for i = 1:numel(tys)
    if isempty(tys{i}), continue; end
    mods = {app.files(strcmp({app.files.type}, tys{i})).modality};
    mods = mods(~cellfun(@isempty,mods));
    if isempty(mods), continue; end
    m(tys{i}) = modeStr(mods);
end
end

function reg = regFor(app,type)
%regFor  THE REGISTRY ONE TYPE IS CONFIGURED WITH - the union filtered down to that
%   type's own modality (wbStepRegistry, which also prunes the requirements naming
%   steps the filter removed).
%
%   WHY PER TYPE AND NOT PER SESSION.  Every Constructor question is asked of a
%   type: which rows it has, which step produces each row, which steps a row
%   offers, what a tick pulls in.  Each of those has exactly one answer only within
%   one modality - two modalities can name the same branch, and a step serving both
%   lists a producer for each in requiresAny - so asking them of the union would
%   resolve a row to another modality's producer.  The union is still right for
%   everything that spans the session: the monitor's columns, the run range, the
%   settings panel, and the per-file gating.
reg = app.reg;
ty = char(type);
if isempty(ty) || ~isfield(app,'typeMod') || ~isKey(app.typeMod,ty), return; end
md = app.typeMod(ty);
if isfield(app,'regCache') && isKey(app.regCache,md), reg = app.regCache(md);
else,                                                 reg = wbStepRegistry(md);
end
end

function md = modalityOfType(app,type)
%modalityOfType  What regFor filtered by ('' when the type has no files).
md = '';
ty = char(type);
if ~isempty(ty) && isfield(app,'typeMod') && isKey(app.typeMod,ty), md = app.typeMod(ty); end
end

function tf = typeIsMixed(app,type)
%typeIsMixed  Does this type hold more than one modality?  Legal, and reported in
%   the selection summary, because the type is then configured for its majority.
tf = false;
if isempty(app.files), return; end
mods = {app.files(strcmp({app.files.type}, char(type))).modality};
tf = numel(unique(mods(~cellfun(@isempty,mods)))) > 1;
end

function tf = stepInReg(reg,stepId)
%stepInReg  Is this column one of the steps that type's recordings can run?
tf = ~isempty(reg) && any(strcmp(char(stepId), {reg.id}));
end

function t = wrongModalityTip(app,type,step)
%wrongModalityTip  Why a column is dead for this type - the file's own kind, in the
%   words the Files tab uses for it.
md = modalityOfType(app,type);
if isempty(md), md = 'these'; end
t = sprintf('%s does not run on %s recordings, and %s is %s.', ...
    step.label, strjoin(step.modalities,'/'), type, md);
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
%   Reading the FULL Data, keyed by the 'file' + 'file type' columns, instead of the
%   edited index keeps this correct whatever the table's sort order is - a sorted
%   uitable renumbers what the user sees, not its Data.  The key works because file
%   names are unique across a scanned tree; a name that is NOT unique is reported by
%   fileProblems and its rows are left alone here.
%
%   A NAME EDIT IS HANDLED FIRST AND ENDS THE CALLBACK, because it moved the very
%   key the reconciliation below reads by: every other row would still match, and
%   the edited one would silently look like an unknown file.
D = src.Data;
if isempty(D), return; end
if handleNameEdit(fig, D), return; end
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

%% ---- renaming a recording (author decision D4) --------------------- %%
function tf = handleNameEdit(fig, D)
%handleNameEdit  Spot a rename typed into the 'file' column and carry it out.
%   FOUND BY VALUE, NEVER BY ROW INDEX: the table is sortable, so the edited row's
%   position says nothing about which file it is.  What a name edit does say is that
%   the table's set of names and the working set's differ in exactly one place - one
%   name has gone and one has appeared - and that pair IS the rename.  Anything else
%   (both sets equal) is a label edit and is left to the reconciliation loop.
%   Returns true when the edit was a name edit, refusals included, so the caller
%   stops rather than reconciling rows against a key that has just moved.
tf = false;
app = getApp(fig);
if isempty(app.files), return; end
col = columnIndex();

shown = cell(1,size(D,1));
for i = 1:size(D,1)
    shown{i} = [charOf(D{i,col.file}) charOf(D{i,col.filetype})];
end
gone = setdiff({app.files.name}, shown);
came = setdiff(shown, {app.files.name});
if isempty(gone) && isempty(came), return; end
tf = true;

if isscalar(gone) && isempty(came)
    setCurStatus(fig,'That name is already used by another file in the list - nothing was renamed.');
    refreshFileTable(fig); return
end
if ~isscalar(gone) || ~isscalar(came)
    setCurStatus(fig,'The file names could not be matched up - nothing was renamed.');
    refreshFileTable(fig); return
end

k = find(strcmp({app.files.name}, gone{1}), 1);
if isempty(k), refreshFileTable(fig); return; end
newBare = '';
for i = 1:size(D,1)
    if strcmp([charOf(D{i,col.file}) charOf(D{i,col.filetype})], came{1})
        newBare = charOf(D{i,col.file}); break
    end
end
renameFromTable(fig, app.files(k).path, newBare);
end

function renameFromTable(fig, oldPath, newBare)
%renameFromTable  Reduce an edited table name to a STEM, then rename the recording.
%   The 'file' column shows the whole bare name - the crop prefix, the stem, the
%   stage/product flags and the role, all at once - but only the stem is the
%   recording's name.  The prefix and the tail belong to the naming GRAMMAR
%   (wbFileModel), and a rename that changed them would rename one product out of
%   its own set of siblings, so both are required to survive the edit.
m = wbFileModel(oldPath);
[~, oldBare] = fileparts(oldPath);
pre  = m.roiPrefix;
tail = oldBare(numel(pre)+numel(m.stem)+1 : end);       % '_t_K_d', '' for a raw file

if ~strncmp(newBare, pre, numel(pre))
    setCurStatus(fig, sprintf(['Keep the "%s" at the front - it is part of the ' ...
        'recording name.  Nothing was renamed.'], pre));
    refreshFileTable(fig); return
end
if ~isempty(tail) && ~endsWith(newBare, tail)
    setCurStatus(fig, sprintf(['Keep the "%s" at the end - you are renaming the ' ...
        'recording, not one of its result files.  Nothing was renamed.'], tail));
    refreshFileTable(fig); return
end
newStem = newBare(numel(pre)+1 : numel(newBare)-numel(tail));
renameRecording(fig, oldPath, newStem);
end

function ok = renameRecording(fig, path, newStem)
%renameRecording  Rename a WHOLE recording on disk (author decision D4): every file
%   named after it moves together - the _d/_r/_s set of every branch, the report
%   images, the raw recording and any workbook written beside them.  Renaming one
%   member is not on offer: every wrapper finds its siblings by name, so a single
%   file renamed out of the set breaks the next step that looks for them.
%
%   The file work is wbRename's (pure, headless, testable); what lives here is the
%   half only the window knows - the refusals, the confirmation that lists every
%   file before any of them moves, and moving the session state onto the new names.
ok = false;
app = getApp(fig);
if isfield(app,'running') && app.running
    setCurStatus(fig,'A run is in progress - stop it before renaming a recording.');
    refreshFileTable(fig); return
end
k = find(strcmp({app.files.path}, char(path)), 1);
if isempty(k)
    setCurStatus(fig,'That file is not in the list - nothing was renamed.');
    refreshFileTable(fig); return
end
model   = app.files(k).model;
newStem = strtrim(char(newStem));
if strcmp(newStem, model.stem), refreshFileTable(fig); return; end   % nothing typed

list = wbRename('plan', model, newStem);
why  = renameProblem(app, model, newStem, list);
if isempty(why)
    [fine, whyNot] = wbRename('check', list);
    if ~fine, why = whyNot; end
end
if ~isempty(why)
    setCurStatus(fig,['Nothing was renamed - ' why]);
    refreshFileTable(fig); return
end
if ~confirmRename(fig, model, newStem, list)
    setCurStatus(fig,'Rename cancelled - nothing was moved.');
    refreshFileTable(fig); return
end

out = wbRename('apply', list);
if ~out.ok
    setCurStatus(fig,['The rename failed - ' out.why]);
    wbLog(fig,['rename failed: ' out.why ternary(out.rolledBack, ...
        '  (the files that had already moved were put back)', ...
        '  (SOME FILES MAY HAVE MOVED - check the folder)')]);
    refreshFileTable(fig); return
end
adoptRename(fig, list, model, newStem);
setCurStatus(fig, sprintf('Renamed %d %s to %s.', numel(list), ...
    plural(numel(list),'file'), [model.roiPrefix newStem]));
wbLog(fig, sprintf('renamed %s -> %s (%d file(s))', [model.roiPrefix model.stem], ...
    [model.roiPrefix newStem], numel(list)));
ok = true;
end

function why = renameProblem(app, model, newStem, list)
%renameProblem  Why the WORKING SET cannot accept this rename ('' when it can).
%   Only what the LIST knows is judged here.  Whether a name is usable at all, and
%   whether anything on disk is in the way, is wbRename's call on the names it
%   composes, so the character rule has one definition and lives with the code that
%   writes to disk.
%
%   THE DUPLICATE-NAME RULE IS ENFORCED BEFORE ANYTHING MOVES.  fileProblems already
%   refuses to let a set with two identically named files go any further - the
%   workbench and the pipeline both identify a recording by name - but by the time
%   it spoke the files would already have been renamed.  So a rename that would
%   create the clash is simply not offered.
why = '';
newStem = char(newStem);
if isempty(newStem), why = 'type a name for the recording.'; return; end
newId  = fullfile(model.folder,[model.roiPrefix newStem]);
lands  = {list.newName};
for i = 1:numel(app.files)
    f = app.files(i);
    if strcmp(f.model.identity, model.identity), continue; end   % the recording itself
    if strcmpi(f.model.identity, newId)
        why = sprintf('"%s" is already another recording in the list.', ...
            [model.roiPrefix newStem]); return
    end
    if any(strcmpi(f.name, lands))
        why = sprintf(['another file in the list is already called "%s", and no two ' ...
            'may share a name.'], f.name); return
    end
end
end

function tf = confirmRename(fig, model, newStem, list)
%confirmRename  Show EVERY file that will move, before any of them does (D4).  A
%   headless window (the tests, an API caller) has nobody to ask and proceeds.
%
%   IT IS ALSO THE ONLY WARNING THERE IS when the results are kept in a folder of
%   their own.  A rename then covers the recording's own folder and nothing else
%   (wbRename), so the closing sentence has to say what stays behind rather than
%   what travels - the two are opposite answers and only one of them can be shown.
%   There is no second dialog and no banner afterwards: this is where it is said.
tf = true;
if ~isvalid(fig) || ~strcmp(fig.Visible,'on'), return; end
lines = cell(1,numel(list));
for i = 1:numel(list)
    lines{i} = ['    ' list(i).oldName '   ->   ' list(i).newName];
end
if resultsApart(getApp(fig))
    note = ['Results and reports already computed for this recording keep their old ' ...
        'name and will no longer be linked to it.'];
else
    note = ['A recording is renamed as a whole - everything computed from it, and its ' ...
        'report images, move with it, because the processing steps find them by name.'];
end
msg = sprintf('Rename the recording %s to %s?\n\n%d %s move together:\n\n%s\n\n%s', ...
    [model.roiPrefix model.stem], [model.roiPrefix newStem], ...
    numel(list), plural(numel(list),'file'), strjoin(lines,newline), note);
sel = uiconfirm(fig,msg,'Rename recording','Options',{'Rename','Cancel'}, ...
    'DefaultOption',2,'CancelOption',2,'Icon','question');
tf = strcmp(sel,'Rename');
end

function adoptRename(fig, list, model, newStem)
%adoptRename  Move every piece of state that is keyed by a PATH or by a recording
%   IDENTITY onto the new names, then rebuild the working set from them.
%
%   THE POINT OF DOING IT IN ONE PLACE.  Nine containers are keyed one of those two
%   ways, and a rename that missed one would leave a label override, a completion
%   record or a per-animal reference pointing at a file that no longer exists - a
%   loss that would only surface on the next run, or after the next session load.
%   So there is one dictionary per kind of key and one two-line helper (remapMap),
%   and every container is routed through it rather than getting a loop of its own.
app = getApp(fig);
pMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(list), pMap(list(i).old) = list(i).new; end
iMap = containers.Map('KeyType','char','ValueType','char');
iMap(model.identity) = fullfile(model.folder,[model.roiPrefix newStem]);

% ---- keyed by FILE PATH -----------------------------------------------------
ax = fieldnames(app.overrides);
for i = 1:numel(ax)
    app.overrides.(ax{i}) = remapMap(app.overrides.(ax{i}), pMap, 'key');
end
app.modalityOvr = remapMap(app.modalityOvr, pMap, 'key');
app.completed   = remapMap(app.completed,   pMap, 'head');   % 'path||stepId'

% ---- keyed by recording IDENTITY -------------------------------------------
app.stale            = remapMap(app.stale,            iMap, 'head');  % 'identity||stepId'
app.runState         = remapMap(app.runState,         iMap, 'head');
app.cellMsg          = remapMap(app.cellMsg,          iMap, 'head');
app.sm.fileOverrides = remapMap(app.sm.fileOverrides, iMap, 'head');  % settings model
% .base / .fileState / .branchState are re-derived from disk by recomputeBase a few
% lines below; they are moved anyway so the app struct is never briefly inconsistent
app.base        = remapMap(app.base,        iMap, 'key');
app.fileState   = remapMap(app.fileState,   pMap, 'key');
app.branchState = remapMap(app.branchState, iMap, 'head');   % 'identity||branch'

% ---- the per-animal references hold an identity as their VALUE --------------
app.animalRef    = remapMap(app.animalRef,    iMap, 'value');
app.autoRef      = remapMap(app.autoRef,      iMap, 'value');
app.animalRefMan = remapMap(app.animalRefMan, iMap, 'value');

% ---- the report links this session collected point at files that moved ------
for i = 1:numel(app.results)
    app.results(i).path = remapKey(app.results(i).path, pMap);
end

paths = cellfun(@(p) remapKey(p,pMap), {app.files.path}, 'UniformOutput', false);
setApp(fig,app);
rebuildWorkingSet(fig, paths, true);   % relabels, re-refs, recomputes, repaints, autosaves
end

function m = remapMap(m, dict, mode)
%remapMap  Move one containers.Map onto renamed keys or values.  THE single place a
%   rename touches keyed state:
%     'key'   the whole key is a path or an identity;
%     'head'  the key is '<path-or-identity>||<something>' and only the head moves;
%     'value' the key is something else and the VALUE is an identity.
if ~isa(m,'containers.Map') || m.Count==0 || dict.Count==0, return; end
k = keys(m); v = values(m);
out = containers.Map('KeyType','char','ValueType',m.ValueType);
for i = 1:numel(k)
    key = k{i}; val = v{i};
    switch mode
        case 'key',   key = remapKey(key, dict);
        case 'head',  key = remapHead(key, dict);
        case 'value', val = remapKey(val, dict);
    end
    out(key) = val;
end
m = out;
end
function s = remapKey(s, dict)
%remapKey  One path or identity, moved when the rename touched it.
s = char(s);
if isKey(dict,s), s = dict(s); end
end
function s = remapHead(s, dict)
%remapHead  A composite '<path-or-identity>||<rest>' key - only its head can move.
s = char(s);
j = strfind(s,'||');
if isempty(j), s = remapKey(s,dict); return; end
s = [remapKey(s(1:j(1)-1),dict) s(j(1):end)];
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
msg = sprintf('%d of %d selected rows set to %s "%s".', n, numel(sel), field, v);
if refused > 0
    msg = [msg sprintf('  %d were left alone - their file format does not allow it.', refused)];
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
    problems{end+1} = sprintf(['%d file names appear more than once: %s. ' ...
        'Recordings are identified by name, here and in the pipeline, so no two ' ...
        'may share one - rename them, or search a narrower folder.'], ...
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
        problems{end+1} = sprintf('%d files have no %s assigned: %s.', ...
            nnz(bad), nm{i}, strjoin(shortList({app.files(bad).name}),', ')); %#ok<AGROW>
    end
end

noRef = {};
animalsNow = wbTypeModel('values', app.labels, 'animal');
for i = 1:numel(animalsNow)
    if ~isKey(app.animalRef, animalsNow{i}), noRef{end+1} = animalsNow{i}; end %#ok<AGROW>
end
if ~isempty(noRef)
    warnings{end+1} = sprintf(['%d animals have no reference recording: %s. ' ...
        'That is allowed, but registration and vessel typing need one.'], ...
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

function s = warnText(w)
if isempty(w), s = ''; else, s = ['   ' strjoin(w,'   ')]; end
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
%guardTabSwitchTo  Move to a tab, or refuse and say what is still missing.  This is
%   what every 'Next' button does, and setting SelectedTab in code does NOT fire the
%   tab group's callback - so the monitor is repainted here as well, or arriving by
%   Next would show a Process tab still holding the configuration's previous shape.
tf = false;
app = getApp(fig);
if ~isfield(app,'tg') || ~isgraphics(app.tg), return; end
if ~filesValid(app), alertIncomplete(fig,app); return; end
app.tg.SelectedTab = tab;
tf = true;
if isequal(tab, app.tabs.process), renderProgress(fig); end
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
%   matches more than one branch (*_K_r.mat) would otherwise emit duplicate rows.
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
%syncFilesControls  Push root/results/glob/patterns back into their boxes (session load).
app = getApp(fig);
if ~isfield(app,'c') || ~isfield(app.c,'files'), return; end
c = app.c.files;
if isgraphics(c.root), c.root.Value = app.root; end
if isgraphics(c.resultsRoot), c.resultsRoot.Value = app.resultsRoot; end
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
    nTy = numel(wbTypeModel('values', app.labels, 'type'));
    nGr = numel(wbTypeModel('values', app.labels, 'expGroup'));
    mods = uniqueStable({app.files.modality});
    mods = mods(~cellfun(@isempty,mods));
    if isempty(mods), mods = {app.modality}; end
    txt = sprintf(['%d files - %d %s - %d recording %s - %d experimental %s - ' ...
        '%d untyped - %d %s without a reference recording.  %s %s.'], ...
        numel(app.files), numel(animalsNow), plural(numel(animalsNow),'animal'), ...
        nTy, plural(nTy,'type'), nGr, plural(nGr,'group'), ...
        nUntyped, nNoRef, plural(nNoRef,'animal'), ...
        ternary(isscalar(mods),'Modality','Modalities'), strjoin(mods,', '));
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
%
%   IT IS RE-RUN WHEN THE CONFIGURATION CHANGES TOO, not only when the disk does:
%   since the look-ahead landed (statesOf), what a type is TICKED for is one of the
%   inputs, so a Constructor tick makes these maps stale exactly as a settings edit
%   does - and both now call this.
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
        bm = branchModelOf(app, app.rows(i).model, brs{b});
        if isempty(bm), continue; end                % that pipeline has no file yet
        app.branchState([id '||' brs{b}]) = statesOf(app, bm, ty);
    end
end
setApp(fig,app);
end

function bs = statesOf(app, model, type)
%statesOf  wbStateEngine's answer for one model, folded into a stepId->state struct.
%
%   THE DISK IS NOT THE ONLY SOURCE (author, 2026-07-31).  A step is also satisfied
%   when the step that will produce its input is selected to run in this sequence -
%   the run order is step-major, so the producer is genuinely there by the time the
%   consumer is reached.  What that configuration will produce is plannedStepIds,
%   the same answer the Constructor and the monitor read; wbStateEngine folds it
%   into its prerequisite test and reports such a cell as READY (never done), which
%   is why only .state is kept here and the distinction lives in .reason.
if nargin<3, type = modelType(model); end
cs = curSettingsFor(app, model, type);
opts = struct('plannedIds', {plannedStepIds(app, type)});
st = wbStateEngine(model, app.reg, cs, opts);
bs = struct();
for k = 1:numel(st), bs.(st(k).id) = st(k).state; end
end

function m = branchModelOf(app, model, branch)
%branchModelOf  ANY data file of one recording sitting on one pipeline ([] if the
%   pipeline has produced nothing yet).  Which file does not matter: wbStateEngine
%   unions the settings ALONG the pipeline, so '_t_K_r' and '_t_BFI_r' give the
%   same answer - and asking for a fixed one would break the moment a step deleted
%   its original (runBFI's deleteOriginal).
%
%   IT TAKES THE ROW'S MODEL, NOT ITS IDENTITY.  The files it is looking for are
%   PRODUCTS, so they are in the results folder, and an identity only ever names
%   the tree the recording was scanned from.  The one it hands back is given the
%   two project folders for the same reason: a branch model with no raw folder
%   sends every needsRaw step of that pipeline looking for the recording among the
%   results.
m = [];
if isempty(model.resultsFolder) || ~isfolder(model.resultsFolder), return; end
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.resultsFolder,[base '*_r.mat']));
for i = 1:numel(d)
    cm = wbFileModel(fullfile(d(i).folder, d(i).name), app.root, app.resultsRoot);
    if strcmp([cm.roiPrefix cm.stem], base) && strcmp(cm.branch, branch), m = cm; return; end
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
%resolveCellState  What one RECORDING's step reads as: the disk baseline with the
%   two live overlays on top - the transient run state, and the in-session
%   invalidation overlay (app.stale).
%
%   THERE IS NO 'checked' STATE ANY MORE (Phase 6).  Selection is a property of the
%   recording TYPE and lives in the Constructor, and the per-cell queue that used to
%   promote 'ready' to 'checked' went with it.
%   THE LOOK-AHEAD IS BACK, BUT ONE LAYER DOWN (2026-07-31).  Phase 6 also dropped
%   the projection that promoted a cell whose prerequisites were merely queued, and
%   that turned out to be a step too far: it made the workbench say a file was
%   missing while the same run was about to write it.  It now lives inside
%   wbStateEngine, which is handed the configuration's planned step ids (statesOf)
%   and answers 'ready' with a reason - so readiness is still exactly what
%   wbStateEngine says it is, and nothing here re-derives it.
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
end

%% ===================== the progress tables (read-only) ============== %%
function rows = progressRows(app, kind)
%progressRows  A monitor table's rows, ordered the way the run goes - animal,
%   reference first, then the file's position in its animal - so watching reads top
%   to bottom.  The two tables mirror the Constructor's two questions (spec D6):
%
%     'raw'     ONE PER FILE of the working set: what happens to the recording
%               itself.  A recording with two product files in the working set
%               contributes two rows, each gated on its OWN pipeline (D8).
%     'derived' ONE PER (recording, PRODUCT), taken from the CONFIGURATION rather
%               than from disk (D7): the branches this file's type actually runs,
%               so a '_c' row appears the moment its raw producer is ticked, long
%               before the '_c' file exists.  Deduped on (identity, branch), since
%               several files of one recording answer to the same product row.
if nargin<2, kind = 'raw'; end
rows = emptyProgressRows();
brsAll = wbTypeSelection('branches', app.reg);
seen = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(app.files)
    f = app.files(i);
    k = find(strcmp({app.rows.identity}, f.model.identity),1);
    if isempty(k), ai = inf; ria = inf; isRef = f.isRef; else
        ai = app.rows(k).animalIdx; ria = app.rows(k).rowInAnimal; isRef = app.rows(k).isRef;
    end
    if strcmp(kind,'raw')
        rows(end+1) = mkProgressRow(app,'raw',f,f.model.branch,ai,ria,isRef,0); %#ok<AGROW>
        continue
    end
    if isempty(f.type), continue; end
    brs = wbTypeSelection('rows', app.typeSel, regFor(app,f.type), f.type);
    for b = 1:numel(brs)
        key = [f.model.identity '||' brs{b}];
        if isKey(seen,key), continue; end              % one row per (recording, product)
        seen(key) = true;
        bi = find(strcmp(brs{b}, brsAll),1);
        if isempty(bi), bi = numel(brsAll)+1; end
        rows(end+1) = mkProgressRow(app,'derived',f,brs{b},ai,ria,isRef,bi); %#ok<AGROW>
    end
end
if isempty(rows), return; end
[~,ord] = sortrows([[rows.animalIdx]', [rows.rowInAnimal]', [rows.branchIdx]', (1:numel(rows))']);
rows = rows(ord);
end

function r = mkProgressRow(app,kind,f,branch,ai,ria,isRef,branchIdx)
%mkProgressRow  One row of either table: the FILE NAME WITHOUT ITS EXTENSION, then
%   the product this row stands for in brackets - the same flag rowFlagFor gives the
%   Constructor (spec D8), so the two surfaces name a product identically.  The
%   extension is the one part of the name that never distinguishes two rows here,
%   so it only costs column width.  The reference recording says so in words.
flag = '';
if strcmp(kind,'derived') && ~isempty(f.type), flag = rowFlagFor(app,f.type,branch); end
[~,lbl] = fileparts(char(f.name));
if ~isempty(flag), lbl = sprintf('%s (%s)', lbl, flag); end
if isRef, lbl = [lbl '  reference']; end
r = struct('kind',kind,'path',f.path,'model',f.model,'identity',f.model.identity, ...
    'label',f.name,'rowLabel',lbl,'animal',f.animal,'type',f.type,'expGroup',f.expGroup, ...
    'branch',char(branch),'flag',flag,'isRef',logical(isRef), ...
    'animalIdx',ai,'rowInAnimal',ria,'branchIdx',branchIdx);
end

function ids = monitorColumns(reg, kind)
%monitorColumns  A table's columns, DERIVED from the registry the same way the
%   Constructor derives its two panels - nothing here lists a step.  The per-animal
%   steps ride with the derived half: they land ON a product row (registration
%   templates the '_t' product, vessel typing paints the '_c' one) and
%   plannedStepsFor already adds them to every row, so between the two tables every
%   registry step has exactly one column.
if strcmp(kind,'raw')
    ids = wbTypeSelection('rawSteps', reg);
else
    ids = orderIds(reg, [wbTypeSelection('derivedSteps',reg), wbTypeSelection('animalSteps',reg)]);
end
end

function renderProgress(fig)
%renderProgress  (Re)build BOTH state tables from scratch: rows, columns, contents,
%   and the From/To items - all of which follow the configuration, so this is what
%   a Constructor edit calls.  refreshCells is the cheap repaint that does not.
app = getApp(fig); c = app.c.process;
if ~isfield(c,'rawTable') || ~isgraphics(c.rawTable), return; end
app.pRows  = progressRows(app,'raw');
app.dRows  = progressRows(app,'derived');
app.pRowOf = containers.Map('KeyType','char','ValueType','double');
app.dRowOf = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(app.pRows), app.pRowOf(app.pRows(i).path) = i; end
for i = 1:numel(app.dRows)
    app.dRowOf([app.dRows(i).identity '||' app.dRows(i).branch]) = i;
end
setApp(fig,app);

setTableColumns(c.rawTable,     'Recording',           monitorColumns(app.reg,'raw'),     app.reg);
setTableColumns(c.derivedTable, 'Recording and result', monitorColumns(app.reg,'derived'), app.reg);
refreshCells(fig,{});
refreshRangeControls(fig);
refreshProgressDetail(fig);
end

function setTableColumns(h, firstName, ids, reg)
%setTableColumns  Column headers of one monitor table, from the registry's labels,
%   and the width policy both of them run under (spec D5).
%
%   THE WIDTHS ARE THE ENGINE'S JOB, NOT OURS.  R2025b's ColumnWidth takes 'fit'
%   (size to this column's own content, header included) and '1x' (take what is
%   left, the way a uigridlayout's '1x' or a WPF star does).  Both were fixed
%   numbers here - 240 for the name and 94 for every step - and both cropped:
%   'Internal cycle' alone measures 95 px and 'Dynamic segmentation' 163, while a
%   row label is a file name plus a product flag plus, sometimes, the word
%   reference.  An earlier attempt estimated the widths in MATLAB instead
%   (px-per-character, measured back from the table's own contents); it was dropped
%   because the estimate is necessarily coarser than the metrics the renderer
%   already has - it wanted 1725 px for the results table where 'fit' needs about
%   1100 - and because it needed a resize hook the layout engine does for free.
%
%   EVERY REAL COLUMN 'fit', PLUS AN EMPTY TRAILING COLUMN THAT CARRIES THE '1x'.
%   The spacer is the whole trick, and it is there because a '1x' has NO MINIMUM: it
%   absorbs the shortfall by shrinking, so whichever real column carries it gets
%   clipped instead of the table scrolling.  Measured, on the recordings table:
%   the star on the NAME column clips the file name, and moving it to the last
%   column just clips 'Internal cycle' down to 'I..' - in both cases with no
%   scrollbar, which is precisely what the fixed widths were being replaced for.
%   An empty column has nothing to clip, so:
%     room to spare  - the spacer takes the leftover and the table covers its block,
%                      header band and all, instead of stopping halfway across it;
%     content wider  - the spacer collapses to nothing and the table SCROLLS
%                      sideways with every real column at its own content width.
%   It costs one column of '' in ColumnName and in Data (paintProgressTable), which
%   is why the two tables have one more column than they have steps.
%
%   Known and accepted: 'fit' re-fits when content changes, so a column whose header
%   is narrower than 'running' (BFI, Guided, Regions) shifts by 8-18 px as a run
%   walks its cells - twice per run, since the column takes the widest cell of all
%   its rows.  Pin those columns to a number if it ever reads as noise.
if ~isgraphics(h), return; end
lbls = cell(1,numel(ids));
for i = 1:numel(ids), lbls{i} = reg(strcmp({reg.id},ids{i})).label; end
h.ColumnName  = [{firstName}, lbls, {''}];
h.ColumnWidth = [repmat({'fit'},1,numel(ids)+1), {'1x'}];
end

function refreshCells(fig,identities)
%refreshCells  Repaint the state tables.  They carry table DATA, never components:
%   the whole point of the read-only monitor is that a 200-file project costs two
%   uitables instead of 3000 live widgets.  'identities' narrows the repaint to the
%   rows of those recordings ({} = all), which is what keeps a per-cell progress
%   tick cheap during a run - it must stay a NARROW repaint, not a rebuild.
if nargin<2, identities = {}; end
paintProgressTable(fig,'raw',identities);
paintProgressTable(fig,'derived',identities);
end

function paintProgressTable(fig,kind,identities)
%paintProgressTable  One table's worth of refreshCells.
app = getApp(fig); c = app.c.process;
if strcmp(kind,'raw')
    fld = 'rawTable'; rows = app.pRows;
else
    fld = 'derivedTable'; rows = app.dRows;
end
if ~isfield(c,fld) || ~isgraphics(c.(fld)), return; end
h = c.(fld);
if isempty(rows), h.Data = {}; return; end
ids  = monitorColumns(app.reg, kind);
D    = h.Data;
nCol = numel(ids)+2;         % the name, one per step, and the '1x' spacer (D5)
if ~iscell(D) || size(D,1)~=numel(rows) || size(D,2)~=nCol
    D = repmat({''}, numel(rows), nCol);               % '' so the spacer draws blank
    identities = {};                                   % a fresh grid: fill it all
end
if isempty(identities)
    which = 1:numel(rows);
else
    which = find(ismember({rows.identity}, identities));
end
for i = which
    r = rows(i);
    D{i,1} = r.rowLabel;
    planned = plannedStepsFor(app, r);
    for s = 1:numel(ids)
        D{i,s+1} = progressCellText(app, r, stepById(app.reg,ids{s}), planned);
    end
end
h.Data = D;
end

function w = monitorWidthPolicy(fig, kind)
%monitorWidthPolicy  One monitor table's ColumnWidth, for a headless check of D5.
%   The widths themselves are the layout engine's and are not readable from MATLAB;
%   what IS checkable is that this window asks for content sizing rather than the
%   fixed numbers that cropped.
c = getApp(fig).c.process;
if strcmp(kind,'raw'), w = c.rawTable.ColumnWidth; else, w = c.derivedTable.ColumnWidth; end
if ~iscell(w), w = {w}; end
end

function ids = plannedStepsFor(app, r)
%plannedStepsFor  Which steps THIS ROW's configuration runs on it, in registry
%   order.  The row hands over its BRANCH rather than having one guessed from a file
%   name: a DERIVED row IS one branch, so it gets exactly that (type,branch) row's
%   selection - what it runs plus what it inherits from the anchor row - and never
%   the other pipeline's steps.  A RAW row stands for the whole recording and gets
%   the union of its type's rows, unless the file itself carries a stage flag (a
%   product file in the working set), in which case it stands for its own pipeline.
%   The per-animal steps are added for the animal, not the row: they span every
%   type by definition.
if isempty(r.type), ids = plannedStepIds(app, ''); return; end
brs = wbTypeSelection('rows', app.typeSel, regFor(app,r.type), r.type);
if strcmp(r.kind,'derived')
    brs = brs(strcmp(brs, r.branch));                  % the row IS a branch
elseif ~isempty(r.branch) && any(strcmp(r.branch, brs))
    brs = {r.branch};
end
ids = plannedStepIds(app, r.type, brs);
end

function out = orderIds(reg, ids)
%orderIds  Unique step ids in registry (dependency) order.
out = {reg(ismember({reg.id}, ids)).id};
end

function txt = progressCellText(app, r, step, planned)
%progressCellText  ONE cell of the monitor (spec D4).  The words are the states a
%   watcher cares about, in this order of authority: what THIS RUN is doing with
%   the cell - or is about to - then what the disk says, then what the
%   configuration intends.
%
%   THE OVERLAY OUTRANKS THE DISK, INCLUDING FOR 'queued' (spec D2).  That is the
%   whole point of the pre-run pass: a step being re-processed has a finished result
%   on disk, and reading 'done' off it while the run is on its way to overwrite it
%   described the last run rather than this one.  Once the run has marked a cell,
%   the run owns the word until it ends.
key = cellKey(r.identity, step.id);
if isKey(app.runState,key) && stepTouchesFile(step, r.branch)
    switch app.runState(key)
        % 'running' carries NO percentage.  A wrapper says what it started and what
        % it finished and nothing in between, so there is no fraction to show and a
        % made-up one would be worse than none.
        case 'queued',  txt = 'queued';  return
        case 'running', txt = 'running'; return
        case 'done',  txt = 'done';  return
        case 'error', txt = 'error'; return
    end
end
st = rowStateOf(app, r, step.id);
if strcmp(st,'done'),                        txt = 'done'; return; end
if ~any(strcmp(step.id, planned)),           txt = char(183); return; end   % not configured
if any(strcmp(st,{'ready','stale'})),        txt = 'queued'; return; end
% configured, but blocked: queued only if everything it needs is itself queued
have = {};
req = wbPrereqs('all', step);
for i = 1:numel(req)
    rs = rowStateOf(app, r, req{i});
    if strcmp(rs,'done') || any(strcmp(req{i}, planned)), have{end+1} = req{i}; end %#ok<AGROW>
end
if wbPrereqs('met', step, have), txt = 'queued'; else, txt = 'skipped'; end
end

function s = rowStateOf(app, r, stepId)
%rowStateOf  The disk state a monitor row reads for one step.  A RAW row is a file,
%   so it asks the per-FILE picture (D6's whole point).  A DERIVED row is a
%   (recording, PIPELINE) pair, so it asks the per-branch picture - which is also
%   what makes D7 work: a row whose product is not on disk yet simply has no branch
%   state, reads 'unavailable', and shows what is QUEUED for it instead of nothing.
if ~strcmp(r.kind,'derived'), s = fileStateOf(app, r.path, stepId); return; end
s = 'unavailable';
key = [r.identity '||' r.branch];
if isKey(app.branchState,key)
    bs = app.branchState(key);
    if isfield(bs,stepId), s = bs.(stepId); end
end
if strcmp(s,'done') && isKey(app.stale, cellKey(r.identity,stepId))
    s = 'stale';                                       % an in-session edit re-opened it
end
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

function onProgressSelect(fig, kind, ev)
%onProgressSelect  Remember which row of WHICH table the user clicked, and describe
%   it below - selecting a row still shows its file, from either half.
if isempty(ev) || ~hasMember(ev,'Indices') || isempty(ev.Indices), return; end
app = getApp(fig);
app.pSel = struct('kind',char(kind),'idx',ev.Indices(1,1));
setApp(fig,app);
refreshProgressDetail(fig);
end

function tf = hasMember(v, name)
%hasMember  Does this event carry that field - WITHOUT converting it.  A real
%   CellSelectionChangeData is an OBJECT, and struct(obj) is the deprecated
%   conversion MATLAB warns about once per call; a headless test passes a plain
%   struct instead, and both callers have to keep working.  So ask each of them in
%   its own language rather than making one look like the other.
if isstruct(v), tf = isfield(v,name); else, tf = isprop(v,name); end
end

function refreshProgressDetail(fig)
%refreshProgressDetail  The selected row's identity card: where its file is, how it
%   is labelled, and what went wrong if something did.
app = getApp(fig); c = app.c.process;
if ~isfield(c,'detail') || ~isgraphics(c.detail), return; end
sel = struct('kind','raw','idx',0);
if isfield(app,'pSel') && isstruct(app.pSel), sel = app.pSel; end
if strcmp(sel.kind,'derived'), rows = app.dRows; else, rows = app.pRows; end
i = sel.idx;
if i<1 || i>numel(rows)
    c.detail.Text = 'Select a row to see its file.';
    return
end
r = rows(i);
bits = {sprintf('%s%s   animal %s | type %s | group %s', r.path, ...
    ternary(isempty(r.flag),'',['   product ' r.flag]), ...
    dashIfEmpty(r.animal), dashIfEmpty(r.type), dashIfEmpty(r.expGroup))};
blocked = firstSkipped(app, r);
if ~isempty(blocked)
    bits{end+1} = sprintf('%s skipped - %s', blocked, unavailableReason(app,r,blocked));
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

function s = unavailableReason(app,r,stepId)
%unavailableReason  A short line for a step that cannot run (modality / prereqs).
%
%   THE MODALITY IS THE ROW'S OWN FILE, not the session's.  It used to be app.modality
%   - the working set's most common modality - which was the only modality there was
%   when a session could hold only one.  Now that a session can hold several, saying
%   'not for LSCI' about a video would name the wrong recording's kind entirely; the
%   gate itself has always been per file (wbStateEngine>applicable), so this just
%   quotes the same fact it read.
s = 'not available yet';
step = stepById(app.reg,stepId);
if isempty(step), return; end
md = '';
if isstruct(r) && isfield(r,'model') && isfield(r.model,'modality'), md = r.model.modality; end
if isKey(app.base,r.identity) && ~isempty(md) && ~any(strcmp(md,step.modalities))
    s = ['not for ' md];
elseif ~isempty(wbPrereqs('all',step))
    s = ['needs: ' wbPrereqs('describe',step)];
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
%stepInfoText  The one line under the Type/Step pair: what this step needs of you
%   and of the data, in plain words.  Everything the reader cannot act on - which
%   settings field gates it, what the wrapper is called - stays out.
bits = {};
if strcmp(step.arity,'perAnimal'), bits{end+1} = 'runs once per animal'; end
if ~isequal(step.interactive,false), bits{end+1} = 'asks you to draw or click'; end
if step.needsRaw, bits{end+1} = 'also reads the raw recording'; end
if ~isempty(wbPrereqs('all',step)), bits{end+1} = ['needs ' wbPrereqs('describe',step) ' first']; end
if isempty(bits), bits{end+1} = 'runs on its own, once its input exists'; end
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
    % The row label is the FIELD NAME unless the registry names a friendlier one in
    % step.labels.  Most fields are their own best label here (a protocol is written
    % in terms of vFR and nHarm), so the override is opt-in per field rather than a
    % table everything has to be listed in.
    lbl = uilabel(gg,'Text',getfieldOr(step.labels,f,f));
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
    uidropdown(parent,'Items',items,'Value',v,'ValueChangedFcn',ui(fig,@(o,~)cb(o)));
elseif islogical(defaultType(step,field))
    uicheckbox(parent,'Text','','Value',logical(firstTrue(val)),'ValueChangedFcn',ui(fig,@(o,~)cb(o)));
elseif iscell(val)
    uieditfield(parent,'text','Value',cellToStr(val),'ValueChangedFcn',ui(fig,@(o,~)cb(o)), ...
        'Tooltip','comma-separated list');
else
    uieditfield(parent,'text','Value',val2str(val),'ValueChangedFcn',ui(fig,@(o,~)cb(o)));
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
renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
wbLog(fig,sprintf('edit %s.%s -> %s',stepId,field,val2str(value)));
end

function applyInvalidation(fig,seed,type)
%applyInvalidation  Mark every DONE cell in the forward set of the edit stale.
%   With a type, only that type's recordings are touched - which is the whole
%   point of keying the settings by type: a BP edit must leave BV's done cells
%   alone.  The seed is a STEP ID (or a shared-key name), never a bare field, so
%   'method' edited on BFI cannot invalidate the internal cycle.
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
renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
renderConstructor(fig);            % a preset carries the per-type layers too
refreshPresetDrop(fig);
wbLog(fig,['loaded preset ' pth]);
autosaveSession(fig);              % the whole settings bag changed
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
    renderProgress(fig);   % the monitor's ROWS, columns and range follow the configuration
    renderConstructor(fig);        % every per-type layer went with the bag
    wbLog(fig,'reset to launcher defaults (global AND per-type settings)');
    autosaveSession(fig);          % the whole settings bag changed
    return
end
loadPreset(fig,fullfile(getApp(fig).presetDir,name));
end
function ensureDir(d), if ~isfolder(d), mkdir(d); end, end

%% ===================== session ===================================== %%
function uiSaveSession(fig)
%uiSaveSession  Save a NAMED session.  It opens in the workbench sessions folder,
%   which is also where the last session lives, so the two sit side by side.
[f,p] = uiputfile('*.mat','Save session',fullfile(sessionDir(),'workbench_session.mat'));
if isequal(f,0), return; end
saveSessionTo(fig,fullfile(p,f));
end
function uiLoadSession(fig)
[f,p] = uigetfile('*.mat','Load session',sessionDir());
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
%   A run that is killed, cancelled or crashes must still leave something resumable
%   behind, so the session is rewritten on EVERY state change - a scan, a label, a
%   reference, a Constructor tick, a settings edit, each finished cell, the end of a
%   run and the exit itself - rather than only when the user remembers to save.
%
%   IT ALWAYS WRITES THE LAST SESSION (lastSessionPath), whatever else it writes.
%   That file is the workbench's own resume point and does not depend on the user
%   having chosen a location; when a session HAS been saved or loaded somewhere
%   else, that copy is kept in step too, so 'Save session...' does not quietly stop
%   tracking the moment it is used.
if ~isvalid(fig), return; end
app = getApp(fig);
if isempty(app.files) && isempty(app.root), return; end     % nothing to resume yet
try
    targets = {lastSessionPath()};
    if ~isempty(app.sessionPath) && ~strcmp(app.sessionPath,targets{1})
        targets{end+1} = app.sessionPath;
    end
    session = wbSession('empty');
    for i = 1:numel(targets), saveSessionInto(app, session, targets{i}); end
    if isempty(app.sessionPath), app.sessionPath = targets{1}; setApp(fig,app); end
catch ME
    wbLog(fig,['session autosave failed: ' ME.message]);
end
end

function saveSessionInto(app, session, pth)
%saveSessionInto  Fill a blank session from the app struct and write it.  ONE
%   definition of what a session contains, used by both the explicit Save and the
%   autosave, so the two can never drift apart.
session.root          = app.root;
session.resultsRoot   = app.resultsRoot;
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
session.staleOverlay  = app.stale;
session.presetRef     = app.presetRef;
session.reprocess     = logical(app.reprocess);
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
app.resultsRoot  = session.resultsRoot;
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
% (type,branch) rows existed is upgraded to them instead of being read as empty.
% AGAINST THE WHOLE REGISTRY, not this window's: the files have not been read back
% yet, so which modalities the session holds is not known here - and fromCells
% drops a key naming a step the registry does not have, which against a
% modality-filtered registry would quietly delete another modality's ticks.
app.typeSel   = wbTypeSelection('fromCells', keys(session.typeSel), wbStepRegistry());
app.animalSel = session.animalSel;
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
setReprocess(fig, session.reprocess);    % after the working set: it repaints the range
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
%
%   THE FROM/TO RANGE (spec D1) is applied on TOP of that expansion, in runOrderIn:
%   it filters the entry list by registry column, and nothing else.  Whether the
%   "already on disk" pruning applies is a SEPARATE question, answered by the
%   Re-process tick (spec D6).  The expansion itself knows about neither.
app = getApp(fig);
entries = runOrderIn(fig, app.rangeFrom, app.rangeTo);
end

function entries = expandEntries(app, prune)
%expandEntries  THE expansion proper (see buildRunOrder for what it asks and why).
%   'prune' is the ONE thing the Process tab changes about it: with it false the
%   already-on-disk test is skipped, so a finished step runs again with new settings
%   (spec D6, the Re-process tick).  Everything else stands - one work item per
%   (file,step), the per-animal steps once per animal - because those are not resume
%   rules, they are what makes the list correct.
%
%   THE PER-ANIMAL STEPS ARE PRUNED TOO (Phase 6).  Until this phase they were
%   emitted unconditionally, so registration and vessel typing re-ran on every Run,
%   their column could never read as finished, and 'Last valid' was pinned to it
%   forever.  That was wrong twice over: vessel typing is an interactive paint step,
%   so a resume aimed at a late column re-opened its window, and the resume frontier
%   could never move past it.  They now answer the same question every other step
%   does - see animalStepAlreadyDone for what 'done' means when the scope is a whole
%   animal.
reg = app.reg;
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
        if prune && animalStepAlreadyDone(app, r.animalIdx, step), continue; end
        entries(end+1) = mkEntry(r.animalIdx, 1, stepIndexOf(reg,step.id), step, ...
            animalRefIdentity(app,r.animalIdx), r.animal, true, ty); %#ok<AGROW>
    end

    % ---- the type's rows, and what each of them runs on this file ------------
    % asked of the TYPE's registry (its own modality), while the step's COLUMN
    % INDEX stays the session's, so one step-major order covers every modality
    if isempty(ty), continue; end
    treg = regFor(app, ty);
    brs = wbTypeSelection('rows', app.typeSel, treg, ty);
    for b = 1:numel(brs)
        ids = wbTypeSelection('steps', app.typeSel, treg, ty, brs{b});
        for j = 1:numel(ids)
            step = stepById(treg, ids{j});
            if isempty(step) || strcmp(step.arity,'perAnimal'), continue; end
            k = cellKey(r.identity, ids{j});
            if isKey(seen,k), continue; end             % one work item per (file,step)
            if prune && stepAlreadyDone(app, r.identity, ty, brs{b}, step), continue; end
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

%% ---- the FROM/TO run range (spec D1-D5) ---------------------------- %%
function [entries, info] = runOrderIn(fig, fromSel, toSel)
%runOrderIn  The run order for ONE range.  Two readings of the same configuration:
%   the PLANNED list (everything the configuration would ever run) says which
%   columns exist and is what a Re-process run executes; the PENDING list (today's
%   expansion, with the already-done steps pruned) is what an ordinary run executes.
%
%   THE RESUME PATH IS UNCHANGED.  With From = 'Last valid' the slice starts at the
%   frontier, and every column before the frontier is by definition finished - so
%   it has no pending entries and the slice removes nothing.  'Last valid' to the
%   last column therefore returns exactly what buildRunOrder returned before the
%   range existed, entry for entry, in the same order (the regression gate).
app = getApp(fig);
[info, planned, pending] = runRangeInfo(app, fromSel, toSel);
if info.forced, src = planned; else, src = pending; end
entries = wbRunRange('slice', app.reg, src, info.fromId, info.toId);
end

function [info, planned, pending] = runRangeInfo(app, fromSel, toSel)
%runRangeInfo  Everything the toolbar and the run headline need about one range:
%   the configured columns, the resume frontier, the two resolved step ids, and
%   whether finished work is to be done again.  The columns come from the PLANNED
%   expansion on purpose - a finished column is still configured, and dropping it
%   would take away the only way to ask for it again.
%
%   'forced' IS THE CHECKBOX, AND ONLY THE CHECKBOX (spec D6).  It used to be a side
%   effect of naming a From column, which is a coupling nobody could see in the
%   window: the same gesture asked for two unrelated things, and there was no way to
%   start at segmentation without also overwriting it.  From/To is now a pure range.
planned = expandEntries(app, false);
pending = expandEntries(app, true);
pendIds = {};
if ~isempty(pending), pendIds = unique({pending.stepId}); end
cols = wbRunRange('columns', app.reg, planned);
fron = wbRunRange('frontier', app.reg, cols, @(id) ~any(strcmp(id, pendIds)));
[fromId, toId] = wbRunRange('resolve', app.reg, cols, fromSel, toSel, fron);
info = struct('cols',{cols},'frontier',fron,'fromId',fromId,'toId',toId, ...
    'forced',logical(getfieldOr(app,'reprocess',false)));
end

function info = currentRangeInfo(fig)
app = getApp(fig);
info = runRangeInfo(app, app.rangeFrom, app.rangeTo);
end
function cols = runRangeColumns(fig)
%runRangeColumns  The configured columns, in registry order (the dropdowns' items).
info = currentRangeInfo(fig); cols = info.cols;
end
function id = runRangeFrontier(fig)
%runRangeFrontier  Where 'Last valid' currently resolves to (spec D3).
info = currentRangeInfo(fig); id = info.frontier;
end

function refreshRangeControls(fig)
%refreshRangeControls  Repopulate From/To from buildRunOrder's OWN expansion, so a
%   Constructor edit is reflected without a second source of truth (spec D5).  To
%   offers only the columns at or after the resolved From and is clamped when From
%   moves past it (spec D4).  With nothing configured both read '(nothing
%   configured)' and Run stays enabled - it says so in the log, as it always did.
app = getApp(fig); c = app.c.process;
if ~isfield(c,'fromDrop') || ~isgraphics(c.fromDrop), return; end
info = runRangeInfo(app, app.rangeFrom, app.rangeTo);
cols = info.cols;
if isempty(cols)
    setDropItems(c.fromDrop, {noRangeItem()}, {''}, '');
    setDropItems(c.toDrop,   {noRangeItem()}, {''}, '');
    return
end
lastTok = wbRunRange('lastValid');
fromItems = [{sprintf('Last valid (%s)', wbRunRange('label',cols,info.frontier))}, {cols.label}];
fromData  = [{lastTok}, {cols.id}];
% show what is STORED, not what it resolved to: the sentinel is a standing
% instruction ("wherever the frontier is now") and replacing it with today's answer
% would quietly freeze it the first time the toolbar was painted
fromVal = char(app.rangeFrom);
if isempty(fromVal), fromVal = lastTok; end
setDropItems(c.fromDrop, fromItems, fromData, fromVal);

k = wbRunRange('index', cols, info.fromId);
if k<1, k = 1; end
setDropItems(c.toDrop, {cols(k:end).label}, {cols(k:end).id}, info.toId);
% NOTE the resolved ids are shown but NOT stored back.  An unset To means "to the
% end", and freezing it to whatever the last column happened to be would pin the
% range the first time it was painted - the next Constructor tick would then add a
% column the range silently refused to reach.
end

function setDropItems(h, items, data, value)
%setDropItems  Set a dropdown's items and value together, without a transient state
%   in which the old Value is not in the new Items (which uidropdown rejects).
if ~isgraphics(h), return; end
h.Items = items;
h.ItemsData = data;
k = find(strcmp(data, char(value)),1);
if isempty(k), k = 1; end
h.Value = data{k};
end

function onRangeEdit(fig, which, value)
%onRangeEdit  One dropdown moved: store it, clamp the other, repaint.  The rule
%   itself lives in wbRunRange; this only decides what is remembered.
app = getApp(fig);
if strcmp(which,'from'), app.rangeFrom = char(value); else, app.rangeTo = char(value); end
setApp(fig,app);
clampStoredTo(fig);
refreshRangeControls(fig);
wbLog(fig, rangeHeadline(currentRangeInfo(fig)));
end

function setRunRange(fig, fromSel, toSel)
%setRunRange  The programmatic twin of the two dropdowns (headless-testable).  An
%   empty From means 'Last valid'; an empty To means "to the last column".
app = getApp(fig);
if nargin>=2, app.rangeFrom = char(fromSel); end
if nargin>=3, app.rangeTo   = char(toSel);   end
if isempty(app.rangeFrom), app.rangeFrom = wbRunRange('lastValid'); end
setApp(fig,app);
clampStoredTo(fig);
refreshRangeControls(fig);
end

function setReprocess(fig, tf)
%setReprocess  The 'Re-process finished files' switch (spec D6) - THE only thing
%   that decides whether a run repeats work that is already on disk.  It is not a
%   window preference like the PDF switches: it changes what the next Run does, it
%   is part of what the session remembers, and the pruning rule it turns off is
%   stepAlreadyDone / animalStepAlreadyDone, which stay the single definition of
%   "already there" (expandEntries consults them or does not).
app = getApp(fig); app.reprocess = logical(tf); setApp(fig,app);
c = app.c.process;
if isfield(c,'reproCheck') && isgraphics(c.reproCheck), c.reproCheck.Value = app.reprocess; end
refreshRangeControls(fig);
end

function onReprocessEdit(fig, tf)
%onReprocessEdit  The checkbox itself: set it, then say in the log what the run will
%   now do - the same courtesy the two range dropdowns pay.
setReprocess(fig, tf);
wbLog(fig, rangeHeadline(currentRangeInfo(fig)));
end

function clampStoredTo(fig)
%clampStoredTo  Spec D4, made durable: an EXPLICIT To that now sits before From is
%   moved up to From rather than left to be clamped afresh on every read - the
%   dropdown would otherwise show one value and remember another.  An unset To is
%   left alone: it means "to the end" and has nothing to clamp.
app = getApp(fig);
if isempty(app.rangeTo), return; end
info = runRangeInfo(app, app.rangeFrom, app.rangeTo);
if isempty(info.cols), return; end
if wbRunRange('index',info.cols,app.rangeTo) < wbRunRange('index',info.cols,info.fromId)
    app.rangeTo = info.fromId;
    setApp(fig,app);
end
end

function v = runRangeValue(fig)
%runRangeValue  What the range currently resolves to - the headless view of the two
%   dropdowns, the frontier, and the Re-process verdict that now stands beside them.
info = currentRangeInfo(fig);
v = struct('from',getApp(fig).rangeFrom,'to',getApp(fig).rangeTo, ...
    'fromId',info.fromId,'toId',info.toId,'frontier',info.frontier,'forced',info.forced);
end

function s = rangeHeadline(info)
%rangeHeadline  One line naming the range and what will happen to work that is
%   already finished inside it.  The two are now INDEPENDENT (spec D6): From/To says
%   how much of the protocol to walk, the Re-process tick says whether finished
%   recordings are done again - so the sentence names them separately instead of
%   implying that picking a start column also means "redo it".
if isempty(info.cols)
    s = 'range: nothing set up yet - tick some steps on the Constructor tab.';
    return
end
s = sprintf('range: %s -> %s.  %s', ...
    wbRunRange('label',info.cols,info.fromId), wbRunRange('label',info.cols,info.toId), ...
    ternary(info.forced, ...
        'Everything in that range will run again, finished or not.', ...
        'Anything already finished is left as it is.'));
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
    brs = wbTypeSelection('rows', app.typeSel, regFor(app,type), type);
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

function tf = animalStepAlreadyDone(app, animalIdx, step)
%animalStepAlreadyDone  The same question as stepAlreadyDone, asked of a step whose
%   scope is a whole ANIMAL (spec D3) - resolved in Phase 6, where they used to be
%   emitted unconditionally.
%
%   WHAT 'DONE' MEANS HERE.  A per-animal step is handed ONE call covering every
%   configured product of every recording of the animal (its branchScope fans out
%   inside the animal's column - wbExecutor>buildFNames), so it is done only when
%   ALL of them carry its gating field.  Anything less and the call has work left:
%   a recording added to the animal after the last run, a product the type only
%   started producing later, or a settings edit that re-opened the step, each puts
%   the whole animal back in the queue.  That is the same conservative reading
%   stepAlreadyDone applies to a recording-level per-file step, one scope up.
%
%   An animal with nothing configured is never 'done' - there is no product to
%   judge it on, so it stays queued rather than being silently skipped.
tf = false;
seenAny = false;
for i = 1:numel(app.rows)
    if app.rows(i).animalIdx ~= animalIdx, continue; end
    id  = app.rows(i).identity;
    ty  = typeOfIdentity(app,id);
    brs = wbTypeSelection('rows', app.typeSel, regFor(app,ty), ty);
    for b = 1:numel(brs)
        key = [id '||' brs{b}];
        if ~isKey(app.branchState,key), return; end          % nothing produced there yet
        bs = app.branchState(key);
        if ~isfield(bs,step.id) || ~strcmp(bs.(step.id),'done'), return; end
        if isKey(app.stale, cellKey(id,step.id)), return; end % an edit re-opened it
        seenAny = true;
    end
end
tf = seenAny;
end

function lines = dryRun(fig)
%dryRun  Preview the SLICE the current range would execute, named first so it is
%   never a mystery why a preview is shorter than the configuration.
app = getApp(fig);
[entries, info] = runOrderIn(fig, app.rangeFrom, app.rangeTo);
lines = {rangeHeadline(info), ...
    sprintf('=== PREVIEW: %d steps would run.  Nothing is processed. ===',numel(entries))};
for i = 1:numel(entries)
    e = entries(i);
    if strcmp(e.arity,'perAnimal'), who = ['ANIMAL ' e.animal]; else, who = e.label; end
    lines{end+1} = sprintf('  %2d. [%s] %-26s :: %s', i, e.animal, who, e.stepLabel); %#ok<AGROW>
end
if numel(entries)==0
    lines{end+1} = '  Nothing to do - set the types up on the Constructor tab, or widen the range.';
end
setLog(fig,lines);
end

%% ===================== real Run (executor) ========================= %%
function runChecked(fig)
%runChecked  Run what the per-type configuration says, through wbExecutor.
app = getApp(fig);
if isfield(app,'running') && app.running, return; end       % already running
% ONE USER ACTION FOR THE WHOLE RUN (spec D4).  The Run button's own callback is
% wrapped, but the latch has to outlive the click: the confirmation dialog, the
% executor loop and the unwind below all drawnow, and each of those dispatches
% whatever the user clicked while the window was frozen.  Taken here, released by
% the token when this function unwinds - alongside finishRun, not inside a try.
% IT IS TAKEN BEFORE the finishRun cleaner below ON PURPOSE: cleanup objects are
% cleared newest-first, so this one goes LAST and finishRun still repaints - and
% still dispatches - with the latch held.  Do not swap the two lines.
latch = wbUiGuard('hold', fig); %#ok<NASGU>
[entries, info] = runOrderIn(fig, app.rangeFrom, app.rangeTo);
if isempty(entries)
    setLog(fig,{rangeHeadline(info), ...
        'Nothing to do - set the types up on the Constructor tab, or widen the range.'});
    return
end
if ~confirmRun(fig, entries, info), return; end
% a fresh run starts with a clean transient overlay
app.runState = containers.Map('KeyType','char','ValueType','any');
app.cellMsg  = containers.Map('KeyType','char','ValueType','any');
app.running  = true;  app.cancel = false;
setApp(fig,app);
% how many entries each COLUMN of THIS run has - counted from the list about to be
% executed (a From/To slice), never re-expanded mid-run.  See beginPdfColumns.
beginPdfColumns(fig, entries);
setRunningUI(fig,true);
autosaveSession(fig);                    % a run that never returns still leaves a session
% The pool token has to be taken BEFORE anything runs a parfor, or "was there a pool
% before this run" is no longer knowable.  It starts nothing; finishRun releases.  It
% rides in the cleanup closure BY VALUE rather than in app, so the workers still go
% back even if the figure itself is gone by the time this unwinds.
poolTok = wbPool('open');
cleaner = onCleanup(@() finishRun(fig,poolTok)); % restore UI + fold results even on error
ctx = buildExecContext(fig);
wbExecutor(entries, ctx);
end

function ok = confirmRun(fig, entries, info)
%confirmRun  The pre-run summary.  A type-level tick can expand into hundreds of
%   file-steps, and the number is the only honest warning about how long this will
%   take - so it is shown before anything runs, not discovered afterwards.  The
%   RANGE is on the same card, and so is the Re-process tick: overwriting finished
%   work has to be visible at the moment of consent, not only in the toolbar, and
%   the tick is easy to leave on from the last time.
%
%   AND SO IS THE DELETION (item 9, spec D10).  Recomputing a source file removes
%   the results derived from it, which is destructive and irreversible, so it is
%   STATED BEFORE the user says Run rather than discovered in the log afterwards.
%   The count is taken from the very list this dialog is describing - a dry run of
%   the same wbProducts question the executor will ask - so the number on the card
%   is the number of files that will actually go.
ok = true;
n  = numel(entries);
tys = unique({entries.type},'stable');
tys = tys(~cellfun(@isempty,tys));
sts = unique({entries.stepId},'stable');
msg = sprintf('%s\nAbout to run %d %s across %d recording %s - %d jobs in all.%s%s  Continue?', ...
    rangeHeadline(info), numel(sts), plural(numel(sts),'step'), ...
    numel(tys), plural(numel(tys),'type'), n, ...
    ternary(info.forced,'  Results that are already saved will be overwritten.',''), ...
    supersededSentence(numel(supersededFiles(getApp(fig), entries))));
wbLog(fig,msg);
if ~isvalid(fig) || strcmp(fig.Visible,'off'), return; end   % headless: never block
sel = uiconfirm(fig, msg, 'Run', 'Options',{'Run','Cancel'}, ...
    'DefaultOption',1, 'CancelOption',2, 'Icon','question');
ok = strcmp(sel,'Run');
if ~ok, wbLog(fig,'Run cancelled at the summary.'); end
end

function s = supersededSentence(n)
%supersededSentence  The author's own line on the consent card, or nothing at all
%   when this run deletes nothing - which is the ordinary case, and a warning that
%   fires every time is a warning nobody reads.
s = '';
if n <= 0, return; end
s = sprintf('  Recomputing source files will delete %d previously computed %s.', ...
    n, plural(n,'result'));
end

function files = supersededFiles(app, entries)
%supersededFiles  A DRY run of the executor's own deletion question over the
%   entries about to run: which results this run will remove.  Deduped across
%   entries, because two producers of one recording can only ever name disjoint
%   sets but the same file must never be counted twice if they ever did.
%
%   It reads the disk, once per producer entry - the only entries that survive the
%   outKind test - and it is asked at the moment a modal dialog is being built, so
%   the cost is paid where the user is already waiting for an answer.
files = cell(1,0);
if isempty(entries), return; end
for i = 1:numel(entries)
    step = stepById(app.reg, entries(i).stepId);
    if isempty(step) || ~strcmp(step.outKind,'new'), continue; end
    ms = entryModels(app, entries(i));
    for k = 1:numel(ms)
        stage = wbProducts('writes', step, contrastStageForModel(app, ms(k)));
        files = [files, wbProducts('below', app.reg, ms(k), step.id, stage)]; %#ok<AGROW>
    end
end
files = unique(files,'stable');
end

function ms = entryModels(app, e)
%entryModels  The recordings one run entry touches: the animal's members for a
%   per-animal step, its own recording otherwise.
if strcmp(e.arity,'perAnimal')
    sel = app.rows([app.rows.animalIdx]==e.animalIdx);
    if isempty(sel), ms = []; return; end
    [~,ord] = sort([sel.rowInAnimal]);
    ms = [sel(ord).model];
    return
end
ms = [];
k = find(strcmp({app.rows.identity}, e.identity),1);
if ~isempty(k), ms = app.rows(k).model; end
end

function ctx = buildExecContext(fig)
%buildExecContext  The callback bundle wbExecutor drives the figure through.
ctx = struct();
ctx.reg           = getApp(fig).reg;
% THE TWO PROJECT FOLDERS, and this is the only place the run is told about them:
% wbExecutor puts them on every wrapper's s, which is the GUI's equivalent of the
% two lines a launcher carries in its STEP 0.  The report pages follow for free -
% reportOpen reads the same two fields (results-folder plan, session 1).
ctx.rootFolder    = getApp(fig).root;
ctx.resultsFolder = getApp(fig).resultsRoot;
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
% THE WORKING SET IS A FENCE (round-2 item 8).  The input glob that fills a
% wrapper's fNames is a directory listing, so a product an EARLIER session left in
% the folder used to join it; this is the answer that stops it, taken from the same
% configuration the derived monitor table is built from.  It is resolved ONCE, here,
% for every recording of the run - the configuration cannot change while the latch
% is held, and the alternative is a settings resolve per candidate file.
ctx.admits        = admitsFcn(getApp(fig));
ctx.setState      = @(id,sid,state,msg) execSetState(fig,id,sid,state,msg);
ctx.log           = @(msg) wbLog(fig,msg);
ctx.isCancelled   = @() execIsCancelled(fig);
ctx.modalGuard    = @(fcn) wbModalGuard(fig,fcn);
ctx.afterDone     = @(id,sid,mdl) execAfterDone(fig,id,sid,mdl);
end

function f = admitsFcn(app)
%admitsFcn  ctx.admits, with the per-recording answer resolved up front.  The map
%   is keyed by recording identity, and a recording the run has never heard of
%   falls through UNFENCED - this hook may narrow what a session uses, never what
%   it can see.
m = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(app.rows)
    id = app.rows(i).identity;
    if ~isKey(m,id), m(id) = admissibleFlagsFor(app,id); end
end
f = @(id,p) admitsWith(m,id,p);
end
function tf = admitsWith(m, identity, path)
id = char(identity);
if ~isKey(m,id), tf = true; return; end
tf = wbProducts('admits', m(id), path);
end

function fl = admissibleFlagsFor(app, identity)
%admissibleFlagsFor  WHICH STAGE FLAGS THIS RECORDING'S CONFIGURATION PRODUCES
%   (spec D7) - what the run may consume and what the reports may list.  It is the
%   configuration that answers, not the disk: the same wbTypeSelection rows the
%   derived monitor table is built from, so "listed in the table" and "allowed into
%   the run" are one fact by construction.
%
%   NO CONFIGURATION MEANS NO FENCE.  A recording whose type has no ticked row
%   falls back to the branches its own files on the Files tab carry, and if those
%   name no branch at all (a working set of raw recordings) the answer is EMPTY,
%   which wbProducts reads as "no answer" and admits everything.  A recording is
%   never fenced down to nothing.
ty  = typeOfIdentity(app, identity);
reg = regFor(app, ty);
if ~isempty(ty) && ~isempty(wbTypeSelection('rows', app.typeSel, reg, ty))
    fl = wbProducts('flags', reg, app.typeSel, ty, contrastStageForType(app,ty));
    return
end
fl = filesTabFlags(app, identity);
end

function fl = filesTabFlags(app, identity)
%filesTabFlags  The fall-back admissible set: the branches this recording's files
%   in the WORKING SET actually carry.  {} when none of them carries one, which is
%   how "no fence" is said.
fl = {};
for i = 1:numel(app.files)
    if ~strcmp(app.files(i).model.identity, identity), continue; end
    st = app.files(i).model.stage;
    if isempty(st) || any(strcmp(st,fl)), continue; end
    fl{end+1} = st; %#ok<AGROW>
end
if ~isempty(fl), fl = [{''}, fl]; end       % the raw recording rides along
end

function st = contrastStageForType(app, type)
%contrastStageForType  The stage flag a settings-driven producer writes for one
%   TYPE - the type-level twin of contrastStageForModel, and the same rule
%   (rowFlagFor > settingStage).  WHICH producer is settings-driven is derived,
%   not named: it is the one whose resolved settings choose a stage at all.
st = '';
reg = regFor(app, type);
brs = wbTypeSelection('branches', reg);
for b = 1:numel(brs)
    pid = wbTypeSelection('producer', reg, brs{b});
    if isempty(pid), continue; end
    alt = settingStage(app, type, stepById(reg,pid));
    if ~isempty(alt), st = alt; return; end
end
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
%
%   IT TAKES A LIST (spec D3).  'identity' is one recording or a cellstr of them,
%   because the executor marks a whole call at once - the pre-run 'queued' pass over
%   a 200-file project is 3000 cells, and one repaint per cell would be 3000 full
%   table repaints before the first wrapper was even called.  One map pass, one
%   refreshCells, one drawnow, whatever the length of the list.
if ~isvalid(fig), return; end
if nargin<5, msg = ''; end
ids = identity; if ~iscell(ids), ids = {char(ids)}; end
if isempty(ids), return; end
app = getApp(fig);
for k = 1:numel(ids)
    key = cellKey(ids{k},stepId);
    if isempty(state)
        if isKey(app.runState,key), remove(app.runState,key); end
        if isKey(app.cellMsg,key),  remove(app.cellMsg,key);  end
    else
        app.runState(key) = state;
        if ~isempty(msg), app.cellMsg(key) = msg; end
    end
end
setApp(fig,app);
if any(strcmp(state,{'done','error'}))
    for k = 1:numel(ids), recordCompletion(fig,ids{k},stepId,state,msg); end
end
refreshCells(fig,ids);
refreshProgressDetail(fig);
% PLAIN drawnow, NOT limitrate (spec D1).  limitrate DROPS an update issued less
% than ~50 ms after the previous one, and the run path no longer has anything
% high-frequency to throttle: the progress axis is gone, and what is left is one
% update per call.  'queued' -> setRunningUI -> 'running' -> the wrapper's first
% line all land inside one 50 ms window, so limitrate dropped exactly the three
% updates that mattered and the monitor stayed a step behind for the whole run.
drawnow;
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
%   downstream done cells to stale, add whatever report images it wrote to the
%   result list - and tell the column it lost an entry, so the LAST one triggers
%   that column's PDF (spec D9).
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
setApp(fig,app);
files = addResults(fig,identity,stepId,model);
notePdfColumn(fig,stepId,files);
autosaveSession(fig);
end

function finishRun(fig, poolTok)
%finishRun  Reset the run UI, fold successful cells into the disk baseline, and
%   keep only error overlays (errors are not an on-disk state).
if nargin<2, poolTok = []; end
if ~isvalid(fig)
    % The window is already gone, so there is no UI to restore and no PDF to
    % assemble - but a parallel pool is a PROCESS-level resource, not a figure's, so
    % the workers still have to be handed back.
    wbPool('close', poolTok);
    return
end
app = getApp(fig);
wasCancel = isfield(app,'cancel') && app.cancel;
app.running = false; app.cancel = false;
setApp(fig,app);
% any column that never reached its last entry - a Stop, an error, a run that ended
% mid-column - is assembled here, so a PDF is never lost because the run was cut
% short.  It happens BEFORE the repaint below, which is what puts it in the list.
flushPdfColumns(fig);
% The workers go now: AFTER the PDFs (assembling a column must not race a pool
% shutdown) and BEFORE the repaint, since there is no reason to hold ~18 GB of idle
% worker heap across it.  finishRun is runChecked's onCleanup target, so a Stop, an
% error and a mid-run Exit all reach this line.  A pool the author opened himself is
% left alone - wbPool decides, from the token taken before the run.
releasePool(fig, poolTok);
recomputeBase(fig);
app = getApp(fig);
ks = keys(app.runState);
for i = 1:numel(ks)
    if ~strcmp(app.runState(ks{i}),'error'), remove(app.runState,ks{i}); end
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
drawnow;
% Exit pressed mid-run: the stop it asked for has now happened and the session is
% written, so this is the first safe moment to close (see requestExit).
if getfieldOr(getApp(fig),'exitAfterStop',false), onClose(fig); end
end

function releasePool(fig, poolTok)
%releasePool  Hand the run's parallel workers back, and SAY SO.
%   A silent fix to an invisible leak is indistinguishable from no fix - the 18 GB
%   never showed up in the GUI, so its disappearance has to be stated.  The line
%   names the worker count and NOT a byte figure: reading a pool's real footprint
%   means walking the OS process table from inside MATLAB, and an unmeasured number
%   is worse than no number.
n = wbPool('close', poolTok);
if n>0
    wbLog(fig, sprintf('released %d parallel %s.', n, plural(n,'worker')));
end
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

%% ===================== wipe (start the protocol over) ============== %%
function files = wipeTargets(app)
%wipeTargets  Every PROCESSED file of the listed recordings - the _d/_r/_s.mat
%   triplet members of each file on the Files tab, whatever stage/product flags sit
%   between the identity and the role.
%
%   THE IDENTITY IS THE FENCE, and it is matched as a whole rather than globbed:
%   'Foo' must not carry away 'Foo2_t_K_r.mat', so the directory listing is filtered
%   by ^<identity>(_...)?_[drs]\.mat$ .  A raw recording (.rls/.avi/...) can never
%   match that shape, which is why the raw data is safe by construction rather than
%   by a list of extensions to spare.
files = cell(1,0);
seen = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(app.files)
    m = app.files(i).model;
    if isempty(m.resultsFolder) || ~isfolder(m.resultsFolder), continue; end
    base = [m.roiPrefix m.stem];
    rx = ['^' regexptranslate('escape',base) '(_.*)?_[drs]\.mat$'];
    d = dir(fullfile(m.resultsFolder,[base '*.mat']));
    for k = 1:numel(d)
        if d(k).isdir || isempty(regexp(d(k).name,rx,'once')), continue; end
        p = fullfile(m.resultsFolder,d(k).name);
        if isKey(seen,p), continue; end
        seen(p) = true; files{end+1} = p; %#ok<AGROW>
    end
end
end

function uiWipeAll(fig)
%uiWipeAll  'Wipe all': throw away the processing and keep the recordings.  A
%   protocol whose settings turned out wrong is otherwise reprocessed one deleted
%   folder at a time, which is exactly when the wrong thing gets deleted - so the
%   list is derived from the curated file set (wipeTargets), COUNTED, and shown
%   before anything goes.  It refuses outright while a run is in flight.
app = getApp(fig);
if isfield(app,'running') && app.running
    alert(fig,'A run is in progress - stop it before wiping.'); return
end
if isempty(app.files)
    setLog(fig,{'Wipe all: no files are listed - nothing to wipe.'}); return
end
victims = wipeTargets(app);
if isempty(victims)
    setLog(fig,{'Wipe all: the listed recordings have no processed files on disk.'}); return
end
msg = sprintf(['Delete %d processed %s produced from the %d %s listed on the Files tab?\n\n' ...
    'This removes every processed file of those recordings - including any ' ...
    'listed .mat product itself.  The raw recordings are NOT touched.\n\n' ...
    'It cannot be undone.'], numel(victims), plural(numel(victims),'file'), ...
    numel(app.files), plural(numel(app.files),'recording'));
if isvalid(fig) && strcmp(fig.Visible,'on')
    sel = uiconfirm(fig,msg,'Wipe all','Options',{'Delete','Cancel'}, ...
        'DefaultOption',2,'CancelOption',2,'Icon','warning');
    if ~strcmp(sel,'Delete'), wbLog(fig,'Wipe all: cancelled.'); return; end
end
wipeAll(fig,victims);
end

function [nGone,failed] = wipeAll(fig,victims)
%wipeAll  Do the deletion and put the workbench back to "nothing has run": the
%   completion record, the staleness overlay and the run overlay all described files
%   that no longer exist, so they go with them and the state is re-read from disk.
app = getApp(fig);
if nargin<2 || isempty(victims), victims = wipeTargets(app); end
nGone = 0; failed = cell(1,0);
for i = 1:numel(victims)
    try
        delete(victims{i});                       % a FILE path, not a handle
        if isfile(victims{i}), failed{end+1} = victims{i}; else, nGone = nGone + 1; end %#ok<AGROW>
    catch
        failed{end+1} = victims{i}; %#ok<AGROW>
    end
end
app = getApp(fig);
app.completed = containers.Map('KeyType','char','ValueType','any');
app.stale     = containers.Map('KeyType','char','ValueType','any');
app.runState  = containers.Map('KeyType','char','ValueType','any');
app.cellMsg   = containers.Map('KeyType','char','ValueType','any');
setApp(fig,app);
recomputeBase(fig);
renderProgress(fig);
renderConstructor(fig);        % a row's flag can change once its products are gone
lines = {sprintf('Wipe all: deleted %d processed file(s).',nGone)};
if ~isempty(failed)
    lines{end+1} = sprintf('  %d could not be deleted (open elsewhere?):',numel(failed));
    for i = 1:min(10,numel(failed)), lines{end+1} = ['    ' failed{i}]; end %#ok<AGROW>
end
setLog(fig,lines);
autosaveSession(fig);
end

function setRunningUI(fig,running)
%setRunningUI  Toggle the Run/Preview/Stop buttons for a run in progress.
if ~isvalid(fig), return; end
c = getApp(fig).c.process;
onoff = ternary(running,'off','on');
if isfield(c,'runBtn')     && isgraphics(c.runBtn),     c.runBtn.Enable     = onoff; end
if isfield(c,'previewBtn') && isgraphics(c.previewBtn), c.previewBtn.Enable = onoff; end
if isfield(c,'fromDrop')   && isgraphics(c.fromDrop),   c.fromDrop.Enable   = onoff; end
if isfield(c,'toDrop')     && isgraphics(c.toDrop),     c.toDrop.Enable     = onoff; end
if isfield(c,'stopBtn')    && isgraphics(c.stopBtn),    c.stopBtn.Enable    = ternary(running,'on','off'); end
if isfield(c,'wipeBtn')    && isgraphics(c.wipeBtn),    c.wipeBtn.Enable    = onoff; end
drawnow;
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
function files = addResults(fig,identity,stepId,model)
%addResults  Fold a finished step's report images into the result list (spec D5).
%   A LIST, NOT THUMBNAILS.  The old panel drew a live uiimage per report of the
%   cell you happened to click; across a 200-file run that is both unaffordable
%   and useless - you want the last thing that finished, and you want it in a real
%   viewer.  So the run accumulates links, newest first, and a double-click hands
%   the file to the desktop (openArtifactViewer), which keeps it zoomable while
%   the next step is already running.
%
%   It RETURNS what it resolved, whether or not the entry was new, so the per-column
%   PDF is assembled from the same answer the list was painted from - wbArtifacts
%   stays the ONE place that knows how a report image is named.
%
%   AND THE SAME FENCE THE RUN USES (item 8).  wbArtifacts lists by base name, so a
%   '_rep_*.jpg' an earlier session left beside this recording matched a tail like
%   any other and joined both the list and the per-column PDF.  Handing it this
%   recording's admissible flags is the whole fix - the PDF inherits it for free,
%   since it is assembled from exactly these files.
files = cell(1,0);
if ~isvalid(fig), return; end
app = getApp(fig);
step = stepById(app.reg,stepId);
if nargin<4 || isempty(model), model = modelByIdentity(fig,identity); end
if isempty(model) || isempty(step), return; end
files = wbArtifacts(model,step,admissibleFlagsFor(app,identity));
lbl = @(p) resultLabel(shortId(identity), step.label, p);
if ~pushResults(fig, files, stepId, lbl), return; end
refreshResultList(fig);
end

function added = pushResults(fig, files, stepId, labelFcn)
%pushResults  Append paths to the result STORE (no repaint), skipping what is
%   already there.  One door for both kinds of result - the images a step wrote and
%   the PDF a column was assembled into - so an entry can never exist without its
%   stepId and its kind.
added = false;
if isempty(files), return; end
app = getApp(fig);
have = {};
if ~isempty(app.results), have = {app.results.path}; end
stamp = datetime('now');
for i = 1:numel(files)
    if any(strcmp(files{i},have)), continue; end
    app.results(end+1) = struct('path',files{i},'label',labelFcn(files{i}), ...
        'when',stamp,'stepId',char(stepId),'kind',artifactKind(files{i}));
    have{end+1} = files{i}; %#ok<AGROW>
    added = true;
end
if added, setApp(fig,app); end
end

function r = shownResults(app)
%shownResults  The result entries the panel currently shows: the whole store,
%   newest first, minus what the kind filter hides.  The store itself is never
%   touched - a filter is a view (spec D10).
r = app.results;
if isempty(r), return; end
[~,ord] = sort([r.when],'descend');
r = r(ord);
switch app.resultFilter
    case 'images', r = r(strcmp({r.kind},'image'));
    case 'pdfs',   r = r(strcmp({r.kind},'pdf'));
end
end

function setResultFilter(fig,kind)
%setResultFilter  Pick what the list shows ('all' | 'images' | 'pdfs').
kind = char(kind);
if ~any(strcmp(kind,resultKindIds())), kind = 'all'; end
app = getApp(fig); app.resultFilter = kind; setApp(fig,app);
c = app.c.process;
if isfield(c,'kindDrop') && isgraphics(c.kindDrop), c.kindDrop.Value = kind; end
refreshResultList(fig);
end

function refreshResultList(fig)
%refreshResultList  Repaint the link list, newest first, through the kind filter.
if ~isvalid(fig), return; end
app = getApp(fig); c = app.c.process;
if ~isfield(c,'resultList') || ~isgraphics(c.resultList), return; end
r = shownResults(app);
if isempty(r)
    c.resultList.Items = {emptyResultItem()};
    c.resultList.ItemsData = {''};
    return
end
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

function txt = derivedCellOf(fig,identity,branch,stepId)
%derivedCellOf  The bottom table's word for one ((recording,product), step) - the
%   headless view of the derived monitor, so a test can assert what a watcher reads
%   for a product that does not exist on disk yet (spec D7).
txt = '';
app = getApp(fig);
key = [char(identity) '||' char(branch)];
if ~isKey(app.dRowOf,key), return; end
r = app.dRows(app.dRowOf(key));
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
%resultLinks  The result list as {label, path} - the programmatic view of D5, and
%   the SHOWN one: it answers with what the panel is displaying, filter and all.
r = shownResults(getApp(fig));
L = cell(0,2);
if isempty(r), return; end
L = [{r.label}', {r.path}'];
end

function k = resultKinds(fig)
%resultKinds  The kind of every SHOWN entry, parallel to resultLinks - the headless
%   view of the filter.
r = shownResults(getApp(fig));
k = cell(1,0);
if ~isempty(r), k = {r.kind}; end
end

function openArtifactViewer(fig, pth)
%openArtifactViewer  Open one report full-size, preferring the DESKTOP's own
%   viewer (winopen on Windows, open/xdg-open elsewhere).  That viewer is a
%   separate process, so the report stays resizable, zoomable and scrollable while
%   the workbench is busy inside a wrapper - a MATLAB window would only repaint
%   when the main thread yields, and uiimage offers no zoom at all.
%
%   A PDF takes exactly the same road (spec D10): the hand-off is by file
%   association and does not care what the file is.  Only the in-MATLAB FALLBACK
%   is image-only - so when the desktop refuses a PDF the log says so, rather than
%   imread throwing behind an empty figure.
if isempty(pth) || ~isfile(pth), return; end
if openInDesktopViewer(pth), return; end
if strcmp(artifactKind(pth),'pdf')
    wbLog(fig,sprintf('  cannot open %s here - no PDF viewer is associated with it', pth));
    return
end
showArtifactInMatlab(pth);
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
%   the axes toolbar's zoom/pan are available (uiimage has neither).  It is only
%   ever handed an IMAGE - openArtifactViewer diverts a PDF to the log (D10).
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

%% ===================== per-column PDF reports (spec D9) ============= %%
function setPdfReports(fig,tf)
%setPdfReports  The 'Create PDF reports' switch.  A WINDOW PREFERENCE: it changes
%   nothing about what runs, what is written by a wrapper or what the session
%   carries (spec D9 - no schema bump), only whether a column's report images are
%   also collected into one document when the column ends.
app = getApp(fig); app.pdfReports = logical(tf); setApp(fig,app);
c = app.c.process;
if isfield(c,'pdfCheck') && isgraphics(c.pdfCheck), c.pdfCheck.Value = app.pdfReports; end
% deleting the images is a consequence of assembling them, so it is only offerable
% while there is a document to assemble
if isfield(c,'delCheck') && isgraphics(c.delCheck)
    c.delCheck.Enable = ternary(app.pdfReports,'on','off');
end
end

function setDeleteAfterPdf(fig,tf)
%setDeleteAfterPdf  The 'Delete the images once the PDF is written' switch (Q3).
%   A WINDOW PREFERENCE like the PDF switch itself, and OFF by default: every entry
%   in the result list is an image resolved from disk (wbArtifacts), so removing
%   the images removes the links, and the PDF is then the only copy of a page there
%   is - nothing rebuilds a deleted JPG except running the step again.  That is
%   what makes it opt-in.  Nothing is ever deleted by a WRAPPER - wbExecutor pins
%   s.reportKeepImages true - so this switch is the only door.
app = getApp(fig); app.deleteAfterPdf = logical(tf); setApp(fig,app);
c = app.c.process;
if isfield(c,'delCheck') && isgraphics(c.delCheck), c.delCheck.Value = app.deleteAfterPdf; end
end

function beginPdfColumns(fig, entries)
%beginPdfColumns  Stash, once, how many entries each COLUMN of THIS run has.
%
%   THE COUNT MUST COME FROM THE LIST THE RUN IS EXECUTING.  The run order is
%   step-major, so a column ends when its last entry lands - but since Phase 7 the
%   list is a From/To SLICE, and re-expanding it mid-run would answer differently
%   as steps complete on disk, so the column would either never reach its end or
%   reach it twice.  It is therefore counted here from the entries runChecked
%   already holds, and never asked for again.  A column outside the range simply
%   never appears in the map, and so never fires.
app = getApp(fig);
app.colTotal = containers.Map('KeyType','char','ValueType','double');
app.colDone  = containers.Map('KeyType','char','ValueType','double');
app.colArt   = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(entries)
    id = entries(i).stepId;
    if isKey(app.colTotal,id)
        app.colTotal(id) = app.colTotal(id) + 1;
    else
        app.colTotal(id) = 1; app.colDone(id) = 0; app.colArt(id) = cell(1,0);
    end
end
setApp(fig,app);
end

function notePdfColumn(fig, stepId, files)
%notePdfColumn  One entry of a column finished: keep its report images and tick the
%   column's counter; when the last entry lands, assemble.  A report is a
%   BY-PRODUCT - nothing in here may reach the run loop, so the whole thing is
%   defended and a failure costs a log line, never the step's state.
try
    app = getApp(fig);
    if ~isKey(app.colTotal,stepId), return; end     % not a column of this run
    if ~isempty(files)
        app.colArt(stepId) = uniqueStable([app.colArt(stepId), files(:)']);
    end
    app.colDone(stepId) = app.colDone(stepId) + 1;
    ended = app.colDone(stepId) >= app.colTotal(stepId);
    setApp(fig,app);
    if ended, flushPdfColumn(fig, stepId); end
catch ME
    wbLog(fig,['  PDF report bookkeeping failed: ' ME.message]);
end
end

function flushPdfColumns(fig)
%flushPdfColumns  Assemble every column still pending - what finishRun calls, so a
%   Stop, an error or a run that ended mid-column still leaves the PDFs it earned.
%   It runs BEFORE the tables are repainted and does not depend on them.
app = getApp(fig);
if ~isa(app.colArt,'containers.Map') || app.colArt.Count==0, return; end
for id = orderIds(app.reg, keys(app.colArt))
    flushPdfColumn(fig, id{1});
end
end

function flushPdfColumn(fig, stepId)
%flushPdfColumn  ONE column -> ONE PDF, then the column is forgotten (so it can
%   never be written twice).  The document joins the result list like any other
%   artifact and the log records what was written or why it was not.
app = getApp(fig);
if ~isKey(app.colArt,stepId), return; end
files = app.colArt(stepId);
remove(app.colArt,stepId);
if isKey(app.colTotal,stepId), remove(app.colTotal,stepId); end
if isKey(app.colDone,stepId),  remove(app.colDone,stepId);  end
setApp(fig,app);
if ~app.pdfReports || isempty(files), return; end

step = stepById(app.reg,stepId);
lbl  = char(stepId); if ~isempty(step), lbl = step.label; end
try
    root = pdfReportRoot(app, files);
    if isempty(root)
        wbLog(fig,sprintf('  no PDF for %s - nowhere to put it', lbl));
        return
    end
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    out = makeReportPdf(files, fullfile(root,'workbenchReports',[char(stepId) '_' stamp '.pdf']));
    for i = 1:numel(out.skipped)
        wbLog(fig,sprintf('  skipped %s (not a readable report image)', out.skipped{i}));
    end
    if isempty(out.path)
        wbLog(fig,sprintf('  no PDF for %s - none of its %d report image(s) could be read', ...
            lbl, numel(files)));
        return
    end
    wbLog(fig,sprintf('  wrote %s (%d page(s) of %s)', out.path, out.pages, lbl));
    app = getApp(fig); app.reportPdfs{end+1} = out.path; setApp(fig,app);
    % 'PDF : Contrast : 18:15:03' - the kind first because that is what the reader
    % is scanning for, and the clock because a column can be re-run in one session
    % and two documents of the same step must be told apart.
    pdfLbl = strjoin({'PDF', lbl, stampClock(stamp)}, sep());
    pushResults(fig, {out.path}, stepId, @(~) pdfLbl);
    dropSourceImages(fig, files, out, lbl);
    refreshResultList(fig);
catch ME
    wbLog(fig,sprintf('  PDF report for %s failed: %s', lbl, ME.message));
end
end

function dropSourceImages(fig, files, out, lbl)
%dropSourceImages  The optional cleanup (Q3), and it is TWO removals, not one: the
%   file goes from the disk AND its entry goes from the result list, in the same
%   step.  Leaving the entry behind would hand the user a link to a file that is no
%   longer there - and since the list is repainted from disk anyway, a stale entry
%   would not even survive the next repaint.
%
%   Only images that really reached a page are removed: makeReportPdf names what it
%   could not read, and an image that was skipped is the one image the PDF cannot
%   stand in for.  Like every other part of the report path this is a by-product -
%   a failure costs a log line and never the column's state.
if isempty(files) || out.pages<=0, return; end
app = getApp(fig);
if ~app.deleteAfterPdf, return; end
gone = cell(1,0);
for i = 1:numel(files)
    f = files{i};
    if any(strcmp(f,out.skipped)), continue; end        % not in the document
    try
        if isfile(f), delete(f); gone{end+1} = f; end %#ok<AGROW>
    catch ME
        wbLog(fig,sprintf('  could not remove %s (%s)', shortName(f), ME.message));
    end
end
if isempty(gone), return; end
dropResults(fig, gone);
wbLog(fig,sprintf('  removed %d report image(s) of %s - they are in the PDF', ...
    numel(gone), lbl));
end

function dropResults(fig, paths)
%dropResults  Take entries out of the result STORE by path (no repaint).  The only
%   remover there is: an entry exists because its file does, so the two are dropped
%   together (see dropSourceImages).
if isempty(paths), return; end
app = getApp(fig);
if isempty(app.results), return; end
keep = ~ismember({app.results.path}, paths);
if all(keep), return; end
app.results = app.results(keep);
setApp(fig,app);
end

function root = pdfReportRoot(app, files)
%pdfReportRoot  Where a column's PDF lands: THE RESULTS ROOT when there is one, so
%   a project keeps all of its reports together and its recordings untouched;
%   otherwise the folder of that column's first artifact - which is what keeps a
%   working set assembled by hand out of several folders writing beside the images
%   each report is made of.  A document is a result like any other, so it goes
%   where the pages it is made of went.
root = '';
r = resultsDirOf(app);
if ~isempty(r) && isfolder(r), root = r; return; end
if ~isempty(files), root = fileparts(files{1}); end
end

function p = reportPdfList(fig)
%reportPdfList  The per-column PDFs assembled this session (headless view).
p = getApp(fig).reportPdfs;
end

function s = shortId(identity)
[~,s] = fileparts(char(identity));
if isempty(s), s = char(identity); end
end
function s = shortName(p)
[~,n,e] = fileparts(char(p));
s = [n e];
end

function lbl = resultLabel(recording, stepLabel, pth)
%resultLabel  How one report reads in the list: the recording, the step, the page.
%
%   It used to read 'Mouse1 - Segmentation - Mouse1_t_K_rep_segments.jpg', which
%   says the recording twice and then spells out a file name whose every token the
%   line has already given.  What distinguishes one entry from another is the
%   RECORDING and the PAGE, so those are what is left: 'Mouse1 : Segmentation :
%   segments'.  The full path stays in ItemsData, so a double-click still opens
%   exactly the file it always did.
lbl = strjoin({char(recording), char(stepLabel), reportTail(pth)}, sep());
end

function t = reportTail(pth)
%reportTail  The page's own name: whatever follows '_rep_' in a report written by
%   Core/Reporting ('Mouse1_t_K_rep_segments.jpg' -> 'segments').  A page from
%   before the rename has no such marker, so its last name token is used instead
%   ('Mouse1_t_K_cm.jpg' -> 'cm') - cryptic, but it is what that file is called.
[~,n] = fileparts(char(pth));
k = strfind(n,'_rep_');
if ~isempty(k)
    t = n(k(end)+5:end);
    return
end
p = strfind(n,'_');
if isempty(p), t = n; else, t = n(p(end)+1:end); end
end

function t = stampClock(stamp)
%stampClock  'yyyymmdd_HHMMSS' -> 'HH:MM:SS', the only half of it a reader needs
%   inside one session.  Anything unexpected is handed back untouched.
t = char(stamp);
if numel(t)==15 && t(9)=='_'
    t = [t(10:11) ':' t(12:13) ':' t(14:15)];
end
end

function s = sep()
%sep  The separator between the parts of a result label.  A middle dot, built from
%   its code point so the source file stays plain ASCII (as progressLegend does).
s = [' ' char(183) ' '];
end

%% ===================== log ========================================= %%
function wbLog(fig,msg)
%wbLog  Append one line and SHOW it.  Plain drawnow, not limitrate (spec D1): the
%   whole log of a run is three lines per recording (Core/Reporting emits Starting,
%   Writing results, Finished), which is nowhere near limitrate's ~50 ms threshold -
%   so all it could ever do here is drop the "Starting ..." that arrives in the same
%   instant as the cell going 'running', which is precisely the line the operator is
%   waiting for while a wrapper blocks the thread for the next few minutes.
app = getApp(fig); c = app.c.process;
v = c.log.Value; if ischar(v), v = {v}; end
v{end+1} = msg;
if numel(v)>400, v = v(end-400:end); end
c.log.Value = v; drawnow;
end
function setLog(fig,lines)
app = getApp(fig);
if ischar(lines), lines = {lines}; end
if numel(lines)>400, lines = lines(end-400:end); end
app.c.process.log.Value = lines; drawnow;
end
function v = getLog(fig)
v = getApp(fig).c.process.log.Value; if ischar(v), v = {v}; end
end

%% ===================== exit ======================================== %%
function requestExit(fig)
%requestExit  Leaving is ONE decision, not three.  Exit stops the run itself - the
%   user should never have to press Stop first, and pressing it twice was only ever
%   a way of saying "yes, I meant it" - saves the session to the last session, and
%   closes.  A run cannot be torn down mid-step, so while one is in flight this
%   requests the same cooperative stop the Stop button does and hands the close to
%   finishRun, which is the moment the executor has actually let go.
if ~isvalid(fig), return; end
app = getApp(fig);
if isfield(app,'running') && app.running
    app.exitAfterStop = true; setApp(fig,app);
    cancelRun(fig);
    wbLog(fig,'Exit: stopping after the current step - the window closes when it does.');
    return
end
autosaveSession(fig);
onClose(fig);
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
req = wbPrereqs('describe',step);
if ~isempty(req), t = [t '.  It needs ' req ' first.']; end
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
