%guiMyograph  Unified GUI for myograph diameter, vasomotion & propagation
%
%   guiMyograph opens one programmatic uifigure app that replaces the
%   Launcher_Myograph script and pickMyographIntervals as the user-facing
%   entry points.  It is a thin layer over the headless functions
%   (getMyographDiameter, cutMyographIntervals, getMyographVasomotion,
%   getMyographPropagation, runMyographFile, loadMyographResults,
%   exportMyographExcel) - every action it performs is also callable from code.
%
%   Five stages (tabs), in workflow order:
%     Setup      - build the file queue, edit the s parameters (grouped as in the
%                  launcher), per-file time crop / row range / pixel size / fps
%                  override / predefined intervals, and "re-detect per interval".
%     Run        - measure the DIAMETER over the whole recording (or crop) only,
%                  batch, per-file/per-stage progress, cancellable, failure isolation.
%     Intervals  - load that diameter, then draw / edit / name / validate intervals
%                  on the median trace and Y-vs-time map (ported pickMyographIntervals,
%                  no inputdlg), with a side video preview of the selected interval
%                  showing the detected walls overlaid (green), exportable as a frame
%                  or an overlaid video.
%     Vasomotion - vasomotion + propagation analysis applied to every defined
%                  interval (or the whole recording if none were defined).
%     Explore    - load one or many *_myograph.mat and plot diameter / spectra /
%                  wavelet spectrograms / vasomotion markers / propagation, per file or as
%                  group comparisons; export figures and an Excel workbook.
%
%   h = guiMyograph('Visible','off') creates the app without showing it and
%   returns the figure; getappdata(h,'workbenchAPI') exposes a struct of function
%   handles (addFiles, addFolder, setParams, setFileOption, setFileIntervals,
%   runDiameter/runDiameterAll, loadDiameter, addInterval, saveIntervals,
%   previewInterval, exportPreviewFrame/Video, runVasomotion/runVasomotionAll,
%   exploreLoad, render, exportFigure, exportExcel, getApp) so the whole workflow
%   can be driven and tested headlessly.
%
%   DEPENDS ON  the Core/Myograph/ headless functions above; base MATLAB + uifigure
%   (R2020b+).  Self-contained.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

function h = guiMyograph(varargin)

vis='on';
for i=1:2:numel(varargin)
    if strcmpi(varargin{i},'Visible'), vis=varargin{i+1}; end
end

app=struct();
app.queue=emptyQueue();          % struct array of files + per-file settings
app.s=launcherDefaults();        % global s parameters
app.cancel=false;
app.explore=struct('path','','intervals',[],'s',struct(),'meta',struct());  % single-file "focus"
app.exploreFiles=struct('path',{},'name',{},'group',{},'rec',{},'recnum',{},'resultPath',{});  % multi-file/group model
app.exploreFocus=1;              % index into exploreFiles used for the single-file plots
app.exploreRoot='';              % folder last scanned (for Scan/apply)
app.exploreGroupOverride=containers.Map('KeyType','char','ValueType','char'); % manual groups (path->name), survive re-scan
app.exploreCache=containers.Map('KeyType','char','ValueType','any');          % LRU cache of loaded {intervals,s,meta}
app.exploreOrder={}; app.exploreCacheLimit=8;
app.detIntervals=[];             % Intervals-tab: current file's whole-recording diameter
app.roiList={}; app.selRoi=[];

% single instance: replace any stale window (by Tag AND by Name, so windows from an
% older code version - which lack the Tag - are also cleared)
delete(findall(groot,'Type','figure','Tag','guiMyograph'));
delete(findall(groot,'Type','figure','Name','Myograph workbench'));
try, delete(timerfind('Name','myographPreviewTimer')); catch, end  % clear any leaked preview timer
fig=uifigure('Name','Myograph workbench','Position',[60 60 1400 860],'Visible',vis,'Tag','guiMyograph');
app.fig=fig;
g=uigridlayout(fig,[1 1],'Padding',[4 4 4 4]);
tg=uitabgroup(g);
app.tabs.setup     =uitab(tg,'Title','1 - File');
app.tabs.run       =uitab(tg,'Title','2 - Diameter');
app.tabs.intervals =uitab(tg,'Title','3 - Intervals');
app.tabs.vasomotion=uitab(tg,'Title','4 - Vasomotion');
app.tabs.explore   =uitab(tg,'Title','5 - Explore');

setappdata(fig,'app',app);
buildSetup(fig); buildRun(fig); buildIntervals(fig); buildVasomotion(fig); buildExplore(fig);

% ---- programmatic API (drives the same internal logic as the UI) ----
api=struct('addFiles',@(p)addFiles(fig,p), 'addFolder',@(d,r)addFolder(fig,d,r), ...
    'setParams',@(s)setParams(fig,s), 'getParams',@()getApp(fig).s, ...
    'setFileOption',@(i,f,v)setFileOption(fig,i,f,v), ...
    'setFileIntervals',@(i,t,n)setFileIntervals(fig,i,t,n), ...
    'runDiameter',@(i)runStage(fig,i,'diameter'), 'runDiameterAll',@()runList(fig,1:numel(getApp(fig).queue),'diameter'), ...
    'runVasomotion',@(i)runStage(fig,i,'analysis'), 'runVasomotionAll',@()runList(fig,1:numel(getApp(fig).queue),'analysis'), ...
    'loadDiameter',@(i)loadDiameterForIntervals(fig,i), 'showIntervals',@()refreshIntervalFileDrop(fig), ...
    'detectForIntervals',@(i)detectForIntervals(fig,i), ...
    'addInterval',@(t0,t1,nm)apiAddInterval(fig,t0,t1,nm), ...
    'getIntervals',@()collectIntervals(fig), 'saveIntervals',@()saveIntervalsToFile(fig), ...
    'previewInterval',@(a,b)setupPreview(fig,a,b), 'exportPreviewFrame',@(p)exportPreviewFrame(fig,p), ...
    'exportPreviewVideo',@(p)exportPreviewVideo(fig,p), ...
    'exploreLoad',@(p)exploreLoad(fig,p), 'render',@()exploreRender(fig), ...
    'exploreLoadGroup',@(paths,grpPat)exploreSetFiles(fig,paths,grpPat), ...
    'exploreSelect',@(idx)apiExploreSelect(fig,idx), ...
    'exploreCreateGroup',@(name,idx)apiExploreCreateGroup(fig,name,idx), ...
    'setExplore',@(field,value)apiSetExplore(fig,field,value), ...
    'exportFigure',@(p,dpi,fmt)exploreExportFigure(fig,p,dpi,fmt), ...
    'exportExcel',@(p)exploreExportExcel(fig,p), ...
    'stop',@()cancelRun(fig), 'exit',@()onCloseWorkbench(fig), 'getApp',@()getApp(fig));
setappdata(fig,'workbenchAPI',api);
refreshQueue(fig);
tg.SelectionChangedFcn=@(~,e)onTabChanged(fig,e);      % refresh the Intervals file list on entry
if nargout>0, h=fig; end
end

function onTabChanged(fig,e)
%ONTABCHANGED  when the Intervals tab is opened, list the processed files and show one
app=getApp(fig);
if isfield(app,'tabs') && isfield(app.tabs,'intervals') && e.NewValue==app.tabs.intervals
    refreshIntervalFileDrop(fig);
end
end

%% ===================== app-state helpers ============================ %%
function app=getApp(fig), app=getappdata(fig,'app'); end
function setApp(fig,app), setappdata(fig,'app',app); end
function q=emptyQueue()
q=struct('path',{},'name',{},'state',{},'timeCrop',{},'rowRange',{}, ...
    'pixelSize',{},'fps',{},'intervals',{},'intervalNames',{},'resultPath',{});
end
function s=launcherDefaults()
s=struct('edgeMode','center','tau',0.85,'rowRange',[1 Inf],'smoothSigma',1.2,'dustRadius',8, ...
    'tSpan',25,'ySpan',31,'outlierK',3,'subpixel',true,'smoothSpan',15, ...
    'vFR',[0.05 0.25],'cFR',[0.4 0.6],'wFR',[0.01 1],'wVPO',10,'normalisation','median', ...
    'normsize',inf,'tgtFS',1,'pcts',0:10:100,'otsuMaxN',5,'otsuElbow',0.05, ...
    'propNShuffle',200,'propSignal','diameter','detectPerInterval',false);
s.segVsmReturn={'bands','moments','series','clustering','reconstruction','spectrum'};   %analysis levels to compute/store ('moments' stores the per-frequency fVectors; 'spectrum' stores the spectrum.amp/.phase grid for the spectrogram)
end

%% ===================== SETUP tab ==================================== %%
function buildSetup(fig)
app=getApp(fig); t=app.tabs.setup;
gl=uigridlayout(t,[1 1],'Padding',[6 6 6 6]);

% -- file queue + per-file settings (analysis parameters live on the Diameter / Vasomotion tabs) --
lp=uigridlayout(gl,[6 1],'RowHeight',{'fit','1x','fit','fit','fit','fit'},'RowSpacing',6);
brow=uigridlayout(lp,[1 4],'Padding',[0 0 0 0]);
uibutton(brow,'Text','Add files','ButtonPushedFcn',@(~,~)uiAddFiles(fig));
uibutton(brow,'Text','Add folder','ButtonPushedFcn',@(~,~)uiAddFolder(fig));
uibutton(brow,'Text','Remove','ButtonPushedFcn',@(~,~)uiRemove(fig));
uibutton(brow,'Text','Clear','ButtonPushedFcn',@(~,~)uiClear(fig));
c.queueTbl=uitable(lp,'ColumnName',{'file','state','crop t0','crop t1','row lo','row hi','px µm','fps'}, ...
    'ColumnEditable',[false false true true true true true true],'CellEditCallback',@(s,e)onQueueEdit(fig,s,e), ...
    'CellSelectionCallback',@(s,e)onQueueSelect(fig,e));
uilabel(lp,'Text','Predefined intervals for the selected file (start, end, name) - optional, act as labels:');
c.ivTbl=uitable(lp,'ColumnName',{'start_s','end_s','name'},'ColumnEditable',[true true true], ...
    'Data',cell(0,3),'CellEditCallback',@(s,e)onIvTblEdit(fig));
ivb=uigridlayout(lp,[1 4],'Padding',[0 0 0 0]);
uibutton(ivb,'Text','Add row','ButtonPushedFcn',@(~,~)ivAddRow(fig));
uibutton(ivb,'Text','Delete row','ButtonPushedFcn',@(~,~)ivDelRow(fig));
uibutton(ivb,'Text','Import...','ButtonPushedFcn',@(~,~)ivImport(fig));
uibutton(ivb,'Text','Export...','ButtonPushedFcn',@(~,~)ivExport(fig));
sess=uigridlayout(lp,[1 3],'Padding',[0 0 0 0]);
uibutton(sess,'Text','Save session...','ButtonPushedFcn',@(~,~)saveSession(fig));
uibutton(sess,'Text','Load session...','ButtonPushedFcn',@(~,~)loadSession(fig));
uibutton(sess,'Text','Exit workbench','ButtonPushedFcn',@(~,~)requestExit(fig),'BackgroundColor',[1 0.82 0.82]);

app.c.setup=c; setApp(fig,app);
end

function g=diameterGroups()
g={ 'Wall detection', {'edgeMode','tau','rowRange'};
    'Robustness',     {'smoothSigma','dustRadius','tSpan','ySpan','outlierK','subpixel','smoothSpan'};
    'Detection mode', {'detectPerInterval'} };
end
function g=vasoGroups()
g={ 'Vasomotion bands',   {'vFR','cFR','wFR','wVPO'};
    'Vasomotion options', {'normalisation','normsize','tgtFS'};
    'Analysis levels',    {'segVsmReturn'};
    'Propagation',        {'propNShuffle','propSignal'} };
end
function ctrls=buildParamEditor(fig,parent,groups)
tips=paramTips();
stack=uigridlayout(parent,[size(groups,1) 1],'RowHeight',repmat({'fit'},1,size(groups,1)),'RowSpacing',8);
ctrls=struct();
app=getApp(fig);
for gi=1:size(groups,1)
    p=uipanel(stack,'Title',groups{gi,1},'FontWeight','bold');
    flds=groups{gi,2};
    gg=uigridlayout(p,[numel(flds) 2],'ColumnWidth',{'fit','1x'},'RowHeight',repmat({'fit'},1,numel(flds)));
    for fi=1:numel(flds)
        f=flds{fi};
        uilabel(gg,'Text',f,'Tooltip',getTip(tips,f));
        val=paramValue(app.s,f);
        if islogical(val)
            ctrls.(f)=uicheckbox(gg,'Value',logical(val),'Text','','ValueChangedFcn',@(s,~)onParamEdit(fig,f,s.Value));
        elseif strcmp(f,'edgeMode')
            ctrls.(f)=uidropdown(gg,'Items',{'center','min','outer','inner'},'Value',valOr(val,'center'),'ValueChangedFcn',@(s,~)onParamEdit(fig,f,s.Value));
        elseif strcmp(f,'normalisation')
            ctrls.(f)=uidropdown(gg,'Items',{'mean','median','mmean','mmedian'},'Value',valOr(val,'median'),'ValueChangedFcn',@(s,~)onParamEdit(fig,f,s.Value));
        elseif strcmp(f,'propSignal')
            ctrls.(f)=uidropdown(gg,'Items',{'diameter','rData'},'Value',valOr(val,'diameter'),'ValueChangedFcn',@(s,~)onParamEdit(fig,f,s.Value));
        elseif strcmp(f,'segVsmReturn')
            % analysis levels are a SET (s.segVsmReturn); render one checkbox per level
            levels={'bands','moments','series','clustering','reconstruction','spectrum'};
            if iscell(val), curVR=val; else, curVR={}; end
            lg=uigridlayout(gg,[numel(levels) 1],'Padding',[0 0 0 0],'RowSpacing',1,'RowHeight',repmat({'fit'},1,numel(levels)));
            for li=1:numel(levels)
                L=levels{li};
                ctrls.(['segVsmReturn_' L])=uicheckbox(lg,'Text',L,'Value',ismember(L,curVR), ...
                    'ValueChangedFcn',@(s,~)onVsmReturnEdit(fig,L,s.Value),'Tooltip',getTip(tips,['segVsmReturn_' L]));
            end
        else
            ctrls.(f)=uieditfield(gg,'text','Value',val2str(val),'ValueChangedFcn',@(s,~)onParamEdit(fig,f,s.Value),'Tooltip',getTip(tips,f));
        end
    end
end
end
function v=valOr(x,d), if isempty(x)||~(ischar(x)||isstring(x)), v=d; else, v=char(x); end, end

%% ===================== RUN tab ====================================== %%
function buildRun(fig)
app=getApp(fig); t=app.tabs.run;
gl=uigridlayout(t,[1 2],'ColumnWidth',{'1x','1.5x'},'Padding',[6 6 6 6],'ColumnSpacing',8);
% -- left: diameter parameters --
pp=uipanel(gl,'Scrollable','on','Title','Diameter parameters (applied to all files)');
pg=uigridlayout(pp,[1 1],'RowHeight',{'fit'},'Padding',[4 4 4 4],'Scrollable','on');  % the grid owns the scroll
buildParamEditor(fig,pg,diameterGroups());
% -- right: run controls --
rc=uigridlayout(gl,[5 1],'RowHeight',{'fit','fit','fit','1x','fit'},'RowSpacing',6);
uilabel(rc,'Text',['Step 2 - measure the wall-to-wall DIAMETER over the whole recording ' ...
    '(or the time crop). Diameter only; define intervals in tab 3, run vasomotion in tab 4.'],'WordWrap','on');
br=uigridlayout(rc,[1 4],'Padding',[0 0 0 0]);
c.runBtn=uibutton(br,'Text','Measure diameter (all)','ButtonPushedFcn',@(~,~)runList(fig,1:numel(getApp(fig).queue),'diameter'));
c.runSelBtn=uibutton(br,'Text','Measure diameter (selected)','ButtonPushedFcn',@(~,~)runSelectedStage(fig,'diameter'));
c.stopBtn=uibutton(br,'Text','Stop','ButtonPushedFcn',@(~,~)cancelRun(fig),'BackgroundColor',[1 0.9 0.7]);
c.exitBtn=uibutton(br,'Text','Exit','ButtonPushedFcn',@(~,~)requestExit(fig),'BackgroundColor',[1 0.82 0.82]);
c.prog=uilabel(rc,'Text','Idle.','FontWeight','bold');
c.log=uitextarea(rc,'Value',{''},'Editable','off');
c.resTbl=uitable(rc,'ColumnName',{'file','stage','status','intervals','detect s','vsm s','prop s','message'});
app=getApp(fig); app.c.run=c; setApp(fig,app);
end

function buildVasomotion(fig)
app=getApp(fig); t=app.tabs.vasomotion;
gl=uigridlayout(t,[1 2],'ColumnWidth',{'1x','1.5x'},'Padding',[6 6 6 6],'ColumnSpacing',8);
% -- left: vasomotion + propagation parameters --
pp=uipanel(gl,'Scrollable','on','Title','Vasomotion & propagation parameters');
pg=uigridlayout(pp,[1 1],'RowHeight',{'fit'},'Padding',[4 4 4 4],'Scrollable','on');  % the grid owns the scroll
buildParamEditor(fig,pg,vasoGroups());
% -- right: run controls --
rc=uigridlayout(gl,[5 1],'RowHeight',{'fit','fit','fit','1x','fit'},'RowSpacing',6);
uilabel(rc,'Text',['Step 4 - vasomotion + propagation, applied to every defined interval ' ...
    '(or the whole recording if none). Needs the diameter from tab 2; define/correct intervals in tab 3 first.'],'WordWrap','on');
br=uigridlayout(rc,[1 4],'Padding',[0 0 0 0]);
c.runBtn=uibutton(br,'Text','Analyse (all)','ButtonPushedFcn',@(~,~)runList(fig,1:numel(getApp(fig).queue),'analysis'));
c.runSelBtn=uibutton(br,'Text','Analyse (selected)','ButtonPushedFcn',@(~,~)runSelectedStage(fig,'analysis'));
c.stopBtn=uibutton(br,'Text','Stop','ButtonPushedFcn',@(~,~)cancelRun(fig),'BackgroundColor',[1 0.9 0.7]);
c.exitBtn=uibutton(br,'Text','Exit','ButtonPushedFcn',@(~,~)requestExit(fig),'BackgroundColor',[1 0.82 0.82]);
c.prog=uilabel(rc,'Text','Idle.','FontWeight','bold');
c.log=uitextarea(rc,'Value',{''},'Editable','off');
c.resTbl=uitable(rc,'ColumnName',{'file','stage','status','intervals','detect s','vsm s','prop s','message'});
app=getApp(fig); app.c.vaso=c; setApp(fig,app);
end

%% ===================== INTERVALS tab ================================ %%
function buildIntervals(fig)
app=getApp(fig); t=app.tabs.intervals;
gl=uigridlayout(t,[1 3],'ColumnWidth',{'1.9x','0.95x','1.15x'},'Padding',[6 6 6 6],'ColumnSpacing',8);

% -- col 1: editable diameter trace + Y-vs-time map --
lp=uigridlayout(gl,[2 1],'RowHeight',{'1x','1x'},'RowSpacing',6);
c.axTrace=uiaxes(lp); title(c.axTrace,'median-over-Y diameter'); xlabel(c.axTrace,'time (s)'); ylabel(c.axTrace,'diameter (px)');
c.axMap=uiaxes(lp); title(c.axMap,'Y vs time diameter map'); xlabel(c.axMap,'time (s)'); ylabel(c.axMap,'Y');

% -- col 2: interval controls (pick a processed file -> its diameter plots automatically) --
rp=uigridlayout(gl,[6 1],'RowHeight',{'fit','fit','fit','1x','fit','fit'},'RowSpacing',6);
c.fileDrop=uidropdown(rp,'Items',{'(no processed files)'},'ValueChangedFcn',@(~,~)onIntervalFileChange(fig));
uilabel(rp,'Text',['Pick a processed file above - its median-over-Y diameter plots below. ' ...
    'Drag on the trace to add an interval (name it, then Add); drag a band''s edge or nudge with ' ...
    'arrow keys; edit start/end/name directly in the table (precise, even when zoomed); ' ...
    'select a table row to preview it.'],'WordWrap','on');
nmrow=uigridlayout(rp,[1 3],'Padding',[0 0 0 0]);
c.ivName=uieditfield(nmrow,'text','Value','');
uibutton(nmrow,'Text','Add','ButtonPushedFcn',@(~,~)uiAddIntervalFromBox(fig));
c.delToggle=uibutton(nmrow,'Text','Delete: off','ButtonPushedFcn',@(s,~)toggleDelete(fig,s));
c.ivListTbl=uitable(rp,'ColumnName',{'start_s','end_s','name'},'ColumnEditable',[true true true], ...
    'CellSelectionCallback',@(s,e)onIvPreviewSelect(fig,e),'CellEditCallback',@(s,e)onIvTableEdit(fig,e));
uibutton(rp,'Text','Validate','ButtonPushedFcn',@(~,~)validateIntervalsUI(fig));
uibutton(rp,'Text','Save intervals to file','ButtonPushedFcn',@(~,~)saveIntervalsToFile(fig));

% -- col 3: video preview with detected walls overlaid --
pp=uigridlayout(gl,[6 1],'RowHeight',{'fit','1x','fit','fit','fit','fit'},'RowSpacing',5);
uilabel(pp,'Text','Interval preview - video with detected walls (green). Select an interval row to load it.','WordWrap','on');
c.previewAx=uiaxes(pp); set(c.previewAx,'XTick',[],'YTick',[]); title(c.previewAx,'(no interval selected)');
c.previewSlider=uislider(pp,'Limits',[1 2],'Value',1,'ValueChangedFcn',@(s,~)onPreviewSlide(fig,s.Value));
prow=uigridlayout(pp,[1 2],'Padding',[0 0 0 0]);
uibutton(prow,'Text','Play','ButtonPushedFcn',@(~,~)previewPlay(fig));
uibutton(prow,'Text','Stop','ButtonPushedFcn',@(~,~)previewStop(fig));
srow=uigridlayout(pp,[1 2],'ColumnWidth',{'fit','1x'},'Padding',[0 0 0 0]);
uilabel(srow,'Text','Frame step (skip):');
c.previewStep=uispinner(srow,'Limits',[1 200],'Value',2,'Step',1,'RoundFractionalValues','on', ...
    'Tooltip','frames advanced per playback tick - higher = faster, skips more frames');
xrow=uigridlayout(pp,[1 2],'Padding',[0 0 0 0]);
uibutton(xrow,'Text','Export frame...','ButtonPushedFcn',@(~,~)uiExportPreviewFrame(fig));
uibutton(xrow,'Text','Export video...','ButtonPushedFcn',@(~,~)uiExportPreviewVideo(fig));

app.fig.WindowKeyPressFcn=@(~,e)onIntervalKey(fig,e);   % arrow-key nudge of the selected band
app.fig.CloseRequestFcn=@(~,~)requestExit(fig);        % window-X stops a run first, then closes
app.c.intervals=c; app.preview=struct('vr',[],'videoPath','','times',[],'curIdx',1,'timer',[]);
setApp(fig,app);
end

function onIntervalKey(fig,e)
%ONINTERVALKEY  nudge/resize the currently selected interval band with the arrow keys
app=getApp(fig); sel=[];
for i=1:numel(app.roiList)
    r=app.roiList{i};
    if isvalid(r) && isprop(r,'Selected') && r.Selected, sel=r; break; end
end
if isempty(sel), return; end
step=1; if any(strcmp(e.Modifier,'control')), step=5; end     % snap unit = 1 s (5 s with Ctrl)
p=sel.Position; resize=any(strcmp(e.Modifier,'shift'));
switch e.Key
    case 'leftarrow',  if resize, p(3)=max(step,p(3)-step); else, p(1)=p(1)-step; end
    case 'rightarrow', if resize, p(3)=p(3)+step;           else, p(1)=p(1)+step; end
    otherwise, return;
end
p(1)=round(p(1)); p(3)=round(p(3));                           % snap to whole seconds
sel.Position=p; refreshIvList(fig);
end

%% ===================== EXPLORE tab ================================== %%
function buildExplore(fig)
app=getApp(fig); t=app.tabs.explore;
gl=uigridlayout(t,[1 2],'ColumnWidth',{'3x','1x'},'Padding',[6 6 6 6],'ColumnSpacing',8);
c.ax=uiaxes(gl); title(c.ax,'Load one or more *_myograph.mat to begin');

% scrollable, sectioned control column
rpPanel=uipanel(gl,'BorderType','none','Scrollable','on');
stack=uigridlayout(rpPanel,[4 1],'RowHeight',repmat({'fit'},1,4),'RowSpacing',8,'Padding',[2 2 2 2],'Scrollable','on');  % the grid owns the scroll

% --- section 1: data source (files + groups) ---
s1=exSection(stack,'1 - Data source (one file, or many + groups)',12);
b1=uigridlayout(s1,[1 2],'Padding',[0 0 0 0]); b1.Layout.Row=1; b1.Layout.Column=[1 2];
uibutton(b1,'Text','Load file(s)...','ButtonPushedFcn',@(~,~)uiExploreLoad(fig));
uibutton(b1,'Text','Load folder...','ButtonPushedFcn',@(~,~)uiExploreLoadFolder(fig));
c.srcLbl=uilabel(s1,'Text','(nothing loaded)','WordWrap','on'); c.srcLbl.Layout.Row=2; c.srcLbl.Layout.Column=[1 2];
c.includeF=exField(s1,3,'Include','_myograph\.mat$');
c.groupF=exField(s1,4,'Group','');
c.includeF.Tooltip=sprintf(['Regexp matched against each file NAME; only matching files are kept.\n' ...
    'Default keeps result files.']);
c.groupF.Tooltip=sprintf(['Regexp whose MATCH labels each file''s experimental group\n' ...
    '(e.g. WT|KO). Files sharing a match form one group. Empty = one group.\n' ...
    'Or: select files below, type a name and press "Create group".']);
c.scanBtn=exBtnFull(s1,5,'Scan / apply',@(~,~)exploreScan(fig));
c.fileList=uilistbox(s1,'Items',{},'Multiselect','on','ValueChangedFcn',@(~,~)onExploreFileSel(fig), ...
    'Tooltip','Ctrl/Shift-click to select files, then "Create group" to group them (overrides the Group pattern).');
c.fileList.Layout.Row=[6 9]; c.fileList.Layout.Column=[1 2];
c.groupNameF=exField(s1,10,'New group','');
c.createBtn=exBtnFull(s1,11,'Create group from selected',@(~,~)exploreCreateGroup(fig));
c.matchLbl=uilabel(s1,'Text','','WordWrap','on','FontColor',[0.5 0.5 0.5]); c.matchLbl.Layout.Row=12; c.matchLbl.Layout.Column=[1 2];

% --- section 2: what to plot ---
s2=exSection(stack,'2 - What to plot',6);
c.plotType=exDrop(s2,1,'Plot',explorePlotTypes(),@(~,~)onExplorePlotType(fig));
c.interval=exDrop(s2,2,'Interval',{'(all)'},@(~,~)exploreRender(fig));
c.vsmMarker=exDrop(s2,3,'vsmMarker',{'ampMeanVB','ampStdVB','ampSkewVB','ampMeanCB','ampStdCB','ampSkewCB', ...
    'fCentMeanVB','fCentStdVB','fSprdMeanVB','fSprdStdVB','shapePeakVB','nPeakMeanVB','nPeakStdVB', ...
    'durFlareMeanVB','durFlareStdVB','durSilenceMeanVB','durSilenceStdVB','ampFlareMeanVB','ampFlareStdVB','ampSilenceMeanVB','ampSilenceStdVB', ...
    'durFlareMeanCB','durFlareStdCB','durSilenceMeanCB','durSilenceStdCB','ampFlareMeanCB','ampFlareStdCB','ampSilenceMeanCB','ampSilenceStdCB'},@(~,~)exploreRender(fig));
c.vsmMarker.Tooltip='Per-Y vasomotion scalar for "group: vsmMarker (box)". The band is the name suffix VB (s.vFR) / CB (s.cFR); read from scalars.<band>.* in the results tree. ampMeanVB/CB = per-Y band-mean amplitude.';
c.propMetric=exDrop(s2,4,'Prop. scalar',{'speed','confidence','R2','domFreq','phaseSpeed'},@(~,~)exploreRender(fig));
c.propMetric.Tooltip='Scalar used by "group: propagation (box)" - one value per file x interval.';
hHint=uilabel(s2,'Text','"group: ..." plots compare the selected files; single-file plots use the focus file.', ...
    'FontColor',[0.5 0.5 0.5],'WordWrap','on','FontSize',10);
hHint.Layout.Row=5; hHint.Layout.Column=[1 2];
c.plotBtn=exBtnFull(s2,6,'Plot / refresh',@(~,~)exploreRender(fig));

% --- section 3: comparison (used by "group: ..." plots) ---
s3=exSection(stack,'3 - Group comparison',3);
c.organize=exDrop(s3,1,'Organise by',{'Auto','Group','Interval','Group x Interval','Pool all'},@(~,~)exploreRender(fig));
c.points=exDrop(s3,2,'Points are',{'Files (mean over Y)','Y-locations'},@(~,~)exploreRender(fig));
c.stat=exDrop(s3,3,'Centre / error',{'Mean +/- SD','Mean +/- SEM','Median +/- IQR'},@(~,~)exploreRender(fig));

% --- section 4: export ---
s4=exSection(stack,'4 - Export',5);
c.dpi=exField(s4,1,'DPI','300');
c.fmt=exDrop(s4,2,'Format',{'PNG','TIFF','PDF','EPS','JPEG'},[]);
c.exportFigBtn=exBtnFull(s4,3,'Export figure...',@(~,~)uiExploreExportFig(fig));
c.exportXlsBtn=exBtnFull(s4,4,'Export Excel...',@(~,~)uiExploreExportXlsx(fig));
c.status=uilabel(s4,'Text','','WordWrap','on'); c.status.Layout.Row=5; c.status.Layout.Column=[1 2];

app.c.explore=c; setApp(fig,app); updateExploreEnable(fig);   % grey out controls that don't apply to the default plot
end
function items=explorePlotTypes()
%EXPLOREPLOTTYPES  single-file types, then the group-comparison variants
items={'median diameter','per-Y diameter','per-Y normalised diameter','normalised diameter', ...
    'Y-vs-time map','spectrum vs f (avg)','percentile spectra','amplitude percentiles','wavelet spectrogram', ...
    'propagation: lag vs Y','propagation: Y-vs-time+fit', ...
    'group: diameter','group: normalised diameter','group: spectrum','group: wavelet spectrogram', ...
    'group: vsmMarker (box)','group: propagation (box)'};
end

%% ===================== queue management ============================= %%
function addFiles(fig,paths)
if ischar(paths)||isstring(paths), paths=cellstr(paths); end
app=getApp(fig);
for i=1:numel(paths)
    app.queue(end+1)=newQueueEntry(paths{i}); %#ok<AGROW>
end
setApp(fig,app); refreshQueue(fig);
end
function addFolder(fig,folder,recursive)
if nargin<3, recursive=true; end
pat=ternary(recursive,fullfile(folder,'**','*.*'),fullfile(folder,'*.*'));
d=dir(pat); d=d(~[d.isdir]);
keep=false(numel(d),1); exts={'.avi','.mp4','.mov','.mj2','.mkv','.wmv'};
for i=1:numel(d), [~,~,e]=fileparts(d(i).name); keep(i)=any(strcmpi(e,exts)); end
d=d(keep); paths=arrayfun(@(x)fullfile(x.folder,x.name),d,'uni',0);
addFiles(fig,paths);
end
function e=newQueueEntry(p)
[~,nn,ex]=fileparts(p);
e=struct('path',p,'name',[nn ex],'state','raw','timeCrop',[],'rowRange',[1 Inf], ...
    'pixelSize',[],'fps',[],'intervals',[],'intervalNames',{{}},'resultPath','');
[pp,nn2]=fileparts(p); rp=fullfile(pp,[nn2 '_myograph.mat']);
if exist(rp,'file'), e.state='result exists'; e.resultPath=rp; end
end
function refreshQueue(fig)
app=getApp(fig); q=app.queue; n=numel(q);
D=cell(n,8);
for i=1:n
    D(i,:)={q(i).name,q(i).state,cropv(q(i).timeCrop,1),cropv(q(i).timeCrop,2), ...
        q(i).rowRange(1),rr2(q(i).rowRange),pxv(q(i).pixelSize),fpsv(q(i).fps)};
end
app.c.setup.queueTbl.Data=D;
setApp(fig,app);
end

function onQueueSelect(fig,e)
if isempty(e.Indices), return; end
app=getApp(fig); app.selFile=e.Indices(1); setApp(fig,app);
renderIvTable(fig,e.Indices(1));
end
function onQueueEdit(fig,~,e)
app=getApp(fig); r=e.Indices(1); col=e.Indices(2); v=e.NewData; q=app.queue(r);
switch col
    case 3, q.timeCrop=setcrop(q.timeCrop,1,v);
    case 4, q.timeCrop=setcrop(q.timeCrop,2,v);
    case 5, q.rowRange(1)=num(v,1);
    case 6, q.rowRange(2)=num(v,Inf);
    case 7, q.pixelSize=numEmpty(v);
    case 8, q.fps=numEmpty(v);
end
app.queue(r)=q; setApp(fig,app);
end
function setFileOption(fig,idx,field,value)
app=getApp(fig); app.queue(idx).(field)=value; setApp(fig,app); refreshQueue(fig);
end
function setFileIntervals(fig,idx,ivT,names)
app=getApp(fig); app.queue(idx).intervals=ivT; app.queue(idx).intervalNames=names;
if ~strcmp(app.queue(idx).state,'result exists'), app.queue(idx).state='intervals set'; end
setApp(fig,app); refreshQueue(fig); renderIvTable(fig,idx);
end
function renderIvTable(fig,idx)
app=getApp(fig); q=app.queue(idx); ivT=q.intervals; nm=q.intervalNames;
D=cell(size(ivT,1),3);
for k=1:size(ivT,1), D(k,:)={ivT(k,1),ivT(k,2),nmk(nm,k)}; end
app.c.setup.ivTbl.Data=D; setApp(fig,app);
end

%% ===================== parameters =================================== %%
function setParams(fig,s)
app=getApp(fig); fn=fieldnames(s);
for i=1:numel(fn), app.s.(fn{i})=s.(fn{i}); end
setApp(fig,app);
end
function onParamEdit(fig,field,rawval)
app=getApp(fig);
if islogical(app.s.(field)) || ismember(field,{'subpixel','detectPerInterval'})
    app.s.(field)=logical(rawval);
elseif (ischar(rawval)||isstring(rawval)) && ismember(field,{'edgeMode','normalisation','propSignal'})
    app.s.(field)=char(rawval);
else
    app.s.(field)=str2paramValue(rawval);
end
setApp(fig,app);
end
function onVsmReturnEdit(fig,level,tf)
%toggle one analysis level in s.segVsmReturn (order-independent; the shared core
%getVasomotionMetrics treats segVsmReturn as a set)
app=getApp(fig);
if isfield(app.s,'segVsmReturn') && iscell(app.s.segVsmReturn), vr=app.s.segVsmReturn; else, vr={}; end
if tf
    if ~ismember(level,vr), vr{end+1}=level; end
else
    vr=vr(~strcmp(vr,level));
end
app.s.segVsmReturn=vr; setApp(fig,app);
end

%% ===================== RUN / VASOMOTION logic ======================= %%
function out=runStage(fig,idx,stage), outs=runList(fig,idx,stage); if isempty(outs), out=[]; else, out=outs{1}; end, end
function runSelectedStage(fig,stage)
app=getApp(fig); idx=1; if isfield(app,'selFile')&&~isempty(app.selFile), idx=app.selFile; end
runList(fig,idx,stage);
end
function outs=runList(fig,idxs,stage)
if nargin<3, stage='diameter'; end
app=getApp(fig); app.cancel=false; app.running=true; app.runStage=stage; setApp(fig,app);
res=cell(0,8); outs=cell(1,numel(idxs));
for ii=1:numel(idxs)
    if ~isvalid(fig), return; end
    app=getApp(fig); if app.cancel, logStage(fig,stage,'Stopped.'); break; end
    i=idxs(ii); q=app.queue(i);
    setProgStage(fig,stage,sprintf('File %d/%d: %s',ii,numel(idxs),q.name));
    s=fileParams(app.s,q);
    s.cancelFcn=@()safeCancel(fig);                        % stop mid-file (checked at progress points)
    s.stageFcn=@(st,det)logStage(fig,stage,sprintf('  [%s] %s',st,det));
    s.progressFcn=@(f,tot)setProgStage(fig,stage,sprintf('%s: frame %d/%d',q.name,f,tot));
    if strcmp(stage,'diameter')
        s.runVasomotion=false; s.runPropagation=false;
        if ~(s.detectPerInterval && ~isempty(q.intervals)), s.intervals=[]; s.intervalNames={}; end % Path B: one interval
        out=runMyographFile(s,q.path);
    else
        out=analyzeStage(fig,q,s);
    end
    if ~isvalid(fig), return; end
    outs{ii}=out;
    app=getApp(fig);
    app.queue(i).state=stageState(stage,out);
    if strcmp(out.status,'done')&&~isempty(out.path), app.queue(i).resultPath=out.path; end
    setApp(fig,app);
    res(end+1,:)={q.name,stage,out.status,out.nIntervals,rnd(out.timings,'detect'), ...
        rnd(out.timings,'vasomotion'),rnd(out.timings,'propagation'),out.message}; %#ok<AGROW>
    cc=stageCtrls(fig,stage); cc.resTbl.Data=res; refreshQueue(fig);
    if strcmp(out.status,'cancelled'), break; end
end
if isvalid(fig)
    app=getApp(fig); app.running=false; wasCancel=app.cancel; setApp(fig,app);
    outs=outs(~cellfun(@isempty,outs));
    setProgStage(fig,stage,ternary(wasCancel,'Stopped.','Done.'));
end
end
function tf=safeCancel(fig)
tf=false; try, if isvalid(fig), tf=getApp(fig).cancel; end, catch, end
end
function out=analyzeStage(fig,q,s)
out=struct('status','failed','message','','path',q.resultPath,'nIntervals',0,'meta',struct(), ...
    'timings',struct('detect',0,'vasomotion',0,'propagation',0),'intervals',[]);
try
    rp=q.resultPath;
    if isempty(rp)||~exist(rp,'file'), error('No diameter result - run tab 2 (Measure diameter) first.'); end
    [intervals,sSaved,meta]=loadMyographResults(rp);
    if numel(intervals)==1 && ~isempty(q.intervals)      % cut a whole-recording result by the defined intervals
        intervals=cutMyographIntervals(intervals(1),q.intervals,q.intervalNames);
    end
    sA=s; sA.runVasomotion=true; sA.runPropagation=true;
    [intervals,tim]=analyzeMyographIntervals(sA,intervals,s.stageFcn);
    s=sSaved; save(rp,'intervals','s','meta','-v7.3');
    out.status='done'; out.message='ok'; out.path=rp; out.nIntervals=numel(intervals);
    out.timings.vasomotion=tim.vasomotion; out.timings.propagation=tim.propagation;
    out.intervals=intervals; out.meta=meta;
catch ME
    if contains(lower(ME.identifier),'cancel')||contains(lower(ME.message),'stopped by user')
        out.status='cancelled'; out.message='Stopped by user.';
    else
        out.message=ME.message;
    end
end
end
function st=stageState(stage,out)
switch out.status
    case 'done',      st=ternary(strcmp(stage,'diameter'),'diameter done','analysed');
    case 'cancelled', st='stopped';
    otherwise,        st='FAILED';
end
end
function c=stageCtrls(fig,stage)
app=getApp(fig); if strcmp(stage,'analysis'), c=app.c.vaso; else, c=app.c.run; end
end
function s=fileParams(sGlobal,q)
s=sGlobal;
if ~isempty(q.timeCrop), s.timeCrop=q.timeCrop; end
if ~isempty(q.rowRange), s.rowRange=q.rowRange; end
if ~isempty(q.pixelSize), s.pixelSize=q.pixelSize; end
if ~isempty(q.intervals), s.intervals=q.intervals; s.intervalNames=q.intervalNames; else, s.intervals=[]; s.intervalNames={}; end
end
function cancelRun(fig)
app=getApp(fig); app.cancel=true; setApp(fig,app);
% update ONLY the running stage's own label - touching the other tab's
% components would pull that tab to the front.
if isfield(app,'running')&&app.running&&isfield(app,'runStage')
    try, c=stageCtrls(fig,app.runStage); c.prog.Text='Stopping after the current step...'; catch, end
end
end
function requestExit(fig)
%REQUESTEXIT  stop a running batch first (click again to close); close cleanly when idle
if ~isvalid(fig), return; end
app=getApp(fig);
if isfield(app,'running') && app.running
    app.cancel=true; setApp(fig,app);
    if isfield(app,'runStage')
        try, c=stageCtrls(fig,app.runStage); c.prog.Text='Stopping... click Exit again to close.'; catch, end
    end
else
    onCloseWorkbench(fig);
end
end
function setProgStage(fig,stage,msg), c=stageCtrls(fig,stage); c.prog.Text=msg; drawnow limitrate; end
function logStage(fig,stage,msg)
c=stageCtrls(fig,stage); v=c.log.Value; if ischar(v), v={v}; end
v{end+1}=msg; if numel(v)>400, v=v(end-400:end); end
c.log.Value=v; drawnow limitrate;
end

%% ===================== INTERVALS logic ============================= %%
function detectForIntervals(fig,idx)
app=getApp(fig); q=app.queue(idx);
s=fileParams(app.s,q); s.intervals=[]; s.intervalNames={};
s.progressFcn=@(f,tot)setProgStage(fig,'diameter',sprintf('re-measuring frame %d/%d',f,tot));
iv=getMyographDiameter(s,q.path);
setDiameterForEditing(fig,idx,iv(1));
end
function loadDiameterForIntervals(fig,idx)
q=getApp(fig).queue(idx); rp=q.resultPath;
if isempty(rp)||~exist(rp,'file'), error('No diameter result for %s - run tab 2 (Measure diameter) first.',q.name); end
[intervals,~,~]=loadMyographResults(rp);
setDiameterForEditing(fig,idx,stitchIntervals(intervals));
q=getApp(fig).queue(idx);                                   % seed bands from stored defs / existing cuts
if ~isempty(q.intervals)
    for k=1:size(q.intervals,1), addCommittedBand(fig,q.intervals(k,1),q.intervals(k,2),nmk(q.intervalNames,k)); end
elseif numel(intervals)>1
    for k=1:numel(intervals), tt=double(intervals(k).time); addCommittedBand(fig,tt(1),tt(end),intervals(k).name); end
end
end
function refreshIntervalFileDrop(fig)
%REFRESHINTERVALFILEDROP  list only files that HAVE a diameter result; auto-show one
app=getApp(fig); proc=[]; names={};
for i=1:numel(app.queue)
    rp=app.queue(i).resultPath;
    if ~isempty(rp) && exist(rp,'file'), proc(end+1)=i; names{end+1}=app.queue(i).name; end %#ok<AGROW>
end
dd=app.c.intervals.fileDrop;
if isempty(proc)
    dd.Items={'(no processed files)'}; dd.ItemsData=[]; clearIntervalEditing(fig); return;
end
dd.Items=names; dd.ItemsData=proc;
if isfield(app,'detFile') && ~isempty(app.detFile) && ismember(app.detFile,proc)
    dd.Value=app.detFile;                                  % keep the current file if still valid
else
    dd.Value=proc(1);
end
onIntervalFileChange(fig);
end
function onIntervalFileChange(fig)
app=getApp(fig); dd=app.c.intervals.fileDrop;
if isempty(dd.ItemsData), return; end
idx=dd.Value; if isempty(idx), return; end
if isfield(app,'detFile') && isequal(app.detFile,idx) && ~isempty(app.detIntervals), return; end % already shown
try, loadDiameterForIntervals(fig,idx); catch ME, clearIntervalEditing(fig); fprintf('[Intervals] %s\n',ME.message); end
end
function clearIntervalEditing(fig)
previewStop(fig); app=getApp(fig);
app.detIntervals=[]; app.detFile=[]; app.roiList={}; app.selRoi=[];
if isfield(app,'preview'), app.preview.vr=[]; app.preview.times=[]; end
setApp(fig,app); c=app.c.intervals;
try, cla(c.axTrace,'reset'); title(c.axTrace,'(no processed file selected)'); catch, end
try, cla(c.axMap,'reset'); catch, end
try, cla(c.previewAx,'reset'); set(c.previewAx,'XTick',[],'YTick',[]); title(c.previewAx,'(no interval selected)'); catch, end
try, c.ivListTbl.Data=cell(0,3); catch, end
end
function onIvTableEdit(fig,e)
%ONIVTABLEEDIT  edit an interval's start/end/name in the table (precise, zoom-independent)
if isempty(e.Indices), return; end
app=getApp(fig); r=e.Indices(1); col=e.Indices(2);
if ~(isfield(app,'ivListRois') && numel(app.ivListRois)>=r && isvalid(app.ivListRois{r})), return; end
roi=app.ivListRois{r}; p=roi.Position;
switch col
    case 3, roi.Label=char(string(e.NewData));
    case 1                                          % new start (keep end fixed)
        a=num(e.NewData,NaN); b=p(1)+p(3);
        if isfinite(a) && b>a, roi.Position=[a p(2) b-a p(4)]; end
    case 2                                          % new end (keep start fixed)
        b=num(e.NewData,NaN);
        if isfinite(b) && b>p(1), roi.Position=[p(1) p(2) b-p(1) p(4)]; end
end
refreshIvList(fig);
end
function setDiameterForEditing(fig,idx,big)
previewStop(fig);
app=getApp(fig); app.detIntervals=big; app.detFile=idx; app.roiList={};
try, app.preview.vr=VideoReader(app.queue(idx).path); catch, app.preview.vr=[]; end
app.preview.videoPath=app.queue(idx).path; app.preview.times=[]; app.preview.curIdx=1;
setApp(fig,app); drawIntervalTrace(fig);
end
function big=stitchIntervals(intervals)
%STITCHINTERVALS  reconstruct one continuous diameter from a set of interval slices
if numel(intervals)==1, big=intervals(1); if ~isfield(big,'name')||isempty(big.name), big.name='full'; end, return; end
t=[]; L=[]; R=[]; D=[]; M=[]; V=[];
for k=1:numel(intervals)
    t=[t; double(intervals(k).time(:))]; %#ok<AGROW>
    L=[L; intervals(k).idxL]; R=[R; intervals(k).idxR]; D=[D; intervals(k).diameter]; %#ok<AGROW>
    nf=size(intervals(k).diameter,1);
    if isfield(intervals,'mask')&&~isempty(intervals(k).mask), M=[M; intervals(k).mask]; else, M=[M; ones(nf,size(intervals(k).diameter,2),'single')]; end %#ok<AGROW>
    if isfield(intervals,'valid')&&~isempty(intervals(k).valid), V=[V; intervals(k).valid(:)]; else, V=[V; true(nf,1)]; end %#ok<AGROW>
end
[tu,iu]=unique(t);                                          % sorted, de-duplicated in time
big=struct('name','full','time',tu,'idxL',L(iu,:),'idxR',R(iu,:),'diameter',D(iu,:),'mask',M(iu,:),'valid',V(iu));
end
function drawIntervalTrace(fig)
app=getApp(fig); iv=app.detIntervals; if isempty(iv), return; end
t=double(iv.time); D=double(iv.diameter); med=median(D,2,'omitnan');
if isfield(iv,'valid') && ~isempty(iv.valid), valid=logical(iv.valid(:)); else, valid=true(numel(t),1); end
ax=app.c.intervals.axTrace; cla(ax); hold(ax,'on');
plot(ax,t,med,'-','Color',[0 0.3 0.8]);                    % valid diameter (blue)
if any(~valid)                                             % off-FOV frames (a wall dilated out of view): lower bound, in RED
    mr=med; mr(valid)=NaN; plot(ax,t,mr,'-','Color',[0.85 0.1 0.1],'LineWidth',1.6);
    title(ax,'median-over-Y diameter  (red = wall off-FOV, invalid)');
else
    title(ax,'median-over-Y diameter (drag to add interval)');
end
grid(ax,'on'); hold(ax,'off');
xlabel(ax,'time (s)'); ylabel(ax,'diameter (px)');
if numel(t)>1, xlim(ax,[t(1) t(end)]); end
app.ivYL=ylim(ax);
axm=app.c.intervals.axMap; cla(axm);
dec=max(1,round(size(D,1)/2000));
imagesc(axm,t(1:dec:end),1:size(D,2),D(1:dec:end,:)'); set(axm,'YDir','normal');
xlabel(axm,'time (s)'); ylabel(axm,'Y'); title(axm,'Y vs time diameter map'); colorbar(axm);
% selection rectangle for adding
yl=app.ivYL; sp=max(t(end)-t(1),eps);
app.selRoi=drawrectangle(ax,'Position',[t(1)+0.02*sp yl(1) 0.1*sp diff(yl)],'Color',[0.95 0.75 0],'FaceAlpha',0.12,'Label','new');
addlistener(app.selRoi,'ROIMoved',@(s,~)snapY(s,yl));
setApp(fig,app);
end
function addCommittedBand(fig,x1,x2,name)
app=getApp(fig); ax=app.c.intervals.axTrace; yl=app.ivYL;
roi=drawrectangle(ax,'Position',[min(x1,x2) yl(1) abs(x2-x1) diff(yl)],'Color',[0.2 0.6 1], ...
    'FaceAlpha',0.12,'Label',char(name),'LabelVisible','on');
addlistener(roi,'ROIMoved',@(s,~)snapY(s,yl));
addlistener(roi,'ROIClicked',@(s,e)onBandClicked(fig,s,e));
app.roiList{end+1}=roi; setApp(fig,app); refreshIvList(fig);
end
function apiAddInterval(fig,t0,t1,nm), addCommittedBand(fig,t0,t1,nm); end
function uiAddIntervalFromBox(fig)
app=getApp(fig); if isempty(app.selRoi)||~isvalid(app.selRoi), return; end
p=app.selRoi.Position; nm=strtrim(app.c.intervals.ivName.Value);
if isempty(nm), nm=sprintf('interval%d',numel(app.roiList)+1); end
addCommittedBand(fig,p(1),p(1)+p(3),nm);             % commits the band (setApp's the new state)
app=getApp(fig); app.c.intervals.ivName.Value='';    % re-fetch the fresh state after the commit
% place the selection box just to the RIGHT of the new band, same duration, kept in-axes
if ~isempty(app.selRoi) && isvalid(app.selRoi)
    ax=app.c.intervals.axTrace; xl=xlim(ax); w=max(p(3),eps); newX=p(1)+p(3);
    if newX+w>xl(2), newX=max(xl(1),xl(2)-w); end     % if it would overflow, tuck it against the right edge
    np=app.selRoi.Position; np(1)=newX; np(3)=w; app.selRoi.Position=np;
end
setApp(fig,app);
end
function toggleDelete(fig,src)
app=getApp(fig); on=~isfield(app,'delMode')||~app.delMode; app.delMode=on;
src.Text=ternary(on,'Delete: ON','Delete: off'); setApp(fig,app);
end
function onBandClicked(fig,src,evt) %#ok<INUSD>
% delete mode removes the band on click; renaming is done by editing the table.
app=getApp(fig);
if isfield(app,'delMode') && app.delMode
    idx=find(cellfun(@(r)isvalid(r)&&r==src,app.roiList)); app.roiList(idx)=[]; delete(src); setApp(fig,app); refreshIvList(fig);
end
end
function refreshIvList(fig)
[ivT,nm,rois]=collectIntervals(fig); app=getApp(fig);
D=cell(size(ivT,1),3); for k=1:size(ivT,1), D(k,:)={ivT(k,1),ivT(k,2),nm{k}}; end
app.c.intervals.ivListTbl.Data=D; app.ivListRois=rois; setApp(fig,app);
end
function [ivT,names,rois]=collectIntervals(fig)
app=getApp(fig); ivT=zeros(0,2); names={}; rois={};
for i=1:numel(app.roiList)
    r=app.roiList{i};
    if isvalid(r), p=r.Position; ivT(end+1,:)=[p(1) p(1)+p(3)]; names{end+1}=char(r.Label); rois{end+1}=r; end %#ok<AGROW>
end
if ~isempty(ivT), [ivT,o]=sortrows(ivT); names=names(o); rois=rois(o); end
end
function ok=validateIntervalsUI(fig)
[ivT,names]=collectIntervals(fig); app=getApp(fig); iv=app.detIntervals;
[err,note]=validateIntervals(ivT,names,iv);
ok=isempty(err);
if ~ok, uialert(fig,err,'Invalid intervals');
elseif ~isempty(note), uialert(fig,['Valid. Note: ' note],'Validate','Icon','warning');
else, uialert(fig,'Intervals valid.','Validate','Icon','success'); end
end
function saveIntervalsToFile(fig)
app=getApp(fig);
if ~isfield(app,'detFile')||isempty(app.detFile)||isempty(app.detIntervals)
    alertUser(fig,'Select a processed file in the dropdown first.','Save intervals'); return;
end
[ivT,names]=collectIntervals(fig);
[err,~]=validateIntervals(ivT,names,app.detIntervals);
if ~isempty(err), alertUser(fig,err,'Cannot save'); return; end
idx=app.detFile; setFileIntervals(fig,idx,ivT,names);
rp=getApp(fig).queue(idx).resultPath;
if isempty(rp)||~exist(rp,'file')
    alertUser(fig,'Interval definitions stored. Run tab 2 to create the diameter result, then tab 4 to analyse.','Saved','info'); return;
end
big=app.detIntervals; intervals=cutMyographIntervals(big,ivT,names); %#ok<NASGU>
[~,sSaved,meta]=loadMyographResults(rp); s=sSaved; %#ok<NASGU>
save(rp,'intervals','s','meta','-v7.3');
alertUser(fig,sprintf('%d interval(s) saved to %s.',size(ivT,1),getApp(fig).queue(idx).name),'Saved','success');
end
function alertUser(fig,msg,titl,icon)
%ALERTUSER  uialert when the window is shown; console fallback (headless/tests)
if nargin<4, icon='info'; end
try
    if isprop(fig,'Visible') && strcmp(fig.Visible,'on'), uialert(fig,msg,titl,'Icon',icon);
    else, fprintf('[%s] %s\n',titl,msg); end
catch
    fprintf('[%s] %s\n',titl,msg);
end
end

%% ===================== INTERVAL VIDEO PREVIEW ======================= %%
function onIvPreviewSelect(fig,e)
if isempty(e.Indices), return; end
app=getApp(fig); r=e.Indices(1); D=app.c.intervals.ivListTbl.Data;
if r>size(D,1), return; end
t0=num(D{r,1},NaN); t1=num(D{r,2},NaN);
if isfinite(t0)&&isfinite(t1)&&t1>t0, setupPreview(fig,t0,t1); end
end
function setupPreview(fig,t0,t1)
previewStop(fig); app=getApp(fig);
if ~isfield(app,'preview')||isempty(app.preview.vr), return; end
fps=app.preview.vr.FrameRate;
app.preview.t0=t0; app.preview.t1=t1; app.preview.times=(t0:1/fps:t1)';
if isempty(app.preview.times), app.preview.times=t0; end
app.preview.curIdx=1;
n=max(2,numel(app.preview.times));
app.c.intervals.previewSlider.Limits=[1 n]; app.c.intervals.previewSlider.Value=1;
setApp(fig,app); renderPreviewFrame(fig);
end
function onPreviewSlide(fig,val)
app=getApp(fig); if isempty(app.preview.times), return; end
app.preview.curIdx=max(1,min(numel(app.preview.times),round(val))); setApp(fig,app);
renderPreviewFrame(fig);
end
function renderPreviewFrame(fig)
if ~isvalid(fig), return; end
app=getApp(fig); pv=app.preview; ax=app.c.intervals.previewAx;
if ~isvalid(ax) || isempty(pv.vr) || isempty(pv.times), return; end   % never draw into a dead axes
t=pv.times(min(pv.curIdx,numel(pv.times))); big=app.detIntervals;
try, pv.vr.CurrentTime=min(max(t,0),max(0,pv.vr.Duration-1/pv.vr.FrameRate)); fr=readFrame(pv.vr); catch, return; end
if ndims(fr)<3, fr=repmat(fr,1,1,3); end
cla(ax); image(ax,fr); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); hold(ax,'on');
[~,fi]=min(abs(double(big.time)-t)); nY=size(big.idxL,2);
plot(ax,big.idxL(fi,:),1:nY,'g-','LineWidth',1.2);
plot(ax,big.idxR(fi,:),1:nY,'g-','LineWidth',1.2);
title(ax,sprintf('t=%.1f s  (frame %d/%d)',t,pv.curIdx,numel(pv.times))); hold(ax,'off');
end
function previewPlay(fig)
app=getApp(fig);
if isempty(app.preview.vr)||isempty(app.preview.times), return; end
if isfield(app.preview,'timer')&&~isempty(app.preview.timer)&&isvalid(app.preview.timer), return; end
tm=timer('ExecutionMode','fixedRate','Period',0.08,'BusyMode','drop','Name','myographPreviewTimer','TimerFcn',@(~,~)previewTick(fig));
app.preview.timer=tm; setApp(fig,app); start(tm);
end
function previewStop(fig)
app=getApp(fig);
if isfield(app,'preview')&&isfield(app.preview,'timer')&&~isempty(app.preview.timer)&&isvalid(app.preview.timer)
    stop(app.preview.timer); delete(app.preview.timer);
end
if isfield(app,'preview'), app.preview.timer=[]; setApp(fig,app); end
end
function previewTick(fig)
if ~isvalid(fig), try, delete(timerfind('Name','myographPreviewTimer')); catch, end, return; end
app=getApp(fig); if isempty(app.preview.times), previewStop(fig); return; end
step=2; try, step=max(1,round(app.c.intervals.previewStep.Value)); catch, end   % frame-skip control
n=numel(app.preview.times); idx=app.preview.curIdx+step; if idx>n, idx=1; end
app.preview.curIdx=idx; setApp(fig,app);
try, app.c.intervals.previewSlider.Value=min(idx,app.c.intervals.previewSlider.Limits(2)); catch, end
renderPreviewFrame(fig);
end
function uiExportPreviewFrame(fig)
[f,p]=uiputfile({'*.png';'*.tif';'*.jpg'},'Export frame'); if isequal(f,0), return; end
try, exportPreviewFrame(fig,fullfile(p,f)); catch ME, uialert(fig,ME.message,'Export frame'); end
end
function exportPreviewFrame(fig,path)
app=getApp(fig); pv=app.preview; big=app.detIntervals;
if isempty(pv.vr)||isempty(pv.times), error('No interval selected for preview.'); end
t=pv.times(min(pv.curIdx,numel(pv.times)));
pv.vr.CurrentTime=min(max(t,0),max(0,pv.vr.Duration-1/pv.vr.FrameRate)); fr=readFrame(pv.vr);
if ndims(fr)<3, fr=repmat(fr,1,1,3); end
[~,fi]=min(abs(double(big.time)-t));
imwrite(drawWalls(im2uint8(fr),big.idxL(fi,:),big.idxR(fi,:)),path);
end
function uiExportPreviewVideo(fig)
[f,p]=uiputfile('*.avi','Export overlaid video'); if isequal(f,0), return; end
try, exportPreviewVideo(fig,fullfile(p,f)); uialert(fig,'Video exported.','Export','Icon','success');
catch ME, uialert(fig,ME.message,'Export video'); end
end
function exportPreviewVideo(fig,path)
app=getApp(fig); pv=app.preview; big=app.detIntervals;
if isempty(pv.vr)||isempty(pv.times), error('No interval selected for preview.'); end
vw=VideoWriter(path,'Motion JPEG AVI'); vw.FrameRate=pv.vr.FrameRate; vw.Quality=90; open(vw);
dec=max(1,round(numel(pv.times)/1500));                     % cap ~1500 frames
for k=1:dec:numel(pv.times)
    t=pv.times(k);
    try, pv.vr.CurrentTime=min(max(t,0),max(0,pv.vr.Duration-1/pv.vr.FrameRate)); fr=readFrame(pv.vr); catch, continue; end
    if ndims(fr)<3, fr=repmat(fr,1,1,3); end
    [~,fi]=min(abs(double(big.time)-t));
    writeVideo(vw,drawWalls(im2uint8(fr),big.idxL(fi,:),big.idxR(fi,:)));
end
close(vw);
end
function I=drawWalls(I,xL,xR)
%DRAWWALLS  paint the left/right wall positions (per row) as green pixels
[H,W,~]=size(I); nY=min(H,numel(xL));
for y=1:nY
    for xx=[round(xL(y)) round(xR(y))]
        if isfinite(xx)&&xx>=1&&xx<=W
            cols=max(1,xx-1):min(W,xx+1);
            I(y,cols,1)=0; I(y,cols,2)=255; I(y,cols,3)=0;
        end
    end
end
end
function onCloseWorkbench(fig)
try, app=getApp(fig); app.cancel=true; app.running=false; setApp(fig,app); catch, end
previewStop(fig);
try, app=getApp(fig); if isfield(app,'preview'), app.preview.vr=[]; setApp(fig,app); end, catch, end
if isvalid(fig), delete(fig); end
end

%% ===================== EXPLORE logic - files & groups ============== %%
function exploreLoad(fig,paths)
%EXPLORELOAD  load one path (single-file, backwards compatible) or a cellstr of paths
exploreSetFiles(fig,paths,'');
end
function exploreSetFiles(fig,paths,groupPat)
%EXPLORESETFILES  register a set of result files, tag each with a group, focus the first
if nargin<3, groupPat=''; end
if ischar(paths)||isstring(paths), paths=cellstr(paths); end
app=getApp(fig);
ef=struct('path',{},'name',{},'group',{},'rec',{},'recnum',{},'resultPath',{});
for i=1:numel(paths)
    e=exploreEntry(char(paths{i}),groupPat);
    if isKey(app.exploreGroupOverride,e.path), e.group=app.exploreGroupOverride(e.path); end % keep manual groups
    ef(end+1)=e; %#ok<AGROW>
end
app.exploreFiles=ef; setApp(fig,app);
refreshExploreFileList(fig);
if ~isempty(ef), setExploreFocus(fig,1); end
refreshExploreIntervalItems(fig);
exploreRender(fig);
end
function e=exploreEntry(path,groupPat)
[~,nn,ex]=fileparts(path); name=[nn ex];
grp=exRegexpMatch(name,groupPat,'all');
e=struct('path',path,'name',name,'group',grp,'rec','1','recnum',1,'resultPath',path);
end
function m=exRegexpMatch(str,pat,dflt)
if isempty(pat), m=dflt; return; end
tok=regexp(str,pat,'match','once'); if isempty(tok), m=dflt; else, m=tok; end
end
function refreshExploreFileList(fig)
app=getApp(fig); ef=app.exploreFiles; L=app.c.explore.fileList;
items=cell(1,numel(ef));
for i=1:numel(ef), items{i}=sprintf('%s   [%s]',ef(i).name,ef(i).group); end
L.Items=items; L.ItemsData=1:numel(ef);
if isempty(ef), L.Value=[]; else, L.Value=1:numel(ef); end       % all selected by default
end
function idx=selectedIdx(fig)
app=getApp(fig); idx=app.c.explore.fileList.Value;
if isempty(idx), idx=1:numel(app.exploreFiles); end
idx=idx(:)';
end
function setExploreFocus(fig,i)
app=getApp(fig); ef=app.exploreFiles; if i<1||i>numel(ef), return; end
rec=getRec(fig,ef(i).path); app=getApp(fig);
app.explore=struct('path',ef(i).path,'intervals',rec.intervals,'s',rec.s,'meta',rec.meta);
app.exploreFocus=i; setApp(fig,app);
app.c.explore.srcLbl.Text=exploreSrcLabel(fig);
end
function txt=exploreSrcLabel(fig)
app=getApp(fig); ef=app.exploreFiles; n=numel(ef);
if n==0, txt='(nothing loaded)'; return; end
ng=numel(unique({ef.group})); foc=ef(app.exploreFocus);
txt=sprintf('%d file(s), %d group(s). Focus: %s - %d intervals, %s', ...
    n,ng,foc.name,numel(app.explore.intervals),pxLabel(app.explore.meta));
end
function onExploreFileSel(fig)
idx=selectedIdx(fig);
if ~isempty(idx), setExploreFocus(fig,idx(1)); end
refreshExploreIntervalItems(fig); exploreRender(fig);
end
function refreshExploreIntervalItems(fig)
%REFRESHEXPLOREINTERVALITEMS  interval filter = '(all)' + the union of level keys across the selection
app=getApp(fig); C=app.c.explore;
[~,levels,anyIndex]=selectedGroupSet(fig);
items=[{'(all)'}, levels(:)'];
keep=C.interval.Value; C.interval.Items=items;
if ismember(keep,items), C.interval.Value=keep; else, C.interval.Value='(all)'; end
if anyIndex, C.matchLbl.Text='Some file(s) have unnamed/duplicate intervals - matched by position (#k).';
else, C.matchLbl.Text=''; end
end
function uiExploreLoadFolder(fig)
d=uigetdir(pwd,'Select a folder of *_myograph.mat'); if isequal(d,0), return; end
app=getApp(fig); app.exploreRoot=d; setApp(fig,app); exploreScan(fig);
end
function exploreScan(fig)
app=getApp(fig);
if isempty(app.exploreRoot), alertUser(fig,'Load a folder first (Load folder...).','Scan'); return; end
inc=app.c.explore.includeF.Value; grp=app.c.explore.groupF.Value;
d=dir(fullfile(app.exploreRoot,'**','*.mat')); d=d(~[d.isdir]); paths={};
for i=1:numel(d)
    if isempty(regexp(d(i).name,inc,'once')), continue; end
    paths{end+1}=fullfile(d(i).folder,d(i).name); %#ok<AGROW>
end
if isempty(paths), alertUser(fig,'No files matched the Include pattern.','Scan'); return; end
exploreSetFiles(fig,paths,grp);
end
function exploreCreateGroup(fig)
app=getApp(fig); name=strtrim(app.c.explore.groupNameF.Value);
if isempty(name), alertUser(fig,'Type a group name first.','Create group'); return; end
idx=selectedIdx(fig);
if isempty(idx)||isempty(app.exploreFiles), alertUser(fig,'Select one or more files in the list first.','Create group'); return; end
for i=idx(:)'
    app.exploreFiles(i).group=name; app.exploreGroupOverride(app.exploreFiles(i).path)=name;
end
setApp(fig,app); refreshExploreFileList(fig); refreshExploreIntervalItems(fig); exploreRender(fig);
end
function apiExploreSelect(fig,idx)
app=getApp(fig); app.c.explore.fileList.Value=idx; setApp(fig,app); onExploreFileSel(fig);
end
function apiExploreCreateGroup(fig,name,idx)
app=getApp(fig); app.c.explore.fileList.Value=idx; app.c.explore.groupNameF.Value=name; setApp(fig,app);
exploreCreateGroup(fig);
end
function apiSetExplore(fig,field,value)
app=getApp(fig); C=app.c.explore;
switch field
    case 'plotType',   C.plotType.Value=value;
    case 'interval',   if ~ismember(value,C.interval.Items), C.interval.Items=[C.interval.Items {value}]; end, C.interval.Value=value;
    case 'organize',   C.organize.Value=value;
    case 'points',     C.points.Value=value;
    case {'vsmMarker','biomarker'}, C.vsmMarker.Value=value;   % 'biomarker' kept as an alias
    case 'propMetric', C.propMetric.Value=value;
    case 'stat',       C.stat.Value=value;
end
updateExploreEnable(fig); exploreRender(fig);
end
function rec=getRec(fig,path)
%GETREC  loaded {intervals,s,meta} for a result path, via a small LRU cache
app=getApp(fig); path=char(path);
if isKey(app.exploreCache,path)
    app.exploreOrder=[{path}, app.exploreOrder(~strcmp(app.exploreOrder,path))];
    rec=app.exploreCache(path); setApp(fig,app); return;
end
[iv,s,meta]=loadMyographResults(path);
rec=struct('intervals',iv,'s',s,'meta',meta);
app.exploreCache(path)=rec; app.exploreOrder=[{path}, app.exploreOrder];
while numel(app.exploreOrder)>app.exploreCacheLimit
    drop=app.exploreOrder{end}; app.exploreOrder(end)=[];
    if isKey(app.exploreCache,drop), remove(app.exploreCache,drop); end
end
setApp(fig,app);
end
function [files,levels,anyIndex]=selectedGroupSet(fig)
%SELECTEDGROUPSET  the selected files + the union of interval level keys across them
app=getApp(fig); ef=app.exploreFiles;
if isempty(ef), files=ef; levels={}; anyIndex=false; return; end
idx=selectedIdx(fig); files=ef(idx);
allKeys={}; anyIndex=false;
for fi=1:numel(files)
    rec=getRec(fig,files(fi).path); [keys,idxUsed]=levelKeys(rec.intervals);
    allKeys=[allKeys keys]; anyIndex=anyIndex||idxUsed; %#ok<AGROW>
end
levels=exOrderedCats(allKeys);
end
function [keys,idxUsed]=levelKeys(iv)
%LEVELKEYS  interval level key = name if non-empty & unique-in-file, else positional '#k'
n=numel(iv); nm=cell(1,n); keys=cell(1,n); idxUsed=false;
for k=1:n
    s=''; if isfield(iv,'name')&&(ischar(iv(k).name)||isstring(iv(k).name)), s=strtrim(char(iv(k).name)); end
    nm{k}=s;
end
for k=1:n
    s=nm{k};
    if isempty(s) || sum(strcmp(nm,s))>1, keys{k}=sprintf('#%d',k); idxUsed=true; else, keys{k}=s; end
end
end

%% ===================== EXPLORE logic - render dispatch ============= %%
function onExplorePlotType(fig)
updateExploreEnable(fig); exploreRender(fig);
end
function updateExploreEnable(fig)
%UPDATEEXPLOREENABLE  grey out the controls that do not apply to the current plot type
C=getApp(fig).c.explore; type=C.plotType.Value; isGroup=startsWith(type,'group:');
setEnable(C.vsmMarker,  strcmp(type,'group: vsmMarker (box)'));
setEnable(C.propMetric, strcmp(type,'group: propagation (box)'));
setEnable(C.organize,   isGroup);
setEnable(C.points,     strcmp(type,'group: vsmMarker (box)'));
setEnable(C.stat,       ismember(type,{'group: diameter','group: normalised diameter','group: spectrum'}));
end
function setEnable(h,on), if on, h.Enable='on'; else, h.Enable='off'; end, end
function exploreRender(fig)
app=getApp(fig);
if isempty(app.exploreFiles) && isempty(app.explore.intervals), return; end
ax=app.c.explore.ax; cla(ax,'reset');
try
    dispatchExplorePlot(fig,ax);
catch ME
    title(ax,'Could not plot'); app.c.explore.status.Text=['Plot error: ' ME.message];
end
end
function dispatchExplorePlot(fig,ax)
%DISPATCHEXPLOREPLOT  single point of truth for both the live axes and the export axes
app=getApp(fig); C=app.c.explore; type=C.plotType.Value;
if startsWith(type,'group:'), dispatchGroupPlot(fig,ax,type); return; end
e=app.explore;
if isempty(e.intervals)
    title(ax,'No focus file - load a file for single-file plots.'); C.status.Text='no focus file'; return;
end
iv=e.intervals; sel=C.interval.Value;
if strcmp(sel,'(all)'), ivIdx=1:numel(iv);
else, ivIdx=exIvIndex(iv,sel); if isempty(ivIdx), ivIdx=1:numel(iv); end
end
switch type
    case 'median diameter',              plotMedianDiam(ax,iv,ivIdx,false);
    case 'per-Y diameter',               plotPerYDiam(ax,iv,ivIdx,false);
    case 'per-Y normalised diameter',    plotPerYDiam(ax,iv,ivIdx,true);
    case 'normalised diameter',          plotMedianDiam(ax,iv,ivIdx,true);
    case 'Y-vs-time map',                plotYtMap(ax,iv,ivIdx(1));
    case 'spectrum vs f (avg)',          plotSpectrum(ax,iv,ivIdx);
    case 'percentile spectra',           plotPctSpectra(ax,iv,ivIdx);
    case 'amplitude percentiles',        plotAmpPct(fig,ax,iv,ivIdx);
    case 'wavelet spectrogram',          plotSpectrogram(ax,iv,ivIdx(1));
    case 'propagation: lag vs Y',        plotPropLag(ax,iv,ivIdx(1));
    case 'propagation: Y-vs-time+fit',   plotPropMap(ax,iv,ivIdx(1));
end
C.status.Text=type;
end
function k=exIvIndex(iv,label)
k=[];
if startsWith(label,'#')
    n=str2double(label(2:end)); if isfinite(n)&&n>=1&&n<=numel(iv), k=n; end, return;
end
m=find(strcmp({iv.name},label),1); if ~isempty(m), k=m; end
end
function dispatchGroupPlot(fig,ax,type)
app=getApp(fig); C=app.c.explore;
if isempty(app.exploreFiles), title(ax,'Load file(s) first for group plots.'); C.status.Text='no files'; return; end
organize=C.organize.Value; stat=C.stat.Value; nSelFiles=numel(selectedIdx(fig));
title(ax,groupTitle(fig,type),'Interpreter','none');       % set first so a renderer's "No data" title wins
switch type
    case 'group: diameter'
        [obs,meta]=gatherGroupCurves(fig,'diam');  renderGroupCurve(ax,obs,meta,organize,stat);
    case 'group: normalised diameter'
        [obs,meta]=gatherGroupCurves(fig,'ndiam'); renderGroupCurve(ax,obs,meta,organize,stat);
    case 'group: spectrum'
        [obs,meta]=gatherGroupCurves(fig,'spct');  renderGroupCurve(ax,obs,meta,organize,stat);
    case 'group: wavelet spectrogram'
        renderGroupSpectrogram(fig,ax);
    case 'group: vsmMarker (box)'
        [vals,tags]=gatherGroupScalars(fig,'vsmMarker'); renderGroupBox(ax,vals,tags,organize,C.vsmMarker.Value);
    case 'group: propagation (box)'
        [vals,tags]=gatherGroupScalars(fig,'prop');      renderGroupBox(ax,vals,tags,organize,['propagation ' C.propMetric.Value]);
end
C.status.Text=sprintf('%s | %d file(s), organise: %s',type,nSelFiles,organize);
end
function t=groupTitle(fig,type)
app=getApp(fig); sel=app.c.explore.interval.Value; what=strrep(type,'group: ','');
if strcmp(sel,'(all)'), t=sprintf('%s - all intervals',what); else, t=sprintf('%s - interval "%s"',what,sel); end
end

%% ===================== EXPLORE logic - group observations ========== %%
function [obs,meta]=gatherGroupCurves(fig,kind)
%GATHERGROUPCURVES  one curve per file x matched-interval, tagged group/interval/file
app=getApp(fig); C=app.c.explore; selInt=C.interval.Value;
[files,~,~]=selectedGroupSet(fig);
obs=struct('x',{},'y',{},'group',{},'interval',{},'file',{});
meta=struct('xlabel','','ylabel','','xlog',false);
for fi=1:numel(files)
    rec=getRec(fig,files(fi).path); iv=rec.intervals; keys=levelKeys(iv);
    for k=1:numel(iv)
        lev=keys{k};
        if ~strcmp(selInt,'(all)') && ~strcmp(lev,selInt), continue; end
        [x,y,ok,xl,yl,xlog]=oneCurve(iv(k),kind);
        if ~ok, continue; end
        meta.xlabel=xl; meta.ylabel=yl; meta.xlog=xlog;
        obs(end+1)=struct('x',x(:),'y',y(:),'group',files(fi).group,'interval',lev,'file',files(fi).name); %#ok<AGROW>
    end
end
end
function [x,y,ok,xl,yl,xlog]=oneCurve(ivk,kind)
%ONECURVE  a single file/interval curve for the group-curve plots (honours the off-FOV valid flag)
x=[]; y=[]; ok=false; xl=''; yl=''; xlog=false;
switch kind
    case {'diam','ndiam'}
        t=double(ivk.time); if isempty(t), return; end
        D=double(ivk.diameter);
        if isfield(ivk,'valid') && ~isempty(ivk.valid), D(~logical(ivk.valid(:)),:)=NaN; end  % drop off-FOV frames
        y=median(D,2,'omitnan'); x=t-t(1);                                  % align to interval start
        if strcmp(kind,'ndiam'), mm=mean(y,'omitnan'); if isfinite(mm)&&mm>0, y=y/mm; end, yl='normalised diameter'; else, yl='diameter (px)'; end % normalise to own interval mean
        xl='time from interval start (s)'; ok=any(isfinite(y));
    case 'spct'
        v=ivk.vasomotion; if isempty(v)||~isstruct(v)||~isfield(v,'fVectors')||~isfield(v.fVectors,'ampMean'), return; end
        x=double(v.f(:)); y=mean(double(v.fVectors.ampMean),1,'omitnan')'; xl='frequency (Hz)'; yl='amplitude (a.u.)'; xlog=true; ok=any(isfinite(y));
end
end
function [vals,tags]=gatherGroupScalars(fig,which)
%GATHERGROUPSCALARS  scalar per file x matched-interval (vsmMarker mean/per-Y, or a prop scalar)
app=getApp(fig); C=app.c.explore; selInt=C.interval.Value;
pointsFiles=strcmp(C.points.Value,'Files (mean over Y)');
[files,~,~]=selectedGroupSet(fig);
vals=[]; tags=struct('group',{},'interval',{},'file',{});
for fi=1:numel(files)
    rec=getRec(fig,files(fi).path); iv=rec.intervals; keys=levelKeys(iv);
    for k=1:numel(iv)
        lev=keys{k};
        if ~strcmp(selInt,'(all)') && ~strcmp(lev,selInt), continue; end
        if strcmp(which,'vsmMarker')
            d=vsmMarkerValue(iv(k).vasomotion,C.vsmMarker.Value); d=d(isfinite(d)); if isempty(d), continue; end
            if pointsFiles
                vals(end+1,1)=mean(d,'omitnan'); tags(end+1)=mkGT(files(fi),lev); %#ok<AGROW>
            else
                for j=1:numel(d), vals(end+1,1)=d(j); tags(end+1)=mkGT(files(fi),lev); end %#ok<AGROW>
            end
        else
            sc=propScalar(iv(k).prop,C.propMetric.Value);
            if isnan(sc), continue; end
            vals(end+1,1)=sc; tags(end+1)=mkGT(files(fi),lev); %#ok<AGROW>
        end
    end
end
end
function t=mkGT(file,lev), t=struct('group',file.group,'interval',lev,'file',file.name); end
function d=vsmMarkerValue(v,name)
%VSMMARKERVALUE  per-Y vector for a band-suffixed vasomotion scalar (scalars.<band>).
%   The tree holds bands as struct branches, so the flat dropdown suffixes each name
%   with the band: 'VB' (s.vFR) or 'CB' (s.cFR).  The remainder is the scalars.<band>
%   field (e.g. ampMeanVB -> v.scalars.VB.ampMean).  [] if absent.
d=[];
if isempty(v)||~isstruct(v)||~isfield(v,'scalars'), return; end
[band,fld]=splitBand(name);
if isempty(band)||~isfield(v.scalars,band)||~isfield(v.scalars.(band),fld), return; end
d=double(v.scalars.(band).(fld)); d=d(:);
end
function [band,fld]=splitBand(name)
%SPLITBAND  split a 'VB'/'CB'-suffixed marker name into its band + scalars field name.
band=''; fld='';
if numel(name)>2 && (endsWith(name,'VB')||endsWith(name,'CB'))
    band=name(end-1:end); fld=name(1:end-2);
end
end
function sc=propScalar(p,name)
sc=NaN; if isempty(p)||~isstruct(p), return; end
switch name
    case 'speed',      sc=gvp(p,'speed');
    case 'confidence', sc=gvp(p,'confidence');
    case 'R2',         sc=gvp(p,'R2');
    case 'domFreq',    sc=gvp(p,'domFreq');
    case 'phaseSpeed', if isfield(p,'phase')&&isstruct(p.phase), sc=gvp(p.phase,'speed'); end
end
end
function y=gvp(s,f), if isfield(s,f)&&~isempty(s.(f))&&isnumeric(s.(f)), y=double(s.(f)); y=y(1); else, y=NaN; end, end

%% ===================== EXPLORE logic - group renderers ============= %%
function renderGroupCurve(ax,obs,meta,organize,statName)
%RENDERGROUPCURVE  mean +/- band per series (series = the active Group/Interval combination)
if isempty(obs), title(ax,'No data for this selection'); return; end
key=exSeriesKeys(obs,organize); uk=unique(key,'stable'); cols=exColours(numel(uk));
[statFun,loFun,hiFun]=exStatFuns(statName); xg=exCommonGrid(obs);
hold(ax,'on'); hLeg=gobjects(1,numel(uk)); names=cell(1,numel(uk));
for m=1:numel(uk)
    sel=strcmp(key,uk{m}); Y=exResample(obs(sel),xg);
    mu=statFun(Y,2); lo=loFun(Y,2); hi=hiFun(Y,2);
    xx=[xg(:);flipud(xg(:))]; band=[lo(:);flipud(hi(:))]; good=isfinite(xx)&isfinite(band);
    if any(good), fill(ax,xx(good),band(good),cols(m,:),'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off'); end
    hLeg(m)=plot(ax,xg,mu,'-','Color',cols(m,:),'LineWidth',1.8); names{m}=uk{m};
end
grid(ax,'on'); xlabel(ax,meta.xlabel); ylabel(ax,meta.ylabel);
if meta.xlog, try, set(ax,'XScale','log'); catch, end, end
if numel(uk)>=2, legend(ax,hLeg,names,'Interpreter','none','Location','best'); else, legend(ax,'off'); end
hold(ax,'off');
end
function renderGroupBox(ax,vals,tags,organize,valueName)
%RENDERGROUPBOX  grouped box plots + jittered points (x/colour = the active Group/Interval dims)
if isempty(vals), title(ax,'No data for this selection'); legend(ax,'off'); return; end
nGrp=numel(unique({tags.group})); nInt=numel(unique({tags.interval}));
dims=exActiveDims(organize,nGrp,nInt); xName=dims{1}; colDims=dims(2:end);
xc=arrayfun(@(t)t.(xName),tags,'uni',0); xcats=exOrderedCats(xc);
hold(ax,'on');
if ~isempty(colDims)
    cc=arrayfun(@(t)strjoin(cellfun(@(d)t.(d),colDims,'uni',0),' | '),tags,'uni',0);
    ccats=exOrderedCats(cc); cols=exColours(numel(ccats)); ax.ColorOrder=cols;
    boxchart(ax,categorical(xc,xcats,'Ordinal',true),vals,'GroupByColor',categorical(cc,ccats),'MarkerStyle','none');
    exOverlayPoints(ax,xc,cc,vals,xcats,ccats,cols,true);
    legend(ax,cellstr(ccats),'Interpreter','none','Location','best');
else
    cols=exColours(numel(xcats));
    for k=1:numel(xcats)
        selk=strcmp(xc,xcats{k}); if ~any(selk), continue; end
        b=boxchart(ax,categorical(xc(selk),xcats,'Ordinal',true),vals(selk),'MarkerStyle','none');
        b.BoxFaceColor=cols(k,:); b.WhiskerLineColor=cols(k,:); b.BoxFaceAlpha=0.45;
    end
    exOverlayPoints(ax,xc,xc,vals,xcats,xcats,cols,false); legend(ax,'off');
end
grid(ax,'on'); ylabel(ax,valueName,'Interpreter','none'); xlabel(ax,exXLabel(xName)); hold(ax,'off');
end
function exOverlayPoints(ax,xc,cg,vals,xcats,cgcats,cols,dodge)
nC=numel(cgcats); bw=0.75/max(nC,1);
for j=1:nC
    sel=strcmp(cg,cgcats{j}); if ~any(sel), continue; end
    xi=double(categorical(xc(sel),xcats,'Ordinal',true));
    if dodge, off=(j-(nC+1)/2)*bw; else, off=0; end
    jit=(rand(numel(xi),1)-0.5)*bw*0.5;
    scatter(ax,xi(:)+off+jit,vals(sel),16,cols(j,:),'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none','HandleVisibility','off');
end
end
function renderGroupSpectrogram(fig,ax)
%RENDERGROUPSPECTROGRAM  average spectrum.amp (mean over Y) across the selected files on a common f/t grid
app=getApp(fig); C=app.c.explore; selInt=C.interval.Value;
[files,~,~]=selectedGroupSet(fig);
S={}; F={}; T={}; grps={};
for fi=1:numel(files)
    rec=getRec(fig,files(fi).path); iv=rec.intervals; keys=levelKeys(iv);
    for k=1:numel(iv)
        lev=keys{k}; if ~strcmp(selInt,'(all)') && ~strcmp(lev,selInt), continue; end
        v=iv(k).vasomotion;
        if isempty(v)||~isstruct(v)||~isfield(v,'spectrum')||~isfield(v.spectrum,'amp')||isempty(v.spectrum.amp)||~isfield(v,'timeDWT'), continue; end
        M=squeeze(mean(double(v.spectrum.amp),1,'omitnan'));  % [nf x nD] mean amplitude over Y (amp is already |WT|)
        f=double(v.f(:)); [f,o]=sort(f); M=M(o,:);        % ascending f
        t=double(v.timeDWT(:)); t=t-t(1);
        if size(M,1)~=numel(f)||size(M,2)~=numel(t)||numel(t)<2, continue; end
        S{end+1}=M; F{end+1}=f; T{end+1}=t; grps{end+1}=files(fi).group; %#ok<AGROW>
    end
end
if isempty(S), title(ax,'No stored wavelet spectrum in the selection (reprocess with ''spectrum'' in s.segVsmReturn)'); return; end
flo=max(cellfun(@min,F)); fhi=min(cellfun(@max,F));                 % common overlap grids
tlo=max(cellfun(@min,T)); thi=min(cellfun(@max,T));
nf=round(median(cellfun(@numel,F))); nt=round(median(cellfun(@numel,T)));
fg=logspace(log10(max(flo,eps)),log10(max(fhi,flo*1.01)),max(nf,2))'; tg=linspace(tlo,max(thi,tlo+eps),max(nt,2));
sumA=zeros(numel(fg),numel(tg)); cntA=zeros(numel(fg),numel(tg));   % NaN-aware average
for i=1:numel(S)
    Mi=interp2(T{i}(:)',F{i}(:),S{i},tg,fg,'linear',NaN);
    good=isfinite(Mi); sumA(good)=sumA(good)+Mi(good); cntA(good)=cntA(good)+1;
end
A=sumA./max(cntA,1); A(cntA==0)=NaN;
him=imagesc(ax,tg,fg,A); set(him,'AlphaData',isfinite(A)); set(ax,'YDir','normal');
try, set(ax,'YScale','log'); catch, end
colormap(ax,'parula'); cb=colorbar(ax); cb.Label.String='amplitude (a.u.)';
xlabel(ax,'time from interval start (s)'); ylabel(ax,'frequency (Hz)');
ng=numel(unique(grps));
if ng>1, C.status.Text=sprintf('averaged %d spectrogram(s) across %d groups - select one group''s files to isolate',numel(S),ng); end
end

%% ===================== EXPLORE export ============================== %%
function exploreExportFigure(fig,path,dpi,fmt), doExportFigure(fig,path,dpi,fmt); end
function exploreExportExcel(fig,path)
app=getApp(fig); ef=app.exploreFiles;
if isempty(ef)
    if isempty(app.explore.path), error('Nothing loaded to export.'); end
    exportMyographExcel(app.explore.path,path); return;
end
idx=selectedIdx(fig); sub=ef(idx);
exportMyographExcel({sub.path},path);
writeGroupSheet(sub,path);
end
function writeGroupSheet(ef,path)
%WRITEGROUPSHEET  append a small file->group mapping sheet next to the standard sheets
if isempty(ef), return; end
C=[{'file','group'}; [{ef.name}' {ef.group}']];
try, writecell(C,path,'Sheet','groups'); catch, end
end

%% ===================== UI glue (dialogs) =========================== %%
function uiAddFiles(fig)
[f,p]=uigetfile({'*.avi;*.mp4;*.mov;*.mj2;*.mkv','Videos'},'Select videos','MultiSelect','on');
if isequal(f,0), return; end; if ischar(f), f={f}; end
addFiles(fig,cellfun(@(x)fullfile(p,x),f,'uni',0));
end
function uiAddFolder(fig)
d=uigetdir(pwd,'Select a folder'); if isequal(d,0), return; end; addFolder(fig,d,true);
end
function uiRemove(fig)
app=getApp(fig); if isfield(app,'selFile')&&~isempty(app.selFile), app.queue(app.selFile)=[]; app.selFile=[]; setApp(fig,app); refreshQueue(fig); end
end
function uiClear(fig), app=getApp(fig); app.queue=emptyQueue(); setApp(fig,app); refreshQueue(fig); end
function ivAddRow(fig)
app=getApp(fig); D=app.c.setup.ivTbl.Data; D(end+1,:)={0,0,sprintf('interval%d',size(D,1)+1)}; app.c.setup.ivTbl.Data=D;
onIvTblEdit(fig);
end
function ivDelRow(fig)
app=getApp(fig); D=app.c.setup.ivTbl.Data; if ~isempty(D), D(end,:)=[]; end; app.c.setup.ivTbl.Data=D; onIvTblEdit(fig);
end
function onIvTblEdit(fig)
app=getApp(fig); if ~isfield(app,'selFile')||isempty(app.selFile), return; end
D=app.c.setup.ivTbl.Data; ivT=zeros(0,2); nm={};
for k=1:size(D,1)
    a=num(D{k,1},NaN); b=num(D{k,2},NaN); if isfinite(a)&&isfinite(b)&&b>a, ivT(end+1,:)=[a b]; nm{end+1}=char(string(D{k,3})); end %#ok<AGROW>
end
setFileIntervals(fig,app.selFile,ivT,nm);
end
function ivImport(fig)
[f,p]=uigetfile('*.mat','Import intervalsPerFile'); if isequal(f,0), return; end
L=load(fullfile(p,f)); app=getApp(fig); idx=app.selFile;
if isfield(L,'intervalsPerFile') && ~isempty(idx)
    ivT=L.intervalsPerFile{min(idx,numel(L.intervalsPerFile))}; nm={};
    if isfield(L,'intervalNamesPerFile'), nm=L.intervalNamesPerFile{min(idx,numel(L.intervalNamesPerFile))}; end
    setFileIntervals(fig,idx,ivT,nm);
end
end
function ivExport(fig)
[f,p]=uiputfile('*.mat','Export intervalsPerFile'); if isequal(f,0), return; end
app=getApp(fig); intervalsPerFile={app.queue.intervals}; intervalNamesPerFile={app.queue.intervalNames}; %#ok<NASGU>
save(fullfile(p,f),'intervalsPerFile','intervalNamesPerFile');
end
function saveSession(fig)
[f,p]=uiputfile('*.mat','Save session'); if isequal(f,0), return; end
app=getApp(fig); session.queue=app.queue; session.s=app.s; %#ok<STRNU>
save(fullfile(p,f),'-struct','session');
end
function loadSession(fig)
[f,p]=uigetfile('*.mat','Load session'); if isequal(f,0), return; end
S=load(fullfile(p,f)); app=getApp(fig);
if isfield(S,'queue'), app.queue=S.queue; end; if isfield(S,'s'), app.s=S.s; end
setApp(fig,app); refreshQueue(fig);
end
function uiExploreLoad(fig)
[f,p]=uigetfile('*_myograph.mat','Load result file(s)','MultiSelect','on'); if isequal(f,0), return; end
if ischar(f), f={f}; end
exploreSetFiles(fig,cellfun(@(x)fullfile(p,x),f,'uni',0),getApp(fig).c.explore.groupF.Value);
end
function uiExploreExportFig(fig)
[f,p]=uiputfile({'*.png';'*.tif';'*.pdf';'*.eps';'*.jpg'},'Export figure'); if isequal(f,0), return; end
app=getApp(fig); doExportFigure(fig,fullfile(p,f),str2double(app.c.explore.dpi.Value),app.c.explore.fmt.Value);
end
function uiExploreExportXlsx(fig)
[f,p]=uiputfile('*.xlsx','Export Excel'); if isequal(f,0), return; end
exploreExportExcel(fig,fullfile(p,f));
end
function doExportFigure(fig,path,dpi,fmt)
%DOEXPORTFIGURE  render the current Explore plot (single-file OR group) into an offscreen figure
if isempty(dpi)||isnan(dpi), dpi=300; end
tmp=figure('Visible','off','Color','w'); axe=axes(tmp); %#ok<LAXES>
dispatchExplorePlot(fig,axe);
switch upper(char(fmt))
    case {'PDF','EPS'}, exportgraphics(axe,path,'ContentType','vector','BackgroundColor','white');
    otherwise,          exportgraphics(axe,path,'Resolution',dpi,'BackgroundColor','white');
end
delete(tmp);
end

%% ===================== plotting primitives ========================= %%
function plotMedianDiam(ax,iv,ivIdx,norm)
hold(ax,'on'); for k=ivIdx(:)'
    t=double(iv(k).time); D=double(iv(k).diameter);
    if isfield(iv,'valid')&&~isempty(iv(k).valid), D(~logical(iv(k).valid(:)),:)=NaN; end   % ignore off-FOV frames
    m=median(D,2,'omitnan');
    if norm, mm=mean(m,'omitnan'); if isfinite(mm)&&mm>0, m=m/mm; end, end                  % normalise to own interval mean
    plot(ax,t,m,'DisplayName',iv(k).name);
end
grid(ax,'on'); xlabel(ax,'time (s)'); ylabel(ax,ternary(norm,'normalised diameter','diameter (px)')); legend(ax,'show');
end
function plotPerYDiam(ax,iv,ivIdx,norm)
k=ivIdx(1); t=double(iv(k).time); D=double(iv(k).diameter);
if isfield(iv,'valid')&&~isempty(iv(k).valid), D(~logical(iv(k).valid(:)),:)=NaN; end
step=max(1,round(size(D,2)/8)); hold(ax,'on');
for y=1:step:size(D,2)
    col=D(:,y);
    if norm, mm=mean(col,'omitnan'); if isfinite(mm)&&mm>0, col=col/mm; end, end            % each Y to its own interval mean
    plot(ax,t,col,'DisplayName',sprintf('Y=%d',y));
end
grid(ax,'on'); xlabel(ax,'time (s)'); ylabel(ax,ternary(norm,'normalised diameter','diameter (px)')); legend(ax,'show'); title(ax,iv(k).name);
end
function plotYtMap(ax,iv,k)
t=double(iv(k).time); D=double(iv(k).diameter); dec=max(1,round(size(D,1)/2000));
imagesc(ax,t(1:dec:end),1:size(D,2),D(1:dec:end,:)'); set(ax,'YDir','normal'); colorbar(ax);
xlabel(ax,'time (s)'); ylabel(ax,'Y'); title(ax,sprintf('%s : diameter (px)',iv(k).name));
end
function plotSpectrum(ax,iv,ivIdx)
hold(ax,'on'); for k=ivIdx(:)'
    v=iv(k).vasomotion; if isempty(v)||~isstruct(v)||~isfield(v,'fVectors')||~isfield(v.fVectors,'ampMean'), continue; end
    plot(ax,double(v.f),mean(double(v.fVectors.ampMean),1,'omitnan'),'DisplayName',iv(k).name);
end
set(ax,'XScale','log'); grid(ax,'on'); xlabel(ax,'frequency (Hz)'); ylabel(ax,'amplitude (a.u.)'); legend(ax,'show');
end
function plotSpectrogram(ax,iv,k)
%PLOTSPECTROGRAM  spectrum.amp (mean over Y) as frequency vs time
v=iv(k).vasomotion;
if isempty(v)||~isstruct(v)||~isfield(v,'spectrum')||~isfield(v.spectrum,'amp')||isempty(v.spectrum.amp)||~isfield(v,'timeDWT')
    title(ax,'no stored wavelet spectrum (reprocess with ''spectrum'' in s.segVsmReturn)'); return;
end
M=squeeze(mean(double(v.spectrum.amp),1,'omitnan'));  % [nf x nD] mean amplitude over Y (amp is already |WT|)
f=double(v.f(:)); [f,o]=sort(f); M=M(o,:);          % ascending f for the log axis
t=double(v.timeDWT(:)); t=t-t(1);
imagesc(ax,t,f,M); set(ax,'YDir','normal');
try, set(ax,'YScale','log'); catch, end
colormap(ax,'parula'); cb=colorbar(ax); cb.Label.String='amplitude (a.u.)';
xlabel(ax,'time from interval start (s)'); ylabel(ax,'frequency (Hz)');
title(ax,sprintf('%s : wavelet spectrogram (mean over Y)',iv(k).name),'Interpreter','none');
end
function plotPctSpectra(ax,iv,ivIdx)
% percentile-resolved spectra: amplitude vs frequency, one line per VB-envelope
% amplitude percentile bin (averaged over Y), legend giving the bin centre of each line.
k=ivIdx(1); v=iv(k).vasomotion;
if isempty(v)||~isstruct(v)||~isfield(v,'fVectors')||~isfield(v.fVectors,'VB') ...
        ||~isfield(v.fVectors.VB,'ampMeanPct')||isempty(v.fVectors.VB.ampMeanPct)
    title(ax,'no percentile spectra (fVectors.VB.ampMeanPct) in this file'); return;
end
f=double(v.f(:)); S=double(v.fVectors.VB.ampMeanPct); % [nY x nf x nBin] mean spectrum per VB amplitude bin
pc=double(v.pctCenters(:));                          % amplitude-percentile bin centres
M=squeeze(mean(S,1,'omitnan'));                      % [nf x nBin] mean amplitude over Y
if isvector(M), M=M(:); end
hold(ax,'on'); cols=lines(numel(pc));
for p=1:min(numel(pc),size(M,2))
    plot(ax,f,M(:,p),'Color',cols(p,:),'LineWidth',1.2,'DisplayName',sprintf('%gth pct',pc(p)));
end
set(ax,'XScale','log'); grid(ax,'on'); xlabel(ax,'frequency (Hz)'); ylabel(ax,'amplitude (a.u.)');
legend(ax,'show','Location','best'); title(ax,sprintf('%s : percentile-resolved spectra (avg over Y)',iv(k).name));
end
function plotAmpPct(fig,ax,iv,ivIdx)
% scalar amplitude percentiles: band-envelope amplitude vs percentile LEVEL, VB & CB
% (mean over Y).  scalars.VB/CB.ampPct is [nY x numel(s.pcts)] = prctile(a(t),s.pcts),
% so it is indexed by the percentile LEVELS s.pcts (0..100 %) - NOT the pctCenters bin
% centres that key the fVector percentile SPECTRA (plotPctSpectra above).  fig is passed
% (not recovered from ax) so app.s.pcts is reachable on the offscreen export axes too.
k=ivIdx(1); v=iv(k).vasomotion;
if isempty(v)||~isstruct(v)||~isfield(v,'scalars')||~isfield(v.scalars,'VB') ...
        ||~isfield(v.scalars.VB,'ampPct')||isempty(v.scalars.VB.ampPct)
    title(ax,'no amplitude percentiles (scalars.VB.ampPct) in this file'); return;
end
aVB=mean(double(v.scalars.VB.ampPct),1,'omitnan'); aVB=aVB(:);   % mean over Y -> [numel(s.pcts) x 1]
app=getApp(fig); pcts=(0:10:100)';
if isfield(app,'s')&&isfield(app.s,'pcts')&&~isempty(app.s.pcts), pcts=double(app.s.pcts(:)); end
if numel(pcts)~=numel(aVB), pcts=(1:numel(aVB))'; end            % axis/level count mismatch -> plain index
hold(ax,'on');
plot(ax,pcts,aVB,'-o','LineWidth',1.3,'DisplayName','VB');
if isfield(v.scalars,'CB')&&isfield(v.scalars.CB,'ampPct')&&~isempty(v.scalars.CB.ampPct)
    aCB=mean(double(v.scalars.CB.ampPct),1,'omitnan'); aCB=aCB(:);
    if numel(aCB)==numel(pcts), plot(ax,pcts,aCB,'-s','LineWidth',1.3,'DisplayName','CB'); end
end
grid(ax,'on'); xlabel(ax,'amplitude percentile (%)'); ylabel(ax,'amplitude (a.u.)');
legend(ax,'show','Location','best');
title(ax,sprintf('%s : amplitude percentiles (avg over Y)',iv(k).name),'Interpreter','none');
end
function plotPropLag(ax,iv,k)
p=iv(k).prop; if isempty(p)||~isfield(p,'lagByRow')||isempty(p.lagByRow), title(ax,'no propagation result'); return; end
yy=p.lagByRow(:,1); lg=p.lagByRow(:,2)*1000;
plot(ax,yy,lg,'.'); hold(ax,'on');
if isfield(p,'diag')&&isfield(p.diag,'fitLag')
    b=p.diag.fitLag; plot(ax,yy,(b(1)+b(2)*yy)*1000,'r-','LineWidth',1.5);
end
grid(ax,'on'); xlabel(ax,'Y (row)'); ylabel(ax,'lag at max correlation (ms)');
title(ax,sprintf('%s : %s %.0f %s | conf %.2f (%s), R2=%.2f p=%.3f', ...
    iv(k).name,dirOr(p),spd(p),unitOr(p),num0(p,'confidence'),confLevelOr(p),num0(p,'R2'),num0(p,'pValue')));
end
function plotPropMap(ax,iv,k)
p=iv(k).prop;
if isempty(p)||~isfield(p,'diag')||~isfield(p.diag,'Xf'), title(ax,'no propagation diagnostics'); return; end
Xf=p.diag.Xf; tt=p.diag.tt; cl=prctile(abs(Xf(:)),98);
imagesc(ax,tt,1:size(Xf,2),Xf'); set(ax,'YDir','normal'); if cl>0, clim(ax,[-cl cl]); end; colorbar(ax);
xlabel(ax,'time (s)'); ylabel(ax,'Y'); title(ax,sprintf('%s : fluctuation (tilt=propagation)',iv(k).name));
end

%% ===================== small helpers =============================== %%
% ---- Explore: sectioned control-column builders (ported from guiExplore.m) ----
function s=exSection(parent,titleText,nBodyRows)
p=uipanel(parent,'Title',titleText,'FontWeight','bold');
s=uigridlayout(p,[nBodyRows 2],'RowHeight',repmat({'fit'},1,nBodyRows),'ColumnWidth',{'fit','1x'}, ...
    'RowSpacing',5,'ColumnSpacing',6,'Padding',[6 6 6 6]);
end
function d=exDrop(g,row,name,items,cb)
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
d=uidropdown(g,'Items',items); d.Layout.Row=row; d.Layout.Column=2;
if nargin>=5 && ~isempty(cb), d.ValueChangedFcn=cb; end
end
function f=exField(g,row,name,val)
lb=uilabel(g,'Text',name); lb.Layout.Row=row; lb.Layout.Column=1;
f=uieditfield(g,'text','Value',val); f.Layout.Row=row; f.Layout.Column=2;
end
function b=exBtnFull(g,row,txt,cb)
b=uibutton(g,'Text',txt,'ButtonPushedFcn',cb); b.Layout.Row=row; b.Layout.Column=[1 2];
end

% ---- Explore: group-comparison math (ported from guiExplore.m) ----
function dims=exActiveDims(organize,nGrp,nInt)
%EXACTIVEDIMS  which tag dims separate the data; dims{1}=box x-axis, the rest=series colour
switch organize
    case 'Group',            prim={'group'};
    case 'Interval',         prim={'interval'};
    case 'Group x Interval', prim={'interval','group'};   % x=interval, colour=group
    case 'Pool all',         prim={};
    otherwise                                             % Auto: use whatever actually varies
        prim={}; if nInt>1, prim{end+1}='interval'; end, if nGrp>1, prim{end+1}='group'; end
end
dims=prim; if isempty(dims), dims={'group'}; end
end
function key=exSeriesKeys(obs,organize)
nGrp=numel(unique({obs.group})); nInt=numel(unique({obs.interval}));
dims=exActiveDims(organize,nGrp,nInt); key=cell(1,numel(obs));
for i=1:numel(obs), key{i}=strjoin(cellfun(@(d)obs(i).(d),dims,'uni',0),' | '); end
end
function cols=exColours(n)
%EXCOLOURS  distinct qualitative (colour-blind-aware) palette, cycling for large n
pal=[0.12 0.47 0.71;0.85 0.37 0.01;0.17 0.63 0.17;0.84 0.15 0.16;0.58 0.40 0.74; ...
     0.55 0.34 0.29;0.89 0.47 0.76;0.45 0.45 0.45;0.74 0.74 0.13;0.09 0.75 0.81];
cols=zeros(max(n,1),3);
for i=1:n
    if i<=size(pal,1), cols(i,:)=pal(i,:); else, cols(i,:)=hsv2rgb([mod((i-1)*0.61803,1) 0.65 0.85]); end
end
end
function [statFun,loFun,hiFun]=exStatFuns(name)
switch name
    case 'Mean +/- SEM'
        statFun=@(Y,d)mean(Y,d,'omitnan');
        sem=@(Y,d)std(Y,0,d,'omitnan')./sqrt(max(sum(isfinite(Y),d),1));
        loFun=@(Y,d)statFun(Y,d)-sem(Y,d); hiFun=@(Y,d)statFun(Y,d)+sem(Y,d);
    case 'Median +/- IQR'
        statFun=@(Y,d)median(Y,d,'omitnan');
        loFun=@(Y,d)exPrctileDim(Y,25,d); hiFun=@(Y,d)exPrctileDim(Y,75,d);
    otherwise
        statFun=@(Y,d)mean(Y,d,'omitnan'); sd=@(Y,d)std(Y,0,d,'omitnan');
        loFun=@(Y,d)statFun(Y,d)-sd(Y,d); hiFun=@(Y,d)statFun(Y,d)+sd(Y,d);
end
end
function p=exPrctileDim(Y,q,dim)
if dim==2, p=nan(size(Y,1),1); for i=1:size(Y,1), p(i)=prctile(Y(i,:),q); end
else,      p=nan(1,size(Y,2)); for i=1:size(Y,2), p(i)=prctile(Y(:,i),q); end
end
end
function xg=exCommonGrid(obs)
%EXCOMMONGRID  shared x-grid for pooling curves (identical grids used as-is)
X={obs.x};
if all(cellfun(@(x)isequal(x,X{1}),X)), xg=X{1}(:); return; end
lo=max(cellfun(@min,X)); hi=min(cellfun(@max,X));
if ~(hi>lo), lo=min(cellfun(@min,X)); hi=max(cellfun(@max,X)); end
nn=round(median(cellfun(@numel,X))); xg=linspace(lo,hi,max(nn,2))';
end
function Y=exResample(obs,xg)
Y=nan(numel(xg),numel(obs));
for i=1:numel(obs)
    x=obs(i).x(:); y=obs(i).y(:);
    if isequal(x,xg), Y(:,i)=y; else, Y(:,i)=interp1(x,y,xg,'linear',NaN); end
end
end
function cats=exOrderedCats(vals)
%EXORDEREDCATS  unique categories: numbered ones ascending (e.g. #1,#2), the rest first-seen
u=unique(vals,'stable');
num=nan(1,numel(u));
for i=1:numel(u), tk=regexp(u{i},'\d+','match','once'); if ~isempty(tk), num(i)=str2double(tk); end, end
if all(isnan(num)), cats=u(:)'; return; end
[~,ord]=sortrows([num(:) (1:numel(u))'],[1 2]); nanmask=isnan(num(ord)); ord=[ord(~nanmask); ord(nanmask)];
cats=u(ord(:)');
end
function xl=exXLabel(name)
switch name, case 'interval', xl='interval'; case 'group', xl='group'; otherwise, xl=name; end
end
function snapY(roi,yl), p=roi.Position; if p(2)~=yl(1)||p(4)~=diff(yl), roi.Position=[p(1) yl(1) p(3) diff(yl)]; end, end
function [err,note]=validateIntervals(ivT,names,iv)
%VALIDATEINTERVALS  err = blocking problems; note = non-blocking overlap/gap warnings
err=''; note='';
if isempty(ivT), err='No intervals defined.'; return; end
if numel(names)~=size(ivT,1), err='Number of names must equal number of intervals.'; return; end
if any(cellfun(@(x)isempty(strtrim(char(x))),names)), err='All intervals must be named.'; return; end
if any(ivT(:,2)<=ivT(:,1)), err='Every interval must have end > start (no zero-length).'; return; end
if ~isempty(iv)
    t=double(iv.time);
    if any(ivT(:,1)<t(1)-eps) || any(ivT(:,2)>t(end)+eps), err='Intervals must lie within the recording bounds.'; return; end
end
so=sortrows(ivT); notes={};
if any(so(2:end,1)<so(1:end-1,2)-eps), notes{end+1}='intervals overlap'; end
gaps=so(2:end,1)-so(1:end-1,2);
if any(gaps>1), notes{end+1}=sprintf('%d gap(s) between intervals (largest %.0f s)',nnz(gaps>1),max(gaps)); end
note=strjoin(notes,'; ');
end
function v=paramValue(s,f), if isfield(s,f), v=s.(f); else, v=[]; end, end
function s=val2str(v)
if ischar(v), s=v; elseif islogical(v), s=mat2str(v); elseif isscalar(v), s=num2str(v);
elseif isempty(v), s=''; else, s=mat2str(v); end
end
function v=str2paramValue(str)
str=strtrim(str);
if isempty(str), v=[]; return; end
n=str2num(str); %#ok<ST2NM>
if ~isempty(n), v=n; else, v=str; end
end
function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
function y=cropv(c,i), if numel(c)>=i, y=c(i); else, y=[]; end, end
function c=setcrop(c,i,v), if numel(c)<2, c=[NaN NaN]; end, c(i)=num(v,NaN); if all(isnan(c)), c=[]; end, end
function y=rr2(rr), if numel(rr)>=2, y=rr(2); else, y=Inf; end, end
function y=pxv(px), if isempty(px), y=[]; else, y=px; end, end
function y=fpsv(f), if isempty(f), y=[]; else, y=f; end, end
function n=num(v,dflt), if isnumeric(v), n=v; elseif ischar(v)||isstring(v), n=str2double(v); else, n=dflt; end, if isempty(n)||isnan(n), n=dflt; end, end
function n=numEmpty(v), if isnumeric(v)&&~isempty(v), n=v; elseif (ischar(v)||isstring(v))&&~isempty(strtrim(char(v))), n=str2double(v); else, n=[]; end, end
function s=nmk(nm,k), if iscell(nm)&&numel(nm)>=k, s=char(nm{k}); else, s=sprintf('interval%d',k); end, end
function r=rnd(tim,f), if isstruct(tim)&&isfield(tim,f), r=round(tim.(f),1); else, r=NaN; end, end
function s=pxLabel(m), if isfield(m,'pixelSize')&&~isempty(m.pixelSize)&&m.pixelSize>0, s=sprintf('%.3g µm/px',m.pixelSize); else, s='uncalibrated (px)'; end, end
function d=dirOr(p), if isfield(p,'direction'), d=p.direction; else, d='?'; end, end
function v=spd(p), if isfield(p,'speed'), v=p.speed; else, v=NaN; end, end
function u=unitOr(p), if isfield(p,'speedUnit'), u=p.speedUnit; else, u=''; end, end
function c=confLevelOr(p), if isfield(p,'confidenceLevel'), c=p.confidenceLevel; else, c=''; end, end
function y=num0(p,f), if isfield(p,f)&&~isempty(p.(f)), y=double(p.(f)); y=y(1); else, y=NaN; end, end
function t=getTip(tips,f), if isKey(tips,f), t=tips(f); else, t=f; end, end
function tips=paramTips()
tips=containers.Map();
tips('edgeMode')='Wall rule: center (darkness-weighted band centre; robust default), min (darkest point), outer/inner edge.';
tips('tau')='Wall threshold on the row-normalised image (<tau = wall).';
tips('rowRange')='[lo hi] rows to measure; others are interpolated.';
tips('smoothSigma')='Gaussian pre-smoothing, px.';
tips('dustRadius')='Vertical half-extent of dust blobs to remove (0=off).';
tips('tSpan')='Temporal outlier window, frames (slow-change prior).';
tips('ySpan')='Along-Y outlier window, rows (straight-vessel prior; set 0 for propagation-safe).';
tips('outlierK')='Outlier threshold factor (scaled MADs).';
tips('subpixel')='Parabolic sub-pixel refinement for min mode.';
tips('smoothSpan')='Final along-Y Savitzky-Golay span (0=off; 0 is propagation-safe).';
tips('vFR')='Vasomotion band [lo hi] Hz.';
tips('cFR')='Comparison band [lo hi] Hz.';
tips('wFR')='Wavelet frequency limits [lo hi] Hz.';
tips('wVPO')='Wavelet voices per octave.';
tips('normalisation')='mean/median (global) or mmean/mmedian (moving).';
tips('normsize')='Window for moving normalisation; inf/0 = global.';
tips('tgtFS')='Target sampling frequency for the kept spectrum, Hz.';
tips('segVsmReturn')='Analysis levels to compute & store (tick what you need).';
tips('segVsmReturn_bands')='(I) Band scalars: ampMean/ampStd/ampSkew per band + VB frequency-shape (fCent/fSprd/shapePeak/nPeak). Keep on.';
tips('segVsmReturn_moments')='(II) Per-frequency spectra (fVectors ampMean/ampStd/ampSkew + VB percentile spectra).';
tips('segVsmReturn_series')='(III) Amplitude / fCent / fSprd / nPeak time series (timeVectors.VB amp/fCent/fSprd/nPeak, timeVectors.CB amp).';
tips('segVsmReturn_clustering')='(IV) Otsu flare/silence markers + amplitude & flare time series.';
tips('segVsmReturn_reconstruction')='(V) Band-limited reconstruction rData (preferred propagation input).';
tips('segVsmReturn_spectrum')='(VI) Decimated amplitude/phase wavelet grid spectrum.amp/.phase (large; needed for the spectrogram).';
tips('propNShuffle')='Surrogate iterations for the propagation p-value.';
tips('propSignal')='Signal for the max-correlation lag: diameter (default, drift-removed) or rData (band-limited).';
tips('detectPerInterval')='OFF = detect once then slice (Path B, default). ON = re-detect per interval (Path A).';
tips('pixelSize')='µm per px; empty or 0 -> report in px.';
end
