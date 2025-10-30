%setTypes  Semi-automatic vessel/ROI labelling for *_BFI_d.mat datasets
%
%   setTypes(s,fNames) opens an interactive GUI for each BFI dataset in
%   fNames, displays the BFI image together with a provisional artery/vein
%   guess (derived from PI(BFI) heuristics), and lets the user paint or
%   relabel regions.  The routine then:
%       • writes the chosen type (“Artery”, “Vein”, “Uncertain”),
%         free-text label, and confidence score into RESULTS.sMetrics
%         (and dvsMetrics, if present)
%       • stores the signed confidence map in RESULTS.mapType (NaN outside
%         mask, −6 … +6 in vessels)
%       • updates *_BFI_r.mat and *_BFI_s.mat in place
%
%   When s.useReference is true the first file in fNames is treated as the
%   reference: subsequent files inherit its type/label definitions via
%   REGID matching, skipping the GUI.
%
%   INPUTS
%     s        parameter struct  
%                • useReference   true / false  
%                • refFName       path to reference file (if used)
%                • prchNSize      grid spacing for parenchyma fill-in
%     fNames   cell array of *_BFI_d.mat paths.
%
%   SIDE-EFFECTS
%       *_BFI_r.mat   RESULTS  – fields mapType, .type, .label, etc.
%       *_BFI_s.mat   SETTINGS – sub-field settings.setTypes added
%
%   EXAMPLE
%     p.useReference = false;
%     D = dir(fullfile(dataRoot,'*_BFI_d.mat'));
%     setTypes(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     Image Processing Toolbox (bwskel, visboundaries, etc.) and core LSCI
%     utility functions.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: Header generated with ChatGPT; minor inconsistencies may remain—
%   please verify before release.
%   ----------------------------------------------------------------------

%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
% %If no reference is used
% % s.useReference=false;
% % s.refFName='';
% %If reference is used without proper file naming
% % s.useReference=ture; %ASSUMES PRE-REGISTERED FILES
% % s.refFName=REFERENCE FILE NAME;
% %If reference is used in an automated setting
% s.useReference=true; %Assumes PRE-registered files
% %Settings tabs with a reference has to be done ROI by ROI if split ROIS are used. It also has
% %to be done ANIMAL by ANIMAL if multiple animals are compared. It can be
% %convienintly done in a loop as below, but REQUIRES proper file naming.
% for idxA=5 %animals index list
%     for idxR=1:2 %ROIs index list
%         files      = dir(fullfile(rootFolder,'**',sprintf('Roi%d*LH%03d*c_BFI_d.mat', idxR, idxA))); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
%         fNames     = fullfile({files.folder}', {files.name}');
%         s.refFName=fNames{1};
%         setTypes(s,fNames)
%     end
% end



function setTypes(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_BFI_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_BFI_d.mat".');
end

for fidx=1:1:numel(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');
    if (fidx==1 &&  s.useReference ) || ~s.useReference
        idxs=results.sMetrics{:,'idx'};
        ctgs=results.sMetrics{:,'category'};
        bfi=results.sMetrics{:,'BFI'};
        guess=zeros(size(idxs));
        cMask=results.cMask;
        if any(strcmp('PI(BFI)',results.sMetrics.Properties.VariableNames))
            var=results.sMetrics{:,"PI(BFI)"};
            thrsh=median(var(bfi<prctile(bfi(ctgs==5),10) & ctgs==5),'omitnan');
            for i=1:1:numel(idxs)
                if ctgs(i)==5
                    guess(idxs(i):idxs(i)+2)=guess(idxs(i):idxs(i)+2)+...
                        (~isnan(idxs(i)+1) && var(idxs(i))>var(idxs(i)+1))+...
                        (~isnan(idxs(i)+1) && ~isnan(idxs(i)+2) && var(idxs(i)+1)>var(idxs(i)+2))+...
                        (~isnan(idxs(i)+2) && var(idxs(i))>var(idxs(i)+2))+...
                        (var(idxs(i))>thrsh)+...
                        (var(idxs(i))>median(var(ctgs==5 & var>thrsh),'omitnan'));
                    guess(idxs(i):idxs(i)+2)=guess(idxs(i):idxs(i)+2)-...
                        (~isnan(idxs(i)+1) && var(idxs(i))<var(idxs(i)+1))-...
                        (~isnan(idxs(i)+1) && ~isnan(idxs(i)+2) && var(idxs(i)+1)<var(idxs(i)+2))-...
                        (~isnan(idxs(i)+2) && var(idxs(i))<var(idxs(i)+2))-...
                        (var(idxs(i))<thrsh)-...
                        (var(idxs(i))<median(var(ctgs==5 & var<thrsh),'omitnan'));
                end
            end
        end
        map=[0;guess];
        map=map(results.sMap+1);
        map(cMask<4)=NaN;
        cmap=zeros(11,3);
        cmap(1:5,3)=(1:5)./10+0.5;
        cmap(6,2)=1;
        cmap(7:11,1)=(5:-1:1)./10+0.5;
        img=sqrt(results.imgBFI);
        img=mat2gray(img,double(prctile(img(cMask(:)>0),[5,99])));
        fSize=floor((min(size(img))./20))*2+1;
        img=(img-imopen(medfilt2(img,[fSize,fSize],"symmetric"),strel('disk',fSize))).*(cMask>0);
        rois=repmat("",1,numel(guess));


        f  = figure('Name','Guided labeling','Color','w','WindowState','maximized', ...
            'CloseRequestFcn',@(~,~)uiresume(gcbf));
        TL = tiledlayout(f,1,5,'TileSpacing','none','Padding','compact');
        axL = nexttile(TL,[1 2]);
        imagesc(axL,img), axis(axL,'image','off'), colormap(axL,'parula');
        title('BFI');
        axR = nexttile(TL,[1 2]);
        hMap = imagesc(axR,map,'AlphaData',~isnan(map));
        axis(axR,'image','off'), clim(axR,[-6 6]), colormap(axR,cmap);
        title('Labels');
        axCtl = nexttile(TL); axis(axCtl,'off');
        axis tight
        panelWH = [150 300];
        p = uipanel(f,'Units','pixels','BorderType','none', 'BackgroundColor','w');
        reposition(f,p,axCtl,panelWH)
        f.SizeChangedFcn = @(~,~)reposition(f,p,axCtl,panelWH);
        uicontrol(p,'Style','text','String','Left panel image', ...
            'Units','normalized','Position',[.05 .90 .90 .08], ...
            'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
        pop = uicontrol(p,'Style','popup','Units','normalized', ...
            'Position',[.05 .81 .90 .08]);
        uicontrol(p,'Style','text','String','Labelling controls', ...
            'Units','normalized','Position',[.05 .7 .90 .06], ...
            'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
        btnY = [.60 .50 .40 .50 .40];
        btnX = [.05 .05 .05 .55 .55];
        lbl   = {'Artery','Vein','Uncertain', 'ROI Vessel','ROI Parench'};
        val = [6 -6 0 7 8];
        btn   = gobjects(1,5);
        for k = 1:5
            btn(k) = uicontrol(p,'String',lbl{k},'Units','normalized', ...
                'Position',[btnX(k) btnY(k) .45 .08], ...
                'Callback',@(src,~)pickCat(src,val(k)));
        end
        txt=uicontrol(p,'Style','edit','String','Label?', ...
            'Units','normalized','Position',[.55 .6 .45 .08], ...
            'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
        finishBtn = uicontrol(p,'String','Finish','Units','normalized', ...
            'Position',[.05 .00 .90 .08], ...
            'Callback',@(~,~)uiresume(gcbf));
        setappdata(f,'btnHandles',btn);
        setappdata(f,'txtHandle',txt);

        % ---------- build selector list (2-D fields == size(map)) ---------------
        names={}; data={};  chk=@(v) isnumeric(v)&&ismatrix(v)&&isequal(size(v),size(map));
        for fn = fieldnames(results).',                  v=results.(fn{1});
            if chk(v), names{end+1}=fn{1}; data{end+1}=v; end, end
        if isfield(results,'extendedMetrics')
            for fn = fieldnames(results.extendedMetrics).', v=results.extendedMetrics.(fn{1});
                if chk(v), names{end+1}=['extendedMetrics.' fn{1}]; data{end+1}=v; end, end
        end
        idx0 = find(strcmp(names,'imgBFI'),1);   % try to start at "imgBFI"
        if isempty(idx0), idx0 = 1; end
        pop.UserData = struct('ax',axL,'data',{data},'names',{names});
        set(pop,'String',names,'Value',idx0,'Callback',@updateImg);   % ← set start

        % ---------- share data for painting routine -----------------------------
        guidata(f,struct('map',map,'cMask',cMask,'guess',guess,'rois',rois,'sMap',results.sMap,'hMap',hMap,'axL',axL,'axR',axR));
        setappdata(f,'paintVal',[]);
        set(f,'WindowButtonDownFcn',@startPaint,'WindowButtonUpFcn',@stopPaint);

        uiwait(f);                               % pause until Finish pressed
        guess = guidata(f).guess;
        rois = guidata(f).rois;
        delete(f);    % retrieve & close

        map=[0;guess];
        map=map(results.sMap+1);
        map(cMask<3)=NaN;
        results.mapType=map;


        type = repmat("Uncertain",size(guess));  % default = 0
        type(guess < 0) = "Vein";                % negative values
        type(guess > 0) = "Artery";
        idxs=results.sMetrics.idx;

        results.sMetrics.type=type;
        results.sMetrics.label=rois';
        results.sMetrics.typeConfidence=guess;
        if isfield(results,"dvsMetrics")
        results.dvsMetrics.label=rois(ismember(results.sMetrics.idx, results.dvsMetrics.idx))';
        results.dvsMetrics.type=type(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
        results.dvsMetrics.typeConfidence=guess(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
        end
    else
        results.sMetrics.type = strings(height(results.sMetrics),1);
        results.sMetrics.label = strings(height(results.sMetrics),1);
        results.sMetrics.typeConfidence = zeros(height(results.sMetrics),1);
        if isfield(results,"dvsMetrics")
        results.dvsMetrics.type = strings(height(results.dvsMetrics),1);
        results.dvsMetrics.label = strings(height(results.dvsMetrics),1);
        results.dvsMetrics.typeConfidence = zeros(height(results.dvsMetrics),1);
        end

        [isHit, loc]   = ismember(results.sMetrics.RegID, idxs);
        results.sMetrics.type(isHit) = type(loc(isHit));
        results.sMetrics.typeConfidence(isHit) = guess(loc(isHit));
        results.sMetrics.label(isHit) = rois(loc(isHit))';

        if isfield(results,"dvsMetrics")
        [isHit, loc]   = ismember(results.dvsMetrics.RegID, idxs);
        results.dvsMetrics.type(isHit) = type(loc(isHit));
        results.dvsMetrics.typeConfidence(isHit) = guess(loc(isHit));
        results.dvsMetrics.label(isHit) = rois(loc(isHit))';
        end
    end

    settings.setTypes=s;
    disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
    save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');
end


%% =========================  LOCAL FUNCTIONS  ========================= %%
function startPaint(fig,~)
set(fig,'WindowButtonMotionFcn',@doPaint);
doPaint(fig);
end
function stopPaint(fig,~)
set(fig,'WindowButtonMotionFcn',[]);
end
function doPaint(fig,~)
val = getappdata(fig,'paintVal');
if isempty(val); return; end

gui = guidata(fig);
if ~isstruct(gui); return; end

% Which axes is the mouse currently over?
axUnder = ancestor(hittest(fig),'axes');
if isempty(axUnder) || ~(axUnder==gui.axL || axUnder==gui.axR)
    return                                   % clicked outside the images
end

% Pixel → matrix indices
cp = get(axUnder,'CurrentPoint');
x  = round(cp(1,1));
y  = round(cp(1,2));
if x<1 || y<1 || x>size(gui.map,2) || y>size(gui.map,1)
    return
end

lbl = gui.sMap(y,x);
if lbl==0, return, end

txt = getappdata(fig,'txtHandle');
txt = get(txt,'String');
if val==7 && lbl>0 && gui.cMask(y,x)==5
    gui.rois(lbl:lbl+2)=txt;
    mask = ismember(gui.sMap,lbl);
    hold(gui.axL,'on')
    visboundaries(gui.axL,mask,'Color','m','LineWidth',0.8);
    hold(gui.axL,'off')
    guidata(fig,gui);
    drawnow limitrate
elseif val==8 && lbl>0 && gui.cMask(y,x)==1
    gui.rois(lbl)=txt;
    mask = (gui.sMap==lbl);
    hold(gui.axL,'on')
    visboundaries(gui.axL,mask,'Color','g','LineWidth',0.8);
    hold(gui.axL,'off')
    guidata(fig,gui);
    drawnow limitrate
elseif val<7 && lbl>0 && gui.cMask(y,x)==5
    switch val
        case  6
            gui.guess(lbl:lbl+2) =  6;
        case -6
            gui.guess(lbl:lbl+2) = -6;
        otherwise
            gui.guess(lbl:lbl+2) =  0;
    end
    gui.map(gui.sMap==lbl) = val;
    gui.map(gui.sMap==lbl+1) = val;
    gui.map(gui.sMap==lbl+2) = val;
    set(gui.hMap,'CData',gui.map);
    guidata(fig,gui);
    drawnow limitrate
else
    return;
end
end

function updateImg(src,~)
u = src.UserData;  sel = src.Value;
mask = u.data{find(strcmp(u.names,'cMask'),1)};   % try to start at "imgBFI"
img=u.data{sel};
imagesc(u.ax,img), axis(u.ax,'image','off')
clim(prctile(img(mask>0),[1,99]));
title(u.ax,u.names{sel},'Interpreter','none')

fig = ancestor(src,'figure');
gui = guidata(fig);
for i = 1:numel(gui.rois)
    if strlength(gui.rois(i))==0, continue, end   % only labelled ROIs

    hold(gui.axL,'on')
    if any(gui.cMask(gui.sMap==i)==5)
        mask  = ismember(gui.sMap,i);
        visboundaries(gui.axL,mask,'Color','m','LineWidth',0.8);
    elseif  any(gui.cMask(gui.sMap==i)==1)
        mask  = (gui.sMap==i);
        visboundaries(gui.axL,mask,'Color','g','LineWidth',0.8);
    end
    hold(gui.axL,'off')
end
end
function pickCat(src,val)
fig = ancestor(src,'figure');
setappdata(fig,'paintVal',val)
bh  = getappdata(fig,'btnHandles');
set(bh,'Enable','on'), set(src,'Enable','off')
end
function reposition(figH,panelH,tileAx,wh)
figPix  = getpixelposition(figH);
tilePos = get(tileAx,'Position');
tilePix = [tilePos(1:2).*figPix(3:4)  tilePos(3:4).*figPix(3:4)];
x = tilePix(1) + max(0,(tilePix(3)-wh(1))/2);
y = tilePix(2) + max(0,(tilePix(4)-wh(2))/2);
set(panelH,'Position',[x y wh]);
end

end