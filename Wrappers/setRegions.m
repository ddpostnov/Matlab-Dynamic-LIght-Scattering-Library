%setRegions  Interactive, per-file selection of segmentation regions (regionsMask)
%
%   setRegions(s,fNames) opens an ROI editor for every file in fNames and writes the
%   labelled region mask (results.regionsMask) that runSegmentation / splitRegions
%   later use to restrict processing to user-chosen parts of the field of view.  It is
%   FULLY INTERACTIVE - there is NO count / headless parameter: the number of regions
%   is simply however many the user draws.  fNames is TWO-DIMENSIONAL (rows = groups,
%   e.g. one animal / FOV per row, exactly as getFileNamesList returns when grouped)
%   and setRegions iterates every file itself - there is no launcher for-loop.
%
%   CARRY-FORWARD WITHIN A GROUP, RESET ACROSS GROUPS.  For each row, files are visited
%   left to right.  On the first file the user draws the region(s); on every subsequent
%   file the SAME ROIs are re-offered as editable objects (nudge / add / delete / reset
%   / draw more).  At the next row the ROIs reset to empty.
%
%   THE EDITOR (one window per file).  ROIs are drawn on the display-enhanced image
%   (enhanceForDisplay), in the native image orientation so createMask aligns with the
%   stored frame.  Nothing advances until Done is pressed, so regions can be freely
%   adjusted first.  Controls:
%     * a shape selector - polygon / rectangle / square / ellipse / circle
%       (drawpolygon / drawrectangle (+FixedAspectRatio for square) / drawellipse /
%       drawcircle);
%     * Add ROI     - draw one more region of the selected shape;
%     * Delete ROI  - remove the selected ROI (or, if none is selected, the most recent
%                     one).  The Delete key also removes the currently selected ROI(s);
%     * Reset ROIs  - clear every ROI on this file;
%     * Done        - accept the current ROIs, save, and advance to the next file.
%
%   WHOLE WINDOW == NO MASK.  If the user draws NO ROIs (or skips setRegions entirely),
%   NO results.regionsMask is written - any stale one is removed.  Downstream readers
%   (runSegmentation, splitRegions) treat a MISSING regionsMask as the whole window (an
%   all-true mask the size of the image), so "draw nothing" == "use everything".  With
%   k ROIs drawn, regionsMask = sum_k createMask(roi_k).*k (k = the ROI index), identical
%   to the old runCategories math, so downstream (splitRegions, getPixelCategories) is
%   unchanged.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.setRegions untouched
%              (setRegions itself reads no numeric fields; the region count is drawn,
%              not configured).
%     fNames   2-D cell array of *_K_d.mat / *_I_d.mat paths.  Rows = groups; each file
%              must have matching *_s.mat and *_r.mat siblings.  Empty cells are skipped
%              (ragged rows from getFileNamesList are fine).
%
%   SIDE-EFFECTS (per file)
%     <name>_r.mat   results.regionsMask added/overwritten when >=1 ROI is drawn, or
%                    REMOVED when none is drawn (double, 0 = excluded, 1..N = labels)
%     <name>_s.mat   settings.setRegions = s
%
%   EXAMPLE
%     fNames = getFileNamesList(root,'*_K_d.mat','[A-Z]+\d+');  % grouped, rows=animals
%     setRegions(s,fNames);                                    % draw regions per file
%     % (skip this call entirely, or draw nothing, to segment the whole window)
%
%   DEPENDS ON
%     enhanceForDisplay (display background), and MATLAB's Image Processing Toolbox
%     ROI tools (drawpolygon/drawrectangle/drawellipse/drawcircle, createMask).
%
% See also: runSegmentation, splitRegions, getPixelCategories, enhanceForDisplay,
%           getFileNamesList
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026

function setRegions(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_K_d.mat')|| contains(x,'_I_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat" or "_I_d.mat".');
end

nGroups=size(fNames,1);
for g=1:1:nGroups
    carried=emptyROISpec();          % ROI geometry carried within this group (reset per group)
    for c=1:1:size(fNames,2)
        fName=fNames{g,c};
        if isempty(fName), continue; end
        disp(['setRegions: group ',num2str(g),'/',num2str(nGroups),', file ',fName])
        s.fName=fName;
        clearvars results settings source
        load(strrep(fName,'_d.mat','_s.mat'),'settings');
        load(strrep(fName,'_d.mat','_r.mat'),'results');

        % --- mean image + modality (only needed for the working frame size / display) ---
        isK=contains(fName,'_K_d.mat');
        if isK
            if isfield(results,'imgK')
                imgIni=results.imgK;
            else
                load(fName,'source');
                imgIni=mean(source.data,3,'omitmissing');
                clearvars source
            end
        else
            imgIni=results.imgI;
        end

        % --- interactive ROI editor (always opens); empty mask => whole window ---
        [regionsMask,carried]=editRegions(imgIni,isK,carried);

        if isempty(regionsMask)
            % no regions drawn: whole window - remove any stale mask, write none, so
            % downstream falls back to an all-true mask (see runSegmentation/splitRegions).
            if isfield(results,'regionsMask'), results=rmfield(results,'regionsMask'); end
        else
            results.regionsMask=regionsMask;
        end
        settings.setRegions=s;
        save(strrep(fName,'_d.mat','_s.mat'),'settings','-v7.3');
        save(strrep(fName,'_d.mat','_r.mat'),'results','-v7.3');
    end
end
end

% =====================================================================
function [regionsMask,carried]=editRegions(imgIni,isK,carried)
%editRegions  ROI editor for one file: draw/edit labelled regions on the enhanced
%   image, return the labelled mask (EMPTY when no ROI is drawn = whole window) plus
%   the ROI geometry to carry to the next file.

% Display-enhanced background the user draws on (mirrors the old runCategories prep).
if isK
    img=imgIni;
    img(img(:)>prctile(img(:),99))=prctile(img(:),99);
    img(img(:)<prctile(img(:),1))=prctile(img(:),1);
    img=imcomplement(img);
else
    img=imgIni;
end
fSize=floor((min(size(img))./20))*2+1;
img(isnan(img))=0;
imgDraw=enhanceForDisplay(img,fSize,min(15,fSize));

rois={};   % live images.roi handles (shared with the nested callbacks)

fig=figure('Name','setRegions - define segmentation regions','NumberTitle','off', ...
    'Color','w','Units','normalized','Position',[0.05 0.08 0.9 0.84]);

axOrig=axes(fig,'Units','normalized','Position',[0.04 0.20 0.44 0.74]);
imagesc(axOrig,imgIni); axis(axOrig,'image'); colormap(axOrig,'gray')
clim(axOrig,prctile(imgIni(:),[1,99])); title(axOrig,'Original image')

axDraw=axes(fig,'Units','normalized','Position',[0.52 0.20 0.44 0.74]);
imagesc(axDraw,imgDraw); axis(axDraw,'image'); colormap(axDraw,'gray')
clim(axDraw,prctile(imgDraw(:),[10,99]))
title(axDraw,'Add / edit regions, then press Done')

uicontrol(fig,'Style','text','Units','normalized','Position',[0.04 0.105 0.09 0.03], ...
    'String','ROI shape:','HorizontalAlignment','left','BackgroundColor','w');
shapePop=uicontrol(fig,'Style','popupmenu','Units','normalized', ...
    'Position',[0.13 0.11 0.14 0.035], ...
    'String',{'polygon','rectangle','square','ellipse','circle'});
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.29 0.105 0.10 0.045], ...
    'String','Add ROI','Callback',@onAdd);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.40 0.105 0.10 0.045], ...
    'String','Delete ROI','Callback',@onDelete);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.51 0.105 0.10 0.045], ...
    'String','Reset ROIs','Callback',@onReset);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.76 0.10 0.13 0.055], ...
    'String','Done (next file)','FontWeight','bold','Callback',@onDone);

hint=uicontrol(fig,'Style','text','Units','normalized','Position',[0.04 0.02 0.68 0.055], ...
    'String','','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);

% Pre-populate the carried ROIs (editable) - empty on the first file of a group.
for i=1:1:numel(carried)
    rois{end+1}=recreateROI(axDraw,carried(i)); %#ok<AGROW>
end

fig.UserData='editing';                       % private wait sentinel (see waitfor below)
fig.WindowKeyPressFcn=@onKey;                 % Delete key removes the selected ROI(s)
fig.CloseRequestFcn=@(~,~) set(fig,'UserData','done');   % closing the window == Done
updateHint();
% Block until Done / close.  We deliberately do NOT use uiwait/uiresume here: the
% interactive drawpolygon / drawrectangle / drawellipse / drawcircle calls run their
% OWN uiwait/uiresume on THIS figure while a shape is being drawn, and that internal
% uiresume also releases an outer uiwait(fig) - which made setRegions jump to the next
% file the instant a shape was finished.  waitfor on a private sentinel property is
% immune to that (the draw* tools never touch UserData), so we advance only when Done
% or window-close sets UserData to 'done'.
waitfor(fig,'UserData','done');

% --- after Done: assemble the labelled region mask + carry-forward geometry ---
% An EMPTY regionsMask signals "whole window" to the caller (no field is written).
regionsMask=[];
carried=emptyROISpec();
if isvalid(fig)
    k=0;
    for i=1:1:numel(rois)
        if isvalid(rois{i})
            if k==0, regionsMask=zeros(size(imgIni)); end
            k=k+1;
            regionsMask=regionsMask+createMask(rois{i}).*k;
            carried(k)=captureROI(rois{i});
        end
    end
    delete(fig);
end

    % ---- nested callbacks (share rois / axDraw / shapePop / hint) ----
    function onAdd(~,~)
        shapes={'polygon','rectangle','square','ellipse','circle'};
        r=drawNewROI(axDraw,shapes{shapePop.Value});
        if ~isempty(r) && isvalid(r)
            rois{end+1}=r;
        end
        updateHint();
    end
    function onDelete(~,~)
        idx=[];
        for jj=numel(rois):-1:1                    % prefer a user-selected ROI
            if isvalid(rois{jj}) && rois{jj}.Selected, idx=jj; break; end
        end
        if isempty(idx)
            for jj=numel(rois):-1:1                % else the most recently added
                if isvalid(rois{jj}), idx=jj; break; end
            end
        end
        if ~isempty(idx)
            if isvalid(rois{idx}), delete(rois{idx}); end
            rois(idx)=[];
        end
        updateHint();
    end
    function onReset(~,~)
        for jj=1:1:numel(rois)
            if isvalid(rois{jj}), delete(rois{jj}); end
        end
        rois={};
        updateHint();
    end
    function onKey(~,evt)
        if strcmp(evt.Key,'delete')                % Delete key = remove selected ROI(s)
            deleted=false;
            for jj=numel(rois):-1:1
                if isvalid(rois{jj}) && rois{jj}.Selected
                    delete(rois{jj}); rois(jj)=[]; deleted=true;
                end
            end
            if deleted, updateHint(); end
        end
    end
    function onDone(~,~)
        if isvalid(fig), fig.UserData='done'; end
    end
    function updateHint()
        n=0;
        for jj=1:1:numel(rois), if isvalid(rois{jj}), n=n+1; end, end
        if n==0
            msg=['No ROIs drawn: the WHOLE WINDOW will be used (no region mask is ', ...
                 'written).  Add ROIs to restrict segmentation, then press Done.'];
        else
            msg=[num2str(n) ' ROI(s).  Add / drag / reshape freely; select an ROI and ', ...
                 'press Delete (or "Delete ROI") to remove it.  Press Done when finished.'];
        end
        if isvalid(hint), set(hint,'String',msg); end
    end
end

% =====================================================================
function r=drawNewROI(ax,shape)
%drawNewROI  Interactive draw of one ROI of the chosen shape on ax.
switch shape
    case 'rectangle'
        r=drawrectangle(ax,'Color','w','FaceAlpha',0);
    case 'square'
        r=drawrectangle(ax,'Color','w','FaceAlpha',0,'FixedAspectRatio',true);
        if isvalid(r)
            p=r.Position; side=max(p(3),p(4)); r.Position=[p(1),p(2),side,side];
        end
    case 'ellipse'
        r=drawellipse(ax,'Color','w','FaceAlpha',0);
    case 'circle'
        r=drawcircle(ax,'Color','w','FaceAlpha',0);
    otherwise   % polygon (default)
        r=drawpolygon(ax,'Color','w','FaceAlpha',0);
end
end

% =====================================================================
function spec=captureROI(r)
%captureROI  Snapshot one live ROI's geometry into a carry-forward spec struct.
spec=emptyROISpec(); spec(1).type='polygon';
switch class(r)
    case 'images.roi.Rectangle'
        if r.FixedAspectRatio, spec.type='square'; else, spec.type='rectangle'; end
        spec.a=r.Position;
    case 'images.roi.Ellipse'
        spec.type='ellipse'; spec.a=r.Center; spec.b=r.SemiAxes; spec.c=r.RotationAngle;
    case 'images.roi.Circle'
        spec.type='circle';  spec.a=r.Center; spec.b=r.Radius;
    otherwise                                  % images.roi.Polygon (default)
        spec.type='polygon'; spec.a=r.Position;
end
end

% =====================================================================
function r=recreateROI(ax,spec)
%recreateROI  Recreate an editable ROI on ax from a carry-forward spec.
switch spec.type
    case 'rectangle'
        r=drawrectangle(ax,'Position',spec.a,'Color','w','StripeColor','b','FaceAlpha',0);
    case 'square'
        r=drawrectangle(ax,'Position',spec.a,'Color','w','StripeColor','b', ...
            'FaceAlpha',0,'FixedAspectRatio',true);
    case 'ellipse'
        r=drawellipse(ax,'Center',spec.a,'SemiAxes',spec.b,'RotationAngle',spec.c, ...
            'Color','w','StripeColor','b','FaceAlpha',0);
    case 'circle'
        r=drawcircle(ax,'Center',spec.a,'Radius',spec.b,'Color','w','StripeColor','b','FaceAlpha',0);
    otherwise   % polygon
        r=drawpolygon(ax,'Position',spec.a,'Color','w','StripeColor','b','FaceAlpha',0);
end
end

% =====================================================================
function spec=emptyROISpec()
%emptyROISpec  0x0 struct array with the ROI-spec fields (type + geometry a/b/c).
spec=struct('type',{},'a',{},'b',{},'c',{});
end
