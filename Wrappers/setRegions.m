%setRegions  Interactive, per-file selection of segmentation regions (regionsMask)
%
%   setRegions(s,fNames) opens an ROI editor for every file in fNames and writes the
%   labelled region mask (results.regionsMask) that runSegmentation / splitRegions
%   later use to restrict processing to user-chosen parts of the field of view.  The
%   regions are DRAWN, never computed - the only thing a protocol fixes in advance is
%   HOW MANY of them there may be (s.nRegions, default unlimited).  fNames is
%   TWO-DIMENSIONAL (rows = groups, e.g. one animal / FOV per row, exactly as
%   getFileNamesList returns when grouped) and setRegions iterates every file itself -
%   there is no launcher for-loop.
%
%   CARRY-FORWARD WITHIN A GROUP, RESET ACROSS GROUPS.  For each row, files are visited
%   left to right.  On the first file the user draws the region(s); on every subsequent
%   file the SAME ROIs are re-offered as editable objects (nudge / add / delete / draw
%   more).  At the next row the ROIs reset to empty.  A carried ROI that does not fit
%   entirely inside the next file's field of view (files of a group need not share an
%   image size) is DROPPED rather than re-offered off-image, and no more than
%   s.nRegions of the survivors are re-offered, so the editor never opens over budget.
%
%   THE EDITOR (one window per file).  ROIs are drawn on the display-enhanced image
%   (enhanceForDisplay), in the native image orientation so createMask aligns with the
%   stored frame.  Nothing advances until Done is pressed, so regions can be freely
%   adjusted first.  Controls:
%     * a shape selector - polygon / rectangle / square / ellipse / circle
%       (drawpolygon / drawrectangle (+FixedAspectRatio for square) / drawellipse /
%       drawcircle);
%     * Add ROI     - draw one more region of the selected shape.  DISABLED once
%                     s.nRegions regions are on the image, and enabled again as soon
%                     as one is deleted;
%     * Delete ROI  - remove the selected ROI (or, if none is selected, the most recent
%                     one).  The Delete key also removes the currently selected ROI(s);
%     * Done        - accept the current ROIs, save, and advance to the next file.
%
%   WHOLE WINDOW == NO MASK.  If the user draws NO ROIs (or skips setRegions entirely),
%   NO results.regionsMask is written - any stale one is removed.  Downstream readers
%   (runSegmentation, splitRegions) treat a MISSING regionsMask as the whole window (an
%   all-true mask the size of the image), so "draw nothing" == "use everything".  This
%   is true whatever s.nRegions says: the field is a CEILING, never a quota.  With
%   k ROIs drawn, regionsMask = sum_k createMask(roi_k).*k (k = the ROI index), identical
%   to the old runCategories math, so downstream (splitRegions, getPixelCategories) is
%   unchanged.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.setRegions untouched.
%       .nRegions      (optional, default Inf) the largest number of regions that may be
%                      drawn on one file.  A positive integer scalar, or Inf / [] for
%                      unlimited - which is what a launcher that never sets it gets, so
%                      the field changes nothing until it is asked for.
%       .fNamesCopyTo  (optional, default {}) copies the drawn mask onto co-registered
%                      siblings - see below.
%     fNames   2-D cell array of *_K_d.mat / *_I_d.mat paths.  Rows = groups; each file
%              must have matching *_s.mat and *_r.mat siblings.  Empty cells are skipped
%              (ragged rows from getFileNamesList are fine).
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail) and
%     s.cancelFcn()->tf (checked between files).
%
%   s.fNamesCopyTo - draw ONCE, apply to every branch of the same recording
%     A recording usually exists as several co-registered products of ONE raw file -
%     the contrast '_t' (or '_s'), the internal cycle '_c', the external cycle '_e'.
%     They share the field of view, so the regions should be drawn once and inherited,
%     exactly as runSegmentation inherits its cMask.  Two shapes are accepted:
%       * ELEMENTWISE - a cell array THE SAME SIZE AS fNames, where element (g,c) holds
%         the target(s) inheriting the mask drawn on fNames{g,c}, as one path (char) or
%         several (nested cellstr).  This is the natural form for a grouped fNames, e.g.
%         s.fNamesCopyTo = regexprep(fNames,'_t_K_d.mat$','_c_K_d.mat').
%       * ROW-PER-SOURCE - the runSegmentation convention: row i lists the targets for
%         the i-th file of fNames in its own (column-major) order.
%     Empty entries copy nowhere.  A target inherits results.regionsMask verbatim (and
%     has a stale one REMOVED when nothing was drawn), plus settings.setRegions; nothing
%     else on the target is touched.  Targets are NOT opened in the editor.
%
%   SIDE-EFFECTS (per file, and per copy target)
%     <name>_r.mat   results.regionsMask added/overwritten when >=1 ROI is drawn, or
%                    REMOVED when none is drawn (double, 0 = excluded, 1..N = labels)
%     <name>_s.mat   settings.setRegions = s
%     <name>_rep_regions.jpg  the record of what was drawn: the displayed image with
%                    the region boundaries outlined and numbered by label index.
%                    Written for EVERY file that receives a mask, drawn or copied,
%                    and also when nothing was drawn ("the whole window is used" is
%                    a decision worth recording).
%
%   EXAMPLE
%     % draw on the temporal contrast, inherit onto the paired internal-cycle files
%     fNames = getFileNamesList(root,'*_t_K_d.mat','[A-Z]+\d+');   % grouped, rows=animals
%     s.fNamesCopyTo = regexprep(fNames,'_t_K_d.mat$','_c_K_d.mat');
%     s.nRegions     = 1;      % one region per recording (omit for as many as you like)
%     setRegions(s,fNames);
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
% Last revision: 31-July-2026

function setRegions(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_K_d.mat')|| contains(x,'_I_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat" or "_I_d.mat".');
end
if ~isfield(s,'fNamesCopyTo'), s.fNamesCopyTo={}; end
s.nRegions=validNRegions(s);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Regions',fNames(~cellfun(@isempty,fNames(:))));

nGroups=size(fNames,1);
fIdx=0;
for g=1:1:nGroups
    carried=emptyROISpec();          % ROI geometry carried within this group (reset per group)
    for c=1:1:size(fNames,2)
        if reportCancelled(rep), return; end  % cooperative cancel between files (across groups)
        fName=fNames{g,c};
        if isempty(fName), continue; end
        fIdx=fIdx+1;
        reportFile(rep,fIdx,fName);
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
        [regionsMask,carried]=editRegions(imgIni,isK,carried,s.nRegions);

        if isempty(regionsMask)
            % no regions drawn: whole window - remove any stale mask, write none, so
            % downstream falls back to an all-true mask (see runSegmentation/splitRegions).
            if isfield(results,'regionsMask'), results=rmfield(results,'regionsMask'); end
        else
            results.regionsMask=regionsMask;
        end
        % The record of what was drawn and confirmed - written for EVERY file,
        % including the ones where nothing was drawn: "the whole window was used"
        % is as much a decision as any ROI.
        writeRegionsReport(rep,fName,imgIni,regionsMask);

        settings.setRegions=reportSettings(s);
        reportWriting(rep);
        save(strrep(fName,'_d.mat','_s.mat'),'settings','-v7.3','-nocompression');
        save(strrep(fName,'_d.mat','_r.mat'),'results','-v7.3','-nocompression');
        reportSaved(rep);

        % --- inherit the same regions on the co-registered siblings (s.fNamesCopyTo) ---
        tgts=copyTargets(s,fNames,g,c);
        for t=1:1:numel(tgts)
            if ~isempty(tgts{t})
                copyRegionsOnto(s,tgts{t},regionsMask,rep,fIdx);
            end
        end
    end
end
reportClose(rep);
end

% =====================================================================
function n=validNRegions(s)
%validNRegions  s.nRegions as a usable ceiling: a positive integer, or Inf for
%   unlimited.  ABSENT OR EMPTY MEANS UNLIMITED, which is what every launcher that
%   predates the field gets, and Inf rather than [] is the point - the editor then
%   compares with a plain >= and needs no special case for "no limit".
if ~isfield(s,'nRegions') || isempty(s.nRegions), n=Inf; return; end
n=s.nRegions;
% Stated as what a ceiling IS, not as the ways of being wrong: NaN passes every
% negative test ever written for it (NaN<1 is false, so is NaN>1), and a NaN ceiling
% would then compare false against every count and silently mean "unlimited".
ok=isnumeric(n) && isscalar(n) && isreal(n) && n>=1 && (isinf(n) || mod(n,1)==0);
if ~ok
    error(['s.nRegions must be a positive whole number (the largest number of ' ...
           'regions that may be drawn on one file), or Inf / [] for unlimited.']);
end
n=double(n);
end

% =====================================================================
function tgts=copyTargets(s,fNames,g,c)
%copyTargets  The copy targets for source file (g,c) as a cellstr, or {}.  Two accepted
%   shapes (see the header): ELEMENTWISE, s.fNamesCopyTo the same size as the 2-D fNames
%   with one path or a nested cellstr per element - the natural form for the grouped
%   launcher list; or the runSegmentation ROW form, one row of targets per source file
%   in fNames' own (column-major) order, which is what a flat 1-file call produces.
tgts={};
if ~isfield(s,'fNamesCopyTo') || isempty(s.fNamesCopyTo), return; end
ct=s.fNamesCopyTo;
if isequal(size(ct),size(fNames))
    e=ct{g,c};
    if isempty(e), return; end
    if ischar(e) || isstring(e), tgts={char(e)}; else, tgts=e(:)'; end
elseif size(ct,1)==numel(fNames)
    row=ct(sub2ind(size(fNames),g,c),:);
    tgts=row(~cellfun(@isempty,row));
end
end

% =====================================================================
function copyRegionsOnto(s,targetName,regionsMask,rep,fIdx)
%copyRegionsOnto  Give a co-registered sibling the SAME regions (verbatim mask, or the
%   removal of a stale one when nothing was drawn) plus the settings stamp.  Nothing
%   else on the target is touched - it is a different recording of the same FOV.
%   It is a RECORDING OF ITS OWN, so it gets its own three lines.
reportFile(rep,fIdx,targetName);
clearvars results settings
load(strrep(targetName,'_d.mat','_s.mat'),'settings');
load(strrep(targetName,'_d.mat','_r.mat'),'results');

if isempty(regionsMask)
    if isfield(results,'regionsMask'), results=rmfield(results,'regionsMask'); end
else
    results.regionsMask=regionsMask;
end

% The sibling gets the mask verbatim, but it is a DIFFERENT recording: its own
% report page is the only way to see whether the mask still fits its field of view.
% The background is the target's own mean image, built exactly as the editor built
% the source's.
writeRegionsReport(rep,targetName,targetImage(targetName,results),regionsMask);

sT=s; sT.fName=targetName;
settings.setRegions=reportSettings(sT);
reportWriting(rep);
save(strrep(targetName,'_d.mat','_s.mat'),'settings','-v7.3','-nocompression');
save(strrep(targetName,'_d.mat','_r.mat'),'results','-v7.3','-nocompression');
reportSaved(rep);
end

% =====================================================================
function imgIni=targetImage(targetName,results)
%targetImage  The copy target's own mean image - results.imgK / results.imgI when
%   they are there, the source cube's time-mean when they are not.  Same rule the
%   main loop uses for the file being edited.
if contains(targetName,'_K_d.mat')
    if isfield(results,'imgK')
        imgIni=results.imgK;
    else
        src=load(targetName,'source');
        imgIni=mean(src.source.data,3,'omitmissing');
    end
else
    imgIni=results.imgI;
end
end

% =====================================================================
function writeRegionsReport(rep,fName,imgIni,regionsMask)
%writeRegionsReport  Save <name>_rep_regions.jpg: the displayed image with the drawn
%   regions outlined and numbered by label index.
%
%   BOUNDARIES, NOT A FILLED OVERLAY, so the image underneath stays readable - the
%   point of the page is to check that the region sits where the operator meant it
%   to, which needs the vessels under it to be visible.  The number is the label
%   index, so the page and results.regionsMask can be compared directly.
%   A FAILURE HERE COSTS NOTHING BUT THE PAGE.  This page is new: before it existed
%   setRegions saved its mask and moved on, and it must still do that - silently -
%   if the drawing goes wrong on some image nobody anticipated.
f=reportFigure(rep,'regions','single');
try
    t=tiledlayout(f,1,1,'TileSpacing','compact','Padding','compact');
    ax=nexttile(t);
    imagesc(ax,imgIni); axis(ax,'image'); colormap(ax,'gray')
    setClim(ax,prctile(imgIni(:),[1,99]))
    if isempty(regionsMask)
        title(ax,'No regions drawn - the whole window is used')
    else
        hold(ax,'on')
        for k=1:1:max(regionsMask(:))
            bw=(regionsMask==k);
            if ~any(bw(:)), continue; end
            visboundaries(ax,bw,'Color','r','LineWidth',1);
            [y,x]=find(bw);
            text(ax,mean(x),mean(y),num2str(k),'Color','y','FontWeight','bold', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
        hold(ax,'off')
        title(ax,[num2str(max(regionsMask(:))) ' region(s)'])
    end
catch
    delete(f);              % no page, no line: a report is a by-product
    return
end
reportSave(rep,f,'regions',fName);
end

% =====================================================================
function setClim(ax,lims)
%setClim  clim() with the degenerate case handled: a flat image gives equal
%   percentiles, and clim rejects a non-increasing pair.
lims=double(lims);
if any(~isfinite(lims)) || lims(2)<=lims(1), return; end
clim(ax,lims);
end

% =====================================================================
function [regionsMask,carried]=editRegions(imgIni,isK,carried,nRegions)
%editRegions  ROI editor for one file: draw/edit labelled regions on the enhanced
%   image, return the labelled mask (EMPTY when no ROI is drawn = whole window) plus
%   the ROI geometry to carry to the next file.  nRegions is the ceiling (Inf =
%   unlimited): Add ROI is greyed out at the limit, and drawing FEWER is always
%   allowed - the limit says how many the operator MAY draw, not how many they must.

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
% Two buttons now that Reset ROIs is gone, sharing the width the three used to have -
% Done keeps its own place on the right, where it has always been held apart.
addBtn=uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.29 0.105 0.13 0.045], ...
    'String','Add ROI','Callback',@onAdd);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.43 0.105 0.13 0.045], ...
    'String','Delete ROI','Callback',@onDelete);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.76 0.10 0.13 0.055], ...
    'String','Done (next file)','FontWeight','bold','Callback',@onDone);

hint=uicontrol(fig,'Style','text','Units','normalized','Position',[0.04 0.02 0.68 0.055], ...
    'String','','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);

% Pre-populate the carried ROIs (editable) - empty on the first file of a group.
% Files within a group need not share a field of view, so an ROI drawn on a larger
% image can fall (partly) outside a smaller one; such an ROI is dropped rather than
% recreated off-image, where it would silently mask nothing and break downstream use.
% The ceiling applies to what is carried in as much as to what is drawn: a group whose
% first file was given three ROIs must not re-open at a limit of one already over
% budget, with an Add button greyed out and no explanation.  Truncation is AFTER the
% out-of-FOV drop, so the ROIs that survive are the first nRegions that still fit.
carried=dropOutsideFOV(carried,size(imgIni));
if numel(carried)>nRegions, carried=carried(1:nRegions); end
for i=1:1:numel(carried)
    rois{end+1}=recreateROI(axDraw,carried(i)); %#ok<AGROW>
end

fig.UserData='editing';                       % private wait sentinel (see waitfor below)
fig.WindowKeyPressFcn=@onKey;                 % Delete key removes the selected ROI(s)
fig.CloseRequestFcn=@(~,~) set(fig,'UserData','done');   % closing the window == Done
refreshControls();                            % hint + Add's enable state, from the ROI count
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

    % ---- nested callbacks (share rois / axDraw / shapePop / hint / addBtn) ----
    function onAdd(~,~)
        if nLive()>=nRegions, return; end   % the button is greyed, but never trust that alone
        shapes={'polygon','rectangle','square','ellipse','circle'};
        r=drawNewROI(axDraw,shapes{shapePop.Value});
        if ~isempty(r) && isvalid(r)
            rois{end+1}=r;
        end
        refreshControls();
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
        refreshControls();
    end
    function onKey(~,evt)
        if strcmp(evt.Key,'delete')                % Delete key = remove selected ROI(s)
            deleted=false;
            for jj=numel(rois):-1:1
                if isvalid(rois{jj}) && rois{jj}.Selected
                    delete(rois{jj}); rois(jj)=[]; deleted=true;
                end
            end
            if deleted, refreshControls(); end
        end
    end
    function onDone(~,~)
        if isvalid(fig), fig.UserData='done'; end
    end
    function n=nLive()
        %nLive  THE ROI count - the number of live handles, computed in ONE place so
        %   the hint and the Add button can never disagree about whether the file is
        %   full.  Deleting an ROI leaves an invalid handle behind until the callback
        %   that removed it also drops the cell, so counting cells is not the same
        %   thing as counting regions.
        n=0;
        for jj=1:1:numel(rois), if isvalid(rois{jj}), n=n+1; end, end
    end
    function refreshControls()
        %refreshControls  Everything the ROI count drives: the hint and Add's enable
        %   state.  Called wherever an ROI appears or disappears.
        n=nLive();
        full=n>=nRegions;
        if isvalid(addBtn), set(addBtn,'Enable',onOff(~full)); end
        if n==0
            msg=['No ROIs drawn: the WHOLE WINDOW will be used (no region mask is ', ...
                 'written).  Add ROIs to restrict segmentation, then press Done.'];
        elseif full
            msg=[num2str(n) ' ROI(s) - as many as this file may have.  Drag / reshape ', ...
                 'freely; delete one to draw a different one.  Press Done when finished.'];
        else
            msg=[num2str(n) ' ROI(s).  Add / drag / reshape freely; select an ROI and ', ...
                 'press Delete (or "Delete ROI") to remove it.  Press Done when finished.'];
        end
        if isvalid(hint), set(hint,'String',msg); end
    end
end

% =====================================================================
function v=onOff(tf)
%onOff  A logical as the 'on'/'off' a uicontrol Enable expects.
if tf, v='on'; else, v='off'; end
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
    otherwise                                      % polygon (default)
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
    otherwise                                      % images.roi.Polygon (default)
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
    otherwise                                      % polygon
        r=drawpolygon(ax,'Position',spec.a,'Color','w','StripeColor','b','FaceAlpha',0);
end
end

% =====================================================================
function carried=dropOutsideFOV(carried,sz)
%dropOutsideFOV  Keep only the carried ROI specs that fit entirely inside an sz =
%   [rows cols] field of view, discarding the ones that do not.  ROIs are carried
%   across files of a group, which may differ in size; an ROI that does not fit is
%   deleted here so nothing downstream sees a region hanging off the image.  The
%   operator sees the outcome on the image in front of them, which is why this is
%   silent.
keep=true(1,numel(carried));
for i=1:1:numel(carried)
    [xLim,yLim]=roiExtent(carried(i));
    keep(i)= xLim(1)>=0.5 && xLim(2)<=sz(2)+0.5 && yLim(1)>=0.5 && yLim(2)<=sz(1)+0.5;
end
carried=carried(keep);
end

% =====================================================================
function [xLim,yLim]=roiExtent(spec)
%roiExtent  Bounding box [min max] in x (columns) and y (rows) of one ROI spec, in the
%   image coordinates ROIs are drawn in (pixel centres at 1..N, so the image spans
%   0.5..N+0.5).  Ellipses use the bounding box of the ROTATED ellipse.
switch spec.type
    case {'rectangle','square'}
        p=spec.a;                                  % [x y w h]
        xLim=[p(1),p(1)+p(3)]; yLim=[p(2),p(2)+p(4)];
    case 'ellipse'
        th=deg2rad(spec.c); a=spec.b(1); b=spec.b(2);
        hx=hypot(a*cos(th),b*sin(th)); hy=hypot(a*sin(th),b*cos(th));
        xLim=spec.a(1)+[-hx,hx]; yLim=spec.a(2)+[-hy,hy];
    case 'circle'
        xLim=spec.a(1)+[-spec.b,spec.b]; yLim=spec.a(2)+[-spec.b,spec.b];
    otherwise                                      % polygon: [x y] vertices
        xLim=[min(spec.a(:,1)),max(spec.a(:,1))];
        yLim=[min(spec.a(:,2)),max(spec.a(:,2))];
end
end

% =====================================================================
function spec=emptyROISpec()
%emptyROISpec  0x0 struct array with the ROI-spec fields (type + geometry a/b/c).
spec=struct('type',{},'a',{},'b',{},'c',{});
end
