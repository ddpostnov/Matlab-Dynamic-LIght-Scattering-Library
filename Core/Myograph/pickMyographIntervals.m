%pickMyographIntervals  Interactive GUI to define / edit time intervals on a diameter trace
%
%   [ivT,names] = pickMyographIntervals(time,d)
%   [ivT,names] = pickMyographIntervals(time,d,ivT0,names0)
%
%   Shows the diameter (median over Y if d is 2-D) versus time on a full-width,
%   zoomable / pannable plot.  Intervals are drawn as editable shaded bands:
%     * ADD    - drag the yellow box over a range, type a title, press
%                "Add interval".
%     * EDIT   - drag a band's left/right edge (or its body) to move its
%                start/end - the band stays full height.
%     * RENAME - double-click a band and type a new title.
%     * DELETE - toggle "Delete interval", then click a band to remove it.
%
%   Pass ivT0 = [start end; ...] (seconds) and names0 (cell of titles) to
%   PRE-LOAD intervals - e.g. the ones defined manually as in OPTION A - so you
%   only have to fine-tune their edges / titles.  Press "Done" (or close the
%   window) to return.
%
%   Navigation: the axes toolbar (zoom / pan / restore) or the scroll wheel.
%
%   OUTPUT
%     ivT    [start end; ...] intervals in seconds, sorted by start time
%     names  1xN cell of the matching titles
%   (same format as the manual intervalsPerFile / intervalNamesPerFile entries).
%
%   DEPENDS ON
%     Image Processing Toolbox (drawrectangle) and base MATLAB (R2023a+).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

function [ivT,names] = pickMyographIntervals(time,d,ivT0,names0)

time=double(time(:));
if ~isvector(d), d=median(double(d),2,'omitnan'); end
d=double(d(:));
if nargin<3 || isempty(ivT0), ivT0=zeros(0,2); end
if nargin<4, names0={}; end

fig=figure('Name','Select intervals','NumberTitle','off','Color','w', ...
    'Units','normalized','Position',[0.08 0.20 0.84 0.62]);
ax=axes('Parent',fig,'Units','normalized','Position',[0.07 0.24 0.90 0.70]);
plot(ax,time,d,'-','Color',[0 0.3 0.8]); grid(ax,'on'); box(ax,'on');
xlabel(ax,'time (s)'); ylabel(ax,'diameter (px)');
title(ax,helpStr);
if time(end)>time(1), xlim(ax,[time(1) time(end)]); end
yl=ylim(ax); span=max(time(end)-time(1),eps);

% ---- state ----
S.ax=ax; S.yl=yl; S.roiList={}; S.deleteMode=false; S.faceColor=[0.2 0.6 1];

% ---- yellow selection box (for adding new intervals) ----
S.selRoi=drawrectangle(ax,'Position',[time(1)+0.02*span yl(1) 0.10*span diff(yl)], ...
    'Color',[0.95 0.75 0],'FaceAlpha',0.12,'LineWidth',1,'Label','new');
addlistener(S.selRoi,'ROIMoved',@(s,~) snapY(s,yl));

% ---- controls ----
uicontrol(fig,'Style','text','Units','normalized','Position',[0.07 0.075 0.06 0.045], ...
    'String','Title:','BackgroundColor','w','HorizontalAlignment','right','FontSize',10);
S.edit=uicontrol(fig,'Style','edit','Units','normalized','Position',[0.14 0.075 0.28 0.05], ...
    'String','','FontSize',10,'HorizontalAlignment','left');
S.addBtn=uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.44 0.075 0.15 0.05], ...
    'String','Add interval','FontSize',10,'Callback',@(~,~) localAdd(fig));
S.delBtn=uicontrol(fig,'Style','togglebutton','Units','normalized','Position',[0.61 0.075 0.15 0.05], ...
    'String','Delete interval','FontSize',10,'Callback',@(src,~) localToggleDelete(fig,src));
S.doneBtn=uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.84 0.075 0.13 0.05], ...
    'String','Done','FontSize',10,'FontWeight','bold','Callback',@(~,~) uiresume(fig));

fig.CloseRequestFcn=@(~,~) uiresume(fig);
guidata(fig,S);

% ---- pre-load user-defined intervals as editable bands ----
if ~isempty(ivT0)
    if isempty(names0)
        names0=arrayfun(@(k)sprintf('interval%d',k),1:size(ivT0,1),'UniformOutput',false);
    end
    for k=1:size(ivT0,1)
        addCommitted(fig,ivT0(k,1),ivT0(k,2),names0{k});
    end
end

uiwait(fig);

% ---- collect results from the committed bands (sorted by start) ----
ivT=zeros(0,2); names={};
if isvalid(fig)
    S=guidata(fig);
    for i=1:numel(S.roiList)
        r=S.roiList{i};
        if isvalid(r)
            p=r.Position;
            ivT(end+1,:)=[p(1) p(1)+p(3)]; %#ok<AGROW>
            names{end+1}=char(r.Label);    %#ok<AGROW>
        end
    end
    if ~isempty(ivT), [ivT,ord]=sortrows(ivT); names=names(ord); end
    delete(fig);
end
end

% =====================================================================
function roi=addCommitted(fig,x1,x2,name)
S=guidata(fig);
yl=S.yl; xlo=min(x1,x2); w=abs(x2-x1);
roi=drawrectangle(S.ax,'Position',[xlo yl(1) w diff(yl)], ...
    'Color',S.faceColor,'FaceAlpha',0.12,'LineWidth',1, ...
    'Label',char(name),'LabelVisible','on');
addlistener(roi,'ROIMoved',@(s,~) snapY(s,yl));
addlistener(roi,'ROIClicked',@(s,e) onRoiClicked(fig,s,e));
S.roiList{end+1}=roi;
guidata(fig,S);
end

% =====================================================================
function localAdd(fig)
S=guidata(fig);
if ~isvalid(S.selRoi), return; end
p=S.selRoi.Position; x1=p(1); x2=p(1)+p(3);
if x2<=x1, return; end
nm=strtrim(S.edit.String);
if isempty(nm), nm=sprintf('interval%d',numel(S.roiList)+1); end
addCommitted(fig,x1,x2,nm);
S.edit.String='';
end

% =====================================================================
function localToggleDelete(fig,src)
S=guidata(fig);
S.deleteMode=logical(src.Value);
if S.deleteMode
    src.BackgroundColor=[1 0.5 0.5];
    S.ax.Title.String='DELETE MODE: click a band to remove it.';
    if isvalid(S.selRoi), S.selRoi.Visible='off'; end
else
    src.BackgroundColor=[0.94 0.94 0.94];
    S.ax.Title.String=helpStr;
    if isvalid(S.selRoi), S.selRoi.Visible='on'; end
end
guidata(fig,S);
end

% =====================================================================
function onRoiClicked(fig,src,evt)
S=guidata(fig);
if S.deleteMode                                  % click removes the band
    idx=find(cellfun(@(r) isvalid(r)&&r==src,S.roiList));
    if ~isempty(idx), S.roiList(idx)=[]; end
    delete(src);
    guidata(fig,S);
    return;
end
st='';
try
    st=evt.SelectionType;
catch %#ok<CTCH>
end
if any(strcmpi(st,{'double','open'}))            % double-click renames
    a=inputdlg('Interval title:','Rename interval',1,{char(src.Label)});
    if ~isempty(a), src.Label=char(a{1}); end
end
end

% =====================================================================
function snapY(roi,yl)                            % keep bands full height (x-only edit)
p=roi.Position;
if numel(p)==4 && (p(2)~=yl(1) || p(4)~=diff(yl))
    roi.Position=[p(1) yl(1) p(3) diff(yl)];
end
end

% =====================================================================
function s=helpStr()
s=['Drag box + "Add interval" to create | drag a band''s edge to adjust | ' ...
   'double-click to rename | "Delete interval" toggle + click to remove.'];
end
