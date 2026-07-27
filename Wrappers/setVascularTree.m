%setVascularTree  Derive & edit the vascular parent-daughter hierarchy
%
%   setVascularTree(s,fNames) is the logical post-processing step that
%   follows setVesselTypes.  For every *_BFI_d.mat dataset in fNames it
%   loads the companion *_BFI_r.mat, automatically derives a
%   parent->daughter vascular hierarchy over the segmented vessels and
%   parenchyma via the Core routine getVascularTree, stores the result in a
%   self-describing schema, and - unless
%   s.autoOnly is true - opens an interactive setVesselTypes-style GUI so
%   the derived hierarchy can be inspected and corrected before it is
%   written back.
%
%   METHOD (staged, type-constrained connection search)
%     A global flow potential phi (from pulse arrival psTimeMin/psTimeMax and
%     pulsatility psPI) orders every node from the arterial inlet (low) to
%     the venous outlet (high).  Connections are then found in stages:
%       (0) FOV BRIDGING - because the image is a 2D projection, a crossing
%           vessel of the opposite type can split one vessel in two so its
%           halves do not touch.  Same-type vessels whose pixels come within
%           s.bridgeTipRadius px of each other's TIP (the two geometric ends,
%           found by PCA) - or s.bridgeWallRadius px of each other's side wall
%           (default 0) - are added as candidate neighbours before the search.
%       (1) ARTERIAL tree - each artery links to the closest touching artery
%           just up-gradient in a side-specific potential that also uses the
%           caliber cues BFI/diameter (larger = more upstream on the arterial
%           side).  Arteries connect only to arteries here.
%       (2) VENOUS tree - each vein links to the closest touching vein just
%           down-gradient (larger vein = more downstream); veins only to veins.
%       (3) PARENCHYMA mesh - only TERMINAL vessels (an artery/vein whose only
%           same-type link is its parent, i.e. no smaller same-type child) may
%           touch parenchyma, one parenchyma each.  Arteries never connect
%           directly to veins - every arterio-venous path runs through the
%           parenchyma.  A parenchyma node links to up to parenchArteryDeg
%           arteries and parenchVeinDeg veins plus parenchParenchDeg other
%           parenchyma (parenchMaxDeg total), forming a capillary mesh from
%           terminal arterioles to terminal venules.  The tipParench parameter
%           (on by default) prefers a parenchyma touching the terminal vessel
%           near its FREE tip (opposite the tip that carries the parent/daughter
%           vessel); tipVessel (off) applies the same opposite-tip bias to the
%           vessel-tree links in stages 1-2.
%     Which parameters are used, their weights and roles (arrival / pulsatility
%     / caliber / tip) are configurable and editable live in the GUI.
%     Output: a flow-directed DAG with Horton-Strahler order and generation.
%
%   INPUTS
%     s        parameter struct
%                • autoOnly       true = derive & save without the GUI
%                                 (batch / headless).  Default false.
%                • flowParams     struct array of the monotone parameters
%                                 used to order the tree; fields name / dirn
%                                 (+1 increases, -1 decreases downstream) /
%                                 weight (0 = off) / enabled / label.  Default
%                                 (defaultFlowParams): psTimeMax, psTimeMin,
%                                 psPI enabled at weight 1.  Editable
%                                 live in the GUI.  BFI/diameter are U-shaped
%                                 and are not used to order the tree.  Two
%                                 role='tip' parameters (tipVessel off,
%                                 tipParench on) are geometric, not part of the
%                                 potential: they reward a connection landing
%                                 near the vessel tip OPPOSITE its committed
%                                 vessel link (see METHOD stage 3).
%                • connectivity   4 or 8 pixel adjacency.  Default 8.
%                • minBorder      min shared-border pixels for adjacency.
%                                 Default 1.
%                • useReference   true = do NOT derive; inherit the
%                                 hierarchy from s.refFName and map it onto
%                                 this file's (matching) segmentation.  The
%                                 reference file itself (s.fName==s.refFName)
%                                 is still derived + editable.  Use this to
%                                 carry the pulsatility-derived hierarchy to
%                                 a registered vasomotion "_t" partner.
%                • refFName       path to the reference *_BFI_d.mat (the
%                                 pulsatility "_c" file, whose *_BFI_r.mat
%                                 already holds results.hierarchy).
%                • propagatePartners  cellstr of recording-variant letters;
%                                 after a "_c" file is derived, the hierarchy
%                                 is auto-copied to sibling recordings whose
%                                 name has _c_BFI replaced by _<L>_BFI, when
%                                 they exist.  Default {'t','s'}; {} disables.
%     fNames   cell array of *_BFI_d.mat paths.
%
%   SIDE-EFFECTS (per file)
%     *_BFI_r.mat   RESULTS – new sMetrics columns parentIdx, daughterIdx
%                   (cell), treeID, generation, strahlerOrder,
%                   hierarchyConfidence, flowPotential, isRoot, isOutlet;
%                   plus RESULTS.hierarchy (edge list / adjacency / roots /
%                   params / version) and RESULTS.mapTree overlay image.
%     *_BFI_s.mat   SETTINGS – field settings.setVascularTree added.
%
%   EXAMPLE
%     s.autoOnly = false;
%     D = dir(fullfile(dataRoot,'*_c_BFI_d.mat'));
%     setVascularTree(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     getVascularTree (Core: automatic derivation), orderForest, getMetric &
%     defaultFlowParams (Core); enhanceForDisplay (background-subtracted GUI
%     preview); Image Processing Toolbox (imdilate/medfilt2), Statistics Toolbox
%     (tiedrank) and MATLAB graph objects (digraph/conncomp/toposort/isdag); core
%     LSCI library utilities.  Consumes the output of runSegmentation,
%     runPulsatility and setVesselTypes.
%
% See also: getVascularTree, setVesselTypes, runPulsatility, runSegmentation, runVasomotion
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - HIERARCHY DERIVATION
% s.autoOnly=false; % true to derive & save without opening the GUI
% %DETECTION PARAMETERS - easiest to tune live in the GUI. Leave s.flowParams
% %unset to use the defaults (psTimeMax,psTimeMin,psPI at weight 1).
% %To trust peak timing more and PI less (build the struct explicitly):
% % s.flowParams=struct('name',{'psTimeMax','psTimeMin','psPI'},'label',{'peak','foot','PI'},'dirn',{1,1,-1},'weight',{2,1,0.5},'enabled',{true,true,true});
% s.connectivity=8; % 4 or 8 pixel adjacency between segments
% s.minBorder=1; % minimum shared-border pixels to treat two segments as touching
% %FOV BRIDGING - reconnect same-type vessels split by a crossing vessel in the
% %2D projection.  Search this many pixels from the vessel tips / walls:
% s.bridgeTipRadius=50; % px gap searched from the vessel ends (0 = off)
% s.bridgeWallRadius=0; % px gap searched from the vessel sides (usually 0)
% %SET FILE NAMES HERE
% fNames=getFileNamesList(rootFolder,'*_c_BFI_d.mat');
% setVascularTree(s,fNames(:));
% %OPTIONAL - carry the hierarchy to the registered vasomotion "_t" partner
% %(no pulse timing there, so inherit instead of deriving):
% % s.useReference=true; s.refFName=fNames{1}; % the matching _c_BFI_d.mat
% % setVascularTree(s, getFileNamesList(rootFolder,'*_t_BFI_d.mat'));

function setVascularTree(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_BFI_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_BFI_d.mat".');
end

% ---- defaults ----------------------------------------------------------
if ~isfield(s,'autoOnly'),         s.autoOnly=false;        end
if ~isfield(s,'flowParams')||isempty(s.flowParams), s.flowParams=defaultFlowParams(); end
if ~isfield(s,'connectivity'),     s.connectivity=8;        end
if ~isfield(s,'minBorder'),        s.minBorder=1;           end
% Parenchyma connection degree caps (adjustable, split by vessel type):
% a parenchyma links to <=parenchArteryDeg arteries and <=parenchVeinDeg
% veins (its capillary feed/drain) plus <=parenchParenchDeg other parenchyma,
% with <=parenchMaxDeg connections total.
if isfield(s,'parenchVesselDeg') && ~isfield(s,'parenchArteryDeg') && ~isfield(s,'parenchVeinDeg')
    s.parenchArteryDeg=max(1,floor(s.parenchVesselDeg/2));   % split a legacy combined vessel cap
    s.parenchVeinDeg  =max(1,floor(s.parenchVesselDeg/2));
end
if ~isfield(s,'parenchArteryDeg'), s.parenchArteryDeg=1;    end  % max arteries a parenchyma links to
if ~isfield(s,'parenchVeinDeg'),   s.parenchVeinDeg=1;      end  % max veins a parenchyma links to
if ~isfield(s,'parenchParenchDeg'),s.parenchParenchDeg=4;   end  % max other parenchyma it links to
if ~isfield(s,'parenchMaxDeg'),    s.parenchMaxDeg=6;       end  % max total connections of a parenchyma
% FOV-bridging: a crossing vessel of the opposite type can split one vessel in
% two in the 2D projection.  Link same-type vessels whose pixels come within
% bridgeTipRadius px of each other's TIP (the geometric ends) or bridgeWallRadius
% px of each other's side wall.  A real split usually happens at the ends, so
% the wall search is off by default.
if ~isfield(s,'bridgeTipRadius'),  s.bridgeTipRadius=50;    end  % px gap searched from vessel tips
if ~isfield(s,'bridgeWallRadius'), s.bridgeWallRadius=0;    end  % px gap searched from vessel walls
if ~isfield(s,'bridgeTipBand'),    s.bridgeTipBand=0.2;     end  % border fraction near each end = "tip"
if ~isfield(s,'useReference'),     s.useReference=false;    end
if ~isfield(s,'refFName'),         s.refFName='';           end
if ~isfield(s,'propagatePartners'),s.propagatePartners={'t','s'}; end

for fidx=1:1:numel(fNames)
    if ~isempty(fNames{fidx})
        tic
        disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        s.fName=fNames{fidx};
        clearvars results settings
        load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
        load(strrep(s.fName,'_d.mat','_r.mat'),'results');

        if s.useReference && ~isempty(s.refFName) && ~strcmp(s.fName,s.refFName)
            % ---- reference mode: inherit hierarchy from the reference ----
            % (the reference must share the segmentation - guaranteed within a
            % registered group; parent/daughter ids are segment idx, valid as-is)
            disp('  Reference mode: propagating hierarchy from reference file.');
            ref=load(strrep(s.refFName,'_d.mat','_r.mat'),'results');
            if ~isfield(ref.results,'hierarchy')
                error('Reference file has no results.hierarchy - run setVascularTree on the reference first.');
            end
            Href=ref.results.hierarchy;
            if ~isequal(size(results.sMap),Href.imgSize)
                error('Reference segmentation size differs from this file - propagation needs a matched segmentation.');
            end
            results=applyHierarchy(results,Href);
            disp(['  Inherited ',num2str(numel(Href.nodeIds)),' nodes, ', ...
                num2str(size(Href.edges,1)),' edges from reference. Elapsed ',num2str(round(toc)),'s'])
        else
            % ---- automatic derivation (headless-testable) ---------------
            H = getVascularTree(results,s);
            results = applyHierarchy(results,H);
            disp(['  Auto-derived ',num2str(numel(H.nodeIds)),' nodes, ', ...
                num2str(size(H.edges,1)),' edges, ',num2str(numel(H.roots)), ...
                ' roots, ',num2str(numel(H.outlets)),' outlets. Elapsed ', ...
                num2str(round(toc)),'s'])

            % ---- interactive correction ---------------------------------
            if ~s.autoOnly
                results = vascularTreeGUI(results,s);
            end
        end

        settings.setVascularTree=s;
        disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
        save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');

        % auto-propagate the hierarchy to registered partner recordings
        % (same base name, different variant letter: _c_BFI -> _t_BFI/_s_BFI)
        if ~s.useReference && isfield(results,'hierarchy') && contains(s.fName,'_c_BFI')
            propagateToPartners(results.hierarchy,s);
        end
    end
end
end


%% ====================  PARTNER PROPAGATION  ========================= %%
function propagateToPartners(H,s)
% Copy the hierarchy onto sibling recordings that share the segmentation
% (same base name, different variant letter, e.g. _c_BFI -> _t_BFI / _s_BFI).
parts=s.propagatePartners; if ischar(parts), parts={parts}; end
for i=1:numel(parts)
    L=parts{i};
    partnerD=strrep(s.fName,'_c_BFI',['_' L '_BFI']);
    if strcmp(partnerD,s.fName), continue; end
    partnerR=strrep(partnerD,'_d.mat','_r.mat');
    partnerS=strrep(partnerD,'_d.mat','_s.mat');
    if exist(partnerR,'file')~=2, continue; end
    disp(['  Auto-propagating hierarchy to partner: ',partnerR]);
    pr=load(partnerR,'results'); presults=pr.results;
    if ~isequal(size(presults.sMap),H.imgSize)
        warning('setVascularTree:partnerMismatch', ...
            'Partner segmentation differs (%s) - skipped.',partnerR); continue;
    end
    presults=applyHierarchy(presults,H);
    tmp.results=presults; save(partnerR,'-struct','tmp','-v7.3'); clear tmp
    if exist(partnerS,'file')==2
        ps=load(partnerS,'settings'); psettings=ps.settings;
        psettings.setVascularTree=s; psettings.setVascularTree.propagatedFrom=s.fName;
        tmp.settings=psettings; save(partnerS,'-struct','tmp','-v7.3'); clear tmp
    end
end
end


%% =====================  WRITE-BACK / SCHEMA  ======================== %%
function results = applyHierarchy(results,H)
% Writes the derived hierarchy into RESULTS.sMetrics columns, the
% RESULTS.hierarchy struct and a RESULTS.mapTree overlay image.  Vessel
% node values are mirrored onto the wall row (lumen idx k and wall idx
% k+1) so the map(sMap+1) coloring works for lumen and wall pixels.

M = results.sMetrics;
N = height(M);
cat = M.category;
nvi = M.nearestVesIdx;

% per-segment scalar/cell columns
parentIdx   = cell(N,1);
daughterIdx = cell(N,1);
treeID      = nan(N,1);
generation  = nan(N,1);
strahlerOrder = nan(N,1);
hierarchyConfidence = nan(N,1);
flowPotential = nan(N,1);
isRoot   = false(N,1);
isOutlet = false(N,1);

nid = H.nodeIds;
for c=1:numel(nid)
    k=nid(c);
    parentIdx{k}   = H.parentList{c};
    daughterIdx{k} = H.daughterList{c};
    treeID(k)      = H.treeID(c);
    generation(k)  = H.generation(c);
    strahlerOrder(k)= H.strahler(c);
    hierarchyConfidence(k)=H.confidence(c);
    flowPotential(k)= H.phi(c);
    isRoot(k)      = H.isRoot(c);
    isOutlet(k)    = H.isOutlet(c);
end

% mirror vessel-node values onto their wall rows (wall idx -> lumen via nvi)
wall = find(cat==3);
for w=wall'
    k = nvi(w);
    if k>=1 && k<=N
        parentIdx{w}=parentIdx{k};
        daughterIdx{w}=daughterIdx{k};
        treeID(w)=treeID(k);
        generation(w)=generation(k);
        strahlerOrder(w)=strahlerOrder(k);
        hierarchyConfidence(w)=hierarchyConfidence(k);
        flowPotential(w)=flowPotential(k);
        isRoot(w)=isRoot(k);
        isOutlet(w)=isOutlet(k);
    end
end

M.parentIdx=parentIdx;
M.daughterIdx=daughterIdx;
M.treeID=treeID;
M.generation=generation;
M.strahlerOrder=strahlerOrder;
M.hierarchyConfidence=hierarchyConfidence;
M.flowPotential=flowPotential;
M.isRoot=isRoot;
M.isOutlet=isOutlet;
results.sMetrics=M;
results.hierarchy=H;

% overlay image: default = Horton-Strahler order (like results.mapType)
val = zeros(N,1);
val(nid)=H.strahler;
for w=wall', if nvi(w)>=1&&nvi(w)<=N, val(w)=val(nvi(w)); end, end
map=[0;val];
map=map(double(results.sMap)+1);
map(results.cMask==0 | results.cMask==2)=NaN;   % show vessels + parenchyma (hide bg/unsegmented)
results.mapTree=map;
end


%% ============================  GUI  ================================= %%
function results = vascularTreeGUI(results,s)
% Interactive inspector/editor for the derived hierarchy, modelled on
% setVesselTypes.  Left panel: a background image (popup-selectable).  Right
% panel: the hierarchy coloured by a popup-selectable attribute (generation
% / Strahler / treeID / confidence / flowPotential / type).  Hover a segment
% to inspect its metrics, parents and daughters (and see its links drawn).
% Each node carries an up-link (parent, for arteries/parenchyma) and/or a
% down-link (daughter, for veins/parenchyma); editing sets one of these by
% clicking a segment then a touching up/down neighbour.  Finish saves.

H=results.hierarchy;
sMap=results.sMap; cMask=results.cMask;
nodeIds=H.nodeIds; nNodes=numel(nodeIds);
mx=double(max(sMap(:))); comp=zeros(mx,1); comp(nodeIds)=(1:nNodes)';

% pixel -> node id and -> compact node index
lut=[0; H.seg2node]; nodeImg=lut(double(sMap)+1);
compImg=zeros(size(nodeImg)); vpix=nodeImg>0; compImg(vpix)=comp(nodeImg(vpix));

% node centroids for drawing links
[yy,xx]=find(vpix); cc=compImg(vpix);
cx=accumarray(cc,xx,[nNodes 1],@mean,NaN);
cy=accumarray(cc,yy,[nNodes 1],@mean,NaN);

% phi-dependent editable state (also rebuilt on Re-derive)
[linkUp,linkDn,pEdges,adjC,autoConf]=linkStateFromH(H,comp,nNodes);
userSet=false(nNodes,1);

% per-node metrics for the inspector
M=results.sMetrics;
bfiN=getNode(M,'BFI',nodeIds); piN=getNode(M,'psPI',nodeIds);
footN=getNode(M,'psTimeMin',nodeIds); peakN=getNode(M,'psTimeMax',nodeIds);
diamN=getNode(M,'diameter',nodeIds);

% background image (as in setVesselTypes)
if isfield(results,'imgBFI'), img=sqrt(double(results.imgBFI)); else, img=double(cMask>0); end
lo=double(prctile(img(cMask(:)>0),[5,99])); if lo(1)==lo(2), lo=[min(img(:)) max(img(:))+eps]; end
img=mat2gray(img,lo);
fSize=floor((min(size(img))./20))*2+1;
img=enhanceForDisplay(img,fSize).*(cMask>0);

% ---- figure & layout ---------------------------------------------------
f=figure('Name','Vascular hierarchy','Color','w','WindowState','maximized', ...
    'CloseRequestFcn',@(~,~)uiresume(gcbf));
TL=tiledlayout(f,1,5,'TileSpacing','none','Padding','compact');
axL=nexttile(TL,[1 2]); imagesc(axL,img); axis(axL,'image','off'); colormap(axL,'parula');
try, clim(axL,prctile(img(cMask(:)>0),[1,99])); catch, end; title(axL,'BFI');
axR=nexttile(TL,[1 2]);
hMapR=imagesc(axR,nan(size(sMap))); axis(axR,'image','off'); title(axR,'Hierarchy');
axCtl=nexttile(TL); axis(axCtl,'off');

p=uipanel(f,'Units','pixels','BorderType','none','BackgroundColor','w');
guiReposition(f,axCtl,p,[260 700]);
f.SizeChangedFcn=@(~,~)guiReposition(f,axCtl,p,[260 700]);
lbl=@(y,str)uicontrol(p,'Style','text','String',str,'Units','normalized', ...
    'Position',[.04 y .92 .026],'HorizontalAlignment','left','BackgroundColor','w','FontWeight','bold');
btn=@(y,str,cb)uicontrol(p,'Style','pushbutton','String',str,'Units','normalized', ...
    'Position',[.05 y .9 .034],'Callback',cb);

lbl(.965,'Right overlay');
overlayPop=uicontrol(p,'Style','popup','Units','normalized','Position',[.05 .938 .9 .026], ...
    'String',{'strahler','generation','treeID','confidence','flowPotential','type'}, ...
    'Callback',@(~,~)guiDrawOverlay(f));
lbl(.905,'Left image');
names={}; data={}; chk=@(v) isnumeric(v)&&ismatrix(v)&&isequal(size(v),size(sMap));
for fn=fieldnames(results)', v=results.(fn{1}); if chk(v), names{end+1}=fn{1}; data{end+1}=v; end, end %#ok<AGROW>
% sMap-based reconstruction of every detection parameter (each segment filled
% with its node value) + the combined flow potential
vpx=compImg>0;
for pk=1:numel(s.flowParams)
    if isfield(s.flowParams,'role') && strcmp(s.flowParams(pk).role,'tip'), continue; end  % geometric, no map
    pnm=s.flowParams(pk).name; nv=getMetric(M,pnm); nv=nv(nodeIds);
    pimg=nan(size(sMap)); pimg(vpx)=nv(compImg(vpx));
    names{end+1}=['recon: ' pnm]; data{end+1}=pimg; %#ok<AGROW>
end
pimg=nan(size(sMap)); pimg(vpx)=H.phi(compImg(vpx));
names{end+1}='recon: flowPotential'; data{end+1}=pimg; %#ok<AGROW>
idx0=find(strcmp(names,'imgBFI'),1); if isempty(idx0), idx0=1; end
imgPop=uicontrol(p,'Style','popup','Units','normalized','Position',[.05 .878 .9 .026], ...
    'String',names,'Value',idx0,'Callback',@(src,~)guiLeftImage(f));
imgPop.UserData=struct('data',{data},'names',{names});

% ---- detection parameters (enable + weight) ---------------------------
yc=.845; lbl(yc,'Detection parameters  (on / weight)'); yc=yc-.028;
paramCtl=struct('chk',{},'wt',{}); np=numel(s.flowParams);
for k=1:np
    nm=s.flowParams(k).name;
    if isfield(s.flowParams,'label')&&~isempty(s.flowParams(k).label), nm=s.flowParams(k).label; end
    paramCtl(k).chk=uicontrol(p,'Style','checkbox','String',nm,'Units','normalized', ...
        'Position',[.05 yc .62 .022],'Value',s.flowParams(k).enabled,'BackgroundColor','w'); %#ok<AGROW>
    paramCtl(k).wt=uicontrol(p,'Style','edit','String',num2str(s.flowParams(k).weight), ...
        'Units','normalized','Position',[.70 yc .25 .022]); %#ok<AGROW>
    yc=yc-.0235;
end
lbl(yc,'Parench degree: artery / vein / parench / max'); yc=yc-.028;
degFld={'parenchArteryDeg','parenchVeinDeg','parenchParenchDeg','parenchMaxDeg'}; degCtl=gobjects(1,4);
for k=1:4
    degCtl(k)=uicontrol(p,'Style','edit','String',num2str(s.(degFld{k})),'Units','normalized', ...
        'Position',[.05+(k-1)*.235 yc .20 .026]); %#ok<AGROW>
end
yc=yc-.030;
lbl(yc,'Bridge gap px  (tip / wall) - splits in 2D FOV'); yc=yc-.027;
brFld={'bridgeTipRadius','bridgeWallRadius'}; brCtl=gobjects(1,2);
for k=1:2
    brCtl(k)=uicontrol(p,'Style','edit','String',num2str(s.(brFld{k})),'Units','normalized', ...
        'Position',[.05+(k-1)*.47 yc .42 .026]); %#ok<AGROW>
end
yc=yc-.032;
btn(yc,'Re-derive with these parameters',@(~,~)guiRederive(f)); yc=yc-.030;
wallsChk=uicontrol(p,'Style','checkbox','String','Show vessel walls','Units','normalized', ...
    'Position',[.05 yc .9 .026],'Value',1,'BackgroundColor','w','Callback',@(src,~)guiToggleWalls(f,src)); yc=yc-.028;
edgeChk=uicontrol(p,'Style','checkbox','String','Show all edges','Units','normalized', ...
    'Position',[.05 yc .9 .026],'BackgroundColor','w','Callback',@(src,~)guiToggleEdges(f,src)); yc=yc-.034;
lbl(yc,'Inspect / edit'); yc=yc-.030;
btn(yc,'Inspect branch (lineage)',@(~,~)guiSetMode(f,'branch')); yc=yc-.032;
btn(yc,'Inspect neighbours',@(~,~)guiSetMode(f,'neighbours')); yc=yc-.032;
btn(yc,'Set parent/daughter',@(~,~)guiSetMode(f,'reparent')); yc=yc-.032;
btn(yc,'Mark root/outlet',@(~,~)guiSetMode(f,'root')); yc=yc-.032;
btn(yc,'Inspect (no edit)',@(~,~)guiSetMode(f,'inspect')); yc=yc-.032;
btn(yc,'Reset to auto',@(~,~)guiResetAuto(f)); yc=yc-.030;
modeTxt=uicontrol(p,'Style','text','String','Mode: inspect','Units','normalized', ...
    'Position',[.05 yc .9 .026],'HorizontalAlignment','left','ForegroundColor',[.2 .2 .7],'BackgroundColor','w');
infoTxt=uicontrol(p,'Style','text','String','(hover a segment)','Units','normalized', ...
    'Position',[.03 .05 .94 .12],'HorizontalAlignment','left','BackgroundColor',[.96 .96 .98]);
uicontrol(p,'Style','pushbutton','String','Finish','Units','normalized', ...
    'Position',[.05 .008 .9 .038],'FontWeight','bold','Callback',@(~,~)uiresume(gcbf));

setappdata(f,'overlayPop',overlayPop); setappdata(f,'imgPop',imgPop);
setappdata(f,'infoTxt',infoTxt); setappdata(f,'edgeChk',edgeChk); setappdata(f,'modeTxt',modeTxt);
setappdata(f,'wallsChk',wallsChk);
setappdata(f,'branchFcn',@guiInspectBranch); setappdata(f,'neighFcn',@guiInspectNeighbours);
setappdata(f,'rederiveFcn',@guiRederive);   % test hooks (harmless)

st=struct('s',s,'results',results,'linkUp',linkUp,'linkDn',linkDn,'memb',H.memb,'adjC',{adjC}, ...
    'nodeIds',nodeIds,'comp',comp,'compImg',compImg,'cx',cx,'cy',cy,'autoConf',autoConf, ...
    'userSet',userSet,'nNodes',nNodes,'imsz',size(sMap),'phi',H.phi,'nodeType',H.nodeType, ...
    'typeCode',typeCodeFromH(H,nNodes),'bfiN',bfiN,'piN',piN,'footN',footN,'peakN',peakN, ...
    'diamN',diamN,'cMask',cMask,'isWallPix',(cMask==3|cMask==4),'showWalls',true, ...
    'paramCtl',{paramCtl},'axL',axL,'axR',axR,'hMapR',hMapR,'sel',[],'mode','inspect', ...
    'showEdges',false,'gen',[],'strahler',[],'treeID',[],'conf',[],'kUp',[],'kDn',[],'pEdges',pEdges, ...
    'degCtl',{degCtl},'brCtl',{brCtl},'linkUpAuto',linkUp,'linkDnAuto',linkDn,'pEdgesAuto',pEdges);
guidata(f,st);
set(f,'WindowButtonMotionFcn',@(o,~)guiHover(o), 'WindowButtonDownFcn',@(o,~)guiClick(o));

guiRefresh(f);                 % compute orders + first overlay
uiwait(f);

% ---- harvest edited state & write back --------------------------------
st=guidata(f); results=st.results; linkUp=st.linkUp; linkDn=st.linkDn; pEdges=st.pEdges; userSet=st.userSet;
H=results.hierarchy; autoConf=st.autoConf; memb=st.memb;
nodeIds=st.nodeIds; nNodes=st.nNodes; delete(f);
H2=rebuildHierarchy(H,linkUp,linkDn,pEdges,userSet,autoConf,nodeIds,nNodes,memb);
results=applyHierarchy(results,H2);
end

function [linkUp,linkDn,pEdges,adjC,autoConf]=linkStateFromH(H,comp,nNodes)
% Editable state from a hierarchy struct: linkUp (an artery's parent), linkDn
% (a vein's daughter), the parenchyma-mesh edge set pEdges (compact), and the
% compact undirected adjacency.
linkUp=zeros(nNodes,1); linkDn=zeros(nNodes,1);
u=find(H.linkUpId>0); linkUp(u)=comp(H.linkUpId(u));
d=find(H.linkDnId>0); linkDn(d)=comp(H.linkDnId(d));
pEdges=zeros(0,2);
if isfield(H,'pEdgesId') && ~isempty(H.pEdgesId)
    pEdges=[comp(H.pEdgesId(:,1)), comp(H.pEdgesId(:,2))];
end
adjC=cell(nNodes,1);
for e=1:size(H.adjacency,1)
    a=comp(H.adjacency(e,1)); b=comp(H.adjacency(e,2));
    if a>0&&b>0, adjC{a}(end+1)=b; adjC{b}(end+1)=a; end
end
autoConf=H.confidence;
end

function tc=typeCodeFromH(H,nNodes)
tc=4*ones(nNodes,1);
tc(H.nodeType=="Artery")=1; tc(H.nodeType=="Vein")=2; tc(H.nodeType=="Parench")=3;
end

function v=getNode(M,name,nodeIds)
if any(strcmp(name,M.Properties.VariableNames)), col=double(M.(name)); else, col=nan(height(M),1); end
v=col(nodeIds);
end

function guiReposition(f,tileAx,panel,WH)
figPix=getpixelposition(f); tp=get(tileAx,'Position');
tpx=[tp(1:2).*figPix(3:4) tp(3:4).*figPix(3:4)];
w=WH(1); h=min(WH(2),max(120,tpx(4)));
x=tpx(1)+max(0,(tpx(3)-w)/2); y=tpx(2)+max(0,(tpx(4)-h)/2);
set(panel,'Position',[x y w h]);
end

function guiSetMode(f,mode)
st=guidata(f); st.mode=mode; st.sel=[]; guidata(f,st);
set(getappdata(f,'modeTxt'),'String',['Mode: ' mode]);
guiDrawOverlay(f);   % restore the normal overlay (also clears any branch colourbar)
end

function guiToggleEdges(f,src)
st=guidata(f); st.showEdges=logical(get(src,'Value')); guidata(f,st); guiRedrawCurrent(f);
end

function guiToggleWalls(f,src)
st=guidata(f); st.showWalls=logical(get(src,'Value')); guidata(f,st); guiRedrawCurrent(f);
end

function guiRedrawCurrent(f)
st=guidata(f);
if strcmp(st.mode,'branch') && ~isempty(st.sel), guiInspectBranch(f,st.sel);
elseif strcmp(st.mode,'neighbours') && ~isempty(st.sel), guiInspectNeighbours(f,st.sel);
else, guiDrawOverlay(f); end
end

function guiRederive(f)
% Read the parameter controls, re-run the auto-derivation, and rebuild the
% editable state (this discards manual link edits - it is a fresh derive).
st=guidata(f); pc=st.paramCtl;
for k=1:numel(pc)
    st.s.flowParams(k).enabled=logical(get(pc(k).chk,'Value'));
    wv=str2double(get(pc(k).wt,'String')); if ~isfinite(wv)||wv<0, wv=0; end
    st.s.flowParams(k).weight=wv;
end
dg=st.degCtl; df={'parenchArteryDeg','parenchVeinDeg','parenchParenchDeg','parenchMaxDeg'};
for k=1:numel(dg)
    dv=round(str2double(get(dg(k),'String')));
    if isfinite(dv)&&dv>=0, st.s.(df{k})=dv; end
end
bg=st.brCtl; bf={'bridgeTipRadius','bridgeWallRadius'};
for k=1:numel(bg)
    bv=round(str2double(get(bg(k),'String')));
    if isfinite(bv)&&bv>=0, st.s.(bf{k})=bv; end
end
guiSetInfo(f,'Re-deriving with the current parameters...'); drawnow;
H=getVascularTree(st.results,st.s);
st.results=applyHierarchy(st.results,H);
[st.linkUp,st.linkDn,st.pEdges,st.adjC,st.autoConf]=linkStateFromH(H,st.comp,st.nNodes);
st.memb=H.memb; st.phi=H.phi; st.nodeType=H.nodeType; st.typeCode=typeCodeFromH(H,st.nNodes);
st.userSet=false(st.nNodes,1); st.linkUpAuto=st.linkUp; st.linkDnAuto=st.linkDn; st.pEdgesAuto=st.pEdges; st.sel=[];
guidata(f,st); guiRefresh(f);
guiSetInfo(f,sprintf('Re-derived: %d edges, %d roots, %d outlets.', ...
    size(H.edges,1), numel(H.roots), numel(H.outlets)));
end

function a=guiAlpha(st,ov)
% Transparency for an overlay: hide NaNs, and hide vessel walls/externals
% (cMask 3/4) when the 'Show vessel walls' box is off.
a=~isnan(ov);
if ~st.showWalls, a=a & ~st.isWallPix; end
end

function guiResetAuto(f)
st=guidata(f);
% restore the initial auto links captured at open time (undoes manual edits).
if isfield(st,'linkUpAuto'), st.linkUp=st.linkUpAuto; st.linkDn=st.linkDnAuto; st.pEdges=st.pEdgesAuto; end
st.userSet(:)=false; st.sel=[]; guidata(f,st); guiRefresh(f);
end

function guiRefresh(f)
st=guidata(f);
[kUp,kDn]=edgesFromLink(st.linkUp,st.linkDn,st.pEdges);
st.kUp=kUp; st.kDn=kDn;                       % cache the directed edge list
[treeID,gen,strahler]=orderForest(st.linkUp,st.linkDn,kUp,kDn,st.nNodes);
st.treeID=treeID; st.gen=gen; st.strahler=strahler;
conf=st.autoConf; conf(st.userSet)=1; st.conf=conf;
guidata(f,st);
guiDrawOverlay(f);
end

function guiDrawOverlay(f)
st=guidata(f);
delete(findall(f,'Type','ColorBar'));       % clear any branch-inspection colourbar
set(st.axR,'Color','w');                     % restore white background (branch mode uses black)
pop=getappdata(f,'overlayPop'); modes=get(pop,'String'); mode=modes{get(pop,'Value')};
switch mode
    case 'generation',   val=st.gen;      cmap=parula(256); cl=[0 max(1,max(st.gen))];
    case 'strahler',     val=st.strahler; cmap=parula(256); cl=[1 max(1,max(st.strahler))];
    case 'treeID',       val=st.treeID;   m=max(1,max(st.treeID)); cmap=hsv(m); cl=[1 m];
    case 'confidence',   val=st.conf;     cmap=parula(256); cl=[0 1];
    case 'flowPotential',val=st.phi;      cmap=parula(256); cl=[0 1];
    case 'type',         val=st.typeCode; cmap=[0.85 .2 .2;.2 .4 .85;.2 .7 .3;.6 .6 .6]; cl=[.5 4.5];
    otherwise,           val=st.gen;      cmap=parula(256); cl=[0 max(1,max(st.gen))];
end
ov=nan(st.imsz); v=st.compImg>0; ov(v)=val(st.compImg(v));
set(st.hMapR,'CData',ov,'AlphaData',guiAlpha(st,ov));
colormap(st.axR,cmap); try, clim(st.axR,cl); catch, end
title(st.axR,['Hierarchy: ' mode],'Interpreter','none');
delete(findobj(st.axR,'Tag','edge'));
if st.showEdges, guiDrawEdges(f); end
guiDrawSelection(f);
end

function guiDrawEdges(f)
st=guidata(f);
kUp=st.kUp; kDn=st.kDn;
if isempty(kUp), return; end
X=[st.cx(kUp)'; st.cx(kDn)'; nan(1,numel(kUp))];
Y=[st.cy(kUp)'; st.cy(kDn)'; nan(1,numel(kUp))];
hold(st.axR,'on');
plot(st.axR,X(:),Y(:),'-','Color',[0 0 0 .25],'LineWidth',.5,'Tag','edge');
hold(st.axR,'off');
end

function guiLeftImage(f)
st=guidata(f); pop=getappdata(f,'imgPop'); u=pop.UserData; sel=get(pop,'Value');
im=double(u.data{sel});
imagesc(st.axL,im); axis(st.axL,'image','off');
try, clim(st.axL,prctile(im(isfinite(im)),[1,99])); catch, end
title(st.axL,u.names{sel},'Interpreter','none');
end

function guiDrawSelection(f)
st=guidata(f);
delete(findobj(st.axR,'Tag','sel'));
if isempty(st.sel), return; end
h=st.sel; [par,dau]=linksOf(h,st.kUp,st.kDn);
hold(st.axR,'on');
for q=par(:)', plot(st.axR,[st.cx(h) st.cx(q)],[st.cy(h) st.cy(q)],'-','Color',[.9 .1 .1],'LineWidth',2,'Tag','sel'); end
for q=dau(:)', plot(st.axR,[st.cx(h) st.cx(q)],[st.cy(h) st.cy(q)],'-','Color',[.1 .1 .9],'LineWidth',2,'Tag','sel'); end
plot(st.axR,st.cx(h),st.cy(h),'o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','y','Tag','sel');
hold(st.axR,'off');
end

function guiHover(f)
st=guidata(f);
if any(strcmp(st.mode,{'branch','neighbours'})), return; end   % keep the inspection summary visible
ax=ancestor(hittest(f),'axes');
if isempty(ax)||~(ax==st.axL||ax==st.axR), return; end
cp=get(ax,'CurrentPoint'); x=round(cp(1,1)); y=round(cp(1,2));
if x<1||y<1||x>st.imsz(2)||y>st.imsz(1), return; end
h=st.compImg(y,x); if h==0, guiSetInfo(f,'(background)'); return; end
[par,dau]=linksOf(h,st.kUp,st.kDn);
pid=st.nodeIds(par); did=st.nodeIds(dau);
txt=sprintf(['seg %d  (%s)\nBFI %.0f  PI %.3f\nfoot %.4f  peak %.4f\ndiam %.1f  phi %.2f\n' ...
    'gen %d  strahler %d  tree %d\nconf %.2f%s\nparents: %s\ndaughters: %s'], ...
    st.nodeIds(h), char(st.nodeType(h)), st.bfiN(h), st.piN(h), st.footN(h), st.peakN(h), ...
    st.diamN(h), st.phi(h), st.gen(h), st.strahler(h), st.treeID(h), st.conf(h), ...
    tern(st.userSet(h),' (edited)',''), numstr(pid), numstr(did));
guiSetInfo(f,txt);
end

function guiClick(f)
st=guidata(f); ax=ancestor(hittest(f),'axes');
if isempty(ax)||~(ax==st.axL||ax==st.axR), return; end
cp=get(ax,'CurrentPoint'); x=round(cp(1,1)); y=round(cp(1,2));
if x<1||y<1||x>st.imsz(2)||y>st.imsz(1), return; end
h=st.compImg(y,x); if h==0, return; end
switch st.mode
    case 'reparent'
        if isempty(st.sel)
            st.sel=h; guidata(f,st); guiDrawSelection(f);
            guiSetInfo(f,sprintf(['seg %d selected.\nClick a touching neighbour to (re)connect:\n' ...
                'same-type vessels set a parent/daughter;\nanything with parenchyma toggles a mesh link.'], st.nodeIds(h)));
        else
            N=st.sel;
            if h==N, st.sel=[]; guidata(f,st); guiDrawSelection(f); return; end
            if ~ismember(h,st.adjC{N})
                guiSetInfo(f,sprintf('seg %d is not adjacent to seg %d - edit ignored.', st.nodeIds(h), st.nodeIds(N)));
                st.sel=[]; guidata(f,st); guiDrawSelection(f); return;
            end
            if st.phi(N)<=st.phi(h), u=N; v=h; else, u=h; v=N; end       % u upstream, v downstream
            snap=struct('linkUp',st.linkUp,'linkDn',st.linkDn,'pEdges',st.pEdges);
            if st.memb(u)==1 && st.memb(v)==1                            % both arteries
                st.linkUp(v)=u; msg=sprintf('seg %d parent set to seg %d.', st.nodeIds(v), st.nodeIds(u));
            elseif st.memb(u)==2 && st.memb(v)==2                        % both veins
                st.linkDn(u)=v; msg=sprintf('seg %d daughter set to seg %d.', st.nodeIds(u), st.nodeIds(v));
            elseif (st.memb(u)==1&&st.memb(v)==2)||(st.memb(u)==2&&st.memb(v)==1)  % artery<->vein: forbidden
                guiSetInfo(f,'Arteries cannot connect directly to veins - connect each to the parenchyma instead. Edit ignored.');
                st.sel=[]; guidata(f,st); guiDrawSelection(f); return;
            else                                                         % parenchyma involved: toggle a mesh edge
                ex=find(st.pEdges(:,1)==u & st.pEdges(:,2)==v,1);
                if ~isempty(ex), st.pEdges(ex,:)=[]; msg=sprintf('mesh link seg %d -> %d removed.', st.nodeIds(u), st.nodeIds(v));
                else, st.pEdges(end+1,:)=[u v]; msg=sprintf('mesh link seg %d -> %d added.', st.nodeIds(u), st.nodeIds(v)); end
            end
            [kUp,kDn]=edgesFromLink(st.linkUp,st.linkDn,st.pEdges);
            if ~isempty(kUp) && ~isdag(digraph(kUp,kDn,[],st.nNodes))
                st.linkUp=snap.linkUp; st.linkDn=snap.linkDn; st.pEdges=snap.pEdges;
                guiSetInfo(f,'Edit rejected: would create a cycle.');
            else
                st.userSet(u)=true; st.userSet(v)=true; guiSetInfo(f,msg);
            end
            st.sel=[]; guidata(f,st); guiRefresh(f);
        end
    case 'root'
        mc=st.memb(h);
        if mc==1,     st.linkUp(h)=0; lbl='root';
        elseif mc==2, st.linkDn(h)=0; lbl='outlet';
        else,         st.pEdges(st.pEdges(:,1)==h | st.pEdges(:,2)==h,:)=[]; lbl='detached'; end
        st.userSet(h)=true; st.sel=[]; guidata(f,st);
        guiSetInfo(f,sprintf('seg %d marked as %s.', st.nodeIds(h), lbl));
        guiRefresh(f);
    case 'branch'
        st.sel=h; guidata(f,st); guiInspectBranch(f,h);
    case 'neighbours'
        st.sel=h; guidata(f,st); guiInspectNeighbours(f,h);
    otherwise
        st.sel=h; guidata(f,st); guiDrawSelection(f); guiHover(f);
end
end

function guiSetInfo(f,txt), set(getappdata(f,'infoTxt'),'String',txt); end

function guiInspectBranch(f,clk)
% Highlight the full flow lineage through the clicked segment: all ancestors
% (upstream, toward the arterial root) and all descendants (downstream,
% toward the venous outlet), coloured red(upstream) -> white(clicked) ->
% blue(downstream) so the clicked segment sits at the middle of the scale and
% the largest connected artery / vein sit at the red / blue extremes.  If the
% clicked segment is terminal (itself the largest connected artery or vein) it
% takes the corresponding extreme colour instead of the middle.
st=guidata(f);
kUp=st.kUp; kDn=st.kDn;
[t,anc,des,artId,veinId]=branchColorValues(kUp,kDn,st.phi,clk,st.nNodes);

% Truecolor image: every segment that is NOT in the lineage (and the whole
% background) is painted BLACK; the lineage is red(artery)->white(clicked)->
% blue(vein).  Opaque, so it does not depend on transparency/axes colour.
cmap=redWhiteBlue(256);
tv=nan(st.imsz); v=st.compImg>0; tv(v)=t(st.compImg(v));      % per-pixel lineage value (NaN elsewhere)
if ~st.showWalls, tv(st.isWallPix)=NaN; end                  % hidden walls stay black too
lin=~isnan(tv); ix=round(tv*255)+1; ix(ix<1)=1; ix(ix>256)=256;
R=zeros(st.imsz); G=zeros(st.imsz); B=zeros(st.imsz);         % black canvas
R(lin)=cmap(ix(lin),1); G(lin)=cmap(ix(lin),2); B(lin)=cmap(ix(lin),3);
set(st.hMapR,'CData',cat(3,R,G,B),'AlphaData',1);
colormap(st.axR,cmap); set(st.axR,'Color','k');
title(st.axR,sprintf('Branch through seg %d   (red = upstream/artery, blue = downstream/vein)', ...
    st.nodeIds(clk)),'Interpreter','none');
delete(findall(f,'Type','ColorBar'));
cb=colorbar(st.axR); cb.Ticks=[0 .5 1];
hasUp=numel(anc)>1; hasDn=numel(des)>1;
if hasUp&&hasDn,   cb.TickLabels={'artery','clicked','vein'};
elseif hasDn,      cb.TickLabels={'','clicked','vein'};      % clicked is the arterial end
elseif hasUp,      cb.TickLabels={'artery','clicked',''};    % clicked is the venous end
else,              cb.TickLabels={'','clicked',''}; end

% branch edges + clicked marker
delete(findobj(st.axR,'Tag','edge')); delete(findobj(st.axR,'Tag','sel'));
S=unique([anc(:);des(:)]); inS=false(st.nNodes,1); inS(S)=true;
if ~isempty(kUp)
    em=inS(kUp)&inS(kDn);
    if any(em)
        X=[st.cx(kUp(em))';st.cx(kDn(em))';nan(1,nnz(em))];
        Y=[st.cy(kUp(em))';st.cy(kDn(em))';nan(1,nnz(em))];
        hold(st.axR,'on'); plot(st.axR,X(:),Y(:),'-','Color',[.55 .55 .55 .4],'LineWidth',.5,'Tag','edge'); hold(st.axR,'off');
    end
end
hold(st.axR,'on'); plot(st.axR,st.cx(clk),st.cy(clk),'o','MarkerSize',11, ...
    'MarkerEdgeColor','w','MarkerFaceColor','y','Tag','sel'); hold(st.axR,'off');

guiSetInfo(f,sprintf(['branch through seg %d\nancestors (upstream): %d\n' ...
    'descendants (downstream): %d\nlargest artery: seg %d\nlargest vein: seg %d'], ...
    st.nodeIds(clk), numel(anc)-1, numel(des)-1, st.nodeIds(artId), st.nodeIds(veinId)));
end

function [t,anc,des,artId,veinId]=branchColorValues(kUp,kDn,phi,clk,nNodes)
% Pure computation for guiInspectBranch (unit-testable).  The colour scale is
% anchored on the CLICKED node, which is ALWAYS white (t=0.5) regardless of its
% level in the tree.  Colouring then spreads out from it: ancestors (upstream,
% toward the arterial root) shade white->RED (t in [0,0.5], reddest at the
% largest connected artery), descendants (downstream, toward the venous outlet)
% shade white->BLUE (t in [0.5,1], bluest at the largest connected vein).  t is
% NaN outside the lineage.  artId/veinId are the compact ids of the arterial
% (min-phi ancestor) and venous (max-phi descendant) extremes.
if isempty(kUp)
    anc=clk; des=clk;
else
    des=bfsearch(digraph(kUp,kDn,[],nNodes),clk);   % clicked + descendants
    anc=bfsearch(digraph(kDn,kUp,[],nNodes),clk);   % clicked + ancestors (reversed graph)
end
anc=anc(:); des=des(:);
phiC=phi(clk); phiUp=min(phi(anc)); phiDn=max(phi(des));
t=nan(nNodes,1);
% ancestors: clicked(white,0.5) -> most-upstream artery(red,0), clamped so an
% ancestor never turns blue even if its phi is (non-monotonically) above phiC.
dU=max(phiC-phiUp,eps);
for n=anc', t(n)=min(max(0.5*(phi(n)-phiUp)/dU,0),0.5); end
% descendants: clicked(white,0.5) -> most-downstream vein(blue,1), clamped so a
% descendant never turns red.
dD=max(phiDn-phiC,eps);
for n=des', t(n)=min(max(0.5+0.5*(phi(n)-phiC)/dD,0.5),1); end
t(clk)=0.5;                                          % clicked node is always white
[~,ia]=min(phi(anc)); artId=anc(ia);
[~,idd]=max(phi(des)); veinId=des(idd);
end

function m=redWhiteBlue(n)
x=linspace(0,1,n)';
R=interp1([0 .5 1],[0.75 0.97 0.10],x);
G=interp1([0 .5 1],[0.10 0.97 0.10],x);
B=interp1([0 .5 1],[0.10 0.97 0.75],x);
m=[R G B];
end

function guiInspectNeighbours(f,clk)
% Colour the clicked segment's adjacent neighbours by how likely each is its
% PARENT (red, strongly up-gradient in the flow potential) vs its DAUGHTER
% (blue, strongly down-gradient).  Ambiguous neighbours (phi close to the
% clicked segment) get near-white, similar colours - so a coin-flip parent
% choice shows up as similar reds/whites among the top candidates.
st=guidata(f);
nb=unique(st.adjC{clk}(:));
if isempty(nb), guiSetInfo(f,sprintf('seg %d has no neighbours.',st.nodeIds(clk))); return; end
phiC=st.phi(clk);
d=st.phi(nb)-phiC;                    % <0 upstream (parent-like), >0 downstream (daughter-like)
ma=max([abs(d);eps]);
t=nan(st.nNodes,1);
t(nb)=0.5+0.5*d./ma;                  % 0 = red (parent) .. 1 = blue (daughter)
t(clk)=0.5;
ov=nan(st.imsz); v=st.compImg>0; ov(v)=t(st.compImg(v));
set(st.hMapR,'CData',ov,'AlphaData',guiAlpha(st,ov));
colormap(st.axR,redWhiteBlue(256)); clim(st.axR,[0 1]); set(st.axR,'Color','w');
title(st.axR,sprintf('Neighbours of seg %d   (red = likely parent, blue = likely daughter)', ...
    st.nodeIds(clk)),'Interpreter','none');
delete(findall(f,'Type','ColorBar'));
cb=colorbar(st.axR); cb.Ticks=[0 .5 1]; cb.TickLabels={'parent','ambiguous','daughter'};
delete(findobj(st.axR,'Tag','edge')); delete(findobj(st.axR,'Tag','sel'));
hold(st.axR,'on');
for q=nb', plot(st.axR,[st.cx(clk) st.cx(q)],[st.cy(clk) st.cy(q)],'-','Color',[.5 .5 .5 .5],'Tag','sel'); end
plot(st.axR,st.cx(clk),st.cy(clk),'o','MarkerSize',11,'MarkerEdgeColor','k','MarkerFaceColor','y','Tag','sel');
hold(st.axR,'off');
% info: neighbours sorted most-parent-like first, current parent(s) flagged
[parC,~]=linksOf(clk,st.kUp,st.kDn); curPar=st.nodeIds(parC);
[~,ord]=sort(d);
lines=cell(0);
for jj=ord(1:min(numel(ord),9))'
    nId=st.nodeIds(nb(jj)); tg='';
    if ismember(nId,curPar), tg=' [current parent]'; end
    lines{end+1}=sprintf('seg %d: dphi=%+.3f%s',nId,d(jj),tg); %#ok<AGROW>
end
guiSetInfo(f,sprintf('neighbours of seg %d (parent-like first):\n%s', ...
    st.nodeIds(clk), strjoin(lines,newline)));
end

%% ---- GUI hierarchy helpers -------------------------------------------
function [kUp,kDn]=edgesFromLink(linkUp,linkDn,pEdges)
% Directed flow edges (parent->child): vessel-tree links linkUp(n)->n and
% n->linkDn(n), plus the parenchyma-mesh edges pEdges.
up=find(linkUp>0); dn=find(linkDn>0);
kUp=[linkUp(up); dn]; kDn=[up; linkDn(dn)];
if nargin>=3 && ~isempty(pEdges), kUp=[kUp; pEdges(:,1)]; kDn=[kDn; pEdges(:,2)]; end
if ~isempty(kUp), E=unique([kUp kDn],'rows'); kUp=E(:,1); kDn=E(:,2); end
end

function [par,dau]=linksOf(h,kUp,kDn)
% Compact parents & daughters of node h from the directed edge list.
par=unique(kUp(kDn==h)); dau=unique(kDn(kUp==h));
end

function H2=rebuildHierarchy(H,linkUp,linkDn,pEdges,userSet,autoConf,nodeIds,nNodes,memb)
% Reassemble the hierarchy struct from the (possibly edited) links + mesh.
[kUp,kDn]=edgesFromLink(linkUp,linkDn,pEdges);
[treeID,gen,strahler]=orderForest(linkUp,linkDn,kUp,kDn,nNodes);
conf=autoConf; conf(userSet)=1;
parentList=cell(nNodes,1); daughterList=cell(nNodes,1);
for e=1:numel(kUp)
    daughterList{kUp(e)}(end+1)=nodeIds(kDn(e));
    parentList{kDn(e)}(end+1)  =nodeIds(kUp(e));
end
isRoot=(memb==1)&(linkUp==0); isOutlet=(memb==2)&(linkDn==0);
H2=H;
up=find(linkUp>0); dn=find(linkDn>0);
H2.linkUpId=zeros(nNodes,1); H2.linkUpId(up)=nodeIds(linkUp(up));
H2.linkDnId=zeros(nNodes,1); H2.linkDnId(dn)=nodeIds(linkDn(dn));
H2.pEdgesId=zeros(0,2); if ~isempty(pEdges), H2.pEdgesId=[nodeIds(pEdges(:,1)), nodeIds(pEdges(:,2))]; end
H2.edges=[nodeIds(kUp), nodeIds(kDn), ones(numel(kUp),1)];
H2.parentList=parentList; H2.daughterList=daughterList;
H2.treeID=treeID; H2.generation=gen; H2.strahler=strahler; H2.confidence=conf;
H2.isRoot=isRoot; H2.isOutlet=isOutlet; H2.roots=nodeIds(isRoot); H2.outlets=nodeIds(isOutlet);
H2.edited=any(userSet);
end

function s=numstr(v)
if isempty(v), s='(none)'; else, s=strtrim(sprintf('%d ',v(:))); end
end
function o=tern(c,a,b), if c, o=a; else, o=b; end, end
