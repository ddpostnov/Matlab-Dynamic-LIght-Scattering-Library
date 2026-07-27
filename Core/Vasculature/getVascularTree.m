%getVascularTree - Derive the vascular parent->daughter hierarchy (pure core)
%
%   H = getVascularTree(results,s) derives a flow-directed vascular hierarchy
%   over the segmented vessels and parenchyma described by RESULTS, with NO file
%   or GUI side effects.  It is the computational core extracted from
%   setVascularTree: the wrapper loads/saves the _r/_s files and runs the
%   interactive editor, while this function only computes and returns the
%   hierarchy struct H.
%
%   METHOD (staged, type-constrained connection search)
%     A global flow potential phi (from pulse arrival psTimeMin/psTimeMax and
%     pulsatility psPI) orders every node from the arterial inlet (low) to the
%     venous outlet (high).  Connections are found in stages:
%       (0) FOV BRIDGING - same-type vessels split by a crossing vessel in the
%           2D projection are re-linked when their tips (or side walls) come
%           within s.bridgeTipRadius (s.bridgeWallRadius) px.
%       (1) ARTERIAL tree - each artery links to the closest touching artery
%           just up-gradient in a side-specific potential (arteries only).
%       (2) VENOUS tree - each vein links to the closest touching vein just
%           down-gradient (veins only).
%       (3) PARENCHYMA mesh - only terminal vessels touch parenchyma; a
%           parenchyma links up to parenchArteryDeg arteries, parenchVeinDeg
%           veins and parenchParenchDeg parenchyma (parenchMaxDeg total).
%     Output: a flow-directed DAG with Horton-Strahler order and generation.
%     See setVascularTree for the full narrative and the live-editable weights.
%
% Syntax:
%    H = getVascularTree(results,s)
%
% Inputs:
%    results - RESULTS struct of a segmented + vessel-typed + pulsatility-
%              analysed *_BFI_r.mat.  Needs sMetrics (columns category, idx,
%              nearestVesIdx, optional type, and the flow metrics
%              psTimeMin/psTimeMax/psPI/BFI/diameter), sMap and cMask.
%    s       - parameter struct, FULLY POPULATED with the fields set by
%              setVascularTree's default block: flowParams (see
%              defaultFlowParams), connectivity, minBorder,
%              parenchArteryDeg/parenchVeinDeg/parenchParenchDeg/parenchMaxDeg,
%              bridgeTipRadius/bridgeWallRadius/bridgeTipBand.  Call via
%              setVascularTree, or replicate those defaults, before calling.
%
% Outputs:
%    H - hierarchy struct: seg2node, nodeIds, nodeType, memb, directed edges,
%        per-node phi/phiA/phiV, treeID, generation, strahler, confidence,
%        roots/outlets, candidate adjacency and params (consumed by
%        applyHierarchy and the setVascularTree GUI).
%
% Depends on: orderForest, getMetric (Core siblings); Image Processing Toolbox
%    (imdilate, bwdist, strel), Statistics Toolbox (tiedrank), MATLAB graph
%    objects (graph/digraph/conncomp/toposort/isdag).
% See also: setVascularTree, orderForest, getMetric, defaultFlowParams
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026

%------------- BEGIN CODE --------------
function H = getVascularTree(results,s)
% Pure (GUI-free) derivation of the vascular hierarchy from RESULTS.
% Returns a struct H with node list, directed edge list, per-node
% attributes and the undirected candidate adjacency (for the GUI).

M   = results.sMetrics;
N   = height(M);
sMap= results.sMap;
cMask=results.cMask;
cat = M.category;
idx = M.idx;
if iscell(idx), idx=cell2mat(idx); end
tp  = strings(N,1);
if any(strcmp('type',M.Properties.VariableNames))
    tp = string(M.type);
end

% ---- node model: collapse lumen(5)+wall(3); parenchyma(1) own node ----
seg2node = zeros(N,1);
isLum = (cat==5);
isWall= (cat==3);
isPar = (cat==1);
nvi   = M.nearestVesIdx;
seg2node(isLum) = idx(isLum);
seg2node(isWall)= nvi(isWall);          % wall -> its lumen node
seg2node(isPar) = idx(isPar);
seg2node(isnan(seg2node)) = 0;          % unsegmented / empty rows excluded

nodeIds = unique(seg2node(seg2node>0));
nodeIds = nodeIds(:);
nNodes  = numel(nodeIds);
comp    = zeros(N,1);  comp(nodeIds)=(1:nNodes)';   % segment idx -> compact node index

% ---- per-node type (node id == representative row) --------------------
nodeType = strings(nNodes,1);
for c=1:nNodes
    k=nodeIds(c);
    if cat(k)==1
        nodeType(c)="Parench";
    else
        t=tp(k);
        if strlength(t)==0, t="Uncertain"; end
        nodeType(c)=t;
    end
end
% ---- membership: A(rtery) / V(ein) / P(arenchyma) ----------------------
% A global monotone flow potential (arrival + pulsatility only) orders every
% node from arterial inlet (low) to venous outlet (high); it drives the
% parenchyma chains and splits uncertain-typed vessels.
phiG = combinePhi(M,s.flowParams,nodeIds,(1:nNodes)',0);
mm = phiG(isfinite(phiG));
if isempty(mm) || max(mm)-min(mm)==0
    warning('setVascularTree:noFlowSignal', ...
        ['No enabled arrival/pulsatility parameter with data (psTimeMax/' ...
         'psTimeMin/psPI) - flow direction is undefined.  Run on a ' ...
         'pulsatility (_c) file, or use s.useReference to propagate.']);
end
phiG(isnan(phiG)) = 0.5;
memb = zeros(nNodes,1);                         % 1=artery 2=vein 3=parenchyma
memb(nodeType=="Parench")=3; memb(nodeType=="Artery")=1; memb(nodeType=="Vein")=2;
unc = (nodeType=="Uncertain");  med = median(phiG,'omitnan');
memb(unc & phiG<=med)=1;  memb(unc & phiG>med)=2;
isA = memb==1;  isV = memb==2;  isP = memb==3;

% Side-specific potentials add the caliber cues (BFI/diameter) with the sign
% that is monotone on that side: on the arterial side caliber DECREASES
% downstream, on the venous side it INCREASES downstream.
phiA = combinePhi(M,s.flowParams,nodeIds,find(isA),-1);
phiV = combinePhi(M,s.flowParams,nodeIds,find(isV),+1);
phiA(isnan(phiA)) = phiG(isnan(phiA));
phiV(isnan(phiV)) = phiG(isnan(phiV));

% ---- undirected pixel adjacency between nodes -------------------------
lut = [0; seg2node];  nodeImg = lut(double(sMap)+1);
[I,J,B] = nodeAdjacency(nodeImg,s.connectivity);
keepB = B>=s.minBorder;  I=I(keepB); J=J(keepB); B=B(keepB);
cI=comp(I); cJ=comp(J);  ok = cI>0 & cJ>0 & cI~=cJ;
cI=cI(ok); cJ=cJ(ok); B=B(ok);
adjList = cell(nNodes,1);
for e=1:numel(cI), adjList{cI(e)}(end+1)=cJ(e); adjList{cJ(e)}(end+1)=cI(e); end

% ---- Stage 0: bridge vessels split by a crossing vessel in the 2D FOV --
% Two same-type vessels can be interrupted by a crossing vessel of the other
% type, so their pixels never touch.  Add candidate links between same-type
% vessels whose pixels come within s.bridgeTipRadius px of each other's TIP
% (the two geometric ends) or s.bridgeWallRadius px of each other's side wall
% (default 0 - a real split usually happens at the ends).  These feed the
% vessel-tree search + the GUI, but NOT the parenchyma mesh (cI/cJ untouched).
[pI,pJ,pGap] = bridgeAdjacency(nodeImg,comp,memb,nNodes, ...
    s.bridgeTipRadius,s.bridgeWallRadius,s.bridgeTipBand);
if ~isempty(pI)
    dup = ismember([pI pJ], sort([cI cJ],2), 'rows');   % drop already-touching pairs
    pI=pI(~dup); pJ=pJ(~dup); pGap=pGap(~dup);
    for e=1:numel(pI), adjList{pI(e)}(end+1)=pJ(e); adjList{pJ(e)}(end+1)=pI(e); end
end

% ---- node centroids + vessel tip anchors (for tip-distance scoring) ----
% A vessel has two ends; the parent (arteries) / daughter (veins) attaches at
% one tip, so a second connection - especially to parenchyma - is preferred at
% the OPPOSITE tip.  tipVessel / tipParench detection parameters weight this.
vpAll=nodeImg>0; ccAll=comp(nodeImg(vpAll)); [ryAll,rxAll]=find(vpAll);
cX=accumarray(ccAll,rxAll,[nNodes 1],@mean,NaN);
cY=accumarray(ccAll,ryAll,[nNodes 1],@mean,NaN);
[T1,T2,vlen]=nodeTips(nodeImg,comp,memb,nNodes,s.bridgeTipBand);
[tipVesW,tipParW]=tipParamWeights(s.flowParams);        % 0 = that tip parameter off

% ---- FOV-edge nodes (touch background) --------------------------------
atEdge = false(nNodes,1);
bg = imdilate(cMask==0, strel('square',3)) & nodeImg>0;
eNodes = comp(unique(nodeImg(bg)));  atEdge(eNodes(eNodes>0)) = true;

gapsG = abs(phiG(cI)-phiG(cJ));  phiScale = median(gapsG(gapsG>0));
if isempty(phiScale)||~isfinite(phiScale)||phiScale==0, phiScale=0.1; end

% ---- Stage 1: ARTERIAL tree (arteries connect only to arteries) -------
% Parent = the touching artery just up-gradient in phiA (higher BFI/diameter
% /PI, earlier arrival) and closest to it.  When the vessel-vessel tip
% parameter is on, the parent is also biased toward the tip OPPOSITE the
% vessel's strongest daughter candidate (an end-to-end join).
linkUp = zeros(nNodes,1);  linkDn = zeros(nNodes,1);  conf = nan(nNodes,1);
for a = find(isA)'
    nb = unique(adjList{a});  nb = nb(isA(nb));
    g  = phiA(a)-phiA(nb);
    up = nb(g>0);  gu = g(g>0);
    if isempty(up), continue; end
    if tipVesW>0
        dn=nb(g<0); gd=-g(g<0); avoid=0;              % strongest daughter candidate -> its tip
        if ~isempty(dn), [~,j]=min(gd); avoid=whichTipIdx(a,cX(dn(j)),cY(dn(j)),T1,T2); end
        [~,ix]=min(tipCost(a,up,cX,cY,T1,T2,vlen,avoid,gu,tipVesW));
    else
        [~,ix]=min(gu);
    end
    linkUp(a)=up(ix); conf(a)=confMargin(sort(gu));
end

% ---- Stage 2: VENOUS tree (veins connect only to veins) ---------------
% Daughter = the touching vein just down-gradient in phiV (larger vein) and
% closest to it; with the tip parameter on, opposite the strongest parent tip.
for v = find(isV)'
    nb = unique(adjList{v});  nb = nb(isV(nb));
    g  = phiV(nb)-phiV(v);
    dn = nb(g>0);  gd = g(g>0);
    if isempty(dn), continue; end
    if tipVesW>0
        up=nb(g<0); gu=-g(g<0); avoid=0;              % strongest parent candidate -> its tip
        if ~isempty(up), [~,j]=min(gu); avoid=whichTipIdx(v,cX(up(j)),cY(up(j)),T1,T2); end
        [~,ix]=min(tipCost(v,dn,cX,cY,T1,T2,vlen,avoid,gd,tipVesW));
    else
        [~,ix]=min(gd);
    end
    linkDn(v)=dn(ix); conf(v)=confMargin(sort(gd));
end

% ---- tree tip: which tip each vessel uses for its parent(A)/daughter(V) --
treeTip=zeros(nNodes,1);
for a=find(isA & linkUp>0)', treeTip(a)=whichTipIdx(a,cX(linkUp(a)),cY(linkUp(a)),T1,T2); end
for v=find(isV & linkDn>0)', treeTip(v)=whichTipIdx(v,cX(linkDn(v)),cY(linkDn(v)),T1,T2); end

% ---- terminal vessels (no smaller same-type vessel connected) ---------
hasArtChild = false(nNodes,1);  ac=find(isA & linkUp>0);  hasArtChild(linkUp(ac))=true;
termA = isA & ~hasArtChild;                     % leaf arterioles -> may feed parenchyma
hasVenTrib = false(nNodes,1);   vc=find(isV & linkDn>0);  hasVenTrib(linkDn(vc))=true;
termV = isV & ~hasVenTrib;                      % leaf venules -> may drain parenchyma

% ---- Stage 3: parenchyma mesh + terminal-vessel links -----------------
% Parenchyma may link to up to s.parenchVesselDeg vessels and
% s.parenchParenchDeg other parenchyma (total <= s.parenchMaxDeg); each
% terminal vessel takes one parenchyma.  Edges oriented low->high phiG,
% closest first.  Carried as pEdges alongside the vessel-tree links so a
% parenchyma can be a mesh node (many up/down neighbours).
% Arteries feed and veins drain the parenchyma; an artery never links to a
% vein directly, so every arterio-venous path runs through parenchyma.  The
% caps are per vessel type: <=paD arteries and <=pvD veins per parenchyma.
paD=s.parenchArteryDeg; pvD=s.parenchVeinDeg; ppD=s.parenchParenchDeg; mxD=s.parenchMaxDeg;
paCnt=zeros(nNodes,1); pvCnt=zeros(nNodes,1); ppCnt=zeros(nNodes,1); taOut=termA; tvIn=termV;
uu=cI; dd=cJ; sw=phiG(cI)>phiG(cJ); tt=uu(sw); uu(sw)=dd(sw); dd(sw)=tt;   % uu=lower phi
wgt=exp(-abs(phiG(cI)-phiG(cJ))./phiScale);
okT=(termA(uu)&isP(dd)) | (isP(uu)&isP(dd)) | (isP(uu)&termV(dd));
uu=uu(okT); dd=dd(okT); wgt=wgt(okT);
% tip parameter: prefer a parenchyma that touches the terminal vessel near its
% FREE tip (opposite the tip that carries its parent/daughter vessel link).
if tipParW>0
    for i=1:numel(uu)
        u=uu(i); d=dd(i);
        if termA(u)&&isP(d),      v=u; par=d;
        elseif isP(u)&&termV(d),  v=d; par=u;
        else,                     continue; end     % parench<->parench: no vessel tip
        ft=3-treeTip(v);                             % free tip index (3 if no tree link)
        if ft==1,     dft=hypot(T1(v,1)-cX(par),T1(v,2)-cY(par));
        elseif ft==2, dft=hypot(T2(v,1)-cX(par),T2(v,2)-cY(par));
        else,         dft=min(hypot(T1(v,1)-cX(par),T1(v,2)-cY(par)), ...
                              hypot(T2(v,1)-cX(par),T2(v,2)-cY(par))); end
        Lv=vlen(v); if ~isfinite(Lv)||Lv<1, Lv=1; end
        wgt(i)=wgt(i)*exp(-tipParW*dft/Lv);          % closer to the free tip -> ordered first
    end
end
[~,ord]=sort(wgt,'descend'); uu=uu(ord); dd=dd(ord);
pEdges=zeros(0,2);
for i=1:numel(uu)
    u=uu(i); d=dd(i); add=false;
    if termA(u) && isP(d)                            % artery -> parenchyma
        if taOut(u) && paCnt(d)<paD && paCnt(d)+pvCnt(d)+ppCnt(d)<mxD
            paCnt(d)=paCnt(d)+1; taOut(u)=false; add=true;
        end
    elseif isP(u) && termV(d)                        % parenchyma -> vein
        if tvIn(d) && pvCnt(u)<pvD && paCnt(u)+pvCnt(u)+ppCnt(u)<mxD
            pvCnt(u)=pvCnt(u)+1; tvIn(d)=false; add=true;
        end
    elseif isP(u) && isP(d)                          % parenchyma -> parenchyma
        if ppCnt(u)<ppD && ppCnt(d)<ppD && ...
           paCnt(u)+pvCnt(u)+ppCnt(u)<mxD && paCnt(d)+pvCnt(d)+ppCnt(d)<mxD
            ppCnt(u)=ppCnt(u)+1; ppCnt(d)=ppCnt(d)+1; add=true;
        end
    end
    if add, pEdges(end+1,:)=[u d]; end %#ok<AGROW>
end
conf(isP & (paCnt+pvCnt+ppCnt)>0) = 0.6;             % nominal confidence for meshed parenchyma

% ---- directed edges: vessel trees + parenchyma mesh -------------------
up = find(linkUp>0);  dn = find(linkDn>0);
kUp=[linkUp(up); dn; pEdges(:,1)];  kDn=[up; linkDn(dn); pEdges(:,2)];
if ~isempty(kUp), E=unique([kUp kDn],'rows'); kUp=E(:,1); kDn=E(:,2); end
parentList=cell(nNodes,1); daughterList=cell(nNodes,1);
for e=1:numel(kUp)
    daughterList{kUp(e)}(end+1)=nodeIds(kDn(e));
    parentList{kDn(e)}(end+1)  =nodeIds(kUp(e));
end

% ---- ordering, roots/outlets -------------------------------------------
[treeID,gen,strahler] = orderForest(linkUp,linkDn,kUp,kDn,nNodes);
isRoot   = isA & (linkUp==0);                   % arterial inlet
isOutlet = isV & (linkDn==0);                   % venous outlet
conf(isnan(conf) & (linkUp>0 | linkDn>0)) = 0.75;

% ---- assemble output ---------------------------------------------------
H = struct();
H.version=2; H.createdBy='setVascularTree'; H.date=datestr(now); H.imgSize=size(sMap);
H.seg2node=seg2node; H.nodeIds=nodeIds; H.nodeType=nodeType; H.memb=memb;
H.isVenous=isV; H.phi=phiG; H.phiA=phiA; H.phiV=phiV; H.atFOVedge=atEdge;
H.linkUpId=zeros(nNodes,1); H.linkUpId(up)=nodeIds(linkUp(up));
H.linkDnId=zeros(nNodes,1); H.linkDnId(dn)=nodeIds(linkDn(dn));
H.pEdgesId=zeros(0,2); if ~isempty(pEdges), H.pEdgesId=[nodeIds(pEdges(:,1)), nodeIds(pEdges(:,2))]; end
H.edges=[nodeIds(kUp), nodeIds(kDn), ones(numel(kUp),1)];
H.adjacency=[nodeIds(cI), nodeIds(cJ), double(B)];         % touching pairs (border length)
if ~isempty(pI)                                            % + FOV-bridge pairs (negative gap)
    H.adjacency=[H.adjacency; nodeIds(pI), nodeIds(pJ), -double(pGap)];
end
H.parentList=parentList; H.daughterList=daughterList;
H.treeID=treeID; H.generation=gen; H.strahler=strahler; H.confidence=conf;
H.isRoot=isRoot; H.isOutlet=isOutlet; H.roots=nodeIds(isRoot); H.outlets=nodeIds(isOutlet);
H.terminalA=nodeIds(termA); H.terminalV=nodeIds(termV);
H.params=struct('flowParams',s.flowParams,'connectivity',s.connectivity, ...
    'minBorder',s.minBorder,'phiScale',phiScale, ...
    'bridgeTipRadius',s.bridgeTipRadius,'bridgeWallRadius',s.bridgeWallRadius, ...
    'bridgeTipBand',s.bridgeTipBand);
end


%% =========================  HELPERS (private)  ====================== %%
function [sig,w] = addSig(sig,w,col,wt)
% Append a rank signal column and its weight only if the column has finite
% values (keeps the signal matrix and weight vector in lock-step).
if any(isfinite(col)), sig=[sig col(:)]; w=[w wt]; end
end

function phi = weightedRankMean(sig,w)
% Weighted mean of per-column rank-normalized signals -> [0,1], NaN safe.
if isempty(sig), phi=nan(0,1); return; end
n=size(sig,1);
R=nan(size(sig));
for j=1:size(sig,2)
    x=sig(:,j); ok=isfinite(x);
    if any(ok)
        rr=nan(n,1); rr(ok)=tiedrank(x(ok))./sum(ok);
        R(:,j)=rr;
    end
end
w=w(:)'; W=repmat(w,n,1); W(isnan(R))=0;
num=sum(W.*R,2,'omitnan'); den=sum(W,2);
phi=num./den; phi(den==0)=NaN;
end

function [I,J,B] = nodeAdjacency(nodeImg,conn)
% Undirected node adjacency (border-length weighted) from a labelled image.
P=[];
h1=nodeImg(:,1:end-1); h2=nodeImg(:,2:end);   P=[P; h1(:) h2(:)];
v1=nodeImg(1:end-1,:); v2=nodeImg(2:end,:);   P=[P; v1(:) v2(:)];
if conn==8
    d1=nodeImg(1:end-1,1:end-1); d2=nodeImg(2:end,2:end); P=[P; d1(:) d2(:)];
    e1=nodeImg(1:end-1,2:end);   e2=nodeImg(2:end,1:end-1); P=[P; e1(:) e2(:)];
end
keep=P(:,1)>0 & P(:,2)>0 & P(:,1)~=P(:,2);
P=P(keep,:);
P=sort(P,2);
[U,~,ic]=unique(P,'rows');
B=accumarray(ic,1);
I=U(:,1); J=U(:,2);
end

function [pI,pJ,pGap] = bridgeAdjacency(nodeImg,comp,memb,nNodes,tipR,wallR,tipBand)
% Candidate links between NON-touching SAME-TYPE vessel segments that a
% crossing vessel splits in the 2D FOV.  Each vessel's outer border is
% classified into TIP pixels (within tipBand of either end of its PCA major
% axis) and WALL pixels; another same-type vessel whose pixels fall within
% tipR px of a tip pixel, or wallR px of a wall pixel, becomes a candidate
% neighbour.  Returns compact node-index pairs pI<pJ and the min pixel gap.
pI=zeros(0,1); pJ=zeros(0,1); pGap=zeros(0,1);
if (tipR<=0 && wallR<=0) || nNodes==0, return; end
lab=double(nodeImg); [Hh,Ww]=size(lab);
bord=borderMask(lab);
vpix=find(lab>0); if isempty(vpix), return; end
lb=lab(vpix); [lbS,ord]=sort(lb); vpixS=vpix(ord);
starts=[1; find(diff(lbS)>0)+1]; stops=[starts(2:end)-1; numel(lbS)];
P=zeros(0,3);
for g=1:numel(starts)
    nid0=lbS(starts(g)); cn=comp(nid0);
    if cn<=0 || memb(cn)==3, continue; end               % vessels only (skip parenchyma)
    pix=vpixS(starts(g):stops(g));
    [py,px]=ind2sub([Hh Ww],pix);
    bsel=bord(pix); bpix=pix(bsel); bx=px(bsel); by=py(bsel);
    if isempty(bpix), continue; end
    isTip=true(numel(bpix),1);
    if numel(pix)>=2
        mu=[mean(px) mean(py)]; C=cov([px py]);
        [V,Dg]=eig(C); [~,mi]=max(diag(Dg)); ax=V(:,mi);
        sPr=([bx by]-mu)*ax; span=max(sPr)-min(sPr);      % project border pts on the major axis
        if span>0
            band=max(2,tipBand*span);
            isTip=(sPr-min(sPr)<=band)|(max(sPr)-sPr<=band);
        end
    end
    if tipR>0
        tp=bpix(isTip);   if ~isempty(tp), P=[P; scanNode(lab,comp,memb,cn,tp,tipR)]; end %#ok<AGROW>
    end
    if wallR>0
        wp=bpix(~isTip);  if ~isempty(wp), P=[P; scanNode(lab,comp,memb,cn,wp,wallR)]; end %#ok<AGROW>
    end
end
if isempty(P), return; end
[U,~,ic]=unique(P(:,1:2),'rows'); G=accumarray(ic,P(:,3),[],@min);
pI=U(:,1); pJ=U(:,2); pGap=G;
end

function P=scanNode(lab,comp,memb,cnSelf,srcPix,R)
% Same-type vessel pixels within R px of this node's source pixels (srcPix),
% as rows [loCompact hiCompact gap].  Uses a cropped distance transform from
% THIS node's sources so a node never shadows its own bridge to a neighbour.
P=zeros(0,3);
[H,W]=size(lab);
[sy,sx]=ind2sub([H W],srcPix);
r0=max(1,min(sy)-R); r1=min(H,max(sy)+R);
c0=max(1,min(sx)-R); c1=min(W,max(sx)+R);
sub=lab(r0:r1,c0:c1);
srcSub=false(size(sub));
srcSub(sub2ind(size(sub),sy-r0+1,sx-c0+1))=true;
D=bwdist(srcSub);
di=find(sub>0 & D<=R); if isempty(di), return; end
cn=comp(sub(di)); gp=D(di);
good=cn>0 & cn~=cnSelf; cn=cn(good); gp=gp(good);
if isempty(cn), return; end
same=memb(cn)==memb(cnSelf); cn=cn(same); gp=gp(same);   % artery<->artery, vein<->vein only
if isempty(cn), return; end
lo=min(cnSelf,cn); hi=max(cnSelf,cn);
P=[lo(:) hi(:) gp(:)];
end

function b=borderMask(lab)
% Node-border pixels: a labelled pixel on the image edge or 4-adjacent to a
% pixel with a different label (incl. background).
[H,W]=size(lab); isN=lab>0; b=false(H,W);
b(1,:)=b(1,:)|isN(1,:); b(end,:)=b(end,:)|isN(end,:);
b(:,1)=b(:,1)|isN(:,1); b(:,end)=b(:,end)|isN(:,end);
b(1:end-1,:)=b(1:end-1,:)|(isN(1:end-1,:)&lab(1:end-1,:)~=lab(2:end,:));
b(2:end,:)  =b(2:end,:)  |(isN(2:end,:)  &lab(2:end,:)  ~=lab(1:end-1,:));
b(:,1:end-1)=b(:,1:end-1)|(isN(:,1:end-1)&lab(:,1:end-1)~=lab(:,2:end));
b(:,2:end)  =b(:,2:end)  |(isN(:,2:end)  &lab(:,2:end)  ~=lab(:,1:end-1));
end

function [wV,wP] = tipParamWeights(fp)
% Weights of the two tip detection parameters (0 if absent/disabled): wV for
% the vessel-vessel scope, wP for the vessel-parenchyma scope.
wV=0; wP=0;
for k=1:numel(fp)
    if ~(isfield(fp,'role') && strcmp(fp(k).role,'tip')), continue; end
    if ~fp(k).enabled || fp(k).weight<=0, continue; end
    sc='parench'; if isfield(fp,'scope') && ~isempty(fp(k).scope), sc=fp(k).scope; end
    if strcmp(sc,'vessel'), wV=fp(k).weight; else, wP=fp(k).weight; end
end
end

function [T1,T2,vlen] = nodeTips(nodeImg,comp,memb,nNodes,tipBand)
% Per vessel node, the [x y] anchor of each end of its PCA major axis (mean of
% the tip-band border pixels there) and the axial length.  NaN for parenchyma.
T1=nan(nNodes,2); T2=nan(nNodes,2); vlen=nan(nNodes,1);
lab=double(nodeImg); [Hh,Ww]=size(lab);
bord=borderMask(lab);
vpix=find(lab>0); if isempty(vpix), return; end
lb=lab(vpix); [lbS,ord]=sort(lb); vpixS=vpix(ord);
starts=[1; find(diff(lbS)>0)+1]; stops=[starts(2:end)-1; numel(lbS)];
for g=1:numel(starts)
    nid0=lbS(starts(g)); cn=comp(nid0);
    if cn<=0 || memb(cn)==3, continue; end
    pix=vpixS(starts(g):stops(g));
    [py,px]=ind2sub([Hh Ww],pix);
    bsel=bord(pix); bx=px(bsel); by=py(bsel);
    if isempty(bx), bx=px; by=py; end
    mu=[mean(px) mean(py)];
    if numel(px)>=2, C=cov([px py]); [V,Dg]=eig(C); [~,mi]=max(diag(Dg)); ax=V(:,mi);
    else,            ax=[1;0]; end
    sB=([bx by]-mu)*ax; span=max(sB)-min(sB); vlen(cn)=max(span,1);
    if span<=0, T1(cn,:)=mu; T2(cn,:)=mu; continue; end
    band=max(2,tipBand*span); lo=sB-min(sB)<=band; hi=max(sB)-sB<=band;
    T1(cn,:)=[mean(bx(lo)) mean(by(lo))];
    T2(cn,:)=[mean(bx(hi)) mean(by(hi))];
end
end

function ti = whichTipIdx(v,px,py,T1,T2)
% Which tip (1 or 2) of vessel v the point (px,py) is nearest.
d1=(T1(v,1)-px)^2+(T1(v,2)-py)^2; d2=(T2(v,1)-px)^2+(T2(v,2)-py)^2;
if d1<=d2, ti=1; else, ti=2; end
end

function cost = tipCost(v,cand,cX,cY,T1,T2,vlen,avoidTip,g,w)
% Combined cost for choosing a vessel neighbour: normalized phi gap g plus w x
% the (length-normalized) distance from the candidate to the tip OPPOSITE
% avoidTip (or to the nearest tip when avoidTip==0).  Lower = better.
gN=g/max([g;eps]); n=numel(cand); tN=zeros(n,1);
Lv=vlen(v); if ~isfinite(Lv)||Lv<1, Lv=1; end
for i=1:n
    c=cand(i);
    d1=hypot(T1(v,1)-cX(c),T1(v,2)-cY(c));
    d2=hypot(T2(v,1)-cX(c),T2(v,2)-cY(c));
    if avoidTip==1,     dt=d2;          % parent occupies tip1 -> want opposite (tip2)
    elseif avoidTip==2, dt=d1;
    else,               dt=min(d1,d2); end
    tN(i)=dt/Lv;
end
cost=gN + w*tN;
end

function phi = combinePhi(M,flowParams,nodeIds,subset,caliberSign)
% Weighted rank-mean of the enabled parameters over a node SUBSET (others are
% NaN), each signed by its role so higher = more downstream.  caliberSign
% picks the caliber direction (-1 arterial, +1 venous, 0 = skip caliber).
nN=numel(nodeIds); sig=[]; w=[];
for k=1:numel(flowParams)
    fp=flowParams(k);
    if ~fp.enabled || fp.weight<=0, continue; end
    role='arrival';
    if isfield(fp,'role') && ~isempty(fp.role), role=fp.role; end
    switch role
        case 'pulsatility', sgn=-1;
        case 'caliber', if caliberSign==0, continue; end;  sgn=caliberSign;
        case 'tip',     continue;              % geometric, not part of the flow potential
        otherwise,      sgn=1;                 % arrival
    end
    col=getMetric(M,fp.name);  v=nan(nN,1);  v(subset)=sgn*col(nodeIds(subset));
    [sig,w]=addSig(sig,w,v,fp.weight);
end
phi=weightedRankMean(sig,w);
end

function c = confMargin(sortedGaps)
% Confidence from the margin between the closest and second-closest candidate.
if numel(sortedGaps)>=2
    g1=sortedGaps(1); g2=sortedGaps(2);
    if g2<=0, c=1; else, c=min(max(1-g1/g2,0),1); end
else
    c=0.75;
end
end
