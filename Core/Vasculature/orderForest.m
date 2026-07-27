%orderForest - Territory, flow generation & Horton-Strahler order of a vessel forest
%
%   [treeID,gen,strahler] = orderForest(linkUp,linkDn,kUp,kDn,nNodes) orders the
%   directed vascular forest:
%     treeID   - weakly-connected territory index per node;
%     gen      - global flow depth (0 at the arterial inlet, increasing to the
%                venous outlet) as the longest path on the flow DAG (kUp->kDn);
%     strahler - max Horton-Strahler order over the arterial (linkUp) and venous
%                (linkDn) single-pointer forests.
%   Pure helper shared by the vascular-tree derivation core (getVascularTree)
%   and the setVascularTree editor GUI (live re-ordering after user edits).
%
% Syntax:
%    [treeID,gen,strahler] = orderForest(linkUp,linkDn,kUp,kDn,nNodes)
%
% Inputs:
%    linkUp,linkDn - per-node single parent/daughter pointers (0 = none).
%    kUp,kDn       - directed edge endpoint lists (compact node indices).
%    nNodes        - number of nodes.
%
% Outputs:
%    treeID,gen,strahler - nNodes x 1 vectors (see above).
%
% Depends on: MATLAB graph objects (graph/digraph/conncomp/toposort/isdag).
% See also: getVascularTree, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026

%------------- BEGIN CODE --------------
function [treeID,gen,strahler] = orderForest(linkUp,linkDn,kUp,kDn,nNodes)
% treeID = weakly-connected territory; gen = global flow depth (0 at the
% arterial inlet, increasing to the venous outlet) via longest path on the
% flow DAG; strahler = max Horton-Strahler order over the arterial (linkUp)
% and venous (linkDn) single-pointer forests.
if isempty(kUp)
    treeID=(1:nNodes)'; gen=zeros(nNodes,1); strahler=ones(nNodes,1); return;
end
treeID=conncomp(graph(kUp,kDn,[],nNodes))';
Gd=digraph(kUp,kDn,[],nNodes); gen=zeros(nNodes,1);
if isdag(Gd)
    for c=toposort(Gd)
        ch=successors(Gd,c);
        if ~isempty(ch), gen(ch)=max(gen(ch),gen(c)+1); end
    end
else
    warning('setVascularTree:cyclic','Flow graph unexpectedly cyclic; generation=0.');
end
strahler=max(strahlerForest(linkUp,nNodes), strahlerForest(linkDn,nNodes));
end

function str=strahlerForest(ptr,nNodes)
% Horton-Strahler order on the forest defined by one pointer per node
% (child -> ptr(child)); ptr is strictly monotone in phi so it is acyclic.
str=ones(nNodes,1);
src=find(ptr>0);
if isempty(src), return; end
G=digraph(src,ptr(src),[],nNodes);
if ~isdag(G), return; end
for c=toposort(G)                      % a node is ordered before its target
    pr=predecessors(G,c);              % nodes pointing at c
    if ~isempty(pr)
        os=str(pr); m=max(os);
        if sum(os==m)>=2, str(c)=m+1; else, str(c)=m; end
    end
end
end
