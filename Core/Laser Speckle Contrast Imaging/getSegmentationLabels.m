%getSegmentationLabels - centerline -> indexed vessel/parenchyma label maps
%
%   [sMap,pMap,sLines,cMask,vsMap,dMask,nodes] = ...
%       getSegmentationLabels(cMask,edgeSize,s)
%   [...] = getSegmentationLabels(cMask,edgeSize,s,centerlines)
%
% DESCRIPTION
%   The indexed-map core lifted verbatim from runSegmentation: it turns the
%   five-level category mask into per-structure integer label maps.  It skeletonizes
%   the lumen, constrains the skeleton to the inner frame, isolates junction nodes,
%   optionally corrects branch nodes (s.correctNodes), geodesically grows each
%   centerline into its surrounding vessel to build the segment map, then tiles the
%   remaining parenchyma into a hex grid of labelled cells.  It is pure mask algebra
%   - it never touches the data cube; the wrapper's single metrics/traces loop fills
%   sMetrics and sData from sMap afterwards, unchanged.
%
%   Internal walls (category 4) are merged into external walls (category 3) up front,
%   as the downstream metrics/dynamic-segmentation code expects; the merged cMask is
%   returned so the caller works on the same categories the labels were built from.
%
% INPUT
%   cMask        int32 five-level category mask (results.cMask): 0 background,
%                1 parenchyma, 2 unsegmented, 3 external wall, 4 internal wall,
%                5 lumen.
%   edgeSize     border width (px) to exclude when constraining the skeleton
%                (settings.runSegmentation.edgeSize, from getPixelCategories).
%   s            parameter struct.  Reads: correctNodes, simR, difR (node
%                correction), sMinL (min segment length), prchNSize (parenchyma
%                cell size), parforSegmentationLabels (logical, default true - run
%                the per-label geodesic loop in parallel; a WORKER BOUND on the
%                parfor, not a branch, so false runs the identical loop body
%                serially in the client and starts no pool).
%   centerlines  OPTIONAL logical skeleton to use instead of the default
%                bwskel(cMask==5) - e.g. from a manual centerline editor.  Pass []
%                or omit for the automatic skeleton.
%
% OUTPUT
%   sMap     indexed label map (int32): vessels get odd ids, external walls the
%            following even id, parenchyma cells the trailing block of ids.
%   pMap     parenchyma-only label map (int32), the hex-grid cells (results.pMap).
%   sLines   labelled, renumbered centerline map (int32) - the metrics loop and the
%            dynamic-segmentation loop both consume it.
%   cMask    the input categories with internal walls (4) merged into external
%            walls (3) - what every downstream consumer uses.
%   vsMap    vessel-only label map (int32) before wall/parenchyma ids were folded
%            into sMap; needed by the dynamic-segmentation overlap tests.
%   dMask    logical wall+lumen mask (original cMask>3), for diameter measurement.
%   nodes    post-correction junction-node mask, consumed by dynamic segmentation.
%
% DEPENDS ON
%   MATLAB's Image Processing Toolbox (bwskel, bwlabel, bwdistgeodesic, bwdist,
%   bwareaopen, imdilate, conv2).  Uses parfor over the segment labels (each page
%   independent, integer-valued geodesic distances - order-independent, bit-exact).
%
% See also: runSegmentation, getPixelCategories, runDynamicSegmentation
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

function [sMap,pMap,sLines,cMask,vsMap,dMask,nodes] = getSegmentationLabels(cMask,edgeSize,s,centerlines)

% Temporary fix - diameter is estimated including inner walls, rest of
% parameters estimated for lumen and merged walls
dMask = cMask > 3;
dtLumen = bwdist(dMask ~= 1);
cMaskOrig = cMask;              % pre-merge categories, for the parenchyma-seed guards
cMask(cMask == 4) = 3;

% Centerlines: default to the lumen skeleton, or use a supplied override.
if nargin < 4 || isempty(centerlines)
    sLines = bwskel(cMask == 5);
else
    sLines = centerlines;
end

% Constrain skeleton to inner region boundaries
validRegionMask = zeros(size(sLines));
validRegionMask(edgeSize*2+1:end-edgeSize*2, edgeSize*2+1:end-edgeSize*2) = ...
    sLines(edgeSize*2+1:end-edgeSize*2, edgeSize*2+1:end-edgeSize*2);
sLines = (validRegionMask == 1);

% Identify and isolate nodes using convolutional adjacency matrix
adjacencyMatrix = [1, 1, 1; 1, 0, 1; 1, 1, 1];
nodes = logical((conv2(single(sLines), adjacencyMatrix, 'same') > 2) .* sLines);
sLinesIni = sLines;
sLines = sLines - nodes;

% Label isolated segments
sLines = int32(bwlabel(sLines));

if s.correctNodes
    numLabels = max(sLines(:));

    % Compute segment lumen distances

    segmentMetrics = zeros(numLabels, 2);
    for labelIndex = 1:numLabels
        segmentMask = (sLines == labelIndex);
        segmentMetrics(labelIndex, 1) = sum(segmentMask(:));
        segmentMetrics(labelIndex, 2) = mean(dtLumen(segmentMask));
    end

    % Define topological criteria
    similarityTolerance = s.simR;
    differenceRatio = s.difR;

    % --- Implement Connected Components Node Clustering ---
    labeledNodes = bwlabel(nodes,8);
    numNodeClusters = max(labeledNodes(:));

    for clusterIndex = 1:numNodeClusters
        % Isolate the current contiguous node cluster
        currentClusterMask = (labeledNodes == clusterIndex);

        % Dilate the cluster by 1 pixel to detect all adjacent segment labels
        dilatedCluster = imdilate(currentClusterMask, strel('square', 3));

        % Extract non-zero labels from sLines that intersect the dilated perimeter
        connectingPixels = sLines(dilatedCluster & (sLines > 0));
        connectedLabels = nonzeros(unique(connectingPixels));

        if numel(connectedLabels) >= 3
            % Retrieve mean radii estimates for connecting segments
            segmentRadii = zeros(numel(connectedLabels), 2);
            for branchIndex = 1:numel(connectedLabels)
                labelIdx = connectedLabels(branchIndex);
                segmentRadii(branchIndex, 1) = double(labelIdx);
                segmentRadii(branchIndex, 2) = segmentMetrics(labelIdx, 2);
            end

            % Sort segments by diameter descendantly
            [~, sortIdx] = sort(segmentRadii(:, 2), 'descend');
            sortedSegments = segmentRadii(sortIdx, :);
            r1 = sortedSegments(1, 2);
            r2 = sortedSegments(2, 2);
            r3 = sortedSegments(3, 2);

            % Evaluate continuity and bifurcation conditions
            isSimilar = r1 < (similarityTolerance * r2 + r2);
            isBifurcation = r2 >= (differenceRatio * r3 + r3);

            if (isSimilar && isBifurcation) || (r1 < (similarityTolerance/2 * r2 + r2)) || (r1>10 && abs(r2 - r1) < 1) || (r2 >= (differenceRatio*2 * r3 + r3))
                % Define truncation limit as average radius of the parent vessel
                truncationDistance = ceil((r1 + r2)/2);

                % Apply geodesic truncation to daughter branches
                for smallBranchIndex = 3:size(sortedSegments, 1)
                    smallLabel = sortedSegments(smallBranchIndex, 1);
                    branchMask = (sLines == smallLabel);

                    % Seed the geodesic transform using the entire junction cluster
                    combinedMask = branchMask | currentClusterMask;

                    % Calculate constrained path distance originating from the cluster edge
                    geoDist = bwdistgeodesic(combinedMask, currentClusterMask, 'quasi-euclidean');

                    % Eliminate intraluminal pixels of the branching vessel
                    pixelsToRemove = (geoDist > 0) & (geoDist <= truncationDistance);
                    sLinesIni(pixelsToRemove) = 0;
                end
                % Remove the modified node cluster from the mask
                nodes(currentClusterMask) = 0;
            end
        end
    end
    sLines = sLinesIni - nodes;
    sLines=bwareaopen(sLines,s.sMinL);


    % Final topological cleanup and reassignment
    sLines = int32(bwlabel(sLines > 0));
end

labels=nonzeros(unique(sLines));
distStack=inf([size(cMask) numel(labels)],'single');
% Parallelism is optional: a BOUND on the parfor (Inf workers or 0), never a branch -
% parfor(...,0) runs the identical body serially in the client and starts no pool.
% Default true when the field is absent, so every existing caller is unaffected.
nwLbl=Inf;
if isfield(s,'parforSegmentationLabels') && ~isempty(s.parforSegmentationLabels) ...
        && ~s.parforSegmentationLabels
    nwLbl=0;
end
parfor (k = 1:numel(labels), nwLbl)   % dtLumen broadcast below is intentional (full-image DT)
    distStack(:,:,k) = bwdistgeodesic(cMask>2, sLines == labels(k), 'quasi-euclidean')-mean(dtLumen(sLines == labels(k))); %#ok<PFBNS>
end
[~,tmp] = min(distStack,[],3);
vsMap     = zeros(size(cMask),'like',sLines);
vsMap(cMask>2) = labels(tmp(cMask>2)); %labeled segments
vsMap(vsMap>0)=(vsMap(vsMap>0)-1).*2+1;
sLines(sLines>0)=(sLines(sLines>0)-1).*2+1;

sMap=vsMap;
sMap(cMask==3)=sMap(cMask==3)+1;
idxs=int32(bwlabel(cMask==2))+max(sMap(:));
sMap(cMask==2)=idxs(cMask==2);


[M,N]       = size(cMaskOrig);                      % image dimensions
step        = double(s.prchNSize);                  % nominal cell width
R           = step/2;                               % half-offset
rows        = 1 : R*sqrt(3) : M;                    % √3·R vertical pitch
cols        = 1 : step         : N;
[C,Rr]      = meshgrid(cols,rows);                  % full grid of centres
C(2:2:end,:)= C(2:2:end,:) + R;                     % shift odd rows by R
C           = round(C);  Rr = round(Rr);            % integer coordinates
inFrame     =  Rr>=1 & Rr<=M & C>=1 & C<=N;         % guard ①
idxAll      = sub2ind([M N], Rr(inFrame), C(inFrame));
idxSeeds    = idxAll(cMaskOrig(idxAll) >= 0);       % guard ②: inside mask
seed        = false(M,N);  seed(idxSeeds) = true;   % binary seed image
[~,lbl]     = bwdist(seed,'euclidean');             % nearest-seed label
valid       = cMaskOrig>=1;     % area to overwrite
lbl(~valid) = 0;
nz                  = lbl>0;
[~,~,lbl(nz)]       = unique(lbl(nz));              % consecutive IDs
lbl                 = int32(lbl) + max(sMap(:));   % avoid clashes
pMap=lbl;
sMap(valid & cMask==1)         = lbl(valid & cMask==1);                   % updated label map
end
