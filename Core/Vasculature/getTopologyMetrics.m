%getTopologyMetrics - Vascular density and shape statistics from a segmented angiogram
%
%   Six published metrics of a vascular field, computed from what runSegmentation
%   already wrote.  It derives NO geometry of its own: the mask, the per-segment
%   lengths, calibres and chord-length ratios all arrive as inputs, and the two
%   labelling transients (sLines, nodes) are re-derived by the caller with
%   getSegmentationLabels exactly as runDynamicSegmentation re-derives them.  It never
%   opens a file and never touches the data cube.
%
%   THE DENOMINATOR IS THE ANALYSED AREA, NOT THE FRAME AND NOT A CONVEX HULL.
%   getPixelCategories ends with cMask = cMask.*mask, so cMask == 0 marks every pixel
%   the segmentation refused to look at - outside the drawn regions, outside the trust
%   mask, inside the edge trim.  Counting only cMask > 0 is the one denominator that
%   cannot be inflated by pixels nothing was measured on.  That follows REAVER (the
%   fraction of the analysed image) and deliberately NOT AngioTool, which normalises by
%   the convex hull of the vessels: the right denominator for an explant floating in a
%   dish, where the specimen's own extent is a measurement, and the wrong one for a
%   cranial window, where the window defines the field and a single vessel found or
%   missed at the rim moves the denominator.
%
%   WHAT IS A VESSEL PIXEL.  cMask > 3, i.e. innerEdge + lumen.  That is the same mask
%   runSegmentation's dMask uses to measure diameter, so the area fraction and the
%   calibre distribution are consistent by construction.  Category 3 (outerEdge) is the
%   band OUTSIDE the vessel boundary and is neither vessel nor tissue; category 2
%   (unsegmented) is vessel-like but unresolved, is NOT counted as vessel, and is
%   reported on its own as a trust number.
%
%   JUNCTIONS ARE CONNECTED COMPONENTS OF THE NODE MASK, NEVER NODE PIXELS.  A three-way
%   meeting is several 8-connected skeleton pixels, and counting them individually is the
%   >100 % over-count VesselVio's clique filter exists to fix (their figure: a mean
%   122.3 % error reduced to 2.6 %).  getSegmentationLabels already clusters adjacent
%   junction pixels and, under s.correctNodes, separates crossings from branchings; this
%   counts what it produced.
%
%   TWO RULES OF THE EXTRAVASCULAR DISTANCE THAT ARE PART OF THE DEFINITION
%     - the distance transform is computed ONCE on the whole frame, never per region.  A
%       region's tissue is allowed to find its nearest vessel outside the region;
%       cropping the transform to a region would invent a large distance at every region
%       edge.
%     - a tissue pixel whose distance to the frame border is smaller than its distance to
%       the nearest vessel is EXCLUDED, because its true nearest vessel may be outside
%       the field.  The retained fraction is reported as evdCoverage; below about 0.8 the
%       field is too sparsely vascularised for the distance statistics to be read.
%
%   AND TWO THINGS THE NUMBERS ARE NOT.  This is a 2-D projection of a 3-D vasculature.
%   A crossing at two depths is indistinguishable from a bifurcation, so junction density
%   is not a branching count; and projection can only shrink a distance, while the
%   capillary bed is not resolved at this pixel size, so the extravascular distance is
%   not an oxygen diffusion distance.  Both are within-study comparators.  The wording
%   that has to reach the reader is on the report page runTopologyAnalysis writes.
%
%   WHAT IS DELIBERATELY ABSENT.  Lacunarity and fractal dimension: measured across a
%   25-point settings sweep, lacunarity is MORE stable than vessel area fraction (1.28x
%   against 1.51x) but correlates with it at Pearson -0.959, so it is that measurement
%   sign-flipped and rescaled rather than a second one.  Branch order: it needs a rooted,
%   directed tree, which is setVascularTree's job.  Nothing here goes into sMetrics - a
%   density is a property of an area, not of a segment.
%
%   UNITS FOLLOW THE LIBRARY'S CONVENTION (Launcher_myograph).  s.pixelSize empty or 0
%   means an uncalibrated recording and everything is reported in pixels; a positive
%   value reports lengths in micrometres, areas in mm^2, densities in mm/mm^2 and mm^-2.
%   m.units records which, and the histogram edges and s.evdBeyond are read in whichever
%   unit is in force.
%
% Syntax:
%    m = getTopologyMetrics(cMask,sMetrics,sLines,nodes,regionsMask,s)
%    [m,maps] = getTopologyMetrics(cMask,sMetrics,sLines,nodes,regionsMask,s)
%
% Inputs:
%    cMask       - int32 five-level category mask (results.cMask, UN-merged: 0
%                  background, 1 parenchyma, 2 unsegmented, 3 outerEdge, 4 innerEdge,
%                  5 lumen).
%    sMetrics    - per-segment table (results.sMetrics).  Rows with category 5 carry
%                  length (px), CLR and diameter (px).
%    sLines      - labelled centre-line map from getSegmentationLabels.  Used ONLY to
%                  attribute each segment's length to a scope; the length itself comes
%                  from sMetrics.
%    nodes       - junction-node mask from getSegmentationLabels (post-correction).
%    regionsMask - integer region label image (results.regionsMask), 0 = excluded.  Empty
%                  or a single label means the whole analysed area only.
%    s           - parameter struct.  Reads:
%                    pixelSize    micrometres per pixel; [] or 0 -> everything in px.
%                    evdBeyond    the distance defining "tissue beyond", working unit.
%                                 Default 25 um, or 10 px uncalibrated.
%                    calibreEdges bin edges of the calibre histogram, working unit.
%                                 Default 0:5:100 um, or 0:2:40 px uncalibrated.
%                    evdEdges     bin edges of the distance histogram, working unit.
%                                 Default 0:2:60 um, or 0:1:24 px uncalibrated.
%                    keepMaps     logical; also return the distance-transform image.
%
% Outputs:
%    m    - struct with
%             .metrics  table, ONE ROW PER SCOPE.  scope 0 is the whole analysed area,
%                       scope k is region k of regionsMask.  Columns: areaAnalysed,
%                       areaFraction, unsegmentedFraction, lengthDensity,
%                       junctionDensity, endpointDensity, calibreMedian, calibreIQR,
%                       tortuosityMedian, evdMean, evdMedian, evdP95,
%                       evdFractionBeyond, evdCoverage, nSegments.
%             .calibre  struct(edges,counts); counts are vessel LENGTH per bin, one row
%                       per scope.
%             .evd      struct(edges,counts,beyond); counts are tissue AREA per bin, one
%                       row per scope.
%             .units    struct(pixelSize,lengthUnit,areaUnit,densityUnit,countUnit).
%    maps - struct with .evd (the distance transform, working unit, single) when
%           s.keepMaps is true, and empty otherwise.
%
% Example:
%    load(fName,'results');  load(getProductPath(fName,'s'),'settings');
%    [~,~,sLines,~,~,~,nodes] = getSegmentationLabels(results.cMask, ...
%                                   settings.runSegmentation.edgeSize,s);
%    s.pixelSize = 2.5;
%    m = getTopologyMetrics(results.cMask,results.sMetrics,sLines,nodes, ...
%                           results.regionsMask,s);
%    disp(m.metrics(1,:));
%
% See also: runTopologyAnalysis, getSegmentationLabels, getPixelCategories,
%           runSegmentation, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function [m,maps] = getTopologyMetrics(cMask,sMetrics,sLines,nodes,regionsMask,s)

if ~isfield(s,'pixelSize'), s.pixelSize=[];   end
if ~isfield(s,'keepMaps'),  s.keepMaps=false; end

% --- unit resolution: the myograph's convention, library-wide ---
calibrated = ~isempty(s.pixelSize) && s.pixelSize>0;
if calibrated
    px2len = s.pixelSize;                       % px -> um
    u = struct('pixelSize',s.pixelSize,'lengthUnit','um','areaUnit','mm^2', ...
               'densityUnit','mm/mm^2','countUnit','1/mm^2');
    area1  = (px2len/1000)^2;                   % one analysed pixel, in mm^2
    len1   = px2len/1000;                       % one pixel of centre line, in mm
else
    px2len = 1;                                 % px -> px
    u = struct('pixelSize',[],'lengthUnit','px','areaUnit','px^2', ...
               'densityUnit','px/px^2','countUnit','1/px^2');
    area1  = 1;
    len1   = 1;
end

% THE BIN EDGES AND THE THRESHOLD ARE IN THE WORKING UNIT, so an uncalibrated recording
% needs its own ladder: a micrometre ladder applied to pixel values silently puts every
% vessel and every tissue pixel in the first bin.  Both ladders are narrower than the
% research prototype's, and that is measured rather than guessed - on real fields all the
% calibre mass is below 60 um and all the distance mass is below 60 um, so the wide
% defaults drew a single bar with most of the axis empty.
if ~isfield(s,'evdBeyond') || isempty(s.evdBeyond)
    if calibrated, s.evdBeyond=25; else, s.evdBeyond=10; end            % um / px
end
if ~isfield(s,'calibreEdges') || isempty(s.calibreEdges)
    if calibrated, s.calibreEdges=0:5:100; else, s.calibreEdges=0:2:40; end
end
if ~isfield(s,'evdEdges') || isempty(s.evdEdges)
    if calibrated, s.evdEdges=0:2:60; else, s.evdEdges=0:1:24; end
end

% --- masks ---
cMask    = int32(cMask);
analysed = cMask>0;                               % the denominator, everywhere
vessel   = cMask>3;                               % innerEdge + lumen
tissue   = analysed & ~(cMask>=3);                % parenchyma + unsegmented
if isempty(regionsMask), regionsMask=double(analysed); end

% --- the distance transform is computed ONCE, on the whole frame ---
dEvd = bwdist(vessel)*px2len;                     % working length unit

% --- edge validity: a pixel closer to the frame than to any vessel may have its true
%     nearest vessel outside the field, so it is EXCLUDED and the loss is counted ---
[H,W]    = size(cMask);
[cc,rr]  = meshgrid(1:W,1:H);
dFrame   = min(min(rr-1,H-rr),min(cc-1,W-cc))*px2len;
evdValid = dEvd<=dFrame;

% --- per-segment quantities, straight out of sMetrics ---
% sMetrics, not dvsMetrics: the dynamic table is filtered by four quality gates that have
% nothing to do with anatomy - the straightest, best-fitting segments survive them - so a
% density computed from it would be biased by settings.
isSeg = sMetrics.category==5 & ~isnan(sMetrics.length);
segId = sMetrics.idx(isSeg);
segL  = sMetrics.length(isSeg);                   % px, geodesic along sLines
segD  = sMetrics.diameter(isSeg);                 % px
segT  = sMetrics.CLR(isSeg);                      % arc / chord, dimensionless

sLines = double(sLines);

% Endpoints follow REAVER and ignore the frame border, "because edge effects cause false
% positives".  The map is the same for every scope, so it is built once and masked below.
endPts = bwmorph(sLines>0,'endpoints');
endPts(1,:)=false; endPts(end,:)=false; endPts(:,1)=false; endPts(:,end)=false;

% --- scopes: 0 = whole analysed area, then every region label ---
% One region is not a region - it is the whole window setRegions leaves behind when
% nothing was drawn, and repeating the same numbers under two scope ids would put one
% measurement in the table twice.
regionIds = double(unique(nonzeros(regionsMask(:))))';
if numel(regionIds)<=1, scopes = 0; else, scopes = [0, regionIds]; end
nS = numel(scopes);

varNames = {'scope','areaAnalysed','areaFraction','unsegmentedFraction', ...
            'lengthDensity','junctionDensity','endpointDensity', ...
            'calibreMedian','calibreIQR','tortuosityMedian', ...
            'evdMean','evdMedian','evdP95','evdFractionBeyond','evdCoverage', ...
            'nSegments'};
T = array2table(nan(nS,numel(varNames)),'VariableNames',varNames);
calCounts = zeros(nS,numel(s.calibreEdges)-1);
evdCounts = zeros(nS,numel(s.evdEdges)-1);

for k=1:1:nS
    if scopes(k)==0, sc = analysed; else, sc = analysed & regionsMask==scopes(k); end
    nA = nnz(sc);
    T.scope(k)=scopes(k);
    T.areaAnalysed(k)=nA*area1;
    if nA==0, continue, end

    % --- 1. area fraction, and the trust companion ---
    T.areaFraction(k)        = nnz(vessel & sc)/nA;
    T.unsegmentedFraction(k) = nnz(cMask==2 & sc)/nA;

    % --- 2. length density: the sum of the LIBRARY'S OWN segment lengths, each weighted
    %        by the fraction of its centre line lying inside this scope ---
    % Our geodesic length is not REAVER's centre-line PIXEL COUNT: measured on a real
    % field ours is 4.6 % larger, because a pixel count scores a diagonal step 1 where the
    % true length is sqrt(2).  Ours is the more accurate of the two and it means a number
    % from this step sits about 5 % above a REAVER number on the same segmentation.
    wSeg = segmentWeight(sLines,segId,sc);
    Ltot = sum(segL(:).*wSeg(:),'omitnan')*len1;
    T.lengthDensity(k) = Ltot/(nA*area1);
    T.nSegments(k)     = sum(wSeg>0.5);

    % --- 3. junction and endpoint density ---
    T.junctionDensity(k) = countComponents(nodes & sc)/(nA*area1);
    T.endpointDensity(k) = nnz(endPts & sc)/(nA*area1);

    % --- 4. calibre, LENGTH-WEIGHTED ---
    % An unweighted median over segments over-weights short stubs, of which a skeleton has
    % many; weighting each segment by its own centre-line length is what makes "the median
    % calibre of this vasculature" mean the calibre of the median millimetre of vessel.
    w = segL(:).*wSeg(:);  d = segD(:)*px2len;
    ok = w>0 & isfinite(d);
    if any(ok)
        T.calibreMedian(k)= weightedPrctile(d(ok),w(ok),50);
        T.calibreIQR(k)   = weightedPrctile(d(ok),w(ok),75)-weightedPrctile(d(ok),w(ok),25);
        calCounts(k,:)    = weightedHist(d(ok),w(ok),s.calibreEdges)*len1;  % length units
    end

    % --- 5. tortuosity: sMetrics.CLR IS VesselVio's arc-chord ratio; no second definition
    t = segT(:); okt = w>0 & isfinite(t);
    if any(okt), T.tortuosityMedian(k)=weightedPrctile(t(okt),w(okt),50); end

    % --- 6. extravascular distance, over the VALID tissue of this scope ---
    tv = tissue & sc & evdValid;
    T.evdCoverage(k) = nnz(tv)/max(nnz(tissue & sc),1);
    if nnz(tv)>0
        dv = dEvd(tv);
        T.evdMean(k)  = mean(dv);
        T.evdMedian(k)= median(dv);
        T.evdP95(k)   = prctile(dv,95);
        T.evdFractionBeyond(k) = nnz(dv>s.evdBeyond)/numel(dv);
        evdCounts(k,:)= histcounts(dv,s.evdEdges)*area1;                   % area units
    end
end

m = struct();
m.metrics = T;
m.calibre = struct('edges',s.calibreEdges,'counts',calCounts);
m.evd     = struct('edges',s.evdEdges,'counts',evdCounts,'beyond',s.evdBeyond);
m.units   = u;

maps = struct();
if s.keepMaps, maps.evd = single(dEvd); end
end

% =====================================================================
function w = segmentWeight(sLines,segId,sc)
%segmentWeight  Fraction of each segment's centre-line pixels lying inside sc.
%   ONE accumarray pass over the label image, rather than a mask compare per segment.
w = zeros(numel(segId),1);
lab = sLines(sLines>0);
inS = sc(sLines>0);
if isempty(lab), return, end
n   = max(lab);
tot = accumarray(lab,1,[n 1]);
ins = accumarray(lab,double(inS),[n 1]);
ok  = segId>=1 & segId<=n;
w(ok) = ins(segId(ok))./max(tot(segId(ok)),1);
end

% =====================================================================
function n = countComponents(bw)
%countComponents  Number of 8-connected components in a logical mask.
if ~any(bw(:)), n=0; return, end
[~,n] = bwlabel(bw,8);
end

% =====================================================================
function v = weightedPrctile(x,w,p)
%weightedPrctile  Percentile p of x with weights w, by linear interpolation on the
%   cumulative weight with the midpoint convention.
x = x(:); w = w(:);
% A scope holding ONE segment is an ordinary case - a small drawn region often does - and
% interp1 needs two sample points, so the one-sample answer is written out rather than
% left to throw.
if isscalar(x), v = x; return, end
[x,i] = sort(x);  w = w(i);
c = cumsum(w) - 0.5*w;
c = c/sum(w);
v = interp1(c,x,p/100,'linear','extrap');
end

% =====================================================================
function h = weightedHist(x,w,edges)
%weightedHist  Sum of w in each bin of edges.
[~,~,b] = histcounts(x,edges);
h = zeros(1,numel(edges)-1);
ok = b>0;
if any(ok), h = accumarray(b(ok),w(ok),[numel(edges)-1 1])'; end
end
%------------- END OF CODE --------------
