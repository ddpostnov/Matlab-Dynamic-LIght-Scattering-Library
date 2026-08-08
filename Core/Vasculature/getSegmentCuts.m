%getSegmentCuts - Where to cut across one segment, and in which direction.
%
%   THE STEP THAT TURNS A SEGMENT INTO A MEASUREMENT.  A wall estimator takes ONE
%   cut; a segment is a CENTRE LINE.  Given the labelled centre line
%   getSegmentationLabels produced, this places perpendicular cuts along it and says
%   how long each one has to be - and nothing else.  It reads no data cube, no
%   intensity and no polarity, so it is the same geometry whatever is measured on it.
%
%   A SEGMENT LABEL IS A TREE, NOT A PATH, AND IGNORING THAT COSTS A FACTOR OF FIVE.
%   Walking the skeleton pixels of one label in bwdistgeodesic order and taking the
%   tangent from their local principal axis was measured on the reference
%   recording's largest artery: the local radius alternated between 12 and 6 pixels
%   every three steps of arc because the ordering was diving into each short spur and
%   back out, the perpendicular came out 36 degrees off, and the wall amplitude 4.8
%   times too small.  A short twig on a bwskel output is not an anomaly; it is what
%   the skeleton of a real vessel looks like.
%
%   AND A CUBIC IS NOT THE ANSWER EITHER, WHICH IS THE MORE USEFUL FINDING.
%   runDynamicSegmentation fits a third-degree polynomial through the label's pixels.
%   On the same artery - 439 skeleton pixels spanning 424 rows - it scores R-squared
%   0.971 and still misses the vessel by a median of 5.9 px: only 17 per cent of the
%   fitted curve comes within three pixels of the skeleton and the longest usable run
%   is 21 px of 424.  A cubic describes a short segment and cannot describe a long
%   one, and R-squared does not notice.  Across the reference field that model left
%   171 of 1 422 segments with any cuts at all, against 1 136 for this one.
%
%   SO THE CENTRE LINE IS THE LONGEST PATH THROUGH THE LABEL, SMOOTHED LOCALLY.  That
%   is the library's OTHER idea about centre lines: runSegmentation measures
%   sMetrics.length as the geodesic distance between the two most distant endpoints,
%   which is exactly a statement that the trunk of the tree IS the segment.  Two
%   geodesic transforms find that pair (farthest from anywhere, then farthest from
%   there - on a tree that pair is the diameter), a steepest descent traces the path
%   between them, and a moving average whose window is one radius gives the position
%   and the tangent.  A spur is off the path and is never visited; a bend of any
%   radius is followed; nothing is fitted globally, so a 424-pixel artery and a
%   20-pixel arteriole are handled by the same code with the same settings.
%
%   THE ENDS OF A SEGMENT ARE DROPPED, BECAUSE THAT IS WHERE THE BRANCHES ARE.
%   getSegmentationLabels already cuts the centre line at the junction nodes, so a
%   segment does not CONTAIN a branch point - but its two ends sit against one, and a
%   perpendicular within a span of the node runs through the vessel that branches off
%   rather than across the one being measured.  The margin is stated in radii because
%   that is the scale on which a neighbour intrudes.  So nothing happens at branch
%   points, because no cut is placed within reach of one, and a segment too short to
%   leave any cuts after the margin reports nothing rather than a guess.
%
%   THE SPAN IS COMPUTED ON THE ARRAY THE CUTS WILL BE MEASURED ON, AND THAT IS THE
%   WHOLE OF A TRAP THAT ONCE COST A THIRD OF AN AMPLITUDE.  The pad beyond three
%   radii is a physical margin of tissue, so a caller working in micrometres converts
%   it ONCE, before this call, and hands over pixels of the array it is about to
%   measure.  Computing the span on one grid and the measurement on another shortens
%   the baseline on one side, which moves the half-way level a wall is found at -
%   measured at 0.0559 px against 0.0830 on a 39 um vessel.  There is one grid here
%   by construction, which is how that stays fixed.
%
%   A CUT THAT WOULD LEAVE THE FRAME IS DROPPED, NEVER CLIPPED.  A clipped cut finds
%   its wall against a shorter baseline on one side, which is a bias on the diameter
%   and therefore indistinguishable from wall motion.
%
% Syntax:
%    [cuts,q] = getSegmentCuts(sLines,id,R,s)
%
% Inputs:
%    sLines - int32 labelled centre-line map from getSegmentationLabels.
%    id     - the label of the segment to cut.
%    R      - [rows columns] the local vessel radius at every pixel, bwdist(~dMask),
%             on the SAME grid as sLines.
%    s      - parameter struct.  Reads:
%             .interpF         points per pixel along the resampled path (4).
%             .smoothRadii     half-window of the tangent smoother, in radii (1).
%             .smoothMin       ...floored here, px (3).
%             .cutSpacing      arc length between cuts, px.  0 takes every point (3).
%             .endMarginRadii  no cut within this many local radii of either end (1.5).
%             .spanRadii       the cut's half-length, in radii (3)...
%             .spanPadPx       ...plus this many pixels of tissue.  ALREADY CONVERTED
%                              to pixels of this array by the caller.
%             .minCuts         a segment offered fewer cuts than this reports nothing (3).
%             .maxCLR          chord-length-ratio gate.  Inf lets every shape through,
%                              which is the default: the detection confidence is what
%                              refuses a segment, and a geometry gate here would turn
%                              the coverage number into a statement about shapes.
%
% Outputs:
%    cuts - [1 n] struct array, EMPTY when the segment cannot carry a measurement.
%           .center [row column] on the smoothed path, pixels of sLines
%           .normal [drow dcolumn] unit, across the vessel
%           .radius the local radius there, px
%           .span   half-length of the cut, px
%           .arc    arc length from one end of the path, px
%    q    - why, and how much of the label was used:
%           .pathPx .clr .nSkel .nPoint .spurFrac .why
%
% Example:
%    [~,~,sLines,cMask,~,dMask] = getSegmentationLabels(results.cMask,edgeSize,s);
%    s.spanPadPx = 30/pixelSize;                 % a 30 um margin, converted ONCE
%    [cuts,q] = getSegmentCuts(sLines,1715,bwdist(~dMask),s);
%    fprintf('%d cuts on %.0f px of path, %.0f%% of the label is spur\n', ...
%        numel(cuts), q.pathPx, 100*q.spurFrac);
%
% Dependencies: Image Processing Toolbox (bwdistgeodesic).
% See also: getWallMotion, runMotionEnhancement, getSegmentationLabels,
%           runSegmentation, meKymograph
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function [cuts,q] = getSegmentCuts(sLines,id,R,s)

cuts = struct('center',{},'normal',{},'radius',{},'span',{},'arc',{});
q    = struct('pathPx',NaN,'clr',NaN,'nSkel',0,'nPoint',0,'spurFrac',NaN,'why','');

[sy,sx] = find(sLines==id);
q.nSkel = numel(sy);
if q.nSkel<4, q.why='too few centre-line pixels'; return; end

% ---- the longest path through the label, on its own bounding box -------------
% Cropping is not a matter of taste: two geodesic transforms per segment over a
% 1300x800 frame is a minute of a sweep over a field, and over a bounding box it is a
% second.  Geodesic distance is a property of the mask and every pixel of the segment
% is inside its own box, so the distances are unchanged.
r0 = min(sy)-1;  c0 = min(sx)-1;
bw = false(max(sy)-r0+1, max(sx)-c0+1);
bw(sub2ind(size(bw), sy-r0, sx-c0)) = true;

[py,px] = longestPath(bw);
if numel(py)<4, q.why='no path through the label'; return; end
py = py+r0;  px = px+c0;
q.nPoint   = numel(py);
q.spurFrac = 1 - numel(py)/q.nSkel;

% ---- resample on arc length, then smooth -------------------------------------
step     = [0; cumsum(hypot(diff(px),diff(py)))];
q.pathPx = step(end);
q.clr    = step(end)/max(hypot(px(end)-px(1),py(end)-py(1)),eps);
if q.clr>s.maxCLR, q.why='more curved than the shape gate allows'; return; end

[nr,nc] = size(sLines);
lin  = sub2ind([nr nc],py,px);
rad0 = max(double(R(lin)),0.5);

a = (0:1/s.interpF:step(end)).';
if numel(a)<5, q.why='path too short to resample'; return; end
yy = interp1(step,double(py),a,'linear');
xx = interp1(step,double(px),a,'linear');
rr = interp1(step,rad0,      a,'linear');

% The smoother's window is the scale over which a vessel is straight: about its own
% calibre.  A wider one cuts corners on a bend; a narrower one follows the skeleton's
% one-pixel staircase and puts a 45-degree wobble on every tangent.
w  = max(s.smoothMin, s.smoothRadii*median(rr));
nW = 2*round(w*s.interpF)+1;
ys = movmean(yy,nW,'Endpoints','shrink');
xs = movmean(xx,nW,'Endpoints','shrink');

ds = [0; cumsum(hypot(diff(xs),diff(ys)))];
ty = gradient(ys,ds);
tx = gradient(xs,ds);

% ---- which points get a cut ---------------------------------------------------
margin = s.endMarginRadii*median(rr);
inside = ds>=margin & ds<=ds(end)-margin;
if ~any(inside), q.why='nothing left after the branch-point margins'; return; end

if s.cutSpacing<=0
    take = find(inside);
else
    take  = zeros(1,0);
    aNext = ds(find(inside,1,'first'));
    for k = find(inside(:)).'
        if ds(k)>=aNext
            take(end+1) = k; %#ok<AGROW> - at most one per fitted point
            aNext = ds(k)+s.cutSpacing;
        end
    end
end
if numel(take)<s.minCuts, q.why='too few cuts fit along the segment'; return; end

% ---- one cut per taken point --------------------------------------------------
for k = take(:).'
    t = [ty(k) tx(k)];
    if ~all(isfinite(t)) || norm(t)==0, continue; end
    t = t./norm(t);
    n = [-t(2) t(1)];

    span = s.spanRadii*rr(k) + s.spanPadPx;
    if ys(k)-span<1 || ys(k)+span>nr || xs(k)-span<1 || xs(k)+span>nc
        continue                                  % dropped, never clipped
    end
    cuts(end+1) = struct('center',[ys(k) xs(k)],'normal',n, ...
        'radius',rr(k),'span',span,'arc',ds(k)); %#ok<AGROW> - one per taken point
end

if numel(cuts)<s.minCuts
    cuts = cuts([]);  q.why = 'too few cuts fit inside the frame';
end
end

% =====================================================================
function [py,px] = longestPath(bw)
%longestPath  The two most distant pixels of the label, and the path between them.
%   Farthest-from-anywhere, then farthest-from-there: on a tree that pair IS the
%   diameter, in two geodesic transforms rather than one per endpoint pair.  The path
%   is then a steepest descent on the distance from the first, which cannot get stuck
%   because every pixel of a connected mask has a neighbour strictly closer to the seed.
py = [];  px = [];
lin = find(bw);
if numel(lin)<2, return; end

d  = bwdistgeodesic(bw,lin(1),'quasi-euclidean');
iA = farthest(d);
if isempty(iA), return; end
dA = bwdistgeodesic(bw,iA,'quasi-euclidean');
iB = farthest(dA);
if isempty(iB), return; end

[nr,nc] = size(bw);
[oy,ox] = ndgrid(-1:1,-1:1);
off     = oy(:)+nr*ox(:);
path    = iB;
cur     = iB;
guard   = numel(lin)+5;
while dA(cur)>0 && numel(path)<guard
    [cy,cx] = ind2sub([nr nc],cur);
    if cy<=1 || cy>=nr || cx<=1 || cx>=nc
        nb = neighboursSafe(cur,nr,nc);
    else
        nb = cur+off;
    end
    nb = nb(bw(nb));
    [~,j] = min(dA(nb));
    if isempty(j) || dA(nb(j))>=dA(cur), break; end
    cur  = nb(j);
    path = [path; cur]; %#ok<AGROW> - bounded by guard
end
[py,px] = ind2sub([nr nc],flipud(path));
end

% =====================================================================
function i = farthest(d)
%farthest  The reachable pixel with the largest geodesic distance.
d(~isfinite(d)) = -1;
[m,i] = max(d(:));
if m<0, i = []; end
end

% =====================================================================
function nb = neighboursSafe(k,nr,nc)
%neighboursSafe  The eight neighbours of a pixel on the border, without wrapping.
[y,x]   = ind2sub([nr nc],k);
[yy,xx] = ndgrid(max(y-1,1):min(y+1,nr), max(x-1,1):min(x+1,nc));
nb = sub2ind([nr nc],yy(:),xx(:));
end
%------------- END OF CODE --------------
