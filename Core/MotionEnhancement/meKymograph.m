%meKymograph - A cut across a vessel, and where its two walls are to a fraction of a pixel.
%
%   THE MEASUREMENT THE WHOLE BLOCK IS JUDGED BY.  A magnified movie is persuasive
%   whether or not it is true, so the deliverable is never the movie - it is the
%   movie beside a number, and this is where the number comes from.  The same
%   function reads the raw stack and the magnified one, which is what makes their
%   ratio mean something: any bias the estimator has cancels between them.
%
%   THE WALL IS THE HALF-WAY CROSSING, READ BY LINEAR INTERPOLATION.  The vessel is
%   bright and the tissue is not, so the wall is where the profile crosses half way
%   between the two.  Interpolating between the two samples that straddle that level
%   is what makes the answer sub-pixel; taking the nearest sample would quantise
%   every displacement to the sample spacing and hide exactly the motion this block
%   exists to see.  On a Gaussian-blurred edge the half-way level sits at the
%   inflection, where the curvature is zero, so the interpolation is not merely
%   convenient - it is unbiased there.
%
%   EACH WALL IS MEASURED AGAINST THE TISSUE ON ITS OWN SIDE.  Real tissue is not
%   equally bright on both sides of a vessel, and one shared half-way level then sits
%   too high for one wall and too low for the other.  That is not extra noise: it is
%   a bias on the diameter, and a bias that drifts with the background is
%   indistinguishable from wall motion.
%
%   TWO NUMBERS COME OUT, AND WHICH ONE IS THE SIGNAL DEPENDS ON THE MOTION.  A
%   pulsating vessel moves its two walls in opposite directions, so the signal is the
%   HALF-WIDTH, (right-left)/2, and the centre is still.  A translating field moves
%   both the same way, so the signal is the CENTRE, (left+right)/2, and the
%   half-width is still.  Both are returned; the caller says which it meant.
%
%   meProbe SAMPLES ITS OWN CUTS.  It reads through a memmap and chooses between two
%   working copies of different scales, which is a decision this function has no
%   business knowing about - so it hands the kymograph it already built straight to
%   the second syntax below and shares the fitting, which is the part that must not
%   have two versions.
%
% Syntax:
%    kymo = meKymograph(stack, cut)
%    kymo = meKymograph(K, pos)
%
% Inputs:
%    stack - [rows columns frames], or one [rows columns] image.
%    cut   - .center [row column] on the vessel's centre line
%            .normal [drow dcolumn] across it; normalised here, so any length does
%            .radius the vessel's half-width in pixels, used to size the cut
%            .span   optional, half-length of the cut in pixels.  Default
%                    3*radius+6, which is meProbe's, so the two agree.
%            .step   optional, sample spacing along the cut.  Default 0.25 px.
%    K     - a kymograph that is already sampled, [nPos frames].
%    pos   - the positions its rows sit at, in pixels from the centre line.
%
% Outputs:
%    kymo - .K        the kymograph, [nPos frames]
%           .pos      [nPos 1] position along the cut, pixels from the centre line
%           .left     [frames 1] the near wall,  pixels along pos
%           .right    [frames 1] the far wall
%           .halfWidth (right-left)/2 - the signal when the vessel pulsates
%           .centre    (left+right)/2 - the signal when the field translates
%           .widthL,.widthR the 10 to 90 per cent edge widths, one per frame
%
% Example:
%    cut  = struct('center',[112 112],'normal',[0 1],'radius',15);
%    kymo = meKymograph(stack, cut);
%    plot(kymo.halfWidth - mean(kymo.halfWidth));
%
% Dependencies: none beyond core MATLAB.
% See also: mePhantom, meValidate, meProbe, meShift
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function kymo = meKymograph(X, cut)

if isstruct(cut)
    [K,pos] = sampleCut(X, cut);
else
    K   = double(X);
    pos = double(cut(:));
end
if size(K,1) ~= numel(pos)
    error('meKymograph:pos', ...
        'The kymograph has %d samples along the cut and %d positions were given.', ...
        size(K,1), numel(pos));
end

n = size(K,2);
kymo           = struct();
kymo.K         = K;
kymo.pos       = pos;
kymo.left      = nan(n,1);
kymo.right     = nan(n,1);
kymo.widthL    = nan(n,1);
kymo.widthR    = nan(n,1);
for t = 1:n
    [kymo.left(t), kymo.right(t), kymo.widthL(t), kymo.widthR(t)] = fitWalls(K(:,t), pos);
end
kymo.halfWidth = (kymo.right - kymo.left)./2;
kymo.centre    = (kymo.right + kymo.left)./2;
end

% =====================================================================
function [K,pos] = sampleCut(X, cut)
%sampleCut  Intensity along the cut for every frame.  The interpolation weights do
%   not change from frame to frame, so they are built once and applied as a single
%   indexed multiply rather than one interp2 call per frame.
nrm = cut.normal(:).'./norm(cut.normal);
if isfield(cut,'span')  && ~isempty(cut.span),  L    = cut.span;  else, L    = 3*cut.radius+6; end
if isfield(cut,'step')  && ~isempty(cut.step),  step = cut.step;  else, step = 0.25;           end

pos = (-L:step:L)';
pt  = cut.center(:).' + pos.*nrm;

[nr,nc,nT] = size(X);
if any(pt(:,1)<1) || any(pt(:,1)>nr) || any(pt(:,2)<1) || any(pt(:,2)>nc)
    error('meKymograph:outside', ...
        'The cut runs off a %dx%d frame; shorten span or move center.',nr,nc);
end

[idx,wgt] = bilinearWeights(pt(:,1), pt(:,2), nr, nc);
Xv = reshape(double(X), nr*nc, nT);
K  = wgt(:,1).*Xv(idx(:,1),:) + wgt(:,2).*Xv(idx(:,2),:) + ...
     wgt(:,3).*Xv(idx(:,3),:) + wgt(:,4).*Xv(idx(:,4),:);
end

% =====================================================================
function [idx,wgt] = bilinearWeights(r, c, nr, nc)
%bilinearWeights  The four pixels each sample point sits between, and their share.
r  = min(max(r,1),nr-1e-6);  c = min(max(c,1),nc-1e-6);
r0 = floor(r); c0 = floor(c);
fr = r-r0;     fc = c-c0;
idx = [sub2ind([nr nc],r0,c0), sub2ind([nr nc],r0+1,c0), ...
       sub2ind([nr nc],r0,c0+1), sub2ind([nr nc],r0+1,c0+1)];
wgt = [(1-fr).*(1-fc), fr.*(1-fc), (1-fr).*fc, fr.*fc];
end

% =====================================================================
function [left,right,widthL,widthR] = fitWalls(p, pos)
%fitWalls  Where the two walls are, to a fraction of a pixel.
left = NaN; right = NaN; widthL = NaN; widthR = NaN;
p = movmean(p,5);
n = numel(p);
m = max(1,round(0.15*n));
baseL = median(p(1:m));
baseR = median(p(end-m+1:end));
[peak,iPk] = max(p(m:end-m+1));
iPk = iPk + m - 1;
if ~(peak > baseL) || ~(peak > baseR), return; end

left   = crossOut(p, iPk, -1, baseL + 0.5*(peak-baseL), pos);
right  = crossOut(p, iPk, +1, baseR + 0.5*(peak-baseR), pos);
l10    = crossOut(p, iPk, -1, baseL + 0.1*(peak-baseL), pos);
l90    = crossOut(p, iPk, -1, baseL + 0.9*(peak-baseL), pos);
r10    = crossOut(p, iPk, +1, baseR + 0.1*(peak-baseR), pos);
r90    = crossOut(p, iPk, +1, baseR + 0.9*(peak-baseR), pos);
widthL = abs(l10-l90);
widthR = abs(r10-r90);
end

% =====================================================================
function x = crossOut(p, iPk, step, level, pos)
%crossOut  Walk outwards from the peak to the first crossing of level, and
%   interpolate between the two samples that straddle it.
x = NaN;
i = iPk;
while i+step>=1 && i+step<=numel(p)
    if p(i+step) <= level
        d = p(i)-p(i+step);
        if d==0, x = pos(i+step); else, x = pos(i) + (pos(i+step)-pos(i))*(p(i)-level)/d; end
        return
    end
    i = i+step;
end
end
%------------- END OF CODE --------------
