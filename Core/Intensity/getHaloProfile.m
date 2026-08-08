%getHaloProfile - How much light sits around the vessels, and how far out it reaches.
%
%   THIS IS THE NUMBER BESIDE THE PICTURE.  A background removal always makes a
%   prettier image, so a report page showing only the before and the after is not
%   evidence of anything; what says the removal took the glow rather than the
%   vasculature is two measurements on the same field, against the same vessels: how
%   much of the halo is left, and how much of the in-vessel signal is left.  Research
%   session R1's whole method is that pair, and this is where it lives in the library.
%
%   THE PROFILE IS A MEDIAN, NOT A MEAN.  At a given distance from a vessel some
%   pixels sit on a second vessel the mask missed, and a mean over that bin follows
%   the vessel rather than the halo.  The median does not.
%
%   THE HEADLINE NUMBER HAS NO MODEL IN IT.  grad is the measured drop from the 10 um
%   ring to the far end of whatever distance range the field actually offers, and the
%   RATIO of two grads is what says how much halo a removal took.  R1's prototype also
%   fitted a sum of two exponentials, and that part is deliberately NOT here: on a
%   field where nothing is more than about 250 um from a vessel the five-parameter
%   model is unidentifiable, and left unbounded it returned amplitudes of 1.8e10
%   counts that cancelled to a plausible-looking curve.  A fit that can do that cannot
%   be what a recommendation rests on.  The two decay lengths are in
%   claude-docs/intensity-branch/research-background-removal.md, measured on a whole
%   field, where they are identifiable.
%
%   THE MASK IS DERIVED ONCE AND HANDED BACK, and that is the point of the two-call
%   shape.  A before and an after measured against masks derived from their own
%   pictures are two different geometries, so their ratio would not mean anything -
%   the removal changes the picture the mask would be derived FROM.  So the caller
%   derives the mask from the raw picture and measures both images against it.
%
%   IT IS NOT A SEGMENTATION AND NOTHING READS IT.  The mask below exists to place
%   distance bins for one measurement on one report page; it is never stored, never
%   compared with results.cMask, and the step that uses it runs before segmentation as
%   often as after.  A generous mask is what this measurement wants and a tight one is
%   what a segmentation wants, which is the other reason these are not the same thing:
%   a mask that clips a vessel wall puts vessel light into the first background bin
%   and the halo then appears to start several hundred counts high, which is the one
%   way this measurement can lie.  So the mask is dilated by a stated margin.
%
%   A BRIGHT-LUMEN VESSEL IS A RIDGE, NOT AN ANNULUS.  The top-hat removes the uneven
%   background the ridge sits on; nothing is hole-filled, because filling the holes of
%   a vascular network welds the tree into one solid region and the distance transform
%   then measures the field rather than the vessel.
%
% Syntax:
%    [p,ref] = getHaloProfile(img,pixelSize)   % derive the vessel mask from img
%    p       = getHaloProfile(img,ref)         % measure another image, SAME mask
%    [p,ref] = getHaloProfile(img,pixelSize,opt)
%
% Inputs:
%    img       - one image in counts, [rows columns].  Whatever image the question is
%                about: a mean picture, or the dye alone (a filled frame minus a
%                pre-bolus one), which is what R1 measured on.  A RATIO of two of
%                these numbers only means something when both came from the same kind
%                of image.
%    pixelSize - micrometres per pixel.  Every length here is in micrometres.
%    ref       - the second output of an earlier call: the mask, the distance map and
%                the pixel size to reuse.
%    opt       - (optional) .tophatUm   top-hat radius, um.        Default 100.
%                           .dilateUm   margin added to the mask.  Default 5.
%                           .minAreaUm2 smallest accepted object.  Default 2000.
%                           .maxUm      furthest distance profiled. Default 400.
%
% Outputs:
%    p   - .d        bin centres, micrometres
%          .value    median counts in each bin
%          .n        pixels in each bin
%          .dMax     furthest distance with enough pixels to be a measurement, um
%          .floor    measured level over the outer fifth of that range, counts
%          .atUm     [5 10 25 50 100 200 300]
%          .atValue  measured excess over .floor at those distances, counts
%          .grad     the halo amplitude: the 10 um ring's excess over .floor, counts.
%                    Floor-free by construction, so a method that shifts the whole
%                    image by a constant scores unchanged - which is correct, a
%                    constant is not haze.
%          .inVessel median over the UNDILATED vessel mask, counts.  The other half of
%                    the pair: a removal that takes the halo and the vessels together
%                    scores well on grad alone.
%          .ok       false when the field offered too few distance bins to measure
%                    anything, in which case every number above is NaN.
%    ref - .mask (dilated) .coreMask (not) .dUm .pixelSize .fraction, plus the three
%          mask lengths as used.
%
% Example:
%    [p0,ref] = getHaloProfile(imgRaw,   2.5);
%     p1      = getHaloProfile(imgClean, ref);
%    fprintf('halo left %.1f %%, vessel kept %.1f %%\n', ...
%        100*p1.grad/p0.grad, 100*p1.inVessel/p0.inVessel);
%
% Dependencies: Image Processing Toolbox (imtophat, imgaussfilt, imbinarize,
%               bwareaopen, imdilate, bwdist, strel).
% See also: getBackground, runBackgroundRemoval
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function [p,ref] = getHaloProfile(img,pixelSizeOrRef,opt)

if nargin<3, opt = struct(); end
if ~isfield(opt,'tophatUm')   || isempty(opt.tophatUm),   opt.tophatUm   = 100;  end
if ~isfield(opt,'dilateUm')   || isempty(opt.dilateUm),   opt.dilateUm   = 5;    end
if ~isfield(opt,'minAreaUm2') || isempty(opt.minAreaUm2), opt.minAreaUm2 = 2000; end
if ~isfield(opt,'maxUm')      || isempty(opt.maxUm),      opt.maxUm      = 400;  end

if isstruct(pixelSizeOrRef)
    ref = pixelSizeOrRef;
    if ~all(isfield(ref,{'mask','coreMask','dUm','pixelSize'}))
        error('getHaloProfile:ref', ...
            ['The second argument is either the micrometres per pixel or the ref ' ...
             'struct an earlier call handed back.']);
    end
    if ~isequal(size(ref.mask),size(img))
        error('getHaloProfile:refSize', ...
            ['The image is %s and the mask it is being measured against is %s.  ' ...
             'Both images have to be the same field for the comparison to mean ' ...
             'anything.'], mat2str(size(img)), mat2str(size(ref.mask)));
    end
else
    px = double(pixelSizeOrRef);
    if ~isscalar(px) || ~isfinite(px) || px<=0
        error('getHaloProfile:pixelSize', ...
            's.pixelSize must be a positive number of micrometres per pixel; it is %s.', ...
            mat2str(pixelSizeOrRef));
    end
    ref = vesselField(img,px,opt);
end

p = profileOf(img,ref,opt);
end

% =====================================================================
function ref = vesselField(img,px,opt)
%vesselField  Where the resolved vessels are, and how far every pixel is from one.
I  = double(img);
I(~isfinite(I)) = 0;
rT = max(1,round(opt.tophatUm/px));
th = imtophat(imgaussfilt(I,1), strel('disk',rT));
% imbinarize wants [0 1].  A top-hat is non-negative by construction (it is the image
% minus its own opening), so dividing by its maximum IS the normalisation - written
% out rather than as mat2gray so that the grep asserting no per-frame rescale ever
% reaches the data has no exception to make for this file.
bw = imbinarize(th./max(max(th(:)),eps));
bw = bwareaopen(bw, max(1,round(opt.minAreaUm2/px^2)));

core = bw;
rD   = max(0,round(opt.dilateUm/px));
if rD>0, bw = imdilate(bw, strel('disk',rD)); end

ref = struct('mask',bw,'coreMask',core,'dUm',bwdist(bw).*px,'pixelSize',px, ...
    'fraction',nnz(bw)/numel(bw),'tophatUm',opt.tophatUm,'dilateUm',opt.dilateUm, ...
    'minAreaUm2',opt.minAreaUm2);
end

% =====================================================================
function p = profileOf(img,ref,opt)
%profileOf  Median counts against distance to the nearest vessel, model free.
px = ref.pixelSize;
% The bin is never finer than one pixel: at 2.5 um/px it is R1's own 2.5 um bin, and
% at a coarser sampling a 2.5 um bin would hold either nothing or one pixel's worth of
% noise and the profile would be a staircase of empty bins.
binUm = max(2.5,px);

I = double(img(:));
d = double(ref.dUm(:));
keep = ~ref.mask(:) & d>0 & d<=opt.maxUm & isfinite(I);
I = I(keep); d = d(keep);

edges = 0:binUm:opt.maxUm;
[~,~,bin] = histcounts(d,edges);
nB  = numel(edges)-1;
val = nan(nB,1); cnt = zeros(nB,1);
for k = 1:1:nB
    v = I(bin==k);
    cnt(k) = numel(v);
    if ~isempty(v), val(k) = median(v); end
end
ctr = edges(1:end-1)'+binUm/2;

p = struct('d',ctr,'value',val,'n',cnt,'binUm',binUm,'maxUm',opt.maxUm);
p.atUm     = [5 10 25 50 100 200 300];
p.atValue  = nan(size(p.atUm));
p.inVessel = median(double(img(ref.coreMask)));

% A bin with a handful of pixels is a median of a handful of pixels; twenty is what
% R1's prototype required and it is what the numbers in that document were measured
% with.  Fewer than eight such bins is not a profile at all - the field is too small
% or the mask has swallowed it - and the honest answer then is NaN with ok false,
% never a number computed from three bins.
solid = cnt>=20 & isfinite(val);
if nnz(solid)<8
    p.dMax=NaN; p.floor=NaN; p.grad=NaN; p.ok=false;
    return
end

p.dMax  = max(ctr(solid));
outer   = solid & ctr>=0.80*p.dMax;
p.floor = median(val(outer));
for i = 1:1:numel(p.atUm)
    if p.atUm(i) > p.dMax, continue; end
    [~,j] = min(abs(ctr-p.atUm(i)));
    if solid(j), p.atValue(i) = val(j)-p.floor; end
end
[~,j10]  = min(abs(ctr-10));
p.grad   = val(j10)-p.floor;
p.ok     = true;
end
%------------- END OF CODE --------------
