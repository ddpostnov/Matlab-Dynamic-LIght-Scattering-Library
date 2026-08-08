%meMask - Where a vessel wall is, and therefore where amplification is allowed.
%
%   THE DESIGN ELEMENT WITH THE LARGEST EFFECT ON WHETHER THE OUTPUT CAN BE TRUSTED.
%   The label is FITC in the plasma, so the lumen is bright and red blood cells are
%   dark objects crossing it.  At 2.5 micrometres a pixel and 10 to 20 mm/s that
%   texture moves forty to eighty pixels between frames - far beyond the phase range
%   of any pyramid level, so it cannot leave a coherent direction behind and lands as
%   broadband phase noise inside the vessel.  Awkward, not misleading.  What IS
%   misleading is the lumen-wall boundary: where a dark cell reaches the wall the
%   local gradient structure genuinely changes, and that is a real phase change that
%   is not wall motion.  Amplifying only where a wall is analyses the lumen and does
%   not magnify it.  Wadhwa et al., SIGGRAPH 2013 do exactly this, masking to the iris.
%
%   THE MASK GOES ON THE SHIFT, NEVER ON THE ANALYSIS.  Everything is measured
%   everywhere - the phase, the amplitude, the spatial mean that meGlobal removes -
%   and only the amplification is confined.  Masking the analysis instead would change
%   the bulk-motion estimate, which is a whole-field quantity, and would make the
%   answer inside the mask depend on where the mask was drawn.
%
%   WHERE IS A DECISION; HOW MUCH IS alpha's JOB.  The weight is one over the whole
%   of the mask and falls to zero over a couple of pixels at its border - it is not a
%   normalised copy of the gradient.  That distinction is the difference between a
%   mask and a distortion: scale the shift by the gradient and a stretch of wall with
%   a strong edge is magnified more than a stretch with a weak one, so a wall moving
%   uniformly comes out moving unevenly, which is a false statement about the very
%   thing the movie is made to show.  Measured on the reference cine, the graded
%   version also halved the amplification on the wall it was supposed to protect.
%
%   THE BORDER IS STILL SOFT, AND THAT IS NOT COSMETIC EITHER.  The phase shift is
%   scaled by this weight before it is exponentiated, so a hard-edged mask puts a step
%   in the shift field, and a step in a shift field is a tear in the picture at
%   exactly the place a reader is looking.  The Gaussian at the end is what keeps it
%   continuous.
%
%   TWO SYNTAXES, TWO DIFFERENT JOBS.  meMask(s, meanImage) BUILDS the weight from a
%   picture of the field; meMask(s, pyr) RESAMPLES an already-built weight onto each
%   pyramid level, which is what meShift calls on every run.  They are one function
%   because a mask that was built one way and resampled another is a mask nobody can
%   reason about.
%
%   EXCLUSION IS BUILT AND IT IS PROBABLY NOT NEEDED HERE.  The plan wanted to keep a
%   large, visibly moving vessel out of the mask so that one α could serve the rest of
%   the field.  Session 0 measured the field and nothing in it moves multiple pixels
%   at the cardiac frequency - the largest vessel's wall moves 0.075 px against a
%   lambda/4 ceiling of 212 at level 5, and what the eye sees moving is bulk
%   preparation motion.  s.maskExclude is therefore cheap rather than central, and it
%   is worth having because it is also how a vessel is held back as a CONTROL.
%
% Syntax:
%    w  = meMask(s, meanImage)
%    wl = meMask(s, pyr)
%
% Inputs:
%    s         - settings.  Fields read when BUILDING:
%                .maskSource  'gradient' or 'segmentation'.
%                .maskSmooth  pixels of Gaussian smoothing.
%                .maskThresh  fraction of the gradient's robust maximum below which
%                             there is no wall.  'gradient' only.
%                .maskRing    half-width of the ring.  'segmentation' only.
%                .maskVessel  a logical vessel mask the ring is built around.
%                             'segmentation' only - mePassA writes one into
%                             <rec>_ME_mean.mat, and the library's own segmentation
%                             can supply a better one for the same field.
%                .maskExclude a logical map of places NOT to amplify, or empty.
%                Fields read when RESAMPLING:
%                .maskMode    'none' amplify everywhere; 'weight' use s.mask.
%                .mask        the weight, the size of the frame being magnified.
%    meanImage - [rows columns] a picture of the field: the time mean, or the cine's
%                own base frame, which is the same thing over the accepted beats.
%    pyr       - the struct from meRieszPyramid.
%
% Outputs:
%    w  - [rows columns] single, in [0,1], one where a wall is.
%    wl - {1 x nLevels}.  An EMPTY entry means no weight at that level, which the
%         shift reads as one everywhere - not as zero.
%
% Example:
%    mn = load(pass.meanPath);
%    s  = meSettings;  s.maskMode = 'weight';
%    s.mask = meMask(s, mn.meanFrame);
%    mag = meShift(stack, fs, s);
%
% Dependencies: Image Processing Toolbox (imgaussfilt, imgradient, imresize,
%               imdilate, imerode, strel).
% See also: meShift, meEnhance, meRieszPyramid, mePassA, meSettings
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function w = meMask(s, x)

if isstruct(x)
    w = resampleToLevels(s, x);
else
    w = buildWeight(s, x);
end
end

% =====================================================================
function w = buildWeight(s, img)
%buildWeight  A wall weight in [0,1] from a picture of the field.

img = double(img);
lo  = prctile(img(:), 1);
hi  = prctile(img(:), 99.9);
img = min(max((img-lo)./max(hi-lo,eps), 0), 1);

switch lower(s.maskSource)

    case 'gradient'
        % THE GRADIENT IS TAKEN OF A SMOOTHED PICTURE, and that ordering is the whole
        % of it: the gradient of a noisy image is noise everywhere, at an amplitude
        % that has nothing to do with whether a wall is there.  The wall edge is 6 to
        % 15 px wide on this recording, so a couple of pixels of blur costs the edge
        % nothing and costs the noise a great deal.
        g = imgradient(imgaussfilt(img, s.maskSmooth));

        % Normalised on a high percentile rather than the maximum.  One saturated
        % pixel or one dust speck sets the maximum, and then the whole mask is scaled
        % by an artefact.  What comes out of the threshold is a DECISION - wall or not
        % - and not a graded copy of the gradient; see the header.
        g = g./max(prctile(g(:), 99.5), eps);
        w = double(g >= s.maskThresh);

    case 'segmentation'
        if isempty(s.maskVessel)
            error('meMask:vessel', ...
                ['maskSource is segmentation and s.maskVessel is empty. mePassA ' ...
                 'writes a vessel mask into <rec>_ME_mean.mat, and the library''s ' ...
                 'own segmentation can supply a better one for the same field.']);
        end
        v = logical(s.maskVessel);
        if ~isequal(size(v), size(img))
            error('meMask:vesselSize','The vessel mask is %dx%d and the frame is %dx%d.', ...
                size(v,1), size(v,2), size(img,1), size(img,2));
        end
        % Dilate minus erode: the shell the wall lives in.  A bright-lumen vessel is a
        % ridge and not an annulus, so there is nothing inside to fill and nothing
        % here tries to - imfill on a vascular network welds the tree into one blob
        % and fills the interstices, which cost Session 0 an afternoon.
        se = strel('disk', max(round(s.maskRing),1));
        w  = double(imdilate(v,se) & ~imerode(v,se));

    otherwise
        error('meMask:source','The mask source is gradient or segmentation, not %s.', ...
            s.maskSource);
end

% ---- soften the border, exclude, normalise ----------------------------------
% A Gaussian of this width leaves the interior of the mask at one and ramps to zero
% over about three sigma at its border, which is what makes the shift field
% continuous without grading the amplification across the wall.
w = imgaussfilt(w, s.maskSmooth);

if ~isempty(s.maskExclude)
    e = double(logical(s.maskExclude));
    if ~isequal(size(e), size(img))
        error('meMask:excludeSize','The exclusion is %dx%d and the frame is %dx%d.', ...
            size(e,1), size(e,2), size(img,1), size(img,2));
    end
    % Softened the same way, so an excluded vessel fades out instead of leaving a
    % hole with a step around it.
    w = w.*max(1 - imgaussfilt(e, s.maskSmooth), 0);
end

w = single(w./max(max(w(:)), eps));
end

% =====================================================================
function wl = resampleToLevels(s, pyr)
%resampleToLevels  The built weight, on each pyramid level's own grid.

wl = cell(1, numel(pyr.lap));

switch lower(s.maskMode)
    case 'none'
        return

    case 'weight'
        if isempty(s.mask)
            error('meMask:mask','maskMode is weight and s.mask is empty.');
        end
        if ~isequal([size(s.mask,1) size(s.mask,2)], pyr.imSize)
            error('meMask:size', ...
                'The mask is %dx%d and the frame is %dx%d.', ...
                size(s.mask,1), size(s.mask,2), pyr.imSize(1), pyr.imSize(2));
        end
        m = single(s.mask);
        for k = pyr.levels
            wk    = imresize(m, pyr.size{k}, 'bilinear');
            wl{k} = min(max(wk,0),1);
        end

    otherwise
        error('meMask:mode','The mask mode is none or weight, not %s.', s.maskMode);
end
end
%------------- END OF CODE --------------
