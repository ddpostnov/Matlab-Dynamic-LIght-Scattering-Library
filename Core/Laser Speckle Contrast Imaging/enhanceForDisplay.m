%enhanceForDisplay - background-subtracted vessel emphasis for display/overlays
%
%   imgOut = enhanceForDisplay(img, fSize)
%   imgOut = enhanceForDisplay(img, fSize, medWin)
%
% DESCRIPTION
%   Removes the slowly varying background of a grayscale frame so vessels stand
%   out for on-screen display and segmentation overlays.  A median-smoothed
%   morphological opening estimates the background, which is subtracted from the
%   image:
%
%       imgOut = img - imopen(medfilt2(img,[medWin medWin],'symmetric'), ...
%                             strel('disk',fSize));
%
%   This is the display-enhancement idiom shared verbatim by the segmentation path
%   (getPixelCategories, setRegions, showSegmentsPreview) and by setVesselTypes and
%   setVascularTree.  Every caller derives the disk radius the same way -
%   fSize = floor(min(size(img))/20)*2+1 - but the median window differs: the
%   categorization path (getPixelCategories / setRegions) caps it at 15
%   (medWin = min(15,fSize)), while the others use the full fSize.  medWin is
%   therefore exposed so each caller reproduces its exact output bit-for-bit.
%
%   Unlike its cousin enhanceForRegistration, this routine does NOT clip, invert,
%   smooth or normalize: callers apply their own percentile scaling / inversion
%   beforehand and their own masking afterwards.  The two are kept separate -
%   enhanceForRegistration has its own contract.
%
% INPUT
%   img     2-D numeric image (already scaled/inverted by the caller as needed).
%   fSize   disk radius (px) of the morphological-opening structuring element,
%           and the default median-filter window.  Positive-integer scalar.
%
% OPTIONAL
%   medWin  side length (px) of the square median-filter window.  Default fSize;
%           getPixelCategories / setRegions pass min(15,fSize).
%
% OUTPUT
%   imgOut  background-subtracted image, same size and class as img.  No masking
%           or normalization is applied here.
%
% DEPENDS ON
%   MATLAB's Image Processing Toolbox (medfilt2, imopen, strel).
%
% See also: enhanceForRegistration, runSegmentation, setVesselTypes, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

function imgOut = enhanceForDisplay(img, fSize, medWin)

if nargin < 3 || isempty(medWin)
    medWin = fSize;
end

imgOut = img - imopen(medfilt2(img,[medWin medWin],'symmetric'),strel('disk',fSize));
end
