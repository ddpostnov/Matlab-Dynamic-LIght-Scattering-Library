%enhanceForRegistration - vessel-emphasis preprocessing for image registration
%
%   e = enhanceForRegistration(img)
%   e = enhanceForRegistration(img, Name,Value, ...)
%
% DESCRIPTION
%   Turns a grayscale frame (LSCI contrast, DLSI, intensity, fluorescence, ...)
%   into a normalized, vessel-emphasized feature image suitable for intensity-
%   or correlation-based registration.  Vessels are usually dark and low
%   contrast; this routine inverts them to bright, removes the slowly varying
%   background and normalizes the result so a registration cost function locks
%   onto the vasculature rather than the illumination profile.
%
%   A single skeleton is shared by every caller:
%     1) clip the image to its [1, 99] percentiles (computed within 'mask');
%     2) invert (imcomplement) so the dark vessels become bright - this happens
%        BEFORE the morphological open, which is the load-bearing ordering:
%        opening a bright-vessel image estimates and removes the background but
%        keeps the vessels, whereas opening the raw (dark-vessel) image would
%        erase them;
%     3) subtract a median-smoothed morphological opening as the background
%        estimate (disk radius = median kernel = floor(min(size)/20)*2+1);
%     4) optionally suppress residual speckle with a Gaussian ('smoothSigma');
%     5) normalize to [0, 1] with mat2gray - optionally within 'mask' with a
%        histogram equalization ('equalizeInMask').
%
%   Two option sets reproduce, to within floating point, the two legacy in-line
%   recipes this primitive replaces: the retina keyframe enhancement
%   (smoothSigma = 1, no mask) and the LSCI-to-LSCI wrapper enhancement
%   (mask + equalizeInMask).
%
% INPUT
%   img   2-D numeric image.  Non-finite pixels are set to 0 before processing.
%
% OPTIONAL Name-Value
%   'mask'           logical (or 0/1) array the size of img.  The percentile
%                    clip is computed within it and, when 'equalizeInMask' is
%                    true, the output is masked and equalized within it.
%                    Default: true(size(img)) (whole image).  Must contain at
%                    least one true element.
%   'smoothSigma'    Gaussian sigma (px) applied after background subtraction to
%                    suppress residual speckle.  0 = no smoothing.  Default 0.
%   'equalizeInMask' logical.  false (default): normalize the whole image with
%                    mat2gray and return single.  true: normalize with mat2gray,
%                    multiply by 'mask', then histogram-equalize the in-mask
%                    pixels and return double (out-of-mask pixels = 0).
%
% OUTPUT
%   e   enhanced feature image in [0, 1].  single when equalizeInMask is false,
%       double (out-of-mask pixels = 0) when it is true.
%
% NOTE
%   mat2gray is a monotone affine rescale and commutes with the median / opening
%   rank filters, so a normalization placed before the inversion would be
%   cosmetic (the final mat2gray cancels it); the invariant that actually
%   matters is that imcomplement precedes the morphological open.
%
% DEPENDS ON
%   MATLAB's Image Processing Toolbox (imcomplement, medfilt2, imopen, strel,
%   imgaussfilt, mat2gray, histeq).
%
% See also: registerRetinaLSCI, registerLSCItoLSCI, imregtform, imregcorr
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 23-July-2026

function e = enhanceForRegistration(img, varargin)

p = inputParser;
p.FunctionName = 'enhanceForRegistration';
addRequired(p,'img',@(x) isnumeric(x) && ismatrix(x));
addParameter(p,'mask',[],@(x) (islogical(x) || isnumeric(x)) && ismatrix(x));
addParameter(p,'smoothSigma',0,@(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'equalizeInMask',false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x) && ismember(x,[0 1])));
parse(p,img,varargin{:});
mask           = p.Results.mask;
smoothSigma    = double(p.Results.smoothSigma);
equalizeInMask = logical(p.Results.equalizeInMask);

img = double(img);
img(~isfinite(img)) = 0;

if isempty(mask)
    mask = true(size(img));
else
    if ~isequal(size(mask),size(img))
        error('enhanceForRegistration:maskSize','mask must be the same size as img.');
    end
    mask = logical(mask);
    if ~any(mask(:))
        error('enhanceForRegistration:emptyMask','mask must contain at least one true element.');
    end
end

% ---- (1) percentile clip within the mask ----------------------------------
inMask = img(mask);
lo = prctile(inMask,1);
hi = prctile(inMask,99);
if hi<=lo, hi = lo+eps; end
img = min(max(img,lo),hi);

% ---- (2) invert BEFORE the open: dark vessels -> bright --------------------
img = imcomplement(img);

% ---- (3) subtract the median-smoothed morphological-opening background -----
fSize = floor(min(size(img))/20)*2+1;
img = img - imopen(medfilt2(img,[fSize fSize],'symmetric'),strel('disk',fSize));

% ---- (4) optional residual-speckle suppression -----------------------------
if smoothSigma>0
    img = imgaussfilt(img,smoothSigma);
end

% ---- (5) normalize ---------------------------------------------------------
if equalizeInMask
    img = mat2gray(img).*mask;
    img(mask) = histeq(img(mask));
    e = img;
else
    e = single(mat2gray(img));
end
end
