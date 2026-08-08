%getBackground - The glow around the vessels in a fluorescence frame, in counts.
%
%   THE BACKGROUND OF A WIDE-FIELD FLUORESCENCE FRAME IS THE VESSELS' OWN LIGHT,
%   SPREAD SIDEWAYS BY THE TISSUE.  Before the dye arrives there is no halo outside
%   the vessels at all (measured: an excess at the wall of 18.7 counts on a 230-count
%   floor, R2 0.14, i.e. nothing), and it bleaches at the vessels' own rate.  Every
%   count of haze is dye.  It falls off with two lengths, about 7 um and 120-160 um,
%   measured the same on two recordings taken seven months apart, and at a vessel wall
%   it is half to six tenths of what is inside the vessel.
%
%   SO REMOVING IT NECESSARILY REMOVES MOST OF A SMALL VESSEL'S SIGNAL, and that is
%   the physics rather than a defect of this filter.  A large part of the light
%   measured inside a vessel IS halo - its neighbours' and its own - so a removal that
%   takes 99.8% of the halo keeps about 41% of the in-vessel first-pass rise, and
%   below 20 um calibre it keeps 14-20%.  Every method measured behaves this way,
%   including the strictly linear ones.  Anything downstream that needs an AMPLITUDE
%   from a small vessel should read the raw frames.
%
%   ONE METHOD, NOT SIX.  Research session R1 measured a Gaussian low pass,
%   Sternberg's rolling ball, a per-pixel temporal baseline, a Fourier deconvolution
%   of the measured scattering kernel and a rank-one split beside this one.  The
%   opening won on the measured trade and the others are in
%   claude-docs/intensity-branch/research-background-removal.md if the question is
%   reopened.  A method switch here would be five untested paths in the library and a
%   settings field with five spellings.
%
%   EVERY LENGTH IS IN MICROMETRES AND THIS IS THE ONLY PLACE THEY BECOME PIXELS.  A
%   structuring element of 151 pixels means nothing at a different magnification, and
%   the halo this removes is a property of the tissue and not of the sensor.
%   pixelSize is therefore REQUIRED and has no default.
%
%   THE RADIUS IS TIED TO THE VASCULATURE, NOT TO THE OPTICS.  200 um leaves 0.2% of
%   the halo where the scratch cell's 151 px (377.5 um at 2.5 um/px) leaves 10.9%, and
%   it costs one percentage point of vessel signal to do it.  But it was measured on a
%   field whose largest vessel is 71-75 um; on a preparation with a 200 um surface
%   artery the opening slides under the artery and the number has to be re-measured.
%
%   IT ASSUMES THE VESSELS ARE THE BRIGHT THINGS, AND THAT IS A REAL LIMIT ON IT.  A
%   morphological OPENING removes small BRIGHT structures and keeps everything else,
%   so on a recording whose vessels are DARKER than the tissue they would survive the
%   opening, land INSIDE the estimate, and the subtraction would take the vasculature
%   rather than the glow - while producing a perfectly ordinary-looking picture.  This
%   library processes bright-vessel fluorescence recordings only (author, 2026-08-08),
%   so that case does not arise here and there is nothing to test for.  The dual is
%   one line - imclose instead of imopen, i.e. invert, estimate, invert back - and
%   whoever writes it owes the measurement R1 made for this arm, the halo left and the
%   vessel kept on a dark-vessel recording, and not a plausible picture.
%
% Syntax:
%    bg        = getBackground(data,s)
%    [bg,info] = getBackground(data,s)
%
% Inputs:
%    data     - [rows columns] or [rows columns frames], counts, any numeric class.
%    s        - .pixelSize  micrometres per pixel.  REQUIRED, no default.
%               .radiusUm   opening radius, micrometres.  Default 200.
%               .medianUm   median-filter window, micrometres.  Default 37.5.
%
% Outputs:
%    bg   - the estimated glow, the same size as data, SINGLE, in counts.  The caller
%           subtracts; this function never returns the cleaned data, so what is
%           stored beside a product is the estimate and the operation stays exactly
%           reversible.
%    info - .radiusPx .medianPx (as used, the median window forced odd),
%           .pixelSize .radiusUm .medianUm .seconds.
%
% Example:
%    s.pixelSize = 2.5;  s.radiusUm = 200;  s.medianUm = 37.5;
%    bg    = getBackground(results.imgI,s);
%    clean = results.imgI - double(bg);
%
% Dependencies: Image Processing Toolbox (medfilt2, imopen, strel).
% See also: runBackgroundRemoval, getHaloProfile, runCTTH
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function [bg,info] = getBackground(data,s)

tCall = tic;

if ~isfield(s,'pixelSize') || isempty(s.pixelSize)
    error('getBackground:pixelSize', ...
        ['s.pixelSize is required: every length this step works in is in ' ...
         'micrometres, and this is the only place they become pixels.']);
end
px = double(s.pixelSize);
if ~isscalar(px) || ~isfinite(px) || px<=0
    error('getBackground:pixelSize', ...
        's.pixelSize must be a positive number of micrometres per pixel; it is %s.', ...
        mat2str(s.pixelSize));
end
if ~isfield(s,'radiusUm') || isempty(s.radiusUm), s.radiusUm = 200;  end
if ~isfield(s,'medianUm') || isempty(s.medianUm), s.medianUm = 37.5; end

% The two lengths in pixels.  The median window is forced ODD because medfilt2 wants
% one and a half-pixel offset in a despeckling window is not something a caller should
% have to think about; the radius is at least one pixel, so a pixelSize larger than
% the radius degenerates to a 1-pixel opening rather than to an empty strel.
mW = 2*floor(double(s.medianUm)/px/2)+1;
rO = max(1,round(double(s.radiusUm)/px));

% strel is built ONCE.  On a 1361-frame cube rebuilding a disk of radius 80 inside the
% loop is the difference between 80 s and several minutes, for the same answer.
se = strel('disk',rO);

[nR,nC,nT] = size(data,1,2,3);
bg = zeros(nR,nC,nT,'single');
for k = 1:1:nT
    % NO mat2gray, EVER, AND THIS IS WHERE IT WOULD GO.  The scratch cell this recipe
    % comes from rescales every frame to its own minimum and maximum before the
    % opening, which turns the frame's own brightness into a signal.  R1 ran the two
    % recipes side by side differing in nothing else: the per-frame rescale adds
    % 1881 ms of arrival bias to a bolus transit time.  It is the obvious
    % "improvement" and it is measured to be wrong.
    bg(:,:,k) = imopen( medfilt2(single(data(:,:,k)),[mW mW],'symmetric'), se );
end

info = struct('radiusPx',rO,'medianPx',mW,'pixelSize',px, ...
    'radiusUm',double(s.radiusUm),'medianUm',double(s.medianUm), ...
    'seconds',toc(tCall));
end
%------------- END OF CODE --------------
