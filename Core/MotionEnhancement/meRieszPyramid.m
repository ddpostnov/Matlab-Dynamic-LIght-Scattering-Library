%meRieszPyramid - Laplacian pyramid plus the 3-tap Riesz transform of each band.
%
%   THE REPRESENTATION THE WHOLE MAGNIFIER RUNS ON.  A Laplacian band and its two
%   Riesz components form a quaternion whose phase is the local position of the
%   structure at that scale; move the structure and the phase moves with it, so a
%   sub-pixel displacement becomes a number rather than a difference of two similar
%   pictures.  Wadhwa et al., ICCP 2014.
%
%   WHY THIS AND NOT A COMPLEX STEERABLE PYRAMID.  Four times over-complete instead
%   of twelve to fifty-six, and ONE orientation-free phase field per level instead
%   of eight oriented ones - which is eight times less temporal filtering, and
%   temporal filtering is where the time goes.  The price is the tighter octave
%   bound, alpha*delta < lambda/4 rather than lambda/2, and on this recording that
%   bound is nowhere near binding.
%
%   THE KERNELS ARE ANTISYMMETRIC, AND THAT IS A TRAP.  conv2 flips a kernel before
%   it slides it; imfilter does not.  These two differ by a sign, so the wrong
%   choice negates both Riesz components.  What that does is worth knowing exactly,
%   because it is not what one first fears: negating BOTH components conjugates the
%   quaternion, and the phase difference and the phase shift then conjugate
%   together, so the real part that comes out of meShift is unchanged.  A
%   consistent flip is harmless.  What is NOT harmless is an INCONSISTENT one - the
%   phase measured one way and the shift applied the other - or a swap of the two
%   kernels, which magnifies motion along x as motion along y.  Both are silent and
%   both produce an entirely plausible movie, which is why the direction is
%   asserted on a phantom rather than reasoned about.
%
%   imfilter IS THE CHOICE HERE, and with these kernels it computes the true Riesz
%   sign: the x component of the Riesz transform is MINUS the derivative along x
%   (its multiplier is -i*w/|w| against the derivative's +i*w), and imfilter with
%   [0 0 0; 0.5 0 -0.5; 0 0 0] returns 0.5*L(c-1) - 0.5*L(c+1), which is exactly
%   that.  Whoever changes this line owes the phantom assertion a re-run.
%
%   TWO RIESZ TRANSFORMS, AND WHICH ONE TO USE WAS DECIDED BY MEASUREMENT.  The
%   three taps are an APPROXIMATION: their response is -i*sin(w) where the Riesz
%   transform's is -i*w/|w|, so they are only right near one frequency of each band
%   and wrong at its edges.  What that costs, measured 2026-08-07, is not small.
%   Amplifying a translating grating with every level active reached 0.80 of the
%   1+alpha the theory promises with the three taps and 0.93 with the exact
%   multiplier; on a vessel wall, 0.73 against 0.81.  A longer lowpass helps
%   neither (a 7- and a 9-tap binomial came out unchanged), so the pyramid's band
%   overlap is a separate term and not the one that moves.
%
%   'exact' IS THE DEFAULT, AND THE PUBLISHED PSEUDOCODE'S THREE TAPS ARE NOT,
%   because on the calibration phantom the exact multiplier gave 24 % more
%   amplification (gain 8.76 against 7.05 at alpha 20) at a detection floor that
%   could not be told apart from it (0.0045 against 0.0051 px, inside the scatter
%   of a thirty-bin estimate) for 17 % more time in the analysis.  Better on the
%   axis that is scarce, no worse on the axis that matters, and cheap.
%
%   ONE RISK IN 'exact' IS UNMEASURED AND SESSION 3 SHOULD CHECK IT.  The exact
%   transform is a multiplier in the Fourier domain, so it treats the band as
%   PERIODIC where the three taps pad it symmetrically.  A Laplacian band is near
%   zero at the frame boundary, which is why this is affordable at all - but a
%   bright vessel running off the edge of a real frame is not near zero, and it can
%   ring across to the opposite side.  The phantom has no vessel at its edge, so
%   the phantom cannot see this.  If a real product shows structure that mirrors
%   across the frame, try '3tap' first.
%
%   'exact' costs two Fourier transforms per level per frame and holds a complex
%   copy of the band while it works; '3tap' is two 3x3 convolutions.
%
%   THE RIESZ IS COMPUTED ONLY WHERE IT IS ASKED FOR.  Three planes per level over
%   every frame is the memory in this block: on a 512x512 window of 1 300 frames
%   the finest level alone is 1.4 GB per plane.  Naming the levels that will be
%   amplified is how the continuous mode fits in memory at all.
%
% Syntax:
%    pyr = meRieszPyramid(I)
%    pyr = meRieszPyramid(I, nLevels)
%    pyr = meRieszPyramid(I, nLevels, rieszLevels)
%    pyr = meRieszPyramid(I, nLevels, rieszLevels, riesz)
%
% Inputs:
%    I           - one image or a stack, [rows columns] or [rows columns frames].
%                  Converted to single unless it is already single or a gpuArray.
%    nLevels     - how many Laplacian bands to build.  Empty or absent takes the
%                  most the frame allows, which leaves the residual at least four
%                  pixels across.  Level 1 is the finest and each level up doubles
%                  the wavelength it carries: level k is centred on 2^(k+1) pixels.
%    rieszLevels - the levels whose Riesz components are wanted.  Empty or absent
%                  takes all of them.
%    riesz       - 'exact' (default) or '3tap'.  See the header.
%
% Outputs:
%    pyr - .lap    {1 x nLevels} the Laplacian bands, each [rows columns frames]
%          .res    the residual lowpass, the coarsest thing that is left
%          .rx     {1 x nLevels} Riesz along columns; empty where not asked for
%          .ry     {1 x nLevels} Riesz along rows;    empty where not asked for
%          .size   {1 x nLevels} [rows columns] of each band, for the collapse
%          .resSize [rows columns] of the residual
%          .imSize [rows columns] of the input frame
%          .nFrames
%          .levels the levels that carry a Riesz transform
%          .class  the class the bands are held in
%
% Example:
%    pyr = meRieszPyramid(single(stack), 5, 3:5);
%    assert(max(abs(meCollapse(pyr)-stack),[],'all') < 1e-4);
%
% Dependencies: Image Processing Toolbox (imfilter); meExpand.
% See also: meCollapse, mePhase, meShift, meExpand
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function pyr = meRieszPyramid(I, nLevels, rieszLevels, riesz)

if ~isa(I,'single') && ~isa(I,'gpuArray'), I = single(I); end
imSize = [size(I,1) size(I,2)];

maxLevels = max(1, floor(log2(min(imSize))) - 2);
if nargin<2 || isempty(nLevels), nLevels = maxLevels; end
if nLevels > maxLevels
    error('meRieszPyramid:levels', ...
        'A %dx%d frame carries %d levels, not %d - the residual would be under four pixels.', ...
        imSize(1), imSize(2), maxLevels, nLevels);
end
if nargin<3 || isempty(rieszLevels), rieszLevels = 1:nLevels; end
rieszLevels = intersect(rieszLevels(:)', 1:nLevels);
if nargin<4 || isempty(riesz), riesz = 'exact'; end
riesz = lower(char(riesz));
if ~any(strcmp(riesz,{'3tap','exact'}))
    error('meRieszPyramid:riesz','The Riesz transform is 3tap or exact, not %s.',riesz);
end

% The two antisymmetric taps.  See the header: imfilter correlates, and with these
% kernels that is the true sign of the Riesz transform.  imfilter takes its kernel
% in double whatever the image is in, and returns the image's class.
kernelX = [0 0 0; 0.5 0 -0.5; 0 0 0];
kernelY = [0 0.5 0; 0 0 0; 0 -0.5 0];

pyr          = struct();
pyr.lap      = cell(1,nLevels);
pyr.rx       = cell(1,nLevels);
pyr.ry       = cell(1,nLevels);
pyr.size     = cell(1,nLevels);
pyr.imSize   = imSize;
pyr.nFrames  = size(I,3);
pyr.levels   = rieszLevels;
pyr.riesz    = riesz;
pyr.class    = class(I);

G = I;
for k = 1:nLevels
    sz          = [size(G,1) size(G,2)];
    coarse      = reduce(G);
    pyr.lap{k}  = G - meExpand(coarse, sz);
    pyr.size{k} = sz;
    if ismember(k, rieszLevels)
        if strcmp(riesz,'exact')
            [pyr.rx{k}, pyr.ry{k}] = exactRiesz(pyr.lap{k});
        else
            pyr.rx{k} = imfilter(pyr.lap{k}, kernelX, 'symmetric', 'same');
            pyr.ry{k} = imfilter(pyr.lap{k}, kernelY, 'symmetric', 'same');
        end
    end
    G = coarse;
end
pyr.res     = G;
pyr.resSize = [size(G,1) size(G,2)];
end

% =====================================================================
function B = reduce(G)
%reduce  Blur with the separable binomial and keep every other sample.  meExpand
%   uses twice this kernel, which is the factor the discarded samples cost.
k = [1 4 6 4 1]./16;
B = imfilter(G, k.', 'symmetric', 'same');
B = imfilter(B, k  , 'symmetric', 'same');
B = B(1:2:end, 1:2:end, :);
end

% =====================================================================
function [rx,ry] = exactRiesz(L)
%exactRiesz  The Riesz transform as its own multiplier, -i*w/|w|, with no
%   approximation - and therefore periodic at the band's edges rather than
%   symmetric.  The band is a Laplacian band, which is already near zero at the
%   frame boundary, so the wrap costs little; that is why this is affordable at all.
[nr,nc,~] = size(L);
wx = reshape(2*pi*ifftshift((0:nc-1)-floor(nc/2))./nc, 1, nc);
wy = reshape(2*pi*ifftshift((0:nr-1)-floor(nr/2))./nr, nr, 1);
W  = hypot(wx,wy);  W(W==0) = 1;
F  = fft2(L);
rx = real(ifft2(F.*cast(-1i.*wx./W,'like',F)));
ry = real(ifft2(F.*cast(-1i.*wy./W,'like',F)));
end
%------------- END OF CODE --------------
