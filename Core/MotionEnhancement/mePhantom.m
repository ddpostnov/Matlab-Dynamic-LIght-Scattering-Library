%mePhantom - A recording whose wall motion is known exactly, to calibrate against.
%
%   THE ONLY DATA IN THIS PROJECT WHOSE TRUTH IS NOT ITSELF A MEASUREMENT.  Every
%   number the block reports about a real recording is an estimate compared with
%   another estimate; here the displacement is a constant in a formula, so the
%   calibration curve has something to be a curve AGAINST.
%
%   BUILT IN CLOSED FORM, NOT BY imwarp, AND THAT IS DELIBERATE.  The plan said
%   warp a frame with a smooth displacement field.  The displacements this recording
%   actually contains are 0.006 to 0.075 px, and at that size a resampling
%   interpolator's own error is not negligible against the signal - the phantom would
%   be calibrating imwarp as much as the magnifier.  A vessel cross-section is a rect
%   convolved with a Gaussian, which has a closed form; evaluating it at a moved wall
%   position is exact at any amplitude, costs nothing, and leaves the truth a
%   constant rather than a thing to be trusted.  Everything that moves - the wall,
%   the texture riding on it, the tissue behind it - is a function of the moved
%   coordinate, so the deformation is smooth and total.
%
%   TWO MOTIONS, BECAUSE TWO DIFFERENT THINGS HAVE TO BE ASSERTED:
%     'dilate'    both walls move outward together, r(t) = r0 + A*sin(2*pi*f*t).
%                 This is pulsatility.  The signal is the half-width; the centre is
%                 still.  It is what the calibration curve is measured on.
%     'translate' the whole field slides, u -> u - A*sin(2*pi*f*t).  This is bulk
%                 preparation motion.  The signal is the centre; the half-width is
%                 still.  It is what the DIRECTION assertion is measured on - and
%                 what meGlobal must remove completely.
%
%   THE NOISE IS SHOT NOISE, MATCHED TO THE RECORDING.  Session 0 measured a still
%   patch at mean 1 072 counts and temporal sigma 46.0, i.e. 1.97 counts per photon
%   and 544 photons - SNR 23.3 per pixel - and confirmed it is independent pixel to
%   pixel by binning 2x2 and getting a gain of 2.01 against a theoretical 2.  So the
%   image is built in photons and drawn from a Poisson distribution.  Every phase
%   argument in this block rests on that independence; putting anything else here
%   would make the phantom easier than the data.
%
%   THE SEED DOES NOT DEPEND ON THE AMPLITUDE.  Every phantom in a sweep therefore
%   carries the SAME noise realisation, so a difference between two points of the
%   calibration curve is the amplitude and nothing else - and the no-motion null is
%   the same recording with A set to zero rather than a different draw of the dice.
%
% Syntax:
%    [stack, truth] = mePhantom(s)
%
% Inputs:
%    s - settings.  Fields read:
%        .phSize       [rows columns] of the frame
%        .phFrames     how many frames
%        .phFs         frames per second
%        .phFreq       the motion frequency, Hz.  Snapped to a whole number of
%                      cycles in the record, so the amplitude can be read off one
%                      Fourier coefficient with no leakage at all.
%        .phAmp        the displacement amplitude, PIXELS PER WALL
%        .phMode       'dilate' or 'translate'
%        .phRadius     the vessel's half-width at rest, pixels
%        .phEdge       sigma of the wall's blur, pixels.  The 10 to 90 per cent
%                      edge width that comes out is 2.563 times this.
%        .phPhotons    photons at the middle of the lumen
%        .phBackground tissue brightness as a fraction of the lumen
%        .phTexture    texture on the wall, as a fraction of the lumen
%        .phTissue     texture in the tissue, as a fraction of the lumen
%        .phSeed       the RNG seed
%
% Outputs:
%    stack - [rows columns frames] single, in photons.
%    truth - .amplitude    the displacement per wall, pixels - EXACT
%            .freq         the snapped motion frequency, Hz
%            .nCycles      whole cycles in the record
%            .fs .nFrames .mode
%            .measure      'halfWidth' for dilate, 'centre' for translate - which
%                          of meKymograph's two outputs carries the motion
%            .cut          the cut struct for meKymograph, across the vessel
%            .band         a passband bracketing the motion, Hz
%            .radius .edgeSigma .edgeWidth1090 .photons .snr .seed
%
% Example:
%    s = meSettings;  s.phAmp = 0.05;
%    [stack, truth] = mePhantom(s);
%    kymo = meKymograph(stack, truth.cut);
%
% Dependencies: Statistics and Machine Learning Toolbox (poissrnd).
% See also: meValidate, meKymograph, meShift, meSettings
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function [stack, truth] = mePhantom(s)

nRow = s.phSize(1);  nCol = s.phSize(2);
nT   = s.phFrames;   fs   = s.phFs;
r0   = s.phRadius;   sig  = s.phEdge;

% A whole number of cycles in the record.  One Fourier coefficient then carries the
% whole motion with no leakage, which is what lets an amplitude of two thousandths
% of a pixel be read off at all.
nCyc = max(1, round(s.phFreq*nT/fs));
f    = nCyc*fs/nT;
t    = (0:nT-1)./fs;
d    = s.phAmp.*sin(2*pi*f.*t);

st  = rng(s.phSeed);
cl  = onCleanup(@() rng(st));

% ---- the fixed random parts, drawn once ------------------------------------
nTex   = 6;
lamR   = 6 + 20*rand(1,nTex);   phiR = 2*pi*rand(1,nTex);   cR = randn(1,nTex);
lamL   = 6 + 20*rand(1,nTex);   phiL = 2*pi*rand(1,nTex);   cL = randn(1,nTex);
cR     = cR./sum(abs(cR));      cL   = cL./sum(abs(cL));

nBg    = 10;
lamB   = 10 + 50*rand(1,nBg);   thB  = pi*rand(1,nBg);      psiB = 2*pi*rand(1,nBg);
cB     = randn(1,nBg);          cB   = cB./sum(abs(cB));
kxB    = 2*pi*cos(thB)./lamB;   kyB  = 2*pi*sin(thB)./lamB;

% ---- geometry ---------------------------------------------------------------
col  = (1:nCol);
row  = (1:nRow)';
c0   = (nCol+1)/2;
u0   = col - c0;

% A slow brightness variation along the vessel, so the frame is not one profile
% repeated.  It multiplies the whole cross-section, so it moves no edge.
M = 1 + 0.15*sin(2*pi*1.5.*row./nRow + 0.7);

% Tissue structure.  When only the walls move it stands still, so it is built once;
% when the whole field slides it has to follow, and only then is it per frame.
moves = strcmpi(s.phMode,'translate');
BG    = tissue(col, row, 0, cB, kxB, kyB, psiB);

stack = zeros(nRow, nCol, nT, 'single');
for k = 1:nT
    if moves
        uu = u0 - d(k);   r = r0;
        BG = tissue(col, row, d(k), cB, kxB, kyB, psiB);
    else
        uu = u0;          r = r0 + d(k);
    end

    % The vessel: a rect of half-width r convolved with a Gaussian of sigma sig.
    P = 0.5*(erf((r+uu)./(sqrt(2)*sig)) + erf((r-uu)./(sqrt(2)*sig)));

    % Texture pinned to the wall, so it travels with it.  v is the distance outside
    % the wall; the envelope keeps it where a wall is.
    v   = abs(uu) - r;
    tR  = sum(cR.'.*cos(2*pi.*v./lamR.' + phiR.'), 1);
    tL  = sum(cL.'.*cos(2*pi.*v./lamL.' + phiL.'), 1);
    tex = ((uu>=0).*tR + (uu<0).*tL).*exp(-(v./(2*sig)).^2);

    cross = s.phBackground + (1-s.phBackground).*P + s.phTexture.*tex;

    clean = M.*cross + s.phTissue.*(1-P).*BG;
    stack(:,:,k) = single(max(clean,0.01).*s.phPhotons);
end

stack = single(poissrnd(double(stack)));

truth = struct();
truth.amplitude     = s.phAmp;
truth.freq          = f;
truth.nCycles       = nCyc;
truth.fs            = fs;
truth.nFrames       = nT;
truth.mode          = lower(s.phMode);
truth.measure       = 'halfWidth';
if strcmpi(s.phMode,'translate'), truth.measure = 'centre'; end
truth.cut           = struct('center',[(nRow+1)/2 c0],'normal',[0 1],'radius',r0);
truth.passband          = f.*[1-s.rangeFrq, 1+s.rangeFrq];
truth.radius        = r0;
truth.edgeSigma     = sig;
truth.edgeWidth1090 = 2*sqrt(2)*erfinv(0.8)*sig;
truth.photons       = s.phPhotons;
truth.snr           = sqrt(s.phPhotons);
truth.seed          = s.phSeed;
end

% =====================================================================
function BG = tissue(col, row, shift, c, kx, ky, psi)
%tissue  A few plane waves, so the field is not a single profile repeated and so
%   there is something outside the vessel for the bulk-motion correction to see.
%   Analytic, so a translation of it is exact rather than interpolated.
BG = sum(reshape(c,1,1,[]).*cos(reshape(kx,1,1,[]).*(col-shift) + ...
                                reshape(ky,1,1,[]).*row + ...
                                reshape(psi,1,1,[])), 3);
end
%------------- END OF CODE --------------
