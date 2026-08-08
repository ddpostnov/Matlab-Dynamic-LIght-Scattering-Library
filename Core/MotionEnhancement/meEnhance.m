%meEnhance - Run one mode of the magnifier over one recording's working copies.
%
%   THE ONE PLACE THE TWO MODES ARE WRITTEN DOWN, and everything they differ in is
%   here rather than spread across the engine.  meShift owns the ORDER the pipeline
%   runs in; this owns WHAT it is run on, and those are two different decisions.
%
%   'cine'        magnify one averaged cardiac cycle.  This is what amplified MRI
%                 does - Holdsworth 2016, Terem 2018 and 2021 all amplify a
%                 retrospectively gated cine, never a raw series - and on this
%                 recording it is not one of two options to compare.  Session 0
%                 measured the same walls both ways: band-passed at the heart rate
%                 they returned 1.08 to 1.25 times an empty band, which is no
%                 detection at all, and folded over the quiet beats they cleared
%                 their own floor by 4.6 to 7.6.
%   'continuous'  magnify real frames in real time.  The honest view: it keeps
%                 beat-to-beat variability, which the cine averages away by
%                 construction, and it is what a viewer would see through the
%                 microscope.  EXPECT IT TO LOOK FAR WORSE, and that is a
%                 measurement rather than a hedge.
%
%   THE CINE'S AXIS IS PHASE, NOT TIME, AND THAT CHANGES THE FILTER.  A cine is one
%   whole cardiac period sampled at nPhaseBins points, so its last bin is adjacent to
%   its first and the wrap-around is the physics.  The band is therefore in cycles per
%   cycle rather than hertz, the sampling rate is bins per cycle, and the filter is
%   'cyclic' - the FFT band-pass with no padding.  meTemporal's default 'ideal'
%   mirrors the series at both ends, which on a genuine period makes it EVEN about
%   both, and an even beat has no systolic upstroke.  Everything derived here is
%   returned in info.settings, so nothing is silently different from what was asked
%   for.
%
%   THE CONTINUOUS RUN STAYS INSIDE ONE UNBROKEN QUIET STRETCH BY DEFAULT, and on the
%   reference recording that is 1 302 frames - 13.0 s - not the 2 000 the plan sized
%   it at.  A bout does not merely look bad; it exceeds the bound the magnifier is
%   valid over and it drags the whole-field phase mean that the differential mode
%   removes, so it can spoil frames that are themselves quiet.  s.contSpanBouts is
%   there for the run whose point is to show what a bout does, and then every bad
%   frame is named in info.badInWindow rather than silently included.
%
%   THE OUTER MARGIN IS CROPPED, BECAUSE IT IS PADDING AND NOT DATA.  The brightest
%   pixels of every Laplacian band are the frame boundary, where the symmetric padding
%   puts a large response - this is measurable on any frame and it caught Session 1's
%   zero-phase test out.  s.edgeCrop pixels come off every side of both the magnified
%   stack and the original that is returned beside it, so a side-by-side comparison is
%   of two pictures of the same thing.
%
%   NOTHING HERE PICKS alpha.  Session 3 does, against Session 1's calibration and the
%   calibre of the vessels a run is meant for.  What this returns is the run at
%   whatever s.alpha says, plus every number needed to judge it.
%
% Syntax:
%    [out, info] = meEnhance(s, pass)
%
% Inputs:
%    s    - settings.  Fields read, beyond everything meShift and meGate read:
%           .mode      'cine' or 'continuous'
%           .band      'puls' or 'vaso' - which rhythm.  'vaso' moves the gate's
%                      search range onto s.vFR, so one gate serves both.
%           .copy      'bin' or 'crop' - which working copy to read.
%           .levels    pyramid levels OF THE RECORDING.  Mapped onto whichever copy
%                      is read, so one value is right for both.
%           .decimate  frames averaged into one before a continuous run.  Empty
%                      derives it from the passband, which is 1 for the cardiac band
%                      and about seventy for the vasomotion one.
%           .edgeCrop  pixels trimmed off every side of the result.
%           .maskMode  'weight' builds a wall mask from the run's own base frame
%                      when s.mask is empty.
%           .null      '' is the product, 'shuffle' is its matched control - the
%                      SAME run with the timing destroyed and nothing else changed.
%           .nullSeed  that control's own random stream.
%    pass - the sidecar struct from mePassA, or the path to its .mat file.
%
% Outputs:
%    out  - the magnified stack, single, [rows columns frames] after the crop.
%    info - .original  the unmagnified input, cropped identically
%           .mode .band .copy .null   what was run
%           .settings  the settings AS DERIVED, including the passband, the filter
%                      and the mask that were actually used
%           .epochs    from meEpochs
%           .gate      from meGate
%           .cine      from meCine, when the mode is 'cine'
%           .window    [first last] frames, when the mode is 'continuous'
%           .decimate  frames averaged into one, for that run
%           .badInWindow  which of the run's frames are inside a bout
%           .base      the picture the run is drawn on, cropped
%           .mask      the wall weight that was used, cropped, or empty
%           .fsRun     the sampling rate of the axis that was filtered
%           .passband  the band that was kept, in that axis's units
%           .levels    the levels of THE COPY that were amplified
%           .lambdaFull the wavelengths those carry, in pixels of the recording -
%                      which is the unit Session 0 measured the wall edge in
%           .seconds   .gate .cine .shift .total
%           .peakGB    the largest MATLAB array footprint seen during the run
%
% Example:
%    s = meSettings;  s.mode = 'cine';  s.band = 'puls';  s.copy = 'bin';
%    s.maskMode = 'weight';  s.alpha = 40;
%    [mag, info] = meEnhance(s, 'C:\Data\CXD\motionEnhancement\rec_ME_passA.mat');
%    fprintf('%d cycles, acceptance %.0f%%\n', info.cine.nCycles, 100*info.gate.acceptanceRate);
%
% Dependencies: meReadRaw, meEpochs, meGate, meCine, meMask, meShift.
% See also: meEpochs, meGate, meCine, meMask, meShift, meProbe, mePassA
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function [out, info] = meEnhance(s, pass)

tAll = tic;
if ischar(pass) || isstring(pass)
    pass = load(char(pass),'pass').pass;
end

probePath = fullfile(pass.outFolder, [pass.stem '_ME_probe.mat']);
if ~isfile(probePath)
    error('meEnhance:noProbe', ...
        ['meProbe has not been run on this recording, and its frame-to-frame ' ...
         'displacement series is what the epoch flag is built from: %s'], probePath);
end
probe = load(probePath,'probe').probe;
tr    = load(pass.tracePath);
fs    = pass.fs;

sd = s;                                   % the settings AS DERIVED, returned in info

% ---- when the animal was still ----------------------------------------------
ep = meEpochs(sd, probe.motion.displacement, fs);

% ---- which rhythm, and the gate for it --------------------------------------
switch lower(sd.band)
    case 'puls'
        % the search range stays where meSettings put it
    case 'vaso'
        sd.minFrqIni = sd.vFR(1);
        sd.maxFrqIni = sd.vFR(2);
    otherwise
        error('meEnhance:band','The band is puls or vaso, not %s.', sd.band);
end
tG   = tic;
gate = meGate(sd, tr.traceVes, fs, ep.bad);
tGate = toc(tG);

% ---- the copy the run reads -------------------------------------------------
[m, geom] = meReadRaw(pass, sd.copy);

% s.levels NAMES A PHYSICAL SCALE, AND THE COPIES DO NOT SHARE ONE.  A pyramid level
% carries a wavelength in pixels OF THE ARRAY IT IS BUILT ON, so level 3 of the
% two-by-two binned copy is 16 binned pixels, which is 32 pixels of the recording -
% twice the scale Session 0 chose it for.  Getting this wrong is not an error anyone
% would see: it amplifies wavelengths above the structure and simply returns less.
% Measured on the reference cine at alpha 20, levels 3:5 taken literally on the binned
% copy delivered 2.77 times the raw wall modulation where the same physical band,
% which is 2:4 there, delivered 5.94.  So s.levels is stated in levels of the
% RECORDING and mapped onto whichever copy is being read.
lvlShift  = round(log2(geom.scale));
sd.levels = s.levels(:)' - lvlShift;
if any(sd.levels < 1)
    error('meEnhance:levels', ...
        ['Levels %s of the recording are levels %s of the %s copy, and there is no ' ...
         'level below 1. Bin less, or name coarser levels.'], ...
        mat2str(s.levels(:)'), mat2str(sd.levels), sd.copy);
end
lambdaFull = 2.^(s.levels(:)'+1);

% ---- what the magnifier is run on -------------------------------------------
tC = tic;  cn = [];  win = [];
switch lower(sd.mode)

    case 'cine'
        cn    = meCine(sd, m, geom.fs, gate);
        stack = cn.cine;
        base  = cn.base;
        fsRun = cn.fsBins;
        sd.passband       = cn.passband;
        sd.temporalFilter = 'cyclic';

    case 'continuous'
        win         = ep.window;
        sd.passband = gate.band;

        % A BAND SETS THE SAMPLING RATE IT NEEDS, AND A SLOW BAND NEEDS A LONG
        % RECORD.  Two cardiac cycles fit in a fifth of a second, so the cardiac run
        % takes the frames as they come and s.decimate derives to 1.  A vasomotion
        % cycle is five or six seconds, so a run that shows any of them has to span
        % minutes - and minutes of full frames do not fit in memory at any resolution
        % this block writes.  Averaging blocks of frames solves both at once: it puts
        % the whole recording inside a few hundred frames AND it improves the photon
        % count per frame by the root of the block, which is the currency the slow
        % band is short of.  A BLOCK MEAN, not a sample every nth: point-sampling
        % folds everything above the new Nyquist - which here is the entire cardiac
        % band - straight down into the vasomotion band.
        dec = sd.decimate;
        if isempty(dec)
            dec = max(1, floor(geom.fs/(6*sd.passband(2))));
        end
        nOut = max(1, floor((win(2)-win(1)+1)/dec));
        if dec == 1
            stack = single(m.Data.d(:,:,win(1):win(1)+nOut-1));
        else
            stack = zeros(geom.size(1), geom.size(2), nOut, 'single');
            bad   = false(nOut,1);
            for j = 1:nOut
                idx = win(1) + (j-1)*dec + (0:dec-1);
                stack(:,:,j) = mean(single(m.Data.d(:,:,idx)),3);
                bad(j)       = any(ep.bad(idx));
            end
            ep.badInWindow = bad;         % a decimated frame is bad if any of its is
        end
        % THE MATCHED CONTROL: the frames the magnifier sees, permuted in time.
        % Real noise, real texture, real bulk motion, no rhythm - and the bout flag
        % travels with its own frame so the video still marks the right ones.
        if strcmpi(sd.null,'shuffle')
            rs  = RandStream('threefry','Seed',sd.nullSeed);
            ord = randperm(rs, size(stack,3));
            stack = stack(:,:,ord);
            if numel(ep.badInWindow) == numel(ord)
                ep.badInWindow = ep.badInWindow(ord);
            end
        end

        base  = mean(stack,3);
        fsRun = geom.fs/dec;
        sd.decimate = dec;

    otherwise
        error('meEnhance:mode','The mode is cine or continuous, not %s.', sd.mode);
end
tBuild = toc(tC);

% ---- where amplification is allowed -----------------------------------------
% Built on the run's OWN base frame, so it is already on the grid of the copy being
% magnified and nothing has to be resampled between building it and using it.
if strcmpi(sd.maskMode,'weight') && isempty(sd.mask)
    if strcmpi(sd.maskSource,'segmentation') && isempty(sd.maskVessel)
        mn = load(pass.meanPath,'vesselMask');
        sd.maskVessel = toCopyGrid(mn.vesselMask, geom) > 0.5;
    end
    sd.mask = meMask(sd, base);
end

% ---- magnify ----------------------------------------------------------------
tS = tic;
[mag, an] = meShift(stack, fsRun, sd);
tShift = toc(tS);
peakGB = analysisGB(an) + numel(stack)*4/2^30;
clear an

% ---- the outer margin is padding --------------------------------------------
c   = max(round(sd.edgeCrop),0);
out = cropMargin(mag, c);

info = struct();
info.original    = cropMargin(stack, c);
info.mode        = sd.mode;
info.band        = sd.band;
info.copy        = sd.copy;
info.null        = sd.null;
info.settings    = sd;
info.epochs      = ep;
info.gate        = gate;
info.cine        = cn;
info.window      = win;
info.decimate    = sd.decimate;
info.badInWindow = ep.badInWindow;
info.base        = cropMargin(base, c);
info.mask        = cropMargin(sd.mask, c);
info.geom        = geom;
info.fsRun       = fsRun;
info.passband    = sd.passband;
info.levels      = sd.levels;
info.lambdaFull  = lambdaFull;
info.edgeCrop    = c;
info.peakGB      = peakGB;
info.seconds     = struct('gate',tGate,'build',tBuild,'shift',tShift,'total',toc(tAll));
end

% =====================================================================
function Y = cropMargin(X, c)
%cropMargin  Take c pixels off every side, of a stack, an image, or nothing at all.
if isempty(X) || c<=0, Y = X; return; end
if size(X,1) <= 2*c || size(X,2) <= 2*c
    error('meEnhance:crop','An edge crop of %d px leaves nothing of a %dx%d frame.', ...
        c, size(X,1), size(X,2));
end
Y = X(c+1:end-c, c+1:end-c, :);
end

% =====================================================================
function gb = analysisGB(an)
%analysisGB  What the analysis is holding, in gigabytes.
%   Not a memory profiler.  The number that decides whether a mode fits is the size
%   of the arrays this block deliberately keeps over every frame at once - the
%   pyramid, its Riesz components and the phase fields - and those are countable.
n = 0;
for k = 1:numel(an.pyr.lap)
    n = n + numel(an.pyr.lap{k}) + numel(an.pyr.rx{k}) + numel(an.pyr.ry{k});
end
for k = an.levels
    n = n + numel(an.pc{k}) + numel(an.ps{k}) + numel(an.amp{k});
end
gb = 4*n/2^30;                 % single throughout
end

% =====================================================================
function Y = toCopyGrid(full, geom)
%toCopyGrid  A full-resolution map, on the grid of one working copy.
%   The copies do not share a coordinate system - meReadRaw's geom carries the map -
%   and a mask built on the wrong grid is a mask of somewhere else.
r0 = geom.origin(1);  nr = geom.size(1);
c0 = geom.origin(2);  nc = geom.size(2);
sc = geom.scale;
sub = double(full(r0:r0+nr*sc-1, c0:c0+nc*sc-1));
if sc > 1
    Y = imresize(sub, [nr nc], 'box');
else
    Y = sub;
end
end
%------------- END OF CODE --------------
