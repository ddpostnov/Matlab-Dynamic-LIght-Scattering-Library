%meSettings - Default settings for the motion-enhancement block.
%
%   One place where every default of the me* block is written down, so a run is
%   described by the difference between this struct and the one it was given.
%   Every me* function documents the fields it reads; none of them invents a
%   default of its own.
%
%   THE BLOCK IS STAND-ALONE.  It may call the library; the library never calls
%   it.  Nothing here is registered with a launcher, a GUI, the workbench or the
%   explorer, and nothing it writes carries a _r.mat or _s.mat tail - those are
%   swept by the workbench and this is not a library product.
%
% Syntax:
%    s = meSettings()
%
% Inputs:
%    None.
%
% Outputs:
%    s - the settings struct.  Fields, grouped by the function that reads them:
%
%        where things go
%          .outFolder    ''  - results go into a motionEnhancement folder beside
%                              the recording.  Name a folder to put the working
%                              copies somewhere else - on a machine where the
%                              recording sits in a synchronised folder that is
%                              worth doing, the copies are tens of gigabytes.
%          .reportLevel  'normal' or 'quiet'
%
%        the recording
%          .pixelSize    2.5 - micrometres per pixel at full resolution.
%
%        mePassA - the one streaming read of the recording
%          .binFactor    2   - the binned working copy is this much smaller.
%          .cropSize     [512 512] - size of the full-resolution working copy,
%                              [height width] in pixels.
%          .roi          []  - [firstRow lastRow; firstCol lastCol] for that
%                              copy.  Empty picks the window carrying the most
%                              pulsatile signal, which is what the probe wants.
%          .pctClip      [0.1 99.9] - intensity percentiles mapped to 0 and 255.
%          .sampleFrames 200 - frames spread over the recording for those limits.
%          .blockFrames  500 - consecutive frames used to measure the noise and
%                              to place the crop.
%          .blockStart   0.2 - where that block starts, as a fraction of the
%                              recording.
%          .bleachBlock  500 - frames averaged per point of the bleaching curve.
%          .writeBlock   250 - frames buffered before each write.
%          .reportEvery  1000 - frames between progress lines.
%          .vesselPct    85  - intensity percentile above which a pixel counts as
%                              vessel, for the second intensity trace.
%
%        meProbe - the report
%          .fftN         8192 - window length of the averaged spectrum.
%          .minFrqIni    1   - lowest frequency searched for the heart rate, Hz.
%          .maxFrqIni    20  - highest frequency searched for the heart rate, Hz.
%                              meEnhance moves this pair onto the vasomotion band
%                              when s.band is 'vaso', so one gate serves both.
%          .rangeFrq     0.3 - the heart-rate band is the found rate plus and
%                              minus this fraction of it.
%          .vFR          [0.05 0.25] - the vasomotion band, Hz.
%          .nCuts        5   - vessel cuts measured across the field.
%          .cutSep       80  - pixels between one cut and the next.
%          .kymoSeconds  5   - seconds of kymograph drawn on the report page.
%          .traceSeconds 30  - seconds of wall position drawn on the page and
%                              band-passed for the amplitude.
%          .nPhaseBins   25  - phase bins one cycle is averaged into, by meCine and
%                              by meProbe's own fold.
%          .motionWin    256 - side of the binned window the frame-to-frame
%                              displacement is measured on.
%          .quietCoef    4   - a frame belongs to a bout when it moves this many
%                              times the recording's own median step.
%
%        meEpochs - when the animal was still
%          .boutPad      0.25 - seconds the bout flag is grown EITHER SIDE of a
%                              threshold crossing.  The displacement crosses in the
%                              MIDDLE of a movement, so without this the shoulders
%                              of every bout stay unflagged.  It is worth a fifth of
%                              the recording: 0.25 s leaves 73.4 % of the reference
%                              recording quiet and 0.5 s leaves 53.0 %.
%          .contFrames   []  - longest window the continuous mode may take.  Empty
%                              takes the whole quiet stretch, which on the reference
%                              recording is 1 302 frames rather than the 2 000 the
%                              plan sized it at.  Inf takes the whole recording,
%                              which is what a vasomotion run needs and what
%                              .decimate makes affordable.
%          .contSpanBouts false - keep the continuous window inside one unbroken
%                              quiet stretch.  True lets it reach contFrames across
%                              a bout, with the bad frames marked for the video.
%
%        meGate - the cardiac gate, recovered from the recording itself
%          .gateDetrend  3   - order of the polynomial removed from the trace before
%                              its spectrum is taken.  The dye bleaches, and a
%                              monotone decay puts power into the search band.
%          .minPromCoef  1/4 - minimum peak prominence, as a fraction of the RHYTHM'S
%                              own standard deviation - not the whole trace's, which
%                              on a fluorescence recording is bleaching, respiration
%                              and vasomotion and comes out several times the pulse.
%          .coeffsSTD    [3 2 2 2 2 3 3 2 2] - one per feature: leave-one-out
%                              standard deviations from the leave-one-out median a
%                              cycle's feature may sit before it is rejected.
%          .coeffsRel    [0.5 0.1] - largest min-to-min jump as a share of the
%                              cycle's excursion, and largest departure of its level
%                              from the recording's median.
%          .coeffsAbs    2   - half-width, in cycles, that a rejection spreads over.
%          .excludeFirstNCycles 0 - rejected outright.
%          THE ACCEPTANCE RATE IS REPORTED, NOT TUNED.  An awake animal fails more
%          cycles than the anaesthetised traces these were set on, and that is
%          information about the preparation.
%
%        the magnifier - meShift and the functions it calls
%          .levels       3:5 - pyramid levels that are amplified, IN PIXELS OF THE
%                              RECORDING.  Level k carries a wavelength of 2^(k+1)
%                              pixels, so this is 16, 32 and 64.  NOT 2:5: Session 0
%                              measured the wall's edge at 6 to 15 px, so level 2
%                              sits below the structure it would be amplifying and
%                              carries mostly noise.
%                              THE COPIES DO NOT SHARE A SCALE, so meEnhance maps
%                              these onto whichever copy it reads - 3:5 becomes 2:4
%                              on the two-by-two binned one.  Taken literally there
%                              it would amplify 32 to 128 recording pixels, above
%                              the structure, and return 2.1x less.  meShift takes
%                              the levels of the array it is handed, so a caller
%                              going straight to it names the copy's own.
%          .nPyrLevels   []  - how many pyramid bands to build.  Empty builds as
%                              many as the highest amplified level needs.
%          .riesz        'exact' - the Riesz multiplier itself, two Fourier
%                              transforms per level per frame.  '3tap' is the
%                              published three-tap approximation.  Measured: the
%                              exact one amplifies 24 per cent more at a detection
%                              floor that cannot be told apart, for 17 per cent
%                              more time.  Try '3tap' if a real product shows
%                              structure mirrored across the frame - the exact
%                              transform treats the band as periodic.
%          .alpha        20  - the amplification factor.  One number, or one per
%                              amplified level - and one per level is worth 1.5x,
%                              because each level's bound scales with its own
%                              wavelength.  What comes out has moved by 1+eta*alpha
%                              times the input, where eta is 0.39 on levels 3:5 and
%                              is measured by meValidate, NOT by 1+alpha: the
%                              pyramid's bands overlap, so a level that is not
%                              amplified holds the edge back.
%          .temporalFilter 'ideal' - 'ideal', 'butter', 'acceleration' or 'cyclic'.
%                              All four are zero-phase; nothing causal belongs here.
%                              'cyclic' is for a series that IS one whole period -
%                              the cine - where the wrap-around is the physics and
%                              padding would be the artefact.
%          .passband     []  - [low high] kept by that filter, in the units of the
%                              series' own axis: hertz for a continuous run, cycles
%                              per cardiac cycle for a cine.  Empty means the caller
%                              has not chosen yet - meShift needs it set, and
%                              meEnhance derives it from s.band.
%          .filterOrder  2   - half the IIR order used by 'butter'.
%          .globalOrder  0   - -1 leaves bulk motion in, which is the absolute
%                              mode; 0 removes a rigid translation, which is the
%                              differential mode and the default; 1 also removes
%                              scale, rotation and parallax.
%          .blurSigma    2   - width of the amplitude-weighted blur of the phase,
%                              in pixels OF EACH LEVEL.
%
%        meMask - where amplification is allowed
%          .maskMode     'none' - 'none' amplifies everywhere; 'weight' uses s.mask.
%          .mask         []  - the weight in [0,1], the size of the frame being
%                              magnified.  meEnhance builds it with meMask when it
%                              is empty and the mode asks for one.
%          .maskSource   'gradient' - how meMask builds that weight.  'gradient'
%                              needs nothing but the time-mean frame; 'segmentation'
%                              calls the library's categorisation and takes the ring
%                              around what it finds.
%          .maskSmooth   2   - pixels of Gaussian smoothing, before the gradient and
%                              again on the weight.  A hard-edged mask rings.
%          .maskThresh   0.25 - fraction of the gradient's robust maximum below
%                              which there is no wall.  It is a DECISION boundary,
%                              not a scaling: what comes out of it is one inside the
%                              mask and zero outside, softened only at the border.
%                              Grading the weight by the gradient magnifies a
%                              strong-edged stretch of wall more than a weak-edged
%                              one, so a wall moving uniformly comes out moving
%                              unevenly.
%          .maskRing     3   - half-width in pixels of the ring 'segmentation'
%                              builds, dilate minus erode.
%          .maskVessel   []  - the logical vessel mask that ring is built around.
%                              mePassA writes one into <rec>_ME_mean.mat and
%                              meEnhance passes it; a caller who has run the
%                              library's own segmentation over the same field can
%                              pass its lumen mask instead and get a better ring.
%          .maskExclude  []  - a logical map, the size of the frame, of places NOT
%                              to amplify.  For excluding one named vessel and using
%                              it as the control instead.
%
%        the two modes - meEnhance
%          .mode         'cine' - 'cine' magnifies one averaged cycle, which is what
%                              amplified MRI does and, on this recording, the only
%                              thing that works.  'continuous' magnifies real frames
%                              in real time: the honest view, keeping beat-to-beat
%                              variability, and expected to look far worse.
%          .band         'puls' - which rhythm.  'puls' takes the gate's own
%                              measured cardiac band, 'vaso' takes s.vFR.  This
%                              names the rhythm; .passband carries the numbers.
%          .copy         'bin' - which working copy to read.  'bin' is the whole
%                              field two-by-two binned, 'crop' the full-resolution
%                              512x512 window.
%          .cineHarmonics 5  - harmonics of the cycle the cine's filter keeps.  Ten
%                              frames a beat samples the fundamental and four
%                              harmonics, so this is what the data support.
%          .decimate     []  - frames averaged into one before a CONTINUOUS run.
%                              Empty derives it from the passband at six samples per
%                              cycle of its fastest component: 1 for the cardiac
%                              band, about seventy for the vasomotion one.  A BLOCK
%                              MEAN, not every nth frame - point-sampling would fold
%                              the whole cardiac band down into the vasomotion one -
%                              and it buys the root of the block in photon count,
%                              which is the currency a slow band is short of.
%          .edgeCrop     16  - pixels trimmed from every side of a magnified frame.
%                              The brightest pixels of every Laplacian band are the
%                              FRAME BOUNDARY, where the padding puts a large
%                              response, so the outer margin is padding and not data.
%                              MEASURED, not guessed: against the same ring of the
%                              unmagnified cine, the magnification runs 4.7 times the
%                              interior at one pixel in, 3.1 at four and 1.0 to 1.6
%                              at sixteen, on both level sets tried.
%          .null         ''  - '' is the product; 'shuffle' is its matched control.
%                              The control destroys the TIMING and nothing else: the
%                              continuous modes permute the frames the magnifier
%                              sees, and the cine gives every accepted cycle its own
%                              random circular time shift before it is folded, so
%                              every bin still averages one interpolated sample of
%                              every cycle and only the phase alignment is gone.
%                              Same frames, same texture, same photon noise, same
%                              alpha - and no rhythm.  Whatever moves in the control
%                              is what the pipeline invents.
%          .nullSeed  20260807 - fixed, and drawn from a private stream, so a control
%                              is reproducible and running one does not move the
%                              random numbers any other part of a session draws.
%
%        meWriteVideo - the deliverable
%          .videoFps     25  - frames per second written.  A cardiac cine played at
%                              the plan's "2x slowed" is 141 fps, which no player
%                              shows; the factor is burnt onto the frame instead, so
%                              a viewer reads the slowdown rather than assuming one.
%          .videoLoops   1   - times the stack is repeated.  A cine IS one cardiac
%                              period, so it is meant to be looped.
%          .videoCodec   'Motion JPEG AVI' - the profile, and it applies to a product
%                              AND its null.  Two codecs are not the identical
%                              treatment the whole comparison rests on.
%          .videoQuality 100
%          .showDifference false - add a third panel of the magnified minus the
%                              original, which is where the magnification actually
%                              is and the fastest way to tell a product from its null.
%          .diffGain     10  - what that difference is multiplied by.
%          .videoPct     [0.5 99.5] - percentiles of the ORIGINAL mapped to black and
%                              white, and the magnified panel is drawn through the
%                              same pair.
%          .videoFontSize 14
%          .scaleBarUm   200 - length of the burnt-in scale bar, micrometres.
%
%        mePhantom - the sequence whose motion is known exactly
%          .phSize       [224 224] - frame size, [rows columns].
%          .phFrames     512 - frames.
%          .phFs         100 - frames per second, near the recording's 99.96.
%          .phFreq       11.32 - the motion frequency, Hz, the measured heart
%                              rate.  Snapped to a whole number of cycles.
%          .phAmp        0.05 - displacement amplitude PER WALL, pixels.  Session
%                              0 measured 0.006 to 0.075 on five real vessels.
%          .phMode       'dilate' - 'dilate' moves the two walls apart and
%                              together, which is pulsatility; 'translate' slides
%                              the whole field, which is bulk motion.
%          .phRadius     15  - the vessel's half-width at rest, pixels, so a
%                              75 micrometre vessel at 2.5 micrometres a pixel.
%          .phEdge       3.9 - sigma of the wall blur, giving a 10 to 90 per cent
%                              edge width of 10 px, inside the 6 to 15 measured.
%          .phPhotons    544 - photons at the middle of the lumen, which is the
%                              measured SNR of 23.3 per pixel.
%          .phBackground 0.25 - tissue brightness as a fraction of the lumen.
%          .phTexture    0.06 - texture riding on the wall, same units.
%          .phTissue     0.08 - texture in the tissue, same units.
%          .phSeed       7919 - fixed, and NOT varied with the amplitude, so two
%                              points of a calibration curve differ by the motion
%                              and by nothing else.
%
%        meValidate - the calibration sweep
%          .calAmp       the displacement amplitudes swept, pixels.  Bracketed on
%                        0.002 to 0.5 because that is where the data live; the
%                        plan's 0.02 to 2 would calibrate a regime this recording
%                        never visits.
%          .calAlpha     the amplification factors swept.
%          .calSnrGain   multipliers on the photon count, i.e. how many frames were
%                        averaged into each frame the magnifier sees.  1 is one raw
%                        frame; 969 is what folding Session 0's 2 202 quiet beats
%                        into 20 phase bins buys per bin.  The detection floor is
%                        reported against this axis rather than as one number,
%                        because that is the axis the whole package turns on.
%          .calCleanGain 400 - the photon multiplier the GAIN grid runs at.  The
%                        gain is a property of the representation, so it is
%                        measured where noise is not competing with it; the floor
%                        is measured separately at the recording's own SNR.
%          .calFloorAmp  0.05 px - the displacement the floor sweep is measured at.
%                        The response is linear in the displacement, so the floor
%                        scales straight off one well-measured point.
%          .calBlur      the phase-blur widths swept.  The blur and the level
%                        choice are both bought with amplification, so the sweep
%                        reports gain and floor together and neither is a decision
%                        on its own.
%
% Example:
%    s = meSettings;
%    s.cropSize = [640 640];
%    mePassA(s,'C:\Data\CXD\EPFL_20241114_2ADTF13BP.cxd');
%
% Dependencies: none.
% See also: mePassA, meProbe, meReadRaw, meShift, mePhantom, meValidate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function s = meSettings()

s = struct();

% ---- where things go -------------------------------------------------------
s.outFolder    = '';
s.reportLevel  = 'normal';

% ---- the recording ---------------------------------------------------------
s.pixelSize    = 2.5;

% ---- mePassA ---------------------------------------------------------------
s.binFactor    = 2;
s.cropSize     = [512 512];
s.roi          = [];
s.pctClip      = [0.1 99.9];
s.sampleFrames = 200;
s.blockFrames  = 500;
s.blockStart   = 0.2;
s.bleachBlock  = 500;
s.writeBlock   = 250;
s.reportEvery  = 1000;
s.vesselPct    = 85;

% ---- meProbe ---------------------------------------------------------------
s.fftN         = 8192;
s.minFrqIni    = 1;
s.maxFrqIni    = 20;
s.rangeFrq     = 0.3;
s.vFR          = [0.05 0.25];
s.nCuts        = 5;
s.cutSep       = 80;
s.kymoSeconds  = 5;
s.traceSeconds = 30;
s.nPhaseBins   = 25;
s.motionWin    = 256;
s.quietCoef    = 4;

% ---- meEpochs --------------------------------------------------------------
s.boutPad       = 0.25;
s.contFrames    = [];
s.contSpanBouts = false;

% ---- meGate ----------------------------------------------------------------
s.gateDetrend         = 3;
s.minPromCoef         = 1/4;
s.coeffsSTD           = [3 2 2 2 2 3 3 2 2];
s.coeffsRel           = [0.5 0.1];
s.coeffsAbs           = 2;
s.excludeFirstNCycles = 0;

% ---- the magnifier ---------------------------------------------------------
s.levels         = 3:5;
s.nPyrLevels     = [];
s.riesz          = 'exact';
s.alpha          = 20;
s.temporalFilter = 'ideal';
s.passband       = [];
s.filterOrder    = 2;
s.globalOrder    = 0;
s.blurSigma      = 2;

% ---- meMask ----------------------------------------------------------------
s.maskMode    = 'none';
s.mask        = [];
s.maskSource  = 'gradient';
s.maskSmooth  = 2;
s.maskThresh  = 0.25;
s.maskRing    = 3;
s.maskVessel  = [];
s.maskExclude = [];

% ---- the two modes ---------------------------------------------------------
s.mode          = 'cine';
s.band          = 'puls';
s.copy          = 'bin';
s.cineHarmonics = 5;
s.decimate      = [];
s.edgeCrop      = 16;
s.null          = '';
s.nullSeed      = 20260807;

% ---- meWriteVideo ----------------------------------------------------------
s.videoFps       = 25;
s.videoLoops     = 1;
s.videoCodec     = 'Motion JPEG AVI';
s.videoQuality   = 100;
s.showDifference = false;
s.diffGain       = 10;
s.videoPct       = [0.5 99.5];
s.videoFontSize  = 14;
s.scaleBarUm     = 200;

% ---- mePhantom -------------------------------------------------------------
s.phSize       = [224 224];
s.phFrames     = 512;
s.phFs         = 100;
s.phFreq       = 11.32;
s.phAmp        = 0.05;
s.phMode       = 'dilate';
s.phRadius     = 15;
s.phEdge       = 3.9;
s.phPhotons    = 544;
s.phBackground = 0.25;
s.phTexture    = 0.06;
s.phTissue     = 0.08;
s.phSeed       = 7919;

% ---- meValidate ------------------------------------------------------------
s.calAmp        = [0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];
s.calAlpha      = [0 2 5 10 20 50 100 200];
s.calSnrGain    = [1 10 100 969];
s.calCleanGain  = 400;
s.calFloorAmp   = 0.05;
s.calBlur       = [0 1 2 4];
end
%------------- END OF CODE --------------
