%meProbe - Say whether the wall motion is there, and how far it can be amplified.
%
%   THE QUESTION THIS ANSWERS IS NOT "DOES THE PIPELINE RUN".  Motion magnification
%   always produces a moving picture; given noise and a large enough factor it
%   produces convincing wall motion from a recording containing none.  So before
%   any magnifier is written, this measures the thing that is supposed to be
%   amplified - the displacement of a vessel wall, in the raw recording, in pixels
%   and micrometres - and the ceiling that displacement puts on the amplification.
%   If the pulsation is already visible in a kymograph, the method will work and
%   the target is known.  If it is not, the number the artifact floor has to beat
%   is known instead.  Either answer is a result.
%
%   IT READS ONLY WHAT mePassA WROTE.  The recording is never opened again.
%
%   THE VESSEL THAT MOVES VISIBLY IS THE CONTROL FOR THE WHOLE PROJECT, so it is
%   measured first and separately.  Its displacement can be read off the raw data
%   without any magnification at all, which is the only way the later validation
%   has anything real to check a magnified movie against.  It is also the tightest
%   constraint: the bound is on the amplified displacement, so the vessel that
%   moves most sets the factor for a uniformly amplified field, and if it and the
%   small vessels differ by an order of magnitude that is the finding, not a
%   nuisance.
%
%   TWO PAGES, NOT ONE.  The library's page is 1400 by 700 and ten panels do not
%   fit on it legibly.  The first page carries the measurement - the field, the
%   spectrum, the kymograph, the wall traces, the amplification ceiling.  The
%   second carries the checks that decide whether the measurement can be believed
%   - the sample-size gate, the cardiac phase filling, the bleaching, the awake
%   animal's bouts.  Both are the canonical size, so they stack.
%
% Syntax:
%    probe = meProbe(s, fName)
%
% Inputs:
%    s     - settings from meSettings.  Fields read: outFolder, pixelSize, fftN,
%            minFrqIni, maxFrqIni, rangeFrq, vFR, nCuts, cutSep, kymoSeconds,
%            traceSeconds, motionWin, quietCoef, levels, reportLevel.
%    fName - path to the recording mePassA was run on.
%
% Outputs:
%    probe - the measurements, also saved as <recording>_ME_probe.mat:
%            .cardiac    heart rate, its harmonics, the beat-to-beat spread
%            .phaseFill  the cardiac phase histogram, and whether it fills
%            .eightBit   the sample-size gate as mePassA measured it
%            .bleach     percent of the starting brightness lost over the record
%            .snr        brightness over noise, per pixel and after binning
%            .cuts       one entry per vessel: calibre, wall edge width, the
%                        band-passed wall displacement in pixels and micrometres,
%                        and the amplification ceiling per pyramid level
%            .motion     frame-to-frame displacement, the bouts, the longest
%                        quiet stretch
%            .shutter    what the acquisition metadata says about readout
%            .verdict    the headline numbers as lines of text
%
% Example:
%    s     = meSettings;
%    probe = meProbe(s,'C:\Data\CXD\EPFL_20241114_2ADTF13BP.cxd');
%
% Dependencies: getFFT; meKymograph; meEpochs; meGate; the reporting module (reportOpen,
%               reportFigure, reportSave, reportClose); Signal Processing Toolbox
%               (findpeaks, butter, filtfilt); Image Processing Toolbox (imtophat,
%               bwskel).
% See also: mePassA, meReadRaw, meSettings, meKymograph, meEpochs, meGate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function probe = meProbe(s, fName)

fName = char(fName);
[recFolder,stem] = fileparts(fName);
outFolder = s.outFolder;
if isempty(outFolder), outFolder = fullfile(recFolder,'motionEnhancement'); end

loaded = load(fullfile(outFolder,[stem '_ME_passA.mat']),'pass');
pass   = loaded.pass;
tr     = load(pass.tracePath);
mn     = load(pass.meanPath);
fs     = pass.fs;
px     = pass.pixelSize;

s.rootFolder    = recFolder;
s.resultsFolder = outFolder;
rep = reportOpen(s,'Motion enhancement probe',{fName});
reportFile(rep,1,fName);

probe          = struct();
probe.stem     = stem;
probe.pass     = pass;
probe.verdict  = cell(0,1);
say = @(varargin) sprintf(varargin{:});

% ===== 1 and 2. the spectrum, and the heart rate =============================
% The trace that follows the bright pixels is the one the pulse lives in: the dye
% is in the plasma, so a wider vessel holds more of it.  Bleaching is removed
% first, because a monotone decay puts power straight into the vasomotion band.
x  = tr.traceVes(:);
xd = x - polyBaseline(x,3);
nFFT = min(s.fftN, 2^floor(log2(numel(xd))));
[amp,~,f] = getFFT(xd, fs, nFFT, 'cpu', 'subtractMean');
probe.spectrum = struct('f',f(:),'amp',amp(:),'fftN',nFFT,'df',f(2)-f(1));

% THE BEATS COME FROM meGate, WHICH IS THE ONE GATE THIS BLOCK SHIPS.  A second
% peak-finder here would be a second answer to the question the whole package turns
% on, and the two would disagree the first time either was touched.  What is wanted
% here is the FIDUCIALS - meGate's cycle rejection is not applied, because the fold
% below is over every quiet beat and the bout flag it uses for that is measured
% further down.  meGate's own acceptance is reported where it is used, in meEnhance,
% so there is one acceptance rate in the package and not two.
gate   = meGate(s, x, fs);
f0     = gate.f0;
minFrq = gate.band(1);
maxFrq = gate.band(2);
beatT  = gate.beatTime;
beatRate = 1./diff(beatT);

probe.cardiac = struct('f0',f0,'band',[minFrq maxFrq], ...
    'nBeats',numel(beatT)-1,'rateMedian',median(beatRate), ...
    'rateIQR',prctile(beatRate,[25 75]),'rateStd',std(beatRate), ...
    'beatTime',beatT,'beatRate',beatRate);

% respiration, and the vasomotion band, taken off the same spectrum
probe.respiration = bandPeak(f,amp,[0.8 4]);
probe.vasomotion  = bandPeak(f,amp,s.vFR);
probe.harmonics   = f0.*(1:4);

probe.verdict{end+1,1} = say('Heart rate %.2f Hz, beat to beat %.2f to %.2f Hz over %d beats', ...
    f0, probe.cardiac.rateIQR(1), probe.cardiac.rateIQR(2), probe.cardiac.nBeats);

% ===== gate 2. does the cardiac phase fill in ================================
% Ten frames per beat only makes a dense average cycle because the beats are not
% locked to the frame clock.  If they were, the sampled phases would pile up in a
% few places and the averaged cycle would come out banded.
nBin = 50;
ph   = nan(numel(tr.time),1);
for k = 1:numel(beatT)-1
    inBeat     = tr.time>=beatT(k) & tr.time<beatT(k+1);
    ph(inBeat) = (tr.time(inBeat)-beatT(k))./(beatT(k+1)-beatT(k));
end
counts = histcounts(ph(~isnan(ph)), linspace(0,1,nBin+1));
expect = mean(counts);
probe.phaseFill = struct('counts',counts(:),'nBin',nBin,'expected',expect, ...
    'emptyBins',sum(counts==0),'minOverMean',min(counts)/expect, ...
    'maxOverMean',max(counts)/expect,'cv',std(counts)/expect);
probe.verdict{end+1,1} = say(['Cardiac phase histogram: %d empty bins of %d, ' ...
    'thinnest bin %.2f and fullest %.2f of the average'], ...
    probe.phaseFill.emptyBins, nBin, probe.phaseFill.minOverMean, probe.phaseFill.maxOverMean);

% ===== gate 1. was eight bits enough ========================================
probe.eightBit = struct('steps',pass.eightBitSteps,'class',pass.rawClass, ...
    'sigmaBack',pass.sigmaBack,'sigmaBright',pass.sigmaBright, ...
    'meanBack',pass.meanBack,'meanBright',pass.meanBright, ...
    'lo',pass.lo,'hi',pass.hi, ...
    'growth',100*(sqrt(1+(0.289/pass.eightBitSteps)^2)-1));
probe.verdict{end+1,1} = say(['Sample size: a still patch carries %.2f steps of noise, ' ...
    'so rounding to eight bits adds %.1f%% and the copies were written as %s'], ...
    pass.eightBitSteps, probe.eightBit.growth, pass.rawClass);

% ===== 3. bleaching, and 9's edge-width drift comes later ====================
bleachPct = 100.*mn.bleach./mn.bleach(1);
probe.bleach = struct('time',mn.bleachTime,'percent',bleachPct, ...
    'lost',100-bleachPct(end),'minutes',tr.time(end)/60);
probe.verdict{end+1,1} = say('Brightness falls to %.1f%% of its start over %.2f minutes', ...
    bleachPct(end), probe.bleach.minutes);

% ===== 4. brightness over noise, per pixel and binned ========================
% THE STILL PATCH IS THE ONLY HONEST PLACE TO MEASURE THIS.  Inside a vessel the
% frame-to-frame change is mostly red cells crossing the lumen, not the camera, so
% a noise figure taken there would be several times too pessimistic.  The pre-pass
% measured the still patch on the recording itself, before any rounding; the
% binned figure is the same patch measured again on the copy every later step
% reads, so the gain from binning is a measurement and not a prediction of two.
[mB,gB] = meReadRaw(pass,'bin');
[mC,gC] = meReadRaw(pass,'crop');

binned = copyNoise(mB, pass, pass.backBox);
probe.snr = struct( ...
    'perPixelBack',   pass.meanBack  /pass.sigmaBack, ...
    'perPixelBright', pass.meanBright/pass.sigmaBright, ...
    'sigmaBackCounts',pass.sigmaBack, ...
    'binnedBack',     binned.mean/binned.sigma, ...
    'binnedSigma',    binned.sigma, ...
    'binningGain',    (binned.mean/binned.sigma)/(pass.meanBack/pass.sigmaBack));
probe.verdict{end+1,1} = say(['Brightness over noise on a still patch is %.0f per pixel ' ...
    'and %.0f after two-by-two binning, a gain of %.2f'], ...
    probe.snr.perPixelBack, probe.snr.binnedBack, probe.snr.binningGain);

% ===== 5. what the awake animal did, and WHEN it was still ===================
% This comes before the vessels on purpose.  The wall displacement measured below
% is the number Session 4 checks a magnified movie against, and measuring it
% across a running bout would report the animal's head instead of its artery.  So
% the quiet stretch is found first and the walls are measured inside it.
fprintf('  measuring frame-to-frame displacement over the whole recording\n');
probe.motion = motionQuality(mB, gB, s, pass, fs);
probe.verdict{end+1,1} = say(['Quiet for %.1f%% of the recording; the longest quiet stretch is ' ...
    '%.1f s starting at %.1f s'], 100*probe.motion.quietFraction, ...
    probe.motion.longestSeconds, probe.motion.longestStart);

% ===== 6, 7. the vessels =====================================================
fprintf('  finding vessel cuts\n');
cuts = findCuts(mn.meanFrame, s, pass.roi);
if isempty(cuts)
    error('meProbe:noVessels','No vessel large enough to cut across was found.');
end

nWant  = min(round(s.traceSeconds*fs), pass.sizeT);
nTrace = min(nWant, max(probe.motion.longestFrames,1));
f0Idx  = round(probe.motion.longestStart*fs) + 1 + floor((probe.motion.longestFrames-nTrace)/2);
f0Idx  = min(max(f0Idx,1), pass.sizeT-nTrace+1);
frames = f0Idx:f0Idx+nTrace-1;
% How many beats the amplitude is averaged over is the precision of every wall
% number below, so it is recorded rather than left to be inferred from the window.
nBeatsUsed = sum(beatT>=(frames(1)-1)/fs & beatT<=(frames(end)-1)/fs);
probe.traceWindow = struct('first',frames(1),'last',frames(end), ...
    'seconds',nTrace/fs,'startSeconds',(frames(1)-1)/fs,'wanted',nWant/fs, ...
    'beats',nBeatsUsed);
fprintf('  measuring the walls over %.1f s of the quietest stretch, from %.1f s, %d beats\n', ...
    nTrace/fs, (frames(1)-1)/fs, nBeatsUsed);

% EVERY POSITION BELOW IS IN PIXELS OF THE FULL-RESOLUTION FRAME, whichever copy
% the profile was sampled from.  The cut is laid out in the original frame's
% coordinates and only the SAMPLING is done in a copy, so a displacement never
% needs converting afterwards.  What the copy does change is how finely the wall
% is resolved, and copy.scale records that.
% THE WALLS ARE FITTED ON EVERY FRAME, not only inside the quiet window, because
% the cycle average below wants every quiet beat in the recording rather than the
% hundred and thirty that happen to sit in the longest unbroken stretch.  That is
% the expensive part of this function and it is where its couple of minutes go.
tAll = (0:pass.sizeT-1)'./fs;
for i = 1:numel(cuts)
    fprintf('  cut %d of %d, %.0f um across\n', i, numel(cuts), 2*cuts(i).radius*px);
    [K, cuts(i).pos, cuts(i).copy] = sampleCut(cuts(i), mB, gB, mC, gC, 1:pass.sizeT);
    w = meKymograph(K, cuts(i).pos);

    cuts(i).K           = K(:,frames);          % only the window is kept, for the page
    cuts(i).wallLeft    = w.left(frames);
    cuts(i).wallRight   = w.right(frames);
    cuts(i).diameter    = cuts(i).wallRight - cuts(i).wallLeft;
    cuts(i).edgeWidthPx = median([w.widthL;w.widthR],'omitnan');
    cuts(i).cycleTrace  = w.halfWidth;                % over the whole record
    cuts(i).cycleTime   = tAll;

    % Every displacement below is the in-band one: band-passed at the heart rate,
    % zero-phase, and reported as the amplitude of the equivalent sinusoid.
    cuts(i).ampLeftPx  = inBandAmplitude(cuts(i).wallLeft , [minFrq maxFrq], fs);
    cuts(i).ampRightPx = inBandAmplitude(cuts(i).wallRight, [minFrq maxFrq], fs);
    cuts(i).ampWallPx  = mean([cuts(i).ampLeftPx cuts(i).ampRightPx]);

    % The difference of the two walls cannot hold a rigid translation of the whole
    % field, so half of it is the wall motion the differential mode will see; the
    % per-wall figure above still holds whatever the brain was doing underneath.
    cuts(i).ampDiamPx = inBandAmplitude(cuts(i).diameter, [minFrq maxFrq], fs);
    cuts(i).ampDiffPx = cuts(i).ampDiamPx/2;

    % AND THE SAME MEASUREMENT WHERE THERE IS NOTHING TO FIND.  Band-passing a
    % noisy trace returns an amplitude whether or not anything is oscillating: a
    % wall position fitted to a fraction of a pixel per frame carries a per-frame
    % scatter several times larger than the motion being looked for, and a slice of
    % that scatter lands inside any band one cares to take.  So the identical
    % measurement is repeated in a band of the SAME WIDTH sitting between the heart
    % rate and its first harmonic, where no physiology puts anything.  What that
    % returns is this measurement's own noise floor, in the same pixels, and the
    % ratio of the two is the only thing that says whether the wall was seen.
    off = offBand([minFrq maxFrq], f0, fs);
    cuts(i).offBand    = off;
    cuts(i).nullWallPx = mean([inBandAmplitude(cuts(i).wallLeft ,off,fs) ...
                               inBandAmplitude(cuts(i).wallRight,off,fs)]);
    cuts(i).nullDiffPx = inBandAmplitude(cuts(i).diameter, off, fs)/2;
    cuts(i).wallOverNull = cuts(i).ampWallPx/max(cuts(i).nullWallPx,eps);
    cuts(i).diffOverNull = cuts(i).ampDiffPx/max(cuts(i).nullDiffPx,eps);

    % ---- and the estimator this whole project is built on -------------------
    % A band-pass does not get better with a longer recording: the noise inside a
    % band is a fixed fraction of the noise everywhere, so a thousand beats leave
    % the floor exactly where a hundred did.  Averaging the beats ON TOP OF EACH
    % OTHER does get better - the wall's motion adds up and the noise adds up as
    % its square root - and that is the whole reason the amplified-MRI work
    % magnifies an averaged cardiac cycle rather than a raw series.  So the same
    % trace is folded onto cardiac phase over every quiet beat in the recording,
    % and folded again on a period that no physiology occupies.  The second fold
    % is what the first one has to beat.
    fold = foldOnPhase(cuts(i).cycleTrace, cuts(i).cycleTime, beatT, s.nPhaseBins, probe.motion.bad, fs);
    null = foldOnPhase(cuts(i).cycleTrace, cuts(i).cycleTime, ...
                       beatT(1):1/(1.37*f0):beatT(end), s.nPhaseBins, probe.motion.bad, fs);
    cuts(i).cycle       = fold.cycle;
    cuts(i).cycleSem    = fold.sem;
    cuts(i).cycleBeats  = fold.n;
    cuts(i).foldPx      = fold.amplitude;
    cuts(i).foldNullPx  = null.amplitude;
    cuts(i).foldOverNull= fold.amplitude/max(null.amplitude,eps);
    cuts(i).foldUm      = fold.amplitude*px;

    % ---- and the same vessel again, through the other copy -------------------
    % A cut that fits inside the full-resolution window is measured a second time
    % from the binned one.  The two answers differ or they do not, and that decides
    % something real: whether the binned copy - which is the only one covering the
    % whole field - can carry a quantitative wall measurement at all, or whether
    % everything numerical has to stay inside the window.  Without this the two
    % groups of cuts are confounded, because Pass A placed that window on the most
    % pulsatile part of the field in the first place.
    cuts(i).foldAltPx = NaN; cuts(i).altCopy = '';
    if strcmp(cuts(i).copy.name,'crop')
        [Kalt,posAlt] = sampleCut(cuts(i), mB, gB, mC, gC, 1:pass.sizeT, 'bin');
        wAlt    = meKymograph(Kalt, posAlt);
        altFold = foldOnPhase(wAlt.halfWidth, tAll, beatT, s.nPhaseBins, probe.motion.bad, fs);
        cuts(i).foldAltPx = altFold.amplitude;
        cuts(i).altCopy   = 'bin';
    end

    % HOW FAR IT CAN BE AMPLIFIED IS SET BY THE DISPLACEMENT THAT IS ACTUALLY
    % THERE, which is the averaged one.  Taking the ceiling from the band-passed
    % figure would set it from this measurement's noise instead, and the noise is
    % larger, so the ceiling would come out too low - conservative by accident
    % rather than by measurement.
    cuts(i).alphaFold = wavelength(s.levels)./(4*max(fold.amplitude,eps));

    cuts(i).diameterPx  = median(cuts(i).diameter,'omitnan');
    cuts(i).diameterUm  = cuts(i).diameterPx *px;
    cuts(i).edgeWidthUm = cuts(i).edgeWidthPx*px;
    cuts(i).ampWallUm   = cuts(i).ampWallPx  *px;
    cuts(i).ampDiffUm   = cuts(i).ampDiffPx  *px;
    cuts(i).pulsatility = 100*2*cuts(i).ampDiffPx/cuts(i).diameterPx;

    % kept for the record: the ceiling the band-passed figure would have implied
    cuts(i).alphaBand = wavelength(s.levels)./(4*max(cuts(i).ampDiffPx,eps));

    % the lumen decorrelates frame to frame if the cells really are moving many
    % pixels between frames; the wall does not
    cuts(i).lumenLag1 = lagOneCorrelation(cuts(i).K);
end

% first and last minute, on the same cuts, to see the dye leaving the vessel
for i = 1:numel(cuts)
    cuts(i).edgeFirst = staticEdge(cuts(i), mn.meanFirst);
    cuts(i).edgeLast  = staticEdge(cuts(i), mn.meanLast );
end
probe.cuts = cuts;

[~,iBig] = max([cuts.diameterPx]);
[~,iSml] = min([cuts.diameterPx]);
probe.iLargest  = iBig;
probe.iSmallest = iSml;
probe.levels    = s.levels;
probe.lambda    = wavelength(s.levels);

probe.verdict{end+1,1} = say(['Largest vessel %.0f um across, wall moves %.3f px (%.2f um) in band, ' ...
    'differential %.3f px (%.2f um), pulsatility %.1f%%'], ...
    cuts(iBig).diameterUm, cuts(iBig).ampWallPx, cuts(iBig).ampWallUm, ...
    cuts(iBig).ampDiffPx, cuts(iBig).ampDiffUm, cuts(iBig).pulsatility);
probe.verdict{end+1,1} = say(['Smallest vessel %.0f um across, wall moves %.3f px (%.2f um) in band, ' ...
    'differential %.3f px (%.2f um), pulsatility %.1f%%'], ...
    cuts(iSml).diameterUm, cuts(iSml).ampWallPx, cuts(iSml).ampWallUm, ...
    cuts(iSml).ampDiffPx, cuts(iSml).ampDiffUm, cuts(iSml).pulsatility);

% THE LINES THAT DECIDE WHETHER ANY OF THE ABOVE IS A MEASUREMENT.
probe.detection = struct('offBand',cuts(1).offBand, ...
    'wallOverNull',[cuts.wallOverNull],'diffOverNull',[cuts.diffOverNull], ...
    'nullWallPx',[cuts.nullWallPx],'nullDiffPx',[cuts.nullDiffPx], ...
    'seenByBand',[cuts.diffOverNull]>=2, ...
    'foldPx',[cuts.foldPx],'foldNullPx',[cuts.foldNullPx], ...
    'foldOverNull',[cuts.foldOverNull],'foldBeats',[cuts.cycleBeats], ...
    'seenByFold',[cuts.foldOverNull]>=2);
probe.verdict{end+1,1} = say(['Band-passed against an empty band of the same width at %.1f to %.1f Hz, ' ...
    'the wall stands %.1f times its own noise on the largest vessel and %.1f on the smallest; ' ...
    '%d of %d vessels clear twice that floor'], ...
    cuts(1).offBand(1), cuts(1).offBand(2), cuts(iBig).diffOverNull, cuts(iSml).diffOverNull, ...
    sum(probe.detection.seenByBand), numel(cuts));
probe.verdict{end+1,1} = say(['Averaged over %d quiet beats instead, the wall moves %.3f px (%.2f um) ' ...
    'on the largest vessel against a floor of %.3f px, which is %.1f times it; ' ...
    '%d of %d vessels clear twice their floor'], ...
    cuts(iBig).cycleBeats, cuts(iBig).foldPx, cuts(iBig).foldUm, cuts(iBig).foldNullPx, ...
    cuts(iBig).foldOverNull, sum(probe.detection.seenByFold), numel(cuts));
probe.verdict{end+1,1} = say(['Amplification ceiling at the coarsest level used: %.0f on the largest ' ...
    'vessel, %.0f on the smallest, a factor of %.1f apart'], ...
    cuts(iBig).alphaFold(end), cuts(iSml).alphaFold(end), ...
    cuts(iSml).alphaFold(end)/max(cuts(iBig).alphaFold(end),eps));

alt = ~isnan([cuts.foldAltPx]);
if any(alt)
    r = [cuts(alt).foldPx]./max([cuts(alt).foldAltPx],eps);
    probe.binPenalty = struct('cuts',find(alt),'full',[cuts(alt).foldPx], ...
        'binned',[cuts(alt).foldAltPx],'ratio',r);
    probe.verdict{end+1,1} = say(['Measured through the binned copy instead, the same %d vessels ' ...
        'give %.2f to %.2f times the full-resolution answer'], nnz(alt), min(r), max(r));
end

% ===== 8. the dye =============================================================
ef = [cuts.edgeFirst]; el = [cuts.edgeLast];
probe.edgeDrift = struct('first',ef,'last',el,'percent',100*(el./ef-1));
probe.verdict{end+1,1} = say('Wall edge width changes by %+.1f%% on average between the first and the last minute', ...
    mean(probe.edgeDrift.percent,'omitnan'));

% The arithmetic behind "the lumen texture cannot be coherent motion": a red cell
% crossing this many pixels between frames is far outside the phase range of any
% pyramid level, so it lands as broadband decorrelation rather than as a
% direction.
probe.rbc = struct('speedMmS',[10 20], ...
    'pxPerFrame',[10 20].*1000./px./fs, ...
    'lumenLag1',[cuts.lumenLag1]);
probe.verdict{end+1,1} = say(['A red cell at 10 to 20 mm/s crosses %.0f to %.0f pixels per frame, ' ...
    'and the lumen holds a frame-to-frame correlation of %.2f'], ...
    probe.rbc.pxPerFrame(1), probe.rbc.pxPerFrame(2), median(probe.rbc.lumenLag1,'omitnan'));

% ===== 9. rolling shutter =====================================================
probe.shutter = shutterLines(pass.metadata);
if isempty(pass.metadata)
    probe.verdict{end+1,1} = ['This recording carries no acquisition metadata at all, so how the ' ...
        'sensor was read has to come from whoever recorded it'];
elseif isempty(probe.shutter)
    probe.verdict{end+1,1} = ['The acquisition metadata says nothing about readout time or ' ...
        'trigger mode, so how the sensor was read has to come from whoever recorded it'];
else
    probe.verdict{end+1,1} = say('The acquisition metadata carries %d readout or trigger entries', ...
        numel(probe.shutter));
end

% ===== the pages ==============================================================
drawMeasurement(reportFigure(rep,'probe'),       probe, mn, s, pass, rep);
drawChecks(     reportFigure(rep,'probeChecks'), probe,     pass, rep);
closed = reportClose(rep);
probe.images = closed.images;

save(fullfile(outFolder,[stem '_ME_probe.mat']),'probe','-v7.3','-nocompression');

fprintf('\n');
for i = 1:numel(probe.verdict), fprintf('  %s\n', probe.verdict{i}); end
end

% =====================================================================
function b = polyBaseline(x,order)
%polyBaseline  The slow decay of the dye, as a low-order polynomial in time.
n = numel(x);
t = ((1:n)'-1)./(n-1);
b = polyval(polyfit(t,x,order),t);
end

% =====================================================================
function y = bandpass0(x, band, fs)
%bandpass0  Zero-phase band-pass.  filtfilt, never filter: a causal lag that
%   varies across the field is the propagating wave this block must not invent.
w = band./(fs/2);
w = min(max(w,1e-6),0.999);
[bb,aa] = butter(2, w, 'bandpass');
y = filtfilt(bb, aa, double(x(:)));
end

% =====================================================================
function a = inBandAmplitude(x, band, fs)
%inBandAmplitude  Amplitude of the equivalent sinusoid inside the band.
%   Root mean square times root two, so a pure sinusoid of amplitude A reports A,
%   and gaps left by a frame whose walls could not be fitted are skipped.
x = x(:);
ok = isfinite(x);
if nnz(ok) < 0.5*numel(x), a = NaN; return; end
x(~ok) = interp1(find(ok), x(ok), find(~ok), 'linear', 'extrap');
a = sqrt(2)*std(bandpass0(x, band, fs));
end

% =====================================================================
function out = foldOnPhase(x, tx, bnd, nBin, bad, fs)
%foldOnPhase  Average one beat out of all of them.
%
%   Every interval between two consecutive boundaries is stretched onto the same
%   phase axis and the stretched copies are averaged.  Anything locked to the
%   beat survives at full size; anything not locked to it falls as the square root
%   of the number of beats.  That is where the sensitivity in this whole approach
%   comes from, and it is why amplified MRI amplifies a gated cine rather than a
%   raw series.
%
%   STRETCHED, NOT WINDOWED, because the heart rate wanders - here between 11.1 and
%   12.5 Hz - and a fixed-length window would smear the late part of the beat
%   across a phase bin or two, which flattens exactly the peak being measured.
%
%   EACH BEAT IS CENTRED ON ITS OWN MEAN, so a slow drift underneath contributes
%   nothing.  What comes out is the SHAPE of the beat.
%
%   A BEAT THE ANIMAL MOVED THROUGH IS DROPPED, not repaired.  Marked, and counted
%   in n, so a low count is visible rather than inferred.
x = x(:); tx = tx(:); bnd = bnd(:);
p   = ((0:nBin-1)'+0.5)./nBin;
acc = zeros(nBin,1); acc2 = zeros(nBin,1); n = 0;
for k = 1:numel(bnd)-1
    t0 = bnd(k); t1 = bnd(k+1);
    if ~isfinite(t0) || ~isfinite(t1) || t1<=t0, continue; end
    i0 = max(1,floor(t0*fs)+1);  i1 = min(numel(x),ceil(t1*fs)+1);
    if i1<=i0 || i1>numel(bad),  continue; end
    if any(bad(i0:i1)),          continue; end
    v = interp1(tx, x, t0 + p.*(t1-t0), 'linear', NaN);
    if any(~isfinite(v)),        continue; end
    v   = v - mean(v);
    acc = acc + v;  acc2 = acc2 + v.^2;  n = n+1;
end
if n < 10
    out = struct('cycle',nan(nBin,1),'sem',nan(nBin,1),'n',n,'amplitude',NaN);
    return
end
cyc = acc./n;
sd  = sqrt(max(acc2./n - cyc.^2, 0));
out = struct('cycle',cyc,'sem',sd./sqrt(n),'n',n,'amplitude',(max(cyc)-min(cyc))/2);
end

% =====================================================================
function b = offBand(band, f0, fs)
%offBand  A band of the same width as the heart-rate band, put where no physiology
%   is: centred half way between the fundamental and its first harmonic.  Same
%   width matters - a wider band would collect more noise and a narrower one less,
%   and the comparison would then be between two different measurements.
w = diff(band);
c = 1.5*f0;
b = [c-w/2, c+w/2];
b = [max(b(1), band(2)+0.1), min(b(2), 0.45*fs)];
end

% =====================================================================
function lam = wavelength(levels)
%wavelength  The wavelength a pyramid level is centred on, in pixels of the
%   full-resolution frame.  Level one holds the finest detail, and each level up
%   doubles.
lam = 2.^(levels(:)'+1);
end

% =====================================================================
function p = bandPeak(f,amp,band)
%bandPeak  The strongest line inside a band, and how far it stands above it.
idx = f(:)>=band(1) & f(:)<=band(2);
if ~any(idx), p = struct('f',NaN,'amp',NaN,'ratio',NaN); return; end
a = amp(:); a(~idx) = -Inf;
[v,i] = max(a);
p = struct('f',f(i),'amp',v,'ratio',v/median(amp(idx)));
end

% =====================================================================
function out = copyNoise(m, pass, boxFull)
%copyNoise  Brightness and noise of a still patch, measured on the working copy
%   and converted back to the recording's own counts.
%
%   The copy is rounded, so its samples carry a rounding noise of 0.289 of a step
%   on top of the recording's.  That is removed here rather than reported as if it
%   were the camera's: it is a known, independent term, so it comes out in
%   quadrature.  Everything is then multiplied back into counts, so the figure
%   stands beside the pre-pass's own measurement without either being rescaled.
b  = pass.binFactor;
r0 = max(1,floor((boxFull(1,1)-1)/b)+1);  r1 = min(pass.binSize(1), ceil(boxFull(1,2)/b));
c0 = max(1,floor((boxFull(2,1)-1)/b)+1);  c1 = min(pass.binSize(2), ceil(boxFull(2,2)/b));

n  = min(300, pass.sizeT);
t0 = max(1, round(pass.sizeT/2) - round(n/2));
X  = single(m.Data.d(r0:r1, c0:c1, t0:t0+n-1));

stepCounts = 1/pass.gain;              % counts per step, whatever the sample size is
sigmaStep  = median(std(diff(X,1,3),0,3)./sqrt(2), 'all');
sigmaStep  = sqrt(max(sigmaStep^2 - 0.289^2, 0));

out.sigma = sigmaStep*stepCounts;
out.mean  = pass.lo + median(X,'all')*stepCounts;
end

% =====================================================================
function cuts = findCuts(meanFrame, s, roi)
%findCuts  Where to cut across a vessel, and in which direction.
%   The dye fills the lumen, so a vessel is a bright ridge.  A top-hat removes the
%   uneven background it sits on, the skeleton of what is left runs along each
%   vessel, and the distance to the nearest edge at a skeleton point is the local
%   radius - which is what makes the calibre available before anything is fitted.
%
%   THE LARGEST VESSEL IS TAKEN FIRST AND ON PURPOSE.  It is the one whose motion
%   can be measured without any magnification at all, so it is the control the
%   whole project is checked against, and it is also the one that sets the
%   amplification ceiling.
%
%   AND NO TWO CUTS MEASURE THE SAME VESSEL.  Distance alone does not prevent it:
%   two points a hundred pixels apart along one artery are one measurement made
%   twice, and a probe that reports it as two has learned nothing about calibre.
%   So one cut is taken per connected piece of vasculature first, widest piece
%   first, and only when the pieces run out are further cuts placed inside a piece
%   - and then only where the vessel is a QUARTER narrower than anything already
%   cut, which is the only way a second cut on one tree says something new.
%   AND THE HOLES IN A VASCULAR NETWORK ARE NOT HOLES.  Filling them - the reflex
%   when a mask comes out speckled - welds the whole tree into one solid region,
%   and the distance to its edge is then a couple of hundred pixels in the middle
%   of the field rather than the radius of a vessel.  Measured here on the
%   reference recording: it produced a "vessel" 651 um across and a cut longer than
%   the frame.  A bright-lumen vessel is a ridge, not an annulus; there is nothing
%   to fill.
I  = double(meanFrame);
th = imtophat(imgaussfilt(I,1), strel('disk',40));
bw = imbinarize(mat2gray(th));
bw = bwareaopen(bw, 300);

D  = bwdist(~bw);
sk = bwskel(bw);
[ry,rx] = find(sk);
if isempty(ry), cuts = []; return; end
lin = sub2ind(size(D),ry,rx);
rad = D(lin);
lab = double(bwlabel(bw));
lab = lab(lin);

cuts = struct('center',{},'normal',{},'radius',{},'piece',{});

% ---- one cut per piece, widest piece first ---------------------------------
pieces = unique(lab);
best   = zeros(numel(pieces),1);
for k = 1:numel(pieces)
    best(k) = max(rad(lab==pieces(k)));
end
[~,order] = sort(best,'descend');
for k = order(:)'
    if numel(cuts)>=s.nCuts, break; end
    inPiece = find(lab==pieces(k));
    [~,iw]  = sort(rad(inPiece),'descend');
    cuts    = tryPoints(cuts, cropFirst(inPiece(iw),ry,rx,rad,roi), ry, rx, rad, pieces(k), s, false, I);
end

% ---- then inside the pieces already used, where the calibre is different ---
if numel(cuts) < s.nCuts
    [~,iw] = sort(rad,'descend');
    cuts   = tryPoints(cuts, cropFirst(iw,ry,rx,rad,roi), ry, rx, rad, [], s, true, I);
end

if isempty(cuts), return; end
[~,order] = sort([cuts.radius],'descend');
cuts = cuts(order);
end

% =====================================================================
function idx = cropFirst(idx, ry, rx, rad, roi)
%cropFirst  Same candidates, but the ones that fit inside the full-resolution copy
%   come first.  A wall measured at half resolution is measured through twice the
%   blur, and mePassA wrote a full-resolution window precisely so that some of this
%   session's numbers could avoid that.  Preference, not a requirement: a vessel
%   outside the window is still worth cutting, and it says which copy it used.
idx = idx(:);
L   = 3*rad(idx) + 6;
in  = ry(idx)-L >= roi(1,1) & ry(idx)+L <= roi(1,2) & ...
      rx(idx)-L >= roi(2,1) & rx(idx)+L <= roi(2,2);
idx = [idx(in); idx(~in)];
end

% =====================================================================
function cuts = tryPoints(cuts, candidates, ry, rx, rad, piece, s, needNewCalibre, meanImg)
%tryPoints  Walk a candidate list and take the first point that is far enough from
%   every cut already made, different enough in calibre when that is being asked
%   for, has enough skeleton around it to say which way the vessel runs, AND
%   actually looks like a vessel when it is cut across.
%
%   The last of those is the one that matters.  A mask, a skeleton and a distance
%   transform will happily nominate a point that is not on a vessel at all - a
%   junction where several vessels meet, a bright patch, a place where the mask
%   welded two things together - and every number measured there afterwards is
%   fiction rather than noise.  So the candidate is cut across on the time-mean
%   image before it is accepted, and it has to come back with two walls, a width
%   that agrees with what the distance transform promised, and edges that are
%   sharp compared with the vessel they bound.
for j = candidates(:)'
    if numel(cuts)>=s.nCuts, return; end
    p = [ry(j) rx(j)];
    if ~isempty(cuts)
        if min(vecnorm(reshape([cuts.center],2,[])-p(:),2,1)) < s.cutSep, continue; end
        if needNewCalibre && min(abs([cuts.radius]./rad(j) - 1)) < 0.25, continue; end
    end
    n = localNormal(ry,rx,p,rad(j));
    if isempty(n), continue; end
    if ~looksLikeAVessel(p, n, rad(j), meanImg), continue; end
    thisPiece = piece;
    if isempty(thisPiece), thisPiece = NaN; end
    cuts(end+1) = struct('center',p,'normal',n,'radius',rad(j),'piece',thisPiece); %#ok<AGROW>
    if ~needNewCalibre, return; end          % one per piece on the first sweep
end
end

% =====================================================================
function tf = looksLikeAVessel(p, n, radius, meanImg)
%looksLikeAVessel  Cut across the candidate on the time-mean image and see.
tf  = false;
pos = (-(3*radius+6):0.25:(3*radius+6))';
pt  = p + pos.*n;
if any(pt(:,1)<1 | pt(:,1)>size(meanImg,1) | pt(:,2)<1 | pt(:,2)>size(meanImg,2)), return; end
prof = interp2(meanImg, pt(:,2), pt(:,1), 'linear', NaN);
if any(isnan(prof)), return; end

w = meKymograph(prof, pos);
if ~isfinite(w.left) || ~isfinite(w.right), return; end
d = w.right - w.left;
if d <= 0, return; end
if d < 0.5*2*radius || d > 2*2*radius, return; end     % the width the mask promised
edge = mean([w.widthL w.widthR],'omitnan');
if ~isfinite(edge) || edge > 0.5*d, return; end        % a wall, not a slope
tf = true;
end

% =====================================================================
function n = localNormal(ry,rx,p,radius)
%localNormal  The direction across the vessel at p, from the skeleton around it.
%   The skeleton runs ALONG the vessel, so the direction its neighbours spread in
%   is the axis and the cut is perpendicular to it.  A point with no neighbours to
%   define an axis is refused rather than guessed at.
w    = max(12, 3*radius);
near = abs(ry-p(1))<=w & abs(rx-p(2))<=w;
if nnz(near) < 5, n = []; return; end
P = [ry(near) rx(near)];
P = P - mean(P,1);
[V,Dg] = eig(cov(P));
[~,iAxis] = max(diag(Dg));
axisDir = V(:,iAxis);
n = [-axisDir(2); axisDir(1)]';
n = n./norm(n);
end

% =====================================================================
function [K,pos,copy] = sampleCut(cut, mB, gB, mC, gC, frames, forceCopy)
%sampleCut  Intensity along the cut, for every frame: the kymograph.
%   The full-resolution copy is used whenever the cut fits inside it, because half
%   a pixel of resolution is half a pixel of the answer.  Otherwise the binned
%   copy of the whole field is used, and the caller is told which by copy.scale.
%
%   The frames are read once, as one block around the cut, and turned into a
%   kymograph by a single sparse multiply: the interpolation weights do not change
%   from frame to frame, so they are built once rather than thirty thousand times.
L   = 3*cut.radius + 6;
pos = (-L:0.25:L)';
ptFull = cut.center + pos.*cut.normal;

if nargin>6 && strcmp(forceCopy,'bin')
    copy = struct('name','bin','origin',gB.origin,'scale',gB.scale,'size',gB.size(1:2));
else
    copy = pickCopy(ptFull, gC, gB);
end
pt   = (ptFull - copy.origin)./copy.scale + 1 - (copy.scale-1)/(2*copy.scale);

r0 = max(1,floor(min(pt(:,1)))-1);  r1 = min(copy.size(1),ceil(max(pt(:,1)))+1);
c0 = max(1,floor(min(pt(:,2)))-1);  c1 = min(copy.size(2),ceil(max(pt(:,2)))+1);
if strcmp(copy.name,'crop'), X = mC.Data.d(r0:r1, c0:c1, frames);
else,                        X = mB.Data.d(r0:r1, c0:c1, frames);
end
nr = r1-r0+1; nc = c1-c0+1;
X  = reshape(X, nr*nc, numel(frames));

% The block stays at its own sample size and only the four rows each sample point
% needs are widened.  Turning the whole block into double first would be four
% hundred megabytes for a cut that needs eight.
[idx,wgt] = bilinearWeights(pt(:,1)-r0+1, pt(:,2)-c0+1, nr, nc);
K = wgt(:,1).*single(X(idx(:,1),:)) + wgt(:,2).*single(X(idx(:,2),:)) + ...
    wgt(:,3).*single(X(idx(:,3),:)) + wgt(:,4).*single(X(idx(:,4),:));
K = double(K);
end

% =====================================================================
function copy = pickCopy(ptFull, gC, gB)
%pickCopy  The full-resolution window if the cut is inside it, the binned field if
%   not.  Nothing is silently truncated: the cut either fits or it moves copy.
inCrop = all(ptFull(:,1)>=gC.origin(1) & ptFull(:,1)<=gC.origin(1)+gC.size(1)-1 & ...
             ptFull(:,2)>=gC.origin(2) & ptFull(:,2)<=gC.origin(2)+gC.size(2)-1);
if inCrop
    copy = struct('name','crop','origin',gC.origin,'scale',gC.scale,'size',gC.size(1:2));
else
    copy = struct('name','bin' ,'origin',gB.origin,'scale',gB.scale,'size',gB.size(1:2));
end
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
function e = staticEdge(cut, img)
%staticEdge  The wall edge width on a single averaged image, for the same cut.
L   = 3*cut.radius + 6;
pos = (-L:0.25:L)';
pt  = cut.center + pos.*cut.normal;
p   = interp2(double(img), pt(:,2), pt(:,1), 'linear', NaN);
if any(isnan(p)), e = NaN; return; end
w = meKymograph(p, pos);
e = mean([w.widthL w.widthR],'omitnan');
end

% =====================================================================
function r = lagOneCorrelation(K)
%lagOneCorrelation  How much a row of the kymograph inside the lumen resembles
%   itself one frame later.  Red cells crossing tens of pixels between frames
%   cannot leave a coherent direction behind, so the lumen should decorrelate
%   almost completely while the wall does not.  The middle of the cut IS the
%   lumen: the cut is laid out from the vessel's own centre line.
mid = round(size(K,1)/2);
band = max(1,mid-3):min(size(K,1),mid+3);
x = mean(K(band,:),1)';
x = x - mean(x);
if numel(x)<3 || std(x)==0, r = NaN; return; end
r = corr(x(1:end-1), x(2:end));
end

% =====================================================================
function mq = motionQuality(mB, gB, s, pass, fs)
%motionQuality  How much the whole field moves between one frame and the next.
%   An awake animal runs, whisks and grooms, and a bout displaces the brain by
%   many pixels.  Amplified, that does not merely look bad: it leaves the range
%   the magnifier is valid over and tears the picture up.  So the bouts are found
%   here and the continuous mode is given the longest stretch between them.
%
%   Phase correlation, not registration: this only has to SAY how far the field
%   moved.  Nothing is resampled, because resampling is what would destroy the
%   fraction of a pixel the rest of the session is measuring.
w  = min(s.motionWin, min(gB.size(1:2)));
r0 = max(1,round((gB.size(1)-w)/2)); c0 = max(1,round((gB.size(2)-w)/2));
win = hann(w)*hann(w)';

nT = pass.sizeT;
d  = zeros(nT,1);
prevF = [];
step  = max(1,round(nT/10));
for t = 1:nT
    A = single(mB.Data.d(r0:r0+w-1, c0:c0+w-1, t)).*win;
    F = fft2(A);
    if ~isempty(prevF)
        R = F.*conj(prevF);
        R = R./max(abs(R),eps);
        c = real(ifft2(R));
        [dy,dx] = peakShift(c,w);
        d(t) = hypot(dy,dx);
    end
    prevF = F;
    if mod(t,step)==0, fprintf('    %3.0f%%\n', 100*t/nT); end
end
d = d.*gB.scale;                       % binned pixels back to frame pixels

% THE CLASSIFICATION IS meEpochs', NOT A SECOND COPY OF IT.  This function's job is
% the expensive one - a phase correlation over every frame of the recording - and
% turning that series into a bout flag is a rule that every later step has to agree
% with.  Two implementations of it would disagree the first time either was touched.
mq = meEpochs(s, d, fs);
mq.displacement = d;
mq.time         = (0:nT-1)'./fs;
mq.corrWindow   = w;
mq.pixelSize    = pass.pixelSize;
end

% =====================================================================
function [dy,dx] = peakShift(c,w)
%peakShift  The correlation peak, wrapped to a signed shift and refined between
%   samples by a parabola through its two neighbours.
[~,i] = max(c(:));
[iy,ix] = ind2sub([w w],i);
dy = subPixel(c, iy, ix, 1, w);
dx = subPixel(c, iy, ix, 2, w);
end

% =====================================================================
function v = subPixel(c, iy, ix, dim, w)
if dim==1
    a = c(mod(iy-2,w)+1, ix); b = c(iy,ix); cc = c(mod(iy,w)+1, ix);
    base = iy-1;
else
    a = c(iy, mod(ix-2,w)+1); b = c(iy,ix); cc = c(iy, mod(ix,w)+1);
    base = ix-1;
end
den = a - 2*b + cc;
delta = 0;
if den ~= 0, delta = 0.5*(a-cc)/den; end
v = base + delta;
if v > w/2, v = v - w; end
end


% =====================================================================
function lines = shutterLines(metadata)
%shutterLines  What the acquisition metadata says about how the sensor is read.
%   A sensor read out row by row samples the top and the bottom of a frame at
%   different times, which is a systematic top-to-bottom phase ramp - and a phase
%   ramp across the field is indistinguishable from a wave travelling across it.
%   That matters here more than it would elsewhere, because this library measures
%   pulse-wave propagation.
%
%   A GLOBAL SHUTTER SETTLES IT AND SOME RECORDINGS WILL NOT SAY.  The reference
%   .cxd carries no metadata table at all, and its camera was confirmed global
%   shutter by the person who recorded it.  So this returns what the file knows
%   and nothing more: an empty answer means ASK, not that there is a problem.
lines = {};
if isempty(metadata), return; end
key = {'readout','shutter','trigger','scan','exposure','line time','sensor','binning','frame'};
for i = 1:numel(metadata)
    l = lower(metadata{i});
    if any(cellfun(@(k) contains(l,k), key))
        lines{end+1,1} = metadata{i}; %#ok<AGROW>
    end
end
end

% =====================================================================
function drawMeasurement(fh, probe, mn, s, pass, rep)
%drawMeasurement  Page one: what was measured.
t = tiledlayout(fh,2,3,'Padding','compact','TileSpacing','compact');

nexttile(t);
imagesc(mn.meanFrame); axis image off; colormap(gca,gray); hold on
rectangle('Position',[pass.roi(2,1) pass.roi(1,1) diff(pass.roi(2,:)) diff(pass.roi(1,:))], ...
    'EdgeColor',[0.2 0.7 1],'LineWidth',1);
for i = 1:numel(probe.cuts)
    cu = probe.cuts(i);
    L  = 3*cu.radius+6;
    p1 = cu.center - L*cu.normal;  p2 = cu.center + L*cu.normal;
    plot([p1(2) p2(2)],[p1(1) p2(1)],'-','Color',[1 0.4 0.1],'LineWidth',1.2);
    text(cu.center(2)+8, cu.center(1), sprintf('%d',i),'Color',[1 0.8 0.2]);
end
title(sprintf('Field, %.0f um per 100 px, cuts 1 to %d', 100*pass.pixelSize, numel(probe.cuts)));

nexttile(t);
loglog(probe.spectrum.f, probe.spectrum.amp,'k'); hold on; grid on
yl = ylim;
patch([s.vFR(1) s.vFR(2) s.vFR(2) s.vFR(1)],[yl(1) yl(1) yl(2) yl(2)], ...
    [0.3 0.6 1],'FaceAlpha',0.15,'EdgeColor','none');
for h = probe.harmonics
    if h < probe.spectrum.f(end), xline(h,'--','Color',[0.9 0.3 0.1]); end
end
if isfinite(probe.respiration.f), xline(probe.respiration.f,':','Color',[0.1 0.6 0.2]); end
xlim([0.02 50]); ylim(yl);
xlabel('Hz'); ylabel('amplitude');
title(sprintf('Spectrum, heart rate %.2f Hz', probe.cardiac.f0));

cu = probe.cuts(probe.iLargest);
nexttile(t);
nk = min(round(s.kymoSeconds*pass.fs), size(cu.K,2));
imagesc((0:nk-1)./pass.fs, cu.pos.*pass.pixelSize, cu.K(:,1:nk));
colormap(gca,gray); xlabel('s'); ylabel('um across the vessel');
title(sprintf('Kymograph, largest vessel, %.0f um across', cu.diameterUm));

% The averaged beat, beside the same average taken on a period nothing occupies.
% One panel, and it is the one that answers the session's question.
nexttile(t);
ph = ((0:numel(cu.cycle)-1)+0.5)./numel(cu.cycle);
errorbar(ph, cu.cycle.*pass.pixelSize, cu.cycleSem.*pass.pixelSize, ...
    '-o','LineWidth',1.4,'Color',[0.85 0.2 0.1],'MarkerSize',3); hold on
yline(0,':k'); grid on
xlabel('phase within the beat'); ylabel('um');
title(sprintf('Averaged beat, %d beats, %.2f um, %.1f times its floor', ...
    cu.cycleBeats, cu.foldUm, cu.foldOverNull));

nexttile(t);
dUm = [probe.cuts.diameterUm];
[dUm,ord] = sort(dUm);
semilogy(dUm, [probe.cuts(ord).ampDiffPx],'s-','LineWidth',1.4,'Color',[0.1 0.4 0.9]); hold on
semilogy(dUm, [probe.cuts(ord).nullDiffPx],'s:','LineWidth',1.2,'Color',[0.5 0.6 0.8]);
semilogy(dUm, [probe.cuts(ord).foldPx],'o-','LineWidth',1.6,'Color',[0.85 0.2 0.1]);
semilogy(dUm, [probe.cuts(ord).foldNullPx],'o:','LineWidth',1.2,'Color',[0.5 0.5 0.5]);
grid on; xlabel('vessel diameter, um'); ylabel('pixels');
legend({'band-passed','its empty band','averaged beat','its empty period'}, ...
    'Location','best','Box','off');
title(sprintf('Wall displacement: %d of %d seen by averaging', ...
    sum(probe.detection.seenByFold), numel(probe.cuts)));

nexttile(t);
for i = 1:numel(probe.cuts)
    semilogy(probe.levels, probe.cuts(i).alphaFold,'o-','LineWidth',1.2); hold on
end
grid on; xlabel('pyramid level'); ylabel('largest usable factor');
legend(arrayfun(@(c) sprintf('%.0f um',c.diameterUm), probe.cuts,'UniformOutput',false), ...
    'Location','best','Box','off');
title('Amplification ceiling, from the averaged beat');

reportSave(rep, fh, 'probe');
end

% =====================================================================
function drawChecks(fh, probe, pass, rep)
%drawChecks  Page two: whether page one can be believed.
t = tiledlayout(fh,2,3,'Padding','compact','TileSpacing','compact');

nexttile(t);
bar(0.5/probe.phaseFill.nBin + (0:probe.phaseFill.nBin-1)./probe.phaseFill.nBin, ...
    probe.phaseFill.counts,1,'FaceColor',[0.3 0.5 0.8],'EdgeColor','none'); hold on
yline(probe.phaseFill.expected,'--k');
xlabel('phase within the beat'); ylabel('frames'); grid on
title(sprintf('Phase filling, %d empty of %d', probe.phaseFill.emptyBins, probe.phaseFill.nBin));

nexttile(t);
bar([probe.eightBit.sigmaBack probe.eightBit.sigmaBright].*255./(probe.eightBit.hi-probe.eightBit.lo), ...
    'FaceColor',[0.3 0.6 0.4],'EdgeColor','none'); hold on
yline(1,'--r','LineWidth',1.2);
set(gca,'XTickLabel',{'dim patch','bright patch'});
ylabel('steps of an eight-bit sample'); grid on
title(sprintf('Sample size, written as %s', probe.eightBit.class));

nexttile(t);
plot(probe.bleach.time./60, probe.bleach.percent,'k','LineWidth',1.2); grid on
xlabel('minutes'); ylabel('% of the start');
title(sprintf('Brightness, %.1f%% lost', probe.bleach.lost));

nexttile(t);
plot(probe.motion.time, probe.motion.displacement,'-','Color',[0.4 0.4 0.4]); hold on
yline(probe.motion.threshold,'--r');
xr = probe.motion.longestStart + [0 probe.motion.longestSeconds];
yl = ylim;
patch([xr(1) xr(2) xr(2) xr(1)],[yl(1) yl(1) yl(2) yl(2)],[0.2 0.8 0.3], ...
    'FaceAlpha',0.2,'EdgeColor','none');
ylim(yl); grid on; xlabel('s'); ylabel('pixels between frames');
title(sprintf('Quiet %.0f%%, longest stretch %.1f s', 100*probe.motion.quietFraction, ...
    probe.motion.longestSeconds));

nexttile(t);
histogram(probe.cardiac.beatRate, 40,'FaceColor',[0.6 0.3 0.5],'EdgeColor','none'); grid on
xlabel('Hz'); ylabel('beats');
title(sprintf('Beat to beat, median %.2f Hz', probe.cardiac.rateMedian));

nexttile(t);
plot([probe.cuts.diameterUm], probe.edgeDrift.first.*pass.pixelSize,'o-'); hold on
plot([probe.cuts.diameterUm], probe.edgeDrift.last .*pass.pixelSize,'s--');
grid on; xlabel('vessel diameter, um'); ylabel('wall edge width, um');
legend({'first minute','last minute'},'Location','best','Box','off');
title(sprintf('Edge width, %+.1f%% over the record', mean(probe.edgeDrift.percent,'omitnan')));

reportSave(rep, fh, 'probeChecks');
end
%------------- END OF CODE --------------
