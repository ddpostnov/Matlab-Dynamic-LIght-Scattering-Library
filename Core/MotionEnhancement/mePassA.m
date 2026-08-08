%mePassA - Read the recording once and write the small working copies.
%
%   THE RECORDING IS OPENED EXACTLY ONCE, EVER.  Bio-Formats indexes a .cxd of
%   this size before it will return a single plane - measured at 1005 s on the
%   64 GB reference recording, which is 17 minutes of nothing happening.  That
%   cost is paid here and never again: everything downstream reads the flat
%   working copies this writes, through meReadRaw, at local disk speed.
%
%   AND IT IS READ A PLANE AT A TIME.  readCXD still allocates the whole stack
%   before the first frame arrives, which for this recording is 129 GB as single,
%   so it cannot be pointed at one.  The loop below holds one frame, accumulates
%   into images the size of one frame, and writes through a buffer of a few
%   hundred.  runIntensity was rewritten the same way and streams too, so it is no
%   longer the cautionary example this line used to name.
%
%   THE INTENSITY SCALING IS FROZEN BEFORE THE FIRST FRAME IS WRITTEN.  A short
%   pre-pass over frames spread across the recording sets the two intensity limits
%   that map to the bottom and the top of the sample, and they are then constant
%   for every frame and recorded beside the data.  Limits that tracked the frame
%   would turn the brightness itself into a signal that was never in the recording
%   - which, in a block whose whole purpose is to amplify small changes, is the
%   worst available bug.
%
%   THE PRE-PASS ALSO DECIDES WHETHER EIGHT BITS ARE ENOUGH, AND IT DECIDES IT BY
%   MEASUREMENT.  Rounding to eight bits adds a noise of 0.289 of a step.  If the
%   noise already in the recording is at least one step after scaling, the total
%   grows by 4% and eight bits are free; below that, rounding becomes the largest
%   noise source at exactly the small amplitudes this block exists to resolve, and
%   the working copies are written at sixteen bits instead.  The measurement, the
%   verdict and the resulting sample size are all written to the sidecar, and
%   meReadRaw reads the sample size from there rather than guessing it.
%
%   IT PRINTS WHILE IT RUNS.  This is a twenty-minute read, so it says how far it
%   has got and how long is left, every so many frames.  That is a deliberate
%   difference from the library's wrappers, which narrate in three lines and stay
%   silent in between: three lines is right for a sixty-file batch and wrong for a
%   single read that outlives the operator's patience.
%
% Syntax:
%    pass = mePassA(s, fName)
%
% Inputs:
%    s     - settings from meSettings.  Fields read: outFolder, binFactor,
%            cropSize, roi, pctClip, sampleFrames, blockFrames, blockStart,
%            bleachBlock, writeBlock, reportEvery, vesselPct, pixelSize.
%    fName - path to the recording.
%
% Outputs:
%    pass - the sidecar struct, also saved as <recording>_ME_passA.mat.  It
%           carries the settings, the frozen intensity limits, the sample size of
%           the working copies, the noise measurement behind that choice, the
%           window the full-resolution copy covers, the acquisition metadata as
%           text, and the timings.
%
%           Written beside it, in the same folder:
%             <recording>_ME_trace.mat - whole-frame and vessel mean intensity
%                                        per frame, and time
%             <recording>_ME_mean.mat  - time-mean frame, the mean of the first
%                                        and of the last minute, the vessel mask,
%                                        per-frame minimum and maximum, and the
%                                        bleaching curve
%             <recording>_ME_bin2.raw  - the binned copy of the whole field
%             <recording>_ME_crop.raw  - the full-resolution copy of one window
%
% Example:
%    s = meSettings;
%    pass = mePassA(s,'C:\Data\CXD\EPFL_20241114_2ADTF13BP.cxd');
%
% Dependencies: the Bio-Formats MATLAB library (bfGetReader, bfGetPlane);
%               Image Processing Toolbox (imgradient); Statistics and Machine
%               Learning Toolbox (prctile).
% See also: meSettings, meReadRaw, meProbe
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function pass = mePassA(s, fName)

tCall = tic;
fName = char(fName);
[recFolder,stem] = fileparts(fName);

outFolder = s.outFolder;
if isempty(outFolder), outFolder = fullfile(recFolder,'motionEnhancement'); end
if ~isfolder(outFolder), mkdir(outFolder); end

fprintf('Starting the single read of %s\n', stem);
fprintf('  working copies go to %s\n', outFolder);

% ---- open, and pay the index once ------------------------------------------
fprintf('  indexing the recording, expect about 17 minutes\n');
tOpen  = tic;
reader = bfGetReader(fName);
tOpen  = toc(tOpen);
guard  = onCleanup(@() closeReader(reader));
fprintf('  indexed in %s\n', hms(tOpen));

omeMeta = reader.getMetadataStore();
sizeX   = double(omeMeta.getPixelsSizeX(0).getValue());
sizeY   = double(omeMeta.getPixelsSizeY(0).getValue());
sizeT   = double(reader.getImageCount());
dt      = double(omeMeta.getPixelsTimeIncrement(0).value()).*1000;   % ms
fs      = 1000/dt;
startT  = char(omeMeta.getImageAcquisitionDate(0).getValue());
little  = reader.isLittleEndian();
fprintf('  %d x %d, %d frames, %.4f ms per frame (%.2f Hz), acquired %s\n', ...
    sizeX, sizeY, sizeT, dt, fs, startT);

% The plane reader below is the raw byte call rather than bfGetPlane, which runs
% an inputParser and five reader queries on every one of thirty thousand frames.
% It is checked against bfGetPlane once, here, so the shortcut is a measured
% equivalence and not an assumption.
if ~isequal(mePlane(reader,1,sizeX,sizeY,little), bfGetPlane(reader,1))
    error('mePassA:planeMismatch','The fast plane reader disagrees with bfGetPlane.');
end

% ---- the pre-pass: intensity limits, noise, and where to crop ---------------
fprintf('  sampling %d frames for the intensity limits and %d for the noise\n', ...
    s.sampleFrames, s.blockFrames);
tPre = tic;
pre  = prePass(reader, s, sizeX, sizeY, sizeT, fs, little);
tPre = toc(tPre);
fprintf('  sampled in %s\n', hms(tPre));

lo = pre.lo; hi = pre.hi;
fprintf('  the %.4g to %.4g intensity percentiles are %.0f and %.0f counts\n', ...
    s.pctClip(1), s.pctClip(2), lo, hi);

% ---- the eight-bit gate, decided here and recorded -------------------------
gain8    = 255/(hi-lo);
sigma8   = pre.sigmaBack*gain8;
useEight = sigma8 >= 1;
if useEight
    rawClass = 'uint8';  rawMax = 255;  nByte = 1;
else
    rawClass = 'uint16'; rawMax = 65535; nByte = 2;
end
fprintf('  a still, dim patch carries %.3g counts of noise, which is %.2f steps of an eight-bit sample\n', ...
    pre.sigmaBack, sigma8);
fprintf('  a still, bright patch carries %.3g counts, which is %.2f steps\n', ...
    pre.sigmaBright, pre.sigmaBright*gain8);
if useEight
    fprintf('  eight bits are enough: rounding adds 0.289 steps, so the total noise grows by %.1f%%\n', ...
        100*(sqrt(1+(0.289/sigma8)^2)-1));
else
    fprintf('  eight bits are NOT enough here, so the working copies are written at sixteen\n');
end
gain = rawMax/(hi-lo);

% ---- the window the full-resolution copy covers -----------------------------
if isempty(s.roi)
    roi = pre.roi;
    fprintf('  the full-resolution copy follows the most pulsatile window\n');
else
    roi = s.roi;
end
roi(1,:) = min(max(roi(1,:),1),sizeY);
roi(2,:) = min(max(roi(2,:),1),sizeX);
cH = roi(1,2)-roi(1,1)+1;  cW = roi(2,2)-roi(2,1)+1;
fprintf('  it covers rows %d to %d and columns %d to %d, %d by %d pixels\n', ...
    roi(1,1),roi(1,2),roi(2,1),roi(2,2), cH, cW);

% ---- geometry of the binned copy --------------------------------------------
b   = s.binFactor;
nyB = floor(sizeY/b);  nxB = floor(sizeX/b);
useY = 1:(nyB*b);      useX = 1:(nxB*b);
fprintf('  the binned copy is %d by %d and %.1f GB; the full-resolution copy is %.1f GB\n', ...
    nyB, nxB, nyB*nxB*sizeT*nByte/2^30, cH*cW*sizeT*nByte/2^30);

% ---- the one pass ------------------------------------------------------------
binPath  = fullfile(outFolder,[stem '_ME_bin' num2str(b) '.raw']);
cropPath = fullfile(outFolder,[stem '_ME_crop.raw']);
fidB = fopen(binPath ,'w');  if fidB<0, error('mePassA:write','Cannot write %s',binPath);  end
fidC = fopen(cropPath,'w');  if fidC<0, error('mePassA:write','Cannot write %s',cropPath); end
closeFiles = onCleanup(@() closeBoth(fidB,fidC));

nb   = s.writeBlock;
bufB = zeros(nyB,nxB,nb,rawClass);
bufC = zeros(cH ,cW ,nb,rawClass);

meanFrame = zeros(sizeY,sizeX);          % double, so the sum over 31 000 frames stays exact
blockSum  = zeros(sizeY,sizeX,'single');
nMinute   = min(round(60000/dt), floor(sizeT/2));
sumFirst  = zeros(sizeY,sizeX,'single');
sumLast   = zeros(sizeY,sizeX,'single');

traceAll  = zeros(sizeT,1);
traceVes  = zeros(sizeT,1);
frameMin  = zeros(sizeT,1,'uint16');
frameMax  = zeros(sizeT,1,'uint16');

vesIdx = pre.vesselIdx;
nPix   = sizeY*sizeX;
k      = 0;
tMain  = tic;

for t = 1:sizeT
    I = mePlane(reader,t,sizeX,sizeY,little);
    A = single(I);

    % ---- what the whole frame says ----
    traceAll(t) = sum(A(:))/nPix;
    traceVes(t) = mean(A(vesIdx));
    frameMin(t) = min(I(:));
    frameMax(t) = max(I(:));
    blockSum    = blockSum + A;
    if t<=nMinute,        sumFirst = sumFirst + A; end
    if t> sizeT-nMinute,  sumLast  = sumLast  + A; end

    % ---- the two working copies ----
    Ab = A(useY,useX);
    Ab = (Ab(1:b:end,1:b:end) + Ab(2:b:end,1:b:end) + Ab(1:b:end,2:b:end) + Ab(2:b:end,2:b:end)).*0.25;
    Ac = A(roi(1,1):roi(1,2), roi(2,1):roi(2,2));

    k = k+1;
    if useEight
        bufB(:,:,k) = uint8((Ab-lo).*gain);
        bufC(:,:,k) = uint8((Ac-lo).*gain);
    else
        bufB(:,:,k) = uint16((Ab-lo).*gain);
        bufC(:,:,k) = uint16((Ac-lo).*gain);
    end

    if k==nb || t==sizeT
        fwrite(fidB, bufB(:,:,1:k), rawClass);
        fwrite(fidC, bufC(:,:,1:k), rawClass);
        meanFrame   = meanFrame + double(blockSum);
        blockSum(:) = 0;
        k = 0;
    end

    if mod(t,s.reportEvery)==0
        el = toc(tMain);
        fprintf('  %6d of %d frames, %4.1f%%, %5.1f frames per second, %s elapsed, %s left\n', ...
            t, sizeT, 100*t/sizeT, t/el, hms(el), hms(el*(sizeT-t)/t));
    end
end
tMainDone = toc(tMain);
fprintf('  read %d frames in %s, %.1f frames per second\n', sizeT, hms(tMainDone), sizeT/tMainDone);

meanFrame = single(meanFrame./sizeT);
meanFirst = sumFirst./nMinute;
meanLast  = sumLast ./nMinute;

% ---- the bleaching curve, and the two small files ---------------------------
nBlk       = floor(sizeT/s.bleachBlock);
bleach     = zeros(nBlk,1);
bleachTime = zeros(nBlk,1);
for i = 1:nBlk
    idx           = (i-1)*s.bleachBlock + (1:s.bleachBlock);
    bleach(i)     = mean(traceAll(idx));
    bleachTime(i) = mean(idx-1).*dt./1000;
end

time = (0:sizeT-1)'.*dt./1000;
save(fullfile(outFolder,[stem '_ME_trace.mat']), ...
    'traceAll','traceVes','time','dt','-v7');
vesselMask = pre.vesselMask;
save(fullfile(outFolder,[stem '_ME_mean.mat']), ...
    'meanFrame','meanFirst','meanLast','vesselMask','frameMin','frameMax', ...
    'bleach','bleachTime','nMinute','-v7.3','-nocompression');

% ---- the sidecar --------------------------------------------------------------
pass                = struct();
pass.settings       = s;
pass.recording      = fName;
pass.stem           = stem;
pass.outFolder      = outFolder;
pass.sizeX          = sizeX;
pass.sizeY          = sizeY;
pass.sizeT          = sizeT;
pass.dt             = dt;
pass.fs             = fs;
pass.startT         = startT;
pass.pixelSize      = s.pixelSize;
pass.lo             = lo;
pass.hi             = hi;
pass.gain           = gain;
pass.rawClass       = rawClass;
pass.eightBitSteps  = sigma8;
pass.sigmaBack      = pre.sigmaBack;
pass.sigmaBright    = pre.sigmaBright;
pass.meanBack       = pre.meanBack;
pass.meanBright     = pre.meanBright;
pass.backBox        = pre.backBox;
pass.brightBox      = pre.brightBox;
pass.binFactor      = b;
pass.binSize        = [nyB nxB];
pass.binPath        = binPath;
pass.roi            = roi;
pass.cropSize       = [cH cW];
pass.cropPath       = cropPath;
pass.tracePath      = fullfile(outFolder,[stem '_ME_trace.mat']);
pass.meanPath       = fullfile(outFolder,[stem '_ME_mean.mat']);
pass.metadata       = readMetadata(reader);
pass.secondsIndex   = tOpen;
pass.secondsPrePass = tPre;
pass.secondsRead    = tMainDone;
pass.secondsTotal   = toc(tCall);

save(fullfile(outFolder,[stem '_ME_passA.mat']),'pass','-v7');
fprintf('Finished the single read of %s, time elapsed: %s\n', stem, hms(pass.secondsTotal));
end

% =====================================================================
function pre = prePass(reader, s, sizeX, sizeY, sizeT, fs, little)
%prePass  Everything that has to be known before the first frame is written.
%   Three questions, one set of reads: what intensity range to map onto the
%   sample, how much noise a still part of the field carries, and which window
%   carries the most pulsatile signal.  The noise and the window both need
%   CONSECUTIVE frames - a strided sample would measure the animal rather than
%   the camera - so the two reads are separate on purpose.

b    = s.binFactor;
nyB  = floor(sizeY/b);  nxB = floor(sizeX/b);
useY = 1:(nyB*b);       useX = 1:(nxB*b);

% ---- frames spread over the whole recording: the intensity limits ----------
idx = unique(round(linspace(1,sizeT,s.sampleFrames)));
sub = zeros(numel(1:4:sizeY)*numel(1:4:sizeX), numel(idx),'single');
sampleMean = zeros(sizeY,sizeX,'single');
for i = 1:numel(idx)
    A = single(mePlane(reader,idx(i),sizeX,sizeY,little));
    p = A(1:4:end,1:4:end);
    sub(:,i)   = p(:);
    sampleMean = sampleMean + A;
end
sampleMean = sampleMean./numel(idx);
lim    = prctile(double(sub(:)), s.pctClip);
pre.lo = lim(1);
pre.hi = lim(2);

% The vessel trace follows the bright pixels, because the dye is in the plasma
% and the lumen is what it fills.
thr            = prctile(double(sampleMean(:)), s.vesselPct);
pre.vesselMask = sampleMean > thr;
pre.vesselIdx  = find(pre.vesselMask);

% ---- consecutive frames: the noise, and where the motion is ----------------
t0     = max(1, min(sizeT-s.blockFrames, round(s.blockStart*sizeT)));
sumD2  = zeros(sizeY,sizeX);            % sum of squared frame-to-frame steps
blockM = zeros(sizeY,sizeX,'single');
blockB = zeros(nyB,nxB,s.blockFrames,'single');
prev   = [];
for i = 1:s.blockFrames
    A      = single(mePlane(reader,t0+i-1,sizeX,sizeY,little));
    blockM = blockM + A;
    Ab     = A(useY,useX);
    blockB(:,:,i) = (Ab(1:b:end,1:b:end) + Ab(2:b:end,1:b:end) + ...
                     Ab(1:b:end,2:b:end) + Ab(2:b:end,2:b:end)).*0.25;
    if ~isempty(prev)
        sumD2 = sumD2 + double(A - prev).^2;
    end
    prev = A;
end
blockM = blockM./s.blockFrames;

% A step between neighbouring frames holds the noise of both, so half its
% variance is the noise of one.  Anything slower than the frame rate - bleaching,
% the animal's own slow drift - cancels in the difference, which is why this is
% the estimator rather than a plain standard deviation over the block.
sigmaPix = sqrt(sumD2./(2*(s.blockFrames-1)));

% ---- the still patch the eight-bit gate is decided on ----------------------
% Still means two things at once: dim, so it is background rather than vessel,
% and flat, so that whatever the animal does moves no intensity across it.  A
% patch chosen on dimness alone can sit on the shoulder of a vessel and report
% the animal's pulsation as camera noise.
w    = 64;
grad = imgradient(blockM);
mBox = boxMean(double(blockM),w);
gBox = boxMean(double(grad  ),w);

[pre.backBox,   pre.sigmaBack,   pre.meanBack  ] = pickPatch( zscore0(mBox)+zscore0(gBox), w, sigmaPix, blockM);
[pre.brightBox, pre.sigmaBright, pre.meanBright] = pickPatch(-zscore0(mBox)+zscore0(gBox), w, sigmaPix, blockM);

% ---- the window the full-resolution copy will cover ------------------------
% What is wanted is the window carrying the most fast, in-band change, so the
% block is high-passed by subtracting a quarter-second running mean and scored on
% what is left.  The vessel that moves visibly scores highest, which is the one
% the probe most wants at full resolution.
blockB = blockB - movmean(blockB, max(3,round(0.25*fs)), 3);
score  = double(std(blockB,0,3));
cH  = min(s.cropSize(1), nyB*b);
cW  = min(s.cropSize(2), nxB*b);
wB  = [floor(cH/b) floor(cW/b)];
box = conv2(ones(wB(1),1), ones(1,wB(2)), score, 'valid');
[~,ib]  = max(box(:));
[rB,cB] = ind2sub(size(box),ib);
r0 = min(max((rB-1)*b + 1, 1), sizeY-cH+1);
c0 = min(max((cB-1)*b + 1, 1), sizeX-cW+1);
pre.roi = [r0 r0+cH-1; c0 c0+cW-1];
end

% =====================================================================
function I = mePlane(reader, iPlane, sizeX, sizeY, little)
%mePlane  One frame, without bfGetPlane's per-call inputParser.
%   Bio-Formats hands back the plane as raw bytes in reading order, so the sample
%   for column x of row y sits at (y-1)*sizeX+x: reshape to [x y] and transpose.
%   The byte order is the recording's, not the machine's, so a big-endian file is
%   swapped here - the same question bfGetPlane asks the reader and answers with
%   DataTools.  mePassA checks this against bfGetPlane once, on the first frame.
I = typecast(reader.openBytes(iPlane-1),'uint16');
if ~little, I = swapbytes(I); end
I = reshape(I, sizeX, sizeY).';
end

% =====================================================================
function [box,sigma,level] = pickPatch(cost, w, sigmaPix, meanImg)
%pickPatch  The lowest-cost w-by-w window, and what it measures.
%   The median rather than the mean over the window, so a single hot pixel or a
%   vessel corner that crept into the box does not set the answer.
[~,i]   = min(cost(:));
[r1,c1] = ind2sub(size(cost),i);
box     = [r1 r1+w-1; c1 c1+w-1];
sigma   = median(sigmaPix(r1:r1+w-1, c1:c1+w-1),'all');
level   = median(double(meanImg(r1:r1+w-1, c1:c1+w-1)),'all');
end

% =====================================================================
function m = boxMean(A,w)
%boxMean  Mean of A over every w-by-w window, one value per top-left corner.
%   Separable, because a full two-dimensional convolution with a 64-by-64 kernel
%   over a megapixel is four thousand times the work for the same answer.
m = conv2(ones(w,1)./w, ones(1,w)./w, A, 'valid');
end

% =====================================================================
function z = zscore0(A)
%zscore0  Standardise, tolerating a constant input.
sd = std(A(:));
if sd==0, z = zeros(size(A)); else, z = (A-mean(A(:)))./sd; end
end

% =====================================================================
function txt = readMetadata(reader)
%readMetadata  The acquisition metadata as lines of text.
%   Kept verbatim so that questions nobody has asked yet - the readout time
%   behind a rolling shutter is this session's - can be answered from the sidecar
%   without opening the recording again.
txt = cell(0,1);
tabs = {};
try
    tabs = {reader.getGlobalMetadata(), reader.getSeriesMetadata()};
catch
end
for i = 1:numel(tabs)
    h = tabs{i};
    if isempty(h), continue; end
    ks = h.keySet().toArray();
    for j = 1:numel(ks)
        try
            txt{end+1,1} = sprintf('%s = %s', char(ks(j)), char(h.get(ks(j)).toString())); %#ok<AGROW>
        catch
        end
    end
end
if ~isempty(txt), txt = unique(txt); end
end

% =====================================================================
function s = hms(sec)
%hms  Seconds as minutes and seconds, for a line an operator reads while waiting.
s = sprintf('%d min %02d s', floor(sec/60), round(mod(sec,60)));
end

% =====================================================================
function closeReader(reader)
try
    reader.close();
catch
end
end

% =====================================================================
function closeBoth(a,b)
try
    fclose(a);
catch
end
try
    fclose(b);
catch
end
end
%------------- END OF CODE --------------
