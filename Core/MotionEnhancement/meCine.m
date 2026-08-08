%meCine - One averaged cardiac cycle, built out of every accepted beat.
%
%   THE PRODUCT THE WHOLE PACKAGE EXISTS TO MAGNIFY.  Amplified MRI does not amplify
%   a raw time series; it amplifies a retrospectively gated cine, so the noise is
%   already down by the root of the number of beats before the magnifier runs.  This
%   is that cine, and on this recording it is not one of two options - Session 0
%   measured the same walls both ways and a band-pass at the heart rate returned
%   1.08 to 1.25 times an empty band on every vessel, while the fold cleared its own
%   floor by 4.6 to 7.6.
%
%   AND IT IS WHAT MAKES THE FULL FIELD AFFORDABLE.  30 975 frames collapse to
%   twenty-five, so the Riesz pyramid over the whole field costs under a gigabyte
%   instead of six hundred and forty, and the expensive path becomes the cheap one.
%   The continuous mode is the one that needs a crop and a time window.
%
%   STRETCHED MIN TO MIN, NEVER WINDOWED AT A FIXED LENGTH.  The heart rate wanders -
%   11.2 to 11.8 Hz interquartile among the accepted cycles here - so a window of
%   constant length would smear the late part of a short beat across a phase bin or
%   two, which flattens exactly the peak being measured.  Each cycle is resampled
%   onto the same phase axis instead, and the fiducials meGate hands over are refined
%   to sub-frame precision, which is what lets twenty-five bins be filled by eight
%   frames a beat without holes.
%
%   EACH CYCLE IS CENTRED ON ITS OWN MEAN, PER PIXEL, so bleaching and any slow drift
%   contribute nothing to the shape.  What is averaged is the SHAPE of the beat; the
%   picture it is drawn on is the mean of the same frames, accumulated alongside, so
%   the cine is a photograph rather than a difference image and nothing outside it
%   has to be loaded to make one.
%
%   THE SEM RIDES ALONG, AND IT IS NOT A GARNISH.  A magnified movie is persuasive
%   whether or not it is true.  The per-bin standard error says which part of the
%   cycle is determined and which is a handful of beats disagreeing, and it is the
%   only thing on the frame that can contradict the movie.  It must reach the report
%   and the video.
%
%   A CYCLE IS IN OR OUT, NEVER PARTLY IN.  meGate has already rejected the cycles
%   the animal moved through, and this function does not second-guess it per bin.
%   Dropping single (cycle, bin) samples would give neighbouring bins different
%   populations of beats, and then the differences BETWEEN beats appear in the cine
%   as structure WITHIN one - which is indistinguishable from the wall motion being
%   looked for.  One population, one beat.
%
%   AND THE MATCHED CONTROL IS BUILT HERE TOO, BECAUSE IT HAS TO BE THE SAME FOLD.
%   s.null = 'shuffle' gives every accepted cycle its own random circular time shift
%   before it is resampled.  Every bin then still receives exactly one interpolated
%   sample from every cycle - same frames, same interpolation weights, same photon
%   noise, same texture - and only the phase alignment is destroyed, so the rhythm
%   falls by the root of the number of cycles and nothing else changes.  That is a
%   permutation of the frames in time, done where it costs nothing: permuting the
%   30 975 frames on disk instead would be a random-access read of eight gigabytes
%   to arrive at the same statistics.
%
% Syntax:
%    cn = meCine(s, X, fs, gate)
%
% Inputs:
%    s    - settings.  Fields read:
%           .nPhaseBins  phase bins one cycle is averaged into.
%           .cineHarmonics  carried through to cn.passband, which is what meShift
%                        should filter the cine with.
%           .null        'shuffle' builds the MATCHED CONTROL instead of the cine -
%                        see below.  '' builds the cine.
%           .nullSeed    the control's own random stream.
%    X    - the frames.  Either a memmapfile whose .Data.d is [rows columns frames],
%           which is how the working copies are reached, or a numeric stack of the
%           same shape.
%    fs   - frames per second.
%    gate - the struct from meGate.  .acceptedTime is what is folded on.
%
% Outputs:
%    cn - .cine       [rows columns nBin] single, one averaged cycle as a picture
%         .shape      the same with the mean frame taken out, which is the motion
%         .sem        [rows columns nBin] single, standard error of the mean per bin
%         .base       [rows columns] the mean of the frames that went in
%         .n          [nBin 1] cycles behind each bin.  A vector because a reader
%                     should be able to see that it is constant rather than be told.
%         .nCycles    how many cycles were folded
%         .nOffered   how many meGate accepted, which is larger only if a cycle ran
%                     off the end of the recording
%         .cycleSeconds  median duration of the cycles that went in
%         .fsBins     bins per cycle - the sampling rate of the cine's own axis
%         .null       '' for the cine, 'shuffle' when this is the matched control
%         .passband   [low high] in cycles per cycle, for meShift
%         .phaseFill  .counts .emptyBins .minOverMean .maxOverMean - the confirmation
%                     that the sampled cardiac phase fills in, on this gate
%         .seconds    what it cost
%
% Example:
%    [m,geom] = meReadRaw(pass,'bin');
%    gate     = meGate(s, trace, geom.fs, ep.bad);
%    cn       = meCine(s, m, geom.fs, gate);
%    s.passband = cn.passband;
%    mag = meShift(cn.cine, cn.fsBins, s);
%
% Dependencies: none beyond core MATLAB.
% See also: meGate, meEpochs, meEnhance, meShift, meReadRaw, meProbe
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function cn = meCine(s, X, fs, gate)

tCall = tic;
[nr,nc,nT] = stackSize(X);
nPix = nr*nc;
nBin = s.nPhaseBins;
if nBin < 4
    error('meCine:bins','A cycle described by %d bins is not a cycle.',nBin);
end

cyc = gate.acceptedTime;
if isempty(cyc)
    error('meCine:noCycles', ...
        ['meGate accepted no cycle at all, so there is nothing to average. On this ' ...
         'recording that is what happens to the vasomotion band: a cycle lasts ' ...
         'several seconds and the animal moves every couple of them.']);
end

p    = (((0:nBin-1)'+0.5)./nBin);
acc  = zeros(nPix, nBin);
acc2 = zeros(nPix, nBin);
accB = zeros(nPix, 1);
n    = 0;

% The control's shifts come from a stream of its own, so a control is reproducible
% and running one does not move any other random number a session draws.
isNull = strcmpi(s.null,'shuffle');
if isNull
    rs = RandStream('threefry','Seed',s.nullSeed);
    u  = rand(rs, size(cyc,1), 1);
else
    u  = zeros(size(cyc,1),1);
end

for k = 1:size(cyc,1)
    t0 = cyc(k,1);  t1 = cyc(k,3);
    tgt = t0 + mod(p+u(k),1).*(t1-t0);

    % One block per cycle, with a frame of margin either side: the refined fiducial
    % can sit up to half a frame before the first sample of the cycle.
    i0 = max(1,  floor(t0*fs)+1 - 1);
    i1 = min(nT, ceil(t1*fs)+1 + 1);
    if i1-i0 < 3, continue; end
    tB = (i0-1:i1-1)'./fs;
    if min(tgt) < tB(1) || max(tgt) > tB(end), continue; end     % tgt is not sorted


    B = single(readBlock(X, i0, i1));
    W = interpWeights(tB, tgt);                    % [nF x nBin], two non-zeros a column
    V = reshape(B, nPix, []) * W;                  % [nPix x nBin]

    b = mean(V,2);
    V = V - b;                                     % this cycle's shape, per pixel
    acc  = acc  + double(V);
    acc2 = acc2 + double(V).^2;
    accB = accB + double(b);
    n    = n + 1;
end

if n < 2
    error('meCine:tooFew','Only %d cycle could be folded; a mean of one is that one.',n);
end

shape = acc./n;
sd    = sqrt(max(acc2./n - shape.^2, 0));
base  = accB./n;

cn = struct();
cn.cine    = single(reshape(base + shape, nr, nc, nBin));
cn.shape   = single(reshape(shape,        nr, nc, nBin));
cn.sem     = single(reshape(sd./sqrt(n),  nr, nc, nBin));
cn.base    = single(reshape(base,         nr, nc));
cn.n       = repmat(n, nBin, 1);
cn.nCycles = n;
cn.nOffered= size(cyc,1);
cn.cycleSeconds = median(cyc(:,3)-cyc(:,1));
cn.fsBins  = nBin;
cn.null    = s.null;

% The filter that belongs on a cine keeps the fundamental and the harmonics the frame
% rate can carry, and drops the static picture.  In cycles per cycle, so the band is
% the same numbers whatever the heart rate was.
cn.passband = [0.5, min(s.cineHarmonics, floor(nBin/2)-0.5) + 0.5];

% ---- the confirmation that the phase axis fills ------------------------------
% Session 0 passed this on all 3 168 beats with the fiducial snapped to a frame, and
% reported a fullest bin holding five times the average - which was the snapping, not
% the recording.  On the refined fiducials that pile-up is gone, and this is where a
% reader sees whether it is gone on THEIR recording.
ph = gate.phase(~isnan(gate.phase));
counts = histcounts(ph, linspace(0,1,nBin+1));
cn.phaseFill = struct('counts',counts(:),'nBin',nBin,'nSample',numel(ph), ...
    'emptyBins',sum(counts==0),'minOverMean',min(counts)/mean(counts), ...
    'maxOverMean',max(counts)/mean(counts));

cn.seconds = toc(tCall);
end

% =====================================================================
function [nr,nc,nT] = stackSize(X)
%stackSize  The shape of either kind of input, without loading any of it.
if isa(X,'memmapfile')
    sz = X.Format{2};
    nr = sz(1); nc = sz(2); nT = sz(3);
else
    nr = size(X,1); nc = size(X,2); nT = size(X,3);
end
end

% =====================================================================
function B = readBlock(X, i0, i1)
%readBlock  Frames i0 to i1 of either kind of input.
if isa(X,'memmapfile')
    B = X.Data.d(:,:,i0:i1);
else
    B = X(:,:,i0:i1);
end
end

% =====================================================================
function W = interpWeights(tB, tgt)
%interpWeights  Linear interpolation onto the phase axis, as a matrix.
%   Two non-zeros per column, so one matrix product resamples every pixel of the
%   block at once.  Doing it with interp1 per pixel is the same arithmetic through a
%   quarter of a million function calls.
nF = numel(tB);  nB = numel(tgt);
dt = tB(2)-tB(1);
u  = (tgt - tB(1))./dt;                    % fractional index, zero based
lo = min(max(floor(u), 0), nF-2);
w  = u - lo;
W  = zeros(nF, nB, 'single');
W(sub2ind([nF nB], lo+1,   (1:nB)')) = single(1-w);
W(sub2ind([nF nB], lo+2,   (1:nB)')) = single(w);
end
%------------- END OF CODE --------------
