% meRun - Motion enhancement of vascular wall motion: the stand-alone entry point.
%
%   Run STEP 0 once per MATLAB session, then the step cells (%%) in order.  Every
%   number in here was measured on the reference recording in Sessions 0 to 3 and is
%   commented with what measured it; edit the recording and the two calibres and the
%   rest follows.
%
%   THE BLOCK IS STAND-ALONE.  It calls the library; the library never calls it.
%   Nothing here is registered with a launcher, a GUI, the workbench step registry or
%   the explorer schema, and the one file it writes beside the videos is
%   <recording>_ME.mat - NOT a _r.mat/_s.mat triplet, because those suffixes are
%   swept by the workbench and the explorer and this is not a library product.
%
%   THE RECORDING IS NEVER OPENED HERE.  mePassA read the CXD once, at a cost of
%   16 min 43 s of Bio-Formats indexing before the first plane, and wrote two flat
%   working copies.  Everything below reads those.
%
%   EVERY PRODUCT SHIPS WITH ITS NULL, and the null is the same run with the timing
%   destroyed and nothing else changed - same frames, same texture, same photon
%   noise, same alpha.  A magnified movie is persuasive whether or not it is true, so
%   the deliverable is never the movie: it is the movie beside the number that says
%   how much of it is real, and that number is the ratio between the product and its
%   own null.
%
%   ONE ALPHA DOES NOT SERVE THE WHOLE FIELD.  The output is proportional to the
%   displacement, so a run is calibrated for ONE calibre of vessel: set alpha from a
%   0.075 px wall and a 0.011 px wall reaches a quarter of a pixel; set it from the
%   small one and the large one is seven times past its bound and deforms.  The
%   cardiac cine is therefore delivered TWICE, at two calibres, and each says on its
%   own frames which vessels it is valid for.
%
% Example:
%    edit meRun            % and press Run Section, cell by cell
%
% Dependencies: meSettings, meEnhance, meWriteVideo, meKymograph.
% See also: mePassA, meProbe, meEnhance, meWriteVideo, meCine, meGate
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library

%% STEP 1 The recording, the shared settings, and where the videos go
close all
clearvars -except libraryFolder
dataFolder = 'C:\Dropbox\Work\Data\CXD\motionEnhancement'; %where mePassA put the working copies
stem       = 'EPFL_20241114_2ADTF13BP';
passPath   = fullfile(dataFolder,[stem '_ME_passA.mat']);
outFolder  = dataFolder; %WHERE THE VIDEOS GO. They come to about a gigabyte - point this off a synchronised folder if that matters

pass  = load(passPath,'pass').pass;
probe = load(fullfile(dataFolder,[stem '_ME_probe.mat']),'probe').probe;
cuts  = probe.cuts; %the five vessel cuts meProbe measured, reused here to check every product

%SHARED - measured in Sessions 0 to 2, not guessed
s0 = meSettings;
s0.levels    = 3:5;      %IN LEVELS OF THE RECORDING - lambda 16, 32 and 64 px. The wall edge is 6 to 15 px wide, so level 2 sits below the structure. meEnhance maps this onto whichever copy it reads
s0.riesz     = 'exact';  %24% more amplification than the published three taps at a floor that cannot be told apart, and it does not wrap at the frame edge
s0.blurSigma = 2;        %the minimum of the detection floor; sigma 4 is worse on BOTH gain and floor
s0.maskMode  = 'weight'; %confine amplification to the vessel walls: 287x in contrast for 31% of the gain, and the 31% comes back with alpha
s0.edgeCrop  = 16;       %the outer margin of a magnified frame is boundary response, measured against the same ring of the unmagnified cine
s0.temporalFilter = 'ideal'; %zero-phase. Acceleration magnification was tried and is WORSE here - it doubles the off-band content of the magnified wall trace
s0.showDifference = true;    %the third panel is where the magnification actually is, and the fastest way to tell a product from its null
s0.videoFps   = 25;
s0.scaleBarUm = 200;

%THE CALIBRES. alpha = lambda/(4*delta) per level, i.e. every level at its own lambda/4
%bound rather than all three held to the finest one's - worth 1.5x, one line.
%
%EVERY delta BELOW IS IN PIXELS OF THE RECORDING, and so is the [4 8 16], which is
%lambda/4 for the three amplified levels at lambda = 16, 32 and 64 recording px. alpha
%itself is scale-free - both delta and lambda halve on the binned copy - but the RATIO
%is only right if the two are in the same units, and a wall amplitude measured through
%meKymograph on the binned copy comes back in BINNED pixels and has to be doubled first.
dLarge = 0.0753;  aLarge = [4 8 16]./dLarge;  %the 38 to 75 um vessels Session 0 measured at 0.064 to 0.075 px
dSmall = 0.0107;  aSmall = [4 8 16]./dSmall;  %the 38 um vessel at 0.011 px. The large vessels are 7x past their bound in this run and say so on the frame
dVaso  = 0.356;   aVaso  = [4 8 16]./dVaso;   %an UPPER BOUND, not a detection: the largest in-band wall excursion on the raw stack (0.178 BINNED px), which is signal plus noise

me = struct('stem',stem,'recording',pass.recording,'written',datetime('now'), ...
    'products',struct([]),'gate',struct(),'notes',{{}});

%% STEP 2 The cardiac cine, calibrated for the 38 to 75 um vessels, and its null
close all
clearvars -except libraryFolder dataFolder stem passPath outFolder pass probe cuts s0 dLarge dSmall dVaso aLarge aSmall aVaso me

s = s0;
s.mode  = 'cine';    %magnify ONE averaged cardiac cycle - what amplified MRI does, and on this recording the only thing that works
s.band  = 'puls';
s.copy  = 'bin';     %the whole field, two-by-two binned. 30 975 frames collapse to 25, so the full field costs 0.1 GB
s.alpha = aLarge;
s.videoLoops = 10;   %a cine IS one cardiac period, so it is meant to be looped

[mag , info ] = meEnhance(s , passPath);
sn = s;  sn.null = 'shuffle';
[magN, infoN] = meEnhance(sn, passPath);

p = mePack(info, infoN, mag, magN, cuts, 'cardiac cine', dLarge);
p.file     = fullfile(outFolder,[stem '_ME_puls_cine.avi']);
p.nullFile = fullfile(outFolder,[stem '_ME_puls_cine_null.avi']);
meta = meMeta(stem, 'cardiac cine - one beat, averaged over every accepted cycle', s, info, p, ...
    sprintf('valid for walls moving about %.3f px - the 38 to 75 um vessels. Output on the largest: %.2f px', ...
    dLarge, p.outPx));
meta.notes = {sprintf('%d of %d cycles accepted (%.0f%%); %d phase bins, %d empty, fullest %.2fx the mean', ...
    info.cine.nCycles, size(info.gate.cycles,1), 100*info.gate.acceptanceRate, ...
    info.cine.fsBins, info.cine.phaseFill.emptyBins, info.cine.phaseFill.maxOverMean)};
p.video     = meWriteVideo(s, struct('original',info.original ,'magnified',mag ), meta);
p.nullVideo = meWriteVideo(s, struct('original',infoN.original,'magnified',magN), meNullMeta(meta, p));

me.products  = pAppend(me.products, p);
me.gate.puls = meTrimGate(info.gate);
clear mag magN

%% STEP 3 The same cine calibrated for the smallest vessels that could be measured
close all
clearvars -except libraryFolder dataFolder stem passPath outFolder pass probe cuts s0 dLarge dSmall dVaso aLarge aSmall aVaso me

s = s0;
s.mode  = 'cine';  s.band = 'puls';  s.copy = 'bin';
s.alpha = aSmall;  %seven times aLarge, because the output is proportional to the displacement
s.videoLoops = 10;

[mag , info ] = meEnhance(s , passPath);
sn = s;  sn.null = 'shuffle';
[magN, infoN] = meEnhance(sn, passPath);

p = mePack(info, infoN, mag, magN, cuts, 'cardiac cine, small vessels', dSmall);
p.file     = fullfile(outFolder,[stem '_ME_puls_cine_smallVessels.avi']);
p.nullFile = fullfile(outFolder,[stem '_ME_puls_cine_smallVessels_null.avi']);
meta = meMeta(stem, 'cardiac cine - calibrated for the smallest vessels', s, info, p, ...
    sprintf('valid for walls moving about %.3f px - vessels near 38 um. LARGER VESSELS ARE PAST THEIR BOUND HERE AND DEFORM', dSmall));
meta.notes = {'read the large vessels in the other cardiac cine, not in this one'};
p.video     = meWriteVideo(s, struct('original',info.original ,'magnified',mag ), meta);
p.nullVideo = meWriteVideo(s, struct('original',infoN.original,'magnified',magN), meNullMeta(meta, p));

me.products = pAppend(me.products, p);
clear mag magN

%% STEP 4 The cardiac continuous run - real frames in real time - and its null
close all
clearvars -except libraryFolder dataFolder stem passPath outFolder pass probe cuts s0 dLarge dSmall dVaso aLarge aSmall aVaso me

s = s0;
s.mode  = 'continuous';  %the honest view: it keeps beat-to-beat variability, which the cine averages away by construction
s.band  = 'puls';
s.copy  = 'crop';        %the full-resolution 512x512 window over the densest vasculature
s.alpha = aLarge;        %THE SAME alpha as the cine, so the only difference between the two products is the gating
s.videoLoops = 1;

[mag , info ] = meEnhance(s , passPath);
sn = s;  sn.null = 'shuffle';
[magN, infoN] = meEnhance(sn, passPath);

p = mePack(info, infoN, mag, magN, cuts, 'cardiac continuous', dLarge);
p.file     = fullfile(outFolder,[stem '_ME_puls_cont.avi']);
p.nullFile = fullfile(outFolder,[stem '_ME_puls_cont_null.avi']);
meta = meMeta(stem, 'cardiac, continuous - single beats, nothing averaged', s, info, p, ...
    'the same alpha as the cardiac cine, so the only difference between them is the gating');
meta.notes = {sprintf('%d frames from frame %d - the longest stretch this awake animal holds still for', ...
        size(mag,3), info.window(1)), ...
    'the mask says where the shift is APPLIED, not where it lands: at level 4 it spreads over 32 to 64 px through the collapse, and this crop is dense'};
p.video     = meWriteVideo(s, struct('original',info.original ,'magnified',mag ), meta);
p.nullVideo = meWriteVideo(s, struct('original',infoN.original,'magnified',magN), meNullMeta(meta, p));

me.products = pAppend(me.products, p);
clear mag magN

%% STEP 5 The vasomotion continuous run over the WHOLE recording, and its null
close all
clearvars -except libraryFolder dataFolder stem passPath outFolder pass probe cuts s0 dLarge dSmall dVaso aLarge aSmall aVaso me

s = s0;
s.mode  = 'continuous';
s.band  = 'vaso';          %moves the gate's search range onto 0.05-0.25 Hz
s.copy  = 'bin';
s.contFrames    = Inf;     %a 5.7 s cycle needs minutes of record, and 13 s is 2.3 cycles
s.contSpanBouts = true;    %so the window has to span bouts. They are MARKED, not dropped
s.decimate      = [];      %derived from the passband at six samples a cycle: 72 here, and a BLOCK MEAN - point-sampling would fold the whole cardiac band into the vasomotion one
s.alpha = aVaso;
s.videoFps   = 12.5;       %a 5.70 s vasomotion cycle then takes 0.9 s on screen
s.videoLoops = 1;

[mag , info ] = meEnhance(s , passPath);
sn = s;  sn.null = 'shuffle';
[magN, infoN] = meEnhance(sn, passPath);

p = mePack(info, infoN, mag, magN, cuts, 'vasomotion continuous', dVaso);
p.file     = fullfile(outFolder,[stem '_ME_vaso_cont.avi']);
p.nullFile = fullfile(outFolder,[stem '_ME_vaso_cont_null.avi']);
meta = meMeta(stem, 'vasomotion, continuous - the whole recording', s, info, p, ...
    sprintf('alpha set from an UPPER BOUND of %.3f px on the wall motion, not from a detection', dVaso));
meta.notes = {sprintf('%.0f%% of these frames touch a movement bout and are marked, not dropped', ...
        100*mean(info.badInWindow)), ...
    'no vasomotion wall dilation was detected above its own null on this recording - read the sidecar'};
p.video     = meWriteVideo(s, struct('original',info.original ,'magnified',mag ), meta);
p.nullVideo = meWriteVideo(s, struct('original',infoN.original,'magnified',magN), meNullMeta(meta, p));

me.products  = pAppend(me.products, p);
me.gate.vaso = meTrimGate(info.gate);
clear mag magN

%% STEP 6 The sidecar - one file, one distinct name, nothing a sweep can mistake for a product
close all
clearvars -except libraryFolder dataFolder stem passPath outFolder pass probe cuts s0 dLarge dSmall dVaso aLarge aSmall aVaso me

%THE VASOMOTION CINE IS NOT DELIVERABLE FROM THIS RECORDING and that is a measurement,
%not an omission: a vasomotion cycle lasts 5.70 s, the animal moves every 2.3 s, and
%0 of 49 cycles survive - 46 of them on the movement rule alone. meCine refuses rather
%than averaging two cycles. A calmer preparation will produce it; the code path is built.
me.notes = {
    'the vasomotion CINE is the one planned product this recording cannot give: 0 of 49 cycles'
    'accepted, 46 of them because the animal moved during them - a 5.70 s cycle against a bout'
    'every 2.3 s. meCine refuses rather than averaging two cycles.'
    'the cardiac cine is delivered TWICE, at two calibres, because one alpha serves one vessel size'
    'the gain is 1 + eta*alpha, and eta here is measured on REAL walls, not on the phantom'
    'every product ships with its null: the same run with the timing destroyed, identical alpha'
    };
save(fullfile(outFolder,[stem '_ME.mat']),'me','-v7');
fprintf('\nwrote %s\n', fullfile(outFolder,[stem '_ME.mat']));
tot = 0;
for k = 1:numel(me.products)
    q = me.products(k);
    fprintf('%-30s alpha %-20s  product/null %5.1f   on-wall %6.2f  off-wall %6.2f  eta %.3f\n', ...
        q.name, mat2str(round(q.alpha,1)), q.productOverNull, q.magOnWall, q.magOffWall, q.etaReal);
    tot = tot + q.video.bytes + q.nullVideo.bytes;
end
fprintf('%d videos, %.2f GB in %s\n', 2*numel(me.products), tot/2^30, outFolder);

%% ===================================================================
% Local functions.  They are here rather than in MotionEnhancement/ because they are
% this script's bookkeeping - packing a product, laying out a caption - and not part
% of the block's interface.
% ====================================================================
function ps = pAppend(ps, p)
%pAppend  One more product.
if isempty(ps), ps = p; else, ps(end+1) = orderfields(p, ps(1)); end
end

% =====================================================================
function p = mePack(info, infoN, mag, magN, cuts, name, dTarget)
%mePack  Everything about one product that belongs in the sidecar, and the wall
%   measurement that says how much of it is real.
%
%   THE PRODUCT IS COMPARED WITH ITS OWN NULL AT EVERY CUT, through the identical
%   estimator, so whatever bias meKymograph has cancels between them.  That ratio is
%   the deliverable; the movie is the illustration.
p = struct();
p.name       = name;
p.mode       = info.mode;      p.band       = info.band;    p.copy = info.copy;
p.alpha      = info.settings.alpha;
p.levels     = info.levels;    p.lambdaFull = info.lambdaFull;
p.passband   = info.passband;  p.fsRun      = info.fsRun;
p.decimate   = info.decimate;  p.window     = info.window;
p.edgeCrop   = info.edgeCrop;  p.pixelSize  = info.geom.pixelSize;
p.frames     = size(mag,3);
p.dTarget    = dTarget;
p.maskFraction   = mean(info.mask(:)>0.5);
p.seconds        = info.seconds;
p.nullSeconds    = infoN.seconds;
p.peakGB         = info.peakGB;
p.acceptanceRate = info.gate.acceptanceRate;
p.nCyclesOffered = size(info.gate.cycles,1);
p.settings   = rmfield(info.settings, intersect({'mask','maskVessel','maskExclude'}, ...
                                                fieldnames(info.settings)));

mk = info.mask > 0.5;
p.magOnWall  = meRegionPtp(mag  - info.original ,  mk);
p.magOffWall = meRegionPtp(mag  - info.original , ~mk);
p.nullOnWall = meRegionPtp(magN - infoN.original,  mk);
p.confine    = p.magOnWall/max(p.magOffWall,eps);

isCine = strcmpi(info.mode,'cine');
if isCine
    p.nCycles    = info.cine.nCycles;
    p.nPhaseBins = info.cine.fsBins;
    p.phaseFill  = info.cine.phaseFill;
    ec = info.edgeCrop;
    sh = info.cine.shape(ec+1:end-ec, ec+1:end-ec, :);
    se = info.cine.sem(  ec+1:end-ec, ec+1:end-ec, :);
    V  = reshape(sh, [], size(sh,3));
    E  = reshape(se, [], size(se,3));
    p.cycleOnWall    = mean(V(mk(:),:),1).';        % the strip meWriteVideo draws
    p.cycleSemWall   = mean(E(mk(:),:),1).';
    p.cycleExcursion = meRegionPtp(sh, mk);
    p.cycleSem       = meRegionMean(se, mk);
    % THE CONTROL GETS ITS OWN STRIP, over the SAME pixels, because that strip is
    % where the two videos differ most legibly: on the product the averaged cycle
    % leaves its own standard error, on the control it does not.
    shN = infoN.cine.shape(ec+1:end-ec, ec+1:end-ec, :);
    seN = infoN.cine.sem(  ec+1:end-ec, ec+1:end-ec, :);
    VN  = reshape(shN, [], size(shN,3));
    EN  = reshape(seN, [], size(seN,3));
    p.cycleOnWallNull  = mean(VN(mk(:),:),1).';
    p.cycleSemWallNull = mean(EN(mk(:),:),1).';
else
    p.nCycles = NaN;  p.nPhaseBins = NaN;  p.phaseFill = [];
    p.cycleOnWall = [];  p.cycleSemWall = [];
    p.cycleOnWallNull = [];  p.cycleSemWallNull = [];
    p.cycleExcursion = NaN;  p.cycleSem = NaN;
end

w = struct([]);  ratios = [];
for i = 1:numel(cuts)
    c = meMapCut(cuts(i), info);
    if isempty(c), continue; end
    r = meKymograph(info.original, c);
    m = meKymograph(mag,           c);
    n = meKymograph(magN,          c);
    if isCine
        aR = meCycleAmp(r.halfWidth);
        aM = meCycleAmp(m.halfWidth);
        aN = meCycleAmp(n.halfWidth);
    else
        aR = meBandAmp(r.halfWidth, info.passband, info.fsRun);
        aM = meBandAmp(m.halfWidth, info.passband, info.fsRun);
        aN = meBandAmp(n.halfWidth, info.passband, info.fsRun);
    end
    % rawPx and magPx are in pixels of THE COPY; the *Full columns are the same in
    % pixels of the recording, which is the unit every displacement in this package is
    % quoted in.  boundRatio is this run's own alpha*delta against lambda/4 at the
    % finest amplified level - taken from the run's OWN raw amplitude rather than from
    % meProbe's cardiac fold, so it means something in the vasomotion band too.
    e = struct('cut',i,'diameterUm',cuts(i).diameterUm,'foldPx',cuts(i).foldPx, ...
        'rawPx',aR,'magPx',aM,'nullPx',aN, ...
        'gain',aM/aR,'overNull',aM/max(aN,eps),'eta',(aM/aR-1)/p.alpha(1), ...
        'rawPxFull',aR*info.geom.scale,'outputPxFull',aM*info.geom.scale, ...
        'boundRatio',p.alpha(1)*aR*info.geom.scale/4, ...
        'edgeWidthRatio',mean(m.widthL+m.widthR,'omitnan')/mean(r.widthL+r.widthR,'omitnan'));
    if isempty(w), w = e; else, w(end+1) = e; end %#ok<AGROW>
    ratios(end+1) = e.overNull; %#ok<AGROW>
end
p.wall            = w;
p.productOverNull = median(ratios);
p.etaReal         = median([w.eta]);
[~,iBig]          = max([w.rawPxFull]);
p.outPx           = w(iBig).outputPxFull;
p.boundRatio      = w(iBig).boundRatio;
p.file = '';  p.nullFile = '';  p.video = [];  p.nullVideo = [];
end

% =====================================================================
function meta = meMeta(stem, product, s, info, p, calibre)
%meMeta  The caption that has to be on every frame of a product.
meta = struct('outPath',p.file,'stem',stem,'product',product, ...
    'alpha',s.alpha,'lambdaFull',info.lambdaFull,'passband',info.passband, ...
    'differential',s.globalOrder>=0,'masked',~isempty(info.mask), ...
    'calibre',calibre,'pixelSize',info.geom.pixelSize, ...
    'bad',info.badInWindow,'isNull',false,'nullText','','notes',{{}}, ...
    'controlName',nameOfFile(p.nullFile), ...
    'passbandUnit','Hz','axis','time','axisValues',[],'rateHz',info.gate.f0,'strip',[]);
if strcmpi(info.mode,'cine')
    meta.passbandUnit = 'cycles/cycle';
    meta.axis         = 'phase';
    meta.axisValues   = (((0:info.cine.fsBins-1)')+0.5)./info.cine.fsBins;
    meta.bad          = false(info.cine.fsBins,1);  %a cycle the animal moved through is not in the cine at all
    meta.strip        = struct('value',p.cycleOnWall,'sem',p.cycleSemWall, ...
        'hi', max([p.cycleOnWall+p.cycleSemWall; -(p.cycleOnWall-p.cycleSemWall)]), ...
        'label',sprintf('averaged cardiac cycle over the wall mask, counts (band = per-pixel SEM, %.1fx smaller than the excursion)', ...
                        p.cycleExcursion/p.cycleSem));
else
    meta.axisValues   = ((info.window(1)-1) + (0:p.frames-1)'.*info.decimate)./info.geom.fs;
end
end

% =====================================================================
function n = nameOfFile(p)
[~,f,e] = fileparts(p);  n = [f e];
end

% =====================================================================
function m = meNullMeta(meta, p)
%meNullMeta  The same caption, saying what it is.
%   THE READOUT LOSES ITS CLOCK, deliberately.  The frames are permuted, so a time in
%   seconds or a cardiac phase would be a caption about a recording this stack no
%   longer is; the readout counts frames instead.
m = meta;
m.outPath = p.nullFile;
m.isNull  = true;
m.nullText= 'frames permuted in time - same frames, same noise, same alpha, no rhythm';
m.product = [meta.product '   -   MATCHED CONTROL'];
m.axis    = 'index';
if ~isempty(p.cycleOnWallNull)
    m.strip = struct('value',p.cycleOnWallNull,'sem',p.cycleSemWallNull, ...
        'hi', meta.strip.hi, ...          %the PRODUCT's scale, so the two strips compare
        'label','the SAME strip on the control - what is left when the timing is destroyed');
end
end

% =====================================================================
function g = meTrimGate(gate)
%meTrimGate  The gate without its two 31 000-sample traces, for the sidecar.
g = rmfield(gate, intersect({'trace','rhythm','phase'}, fieldnames(gate)));
end

% =====================================================================
function c = meMapCut(cut, info)
%meMapCut  A probe cut, on the grid of the stack info holds.
%   meProbe's cut centres are FULL-RESOLUTION frame coordinates; a magnified stack
%   sits on a working copy's grid and has had info.edgeCrop off every side.
%   meReadRaw's map is centre to centre, so the half pixel a binned sample sits at is
%   in here - and in a block that resolves tenths of a pixel that is not a rounding
%   detail.  A cut whose span runs off the stack is dropped, not clipped.
g = info.geom;
r = (cut.center(1) - g.origin(1) - (g.scale-1)/2)/g.scale + 1 - info.edgeCrop;
k = (cut.center(2) - g.origin(2) - (g.scale-1)/2)/g.scale + 1 - info.edgeCrop;
rad  = cut.radius/g.scale;
span = 3*rad + 6;
[nr,nc,~] = size(info.original);
if r-span < 1 || r+span > nr || k-span < 1 || k+span > nc, c = []; return; end
c = struct('center',[r k],'normal',cut.normal,'radius',rad);
end

% =====================================================================
function a = meCycleAmp(v)
%meCycleAmp  Amplitude of the fundamental of a CLOSED cycle - one cycle per cycle.
%   A cine's last bin is adjacent to its first, so the natural amplitude is the first
%   Fourier coefficient rather than a peak-to-peak, which one noisy bin can carry on
%   its own.
v = v(:);  v(~isfinite(v)) = mean(v(isfinite(v)));
n = numel(v);
a = 2*abs(sum(v.*exp(-2i*pi*(0:n-1)'/n)))/n;
end

% =====================================================================
function a = meBandAmp(x, band, fs)
%meBandAmp  Amplitude of the equivalent sinusoid inside a band, zero-phase.
x = x(:);  ok = isfinite(x);
if nnz(ok) < 0.5*numel(x), a = NaN; return; end
x(~ok) = interp1(find(ok), x(ok), find(~ok), 'linear','extrap');
w = band./(fs/2);  w = min(max(w,1e-6),0.999);
[bb,aa] = butter(2, w, 'bandpass');
a = sqrt(2)*std(filtfilt(bb,aa,double(x)));
end

% =====================================================================
function v = meRegionPtp(X, mask)
%meRegionPtp  Mean over the region of EACH PIXEL's peak-to-peak over time.
%   Not the peak-to-peak of the region's mean: the two walls of one vessel move in
%   opposite directions, so a spatial mean cancels the thing being measured.
P = max(X,[],3) - min(X,[],3);
v = double(mean(P(mask)));
end

% =====================================================================
function v = meRegionMean(X, mask)
%meRegionMean  Mean over the region of each pixel's mean over time.
P = mean(X,3);
v = double(mean(P(mask)));
end
