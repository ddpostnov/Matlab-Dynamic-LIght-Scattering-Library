%runIntensityInternalCycle - Collapse a .cxd fluorescence recording to one cardiac cycle.
%
%   THE TWIN OF runContrastInternalCycle, ON THE OTHER MODALITY.  That step folds the
%   CONTRAST of an *.rls speckle acquisition onto one mean heartbeat and writes '_c_K';
%   this one folds the FLUORESCENCE of a *.cxd stack and writes '_c_I'.  Both sit on the
%   'cardiac' branch, which they can share because a recording is a .cxd OR an .rls and
%   never both - the modality filter is what separates the two steps.  The pulse
%   detection, the nine cycle features, the leave-one-out rejection and the phase
%   binning are that step's, deliberately: a second idiom for one idea would be a cost
%   with no benefit.
%
%   THE PRODUCT IS ABSOLUTE INTENSITY, AND THAT IS THE WHOLE POINT.  Each accepted cycle
%   is centred on ITS OWN PER-PIXEL MEAN before it is added in - which is what removes
%   the bleaching, and this preparation bleaches - and the MEAN IMAGE IS ADDED BACK, so
%   what comes out is a picture that beats rather than a fractional change.  A cycle
%   average left as a fractional change would silently make every downstream pulsatility
%   index a different quantity: (max-min)/mean on '_c_I' means exactly what it means on
%   '_c_BFI', so the two branches' numbers can be put in one table.  The identity that
%   says it worked is asserted on every run - the mean of the returned cycle over its
%   phase bins IS results.imgI, to floating-point.
%
%   IT READS THE RECORDING ONCE AND HOLDS ONE PLANE.  A .cxd of the reference protocol
%   is 30 975 frames of 1300 x 800, i.e. 129 GB as single, and the whole point of a cycle
%   average is that the PRODUCT is one cycle.  So nothing here ever allocates the
%   recording: pass 1 streams every plane into a spatially decimated reference cube and a
%   mean image, and pass 2 re-reads the accepted cycles in contiguous runs, THROUGH THE
%   SAME OPEN READER.  Bio-Formats spends seventeen minutes indexing a recording of that
%   size before the first plane arrives, so the one thing this step must never do is open
%   it twice.  Do not model it on readCXD, which still preallocates the stack.
%
%   TWO PASSES, BECAUSE A BEAT IS NOT KNOWN UNTIL THE TRACE IS COMPLETE.  Rejection is
%   leave-one-out against every OTHER cycle, so no cycle can be folded before the last
%   one has been seen.  Pass 1 therefore keeps a spatially decimated cube - the smallest
%   thing from which any block-level mask's reference trace can be formed after the fact
%   - and s.decimationSpace is RAISED by itself when that cube will not fit, exactly as
%   the contrast step raises its own: a coarser sampling of a spatial average is still
%   that average, so nothing the user asked for changes.  s.nPhaseBins is not like that -
%   it is the scientific content of the product - so an output that will not fit is
%   REFUSED, naming the number of bins that would have.
%
%   WHAT DRIVES THE GATE IS A VESSEL-LIMITED MEAN, NOT THE WHOLE FIELD (s.maskPctI).
%   The dye is in the plasma, so the beat is in the pixels the plasma fills, and diluting
%   them with parenchyma costs signal for nothing.  Measured on the small fixture
%   (1300 x 1300 x 668 at 33.33 Hz): the cardiac line at 10.42 Hz stands 3.71 times the
%   band median on the whole-field mean, 7.22 above the median pixel, 7.69 above the 85th
%   percentile and 4.13 above the 95th.  The optimum is INTERIOR - past the 85th there
%   are too few pixels and photon noise takes the gain back - which is what says the
%   choice is a real one and not "more masking is always better".  The percentile is used
%   rather than an absolute intensity range because it is scale-free: these recordings
%   are uint16 and their exposure is not a constant of the rig.  It is the same 85th
%   percentile the motion-enhancement block arrived at independently, on a different
%   protocol and a different animal.
%
%   AND THE TRACE IS CONDITIONED THE WAY meGate MEASURED IT HAS TO BE.  The contrast
%   step smooths with a loess span of a fraction of a cycle and takes its peak prominence
%   from the standard deviation of that smoothed trace.  Both fail on fluorescence, and
%   for reasons that belong to the signal rather than to the rule: there are only a
%   handful of samples in a beat, so the span rounds to one and nothing is smoothed; and
%   the variance of a fluorescence trace is bleaching, respiration and vasomotion rather
%   than pulse, so the threshold comes out several times the cardiac excursion.  Run that
%   way on the reference recording it found 1 167 of 3 168 beats.  So the slow decay is
%   removed as a low-order polynomial, the rhythm is taken with a zero-phase band-pass at
%   its own band, the minima are found on THAT alone - a beat on a rising level is less
%   prominent than the same beat on a falling one - and the fiducials are refined to
%   sub-frame precision, which is what lets the phase axis fill without holes.
%   s.smoothCoef1 is therefore NOT a setting here: it would be one nobody reads.
%
%   THE MATCHED CONTROL IS BUILT HERE, BECAUSE HERE IS THE ONLY PLACE IT CAN BE.  A
%   cycle average is persuasive whether or not it is true, and the only thing that can
%   contradict it is the same fold with the timing destroyed: every accepted beat given
%   its own random circular phase shift before it is resampled, so every bin still
%   receives exactly one interpolated sample from every beat - same frames, same
%   interpolation weights, same photon noise, same texture - and only the phase
%   alignment goes.  The rhythm then falls by the root of the number of beats and
%   nothing else changes, which is what makes a ratio to it a detection.
%
%   IT CANNOT BE RECONSTRUCTED AFTERWARDS AND THAT IS WHY IT IS STORED.  Once the fold
%   has run, the product IS one cycle; the beats it was made of are back in a 60 GB
%   recording that costs seventeen minutes to reopen.  A later step handed only the
%   averaged cycle can permute its bins or add noise of the right size, but neither
%   reproduces the residual this fold leaves - red cells crossing a wall and residual
%   bulk motion are spatially correlated and a surrogate built from the marginal error
%   is not.  So s.nControls of them ride in the source member beside the cycle, and
%   runMotionEnhancement refuses to measure a wall without them.  THREE, because one
%   realisation of the floor is uncertain by about a factor of two.
%
%   ONE RULE OF THE CONTRAST STEP IS ABSENT AND IT IS NAMED RATHER THAN LEFT IMPLIED.
%   There is no movement-bout rule.  Detecting a bout needs a frame-to-frame displacement
%   of the field, which is a measurement of its own and belongs with the motion step.
%   What partly covers it is rules 1, 9 and 13, which reject a cycle sitting at a
%   brightness unlike the recording's, and a bout usually changes brightness; what is NOT
%   covered is a bout that translates the field without changing its mean.  On an awake
%   preparation that is the known gap in this step.
%
% Syntax:
%    runIntensityInternalCycle(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%               .maskPctI      percentile of the mean image above which a pixel is
%                              averaged into the reference trace the beat is found in.
%               .minFrqIni .maxFrqIni  the band the heart rate is looked for in, Hz.
%               .rangeFrq      a cycle's own rate has to fall within this fraction
%                              either side of the rate that was found.
%               .gateDetrend   order of the polynomial removed before the spectrum.
%               .minPromCoef   minimum peak prominence, as a fraction of the RHYTHM'S
%                              own standard deviation.
%               .coeffsSTD     one per feature: how many leave-one-out standard
%                              deviations from the leave-one-out median a feature may
%                              sit before its cycle is rejected.
%               .coeffsRel     (1) largest min-to-min level jump as a share of the
%                              cycle's excursion; (2) largest departure of the cycle's
%                              level from the recording's median.
%               .coeffsAbs     half-width, in cycles, of the run a rejection spreads to.
%               .excludeFirstNCycles  rejected outright.
%               .nPhaseBins    phase bins one cycle is averaged into.
%               .nControls     matched controls stored beside the cycle - the same
%                              beats folded with their timing scrambled.  0 writes
%                              none, and then no wall motion can be measured on this
%                              product.
%               .controlSeed   the controls' own random stream, so a product is
%                              reproducible and building one moves no other draw.
%               .decimationSpace  block size of the pass-1 reference cube.
%    fNames - cell array of .cxd file paths.
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.
%
% Outputs:
%    (none) - writes, through getResultsPath, beside the recording or in the project's
%             results folder:
%               <stem>_c_I_d.mat  SOURCE   source.data [Y x X x nPhaseBins] single, the
%                                          averaged cycle in the recording's own
%                                          intensity units; source.sem, the standard
%                                          error of the mean per pixel per bin;
%                                          source.control [Y x X x nPhaseBins x
%                                          nControls], the same fold with each beat's
%                                          phase scrambled - the floor every wall
%                                          measurement is divided by;
%                                          source.time [nPhaseBins x 1]
%               <stem>_c_I_r.mat  RESULTS  time, timeStamp, fps, imgI (the mean image
%                                          the cycle is drawn on), cycleTrace and
%                                          cycleTraceSEM (the whole field over one
%                                          cycle, and the error OF THAT MEAN - not the
%                                          per-pixel one, which is source.sem and is
%                                          orders of magnitude larger), the reference
%                                          trace and its clock, every candidate cycle,
%                                          which were accepted, why the others were not,
%                                          and the rate statistics
%               <stem>_c_I_s.mat  SETTINGS settings.runIntensityInternalCycle
%               <stem>_rep_cardiac-detect.jpg    the beats, the spectrum, the rejections
%               <stem>_rep_cardiac-average.jpg   the cycle, its modulation, its error
%
%    THE TIME AXIS COVERS ONE CYCLE AND t = 0 IS THE CARDIAC FOOT - the minimum of the
%    reference trace that opens the cycle, refined to a fraction of a frame.  The stored
%    phases are BIN CENTRES: bin k sits at (k-0.5)/nPhaseBins of the way from that foot
%    to the next, so results.time(1) is half a bin and not zero.  Endpoints would put the
%    same phase in the first bin and the last, and a mean over the bins would then
%    double-count the foot - which is exactly the quantity that has to equal the mean
%    image for the product to be in absolute intensity.  The bins wrap, so the cycle is
%    one whole period and can be filtered as one.
%
% Example:
%    s.libraryFolder = libraryFolder;
%    s.maskPctI      = 85;
%    s.nPhaseBins    = 25;
%    D = dir(fullfile(dataRoot,'*.cxd'));
%    runIntensityInternalCycle(s, fullfile({D.folder}',{D.name}'));
%
% Dependencies: the Bio-Formats MATLAB library (bfGetReader, bfGetPlane), getFFT,
%               and the Signal Processing Toolbox (findpeaks, butter, filtfilt).
% See also: runContrastInternalCycle, runIntensity, runIntensityBolus, getResultsPath,
%           wbFileModel, wbStepRegistry
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
%
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - PULSE DETECTION
% s.maskPctI=85; %percentile of the mean image above which a pixel is averaged to detect the pulse
% s.maxFrqIni=20; % initial max frequency of the activity of interest, Hz
% s.minFrqIni=1; % initial min frequency of the activity of interest, Hz
% s.rangeFrq=0.3; % relative frq range around the central frequency
%
% %ADJUSTED IF NECESSARY - EXCLUSION CRITERIA
% s.excludeFirstNCycles=0; %reject given number of cycles
% s.coeffsSTD=[3,2,2,2,2,3,3,2,2]; %pulses rejection coefficients relative to the feature standard deviation
% s.coeffsRel=[0.5,0.1]; %pulses rejection coefficients relative to the feature value
% s.coeffsAbs=2; %pulses rejection coefficients relative to the absolute feature value
% s.gateDetrend=3; %order of the polynomial removed from the trace before its spectrum
% s.minPromCoef=1/4; %in respect to the std of the band-passed rhythm
%
% %ADJUSTED IF NECESSARY - CYCLE CALCULATION
% s.nPhaseBins=25; %phase bins one cycle is averaged into
% s.nControls=3; %matched controls stored beside it - the same beats, timing scrambled
% s.controlSeed=20260807; %the controls' own random stream
% s.decimationSpace=4; %spatial decimation used to conserve memory while the pulse is detected

%------------- BEGIN CODE --------------
function runIntensityInternalCycle(s,fNames)

REPORTEVERYSEC = 30;    % how often a read that is running long says where it is

if ~all( cellfun(@(f) isempty(f) || contains(f,'.cxd'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".cxd".');
end

s = withDefaults(s);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  'Cardiac cycle' is
% the REGISTRY LABEL, not the function name - the house rule runIntensity follows with
% 'Angiogram'.
rep=reportOpen(s,'Cardiac cycle',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=char(fNames{fidx});
        reportFile(rep,fidx,s.fName);
        clearvars results source settings cube

        % Everything judgeable from the settings alone is judged BEFORE the recording is
        % opened, because a Bio-Formats index of a large .cxd is seventeen minutes and
        % nobody should pay it to be told about a typo.
        checkSettings(s);

        reader=bfGetReader(s.fName);
        % An ANONYMOUS handle capturing the reader by value: an onCleanup holding a
        % nested function handle would pin this workspace and never fire.  Re-assigning
        % it on the next file fires this one, which is why the explicit close below is
        % safe to follow it - closeReader swallows a double close.
        guard=onCleanup(@() closeReader(reader));

        omeMeta=reader.getMetadataStore();
        sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
        sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
        sizeT=double(omeMeta.getPlaneCount(0));
        startT=omeMeta.getImageAcquisitionDate(0);
        timeStamp=round(posixtime(datetime(string(startT),'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC')))*1000; %ms
        sampling=double(omeMeta.getPixelsTimeIncrement(0).value()); %s
        checkClock(sampling,s.fName);
        fps=1/sampling;

        %The shape of the whole product, before a single plane has been read.  The two
        %large things are sized here and treated DIFFERENTLY on purpose: the reference
        %cube is adjusted, the fold is refused.  See the header.
        nPix=sizeY*sizeX;
        checkFoldFits(nPix,s.nPhaseBins,s.nControls,sizeY,sizeX,s.fName);
        decimSpace=fitDecimation(s.decimationSpace,sizeY,sizeX,sizeT);
        if decimSpace>s.decimationSpace
            warning(['s.decimationSpace=%d needs %.1f GB for the reference cube and ' ...
                'half the free memory is %.1f GB; raised to %d (%.1f GB).'], ...
                s.decimationSpace, cubeBytes(s.decimationSpace,sizeY,sizeX,sizeT)/2^30, ...
                halfFree()/2^30, decimSpace, cubeBytes(decimSpace,sizeY,sizeX,sizeT)/2^30);
        end
        s.decimationSpaceUsed=decimSpace;

        % ---- PASS 1: the reference trace and the mean image --------------------
        [cube,imgMean]=pass1(reader,sizeY,sizeX,sizeT,decimSpace,rep,REPORTEVERYSEC);

        %Mask A - the pixels averaged into the reference trace.  It decides which pixels
        %the heartbeat is DETECTED from and goes no further; nothing downstream sees it.
        %It is built on the DECIMATED grid, from the block means of the same reduction
        %every frame of the cube went through, so a threshold here is a threshold on a
        %block average and not on one pixel - a single hot pixel cannot carry its block.
        imgD=blockMean(imgMean,decimSpace);
        maskA=imgD>prctile(imgD(:),s.maskPctI) & isfinite(imgD);
        if ~any(maskA(:))
            error('runIntensityInternalCycle:emptyMask', ...
                ['No pixel of %s is above the %g-th percentile of its own mean image, ' ...
                 'so there is nothing to detect the pulse in.'], ...
                shortName(s.fName), s.maskPctI);
        end
        %a matrix-vector product, not sum(cube.*mask,1): the second allocates a copy of
        %the largest array this function holds
        tsI=double((single(maskA(:)).'*reshape(cube,[],sizeT))./nnz(maskA)).';
        tsTime=((0:sizeT-1)./fps).';
        clear cube

        % ---- the gate ----------------------------------------------------------
        g=detectCycles(s,tsI,fps);

        % ---- PASS 2: fold the accepted cycles ----------------------------------
        C=pass2(reader,s,g,fps,sizeY,sizeX,sizeT,rep,REPORTEVERYSEC);
        reader.close();

        %BIN CENTRES, not endpoints - see the header.  This is what makes the mean over
        %the bins equal the mean image rather than double-counting the foot.
        time=(((0:s.nPhaseBins-1).'+0.5)./s.nPhaseBins).*C.cycleSeconds;

        source.data   =C.data;
        source.sem    =C.sem;
        source.control=C.control;
        source.time   =time;

        results.time=time;
        results.timeStamp=timeStamp;
        results.fps=fps;
        results.imgI=C.imgI;                    % the picture the cycle is drawn on
        %the whole field's own cycle and the error OF THAT MEAN.  source.sem is the error
        %of a PIXEL and is the wrong bar to draw round a trace averaged over a million.
        results.cycleTrace=C.trace;
        results.cycleTraceSEM=C.traceSEM;
        results.tsI=tsI;                        % the trace the beats were found in
        results.tsTime=tsTime;
        results.centralFrq=g.f0;
        results.band=g.band;
        results.cycles=g.cycles;
        results.cycleTime=g.cycleTime;
        results.cycleFrq=g.cycleFrq;
        results.accepted=g.accepted;
        results.reject=g.reject;
        results.reasons=g.reasons;              % the phrasebook lives with the rules
        results.rejectCount=g.rejectCount;
        results.acceptanceRate=g.acceptanceRate;
        results.nCycles=C.n;
        results.nOffered=size(g.acceptedCycles,1);
        results.cycleSeconds=C.cycleSeconds;
        results.rateMedian=median(g.cycleFrq(g.accepted));
        results.rateIQR=prctile(g.cycleFrq(g.accepted),[25 75]);
        settings.runIntensityInternalCycle=reportSettings(s);

        %THE IDENTITY [D11] RESTS ON, CHECKED ON EVERY RUN.  Each cycle went in with its
        %own per-pixel mean removed, so the shape sums to zero over the bins and the mean
        %image is what is left.  If this ever fails the product is a fractional change
        %wearing an absolute name, and every pulsatility index taken off it is a
        %different quantity from the one taken off '_c_BFI'.
        checkAbsolute(C.data,C.imgI,s.fName);

        drawDetect(rep,g,tsI,tsTime,fps);
        drawAverage(rep,source,results);

        % Save the settings and results
        reportWriting(rep);
        % This is the ENTRY step, so this is where a raw path becomes a product path -
        % and the only place that has to know a project may keep its results apart from
        % its recordings.  s.fName is NOT reassigned: it stays the raw recording.  With
        % no results folder set the name comes back verbatim.
        outBase=getResultsPath(s.fName,s);
        save(strrep(outBase,'.cxd','_c_I_d.mat'),'source','-v7.3','-nocompression');
        save(strrep(outBase,'.cxd','_c_I_r.mat'),'results','-v7.3','-nocompression');
        save(strrep(outBase,'.cxd','_c_I_s.mat'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
        %Before the next file allocates its own: a plain re-assignment would build the
        %new arrays while the old ones are still resident.
        clear source C
    end
end
reportClose(rep);

end

% =====================================================================
function s = withDefaults(s)
%withDefaults  The settings this step has, defaulted where a hand-written s omits one.
%   The registry preset supplies all of them, so this is for a script that names only
%   what it wants changed.  There is deliberately no dataTypeOut and no saveSource:
%   every value of this product is a mean over hundreds of beats and the quantity of
%   interest is a 0.5-2 per cent modulation on it, so an integer store would quantise
%   away exactly the signal - and the cycle IS the product, so there is nothing here
%   that could be switched off and still leave one.
d = struct('maskPctI',85,'minFrqIni',1,'maxFrqIni',20,'rangeFrq',0.3, ...
    'gateDetrend',3,'minPromCoef',1/4,'coeffsSTD',[3 2 2 2 2 3 3 2 2], ...
    'coeffsRel',[0.5 0.1],'coeffsAbs',2,'excludeFirstNCycles',0, ...
    'nPhaseBins',25,'nControls',3,'controlSeed',20260807,'decimationSpace',4);
f = fieldnames(d);
for i=1:numel(f)
    if ~isfield(s,f{i}) || isempty(s.(f{i})), s.(f{i})=d.(f{i}); end
end
end

% =====================================================================
function checkSettings(s)
%checkSettings  Everything knowable without opening the recording, checked before it is.
if ~isscalar(s.nPhaseBins) || ~isfinite(s.nPhaseBins) || mod(s.nPhaseBins,1)~=0 || s.nPhaseBins<4
    error('runIntensityInternalCycle:nPhaseBins', ...
        ['s.nPhaseBins must be a whole number of 4 or more - a cycle described by ' ...
         'fewer is not a cycle; it is %s.'], mat2str(s.nPhaseBins));
end
if ~isscalar(s.maskPctI) || ~isfinite(s.maskPctI) || s.maskPctI<0 || s.maskPctI>=100
    error('runIntensityInternalCycle:maskPctI', ...
        ['s.maskPctI is a percentile of the mean image, from 0 (every pixel) up to ' ...
         'but not including 100; it is %s.'], mat2str(s.maskPctI));
end
if ~isscalar(s.decimationSpace) || ~isfinite(s.decimationSpace) || ...
        mod(s.decimationSpace,1)~=0 || s.decimationSpace<1
    error('runIntensityInternalCycle:decimationSpace', ...
        's.decimationSpace must be a whole number of pixels, 1 or more; it is %s.', ...
        mat2str(s.decimationSpace));
end
if ~isscalar(s.minFrqIni) || ~isscalar(s.maxFrqIni) || ~(s.maxFrqIni>s.minFrqIni) || s.minFrqIni<=0
    error('runIntensityInternalCycle:band', ...
        ['The heart rate is looked for between s.minFrqIni and s.maxFrqIni Hz, so the ' ...
         'second has to be the larger and both positive; they are %g and %g.'], ...
        s.minFrqIni, s.maxFrqIni);
end
if ~isscalar(s.rangeFrq) || ~(s.rangeFrq>0)
    error('runIntensityInternalCycle:rangeFrq', ...
        's.rangeFrq is a fraction of the rate that was found and must be positive; it is %s.', ...
        mat2str(s.rangeFrq));
end
if numel(s.coeffsSTD)~=9
    error('runIntensityInternalCycle:coeffsSTD', ...
        's.coeffsSTD is one coefficient per cycle feature, so it has nine; it has %d.', ...
        numel(s.coeffsSTD));
end
if ~isscalar(s.nControls) || ~isfinite(s.nControls) || mod(s.nControls,1)~=0 || s.nControls<0
    error('runIntensityInternalCycle:nControls', ...
        ['s.nControls is how many matched controls are stored beside the cycle, so it ' ...
         'is a whole number of 0 or more; it is %s.'], mat2str(s.nControls));
end
end

% =====================================================================
function checkClock(sampling,fName)
%checkClock  A cycle product with no clock is not a cycle product.
if ~isscalar(sampling) || ~isfinite(sampling) || sampling<=0
    error('runIntensityInternalCycle:noClock', ...
        ['%s reports a time increment of %s, so its frame rate is not knowable and a ' ...
         'heart rate cannot be measured against it.'], shortName(fName), mat2str(sampling));
end
end

% =====================================================================
function b = cubeBytes(d,sizeY,sizeX,sizeT)
%cubeBytes  What pass 1's reference cube costs at a given block size, as single.
b = floor(sizeY/d)*floor(sizeX/d)*sizeT*4;
end

% =====================================================================
function a = halfFree()
%halfFree  Half the free physical memory - the bound the contrast cycle applies to its
%   own largest array, and for the same reason: it is not the only thing resident.
[~,mem] = memory;
a = mem.PhysicalMemory.Available/2;
end

% =====================================================================
function d = fitDecimation(d,sizeY,sizeX,sizeT)
%fitDecimation  Raise the block size until the reference cube fits, never lower it.
%   The cube feeds one mask-weighted spatial average per frame, and a coarser sampling
%   of an average is still that average - so adjusting it changes nothing the user asked
%   for.  The loop needs no bound: floor(size/d) falls to 1, so the cube shrinks
%   monotonically to sizeT*4 bytes, which any machine running MATLAB clears.
%
%   The one CLAMP is at the frame itself.  A block wider than the image leaves no whole
%   block at all, and the mask built on it would be empty - which this step does report
%   by name, but as "no pixel is above the percentile", which is true and describes the
%   wrong thing.  One block covering the whole frame is what such a setting means, and
%   s.decimationSpaceUsed records what was actually applied.
d   = max(1, min(d, min(sizeY,sizeX)));
lim = halfFree();
while cubeBytes(d,sizeY,sizeX,sizeT) > lim && (floor(sizeY/d)>1 || floor(sizeX/d)>1)
    d = d+1;
end
end

% =====================================================================
function checkFoldFits(nPix,nBins,nCtrl,sizeY,sizeX,fName)
%checkFoldFits  REFUSE BEFORE THE READ, and say what would have worked.
%
%   Unlike the reference cube, this one is not adjusted.  s.nPhaseBins sets how finely
%   the cycle is described, which is the scientific content of the product, so choosing
%   it is the user's job and not this function's - the same distinction runIntensity
%   draws between decimFactor (refused) and the contrast step's decimationSpace
%   (raised).  The fold holds two double accumulators and a per-pixel offset while it
%   runs, and writes two single arrays; all four are sized here.
%
%   EACH CONTROL COSTS ONE MORE DOUBLE ACCUMULATOR, ONE MORE PER-PIXEL OFFSET AND ONE
%   MORE SINGLE OUTPUT.  It needs no squared accumulator because a control has no
%   standard error to report - it is the floor, not a measurement.  The message names
%   both ways out, because dropping a control and dropping a bin are different losses:
%   fewer bins describes the cycle more coarsely, no controls means no wall motion can
%   be measured on this product at all.
perCopy = nPix*nBins*8 + nPix*8 + nPix*nBins*4;         % accumulator + offset + output
need    = perCopy*(1+nCtrl) + nPix*nBins*8 + nPix*nBins*4;   % the product's acc2 and sem
lim     = halfFree();
if need <= lim, return; end
perBin = nPix*(8+4)*(1+nCtrl) + nPix*(8+4);
nFit   = max(4,floor((lim - nPix*8*(1+nCtrl))/perBin));
error('runIntensityInternalCycle:cycleTooLarge', ...
    ['A %d-bin cycle of %s (%d x %d) with %d matched controls needs %.1f GB while it ' ...
     'is being built and half the free memory is %.1f GB.  s.nPhaseBins=%d would fit ' ...
     'at this many controls; at this frame rate anything past about 25 bins is ' ...
     'interpolation in any case.  Fewer controls also fits, but a product with none ' ...
     'cannot have its wall motion measured.'], ...
    nBins, shortName(fName), sizeY, sizeX, nCtrl, need/2^30, lim/2^30, nFit);
end

% =====================================================================
function [cube,imgMean] = pass1(reader,sizeY,sizeX,sizeT,d,rep,everySec)
%pass1  Every plane, once, into a spatially decimated cube and a full-resolution mean.
%   ONE PLANE IS RESIDENT.  The cube is what makes the mask a decision that can be taken
%   AFTER the whole recording has been seen: it holds a block mean per block per frame,
%   which is the smallest object from which any block-level mask's reference trace can
%   still be formed.
Yb=floor(sizeY/d); Xb=floor(sizeX/d);
cube=zeros(Yb,Xb,sizeT,'single');
imgMean=zeros(sizeY,sizeX);
tRead=tic; tNext=everySec;
for t=1:1:sizeT
    A=double(bfGetPlane(reader,t));
    imgMean=imgMean+A;
    cube(:,:,t)=single(blockMean(A,d));
    if toc(tRead)>=tNext
        el=toc(tRead);
        rep.emit(sprintf('  finding the beats: read %d of %d frames, %.0f%%, about %s left', ...
            t,sizeT,100*t/sizeT,hms(el*(sizeT-t)/t)));
        tNext=el+everySec;
    end
end
imgMean=imgMean./sizeT;
end

% =====================================================================
function B = blockMean(A,d)
%blockMean  Non-overlapping d-by-d block means of one plane, trailing partial dropped.
%   NOT a 'same' convolution followed by subsampling: block means have no edge
%   normalisation to go wrong and no half-pixel offset to cancel, so the mask image and
%   every frame of the cube sit on exactly the same grid by construction.  Dropping the
%   trailing partial block is the same truncation the contrast core applies to its own
%   frame blocks.
if d==1, B=A; return; end
Yb=floor(size(A,1)/d); Xb=floor(size(A,2)/d);
B=reshape(A(1:Yb*d,1:Xb*d),d,Yb,d,Xb);
B=reshape(sum(sum(B,1),3),Yb,Xb)./(d*d);
end

% =====================================================================
function g = detectCycles(s,x,fs)
%detectCycles  The heart rate, the beats between its minima, and which of them to keep.
%
%   The contrast step's structure - central frequency from the most prominent spectral
%   peak in the plausible band, minima by findpeaks with a distance and a prominence, a
%   min-to-max-to-min cycle list, nine features per cycle, outliers rejected against the
%   OTHER cycles - with the conditioning meGate measured has to be different on a
%   fluorescence trace.  See the header for what that measurement was.
x  = double(x(:));
nT = numel(x);
time = (0:nT-1).'./fs;
xd = x - polyBaseline(x,s.gateDetrend);

% ---- the rhythm's own frequency ---------------------------------------------
nShort = min(2.^nextpow2(fs/s.maxFrqIni*50), 2^floor(log2(nT)));
[pw,~,fShort] = getFFT(xd, fs, nShort, 'cpu');
pw = squeeze(pw); pw = pw(:); fShort = fShort(:);
maxF = min(s.maxFrqIni, fs/2 - eps(fs));
lim = [find(fShort>s.minFrqIni,1,'first'), find(fShort>maxF,1,'first')-1];
if isempty(lim(1)), lim(1)=[]; end
if numel(lim)<2 || isempty(lim(2)) || lim(2)<=lim(1)
    lim = [find(fShort>s.minFrqIni,1,'first'), numel(fShort)];
end
if isempty(lim(1)) || lim(2)<=lim(1)
    error('runIntensityInternalCycle:band', ...
        'No spectral line fits between %g and %g Hz at %g Hz sampling.', ...
        s.minFrqIni, s.maxFrqIni, fs);
end
[~,iPk,~,prom] = findpeaks(pw(lim(1):lim(2)));
if isempty(iPk)
    error('runIntensityInternalCycle:noPeak', ...
        ['No peak at all between %g and %g Hz in the reference trace of this ' ...
         'recording, so there is no rhythm to fold on.'], s.minFrqIni, s.maxFrqIni);
end
[~,iBest] = max(prom);
f0     = fShort(iPk(iBest)+lim(1)-1);
minFrq = max(f0*(1-s.rangeFrq), s.minFrqIni);
maxFrq = min(f0*(1+s.rangeFrq), maxF);

% ---- condition the trace ----------------------------------------------------
% filtfilt, never filter: a causal lag varying across a field is a propagating wave
% this step must not invent, and a gate is not exempt because it is a scalar.
xOsc = bandpass0(xd,[minFrq maxFrq],fs);
wLev = 2*floor(fs/minFrq)+1;                  % more than one cycle of the slowest rate
xf   = xOsc + movmean(x,wLev);

% ---- the minima, and the cycles between them --------------------------------
% THE MINIMA COME OFF THE RHYTHM ALONE.  A minimum is found by its PROMINENCE, which is
% measured against the surrounding baseline, so a beat sitting on a rising level is less
% prominent than the same beat on a falling one and the threshold would then reject
% beats for where they sit rather than for what they are.  The FEATURES are read off the
% summed trace, because a level is what several of them are about.
minProm = s.minPromCoef*std(xOsc);
[~,locsMin] = findpeaks(-xOsc,'MinPeakDistance',max(1,floor(fs/maxFrq)), ...
    'MinPeakProminence',minProm);
nC = numel(locsMin)-1;
if nC < 3
    error('runIntensityInternalCycle:noCycles', ...
        'Only %d cycles were found in %.1f s of recording.', max(nC,0), nT/fs);
end

cycles = zeros(nC,3);
for i = 1:nC
    seg = xf(locsMin(i):locsMin(i+1));
    [~,iMax] = max(seg);
    cycles(i,:) = [locsMin(i), locsMin(i)+iMax-1, locsMin(i+1)];
end

% THE FIDUCIAL IS REFINED TO SUB-FRAME PRECISION, and that is what fills the phase axis.
% Snap every beat's phase zero to the nearest frame and a frame's phase inside its beat
% is k/L for whole k and L - with a handful of frames to a beat there are only about a
% dozen values it can take, so the phase histogram comes out with holes in it.  A
% parabola through the three samples around the minimum puts the fiducial where the trace
% actually turns, which dithers the phases off the frame grid and is also simply a better
% answer: the heart does not beat on the frame clock.
tMin = (subSampleMinimum(xOsc,locsMin)-1)./fs;    % a 1-based index; time starts at 0

% ---- nine features per cycle ------------------------------------------------
feat = zeros(nC,9);
for i = 1:nC
    a=cycles(i,1); b=cycles(i,2); c=cycles(i,3);
    three     = xf(cycles(i,:));
    feat(i,1) = mean(three);                                  % level
    feat(i,2) = std(three);
    feat(i,3) = max(three)-min(three);                        % excursion
    feat(i,4) = (tMin(i+1)-tMin(i))*fs;                       % length, frames
    feat(i,5) = abs(xf(a)-xf(c));                             % min-to-min jump
    feat(i,6) = sum(diff(xf(a:b))<0);                         % dips on the upstroke
    feat(i,7) = sum(diff(xf(b:c))>0);                         % rises on the decay
    feat(i,8) = mean(diff(xf(b:c)));                          % slope of the decay
    feat(i,9) = mean(x(a+1:c));                               % the brightness it sat at
end

% ---- one row per rejection rule, one column per cycle ------------------------
nRule  = 14;
reject = false(nRule,nC);
reject(1:9,:) = featureOutliers(feat,s.coeffsSTD);
reject(10,:)  = (feat(:,5)./feat(:,3) > s.coeffsRel(1)).';
cycleFrq      = fs./feat(:,4);
reject(11,:)  = (cycleFrq < minFrq | cycleFrq > maxFrq).';
reject(12,:)  = (abs(1-feat(:,1)./median(feat(:,1))) > s.coeffsRel(2)).';
reject(13,:)  = any(~isfinite(feat),2).';
if s.excludeFirstNCycles > 0
    reject(1:13, 1:min(s.excludeFirstNCycles,nC)) = true;
end
% A single bad cycle in the middle of a good run is usually the detection slipping
% rather than the heart, and its neighbours inherit the doubt.
reject(14,:) = round(movmean(double(max(reject,[],1)),[s.coeffsAbs(1) s.coeffsAbs(1)])) > 0;

accepted = ~any(reject,1).';
if nnz(accepted) < 2
    error('runIntensityInternalCycle:allRejected', ...
        ['%d of %d cycles were rejected, so there is nothing to average.  The ' ...
         'acceptance rate is information about the preparation - loosen s.coeffsSTD ' ...
         'or s.rangeFrq only if the report page says the beats themselves are right.'], ...
        nC-nnz(accepted), nC);
end

g = struct();
g.f0              = f0;
g.band            = [minFrq maxFrq];
g.minProminence   = minProm;
g.trace           = xf;
g.rhythm          = xOsc;
g.time            = time;
g.spectrum        = struct('f',fShort,'power',pw,'iPeak',iPk(iBest)+lim(1)-1);
g.cycles          = cycles;
g.cycleTime       = [tMin(1:end-1), time(cycles(:,2)), tMin(2:end)];
g.cycleFrq        = cycleFrq;
g.features        = feat;
g.reject          = reject;
g.reasons         = ruleNames();
g.rejectCount     = sum(reject,2);
g.accepted        = accepted;
g.acceptedCycles  = cycles(accepted,:);
g.acceptedTime    = g.cycleTime(accepted,:);
g.acceptanceRate  = mean(accepted);
end

% =====================================================================
function r = ruleNames()
%ruleNames  What each row of the rejection matrix means, in a reader's words.
%   The phrasebook lives with the rules that fill it, so a report page and an export
%   column cannot describe a rule the code does not apply.
r = {'level unlike the other cycles'
     'variability unlike the other cycles'
     'excursion unlike the other cycles'
     'length unlike the other cycles'
     'start and end levels unlike the other cycles'
     'dips on the upstroke unlike the other cycles'
     'rises on the decay unlike the other cycles'
     'decay slope unlike the other cycles'
     'brightness unlike the other cycles'
     'starts and ends at very different levels'
     'its own rate is outside the accepted band'
     'sits at a brightness far from the recording''s'
     'a feature could not be computed'
     'next to other rejected cycles'};
end

% =====================================================================
function C = pass2(reader,s,g,fs,sizeY,sizeX,sizeT,rep,everySec)
%pass2  Re-read the accepted cycles and fold them onto one phase axis.
%
%   READ IN CONTIGUOUS RUNS, THROUGH THE READER THAT IS ALREADY OPEN.  Cycles are
%   min-to-min, so two neighbours share a frame whenever nothing between them was
%   rejected; one read per run keeps the pass sequential over a very large file and
%   stops the shared frame being fetched twice.
%
%   A CYCLE IS IN OR OUT, NEVER PARTLY IN.  Dropping single (cycle, bin) samples would
%   give neighbouring bins different populations of beats, and then the differences
%   BETWEEN beats appear in the cycle as structure WITHIN one - which is exactly what
%   the wall motion being looked for also looks like.
%
%   THE PER-CYCLE CENTRING IS FOLDED INTO THE INTERPOLATION WEIGHTS, and that is a real
%   saving rather than a tidiness.  Taking this cycle's own per-pixel mean out is
%   V - mean(V,2), and V is B*W, so it is B*(W - mean(W,2)) - the subtraction moves from
%   an array the size of the whole output to one the size of a handful of frames by a
%   number of bins.  The same identity gives this cycle's own mean image as B*mean(W,2),
%   a matrix-vector product.  Done the obvious way instead, every cycle walks five
%   [pixels x bins] arrays rather than three, which on the reference recording is an
%   hour of memory traffic bought for nothing.
%
%   ONLY THE CYCLE'S OWN FRAMES ENTER ITS PRODUCT, not the whole run's.  A run is read in
%   one go because that is what the disk wants; the matrix product is then taken over the
%   half-dozen columns the cycle actually spans, because that is what the arithmetic
%   wants.  Multiplying every cycle by the whole run's block is the same answer at fifty
%   times the cost.
%
%   THE CONTROLS ARE FOLDED IN THE SAME LOOP, OFF THE SAME BLOCK, AND THAT IS THE ONLY
%   WAY THEY MEAN ANYTHING.  Control j gives cycle k a random circular phase offset
%   u(k,j) drawn once before the read, so its bin b samples the beat at phase
%   mod(p(b)+u(k,j),1).  The offsets rotate a uniform set of phases into another
%   uniform set, so every bin of every control still receives exactly one interpolated
%   sample from every accepted beat, through weights of the same kind - what is
%   destroyed is which beat's phase lands in which bin, and nothing else.  Reading the
%   recording a second time to build them would cost another seventeen minutes and
%   would sample DIFFERENT frames, so the difference between the two would stop being
%   the timing.
nPix=sizeY*sizeX;
nBin=s.nPhaseBins;
nCtrl=s.nControls;
cyc =g.acceptedTime;
p   =(((0:nBin-1).')+0.5)./nBin;

acc =zeros(nPix,nBin);
acc2=zeros(nPix,nBin);
accB=zeros(nPix,1);
%One accumulator and one offset per control; no squared accumulator, because a floor
%does not carry a standard error.  A private stream, so a product is reproducible and
%drawing these moves no other random number this session takes.
accC =cell(1,nCtrl);
accBC=cell(1,nCtrl);
for j=1:nCtrl
    accC{j} =zeros(nPix,nBin);
    accBC{j}=zeros(nPix,1);
end
if nCtrl>0
    rs=RandStream('threefry','Seed',s.controlSeed);
    u =rand(rs,size(cyc,1),nCtrl);
else
    u =zeros(size(cyc,1),0);
end
%THE FIELD'S OWN CYCLE, ACCUMULATED SEPARATELY, AND IT IS NOT A DUPLICATE.  The
%per-pixel standard error below is the right error bar for a pixel and the WRONG one for
%a trace averaged over a million of them - drawn around a spatial mean it says the cycle
%is lost in its own noise when what it is describing is the noise of one pixel.  The
%error of the spatial mean cannot be derived from the per-pixel one either, because
%beat-to-beat variation is spatially correlated, so it is accumulated where it can be
%measured: two [1 x nBin] rows, which cost nothing beside the arrays above.
accM =zeros(1,nBin);
accM2=zeros(1,nBin);
n   =0;

runs=cycleRuns(g.acceptedCycles,maxRunFrames(sizeY,sizeX,sizeT));
tRead=tic; tNext=everySec;
for r=1:1:size(runs,1)
    kk=runs(r,1):runs(r,2);
    % a frame of margin either side: the refined fiducial can sit up to half a frame
    % before the first sample of the cycle and after the last
    i0=max(1,     floor(min(cyc(kk,1))*fs)+1-1);
    i1=min(sizeT, ceil( max(cyc(kk,3))*fs)+1+1);
    if i1-i0 < 3, continue; end
    B=zeros(sizeY,sizeX,i1-i0+1,'single');
    for t=i0:1:i1
        B(:,:,t-i0+1)=single(bfGetPlane(reader,t));
    end
    B=reshape(B,nPix,[]);
    nF=size(B,2);
    tB=((i0-1:i1-1).')./fs;
    for k=kk
        tgt=cyc(k,1)+p.*(cyc(k,3)-cyc(k,1));
        j0=max(1, floor(cyc(k,1)*fs)+1-1-i0+1);
        j1=min(nF,ceil( cyc(k,3)*fs)+1+1-i0+1);
        if j1-j0<2, continue; end
        tBk=tB(j0:j1);
        %ONE ACCEPTANCE TEST FOR THE CYCLE AND ALL ITS CONTROLS.  A rotated phase set
        %reaches marginally closer to both ends of the beat than the product's centred
        %one does, so the span checked here is the whole beat rather than the product's
        %own targets.  Skipping a cycle for the product but keeping it for a control -
        %or the reverse - would let the two average different populations of beats, and
        %then the difference BETWEEN beats appears in the ratio as detection.
        if cyc(k,1)<tBk(1) || cyc(k,3)>tBk(end), continue; end
        Bk=double(B(:,j0:j1));
        W =interpWeights(tBk,tgt);               % [nFk x nBin], double
        wb=mean(W,2);
        V =Bk*(W-wb);                            % the SHAPE, centring folded in
        b =Bk*wb;                                % this cycle's own mean image
        acc =acc +V;
        acc2=acc2+V.*V;
        accB=accB+b;
        mk  =mean(V,1)+mean(b);                  % this cycle's whole-field trace
        accM =accM +mk;
        accM2=accM2+mk.*mk;
        for j=1:nCtrl
            tgtC=cyc(k,1)+mod(p+u(k,j),1).*(cyc(k,3)-cyc(k,1));
            Wc  =interpWeights(tBk,tgtC);
            wbc =mean(Wc,2);
            accC{j} =accC{j} +Bk*(Wc-wbc);
            accBC{j}=accBC{j}+Bk*wbc;
        end
        n=n+1;
    end
    if toc(tRead)>=tNext
        el=toc(tRead);
        rep.emit(sprintf('  averaging the beats: %d of %d, %.0f%%, about %s left', ...
            r,size(runs,1),100*r/size(runs,1),hms(el*(size(runs,1)-r)/r)));
        tNext=el+everySec;
    end
end

if n < 2
    error('runIntensityInternalCycle:tooFewFolded', ...
        ['Only %d accepted cycle could be read back, so there is nothing to average; ' ...
         'a mean of one cycle is that cycle.'], n);
end

shape=acc./n;
sd   =sqrt(max(acc2./n - shape.*shape,0));
base =accB./n;

%THE MEAN IMAGE IS ADDED BACK.  See the header, and see checkAbsolute.
C=struct();
C.data =single(reshape(base+shape,sizeY,sizeX,nBin));
C.sem  =single(reshape(sd./sqrt(n), sizeY,sizeX,nBin));
C.imgI =reshape(base,sizeY,sizeX);
C.n    =n;
%Each control gets its own mean image added back, exactly as the product does, so a
%control is a picture that does not beat rather than a difference image - which is what
%lets the same wall estimator read both without knowing which it has.
C.control=zeros(sizeY,sizeX,nBin,nCtrl,'single');
for j=1:nCtrl
    C.control(:,:,:,j)=single(reshape(accBC{j}./n + accC{j}./n, sizeY,sizeX,nBin));
end
%the whole field's own cycle and the error OF THAT MEAN, from the per-cycle traces
C.trace   =(accM./n).';
C.traceSEM=(sqrt(max(accM2./n - (accM./n).^2,0))./sqrt(n)).';
used   =g.acceptedTime;
C.cycleSeconds=median(used(:,3)-used(:,1));
end

% =====================================================================
function m = maxRunFrames(sizeY,sizeX,sizeT)
%maxRunFrames  How many frames one pass-2 read may hold, from the free memory.
%   An eighth rather than a half: the two fold accumulators are already resident and are
%   the larger pair.  Never a setting - it is a property of the machine and changing it
%   cannot change an answer.
m = max(4, min(sizeT, floor((halfFree()/4)/(sizeY*sizeX*4))));
end

% =====================================================================
function runs = cycleRuns(acceptedCycles,maxFrames)
%cycleRuns  Accepted cycles grouped into contiguous, length-capped runs.
%   Two neighbours in the list share a frame whenever nothing between them was rejected,
%   which is what makes reading the stretch in one go worth doing.  Each row is
%   [firstCycle,lastCycle].
runs=zeros(0,2);
iFirst=1;
for i=1:1:size(acceptedCycles,1)
    if i==size(acceptedCycles,1) ...
            || acceptedCycles(i+1,1)~=acceptedCycles(i,3) ...
            || (acceptedCycles(i+1,3)-acceptedCycles(iFirst,1)+1)>maxFrames
        runs(end+1,:)=[iFirst,i]; %#ok<AGROW>
        iFirst=i+1;
    end
end
end

% =====================================================================
function W = interpWeights(tB,tgt)
%interpWeights  Linear interpolation onto the phase axis, as a matrix.
%   Two non-zeros per column, so ONE matrix product resamples every pixel of the block
%   at once.  Doing it with interp1 per pixel is the same arithmetic through a million
%   function calls.  DOUBLE, not single: the centring is folded into these weights, so
%   how nearly their centred columns sum to zero is exactly how nearly the finished
%   cycle averages to its own mean image - which is the identity the product's units
%   rest on, and it is not worth risking to save a few kilobytes.
nF=numel(tB); nB=numel(tgt);
dt=tB(2)-tB(1);
u =(tgt-tB(1))./dt;                        % fractional index, zero based
lo=min(max(floor(u),0),nF-2);
w =u-lo;
W =zeros(nF,nB);
W(sub2ind([nF nB],lo+1,(1:nB).'))=1-w;
W(sub2ind([nF nB],lo+2,(1:nB).'))=w;
end

% =====================================================================
function checkAbsolute(dataCycle,imgI,fName)
%checkAbsolute  The round trip [D11] rests on: the cycle's mean IS the mean image.
%   Asserted rather than assumed, because the failure is silent - a product that is a
%   fractional change under an absolute name reads perfectly well and every pulsatility
%   index taken off it is a different quantity from the same index on '_c_BFI'.  The
%   tolerance is a floating-point one, scaled by the image: the accumulation is in
%   double and the store is single, so the only difference allowed is the rounding of
%   that store.
m   = mean(double(dataCycle),3);
scl = max(abs(imgI(:)));
if scl==0, scl=1; end
rel = max(abs(m(:)-imgI(:)))/scl;
if ~(rel < 1e-5)
    error('runIntensityInternalCycle:notAbsolute', ...
        ['The averaged cycle of %s does not average to its own mean image (%.3g of ' ...
         'the image scale), so the product is not in absolute intensity.'], ...
        shortName(fName), rel);
end
end

% =====================================================================
function drawDetect(rep,g,tsI,tsTime,fs)
%drawDetect  What the beat detection saw, what it kept, and why it dropped the rest.
fh=reportFigure(rep,'cardiac-detect');
tl=tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

nShow=min(numel(tsI),max(8,round(4*fs/g.f0)));
ax=nexttile(tl,1);
plot(ax,tsTime(1:nShow),tsI(1:nShow))
hold(ax,'on')
plot(ax,tsTime(1:nShow),g.trace(1:nShow))
hold(ax,'off')
legend(ax,{'Measured','Conditioned'},'Location','best')
xlabel(ax,'Time, s'); ylabel(ax,'Intensity')
title(ax,'The first few beats')

ax=nexttile(tl,3);
plot(ax,g.spectrum.f,g.spectrum.power)
hold(ax,'on')
plot(ax,g.spectrum.f(g.spectrum.iPeak),g.spectrum.power(g.spectrum.iPeak),'or')
hold(ax,'off')
xlabel(ax,'Frq, Hz'); ylabel(ax,'Amplitude')
title(ax,sprintf('Heart rate %.2f Hz, kept %.0f%% of %d beats', ...
    g.f0,100*g.acceptanceRate,size(g.cycles,1)))

ax=nexttile(tl,2);
plot(ax,tsTime,g.trace,'Color',[0.7 0.7 0.7])
hold(ax,'on')
for i=1:1:size(g.cycles,1)
    idx=g.cycles(i,1):g.cycles(i,3);
    if g.accepted(i), col='g'; else, col='r'; end
    plot(ax,tsTime(idx),g.trace(idx),col)
end
hold(ax,'off')
xlabel(ax,'Time, s'); ylabel(ax,'Intensity')
title(ax,sprintf('%d beats kept, %d dropped', ...
    nnz(g.accepted),nnz(~g.accepted)))

%THE PHRASEBOOK, DRAWN.  A rejection count nobody can read is a rejection count nobody
%checks - and the acceptance rate is meant to be REPORTED, never tuned up: an awake
%animal fails more cycles than an anaesthetised one, and that is information about the
%preparation rather than a number to improve.
ax=nexttile(tl,4);
barh(ax,g.rejectCount)
set(ax,'YTick',1:numel(g.reasons),'YTickLabel',g.reasons,'YDir','reverse')
ax.FontSize=7;
xlabel(ax,'Beats caught by this rule')
title(ax,'Why beats were dropped')

reportSave(rep,fh,'cardiac-detect');        % reportSave titles the page
end

% =====================================================================
function drawAverage(rep,source,results)
%drawAverage  The cycle, how much of it moves, and how sure of it the recording is.
fh=reportFigure(rep,'cardiac-average');
tl=tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl,1);
imagesc(ax,results.imgI)
clim(ax,[prctile(results.imgI(:),1),prctile(results.imgI(:),99)])
colorbar(ax); axis(ax,'image')
title(ax,'Mean intensity over the cycle')

%THE PANEL THE PRODUCT EXISTS FOR.  (max-min)/mean over the cycle is the per-pixel
%pulsatility index, and it is the SAME quantity as on a '_c_BFI' product precisely
%because the mean image was added back.
%AND THE SAME MAP OF THE MATCHED CONTROLS IS WHAT SAYS HOW MUCH OF IT IS REAL.  The
%vessel tree appearing in this panel is already a check nothing else provides - photon
%noise and a mis-detected rhythm are both spatially unstructured and neither can draw
%an angiogram - but it is a shape, not a number.  The controls are the same beats with
%their timing scrambled, so the modulation they show is what this recording produces
%with no rhythm to find, and the two medians beside each other turn the shape into a
%ratio a reader can quote.
ax=nexttile(tl,2);
d=double(source.data);
pulse=(max(d,[],3)-min(d,[],3))./max(results.imgI,eps);
imagesc(ax,100*pulse)
clim(ax,[0,prctile(100*pulse(:),99)])
colorbar(ax); axis(ax,'image')
nCtrl=size(source.control,4);
if nCtrl>0
    cm=zeros(1,nCtrl);
    for j=1:nCtrl
        dj=double(source.control(:,:,:,j));
        cm(j)=median((max(dj,[],3)-min(dj,[],3))./max(results.imgI,eps),'all');
    end
    title(ax,sprintf(['Modulation over the cycle, %% of the mean\n' ...
        'median %.2f%% against %.2f%% with the timing scrambled (x%.1f)'], ...
        100*median(pulse,'all'), 100*median(cm), median(pulse,'all')/max(median(cm),eps)))
else
    title(ax,sprintf(['Modulation over the cycle, %% of the mean\n' ...
        'median %.2f%% - no matched control was stored'], 100*median(pulse,'all')))
end

%THE ERROR BAND IS THE ERROR OF THIS TRACE, not the typical error of a pixel.  Those
%differ by orders of magnitude - a mean over a million pixels is far better determined
%than any one of them - and drawing the per-pixel one round a whole-field mean would say
%the cycle is lost in noise when what is being described is one pixel's noise.
ax=nexttile(tl,3,[1,2]);
m=results.cycleTrace;
e=results.cycleTraceSEM;
t=[source.time; source.time(end)+(source.time(2)-source.time(1))];
mm=[m; m(1)]; ee=[e; e(1)];              % the cycle wraps - draw it closed
fill(ax,[t;flipud(t)],[mm-ee;flipud(mm+ee)],[0.85 0.85 0.95],'EdgeColor','none')
hold(ax,'on')
plot(ax,t,mm,'k')
hold(ax,'off')
xlabel(ax,'Time from the cardiac foot, s')
ylabel(ax,'Intensity')
title(ax,sprintf(['The whole field over one cycle, from %d beats, shaded by the ' ...
    'standard error OF THIS MEAN'],results.nCycles))

reportSave(rep,fh,'cardiac-average');
end

% =====================================================================
function reject = featureOutliers(feat,coeffsSTD)
%featureOutliers  Rules 1-9: a cycle whose feature stands out from the rest.
%   reject(k,i) is true when |f - median(f without i)| > coeffsSTD(k)*std(f without i).
%   LEAVE-ONE-OUT, so a single wild cycle cannot widen the spread it is then compared
%   against and hide inside it.
%
%   Taken in CLOSED FORM rather than recomputed per (cycle,feature) pair, which builds
%   an N-element index vector 9N times and grows as N^2 - the contrast step measured
%   5.5 s of it at 7 000 cycles, and a .cxd of the reference protocol has 3 000.
%   With mu=mean(x) and Q=sum((x-mu).^2) over all N values of one feature,
%       mu_i = (N*mu-x_i)/(N-1),  Q_i = Q-(x_i-mu).^2*N/(N-1),
%       std_i = sqrt(max(Q_i,0)/(N-2))     std of N-1 samples normalises by N-2.
%   Q is the CENTRED sum: sum(x.^2)-N*mu^2 would cancel two digits on feature 9, whose
%   values are intensities with a small spread.
nF=size(feat,2);
nC=size(feat,1);
reject=false(nF,nC);
for k=1:1:nF
    x=feat(:,k);
    dev=x-mean(x);
    devQ=dev.*dev;
    devSum=sum(devQ);
    if nC<=2 || ~isfinite(devSum)
        %Two exceptions, both rare, neither with a closed form.  With two cycles the
        %leave-one-out set is one sample, whose std MATLAB defines as 0 while the
        %identity reads 0/0; and a non-finite Q covers the rest - median and std return
        %NaN for every cycle EXCEPT the one holding the offending value, whose own
        %leave-one-out set has just dropped it, and no sum reproduces that asymmetry.
        for i=1:1:nC
            xLoo=x([1:i-1,i+1:end]);
            reject(k,i)=abs(x(i)-median(xLoo))>coeffsSTD(k)*std(xLoo);
        end
        continue
    end

    sdLoo=sqrt(max(devSum-devQ.*(nC/(nC-1)),0)./(nC-2));

    [xs,ord]=sort(x);
    rnk=zeros(nC,1);
    rnk(ord)=1:nC;
    if rem(nC,2)==0                             % the leave-one-out set is odd
        medLoo=xs(nC/2+(rnk<=nC/2));
    else                                        % even - two order statistics to combine
        jLo=(nC-1)/2;
        a=xs(jLo+(rnk<=jLo));
        b=xs(jLo+1+(rnk<=jLo+1));
        medLoo=a+(b-a)./2;
        acrossZero=sign(a)~=sign(b);            % median's own guard against overflow
        medLoo(acrossZero)=(a(acrossZero)+b(acrossZero))./2;
    end

    reject(k,:)=(abs(x-medLoo)>coeffsSTD(k).*sdLoo).';
end
end

% =====================================================================
function p = subSampleMinimum(x,loc)
%subSampleMinimum  Where the trace really turns, in fractional samples.
%   A parabola through the sample at the minimum and its two neighbours.  Clamped to
%   half a sample either way: a vertex further out than that means the three points do
%   not describe a turn, and then the sample itself is the best available answer.
loc=loc(:);
p  =double(loc);
in =loc>1 & loc<numel(x);
a=x(loc(in)-1); b=x(loc(in)); c=x(loc(in)+1);
den=a-2*b+c;
d=zeros(size(den));
ok=den~=0;
d(ok)=0.5*(a(ok)-c(ok))./den(ok);
p(in)=p(in)+min(max(d,-0.5),0.5);
end

% =====================================================================
function y = bandpass0(x,band,fs)
%bandpass0  Zero-phase band-pass.  filtfilt, never filter.
w=band./(fs/2);
w=min(max(w,1e-6),0.999);
[bb,aa]=butter(2,w,'bandpass');
y=filtfilt(bb,aa,double(x(:)));
end

% =====================================================================
function b = polyBaseline(x,order)
%polyBaseline  The slow decay of the dye, as a low-order polynomial in time.
%   It comes off BEFORE the spectrum: a monotone decay puts power straight into the low
%   end of the search band.  The raw trace is kept for the level and for feature 9,
%   which is about the brightness a cycle actually sat at.
n=numel(x);
t=((1:n).'-1)./max(n-1,1);
b=polyval(polyfit(t,x,order),t);
end

% =====================================================================
function s = hms(sec)
%hms  Seconds as minutes and seconds, for a line an operator reads while waiting.
s = sprintf('%d min %02d s', floor(sec/60), round(mod(sec,60)));
end

% =====================================================================
function n = shortName(f)
%shortName  Name and extension - what the user sees on disk.
[~,b,e] = fileparts(char(f));
n = [b e];
end

% =====================================================================
function closeReader(reader)
try
    reader.close();
catch
end
end
%------------- END OF CODE --------------
