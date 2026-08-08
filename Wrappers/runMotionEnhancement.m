%runMotionEnhancement - How far every vessel wall moved, and whether it moved at all.
%
%   EIGHT PER-SEGMENT COLUMNS AND THE EVIDENCE BEHIND THEM.  Perpendicular cuts are
%   placed along every segment the library's own segmentation found, both walls are
%   located on every cut to a fraction of a pixel, the cuts are collapsed along the
%   segment, and the result is divided by the same measurement on a copy of the same
%   recording with its timing destroyed.  A segment that does not clear that floor
%   reports NaN - not a small number.
%
%   EVERY MOTION NUMBER SHIPS WITH ITS OWN FLOOR OR IT DOES NOT SHIP.  That is the
%   single sentence Core/MotionEnhancement/ exists to enforce and it is why this step
%   refuses a cardiac product that carries no matched control.  On the reference
%   recording the cycle-averaged product clears its control by 5.3 to 19.1 times while
%   the CONTINUOUS mode of the same block is 1.0 - exactly its own control - and the
%   continuous movie is the most convincing-looking of the four.  A wall amplitude
%   without a floor beside it is not a weak measurement, it is not a measurement.
%
%   WHAT IT DOES DEPENDS ON WHICH PRODUCT IT IS ON, AND BOTH ARMS REPORT.
%
%     '_c_I' the averaged cardiac cycle - the amplitude is the fundamental of the
%            closed cycle, and the floor is the matched controls stored beside it by
%            runIntensityInternalCycle: the same beats folded with their phases
%            scrambled.  This is the arm every number in the findings was measured on.
%
%     '_a_I' / '_b_I' a frame series - the SAME estimator on the SAME cuts, in a
%            band, and the floor is the identical measurement in a band of the same
%            width where no physiology sits.  On the reference recording every row
%            comes back NaN because the confidence reads about 1.0.  THAT IS THE
%            CORRECT OUTPUT AND THE STEP MUST PRODUCE IT: a step that declines the
%            continuous product teaches nobody anything, and one that returns a table
%            of NaN beside the sentence "no wall motion is detectable without cardiac
%            averaging" teaches them what five sessions of that package established.
%
%   THE STEP IS SELF-CONTAINED IN THE SAME SENSE runDynamicSegmentation IS.  The
%   labelled centre lines are not persisted by anything, so they are re-derived here
%   from the two things that are - results.cMask and
%   settings.runSegmentation.edgeSize - through the same getSegmentationLabels call.
%   A product that was never segmented is refused by name, before any array is read.
%
%   THE THREE GATES, AND THE FIRST IS THE ONE THAT DOES THE WORK.  A row carries a
%   number only if the segment is wide enough for this estimator to resolve, was
%   reduced over enough cuts, and cleared its own floor.  The expectation that the
%   confidence would be the honest one is REVERSED by measurement: on the reference
%   field the confidence gate alone passes 609 of 1 422 segments and 585 of those are
%   vessels this estimator cannot resolve.  A scrambled control destroys the TIMING,
%   so it catches anything that is not rhythm-locked - and a rhythm-locked change in
%   the tissue behind a wall is rhythm-locked.  Only the calibre limit separates
%   those two, and no per-segment statistic measured on that data reproduces it.
%
%   s.minDiameterUm IS A STATED LIMIT OF ONE RIG, NOT A DERIVED ONE, and it decides
%   98 per cent of the table.  37 um is where Core/MotionEnhancement/README.md puts
%   the smallest calibre with a signal above its own control (36.9 um) and the
%   smallest the raw wall measurement reaches (38 um), measured once, on one
%   recording, on one animal.  It does NOT improve with sampling - the same segments
%   measured at 2.5 and 5.0 um/px give a wall-amplitude ratio of median 1.00 - so it
%   does not follow the pixel size.  The report page says all of this on the page,
%   because a settings tooltip is not where the reader of a figure looks.
%
%   THE CALIBRE THE SIZE GATE READS IS THE SEGMENTATION'S, NOT THE WALL FIT'S.  It is
%   results.sMetrics.diameter, twice the mean distance transform along the centre
%   line.  A gate that a noisy fit can widen for itself is not a gate.
%
%   IT ASKS NOTHING ABOUT POLARITY AND IT DOES NOT CALL getVesselPolarity.  meKymograph
%   finds the brightest point of a cut and walks outwards to the half-way level, so it
%   assumes a bright vessel - which on this branch is the stated assumption rather than
%   a question put to the operator (author, 2026-08-08).  A dark-vessel recording
%   produces wrong numbers here that look entirely ordinary, exactly as it does at
%   segmentation and at background removal, and the trade is taken knowingly.
%
%   IT WORKS IN MICROMETRES AND HAS NO DEFAULT FOR THEM.  The headline gate is a
%   physical calibre and the cut's margin is a physical distance, so there is nothing
%   here that has a sensible answer in pixels; 2.5 um/px is measured for one rig and a
%   default would be a silent wrong answer at any other magnification.  s.pixelSize is
%   refused by name before the first product is opened, exactly as runBackgroundRemoval
%   refuses it.  results.pixelSize is deliberately NOT read: one calibration, in one
%   place, is the point.
%
%   THE AVI IS OFF BY DEFAULT AND ITS CONTROL SHIPS WITH IT OR NEITHER DOES.  The
%   deliverable is the columns; the movie is an illustration, and a magnified movie is
%   persuasive whether or not it is true.  When s.writeVideo is on, the magnified
%   cycle and the magnified control are written as a pair into a folder the user
%   names, defaulting beside the results tree and never inside it - 0.1 to 0.8 GB
%   each, against the 0.5 MB of a report image, and the workbench and the explorer
%   sweep the results tree by name.
%
% Syntax:
%    runMotionEnhancement(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%             .pixelSize      micrometres per pixel.  NO DEFAULT - refused by name.
%             .cutSpacing     arc length between cuts along a segment, px.
%             .endMarginRadii no cut within this many radii of a segment's end.
%             .spanRadii .spanPadUm   a cut's half-length: this many radii plus this
%                             many micrometres of tissue.  The pad is a PHYSICAL
%                             margin and is converted to pixels once, here.
%             .smoothRadii .smoothMin .interpF .maxCLR   the centre line.
%             .minCuts        fewest cuts a segment may be reduced over.
%             .cutStepPx .widthTol .edgeTol   how a cut is measured and refused.
%             .minDiameterUm .confMin .minCohere   the three gates.
%             .bandHz         the rhythm's band on a CONTINUOUS product, Hz.  NO
%                             DEFAULT - refused by name when such a product is met.
%             .writeVideo .videoFolder .alpha .levels   the optional AVI pair.
%             .sMinL .prchNSize .correctNodes .simR .difR .parforSegmentationLabels
%                             the labelling settings, which must match the ones
%                             runSegmentation ran with.
%    fNames - cell array of '*_I_r.mat' paths.
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.
%
% Outputs:
%    (none) - writes IN PLACE into the product it was given:
%       results.sMetrics gains EIGHT columns, NaN in all of them for a segment that
%       failed any gate and for every row that is not a vessel:
%           wmWallPx    wall displacement amplitude, peak to peak, recording px
%           wmWallUm    the same, micrometres
%           wmDilRel    the diameter's change over the cycle, per cent of itself
%           wmCentrePx  how far the vessel's CENTRE travelled, peak to peak, px
%           wmCentreUm  the same, micrometres
%           wmDilShare  wall / (wall + centre): 1 is all calibre change, 0 is all
%                       sliding.  Measured median across the reference field: 0.25,
%                       and the centre exceeds the wall in 20 of 24 segments, so a
%                       "pulsatile diameter" taken without this split is reporting
%                       the preparation moving
%           wmConf      the measurement divided by the same measurement with the
%                       timing scrambled
%           wmCuts      how many cuts the segment was reduced over
%       results.wallMotion - the same rows, ungated, plus the floor, both
%       coherences, the geometry and WHICH gate refused each row.
%       <stem>_rep_wall-motion.jpg
%       Optionally <folder>/<stem>_wallMotion.avi and its _control.avi.
%
%    EVERY DISPLACEMENT IS PEAK TO PEAK.  Core/MotionEnhancement/README.md quotes
%    zero-to-peak throughout - its 0.003 to 0.083 px and its 0.0042 px floor - so
%    these columns read TWICE those, and they are per SEGMENT rather than per cut:
%    the package's numbers came from cuts a probe SEARCHED for, and averaging a whole
%    vessel removes the hot spot it found.  The two sets are not comparable and
%    nobody should quote one against the other.
%
% Example:
%    s.libraryFolder = libraryFolder;
%    s.pixelSize     = 2.5;
%    s.minDiameterUm = 37;
%    runMotionEnhancement(s, getFileNamesList(rootFolder,'*_c_I_r.mat'));
%
% Dependencies: getSegmentCuts, getWallMotion, getSegmentationLabels, meKymograph,
%               and (only when s.writeVideo is on) meMask, meShift, meWriteVideo,
%               meSettings.  Image Processing Toolbox; Signal Processing Toolbox on a
%               continuous product.
% See also: getWallMotion, getSegmentCuts, runIntensityInternalCycle,
%           runSegmentation, runDynamicSegmentation, runTopologyAnalysis
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
%
% %ADJUSTED (OR VERIFIED) PER RIG - THERE IS NO DEFAULT
% s.pixelSize=2.5; %micrometres per pixel
%
% %ADJUSTED IF NECESSARY - WHERE THE CUTS GO
% s.cutSpacing=3; %arc length between cuts along a segment, px
% s.endMarginRadii=1.5; %no cut within this many vessel radii of a segment's end
% s.spanRadii=3; %a cut reaches this many radii either side of the centre line...
% s.spanPadUm=30; %...plus this much tissue beyond it, micrometres
% s.smoothRadii=1; %tangent smoother half-window, in radii
% s.smoothMin=3; %...floored here, px
% s.interpF=4; %points per pixel along the resampled centre line
% s.maxCLR=Inf; %refuse a segment more curved than this.  Inf lets every shape through
% s.minCuts=3; %fewest cuts a segment may be reduced over
%
% %ADJUSTED IF NECESSARY - HOW A CUT IS MEASURED
% s.cutStepPx=0.25; %sample spacing along a cut, px
% s.widthTol=[0.5,2]; %fitted diameter over the mask's, accepted range
% s.edgeTol=1.0; %10-to-90 wall edge as a share of the diameter, above which a cut is
%                %refused.  NOT meProbe's 0.5 - see the header
%
% %ADJUSTED IF NECESSARY - WHAT IS ALLOWED TO CARRY A NUMBER
% s.minDiameterUm=37; %vessels narrower than this report nothing
% s.confMin=2; %a segment must double its own scrambled-timing floor
% s.minCohere=0; %refuse a segment whose cuts disagree more than this
%
% %CONTINUOUS PRODUCTS ONLY - THERE IS NO DEFAULT
% s.bandHz=[]; %the heartbeat's band, Hz.  The Cardiac cycle step reports the rate
%
% %OPTIONAL - THE ILLUSTRATION, WHICH IS NEVER THE DELIVERABLE
% s.writeVideo=false;
% s.videoFolder=''; %empty writes beside the results folder, never inside it
% s.alpha=[53.1,106.2,212.5]; %amplification per level - one alpha serves ONE calibre
% s.levels=3:5; %pyramid levels OF THE RECORDING

%------------- BEGIN CODE --------------
function runMotionEnhancement(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_I_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_I_r.mat".  ' ...
        'Every step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_I_r.mat'').']);
end

s = withDefaults(s);

% MICROMETRES, WITH NO DEFAULT, REFUSED BEFORE THE FIRST PRODUCT IS OPENED.  The
% headline gate is a calibre and the cut's pad is a margin of tissue, so there is no
% answer to this step in pixels at all - a 37 um gate applied as 37 px would keep and
% drop the wrong vessels with nothing about the result looking wrong.
if ~isfield(s,'pixelSize') || isempty(s.pixelSize) || ~isscalar(s.pixelSize) || ...
        ~isfinite(s.pixelSize) || s.pixelSize<=0
    error('runMotionEnhancement:pixelSizeNotSet', ...
        ['This step measures vessels in micrometres, so it needs to know how many ' ...
         'micrometres a pixel is on this recording.  Set it and run again.']);
end

% The physical margin, converted ONCE.  Everything after this line is in pixels of the
% product's own grid and no function below has to remember which unit it is in.
s.spanPadPx = s.spanPadUm/s.pixelSize;

% What each file recorded about its own segmentation, judged from the sidecars alone
% and BEFORE any array is read: a product that was never segmented is refused by name
% here rather than on the line that reads results.cMask.  The sidecars are kilobytes.
prov = readSegmentationProvenance(fNames);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  'Wall motion' is
% the REGISTRY LABEL, not the function name.
rep = reportOpen(s,'Wall motion',fNames);

for fidx = 1:1:numel(fNames)
    if reportCancelled(rep), break; end          % cooperative cancel between files
    if ~isempty(fNames{fidx})

    s.fName = fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results settings source

    load(getProductPath(s.fName,'s'),'settings');
    load(s.fName,'results');
    dPath = getProductPath(s.fName,'d');
    if ~isfile(dPath)
        error('runMotionEnhancement:noFrames', ...
            ['%s has no frames stored beside it, and a wall is measured in the ' ...
             'pictures, not in a table.  Re-run the step that made this recording ' ...
             'with its frames kept.'], shortName(s.fName));
    end
    load(dPath,'source');

    % ---- which arm, and what the floor is going to be -------------------------
    mode = modeOf(s.fName);
    sM   = s;
    if strcmp(mode,'cycle')
        sM.motionMode = 'cycle';
        if ~isfield(source,'control') || isempty(source.control)
            error('runMotionEnhancement:noControl', ...
                ['%s was made before its comparison beats could be stored, so there ' ...
                 'is nothing to say how much of its movement is real.  Re-run ' ...
                 'Cardiac cycle on the recording and this step will have a floor to ' ...
                 'measure against.'], shortName(s.fName));
        end
        nCtrl = size(source.control,4);
    else
        sM.motionMode = 'band';
        if isempty(s.bandHz) || numel(s.bandHz)~=2 || ~all(isfinite(s.bandHz)) || ...
                s.bandHz(2)<=s.bandHz(1)
            error('runMotionEnhancement:bandNotSet', ...
                ['%s is a continuous recording, so this step has to be told which ' ...
                 'rhythm to look for - the Cardiac cycle step reports the heart rate ' ...
                 'it found on the same animal.  Give it as a band in beats per ' ...
                 'second and run again.'], shortName(s.fName));
        end
        sM.fps       = 1/median(diff(double(results.time)));
        sM.offBandHz = offBandOf(s.bandHz, sM.fps);
        nCtrl        = 1;                    % the off band IS the control, and it is
    end                                      % measured on the same trace

    % ---- the labelling transients, re-derived from the persisted seam ---------
    % The same two arrays runDynamicSegmentation and runTopologyAnalysis re-derive,
    % from the same call.  sLines gives the centre lines the cuts are placed along and
    % dMask the wall+lumen mask their local radius is measured in.
    edgeSize = settings.runSegmentation.edgeSize;
    [~,~,sLines,~,~,dMask] = getSegmentationLabels(results.cMask,edgeSize,s);
    R = bwdist(~dMask);

    % ---- the cuts, once, and they are used on the product and on every control --
    % A control is the same recording with the timing destroyed; measuring it through
    % DIFFERENT cuts would compare two geometries rather than two timings.
    T    = results.sMetrics;
    isVes= T.category==5 & isfinite(T.idx);
    ids  = T.idx;
    n    = height(T);
    cuts = cell(n,1);
    q    = repmat(struct('pathPx',NaN,'clr',NaN,'nSkel',0,'nPoint',0, ...
                         'spurFrac',NaN,'why',''),n,1);
    for i = find(isVes(:)).'
        [cuts{i},q(i)] = getSegmentCuts(sLines,ids(i),R,s);
    end
    nOffered = cellfun(@numel,cuts);

    % ---- refuse a sweep that will not finish, before it starts ----------------
    checkSweepFits(sum(nOffered), size(source.data,3), nCtrl, mode, s);

    % ---- the sweep -------------------------------------------------------------
    % THE CUBE AND EVERY CONTROL ARE CAST TO DOUBLE ONCE, HERE, AND THAT IS NOT A
    % TIDINESS.  meKymograph casts whatever it is handed, so a cast inside the segment
    % loop is one copy of the whole field per segment per control - on a field of
    % fourteen hundred segments with three controls that is four thousand copies of a
    % hundred-megabyte array, which is the entire cost of the sweep and nothing else.
    % The single-precision originals are dropped as each one is converted, so the peak
    % is the doubles alone.
    X  = double(source.data);            source.data    = [];
    Xc = cell(1,nCtrl);
    if strcmp(mode,'cycle')
        for j = 1:nCtrl, Xc{j} = double(source.control(:,:,:,j)); end
        source.control = [];
    end
    tree = blankTree(n,ids);
    reportProgressEvery = 30;  tSweep = tic;  tNext = reportProgressEvery;
    todo = find(nOffered>0).';
    for c = 1:numel(todo)
        i = todo(c);
        m = getWallMotion(X,cuts{i},sM);
        if m.nCut==0, continue; end
        tree.wallPxRaw(i)    = m.wallAmp;
        tree.centrePxRaw(i)  = m.centreAmp;
        tree.cohere(i)       = m.wallCohere;
        tree.centreCohere(i) = m.centreCohere;
        tree.nCut(i)         = m.nCut;
        tree.diameterPx(i)   = 2*m.radiusPx;
        tree.taperPx(i)      = 2*m.radiusSdPx;
        tree.edgeUm(i)       = m.edgePx*s.pixelSize;

        if strcmp(mode,'cycle')
            aN = nan(nCtrl,1);  cN = nan(nCtrl,1);
            for j = 1:nCtrl
                mc = getWallMotion(Xc{j},cuts{i},sM,m.keep);
                aN(j) = mc.wallAmp;  cN(j) = mc.centreAmp;
            end
            tree.nullPx(i)       = median(aN,'omitnan');
            tree.nullSpreadPx(i) = max(aN,[],'omitnan')-min(aN,[],'omitnan');
            tree.centreNullPx(i) = median(cN,'omitnan');
        else
            tree.nullPx(i)       = m.wallNullAmp;
            tree.nullSpreadPx(i) = NaN;      % one band, so there is no spread to show
            tree.centreNullPx(i) = m.centreNullAmp;
        end
        tree.conf(i)       = tree.wallPxRaw(i)  /max(tree.nullPx(i),      eps);
        tree.centreConf(i) = tree.centrePxRaw(i)/max(tree.centreNullPx(i),eps);

        if toc(tSweep)>=tNext
            el = toc(tSweep);
            rep.emit(sprintf(['  measuring the walls: %d of %d segments, %.0f%%, ' ...
                'about %s left'],c,numel(todo),100*c/numel(todo), ...
                hms(el*(numel(todo)-c)/c)));
            tNext = el+reportProgressEvery;
        end
    end

    % ---- the centre line's own geometry, so a segment that reports nothing can be
    % told apart from one that was never offered a cut ---------------------------
    tree.nCutOffered = nOffered(:);
    tree.pathPx      = [q.pathPx].';
    tree.clr         = [q.clr].';
    tree.spurFrac    = [q.spurFrac].';
    tree.why         = {q.why}.';

    % ---- the three gates, kept apart so a reader can see WHICH one refused ------
    % The calibre is the SEGMENTATION's own - sMetrics.diameter, twice the mean
    % distance transform along the centre line - and never the wall fit's.
    tree.diameterMaskUm = double(T.diameter)*s.pixelSize;
    tree.gateSize = tree.diameterMaskUm >= s.minDiameterUm;
    tree.gateCuts = tree.nCut           >= s.minCuts;
    tree.gateConf = tree.conf           >= s.confMin & tree.cohere >= s.minCohere;
    tree.cleared  = tree.gateSize & tree.gateCuts & tree.gateConf & isVes;
    tree.pixelSize= s.pixelSize;
    tree.mode     = mode;
    tree.nControls= nCtrl;
    if strcmp(mode,'band'), tree.bandHz = s.bandHz; tree.offBandHz = sM.offBandHz; end

    % ---- the table: the four columns and their two companions ------------------
    % NaN ACROSS ALL OF THEM, NOT A SMALL NUMBER IN ONE.  A wall amplitude of 0.004 px
    % beside a confidence of 1.3 is a measurement to nine readers out of ten, and
    % almost every row of a real segment table is that kind.  What a refused segment
    % measured is kept in the sub-tree, where a number cannot be mistaken for a result.
    drop = ~tree.cleared;
    wallPx   = tree.wallPxRaw;    wallPx(drop)   = NaN;
    centrePx = tree.centrePxRaw;  centrePx(drop) = NaN;
    conf     = tree.conf;         conf(drop)     = NaN;
    nCutCol  = tree.nCut;         nCutCol(drop)  = NaN;
    dilRel   = 100*wallPx./(tree.diameterPx/2);       % both walls move, two twos cancel
    share    = wallPx./(wallPx+centrePx);

    results.sMetrics.wmWallPx   = wallPx;
    results.sMetrics.wmWallUm   = wallPx*s.pixelSize;
    results.sMetrics.wmDilRel   = dilRel;
    results.sMetrics.wmCentrePx = centrePx;
    results.sMetrics.wmCentreUm = centrePx*s.pixelSize;
    results.sMetrics.wmDilShare = share;
    results.sMetrics.wmConf     = conf;
    results.sMetrics.wmCuts     = nCutCol;
    results.wallMotion          = tree;

    % ---- the report page, and then the optional pair of movies -----------------
    writeWallMotionPage(rep,results,tree,prov,fidx,s);
    if s.writeVideo && strcmp(mode,'cycle')
        writeVideoPair(rep,s,X,Xc{1},results,tree);
    elseif s.writeVideo
        rep.emit(['  no movie: magnifying a continuous recording shows its own noise ' ...
                  'at the same size as its signal']);
    end

    settings.runMotionEnhancement = reportSettings(s);
    reportWriting(rep);
    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
    save(s.fName,'results','-v7.3','-nocompression');
    reportSaved(rep);
    %Before the next file allocates its own: a plain re-assignment would build the new
    %arrays while these are still resident, and on a full field they are the largest
    %pair this step holds.
    clear X Xc source
    end
end
reportClose(rep);

end

% =====================================================================
function s = withDefaults(s)
%withDefaults  Everything this step has, defaulted where a hand-written s omits one.
%   The registry preset supplies all of them, so this is for a script that names only
%   what it wants changed.  pixelSize and bandHz are NOT here: both are facts about a
%   rig or a preparation and a default for either is a silent wrong answer.
d = struct('cutSpacing',3,'endMarginRadii',1.5,'spanRadii',3,'spanPadUm',30, ...
    'smoothRadii',1,'smoothMin',3,'interpF',4,'maxCLR',Inf,'minCuts',3, ...
    'cutStepPx',0.25,'widthTol',[0.5 2],'edgeTol',1.0, ...
    'minDiameterUm',37,'confMin',2,'minCohere',0, ...
    'sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3,'difR',0.4, ...
    'parforSegmentationLabels',true, ...
    'writeVideo',false,'videoFolder','','alpha',[53.1 106.2 212.5],'levels',3:5);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s,f{i}) || isempty(s.(f{i})), s.(f{i}) = d.(f{i}); end
end
if ~isfield(s,'bandHz'), s.bandHz = []; end
end

% =====================================================================
function m = modeOf(fName)
%modeOf  Which arm this product takes, from its stage token.
%   THE CLOSED CYCLE IS THE ONLY THING THAT HAS A FUNDAMENTAL.  '_c_I' is one beat
%   whose last sample is adjacent to its first; '_a_I' and '_b_I' are frame series
%   with a clock.  Reading the name is this function's job, which is why nothing else
%   below does it.
[~,nm] = fileparts(char(fName));
if     contains(nm,'_c_I_'), m = 'cycle';
elseif contains(nm,'_a_I_') || contains(nm,'_b_I_'), m = 'band';
else
    error('runMotionEnhancement:unknownProduct', ...
        ['%s is not an angiogram, a cardiac cycle or a bolus, so this step cannot ' ...
         'tell whether it is one averaged beat or a series of frames.'], nm);
end
end

% =====================================================================
function b = offBandOf(band,fs)
%offBandOf  A band of the same width as the rhythm's, where no physiology sits.
%   Centred half way between the fundamental and its first harmonic.  SAME WIDTH
%   matters: a wider band collects more noise and a narrower one less, and the
%   comparison would then be between two different measurements rather than between
%   a signal and its floor.  This is meProbe's own continuous control, and meProbe is
%   where the 1.0-times-its-own-control verdict on the continuous mode came from.
w  = diff(band);
f0 = mean(band);
b  = [1.5*f0-w/2, 1.5*f0+w/2];
b  = [max(b(1),band(2)+0.1), min(b(2),0.45*fs)];
if b(2)<=b(1)
    error('runMotionEnhancement:noOffBand', ...
        ['A band of %g to %g Hz leaves no room beside it at %g frames per second, ' ...
         'so there is nowhere to measure this recording''s own noise floor.'], ...
        band(1),band(2),fs);
end
end

% =====================================================================
function checkSweepFits(nCut,nSample,nCtrl,mode,s)
%checkSweepFits  REFUSE BEFORE THE SWEEP, and say what would have worked.
%   A wall is fitted once per cut per sample, on the product and on every control, and
%   a fit is a profile smooth and six interpolated crossings.  On one averaged cycle
%   that is tens of bins and the sweep is seconds; on a continuous recording of thirty
%   thousand frames the SAME estimator over the same cuts is hours, and nothing about
%   the settings says so.  The bound below is this library's own: a user processes a
%   whole recording in under twenty minutes.
FITSPERSEC = 2.5e5;                 % measured order of magnitude, not a promise
MAXSECONDS = 20*60;
nFit = nCut*nSample*(1 + nCtrl*strcmp(mode,'cycle'));   % the product AND every control
if nFit/FITSPERSEC <= MAXSECONDS, return; end
nFitOK = MAXSECONDS*FITSPERSEC;
error('runMotionEnhancement:sweepTooLarge', ...
    ['Measuring %d cuts on %d frames would fit walls about %.0f million times, ' ...
     'which is hours rather than minutes.  A larger smallest vessel (s.minDiameterUm ' ...
     'is %g um) or a wider cut spacing (s.cutSpacing is %g px) offers fewer cuts; ' ...
     'about %.0f million fits is what finishes inside twenty minutes.'], ...
    nCut, nSample, nFit/1e6, s.minDiameterUm, s.cutSpacing, nFitOK/1e6);
end

% =====================================================================
function tree = blankTree(n,ids)
%blankTree  The sub-tree, one row per row of sMetrics, so a reader joins by row.
f = {'wallPxRaw','centrePxRaw','nullPx','nullSpreadPx','centreNullPx','conf', ...
     'centreConf','cohere','centreCohere','nCut','diameterPx','diameterMaskUm', ...
     'taperPx','edgeUm','nCutOffered','pathPx','clr','spurFrac'};
tree = struct();
tree.idx = double(ids(:));
for k = 1:numel(f), tree.(f{k}) = nan(n,1); end
tree.cleared = false(n,1);
end

% =====================================================================
function prov = readSegmentationProvenance(fNames)
%readSegmentationProvenance  What each file's sidecar says it was segmented with.
%   Judged from the sidecars alone and before any array is read, so a product that was
%   never segmented is refused by name here rather than on the line that reads cMask,
%   and a batch segmented at two different settings is named on every page.
prov = struct('edgeSize',{},'sMinL',{},'file',{});
for i = 1:numel(fNames)
    if isempty(fNames{i}), prov(i) = struct('edgeSize',NaN,'sMinL',NaN,'file',''); continue; end
    p = getProductPath(fNames{i},'s');
    if ~isfile(p)
        error('runMotionEnhancement:noSettings', ...
            'There is no settings file beside %s.', shortName(fNames{i}));
    end
    S = load(p,'settings');
    if ~isfield(S,'settings') || ~isfield(S.settings,'runSegmentation')
        error('runMotionEnhancement:notSegmented', ...
            ['%s has not been segmented, so there are no vessels to measure a wall ' ...
             'on.  Run Segmentation on it first.'], shortName(fNames{i}));
    end
    prov(i) = struct('edgeSize',S.settings.runSegmentation.edgeSize, ...
                     'sMinL',getfielddef(S.settings.runSegmentation,'sMinL',NaN), ...
                     'file',shortName(fNames{i}));
end
end

% =====================================================================
function v = getfielddef(st,f,d)
%getfielddef  One field of a settings struct, or a stated default.
if isfield(st,f) && ~isempty(st.(f)), v = st.(f); else, v = d; end
end

% =====================================================================
function writeWallMotionPage(rep,results,tree,prov,fidx,s)
%writeWallMotionPage  The field, the four columns against calibre, and the refusals.
%
%   THE COUNT OF ROWS THIS STEP REFUSED IS ON THE PAGE, SPLIT BY GATE.  A page that
%   shows only what cleared is the page that made this measurement look easy: on the
%   reference field 24 of 1 422 segments carry a number and the other 1 398 are the
%   result too.
%
%   THE HEADLINE GOES IN THE SUBTITLE AND NOT IN A TITLE.  reportSave writes the
%   recording's name across every page with sgtitle, which attaches to whichever grid
%   the figure holds - so a title set here would be silently overwritten and the page
%   would export with the numbers simply missing.  And a text object does not wrap
%   itself on a fixed canvas, so the lines are broken here or they lose both ends.
fh = reportFigure(rep,'wall-motion');
tl = tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

cl  = tree.cleared;
dUm = tree.diameterMaskUm;

% ---- 1. the field, with the segments that carry a number drawn on it ----------
ax = nexttile(tl,1);
imagesc(ax,results.imgI)
clim(ax,[prctile(results.imgI(:),1),prctile(results.imgI(:),99)])
colormap(ax,gray); axis(ax,'image'); hold(ax,'on')
if any(cl)
    ok = ismember(double(results.sMap),tree.idx(cl));
    [B,~] = bwboundaries(ok,'noholes');
    for b = 1:numel(B), plot(ax,B{b}(:,2),B{b}(:,1),'-','Color',[1 0.3 0.1],'LineWidth',1); end
end
hold(ax,'off')
title(ax,sprintf('%d of %s carry a wall measurement',nnz(cl), ...
    plural(nnz(isfinite(dUm)),'vessel','vessels')))

% ---- 2. the two displacement columns against calibre --------------------------
% THE PANEL A PAGE WITH NOTHING ON IT STILL HAS TO DRAW, and on this measurement that
% is the COMMON case: on the reference field 24 of 1 422 segments carry a number and on
% a continuous product none of them do.  Empty data makes semilogy return a 0x1 handle
% array, so a legend then has two labels and no series and MATLAB quietly renames them
% - one NaN point keeps the axes empty and the legend honest.  The two series are
% handed to legend BY HANDLE for the same reason: the gate line is a third object on
% these axes and a legend given only labels matches them in draw order.
ax = nexttile(tl,2);
[xW,yW] = orNaN(dUm(cl), tree.wallPxRaw(cl)*s.pixelSize);
[xC,yC] = orNaN(dUm(cl), tree.centrePxRaw(cl)*s.pixelSize);
hW = semilogy(ax,xW,yW,'o','MarkerSize',4, ...
    'MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','none'); hold(ax,'on')
hC = semilogy(ax,xC,yC,'s','MarkerSize',4, ...
    'MarkerFaceColor',[0.85 0.5 0.1],'MarkerEdgeColor','none');
xline(ax,s.minDiameterUm,'k--');
hold(ax,'off'); grid(ax,'on')
padX(ax,dUm(cl),s.minDiameterUm);
xlabel(ax,'Vessel width, micrometres'); ylabel(ax,'Movement, micrometres')
legend(ax,[hW hC],{'each wall','the whole vessel'},'Location','northwest','Box','off')
title(ax,'How far the walls moved, and how far the vessel slid')

% ---- 3. the fractional change and the split -----------------------------------
ax = nexttile(tl,3);
yyaxis(ax,'left')
plot(ax,dUm(cl),tree.wallPxRaw(cl)./max(tree.diameterPx(cl)/2,eps)*100,'o', ...
    'MarkerSize',4,'MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','none')
ylabel(ax,'Width change over the beat, %')
yyaxis(ax,'right')
plot(ax,dUm(cl),tree.wallPxRaw(cl)./max(tree.wallPxRaw(cl)+tree.centrePxRaw(cl),eps), ...
    's','MarkerSize',4,'MarkerFaceColor',[0.85 0.5 0.1],'MarkerEdgeColor','none')
ylim(ax,[0 1]); ylabel(ax,'Share that is width change, not sliding')
xline(ax,s.minDiameterUm,'k--'); grid(ax,'on')
padX(ax,dUm(cl),s.minDiameterUm);
xlabel(ax,'Vessel width, micrometres')
title(ax,'How much the width changed, and how much of it was width at all')

% ---- 4. what was refused, and by which gate -----------------------------------
ax   = nexttile(tl,4);
isV  = isfinite(dUm);
lab  = {'too narrow','too few cuts','no more than its own floor','carries a number'};
cnt  = [ nnz(isV & ~tree.gateSize)
         nnz(isV &  tree.gateSize & ~tree.gateCuts)
         nnz(isV &  tree.gateSize &  tree.gateCuts & ~tree.gateConf)
         nnz(cl) ];
barh(ax,cnt); grid(ax,'on')
set(ax,'YTick',1:4,'YTickLabel',lab)
xlabel(ax,'Vessels')
title(ax,'Why a vessel reports nothing')

notes = wrapLines([{headline(tree,s)}, pageNotes(prov,fidx)],230);
subtitle(tl,notes,'FontSize',8,'Interpreter','none');

reportSave(rep,fh,'wall-motion');
end

% =====================================================================
function [x,y] = orNaN(x,y)
%orNaN  One invisible point instead of nothing, so a series always has a handle.
if isempty(x), x = NaN; y = NaN; end
end

% =====================================================================
function padX(ax,x,gate)
%padX  A margin either side, so a marker at the extreme is not half off the axes.
%   With one cleared vessel MATLAB's automatic limits land exactly on its calibre and
%   the marker is clipped by the box - which reads as a missing point rather than an
%   edge case.
lo = min([x(:); gate]);  hi = max([x(:); gate]);
if ~isfinite(lo) || ~isfinite(hi), return; end
pad = max(0.05*(hi-lo), 1);
xlim(ax,[lo-pad hi+pad]);
end

% =====================================================================
function t = headline(tree,s)
%headline  One line a reader can quote, and it names the refusals as well as the wins.
cl  = tree.cleared;
nV  = nnz(isfinite(tree.diameterMaskUm));
if ~any(cl)
    t = sprintf(['No vessel on this recording cleared all three gates, so every wall ' ...
        'column is blank. That is a result: it says the movement here is no larger ' ...
        'than what the same vessels show with the heartbeat''s timing scrambled ' ...
        'away. %s were measured.'], plural(nV,'vessel','vessels'));
    return
end
wallUm = tree.wallPxRaw(cl)*s.pixelSize;
dilRel = 100*tree.wallPxRaw(cl)./max(tree.diameterPx(cl)/2,eps);
share  = 100*median(tree.wallPxRaw(cl)./max(tree.wallPxRaw(cl)+tree.centrePxRaw(cl),eps));
% A range with one member is not a range, and printing "0.261 to 0.261" reads as an
% arithmetic mistake rather than as a single vessel.
if nnz(cl)==1
    t = sprintf(['1 of %s cleared all three gates: its walls moved %.3f um peak to ' ...
        'peak, its width changed %.2f %%, and %.0f %% of its movement was width ' ...
        'change rather than sliding.'], plural(nV,'vessel','vessels'), ...
        wallUm, dilRel, share);
else
    t = sprintf(['%d of %s cleared all three gates: walls moved %.3f to %.3f um peak ' ...
        'to peak, width changed %.2f to %.2f %%, and the width change was %.0f %% of ' ...
        'the total movement in the median one.'], nnz(cl), plural(nV,'vessel','vessels'), ...
        min(wallUm), max(wallUm), min(dilRel), max(dilRel), share);
end
end

% =====================================================================
function t = plural(n,one,many)
%plural  "1 vessel" rather than "1 vessels", on a page a biologist reads.
if n==1, t = sprintf('%d %s',n,one); else, t = sprintf('%d %s',n,many); end
end

% =====================================================================
function c = pageNotes(prov,fidx)
%pageNotes  The sentences that say what these numbers are NOT.
c = {['The smallest vessel this can see, ' ...
      'and the way the vessels were labelled, are stated limits rather than ' ...
      'measurements of this recording. The smallest calibre with a signal above its ' ...
      'own scrambled-timing copy was 36.9 micrometres, measured once, on one ' ...
      'recording, on one animal, and it does not improve with a finer pixel - the ' ...
      'same vessels measured at half the pixel size gave the same answer.'], ...
     ['A vessel is compared only with itself. The floor each one is divided by is ' ...
      'the same measurement on the same vessel with the heartbeat''s timing ' ...
      'destroyed, so it already contains that vessel''s own texture, its red cells ' ...
      'and whatever the preparation was doing.'], ...
     ['How far the vessel slid is reported beside how far its walls moved because ' ...
      'the two are the same size here. A width change quoted without it can be the ' ...
      'preparation moving.']};
if ~isempty(prov) && fidx<=numel(prov)
    es = [prov.edgeSize];
    if numel(unique(es(isfinite(es))))>1
        c{end+1} = sprintf(['This batch was segmented at more than one setting ' ...
            '(this recording used %g); the vessel outlines are not comparable ' ...
            'across it.'], prov(fidx).edgeSize);
    end
end
end

% =====================================================================
function c = wrapLines(c,n)
%wrapLines  Break each line on a word boundary so nothing runs off the page.
%   The canvas is a fixed 1400 px and a text object does not wrap itself: unwrapped,
%   a long line loses BOTH of its ends, which for a refusal count is a message nobody
%   can read.
out = cell(1,0);
for i = 1:1:numel(c)
    t = c{i};
    while numel(t)>n
        k = find(isspace(t(1:n+1)),1,'last');
        if isempty(k), k = n; end
        out{end+1} = strtrim(t(1:k-1)); %#ok<AGROW> - a handful of lines
        t = strtrim(t(k:end));
    end
    out{end+1} = t;                     %#ok<AGROW>
end
c = out;
end

% =====================================================================
function writeVideoPair(rep,s,X,Xc,results,tree)
%writeVideoPair  The magnified cycle and its matched control, or neither of them.
%
%   THE TWO FILES SHIP TOGETHER.  A magnified movie is persuasive whether or not it
%   is true, so the control - the same beats with their timing scrambled, magnified
%   at the same alpha through the same pyramid - is written first and named on the
%   product's own banner.  If the pair cannot both be written, neither is.
%
%   AND THE FOLDER IS NEVER THE RESULTS TREE.  These are 0.1 to 0.8 GB each; the
%   workbench and the explorer sweep the results tree by name, and a file that large
%   which is not a product has no business being swept.  The default is a sibling of
%   it, and the settings tip says what it costs on a synced drive.
outDir = s.videoFolder;
if isempty(outDir)
    outDir = fullfile(fileparts(fileparts(s.fName)),'motionEnhancement');
end
if ~isfolder(outDir), mkdir(outDir); end
% The role letter is stripped from the END only, so a recording whose own stem ends in
% one keeps it - the same anchoring getProductPath uses.
[~,stem] = fileparts(s.fName);
stem = regexprep(stem,'_r$','');

sv          = meSettings;
sv.alpha    = s.alpha;
sv.levels   = s.levels;
sv.maskMode = 'weight';
sv.maskSource  = 'segmentation';
sv.maskVessel  = results.cMask>3;      % the library's own segmentation of this field
sv.mask     = meMask(sv,double(results.imgI));
nBin        = size(X,3);
sv.passband = [0.5, min(sv.cineHarmonics, floor(nBin/2)-0.5)+0.5];
sv.reportLevel = 'quiet';

cl   = tree.cleared;
if any(cl)
    calibre = sprintf('valid for vessels of about %.0f to %.0f micrometres', ...
        min(tree.diameterMaskUm(cl)), max(tree.diameterMaskUm(cl)));
else
    calibre = sprintf('no vessel on this recording cleared its own floor at %g um', ...
        s.minDiameterUm);
end
meta = struct('stem',stem,'alpha',sv.alpha,'passband',sv.passband, ...
    'passbandUnit','cycles/cycle','axis','phase', ...
    'axisValues',(0:nBin-1).'/nBin,'pixelSize',s.pixelSize,'masked',true, ...
    'calibre',calibre);

ctrlName = [stem '_wallMotion_control.avi'];
prodName = [stem '_wallMotion.avi'];

rep.emit('  magnifying the comparison beats');
mc = meShift(Xc,nBin,sv);
metaC = meta;
metaC.outPath = fullfile(outDir,ctrlName);
metaC.product = 'the same beats with their timing scrambled';
metaC.isNull  = true;
metaC.nullText= 'the heartbeat''s timing has been destroyed; anything that moves here is not the heartbeat';
meWriteVideo(sv,struct('original',Xc,'magnified',mc),metaC);

rep.emit('  magnifying the averaged beat');
mp = meShift(X,nBin,sv);
metaP = meta;
metaP.outPath     = fullfile(outDir,prodName);
metaP.product     = 'the averaged heartbeat';
metaP.isNull      = false;
metaP.controlName = ctrlName;
meWriteVideo(sv,struct('original',X,'magnified',mp),metaP);

rep.emit(sprintf('  wrote %s and its control to %s',prodName,outDir));
end

% =====================================================================
function t = hms(sec)
%hms  Seconds as a human duration, for a progress line.
if ~isfinite(sec) || sec<0, t = 'a while'; return; end
if sec<90, t = sprintf('%.0f s',sec);
elseif sec<5400, t = sprintf('%.0f min',sec/60);
else, t = sprintf('%.1f h',sec/3600);
end
end

% =====================================================================
function n = shortName(p)
%shortName  A file's name without its folder, for a message a person reads.
[~,b,e] = fileparts(char(p));
n = [b e];
end
%------------- END OF CODE --------------
