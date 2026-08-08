%runNVC  Per-epoch neurovascular-coupling response into the results.nvc tree
%
%   runNVC(s,fNames) loads every contrast-branch *_t_BFI_r.mat (or *_s_BFI_r.mat)
%   whole-recording product in fNames (runContrastFromRLS -> runSegmentation ->
%   runBFI: results.sData/gsData/dvsData/dvsDiameter on results.time, plus the flow
%   cube source.data [Y x X x nT]), cuts every stimulus repetition out of the
%   recording, and measures each one ON ITS OWN.
%
%   NOTHING IS AVERAGED BEFORE IT IS MEASURED, and that is the whole point.  The step
%   this replaces built a second product triplet whose only content was the average of
%   twenty repetitions, and every number downstream described that average - so
%   habituation, drift and a single bad repetition were all invisible, by construction.
%   Here the epochs stay separate: results.nvc.esMetrics.<signal>.<name> is
%   [nSeg x nEp], every segment, every epoch, every metric, and the average is
%   something a reader takes afterwards rather than something the pipeline took first.
%   4.7 GB per recording of derived data stops existing with it.
%
%   THE STEP RUNS STRICTLY AFTER THE BFI CONVERSION, and that is load-bearing rather
%   than merely an ordering.  Every trace it reads is already a flow index, so NO
%   1/K^2 APPEARS ANYWHERE IN THIS FILE OR IN ITS CORE, and nothing here asks what
%   format a trace is in - the '_BFI_' token in the product name is the answer.  A
%   '_K_' product is refused by name even though the registry cannot offer one, since
%   a launcher file is written by hand.
%
%   IT NEVER READS THE IMAGE CUBE FOR THE PER-SEGMENT PATH.  The epoch-rejection rules
%   are computed from the recording's own mean trace, so the 21 GB source member and
%   the 31.6 s it takes to load are not on this step's path at all.  Only the optional
%   per-pixel path (s.ppxNvcReturn) opens it, and it is off by default.
%
%   TWO SSIM IMAGE RULES WENT WITH THE AVERAGE THEY PROTECTED (author, 2026-08-05).
%   An epoch that moved must not be added to the other nineteen - that is what the
%   image-similarity tests were for.  Nothing is added to anything any more, so the
%   comparison between epochs has no consumer, and rejectBlSimCoef / rejectSimCoef are
%   gone rather than reimplemented on the traces.  A per-epoch outlier in the metrics
%   themselves is the more direct signal now that every epoch is kept.
%
%   NOTHING IS DELETED BY A REJECTION EITHER.  results.nvc.epochValid is a FLAG:
%   every epoch's metrics are computed and stored whatever the rules said, and the
%   flag is what the aggregations honour.  In the step this replaces a rejected epoch
%   simply never entered the average and left no trace of having existed.
%
%   OUTPUT TREE  (RESULTS.nvc; every float SINGLE)
%     Shared axes - the EPOCH clock every metric's times are quoted against, and the
%     stimulus on it, so a downstream plot is self-describing:
%       .time          [nES x 1]  the epoch clock, 0 .. one epoch, seconds
%       .stimulus      [nES x 1]  the stimulus boxcar on that clock
%       .epochStart    [nEp x 1]  where each epoch was cut, on the RECORDING clock, s
%       .epochValid    [nEp x 1]  logical - no quality rule rejected this epoch
%       .epochReject   [nRule x nEp] logical - ONE ROW PER RULE, in the field order of
%                                 .epochQuality below, which is the only place that
%                                 order is written down
%       .epochQuality  struct, each field [nEp x 1]: blDev finDev meanDev peakDev
%                                 timeLoss.  EACH FIELD IS THE QUANTITY ITS RULE
%                                 TESTS, IN THE UNITS OF ITS OWN COEFFICIENT - the
%                                 first four in standard deviations of the spread
%                                 ACROSS epochs, timeLoss in seconds - so a rule and
%                                 the number it fired on are read on the same scale
%       .esMetrics     <- THE PRODUCT.  One sub-struct per analysed signal, each
%                         field [nSeg x nEp], row order = segment order (aligned to
%                         sMetrics / dvsMetrics rows), column order = epoch order:
%           .sData        .<ns*>      .gsData       .<ng*>
%           .dvsData      .<ns*>      .dvsDiameter  .<nd*>
%       .agg           the AGGREGATE FIT, one fit per segment on that segment's
%                      epoch-average over the VALID epochs, [nSeg x 1] per scalar:
%           .<signal>.<prefix><Name>
%       .ppx           (only when s.ppxNvcReturn is non-empty) [Y x X x nEp] per
%                      metric.  MARKERS ONLY - see below
%     esMetrics IS SPLIT BY SIGNAL because the segment sets differ: sData has 1940
%     rows on the reference recording and dvsDiameter has 67, and they cannot share
%     one array.  The PREFIXES keep the existing vocabulary - ns = NVC flow,
%     nd = NVC diameter, ng = NVC guided contrast - and an nd* column must never be
%     pooled with an ns* one: dvsDiameter is a length in pixels, not a flow index.
%
%   THE EPOCH-CUT TRACES THEMSELVES ARE NOT STORED, and that is a decision rather than
%   an omission.  [nES x nSeg x nEp] is 18.6 GB on the reference recording, the core's
%   levels are {'markers','timing','areas'} so nothing can ask for it, and the two
%   curves a reader actually wants - the epoch-average and its fit - are on the
%   response page.  Flagged for the consumers session.
%
%   METRICS-TABLE DUPLICATION (a small key set, same names as the tree, single)
%     results.sMetrics    (from sData)      : nsPeakRel nsTTP nsAUCf nsAUCblkRel
%                                             nsSNR nsNEpochs
%     results.dvsMetrics  (from dvsData)    : the same six ns*
%                         (from dvsDiameter): the six nd*
%     gsData writes NO table - there is no gMetrics, exactly as in runVasomotion.
%     THE AGGREGATE IS A MEDIAN OVER THE VALID EPOCHS, NEVER A MEAN, and nsSNR ships
%     beside it for a reason that is measured rather than stylistic.  The per-epoch
%     amplitudes are read off the RAW trace, so that a peak and its time name the same
%     sample of the data (getNVCMetrics' header argues this in full); a max over ~80
%     noisy samples is consequently a biased estimator, sitting about 2.4 baseline SDs
%     above baseline on an epoch with no stimulus in it at all against 3.96 on a real
%     one.  What removes that bias is AGGREGATION, not smoothing - hence the median -
%     and what makes it readable rather than hidden is SNR travelling beside every
%     amplitude.  nsNEpochs counts, per segment, how many valid epochs the median was
%     actually taken over: without it a median of one epoch and a median of twenty are
%     the same number.
%     BFI / std(BFI) are runBFI's columns and are LEFT UNTOUCHED.
%
%   THE TIMING METRICS ARE AGGREGATION INPUTS, NOT PER-EPOCH READINGS.  Measured on
%   real 4 Hz data, the fraction of a metric's variance that belongs to the SEGMENT
%   rather than to the trial is 0.61 for Peak, 0.49 for AUCblk and 0.42 for AUCf, but
%   0.18 for TRise, 0.13 for TTP and 0.06 for TDec.  Neither a longer median filter
%   nor fitting the epoch changes that - it is the SNR of one repetition, not a fault
%   in the definitions.  This is why the trials page plots the AMPLITUDE metrics and
%   why the metrics tables aggregate before they report.
%
%   LEVELS  (selector cell arrays; a subset of {'markers','timing','areas'})
%     s.segNvcReturn  per-segment levels; default all three.
%     s.ppxNvcReturn  per-pixel levels + GATE: NON-EMPTY runs the per-pixel path,
%                     [] skips it; DEFAULT [] - OFF.  'markers' is the per-pixel set:
%                     amplitude is spatially structured (rho -0.39 against distance
%                     from the activation focus) and per-epoch timing is not
%                     (rho -0.04).
%     THE LEVELS GATE THE TREE, NOT THE METRICS TABLE.  The table's key set spans the
%     marker and the area blocks, both of which cost no fit, so both are always
%     computed; the user's selection still decides what reaches the saved tree.  Same
%     rule as runFitVasoreactivity and runPulsatility, and it keeps a metrics column
%     from appearing and disappearing with a display setting.
%
%   THE TWO FITS, AND WHAT EACH COSTS
%     THE AGGREGATE FIT IS UNCONDITIONAL.  One fitNVC per segment on that segment's
%     epoch-average - 1.3 min on 8 workers for the 1938 segments of the reference
%     recording - and it is where the two-timescale structure the data really has can
%     be estimated at all, because it is the only place the SNR supports it.
%     THE PER-EPOCH FIT IS NOT, and s.nvcFit (default false) gates it.  IT COSTS
%     26.7 MIN PER RECORDING ON 8 WORKERS for sData alone, twenty times the aggregate,
%     and it buys less than it looks: fitting each epoch individually moves TTP's
%     between-segment share from 0.086 to 0.096, because the signal it is chasing is
%     itself only about 1.75 s wide.  It IS a better estimator of an individual epoch
%     (the p90 TTP outlier falls from 11.1 s to 7.0 s) and that is what it is for.
%     THERE IS NO PER-PIXEL FIT AT ANY SETTING.  23.2 M fits is 266 hours.
%
%   THE SCALAR SET OF EITHER FIT IS THE CORE'S, NOT THIS FILE'S.  Both fits are driven
%   by layout.scalarNames, so fitNVC can add or drop a scalar - Identified, say -
%   without a line changing here.  Both are called at the 'model' level only, which
%   leaves out the reconstruction: a fitted CURVE is [nES x 1] per trace and has no
%   place in a [nSeg x nEp] array of scalars.  The response page asks for it separately
%   for the one curve it draws.
%   AND EITHER FIT CAN REPORT THAT IT MEASURED NOTHING.  fitNVC emits Identified, 0
%   when t0 or a time constant came to rest on a bound, and it emits nothing for
%   Zeta / Omega0 / RingPeriod / Damping when that happens, because all four are
%   functions of the railed constant.  On real 4 Hz data that is 56 % of segments even
%   on the epoch-average, and more per epoch - so those four columns are mostly the
%   aggregate's, and where they carry a number that number is backed by the trace.
%
%   POLARITY IS RESOLVED ONCE PER SEGMENT, FROM ITS EPOCH-AVERAGE.  Constricting
%   vessels and the negative surround response are real, so s.nvcPolarity 'auto'
%   (the default) takes the sign of the mean of the filtered epoch-average minus its
%   baseline over the response window, and hands that one sign to every epoch of that
%   segment.  A sign resolved per EPOCH would flip on noise in a segment that does not
%   respond, which makes Peak the maximum of |noise| - the null distribution the
%   median aggregate exists to keep out.  'positive' / 'negative' override both.
%
%   INPUTS
%     s        parameter struct with fields
%              % where the epochs are
%                • stimStartType   'offset' | 'manual'
%                • stimOffset      seconds from the recording start to the FIRST
%                                  epoch's start ('offset' mode)
%                • stimStartAll    per-file 'HH:mm:ss.SSS' list ('manual' mode)
%                • epochsN         how many epochs
%                • epochDurationSec       one epoch's length, s
%                • epochBaselineSec       [t1 t2] baseline window on the epoch clock
%                • epochStimStartSec      when the stimulus starts on the epoch clock
%                • stimDurationSec        how long it lasts, s.  REQUIRED - this step
%                                  USES it, unlike the epoching step it replaces,
%                                  which merely recorded it for someone else
%                • epochFinaleSec  [-t1 0] finale window, back from the epoch's END
%              % how a response is measured (passed to getNVCMetrics)
%                • nvcMedFiltSec   median-filter length for the CROSSING RULES, in
%                                  SECONDS (default 1.25).  Seconds and not samples
%                                  because sData is 4 Hz and gsData is 100 Hz in the
%                                  same run
%                • nvcPolarity     'auto' (default) | 'positive' | 'negative'
%                • nvcHoldSamples  consecutive samples a crossing must hold (2)
%                • nvcDecayFrac    fraction-of-peak decay level (0.1)
%                • nvcSignals      which traces to analyse (default all four)
%                • segNvcReturn    per-segment level cell (default all three)
%                • ppxNvcReturn    per-pixel level cell / gate (default [] = off)
%              % which epochs to trust
%                • rejectBlCoef rejectFinCoef rejectEpochCoef rejectPeakCoef
%                                  how far an epoch's baseline / finale / whole-epoch
%                                  mean / peak may sit from the median across epochs,
%                                  in SDs of that spread (Inf disables the rule)
%                • rejectTimeLoss  seconds of dropped frames an epoch may carry
%                • rejectFirstEpoch  1 = always reject the first epoch
%                • enablelRejectionModification  1 = open the accept/reject picker
%              % fitting
%                • nvcFit          false (default) | true - ALSO fit every epoch
%                • nvcModel        'secondOrder' (default) | 'doubleGamma'
%                • nvcStarts       multi-start start points (default 16)
%                • parforNvcSegments logical, default true
%                • parforNvcPixels   logical, default true - both are WORKER BOUNDS,
%                                  not branches: false runs the identical loop body
%                                  serially in the client and starts no pool.
%     fNames   cell array of *_t_BFI_r.mat (or *_s_BFI_r.mat) paths, IN THE ORDER
%              s.stimStartAll is written in.
%                • Optional workbench hooks in s (no-op when absent):
%                  s.stageFcn(stage,detail), s.cancelFcn()->tf.  Cancel is checked
%                  between files (never inside a parfor).
%
%   OUTPUT FILES (side-effects) - NON-DESTRUCTIVE: the flow cube is preserved
%       *_t_BFI_d.mat   SOURCE   - NEVER re-saved, and read at all only by the
%                                  optional per-pixel path
%       *_t_BFI_r.mat   RESULTS  - RESULTS.nvc.* + the ns*/nd* metrics columns
%       *_t_BFI_s.mat   SETTINGS - field settings.runNVC added
%       *_rep_nvccuts.jpg     - the recording's mean trace, where every epoch was
%                               cut, and one accept/reject bar per rule
%       *_rep_nvcresponse.jpg - the epoch-average response, the aggregate fit, the
%                               residuals and the fitted parameters
%       *_rep_nvctrials.jpg   - the per-epoch amplitude across the session, so
%                               habituation or drift is visible at a glance.  This is
%                               the page that could not exist before
%
%   EXAMPLE
%     s.stimStartType='offset'; s.stimOffset=0;
%     s.epochsN=20; s.epochDurationSec=30;
%     s.epochBaselineSec=[0 10]; s.epochStimStartSec=10; s.stimDurationSec=5;
%     s.epochFinaleSec=[-5 0];
%     fNames = getFileNamesList(resultsFolder,'*_t_BFI_r.mat');
%     runNVC(s, fNames(:));
%
%   DEPENDS ON
%     Core/NVC/getNVCMetrics (the per-epoch science), Core/NVC/fitNVC (both fits),
%     MATLAB Optimization Toolbox (lsqcurvefit); core LSCI library utilities.
%
% See also: getNVCMetrics, fitNVC, runFitVasoreactivity, runVasomotion,
%           runPulsatility, getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 05-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Where the stimulus repetitions are
% s.stimStartType='offset'; % 'offset' for a stimulus that starts a fixed time after
%                         % the recording begins, 'manual' for a list of clock times.
%                         % Either way this is the start of the FIRST EPOCH, not of
%                         % the stimulus inside it.
% s.stimOffset=0;         % seconds ('offset' mode)
% s.stimStartAll{1}='09:23:31.346';  % 'HH:mm:ss.SSS', one per file ('manual' mode)
% s.epochsN=20;           % how many repetitions
% s.epochDurationSec=30;  % how long one repetition is, seconds
% s.epochBaselineSec=[0,10];  % the quiet stretch at the start of each epoch that the
%                         % response is measured against, seconds from the epoch start
% s.epochStimStartSec=10; % when the stimulus actually starts within the epoch
% s.stimDurationSec=5;    % how long it lasts, seconds
% s.epochFinaleSec=[-5,0];% the stretch at the END of the epoch where flow is expected
%                         % to be back to baseline, measured back from the epoch's end
% %ADJUSTED IF NECESSARY - How a response is measured
% s.nvcMedFiltSec=1.25;   % smoothing used ONLY to decide when the response rose and
%                         % when it came back, in seconds.  The peak and the areas are
%                         % read off the trace as recorded.
% s.nvcPolarity='auto';   % 'auto' lets each vessel respond in its own direction;
%                         % 'positive' or 'negative' forces one for all of them.
% s.nvcHoldSamples=2;     % how many samples in a row a level has to be held before it
%                         % counts as crossed
% s.nvcDecayFrac=0.1;     % the second, fixed decay level: this fraction of the peak
% s.nvcSignals={'sData','gsData','dvsData','dvsDiameter'};  % which traces to analyse
% %ADJUSTED IF NECESSARY - Which levels to store
% s.segNvcReturn={'markers','timing','areas'};  % markers (baseline, peak, SNR),
%                         % timing (rise, decay, duration), areas.  Default = all.
% s.ppxNvcReturn=[];      % [] = off (the default).  Non-empty measures every pixel of
%                         % the field in every epoch and stores a map per epoch.
% %ADJUSTED IF NECESSARY - Which epochs to trust
% s.rejectBlCoef=1;       % use Inf to disable rejection by this parameter
% s.rejectFinCoef=1;      % use Inf to disable rejection by this parameter
% s.rejectEpochCoef=1;    % use Inf to disable rejection by this parameter
% s.rejectPeakCoef=1;     % use Inf to disable rejection by this parameter
% s.rejectTimeLoss=0.5;   % allowed time loss due to grabbing failure, seconds per epoch
% s.rejectFirstEpoch=1;   % always reject the first epoch
% s.enablelRejectionModification=1;  % 1 = look at the epochs and change the decisions
% %ADJUSTED IF NECESSARY - Fitting
% s.nvcFit=false;         % true ALSO fits every epoch of every segment separately.
%                         % Expect about 27 minutes per recording on 8 workers, against
%                         % about 1.3 for the per-segment fit that always runs.
% s.nvcModel='secondOrder';
% s.nvcStarts=16;         % how many start points a fit is tried from.  Fewer is faster
%                         % and more likely to settle on a wrong answer.
% %ADJUSTED IF NECESSARY - Parallel execution
% s.parforNvcSegments=true;
% s.parforNvcPixels=true; % false measures one at a time in this MATLAB and starts no
%                         % parallel pool - slower, but it costs no worker processes.

%------------- BEGIN CODE --------------
function runNVC(s,fNames)

%THE INPUT IS A WHOLE RECORDING ON THE CONTRAST BRANCH.  An epoch average ('_e_') and
%a cardiac average ('_c_') both collapse the recording clock, so there is no timeline
%left to place a stimulus on; a '_K_' product is a contrast, and this step applies no
%conversion because runBFI already did.  The guard is written positively, naming the
%two stages that ARE a whole recording of a flow index.
if ~all( cellfun(@(x) isempty(x) || contains(x,'_t_BFI_r.mat') || ...
        contains(x,'_s_BFI_r.mat'), fNames(:)) )
    error('runNVC:badInput', ...
        ['One or more *non-empty* entries are not a whole-recording BFI product.  ' ...
         'This step cuts stimulus repetitions out of the RECORDING, so it needs the ' ...
         'whole trace, already converted to a flow index: a cardiac average ' ...
         '("_c_") carries no recording clock, and a contrast product ("_K_") has ' ...
         'not been converted.  List the right files with ' ...
         'getFileNamesList(resultsFolder,''*_t_BFI_r.mat'').']);
end

s=nvcDefaults(s);
%THE GEOMETRY IS CHECKED ONCE, BEFORE ANY FILE IS OPENED, so a mistyped protocol is a
%message rather than a partly-processed data set.  It is also checked HERE and nowhere
%upstream: the geometry is this step's own input and there is no producing step to
%inherit it from - the epoching step that used to record it does not exist any more.
requireGeometry(s,numel(fNames));

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  This step can block
% on the epoch picker, so cancel is only checked between files.
rep=reportOpen(s,'NVC response',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        %SOURCE is read only for the per-pixel path (and never written back).
        if ~isempty(s.ppxNvcReturn)
            load(getProductPath(s.fName,'d'),'source')
        end
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        if ~isfield(results,'time') || isempty(results.time)
            error('runNVC:noTime', ...
                '%s carries no results.time, so there is no recording clock to cut.', ...
                s.fName);
        end
        %The ABSOLUTE start of the recording is what places the stimulus, and only a
        %step that read the raw recording clock can supply it.  runInternalCycle writes
        %no timeStamp - its product is one AVERAGED cardiac period - so a '_c_' product
        %cannot be epoched by construction rather than by omission.
        if ~isfield(results,'timeStamp')
            error('runNVC:noTimeStamp', ...
                ['%s carries no results.timeStamp, so the stimulus cannot be placed ' ...
                 'on the recording clock.  Cut a contrast-branch recording ' ...
                 '(*_t_BFI_r.mat or *_s_BFI_r.mat); an averaged cardiac cycle has no ' ...
                 'timeline.'],s.fName);
        end

        sf=s;
        sf.epochStartSec=epochStarts(sf,results.timeStamp,fidx);

        %runNVC owns RESULTS.nvc entirely; rebuild it from scratch so no stale
        %sub-branch survives a re-run with different levels or a different geometry.
        if isfield(results,'nvc'), results=rmfield(results,'nvc'); end

        % =============================================================
        % The cuts, on the PRIMARY clock (results.time).  gsData gets its own below;
        % everything the tree publishes at the root is on this one.
        % =============================================================
        [cut,epochStart]=cutIndex(results.time,sf);
        nEp=size(cut,2);
        %TWO LAYOUTS, ONE TIME BASE.  layoutCalc is what the core is CALLED with - the
        %user's levels plus the two the metrics tables need, because those cost no fit
        %and a table column must not appear and disappear with a display setting.
        %layoutTree is the user's own selection, and its scalarNames (a subsequence of
        %layoutCalc's) is exactly what enters the saved tree.
        sCalc=sf; sCalc.segNvcReturn=withTableLevels(sf.segNvcReturn);
        epochTime=results.time(cut(:,1))-results.time(cut(1,1));
        layoutRoot=getNVCMetrics(epochTime,sCalc);

        results.nvc.time      =single(layoutRoot.time(:));
        results.nvc.stimulus  =single(layoutRoot.u(:));
        results.nvc.epochStart=single(epochStart(:));

        % ---- epoch quality, from the recording's own mean trace ----------------
        [qTrace,qWhat]=meanTrace(results);
        if isempty(qTrace)
            error('runNVC:noTraces', ...
                ['%s carries none of sData / dvsData / dvsDiameter, so there is ' ...
                 'nothing to cut.  Run the segmentation first.'],s.fName);
        end
        [reject,quality]=epochQuality(qTrace,results.time,cut,layoutRoot,sf);

        % ---- the picker, and then the page drawn fresh from the decisions ------
        if isequal(sf.enablelRejectionModification,1)
            reject=pickEpochs(results.time,qTrace,qWhat,epochStart, ...
                epochStart+layoutRoot.time(end),quality.timeLossSeries,reject);
        end
        drawCutsPage(rep,results.time,qTrace,qWhat,epochStart, ...
            epochStart+layoutRoot.time(end),quality.timeLossSeries,reject,quality);

        valid=~any(reject,1);
        if ~any(valid)
            error('runNVC:noEpochsAccepted', ...
                ['Every epoch of %s was rejected, so no aggregate can be taken.  ' ...
                 'Relax the reject* coefficients, or keep an epoch in the picker.'], ...
                s.fName);
        end
        results.nvc.epochValid  =valid(:);
        results.nvc.epochReject =reject;
        results.nvc.epochQuality=rmfield(quality,'timeLossSeries');

        % =============================================================
        % Per-signal analysis.  Each signal brings its OWN clock and therefore its own
        % cut index and its own layout: sData is 4 Hz and gsData is 100 Hz in the same
        % recording, which is why the filter length is a number of SECONDS.
        % =============================================================
        nwSeg=0; if sf.parforNvcSegments, nwSeg=Inf; end   %worker bound, not a branch
        for kSig=1:numel(sf.nvcSignals)
            sigName=sf.nvcSignals{kSig};
            if ~isfield(results,sigName) || isempty(results.(sigName)), continue; end
            prefix=nvcPrefix(sigName);

            [sigCut,L,LTree,LFit]=signalLayout(results,sigName,sf,sCalc,cut,layoutRoot);
            names=L.scalarNames;
            treeNames=LTree.scalarNames;
            fitNames=fitScalarNames(LFit,names);

            sigMat=results.(sigName);
            nSeg=size(sigMat,2);

            % ---- the epoch-average, and the polarity it resolves ---------------
            %Accumulated epoch by epoch rather than reshaped in one go: a 100 Hz signal
            %is [3000 x 20 x 1940], which is a gigabyte nobody needs to hold at once.
            avg=epochAverage(sigMat,sigCut,valid);
            pol=resolvePolarity(avg,L,sf);

            % ---- one parfor over SEGMENTS, epochs serial inside ----------------
            %Segments and not epochs: 1940 against 20, and the layout, the cut index
            %and both fit layouts are shared by every iteration.  BOTH sliced outputs
            %are assigned on EVERY iteration, so nothing about the parfor depends on
            %which branch the core took.
            nFit=1; if sf.nvcFit, nFit=numel(fitNames); end
            M=nan(nSeg,numel(names),nEp,'single');
            A=nan(nSeg,numel(fitNames),'single');
            F=nan(nSeg,nFit,nEp,'single');
            %the workbench hooks are function handles bound to a GUI; they are
            %transport rather than parameters, the core never reads them, and
            %broadcasting them to workers is the one thing in s that can fail to
            %serialise.  reportSettings strips exactly those two, here as at the save.
            sPar=reportSettings(sCalc);
            doFit=sf.nvcFit;
            parfor (i=1:nSeg, nwSeg)
                [M(i,:,:),A(i,:),F(i,:,:)]=nvcSegment(sigMat(:,i),avg(:,i),pol(i), ...
                    sigCut,L,LFit,sPar,names,fitNames,doFit,nFit);
            end

            % ---- assemble the sub-tree from the USER's level selection ---------
            T=struct();
            for k=1:numel(treeNames)
                T.(pfx(prefix,treeNames{k}))= ...
                    reshape(M(:,strcmp(names,treeNames{k}),:),nSeg,nEp);
            end
            if sf.nvcFit
                for k=1:numel(fitNames)
                    T.(pfx(prefix,fitNames{k}))=reshape(F(:,k,:),nSeg,nEp);
                end
            end
            results.nvc.esMetrics.(sigName)=T;

            G=struct();
            for k=1:numel(fitNames), G.(pfx(prefix,fitNames{k}))=A(:,k); end
            results.nvc.agg.(sigName)=G;

            % ---- the metrics tables: the MEDIAN over the valid epochs ----------
            results=duplicateToTable(results,sigName,prefix,M,names,valid,nSeg,nEp);
        end

        % =============================================================
        % Per-pixel markers (gated by a non-empty s.ppxNvcReturn; ns prefix, since
        % source.data is a flow cube).  Every pixel goes through the same core the
        % segments did, so a per-pixel map and a per-segment scalar cannot drift.
        % NEVER A FIT: 23.2 M fits is 266 hours.
        % =============================================================
        if ~isempty(sf.ppxNvcReturn)
            results.nvc.ppx=perPixel(results,source,cut,valid,sf);
        end

        % =============================================================
        % The two remaining pages: the response with its fit, and the trials.
        % =============================================================
        drawResponsePage(rep,results,qTrace,qWhat,cut,valid,sf);
        drawTrialsPage(rep,results,valid,sf);

        settings.runNVC=reportSettings(sf);
        reportWriting(rep);
        %NON-DESTRUCTIVE: SOURCE (_d) is never re-saved - only RESULTS and SETTINGS.
        save(s.fName,'results','-v7.3','-nocompression');
        save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);

end

% =====================================================================
function [rows,aggRow,fitRows]=nvcSegment(y,yAvg,sg,cut,L,LFit,s,names,fitNames,doFit,nFit)
%nvcSegment  One segment: every epoch measured, then the aggregate fit.
%   THE RESOLVED POLARITY IS SET HERE AND NOWHERE ELSE.  s arrives as the broadcast
%   copy every worker shares; adding the sign to a LOCAL copy is what keeps one sign
%   per segment without a per-segment broadcast.  The core falls back to resolving from
%   the single epoch in hand when the field is absent, which is a convenience for a
%   lone caller and exactly what a producer must not rely on.
s.nvcPolarityResolved=sg;
E=y(cut);                                   % [nES x nEp], the epochs of this segment
nEp=size(E,2);

rows=nan(numel(names),nEp,'single');
for k=1:nEp
    m=getNVCMetrics(E(:,k),L,s);
    for j=1:numel(names), rows(j,k)=m.(names{j}); end
end
rows=reshape(rows,1,numel(names),nEp);

% THE AGGREGATE FIT IS OF THE EPOCH-AVERAGE, ONCE.  Averaging the per-epoch fits would
% give a parameter set that is not a fit of anything; fitting the average gives one
% curve whose SNR is nine baseline SDs rather than four, which is the only level this
% kernel's two timescales are identifiable at.
aggRow=fitRow(yAvg,LFit,s,fitNames);

fitRows=nan(nFit,nEp,'single');
if doFit
    for k=1:nEp, fitRows(:,k)=fitRow(E(:,k),LFit,s,fitNames); end
end
fitRows=reshape(fitRows,1,nFit,nEp);
end

% =====================================================================
function row=fitRow(y,LFit,s,fitNames)
%fitRow  One trace through fitNVC, flattened to a fixed-width row.
%   Driven by the layout's own name list, so the core owns which scalars exist and
%   this file never types one.
row=nan(numel(fitNames),1,'single');
if isempty(fitNames) || isempty(LFit), return, end
m=fitNVC(y,LFit,s);
for k=1:numel(fitNames), row(k)=m.(fitNames{k}); end
end

% =====================================================================
function [sigCut,L,LTree,LFit]=signalLayout(results,sigName,sf,sCalc,cut,layoutRoot)
%signalLayout  This signal's cut index and its three layouts.
%   gsData is the per-region guided trace and is stored at the RAW frame rate on
%   results.gsTime - 100 Hz against the segment traces' 4 Hz on the reference
%   recording - so it is cut on its own clock and filtered by its own sample count.
%   Everything else shares the primary cut, which is already built.
%   L is what the core is CALLED with (the user's levels plus the table's), LTree is
%   the user's own selection and supplies the names that reach the saved tree, and
%   LFit serves both fits.
if strcmp(sigName,'gsData')
    if ~isfield(results,'gsTime') || isempty(results.gsTime)
        error('runNVC:noGsTime', ...
            ['results.gsData is present but results.gsTime is not, so the guided ' ...
             'trace has no clock to cut it on.  Run the guided step again.']);
    end
    sigCut=cutIndex(results.gsTime,sf);
    epochT=results.gsTime(sigCut(:,1))-results.gsTime(sigCut(1,1));
    L=getNVCMetrics(epochT,sCalc);
else
    sigCut=cut;
    L=layoutRoot;
end
LTree=getNVCMetrics(L.time,sf);
%BOTH FITS SHARE ONE LAYOUT, AND IT IS BUILT AT THE 'model' LEVEL ONLY - the other
%level is the fitted CURVE, which is [nES x 1] per trace and cannot go into a
%[nSeg x nEp] array of scalars.  The layout is this signal's own, so the fit core's
%marker twins are measured on the same clock, in the same windows and with the same
%hold rule as the per-epoch metrics they are twins of.
sFit=sf;
sFit.segNvcReturn={'model'};
LFit=fitNVC(L.time,sFit);
end

% =====================================================================
function names=fitScalarNames(LFit,metricNames)
%fitScalarNames  The fit's scalar contract, checked against the metrics' one.
%   The two name sets share a namespace inside esMetrics.<signal>, so a collision is
%   an overwrite nobody would see.  It cannot happen today - fitNVC's model-free
%   measurement IS getNVCMetrics and its twins carry a Fit suffix - and this is one
%   line to make sure it never starts.
names=LFit.scalarNames(:)';
clash=intersect(names,metricNames,'stable');
if ~isempty(clash)
    error('runNVC:scalarNameClash', ...
        ['fitNVC and getNVCMetrics both emit %s, and both are written into ' ...
         'results.nvc.esMetrics with the same prefix.  One of them has to be ' ...
         'renamed in its own core.'],strjoin(clash,', '));
end
end

% =====================================================================
function avg=epochAverage(sigMat,cut,valid)
%epochAverage  Each segment's average over the VALID epochs, [nES x nSeg].
%   Accumulated epoch by epoch: the reshape that would do it in one expression is
%   [nES x nEp x nSeg], which for a 100 Hz signal over 1940 segments is a gigabyte.
%   A segment that is non-finite anywhere stays non-finite here, which is what makes
%   its polarity and its aggregate fit NaN rather than an average of the good half.
k=find(valid);
avg=zeros(size(cut,1),size(sigMat,2));
for i=1:numel(k)
    avg=avg+sigMat(cut(:,k(i)),:);
end
avg=avg/numel(k);
end

% =====================================================================
function sg=resolvePolarity(avg,L,s)
%resolvePolarity  The sign every threshold of a segment is mirrored by, [1 x nSeg].
%   'auto' takes it from the segment's EPOCH-AVERAGE - the mean of the filtered
%   average minus its own baseline, over the response window - so one segment has one
%   direction.  Resolved per epoch instead, a segment that does not respond flips sign
%   on noise and its Peak becomes the maximum of |noise|.
nSeg=size(avg,2);
switch L.polarity
    case 'positive', sg=ones(1,nSeg);  return
    case 'negative', sg=-ones(1,nSeg); return
end
if isfield(s,'nvcPolarityResolved') && ~isempty(s.nvcPolarityResolved)
    sg=sign(double(s.nvcPolarityResolved(1)))*ones(1,nSeg); return
end
af=movmedian(avg,L.medN,1);
bl=mean(avg(L.blIdx,:),1);
sg=sign(mean(af(L.respIdx,:)-bl,1));
sg(sg==0|~isfinite(sg))=1;                  % a flat or invalid trace gets a direction
end

% =====================================================================
function [cut,epochStart]=cutIndex(time,s)
%cutIndex  Where every epoch sits in a signal's own samples, [nES x nEp].
%
%   THE NEAREST FRAME TO EACH EDGE, which is the rule the epoching step used and the
%   reason a dropped frame is a time LOSS rather than a time SHIFT: the epoch still
%   starts where the protocol says it starts, and the gap shows up in the quality
%   number instead of moving every later epoch along.
%
%   EVERY EPOCH IS THE SAME NUMBER OF SAMPLES, and it has to be - the epoch clock, the
%   filter length, the windows and the [nSeg x nEp] arrays are all shared.  The count
%   is the TYPICAL span, the median over the epochs, and not the shortest one: the
%   last epoch of a recording that stops on its own final sample is one frame short
%   by arithmetic rather than by any fault, and letting that shorten all twenty would
%   throw away a sample from every epoch to describe the edge of one.  An epoch that
%   really did lose a frame then covers slightly more wall-clock time in the same
%   number of samples, which is exactly what the timeLoss rule is there to say.
time=double(time(:));
nT=numel(time);
startSec=s.epochStartSec(:)';
endSec  =startSec+s.epochDurationSec;

%[nT x 1] against [1 x nEp] expands to [nT x nEp]; the nearest frame to each edge is
%then the min DOWN the columns.
[~,startFrame]=min(abs(time-startSec),[],1);
[~,endFrame]  =min(abs(time-endSec),[],1);

if endSec(end)>time(end)+median(diff(time))
    error('runNVC:protocolOutlastsRecording', ...
        ['The protocol asks for %d epochs of %g s starting %g s in, which ends at ' ...
         '%.1f s - but the recording is %.1f s long.  Either epochsN, ' ...
         'epochDurationSec or the stimulus start is wrong.'], ...
        numel(startSec),s.epochDurationSec,startSec(1),endSec(end),time(end));
end

nES=min(round(median(endFrame-startFrame)), nT-startFrame(end)+1);
if nES<4
    error('runNVC:shortEpoch', ...
        ['An epoch of %g s holds only %d samples of this recording, which is too ' ...
         'few to measure a response in.'],s.epochDurationSec,nES);
end
cut=startFrame+(0:nES-1)';                  % [nES x nEp], column k = epoch k
epochStart=time(startFrame);
end

% =====================================================================
function [reject,q]=epochQuality(y,time,cut,L,s)
%epochQuality  The five trace-only rules, and the numbers they fired on.
%
%   FIVE RULES AND NO IMAGES.  The two image-similarity rules this replaces existed to
%   keep a moved epoch out of an AVERAGE; nothing is averaged before it is measured any
%   more, so they went with the average rather than being ported onto the traces
%   (author, 2026-08-05).  What is left needs no cube, which is what makes this step
%   cost seconds instead of the 31.6 s and 21 GB of loading the flow cube.
%
%   EACH STORED NUMBER IS IN THE UNITS OF ITS OWN COEFFICIENT - four of them in
%   standard deviations of the spread ACROSS epochs, timeLoss in seconds - so the
%   report page can print the rule and the number it fired on side by side, and a
%   coefficient can be read straight off a rejected epoch.
y=double(y(:));
nEp=size(cut,2);

blLevel =mean(y(cut(L.blIdx ,:)),1);
finLevel=mean(y(cut(L.finIdx,:)),1);
epLevel =mean(y(cut),1);
pkLevel =max(y(cut),[],1);

%the sampling gaps: everything beyond one nominal step is time the camera did not
%deliver.  Kept as a SERIES as well as a sum, because it is what the cuts page draws.
step=median(diff(time));
gaps=[0; diff(double(time(:)))-step];
gaps(gaps<0)=0;

q.blDev   =single(spreadSD(blLevel));
q.finDev  =single(spreadSD(finLevel));
q.meanDev =single(spreadSD(epLevel));
q.peakDev =single(spreadSD(pkLevel));
q.timeLoss=single(sum(gaps(cut),1));
q.timeLossSeries=gaps;

%ONE ROW PER RULE, in the field order of q above - which is where that order is
%written down, and the only place it is.
reject=false(5,nEp);
reject(1,:)=double(q.blDev)  >s.rejectBlCoef;
reject(2,:)=double(q.finDev) >s.rejectFinCoef;
reject(3,:)=double(q.meanDev)>s.rejectEpochCoef;
reject(4,:)=double(q.peakDev)>s.rejectPeakCoef;
reject(5,:)=double(q.timeLoss)>=s.rejectTimeLoss;
if isequal(s.rejectFirstEpoch,1), reject(:,1)=true; end

for f={'blDev','finDev','meanDev','peakDev','timeLoss'}
    q.(f{1})=q.(f{1})(:);
end
end

% =====================================================================
function d=spreadSD(x)
%spreadSD  How far each epoch sits from the median of them all, in SDs of that spread.
%   A spread of zero is not a division to guard against afterwards - it means every
%   epoch agrees, so nothing deviates.
x=double(x(:))';
sd=std(x,'omitnan');
d=abs(x-median(x,'omitnan'));
if sd>0, d=d/sd; else, d=zeros(size(x)); end
end

% =====================================================================
function reject=pickEpochs(time,y,what,epochStart,epochEnd,gaps,reject)
%pickEpochs  The operator's accept/reject pass, on a VISIBLE figure.
%   The picker and the report page are two different objects on purpose: the picker is
%   maximised to whatever monitor is in front of the operator, the page is drawn fresh
%   from the FINAL decisions onto the canonical canvas.  Drawn as one figure, a JPEG's
%   pixel size became a function of the screen it was made on and the window was still
%   open when the next file started.
hInt=figure('Name','runNVC - accept or reject epochs','NumberTitle','off');
cleaner=onCleanup(@() delete(hInt(isgraphics(hInt))));   % closes it on every path
hInt.WindowState='maximized';
axInt=axes('Parent',hInt);
drawEpochs(axInt,time,y,what,epochStart,epochEnd,gaps,reject);
title(axInt,['Epochs (green - kept, red - rejected).  Click an epoch to change the ' ...
    'decision, then press Enter'])
[x,~]=ginput();
for i=1:numel(x)
    k=find(x(i)>=epochStart & x(i)<epochEnd,1);
    if ~isempty(k), reject(:,k)=~any(reject(:,k)); end
end
end

% =====================================================================
function drawCutsPage(rep,time,y,what,epochStart,epochEnd,gaps,reject,q)
%drawCutsPage  '_rep_nvccuts.jpg' - where every epoch was cut, and what was kept.
%   The successor to the epochs page of the step this replaces, drawn on the whole
%   recording rather than on a cube that no longer exists.
fh=[];
try
    fh=reportFigure(rep,'nvccuts','single');
    tl=tiledlayout(fh,1,1,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl);
    drawEpochs(ax,time,y,what,epochStart,epochEnd,gaps,reject);
    kept=sum(~any(reject,1));
    title(ax,sprintf('%d of %d epochs kept  (worst deviation %.1f SD, worst loss %.2f s)', ...
        kept,size(reject,2),max([0;double(q.meanDev(:))]),max([0;double(q.timeLoss(:))])))
    reportSave(rep,fh,'nvccuts');   % deletes the figure on every path of its own
    fh=[];
catch
    % a page is a by-product and must never kill a run - reportSave says the same
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function drawEpochs(ax,time,y,what,epochStart,epochEnd,gaps,reject)
%drawEpochs  The trace, the cuts, one accept/reject bar per rule, and the time loss.
%   ONE function so the picker and the report page cannot drift apart - drawn twice
%   they differed by the colour of the epoch lines and nobody noticed for a year.
nEp=numel(epochStart); nRule=size(reject,1);
lo=min(y); hi=max(y);
if ~(isfinite(lo)&&isfinite(hi)&&hi>lo), lo=0; hi=1; end
yyaxis(ax,'left')
plot(ax,time,y,'-','Color',[0.25 0.35 0.55])
hold(ax,'on')
for i=1:nEp
    plot(ax,[epochStart(i) epochStart(i)],[lo hi],'--k')
    for ii=1:nRule
        yb=lo-0.015*ii*(hi-lo);
        if reject(ii,i), c='-r'; else, c='-g'; end
        plot(ax,[epochStart(i) epochEnd(i)],[yb yb],c,'LineWidth',1.5)
    end
end
plot(ax,[epochEnd(end) epochEnd(end)],[lo hi],'--r')
ylabel(ax,what)
hold(ax,'off')
%the bars hang BELOW the trace, so the left limits are set from both of them together
%rather than left to a tight fit that would clip whichever was drawn last
pad=0.015*(nRule+1)*(hi-lo);
ylim(ax,[lo-pad hi+0.04*(hi-lo)])
yyaxis(ax,'right')
plot(ax,time,gaps,'-')
ylabel(ax,'Time loss, s')
xlabel(ax,'Time, s')
xlim(ax,[time(1) time(end)])
end

% =====================================================================
function drawResponsePage(rep,results,y,what,cut,valid,s)
%drawResponsePage  '_rep_nvcresponse.jpg' - the epoch-average, its fit, the residuals
%   and the parameters.  The successor to the epoch-average and the NVC-fit pages in
%   one, which is what they always were: two views of the same curve.
%
%   THE TRACE IS DRAWN RAW.  The median filter this step carries is for the crossing
%   rules alone, and a page that smoothed the curve would be showing a peak the
%   recording does not contain at a time it did not happen.
%
%   AND THE PAGE CANNOT KILL THE RUN.  reportSave already swallows a failed export;
%   this swallows a failed DRAWING, because the results are computed by the time the
%   page is attempted and losing a recording's analysis to a plotting fault would be
%   absurd.
fh=[];
try
    avg=mean(y(cut(:,valid)),2);
    t=double(results.nvc.time);
    %the page needs the CURVE as well as the parameters, so it asks for the
    %reconstruction the analysis levels deliberately do not
    sPage=s; sPage.segNvcReturn={'model','reconstruction'};
    L=fitNVC(t,sPage);
    m=fitNVC(avg,L,sPage);
    fd=double(m.fData);
    if isempty(fd)||numel(fd)~=numel(t), fd=nan(size(t)); end

    fh=reportFigure(rep,'nvcresponse');
    tl=tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');

    % ---- the average response and the fit ----
    ax=nexttile(tl,1,[1 2]);
    hold(ax,'on')
    xregion(ax,L.tStim,L.tStim+L.D,'FaceColor',[0.30 0.55 0.85],'FaceAlpha',0.12);
    plot(ax,t,avg,'-','Color',[0.45 0.45 0.45],'LineWidth',1)
    plot(ax,t,fd,'-','Color',[0.85 0.20 0.15],'LineWidth',1.5)
    hold(ax,'off')
    xlabel(ax,'Time in the epoch, s'); ylabel(ax,what)
    legend(ax,{'stimulus','measured','fitted'},'Location','best','Box','off')
    axis(ax,'tight'); grid(ax,'on')
    title(ax,sprintf('Response averaged over %d of %d epochs',sum(valid),numel(valid)))

    % ---- the residuals ----
    ax=nexttile(tl,4,[1 2]);
    plot(ax,t,avg-fd,'-','Color',[0.20 0.35 0.60])
    yline(ax,0,'-','Color',[0.6 0.6 0.6]);
    xlabel(ax,'Time in the epoch, s'); ylabel(ax,'Measured - fitted')
    axis(ax,'tight'); grid(ax,'on')
    title(ax,'Residuals')

    % ---- the parameters ----
    ax=nexttile(tl,3,[2 1]);
    axis(ax,'off')
    text(ax,0.02,0.98,pageLines(m,L),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName','Courier New','FontSize',10,'Interpreter','none')
    title(ax,'Fitted response')

    reportSave(rep,fh,'nvcresponse');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function txt=pageLines(m,L)
%pageLines  The parameter block, written for a reader rather than for a parser.
%   Every line asks the bag whether it has a NUMBER for the field, because which
%   scalars carry one is the core's decision and it differs by model AND by fit -
%   'doubleGamma' has no damping ratio at all, and the second-order model's damping
%   block is gated by Identified, so it is printed for one segment and not the next.
lines={sprintf('model        %s',L.model)};
lines=[lines,num(m,'Gain',      'gain         %.4g')];
lines=[lines,num(m,'Onset',     'onset        %.3g s')];
lines=[lines,num(m,'TauS',      'tau signal   %.3g s')];
lines=[lines,num(m,'TauF',      'tau feedback %.3g s')];
lines=[lines,num(m,'Zeta',      'damping      %.3f')];
lines=[lines,num(m,'RingPeriod','ring period  %.3g s')];
lines=[lines,num(m,'A1',        'peak shape   %.3g')];
lines=[lines,num(m,'Beta1',     'peak rate    %.3g')];
lines=[lines,num(m,'CRatio',    'dip fraction %.3f')];
lines=[lines,{''}];
lines=[lines,num(m,'PeakFit',   'peak         %.4g')];
lines=[lines,num(m,'TTPFit',    'time to peak %.3g s')];
lines=[lines,num(m,'DurFit',    'duration     %.3g s')];
lines=[lines,{''}];
lines=[lines,num(m,'R2',        'R2           %.4f')];
lines=[lines,num(m,'RMSE',      'RMSE         %.4g')];
lines=[lines,num(m,'StartsAgree',sprintf('starts agree %%d of %d',L.nStarts))];
% IDENTIFIED IS THE ONE NUMBER THAT CHANGES WHAT THE REST MEAN, so it is called out
% rather than filed in the block.  On real data a shape parameter sat on its bound in
% 51.5 % of segments, and a railed parameter reported as a number is the failure this
% design exists to prevent.  The damping lines above are already absent when it fires -
% the core emits no number for them - and this says why the block got shorter.
if isfield(m,'Identified') && isfinite(m.Identified) && m.Identified<=0
    lines=[lines,{'','*** NOT IDENTIFIED ***', ...
        'a shape parameter sat on its', ...
        'bound, so the time constants', ...
        'above describe the bound and', ...
        'not the vessel'}];
end
txt=strjoin(lines,newline);
end

function c=num(m,name,fmt)
%num  One parameter line, or nothing at all when this model does not emit it.
c={};
if isfield(m,name) && isfinite(m.(name)), c={sprintf(fmt,double(m.(name)))}; end
end

% =====================================================================
function drawTrialsPage(rep,results,valid,s)
%drawTrialsPage  '_rep_nvctrials.jpg' - the response epoch by epoch across the session.
%
%   THIS IS THE PAGE THE REDESIGN IS FOR.  Averaging twenty repetitions before
%   measuring them made habituation, drift and a single bad repetition all look
%   identical - like nothing.  Here each epoch is a point: the median across segments
%   with the interquartile band around it, so a response that fades over the session is
%   a slope a reader sees rather than a number nobody computed.
%
%   IT PLOTS THE AMPLITUDES AND NOT THE TIMINGS, and that is measured rather than
%   aesthetic: 0.61 of Peak's variance belongs to the segment against 0.13 of TTP's, so
%   a per-epoch timing plot would be a plot of trial noise.
fh=[];
try
    [sig,mets]=trialsSource(results,s);
    if isempty(sig), return, end
    fh=reportFigure(rep,'nvctrials');
    tl=tiledlayout(fh,2,1,'TileSpacing','compact','Padding','compact');
    ep=1:numel(valid);
    lbl={'Peak, fraction of baseline','Flow through the response / flow expected'};
    for k=1:2
        v=double(results.nvc.esMetrics.(sig).(mets{k}));
        p=prctile(v,[25 50 75],1);
        ax=nexttile(tl,k);
        hold(ax,'on')
        fill(ax,[ep flip(ep)],[p(1,:) flip(p(3,:))],[0.30 0.55 0.85], ...
            'FaceAlpha',0.18,'EdgeColor','none');
        plot(ax,ep,p(2,:),'-o','Color',[0.20 0.35 0.60],'MarkerFaceColor','w', ...
            'LineWidth',1.2,'MarkerSize',4)
        %the legend names what was DRAWN: a session in which every epoch was kept has
        %no rejected marker, and a third entry for it would be a key to nothing
        names={'middle half of the segments','median segment'};
        if any(~valid)
            plot(ax,ep(~valid),p(2,~valid),'x','Color',[0.85 0.20 0.15], ...
                'LineWidth',1.5,'MarkerSize',10)
            names=[names,{'rejected'}];
        end
        hold(ax,'off')
        ylabel(ax,lbl{k}); grid(ax,'on')
        xlim(ax,[0.5 numel(ep)+0.5])
        %the band is the point of the panel, so it gets room rather than being fitted
        %flush to the box where its edges read as the axes rather than as the data
        b=[min(p(1,:)) max(p(3,:))];
        if all(isfinite(b)) && b(2)>b(1), ylim(ax,b+0.06*diff(b)*[-1 1]); end
        if k==1
            legend(ax,names,'Location','best','Box','off')
            title(ax,'Response epoch by epoch')
        end
    end
    xlabel(ax,'Epoch')
    reportSave(rep,fh,'nvctrials');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function [sig,mets]=trialsSource(results,s)
%trialsSource  The first analysed signal that carries both amplitude metrics.
%   A run that stored neither the marker nor the area level has nothing to draw here,
%   and an empty page is better than an invented one.
sig=''; mets={};
if ~isfield(results,'nvc') || ~isfield(results.nvc,'esMetrics'), return, end
for c=s.nvcSignals(:)'
    nm=c{1};
    if ~isfield(results.nvc.esMetrics,nm), continue, end
    p=nvcPrefix(nm);
    want={pfx(p,'PeakRel'),pfx(p,'AUCblkRel')};
    if all(isfield(results.nvc.esMetrics.(nm),want)), sig=nm; mets=want; return, end
end
end

% =====================================================================
function ppx=perPixel(results,source,cut,valid,s)
%perPixel  Per-epoch marker maps, [Y x X x nEp] each.  THE ONE PATH THAT READS SOURCE.
%   Every pixel goes through the core the segments did, resolving its own polarity from
%   its own epoch-average exactly as a segment does, so a map and a scalar cannot mean
%   two different things.  A background pixel is written by no iteration and stays NaN:
%   an unmeasured pixel has no response, and a zero there would be a measurement that
%   was never made.
if ~isfield(results,'sMap')
    error('runNVC:noSMap', ...
        'results.sMap not found; runSegmentation must be run before per-pixel NVC.');
end
sPix=s; sPix.segNvcReturn=s.ppxNvcReturn;
L=getNVCMetrics(results.nvc.time,sPix);
names=L.scalarNames;
nEp=size(cut,2);

sz=size(source.data); Y=sz(1); X=sz(2); npx=Y*X;
D=reshape(source.data,npx,sz(3));
sMapLin=results.sMap(:);

Mp=nan(npx,numel(names),nEp,'single');
nwPix=0; if s.parforNvcPixels, nwPix=Inf; end   %worker bound, not a branch
sPar=reportSettings(sPix);
parfor (p=1:npx, nwPix)
    if sMapLin(p)==0, continue; end                  %background stays NaN
    Mp(p,:,:)=ppxPixel(double(D(p,:)).',cut,valid,L,sPar,names);
end

ppx=struct();
for k=1:numel(names)
    ppx.(pfx('ns',names{k}))=reshape(Mp(:,k,:),Y,X,nEp);
end
end

% =====================================================================
function rows=ppxPixel(y,cut,valid,L,s,names)
%ppxPixel  One pixel: its own polarity from its own epoch-average, then every epoch.
E=y(cut);
s.nvcPolarityResolved=resolvePolarity(mean(E(:,valid),2),L,s);
nEp=size(E,2);
rows=nan(numel(names),nEp,'single');
for k=1:nEp
    m=getNVCMetrics(E(:,k),L,s);
    for j=1:numel(names), rows(j,k)=m.(names{j}); end
end
rows=reshape(rows,1,numel(names),nEp);
end

% =====================================================================
function results=duplicateToTable(results,sigName,prefix,M,names,valid,nSeg,nEp)
%duplicateToTable  The epoch-aggregated key set into sMetrics / dvsMetrics.
%   THE AGGREGATE IS A MEDIAN OVER THE VALID EPOCHS, never a mean - see the header on
%   the measured bias in a per-epoch maximum.  Row order matches the 1:nSeg loop, the
%   same segment invariant runPulsatility relies on.  BFI / std(BFI) are runBFI's
%   columns and are intentionally not touched.
switch sigName
    case 'sData',                    tbl='sMetrics';
    case {'dvsData','dvsDiameter'},  tbl='dvsMetrics';
    otherwise, return                % gsData has no table, exactly as in runVasomotion
end
if ~isfield(results,tbl), return, end

dup={'PeakRel','TTP','AUCf','AUCblkRel','SNR'};
for k=1:numel(dup)
    j=find(strcmp(names,dup{k}),1);
    if isempty(j), continue, end
    v=reshape(M(:,j,:),nSeg,nEp);
    v(:,~valid)=NaN;
    results.(tbl).([prefix dup{k}])=single(median(v,2,'omitnan'));
end
%NEpochs counts, per segment, the valid epochs the medians above were taken over.  It
%is read off Peak because the core makes an unmeasurable epoch NaN in EVERY field at
%once, so one column answers it for all of them - and without it a median of one epoch
%and a median of twenty are the same number.
j=find(strcmp(names,'Peak'),1);
if ~isempty(j)
    v=reshape(M(:,j,:),nSeg,nEp);
    results.(tbl).([prefix 'NEpochs'])=single(sum(isfinite(v(:,valid)),2));
end
end

% =====================================================================
function [y,what]=meanTrace(results)
%meanTrace  The recording's average trace, and what it is a trace of.
%   The DIAMETER is last on purpose: it is a length in pixels rather than a flow index,
%   so it is what the page says it is only when there is no flow trace at all.
y=[]; what='';
for c={'sData','Flow, BFI'; 'dvsData','Flow, BFI'; 'dvsDiameter','Diameter, px'}'
    nm=c{1}; lbl=c{2};
    if isfield(results,nm) && ~isempty(results.(nm))
        v=mean(double(results.(nm)),2,'omitnan');
        if all(isfinite(v)), y=v; what=lbl; return, end
    end
end
end

% =====================================================================
function startSec=epochStarts(s,timeStamp,fidx)
%epochStarts  Where every epoch begins, in seconds from the recording start.
%   The recording carries its absolute start time; the protocol says when the first
%   epoch began, either as an offset from that or as a wall-clock time per file.  The
%   rest tile forward at the epoch duration.
rlsStart=datetime(timeStamp,'ConvertFrom','epochtime', ...
    'Epoch',datetime(1970,1,1),'TicksPerSecond',1e3,'Format','HH:mm:ss.SSS');
switch s.stimStartType
    case 'manual'
        list=s.stimStartAll;
        if ~iscell(list), list={list}; end
        if isscalar(list), item=list{1}; else, item=list{fidx}; end
        stimStart=datetime(item,'InputFormat','HH:mm:ss.SSS');
    case 'offset'
        stimStart=rlsStart+seconds(s.stimOffset);
end
first=((hour(stimStart)-hour(rlsStart))*60 + ...
       (minute(stimStart)-minute(rlsStart)))*60 + ...
       (second(stimStart)-second(rlsStart));
if first<0
    error('runNVC:stimBeforeRecording', ...
        ['The first epoch is timed %.1f s BEFORE the recording started.  Check ' ...
         's.stimStartAll against this recording, or use the offset mode.'],-first);
end
startSec=first+(0:s.epochsN-1)*s.epochDurationSec;
end

% =====================================================================
function s=nvcDefaults(s)
%nvcDefaults  Resolve every default HERE, so the saved settings record what ran.
%   segNvcReturn defaults when absent OR empty; ppxNvcReturn defaults ONLY when the
%   field is ABSENT - an explicit [] must stay empty (per-pixel off), exactly like
%   runPulsatility's s.ppxPulsReturn.
if ~isfield(s,'nvcMedFiltSec')||isempty(s.nvcMedFiltSec), s.nvcMedFiltSec=1.25;   end
if ~isfield(s,'nvcPolarity')  ||isempty(s.nvcPolarity),   s.nvcPolarity='auto';   end
if ~isfield(s,'nvcHoldSamples')||isempty(s.nvcHoldSamples),s.nvcHoldSamples=2;    end
if ~isfield(s,'nvcDecayFrac') ||isempty(s.nvcDecayFrac),  s.nvcDecayFrac=0.1;     end
if ~isfield(s,'nvcModel')     ||isempty(s.nvcModel),      s.nvcModel='secondOrder'; end
if ~isfield(s,'nvcStarts')    ||isempty(s.nvcStarts),     s.nvcStarts=16;         end
if ~isfield(s,'nvcFit')       ||isempty(s.nvcFit),        s.nvcFit=false;         end
if ~isfield(s,'nvcSignals')   ||isempty(s.nvcSignals)
    s.nvcSignals={'sData','gsData','dvsData','dvsDiameter'};
end
if ~isfield(s,'segNvcReturn') ||isempty(s.segNvcReturn)
    s.segNvcReturn={'markers','timing','areas'};
end
if ~isfield(s,'ppxNvcReturn'), s.ppxNvcReturn=[]; end
%the rejection rules.  Inf disables one, which is why the defaults are a number rather
%than a switch: a rule is tuned, not turned off.
if ~isfield(s,'rejectBlCoef')   ||isempty(s.rejectBlCoef),   s.rejectBlCoef=1;    end
if ~isfield(s,'rejectFinCoef')  ||isempty(s.rejectFinCoef),  s.rejectFinCoef=1;   end
if ~isfield(s,'rejectEpochCoef')||isempty(s.rejectEpochCoef),s.rejectEpochCoef=1; end
if ~isfield(s,'rejectPeakCoef') ||isempty(s.rejectPeakCoef), s.rejectPeakCoef=1;  end
if ~isfield(s,'rejectTimeLoss') ||isempty(s.rejectTimeLoss), s.rejectTimeLoss=0.5;end
if ~isfield(s,'rejectFirstEpoch')||isempty(s.rejectFirstEpoch),s.rejectFirstEpoch=1; end
%OFF BY DEFAULT HERE, ON BY DEFAULT IN A PROTOCOL.  A wrapper that blocks on a figure
%because a field was omitted cannot be run from a script; the launcher and the
%workbench preset both switch it on, which is where a person is choosing.
if ~isfield(s,'enablelRejectionModification')||isempty(s.enablelRejectionModification)
    s.enablelRejectionModification=0;
end
%PARALLELISM IS OPTIONAL.  Each switch is a BOUND on its parfor (Inf workers or 0),
%never a branch: parfor(...,0) runs the identical loop body serially IN THE CLIENT and
%starts no pool at all.
if ~isfield(s,'parforNvcSegments')||isempty(s.parforNvcSegments)
    s.parforNvcSegments=true;
end
if ~isfield(s,'parforNvcPixels')||isempty(s.parforNvcPixels)
    s.parforNvcPixels=true;
end

allSignals={'sData','gsData','dvsData','dvsDiameter'};
if ischar(s.nvcSignals)||isstring(s.nvcSignals), s.nvcSignals=cellstr(s.nvcSignals); end
if ~iscellstr(s.nvcSignals)||~all(ismember(s.nvcSignals,allSignals))
    error('runNVC:nvcSignals', ...
        's.nvcSignals must be a cell subset of {%s}.', ...
        strjoin(cellfun(@(x) ['''' x ''''],allSignals,'UniformOutput',false),', '));
end
end

% =====================================================================
function requireGeometry(s,nFiles)
%requireGeometry  The protocol, checked before any file is opened.
%   NOTHING UPSTREAM WRITES THIS GEOMETRY ANY MORE.  The epoching step that used to
%   record it for a later fit to inherit does not exist; the geometry is this step's
%   own input, there is no second place to read it from, and a missing or mis-sized
%   entry has to be a named message rather than a subscript error inside the loop.
need={'epochsN','epochDurationSec','epochStimStartSec','stimDurationSec'};
for k=1:numel(need)
    if ~isfield(s,need{k})||isempty(s.(need{k}))||~isscalar(s.(need{k}))|| ...
            ~isfinite(s.(need{k}))||s.(need{k})<=0
        error('runNVC:missingGeometry', ...
            's.%s is required - a positive number of %s.',need{k}, ...
            ternary(strcmp(need{k},'epochsN'),'epochs','seconds'));
    end
end
for k={'epochBaselineSec','epochFinaleSec'}
    if ~isfield(s,k{1})||numel(s.(k{1}))~=2||any(~isfinite(s.(k{1})))
        error('runNVC:missingWindow', ...
            ['s.%s is required - the [start end] of the window, in seconds.  The ' ...
             'baseline is measured from the epoch start and the finale back from ' ...
             'its end, so a typical pair is [0 10] and [-5 0].'],k{1});
    end
end
%A STIMULUS THAT OUTLASTS ITS EPOCH means every window below is wrong, and hearing it
%now costs a second instead of sixty files.
if s.epochStimStartSec+s.stimDurationSec>s.epochDurationSec
    error('runNVC:stimOutlastsEpoch', ...
        ['The stimulus starts %g s into the epoch and lasts %g s, which runs past ' ...
         'the end of the %g s epoch.  Either the stimulus duration or the epoch ' ...
         'duration is wrong.'],s.epochStimStartSec,s.stimDurationSec,s.epochDurationSec);
end
if ~isfield(s,'stimStartType')||~any(strcmp(s.stimStartType,{'offset','manual'}))
    error('runNVC:stimStartType', ...
        's.stimStartType must be ''offset'' or ''manual''.');
end
switch s.stimStartType
    case 'offset'
        if ~isfield(s,'stimOffset')||~isscalar(s.stimOffset)||~isfinite(s.stimOffset)
            error('runNVC:noStimOffset', ...
                ['s.stimOffset is required in offset mode - how many seconds after ' ...
                 'the recording starts the first epoch begins.']);
        end
    case 'manual'
        if ~isfield(s,'stimStartAll')||isempty(s.stimStartAll)
            error('runNVC:noStimStartAll', ...
                ['s.stimStartAll is required in manual mode - one ''HH:mm:ss.SSS'' ' ...
                 'time per file, in the order fNames lists them.']);
        end
        list=s.stimStartAll; if ~iscell(list), list={list}; end
        if ~(isscalar(list)||numel(list)==nFiles)
            error('runNVC:stimStartCount', ...
                ['s.stimStartAll has %d entries but there are %d files.  Give one ' ...
                 'time for all of them, or exactly one per file in the order fNames ' ...
                 'lists them.'],numel(list),nFiles);
        end
end
end

% =====================================================================
function lv=withTableLevels(lv)
%withTableLevels  The level set the core is CALLED with: the user's plus the two the
%   metrics tables draw from.  Neither costs a fit, and a table column that came and
%   went with a display setting would be worse than the duplication.
if ischar(lv)||isstring(lv), lv=cellstr(lv); end
if isempty(lv), lv={}; end
lv=reshape(lv,1,[]);
for k={'markers','areas'}
    if ~ismember(k{1},lv), lv=[lv k]; end %#ok<AGROW>
end
end

% =====================================================================
function p=nvcPrefix(sigName)
%nvcPrefix  The branch prefix of a signal: ns flow, nd diameter, ng guided contrast.
switch sigName
    case 'dvsDiameter', p='nd';
    case 'gsData',      p='ng';
    otherwise,          p='ns';
end
end

function n=pfx(p,name)
%pfx  The prefixed field name.  Every scalar starts with a capital, so this is a plain
%   concatenation and there is exactly one spelling of each quantity in the tree.
n=[p name];
end

function v=ternary(c,a,b)
if c, v=a; else, v=b; end
end
