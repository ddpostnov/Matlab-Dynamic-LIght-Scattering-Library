%runNVC  Model-free neurovascular response, every stimulus repetition on its own
%
%   runNVC(s,fNames) loads every contrast-branch *_t_BFI_r.mat (or *_s_BFI_r.mat)
%   whole-recording product in fNames (runContrastFromRLS -> runSegmentation ->
%   runBFI: results.sData/gsData/dvsData/dvsDiameter on results.time, plus the flow
%   cube source.data [Y x X x nT]), cuts every stimulus repetition out of the
%   recording, and MEASURES each one on its own.
%
%   THERE IS NO FIT.  The author's ruling of 06-Aug-2026 - "at this stage I do not see
%   value in fitting the data" - retired fitNVC, and with it the aggregate fit, the
%   per-epoch fit and every setting that served them.  The step measures the trace and
%   reports what it measured: fourteen markers per (segment, epoch) from
%   getNVCMetrics, two confidence numbers per (segment, epoch) from getNVCConfidence,
%   and one recording-wide trusted set from getNVCEpochTrust.
%
%   NOTHING IS AVERAGED BEFORE IT IS MEASURED, AND THERE IS NO EXCEPTION LEFT.  The
%   step this replaces averaged each segment's epochs twice - once to resolve a
%   polarity and once to fit - and both are gone: polarity is now sign(Auc) of the
%   epoch in hand.  The ONLY average in this file is the optional representative epoch,
%   it happens AFTER every marker and both confidence numbers exist for every
%   (segment, epoch), and the epochs it averages are chosen by them.  Author's words:
%   "Even in that case all markers are first calculated for all epochs and all segments
%   individually and then based on that the epochs to average are defined."
%
%   THE STEP RUNS STRICTLY AFTER THE BFI CONVERSION, and that is load-bearing rather
%   than merely an ordering.  Every trace it reads is already a flow index, so NO
%   1/K^2 APPEARS ANYWHERE IN THIS FILE OR IN ITS CORES, and nothing here asks what
%   format a trace is in - the '_BFI_' token in the product name is the answer.  A
%   '_K_' product is refused by name even though the registry cannot offer one, since
%   a launcher file is written by hand.
%
%   IT NEVER READS THE IMAGE CUBE FOR THE PER-SEGMENT PATH.  Only the per-pixel path
%   (s.ppxNvcReturn) and the representative-epoch collapse (s.nvcRepresentative) open
%   the _d member, and both are off by default, so the usual run never pays the 21 GB.
%
%   THE RUN, PER FILE
%     guard the input, and REFUSE a collapsed product by name
%     resolve the epoch starts and cut the index                        [nES x nEp]
%     PASS 1 - MEASURE.  per signal: the [nSeg x nMetric x nEp] marker block in a
%              parfor over segments, then ONE getNVCConfidence call on the whole block
%     the responding-area rule on the segmented flow signal -> epochTrust, areaFrac
%     the epoch editor, when s.nvcEditEpochs - the operator's set replaces the rule's
%     PASS 2 - AGGREGATE.  the metrics tables: the MEAN over the GLOBALLY-trusted set
%     four report pages
%     optionally the representative epoch, IN PLACE and destructive
%
%   THE TWO PASSES ARE NOT A STYLE CHOICE.  The tables average over the GLOBALLY
%   trusted epochs; that set does not exist until every signal's confidence has been
%   computed and the area rule has run, and the operator may still change it in the
%   editor.  So the aggregation cannot live inside the measuring loop.  Holding the
%   marker blocks costs 2.2 MB for 1940 segments, 15 metrics and 20 epochs, which is
%   nothing next to measuring twice.
%
%   THE CONFIDENCE IS COMPUTED ONCE PER SIGNAL, OUTSIDE THE parfor, and that is not a
%   preference either.  Every optional deviation factor is a deviation along the EPOCH
%   axis of one unit, so a worker holding one segment does not carry enough of the
%   problem to compute one.  Getting this wrong is not a correctness bug a test would
%   catch - it is a call that silently cannot be made.
%
%   TWO LAYOUTS, ONE TIME BASE, AND THEY ARE NOT OPTIONAL.  getNVCConfidence needs
%   Bl BlStd Fn FnStd Ep St AucRel and every T<p>, so a run with s.segNvcReturn
%   {'times'} would have no confidence at all.  SETUP is therefore called TWICE: once
%   with the COMPLETE marker set, which is what the block is measured with, and once
%   with the user's own levels, which is what NAMES the tree columns.  The producer
%   types no scalar name of its own - layout.blockNames is what M is built from,
%   layout.scalarNames is what is written, and they differ by exactly one name
%   (the lowercase tPeak, which confidence factor fPeak grades and decision 5 forbids
%   storing).
%
%   OUTPUT TREE  (RESULTS.nvc; every float SINGLE)
%       .time        [nES x 1]  the EPOCH clock, seconds - what every stored time is
%                               quoted against, so a downstream plot is self-describing
%       .stimulus    [nES x 1]  the boxcar on that clock
%       .epochStart  [nEp x 1]  where each epoch was cut, on the RECORDING clock, s
%       .epochTrust  [nEp x 1]  logical, the recording-wide trusted set
%       .areaFrac    [nEp x 1]  the responding-area share each decision was made on
%       .timeLoss    [nEp x 1]  seconds of dropped frames per epoch
%       .esMetrics   <- THE PRODUCT.  One sub-struct per analysed signal, each field
%                       [nSeg x nEp], row order = segment order (aligned to sMetrics /
%                       dvsMetrics rows), column order = epoch order:
%           .sData       .<ns*>     .gsData      .<ng*>
%           .dvsData     .<ns*>     .dvsDiameter .<nd*>
%                       The fourteen markers, plus Conf, ConfMin and Trust.
%       .ppx         (only when s.ppxNvcReturn is non-empty) [Y x X x nEp] per
%                       requested marker, plus nsConf and nsConfMin
%     esMetrics IS SPLIT BY SIGNAL because the segment sets differ: sData has 1940 rows
%     on the reference recording and dvsDiameter has 67, and they cannot share one
%     array.  The PREFIXES keep the existing vocabulary - ns = NVC flow, nd = NVC
%     diameter, ng = NVC guided contrast - and an nd* column must never be pooled with
%     an ns* one: dvsDiameter is a length in pixels, not a flow index.
%
%     epochValid, epochReject and epochQuality are GONE.  Their content moved: the five
%     across-epoch deviations became per-segment confidence factors (optional, off by
%     default), and "which rule fired" is answered by ConfMin and the confidence page.
%
%   THE EPOCH-CUT TRACES THEMSELVES ARE NOT STORED.  [nES x nSeg x nEp] is 18.6 GB on
%   the reference recording; the curves a reader wants are on the response page and the
%   marker viewer re-cuts them from results.time against .epochStart.
%
%   METRICS-TABLE DUPLICATION (a small key set, same names as the tree, single)
%     results.sMetrics    (from sData)      : nsStRel nsPkRel nsAucRel nsT50 nsConf
%     results.dvsMetrics  (from dvsData)    : the same five ns*
%                         (from dvsDiameter): the five nd*
%     gsData writes NO table - there is no gMetrics, exactly as in runVasomotion.
%     BFI / std(BFI) are runBFI's columns and are LEFT UNTOUCHED.
%
%     EVERY COLUMN IS THE MEAN OVER THE GLOBALLY-TRUSTED EPOCHS, the same epochs for
%     every row - the author's ruling of 06-Aug-2026.  That is why there is NO
%     nsNEpochs column: a count identical down 1940 rows is a recording-level number
%     and it lives in results.nvc.epochTrust and on the cuts page.  nsConf takes its
%     place as the per-segment trust number - that segment's MEAN confidence over the
%     accepted epochs - so a row whose segment was never trusted still exists and says
%     so, which is exactly the information a reader needs to discard it.
%
%     A MEAN, NOT A MEDIAN, and it costs something worth writing down.  The per-epoch
%     Pk is a positively biased estimator - a maximum over ~80 noisy samples sits about
%     2.4 baseline SDs above baseline on an epoch with no stimulus in it at all - so
%     nsPkRel carries that bias where the old median blunted it.  What changed is that
%     untrustworthy epochs are now excluded BY NAME instead of being averaged over
%     robustly, which is the better instrument and a worse excuse for not saying so.
%
%     nsT50 IS NAMED FROM THE SETTING, like every other time.  If s.nvcAreaPcts does
%     not contain 50 the table takes the requested level nearest to it and the column
%     is named after THAT level, so a column called nsT50 always holds T50.
%
%   THE LEVELS GATE THE TREE, NOT THE TABLE.  s.segNvcReturn decides only what reaches
%   the saved tree; the table's key set spans two levels, both cost nothing, and a
%   table column must not appear and disappear with a display setting.  Same rule as
%   runPulsatility and runFitVasoreactivity.
%
%   EACH SIGNAL BRINGS ITS OWN CLOCK.  gsData is the per-region guided trace stored at
%   the raw frame rate on results.gsTime - 100 Hz against the segment traces' 4 Hz on
%   the reference recording - so it is cut on its own clock, gets its own layout and
%   its own dropped-frame series.
%
%   THE EPOCH-LEVEL TRUST COMES FROM ONE SIGNAL, AND IT IS THE SEGMENTED FLOW TRACE.
%   sData, else dvsData, else dvsDiameter: results.sMap is the only spatial map the
%   product carries, so the responding-area rule is a sData rule.  The cuts page says
%   which signal was used.  An epoch is judged by AREA and never by an average over its
%   segments - the author's answer to "only a part of the segments are responders".
%
%   PARALLELISM IS A BOUND, NEVER A BRANCH.  parfor (i=1:n, 0) runs the identical loop
%   body serially in the client and starts no pool at all, so s.parforNvcSegments and
%   s.parforNvcPixels change what it costs and never what it computes.  reportSettings
%   strips the two workbench hooks before the broadcast - they are the one thing in s
%   that can fail to serialise.
%
%   INPUTS
%     s        parameter struct with fields
%              % where the epochs are
%                . stimStartType   'offset' | 'manual'
%                . stimOffset      seconds from the recording start to the FIRST
%                                  epoch's start ('offset' mode)
%                . stimStartAll    per-file 'HH:mm:ss.SSS' list ('manual' mode)
%                . epochsN         how many epochs
%                . epochDurationSec       one epoch's length, s
%                . epochBaselineSec       [t1 t2] baseline window on the epoch clock
%                . epochStimStartSec      when the stimulus starts on the epoch clock
%                . stimDurationSec        how long it lasts, s.  REQUIRED
%                . epochFinaleSec  [-t1 0] finale window, back from the epoch's END
%              % how a response is measured (passed to getNVCMetrics)
%                . nvcPeakGraceSec how far past the stimulus offset Pk may be found
%                                  (default 2)
%                . nvcAreaPcts     the cumulative-area levels the times are read at,
%                                  integers in (0,100) (default [10 50 90])
%                . nvcSignals      which traces to analyse (default all four)
%                . segNvcReturn    per-segment level cell (default all three)
%                . ppxNvcReturn    per-pixel level cell / gate (default [] = off)
%              % how much to trust a (segment, epoch) pair (passed to getNVCConfidence)
%                . nvcConfThreshold     geometric-mean threshold (default 0.6)
%                . nvcConfMinThreshold  weakest-factor threshold (default 0.2)
%                . nvcReturnScale       fReturn's scale, in noise units (default 5)
%                . nvcDevRules          add the five across-epoch deviation factors
%                                       (default false)
%                . nvcDevScale          their shared scale, in robust SDs (default 3)
%                . rejectFirstEpoch     1 zeroes epoch 1's confidence (default 1)
%              % how much to trust an EPOCH (passed to getNVCEpochTrust)
%                . nvcEpochAreaFrac  responding-area share an epoch needs (default 0.10)
%                . nvcEpochConnected count only the largest connected trusted region
%                                    (default true)
%              % what to do with the trusted set
%                . nvcEditEpochs     open the epoch editor (default false)
%                . nvcRepresentative average the trusted epochs and REPLACE THE
%                                    RECORDING WITH IT (default false)
%                . parforNvcSegments logical, default true
%                . parforNvcPixels   logical, default true - both are WORKER BOUNDS
%     fNames   cell array of *_t_BFI_r.mat (or *_s_BFI_r.mat) paths, IN THE ORDER
%              s.stimStartAll is written in.
%                . Optional workbench hooks in s (no-op when absent):
%                  s.stageFcn(stage,detail), s.cancelFcn()->tf.  Cancel is checked
%                  between files (never inside a parfor) - this step can block on the
%                  epoch editor.
%
%   OUTPUT FILES (side-effects) - NON-DESTRUCTIVE unless s.nvcRepresentative is true
%       *_t_BFI_d.mat   SOURCE   - read by the per-pixel path and by the collapse;
%                                  RE-SAVED ONLY BY THE COLLAPSE
%       *_t_BFI_r.mat   RESULTS  - RESULTS.nvc.* + the ns*/nd* metrics columns
%       *_t_BFI_s.mat   SETTINGS - field settings.runNVC added
%       *_rep_nvccuts.jpg       - where every epoch was cut, its responding area, its
%                                 confidence and whether the recording trusted it
%       *_rep_nvcconfidence.jpg - the [nSeg x nEp] confidence map, and per epoch the
%                                 share of segments each factor was weakest for
%       *_rep_nvcresponse.jpg   - every trusted epoch's field-mean and their average,
%                                 with the four windows and the three times marked
%       *_rep_nvctrials.jpg     - StRel and AucRel epoch by epoch across the session
%
%   THE COLLAPSE IS DESTRUCTIVE AND LOOKS LIKE IT.  With s.nvcRepresentative true the
%   trusted epochs are averaged, every signal and the flow cube are replaced by that
%   average on the epoch clock, the markers and the confidence are computed AGAIN on it
%   and stored as [nSeg x 1], and everything per-epoch is dropped.  results.timeStamp is
%   REMOVED, which arms the "no recording clock" guard so the step cannot run on its own
%   output; the explicit refusal below says so in words as well.  NOTHING PER-EPOCH
%   SURVIVES IN THE PRODUCT, by the author's ruling, so *_rep_nvccuts.jpg is the only
%   record there will ever be of which epochs went in - which is why it is written
%   before the overwrite, on every run, and names them.
%
%   TWO CONSEQUENCES OF THE COLLAPSE THAT THIS FILE CANNOT FIX, stated rather than
%   hidden.  Every OTHER step that globs *_BFI_r.mat - runVasomotion, runPulsatility,
%   runFitVasoreactivity - would accept a collapsed product and produce arithmetic about
%   one averaged 30-second epoch; and any branch those steps already wrote into this
%   product describes a recording that no longer exists.  Both are the author's call and
%   out of this step's scope.
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
%     Core/NVC/getNVCMetrics (the fourteen markers of one epoch of one trace),
%     Core/NVC/getNVCConfidence (the two confidence numbers), Core/NVC/getNVCEpochTrust
%     (the responding-area rule; MATLAB Image Processing Toolbox for bwlabel),
%     Core/NVC/editNVCEpochs (the epoch editor, only when s.nvcEditEpochs);
%     core LSCI library utilities.
%
% See also: getNVCMetrics, getNVCConfidence, getNVCEpochTrust, editNVCEpochs,
%           guiResponse, runFitVasoreactivity, runVasomotion, runPulsatility,
%           getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 06-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Where the stimulus repetitions are
% s.stimStartType='offset'; % 'offset' for a stimulus that starts a fixed time after
%                         % the recording begins, 'manual' for a list of clock times.
%                         % Either way this is the start of the FIRST EPOCH, not of
%                         % the stimulation inside it.
% s.stimOffset=0;         % seconds ('offset' mode)
% s.stimStartAll{1}='09:23:31.346';  % 'HH:mm:ss.SSS', one per file ('manual' mode)
% s.epochsN=20;           % how many repetitions
% s.epochDurationSec=30;  % how long one repetition is, seconds
% s.epochBaselineSec=[0,10];  % the quiet stretch at the start of each epoch that the
%                         % response is measured against, seconds from the epoch start
% s.epochStimStartSec=10; % when the stimulation actually starts within the epoch
% s.stimDurationSec=5;    % how long it lasts, seconds
% s.epochFinaleSec=[-5,0];% the stretch at the END of the epoch where flow is expected
%                         % to be back to baseline, measured back from the epoch's end
% %ADJUSTED IF NECESSARY - How a response is measured
% s.nvcPeakGraceSec=2;    % how long after the stimulation stops the peak may still
%                         % arrive, seconds
% s.nvcAreaPcts=[10,50,90]; % the response is timed by when 10, 50 and 90 % of the
%                         % accumulated flow increase had been delivered
% s.nvcSignals={'sData','gsData','dvsData','dvsDiameter'};  % which traces to analyse
% %ADJUSTED IF NECESSARY - Which levels to store
% s.segNvcReturn={'levels','amplitudes','times'};  % levels (baseline, finale, epoch
%                         % mean), amplitudes (stimulus, peak, area) and times.
% s.ppxNvcReturn=[];      % [] = off (the default).  Non-empty measures every pixel of
%                         % the field in every repetition and stores a map per epoch.
% %ADJUSTED IF NECESSARY - How much to trust one repetition of one segment
% s.nvcConfThreshold=0.6; % overall confidence a repetition needs, 0 to 1
% s.nvcConfMinThreshold=0.2; % the lowest any single check may score
% s.nvcReturnScale=5;     % how far flow may be from its own baseline at the end of a
%                         % repetition, in noise levels, before the score starts falling
% s.nvcDevRules=false;    % true also compares each repetition against the others of the
%                         % same segment
% s.nvcDevScale=3;        % how far it may sit from them, in robust noise levels
% s.rejectFirstEpoch=1;   % 1 never trusts the first repetition
% %ADJUSTED IF NECESSARY - Which repetitions the recording trusts
% s.nvcEpochAreaFrac=0.10;  % how much of the imaged area has to respond
% s.nvcEpochConnected=true; % count only the largest connected responding region
% s.nvcEditEpochs=false;  % true opens the editor and lets you change the decisions
% s.nvcRepresentative=false; % true REPLACES the recording with the average of the
%                         % trusted repetitions.  This cannot be undone.
% %ADJUSTED IF NECESSARY - Parallel execution
% s.parforNvcSegments=true;
% s.parforNvcPixels=true; % false measures one at a time in this MATLAB and starts no
%                         % parallel pool - slower, but it costs no worker processes.

%------------- BEGIN CODE --------------
function runNVC(s,fNames)

%THE INPUT IS A WHOLE RECORDING ON THE CONTRAST BRANCH.  A cardiac average ('_c_')
%collapses the recording clock, so there is no timeline left to place a stimulus on; a
%'_K_' product is a contrast, and this step applies no conversion because runBFI
%already did.  The guard is written positively, naming the two stages that ARE a whole
%recording of a flow index.
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
%inherit it from.
requireGeometry(s,numel(fNames));

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  This step can block
% on the epoch editor, so cancel is only checked between files.
rep=reportOpen(s,'NVC response',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        %SOURCE is opened for the per-pixel path and for the collapse - the collapse
        %averages the flow cube - and for nothing else.  Both are off by default.
        if ~isempty(s.ppxNvcReturn) || s.nvcRepresentative
            load(getProductPath(s.fName,'d'),'source')
        end
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        requireRecording(results,s.fName);

        sf=s;
        sf.epochStartSec=epochStarts(sf,results.timeStamp,fidx);

        %runNVC owns RESULTS.nvc entirely; rebuild it from scratch so no stale
        %sub-branch survives a re-run with different levels or a different geometry.
        if isfield(results,'nvc'), results=rmfield(results,'nvc'); end
        %AND IT OWNS ITS TABLE COLUMNS THE SAME WAY.  ns / nd / ng belong to this step
        %alone - the vasoreactivity fit is rs / rd and pulsatility is ps / pd - so every
        %one of them is dropped before pass 2 writes the current set.  Without this a
        %product processed by an earlier version of the step keeps that version's
        %columns beside the new ones, which is one file disagreeing with itself: a
        %reader would find the retired nsNEpochs sitting next to nsConf and no way to
        %tell which run either came from.
        results=dropOwnColumns(results,'sMetrics');
        results=dropOwnColumns(results,'dvsMetrics');

        % =============================================================
        % The cuts, on the PRIMARY clock (results.time).  gsData gets its own inside
        % measureSignal; everything the tree publishes at the root is on this one.
        % =============================================================
        [cut,epochStart,timeLoss,gaps]=cutIndex(results.time,sf);
        epochTime=results.time(cut(:,1))-results.time(cut(1,1));

        %TWO LAYOUTS, ONE TIME BASE.  layoutCalc is the COMPLETE marker set and it is
        %what every block is measured with, because getNVCConfidence needs all of it.
        %layoutTree is the user's own selection and its scalarNames are exactly the
        %columns that reach the saved tree.  See the header.
        sCalc=sf; sCalc.segNvcReturn={};
        layoutCalc=getNVCMetrics(epochTime,sCalc);
        layoutTree=getNVCMetrics(epochTime,sf);

        results.nvc.time      =single(layoutCalc.time(:));
        results.nvc.stimulus  =single(layoutCalc.u(:));
        results.nvc.epochStart=single(epochStart(:));
        results.nvc.timeLoss  =single(timeLoss(:));

        % =============================================================
        % PASS 1 - MEASURE.  One marker block and one confidence per signal.
        % =============================================================
        A=cell(1,numel(sf.nvcSignals));
        for kSig=1:numel(sf.nvcSignals)
            a=measureSignal(results,sf.nvcSignals{kSig},sf,sCalc, ...
                cut,timeLoss,layoutCalc,layoutTree);
            if isempty(a), continue, end
            results.nvc.esMetrics.(a.name)=a.tree;
            A{kSig}=rmfield(a,'tree');
        end
        A=A(~cellfun(@isempty,A));
        if isempty(A)
            error('runNVC:noTraces', ...
                ['%s carries none of %s, so there is nothing to cut.  Run the ' ...
                 'segmentation first.'],s.fName,strjoin(sf.nvcSignals,' / '));
        end

        % =============================================================
        % The recording-wide trusted set: RESPONDING AREA of ONE signal, and that
        % signal is the segmented flow trace.  sMap is the only spatial map the product
        % carries, so this is a sData rule; the cuts page says which one was used.
        % =============================================================
        kPri=primaryIndex(A);
        if ~isfield(results,'sMap') || isempty(results.sMap)
            error('runNVC:noSMap', ...
                ['%s carries no results.sMap, so the responding area of an epoch ' ...
                 'cannot be measured.  Run the segmentation first.'],s.fName);
        end
        [epochTrust,areaFrac]=getNVCEpochTrust(A{kPri}.C.trust,results.sMap,sf);

        ed=editorData(results,A{kPri},epochStart,layoutCalc,gaps,epochTrust,areaFrac);
        % THE EPOCH EDITOR, AND IT SITS HERE FOR A REASON.  The area rule has proposed a
        % trusted set and nothing has yet been built on it, so the operator's answer is
        % what pass 2, the pages and the collapse all read - one decision, taken once,
        % before anything depends on it.  `ed` is the whole contract, and it is the same
        % struct '_rep_nvccuts.jpg' draws from, so the window and the page cannot drift
        % apart.  OFF BY DEFAULT: a wrapper that blocks on a figure because a field was
        % omitted cannot be run from a script.
        if sf.nvcEditEpochs
            ed.epochTrust=editNVCEpochs(ed,struct('title',editorTitle(s.fName)));
        end

        results.nvc.epochTrust=logical(ed.epochTrust(:));
        results.nvc.areaFrac  =single(areaFrac(:));

        % =============================================================
        % PASS 2 - AGGREGATE.  The metrics tables, over the GLOBAL trusted set.
        % =============================================================
        results=writeTables(results,A,ed.epochTrust);

        % =============================================================
        % Per-pixel markers (gated by a non-empty s.ppxNvcReturn; ns prefix, since
        % source.data is a flow cube).  Every pixel goes through the same core the
        % segments did, so a map and a scalar cannot mean two different things.
        % SKIPPED WHEN THE COLLAPSE IS ON: the collapse recomputes the maps on the
        % averaged epoch, and the per-epoch ones it would throw away are the single
        % most expensive thing this step can be asked for.
        % =============================================================
        if ~isempty(sf.ppxNvcReturn) && ~sf.nvcRepresentative
            results.nvc.ppx=perPixel(reshape(source.data,[],size(source.data,3)), ...
                cut,results.sMap,timeLoss,layoutCalc,sf);
        end

        % =============================================================
        % The four pages.  ALL OF THEM BEFORE THE COLLAPSE - the cuts page is the only
        % provenance a collapsed product will ever have.
        % =============================================================
        drawCutsPage(rep,ed);
        drawConfidencePage(rep,ed);
        drawResponsePage(rep,results,A{kPri},ed,layoutCalc);
        drawTrialsPage(rep,A{kPri},ed);

        % =============================================================
        % The representative epoch - IN PLACE, destructive, and off by default.
        % =============================================================
        collapsed=sf.nvcRepresentative;
        if collapsed
            [results,source]=collapseToEpoch(results,source,A,cut, ...
                ed.epochTrust,timeLoss,layoutCalc,sf);
        end

        settings.runNVC=reportSettings(sf);
        reportWriting(rep);
        save(s.fName,'results','-v7.3','-nocompression');
        save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
        %SOURCE (_d) IS RE-SAVED ONLY BY THE COLLAPSE.  On every other path the flow
        %cube is preserved exactly as it was read.
        if collapsed
            save(getProductPath(s.fName,'d'),'source','-v7.3','-nocompression');
        end
        reportSaved(rep);
    end
end
reportClose(rep);

end

% =====================================================================
function a=measureSignal(results,sigName,sf,sCalc,cut,timeLoss,layoutCalc,layoutTree)
%measureSignal  PASS 1 for one signal: the marker block, then ONE confidence call.
%
%   THE CONFIDENCE CALL IS OUTSIDE THE parfor AND IT HAS TO BE.  Every optional
%   deviation factor measures one epoch against the OTHER epochs of the same unit, so a
%   worker holding one segment cannot compute one.  The block is what carries enough of
%   the problem, and building it is the only thing the workers do.
a=[];
if ~isfield(results,sigName) || isempty(results.(sigName)), return, end

%EACH SIGNAL BRINGS ITS OWN CLOCK.  gsData is the per-region guided trace at the raw
%frame rate, so it is cut on results.gsTime and gets its own layout and its own
%dropped-frame series; everything else shares the primary cut, which is already built.
if strcmp(sigName,'gsData')
    if ~isfield(results,'gsTime') || isempty(results.gsTime)
        error('runNVC:noGsTime', ...
            ['results.gsData is present but results.gsTime is not, so the guided ' ...
             'trace has no clock to cut it on.  Run the guided step again.']);
    end
    [sigCut,~,sigLoss]=cutIndex(results.gsTime,sf);
    epochT=results.gsTime(sigCut(:,1))-results.gsTime(sigCut(1,1));
    LCalc=getNVCMetrics(epochT,sCalc);
    LTree=getNVCMetrics(epochT,sf);
else
    sigCut=cut; sigLoss=timeLoss; LCalc=layoutCalc; LTree=layoutTree;
end

names=LCalc.blockNames;
sigMat=results.(sigName);
nSeg=size(sigMat,2);
nEp =size(sigCut,2);

%ONE parfor OVER SEGMENTS, epochs serial inside: 1940 against 20, and the layout and
%the cut index are shared by every iteration.  The sliced output is assigned on EVERY
%iteration, so nothing about the loop depends on which branch the core took.
%the workbench hooks are function handles bound to a GUI; they are transport rather
%than parameters, the core never reads them, and broadcasting them to workers is the
%one thing in s that can fail to serialise.  reportSettings strips exactly those two.
nw=0; if sf.parforNvcSegments, nw=Inf; end        % a worker BOUND, not a branch
sPar=reportSettings(sCalc);
M=nan(nSeg,numel(names),nEp,'single');
parfor (i=1:nSeg, nw)
    M(i,:,:)=nvcSegment(sigMat(:,i),sigCut,LCalc,sPar,names);
end

C=getNVCConfidence(M,names,sigLoss,LCalc,sf);

%THE TREE IS WRITTEN FROM layoutTree.scalarNames AND THE BLOCK IS INDEXED BY NAME, so
%this file types no scalar name of its own and a gated level cannot silently shift a
%column.  Conf, ConfMin and Trust are ordinary members with the branch prefix, which is
%what gets them into the explorer and the export sheet with no special case.
prefix=nvcPrefix(sigName);
T=struct();
for k=1:numel(LTree.scalarNames)
    nm=LTree.scalarNames{k};
    T.([prefix nm])=blockPlane(M,names,nm,nSeg,nEp);
end
T.([prefix 'Conf'])   =C.conf;
T.([prefix 'ConfMin'])=C.confMin;
T.([prefix 'Trust'])  =C.trust;

b.name    =sigName;
b.prefix  =prefix;
b.cut     =sigCut;
b.L       =LCalc;
b.M       =M;
b.C       =C;
b.nSeg    =nSeg;
b.nEp     =nEp;
b.timeLoss=sigLoss;
b.tree    =T;
a=b;
end

% =====================================================================
function rows=nvcSegment(y,cut,L,s,names)
%nvcSegment  One unit: every epoch of it measured, flattened to [1 x nMetric x nEp].
%   The core is a pure function of (series, layout, s) - the polarity is sign(Auc) of
%   the epoch in hand - so there is nothing per-segment to prepare and nothing a
%   worker has to be told beyond the trace.
E=y(cut);                                   % [nES x nEp], the epochs of this unit
nEp=size(E,2);
rows=nan(numel(names),nEp,'single');
for k=1:nEp
    m=getNVCMetrics(E(:,k),L,s);
    for j=1:numel(names), rows(j,k)=m.(names{j}); end
end
rows=reshape(rows,1,numel(names),nEp);
end

% =====================================================================
function results=writeTables(results,A,epochTrust)
%writeTables  PASS 2: the key set into sMetrics / dvsMetrics, as the MEAN over the
%   GLOBALLY-trusted epochs.
%
%   THE SAME EPOCHS FOR EVERY ROW - the author's ruling, and the reason there is no
%   nsNEpochs column.  nsConf is the per-segment number that says how much to trust the
%   row: that segment's mean confidence over the accepted epochs.  A segment trusted in
%   no epoch still gets a row, its markers are the mean over the global set and its
%   nsConf is low, which is exactly what a reader needs to discard it.
%   Row order matches the 1:nSeg loop, the same segment invariant runPulsatility
%   relies on.  BFI / std(BFI) are runBFI's columns and are intentionally not touched.
%   AN EMPTY TRUSTED SET IS NOT AN ERROR HERE: the per-epoch tree is the product and it
%   is complete either way, so the columns come out NaN and the cuts page says why.
for k=1:numel(A)
    a=A{k};
    switch a.name
        case 'sData',                    tbl='sMetrics';
        case {'dvsData','dvsDiameter'},  tbl='dvsMetrics';
        otherwise, continue              % gsData writes no table, as in runVasomotion
    end
    if ~isfield(results,tbl), continue, end
    keys=tableKeys(a.L);
    for j=1:numel(keys)
        v=double(blockPlane(a.M,a.L.blockNames,keys{j},a.nSeg,a.nEp));
        results.(tbl).([a.prefix keys{j}])= ...
            single(mean(v(:,epochTrust),2,'omitnan'));
    end
    results.(tbl).([a.prefix 'Conf'])= ...
        single(mean(double(a.C.conf(:,epochTrust)),2,'omitnan'));
end
end

% =====================================================================
function results=dropOwnColumns(results,tbl)
%dropOwnColumns  Remove every column this step owns from one metrics table.
%   The branch prefixes ns / nd / ng are runNVC's and nobody else's, and the naming
%   grammar puts a capital straight after them, so the match is exact rather than a
%   prefix sweep that would also take a column called `ndvi`.  BFI / std(BFI) are
%   runBFI's and every ps / pd / rs / rd column belongs to another step; none of them
%   is touched.
if ~isfield(results,tbl), return, end
T=results.(tbl);
if istable(T), have=T.Properties.VariableNames; else, have=fieldnames(T)'; end
mine=have(~cellfun(@isempty,regexp(have,'^n[sdg][A-Z]','once')));
if isempty(mine), return, end
if istable(T), results.(tbl)=removevars(T,mine); else, results.(tbl)=rmfield(T,mine); end
end

% =====================================================================
function keys=tableKeys(L)
%tableKeys  The four marker columns the tables carry, NAMED FROM THE SETTING.
%   The timing column is whichever requested level is nearest 50 %, so a column called
%   nsT50 always holds T50 and a protocol that asked for [20 80] gets nsT80 rather than
%   an nsT50 quietly holding something else.  The relatives are the three markers with
%   validity 0.95, and an area-weighted ROI mean of nsStRel is a meaningful number where
%   a mean of nsSt across two contrast products is not.
pcts=L.pcts;
if any(pcts==50)
    pT=50;
else
    [~,i]=min(abs(pcts-50));
    pT=pcts(i);
end
keys={'StRel','PkRel','AucRel',sprintf('T%d',pT)};
end

% =====================================================================
function P=blockPlane(M,names,want,nUnit,nEp)
%blockPlane  One marker of the block as a [nUnit x nEp] single, indexed BY NAME.
%   The block is always measured with the complete set, so a miss here is a programming
%   error rather than a level the user switched off - and it is named as one.
k=find(strcmp(names,want),1);
if isempty(k)
    error('runNVC:missingMetric', ...
        ['The marker block has no ''%s''.  The block is measured with the COMPLETE ' ...
         'layout; only the tree is gated.'],want);
end
P=reshape(M(:,k,:),nUnit,nEp);
end

% =====================================================================
function ppx=perPixel(D,cut,sMap,timeLoss,layoutCalc,s)
%perPixel  Per-epoch marker maps, [Y x X x nEp] each.  THE ONE PATH THAT READS SOURCE.
%
%   Every pixel goes through the core the segments did, on the same epoch clock and
%   with the same complete layout, so a per-pixel map and a per-segment scalar cannot
%   drift.  A BACKGROUND PIXEL IS WRITTEN BY NO ITERATION AND STAYS NaN, confidence
%   included: an unmeasured pixel has no response, and a zero confidence there would be
%   a measurement that was never made rather than one that failed.
%
%   IT SERVES THE COLLAPSE TOO.  Called with a single-column cut index on the averaged
%   cube it returns [Y x X] leaves, because reshape drops the trailing singleton - one
%   function, two shapes, no second spelling of the per-pixel path.
%
%   THE BLOCK IS THE COST.  [nPixel x 15 x nEp] single is 1.2 GB for a 1000 x 1000
%   field over 20 epochs, and it cannot be trimmed to the requested levels because the
%   confidence needs the complete set.  This is why s.ppxNvcReturn is [] by default.
[Y,X]=size(sMap);
npx=Y*X;
nEp=size(cut,2);
sPix=s; sPix.segNvcReturn=s.ppxNvcReturn;
LTree=getNVCMetrics(layoutCalc.time,sPix);
names=layoutCalc.blockNames;

sMapLin=sMap(:);
Mp=nan(npx,numel(names),nEp,'single');
nw=0; if s.parforNvcPixels, nw=Inf; end           % a worker BOUND, not a branch
sPar=reportSettings(sPix);
parfor (p=1:npx, nw)
    if sMapLin(p)==0, continue; end               % background stays NaN
    Mp(p,:,:)=nvcSegment(double(D(p,:)).',cut,layoutCalc,sPar,names);
end

Cp=getNVCConfidence(Mp,names,timeLoss,layoutCalc,s);
bg=sMapLin==0;
conf=Cp.conf;       conf(bg,:)=NaN;
confMin=Cp.confMin; confMin(bg,:)=NaN;

ppx=struct();
for k=1:numel(LTree.scalarNames)
    nm=LTree.scalarNames{k};
    ppx.(['ns' nm])=reshape(blockPlane(Mp,names,nm,npx,nEp),Y,X,nEp);
end
ppx.nsConf   =reshape(conf,Y,X,nEp);
ppx.nsConfMin=reshape(confMin,Y,X,nEp);
end

% =====================================================================
function [results,source]=collapseToEpoch(results,source,A,cut,epochTrust,timeLoss,L,s)
%collapseToEpoch  The representative epoch, written IN PLACE over the BFI product.
%
%   AFTER every marker and both confidence numbers exist for every (segment, epoch) -
%   which is the author's condition, not a convenience.  The trusted set chooses the
%   epochs, every signal and the flow cube are averaged over them on the epoch clock,
%   and the markers and the confidence are then computed AGAIN on that average and
%   stored as [nSeg x 1].
%
%   THE AVERAGED EPOCH IS NOT "EPOCH 1", so the re-measurement runs with
%   rejectFirstEpoch off; leaving it on would zero fFirst for every segment and make
%   the collapsed confidence identically zero.  Its dropped-frame allowance is the mean
%   loss of the epochs that went in, which is the loss the average actually carries.
%
%   THE METRICS TABLES ARE NOT RECOMPUTED, and the difference is real rather than an
%   oversight: the table is the MEAN OF THE MARKERS over the trusted epochs (pass 2,
%   which has already run), while esMetrics is now the MARKERS OF THE MEAN.  Both are
%   defined, they are not the same number, and the table's definition is the author's
%   ruling for every mode.
k=find(epochTrust);
if isempty(k)
    error('runNVC:noEpochsTrusted', ...
        ['No epoch of %s reached the responding-area threshold, so there is nothing ' ...
         'to average into a representative epoch.  Lower s.nvcEpochAreaFrac, relax ' ...
         'the confidence thresholds, or keep an epoch in the editor.'],s.fName);
end
sAvg=s; sAvg.rejectFirstEpoch=false;
sCalc=sAvg; sCalc.segNvcReturn={};
lossAvg=mean(double(timeLoss(k)));

esMetrics=struct();
for i=1:numel(A)
    a=A{i};
    nES=size(a.cut,1);
    yAvg=zeros(nES,a.nSeg);
    for j=1:numel(k)
        yAvg=yAvg+double(results.(a.name)(a.cut(:,k(j)),:));
    end
    yAvg=yAvg/numel(k);

    %the re-measurement: one epoch per unit, so it is 1/nEp of pass 1 and runs serially
    names=a.L.blockNames;
    LTree=getNVCMetrics(a.L.time,sAvg);
    M1=nan(a.nSeg,numel(names),1,'single');
    for iSeg=1:a.nSeg
        M1(iSeg,:,1)=nvcSegment(yAvg(:,iSeg),(1:nES)',a.L,sCalc,names);
    end
    C1=getNVCConfidence(M1,names,lossAvg,a.L,sAvg);

    T=struct();
    for kk=1:numel(LTree.scalarNames)
        nm=LTree.scalarNames{kk};
        T.([a.prefix nm])=blockPlane(M1,names,nm,a.nSeg,1);
    end
    T.([a.prefix 'Conf'])   =C1.conf;
    T.([a.prefix 'ConfMin'])=C1.confMin;
    T.([a.prefix 'Trust'])  =C1.trust;
    esMetrics.(a.name)=T;

    results.(a.name)=cast(yAvg,'like',results.(a.name));
    if strcmp(a.name,'gsData')
        results.gsTime=cast(a.L.time(:),'like',results.gsTime);
    end
end

%the flow cube, averaged over the same epochs and accumulated in its own class so a
%[Y x X x nES] double copy of a large field is never allocated
nES=size(cut,1);
acc=zeros(size(source.data,1),size(source.data,2),nES,'like',source.data);
for j=1:numel(k)
    acc=acc+source.data(:,:,cut(:,k(j)));
end
source.data=acc/numel(k);
if isfield(source,'time')
    source.time=cast(L.time(:),'like',source.time);
end

ppx=[];
if ~isempty(s.ppxNvcReturn)
    ppx=perPixel(reshape(source.data,[],nES),(1:nES)',results.sMap,lossAvg,L,sAvg);
end

%NOTHING PER-EPOCH SURVIVES.  epochStart, epochTrust, areaFrac and timeLoss are absent
%rather than empty, which is what makes the explorer's epoch axis simply not exist for
%a collapsed product, and what the re-run refusal keys on.
results.time=cast(L.time(:),'like',results.time);
results=rmfield(results,'nvc');
results.nvc.time     =single(L.time(:));
results.nvc.stimulus =single(L.u(:));
results.nvc.esMetrics=esMetrics;
if ~isempty(ppx), results.nvc.ppx=ppx; end

%REMOVING timeStamp ARMS THE "no recording clock" GUARD, so this step cannot run on its
%own output.  The named refusal in requireRecording says the same thing in words.
results=rmfield(results,'timeStamp');
end

% =====================================================================
function ed=editorData(results,a,epochStart,L,gaps,epochTrust,areaFrac)
%editorData  Everything the cuts page draws and the epoch editor will need, once.
%   ONE STRUCT so the page and the editor cannot drift apart - drawn twice, the old
%   picker and the old page differed by the colour of the epoch lines and nobody
%   noticed for a year.
[ed.trace,ed.traceLabel]=fieldMean(results,a.name);
ed.signal          =a.name;
ed.time            =double(results.time(:));
ed.epochStart      =double(epochStart(:));
ed.epochEnd        =ed.epochStart+L.time(end);
ed.timeLossSeries  =double(gaps(:));
ed.epochTrust      =logical(epochTrust(:));
ed.areaFrac        =double(areaFrac(:));
ed.conf            =a.C.conf;
ed.confMin         =a.C.confMin;
ed.trust           =a.C.trust;
ed.factors         =a.C.factors;
ed.factorNames     =a.C.factorNames;
%THE SENTENCES COME WITH THE FACTORS, from the core that built them.  How many checks
%there are is not a constant and the timing family is named from a setting, so a window
%holding its own list of words would go stale the first time a factor was added.
ed.factorWords     =a.C.factorWords;
end

% =====================================================================
function t=editorTitle(fName)
%editorTitle  How the editor names the recording it is open on.  The stem rather than the
%   path: a run of twenty files must never leave any doubt about which one is on screen,
%   and a title bar is not where a folder tree belongs.
[~,stem]=fileparts(regexprep(char(fName),'_[drs]\.mat$',''));
t=['Repetitions - ' stem];
end

% =====================================================================
function [y,what]=fieldMean(results,sigName)
%fieldMean  The recording's average trace of one signal, and what it is a trace of.
y=mean(double(results.(sigName)),2,'omitnan');
if strcmp(sigName,'dvsDiameter'), what='Diameter, px'; else, what='Flow, BFI'; end
end

% =====================================================================
function k=primaryIndex(A)
%primaryIndex  The signal the responding-area rule reads: sData, else dvsData, else
%   dvsDiameter.  results.sMap indexes the segmented flow trace, so the area rule is a
%   sData rule and the fallback order is the one the old mean trace used.  gsData is
%   never it: its rows are guided regions, not segments of sMap.
for want={'sData','dvsData','dvsDiameter'}
    k=find(cellfun(@(a) strcmp(a.name,want{1}),A),1);
    if ~isempty(k), return, end
end
error('runNVC:noAreaSignal', ...
    ['None of sData / dvsData / dvsDiameter was analysed, so the responding area of ' ...
     'an epoch cannot be measured.  s.nvcSignals must include a segmented flow trace.']);
end

% =====================================================================
function [cut,epochStart,timeLoss,gaps]=cutIndex(time,s)
%cutIndex  Where every epoch sits in a signal's own samples, [nES x nEp].
%
%   THE NEAREST FRAME TO EACH EDGE, which is the reason a dropped frame is a time LOSS
%   rather than a time SHIFT: the epoch still starts where the protocol says it starts,
%   and the gap shows up in timeLoss instead of moving every later epoch along.
%
%   EVERY EPOCH IS THE SAME NUMBER OF SAMPLES, and it has to be - the epoch clock, the
%   windows and the [nSeg x nEp] arrays are all shared.  The count is the TYPICAL span,
%   the median over the epochs, and not the shortest one: the last epoch of a recording
%   that stops on its own final sample is one frame short by arithmetic rather than by
%   any fault, and letting that shorten all twenty would throw a sample away from every
%   epoch to describe the edge of one.
time=double(time(:));
nT=numel(time);
startSec=s.epochStartSec(:)';
endSec  =startSec+s.epochDurationSec;

%[nT x 1] against [1 x nEp] expands to [nT x nEp]; the nearest frame to each edge is
%then the min DOWN the columns.
[~,startFrame]=min(abs(time-startSec),[],1);
[~,endFrame]  =min(abs(time-endSec),[],1);

step=median(diff(time));
if endSec(end)>time(end)+step
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

%the sampling gaps: everything beyond one nominal step is time the camera did not
%deliver.  Kept as a SERIES as well as a per-epoch sum, because the cuts page draws it
%and getNVCConfidence grades the sum as a FRACTION of the epoch span.
gaps=[0; diff(time)-step];
gaps(gaps<0)=0;
timeLoss=sum(gaps(cut),1);                  % [1 x nEp]
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
function requireRecording(results,fName)
%requireRecording  The three things a product must be before it can be epoched.
%   A COLLAPSED PRODUCT IS REFUSED IN WORDS.  s.nvcRepresentative writes the average of
%   the trusted epochs back over the recording and drops everything per-epoch, so a
%   results.nvc without an epochStart is a product that has already been collapsed.
%   Removing timeStamp already arms the guard below; this says which of the two it is.
if isfield(results,'nvc') && ~isfield(results.nvc,'epochStart')
    error('runNVC:collapsedProduct', ...
        ['%s is a COLLAPSED product: it already holds the average of its trusted ' ...
         'epochs in place of the recording, so there are no repetitions left to cut. ' ...
         'Re-process the recording from the contrast step if it has to be measured ' ...
         'again.'],fName);
end
if ~isfield(results,'time') || isempty(results.time)
    error('runNVC:noTime', ...
        '%s carries no results.time, so there is no recording clock to cut.',fName);
end
%The ABSOLUTE start of the recording is what places the stimulus, and only a step that
%read the raw recording clock can supply it.  runContrastInternalCycle writes no timeStamp -
%its product is one AVERAGED cardiac period - so a '_c_' product cannot be epoched by
%construction rather than by omission.
if ~isfield(results,'timeStamp') || isempty(results.timeStamp)
    error('runNVC:noTimeStamp', ...
        ['%s carries no results.timeStamp, so the stimulus cannot be placed on the ' ...
         'recording clock.  Cut a contrast-branch recording (*_t_BFI_r.mat or ' ...
         '*_s_BFI_r.mat); an averaged cardiac cycle has no timeline.'],fName);
end
end

% =====================================================================
function drawCutsPage(rep,ed)
%drawCutsPage  '_rep_nvccuts.jpg' - where every epoch was cut, and what the recording
%   made of it.
%
%   IN COLLAPSE MODE THIS PAGE IS THE ONLY PROVENANCE THERE WILL EVER BE.  Nothing
%   per-epoch survives in a collapsed product by the author's ruling, so this page has
%   to name the trusted epochs, their responding area and both confidence summaries -
%   and it is written on every run, before the overwrite.
%
%   A PAGE IS A BY-PRODUCT AND MUST NEVER KILL A RUN.  reportSave already swallows a
%   failed export; this swallows a failed DRAWING, because the analysis is complete by
%   the time the page is attempted.
fh=[];
try
    fh=reportFigure(rep,'nvccuts');
    tl=tiledlayout(fh,1,3,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1,[1 2]);
    drawEpochs(ax,ed);
    title(ax,sprintf('%d of %d epochs trusted  (responding area of %s, best %.0f %%)', ...
        sum(ed.epochTrust),numel(ed.epochTrust),ed.signal,100*max([0;ed.areaFrac])))

    ax=nexttile(tl,3);
    axis(ax,'off')
    nEp=numel(ed.epochTrust);
    yn={'no','yes'};
    lines=cell(1,nEp+2);
    lines{1}='ep   area    Conf  ConfMin  trusted';
    lines{2}='-------------------------------------';
    for k=1:nEp
        lines{k+2}=sprintf('%2d  %5.1f%%   %.2f    %.2f    %s', ...
            k,100*ed.areaFrac(k),median(double(ed.conf(:,k)),'omitnan'), ...
            median(double(ed.confMin(:,k)),'omitnan'),yn{1+ed.epochTrust(k)});
    end
    text(ax,0.0,1.0,strjoin(lines,newline),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName','Courier New','FontSize',max(6,min(10,round(560/(nEp+6)))), ...
        'Interpreter','none')
    title(ax,'Per epoch (median segment)')

    reportSave(rep,fh,'nvccuts');   % deletes the figure on every path of its own
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function drawEpochs(ax,ed)
%drawEpochs  The recording's mean trace, where every epoch was cut, whether it was
%   trusted, and the dropped-frame series beside it.  ONE function, so the page and
%   the session-3 editor draw the same picture.
nEp=numel(ed.epochStart);
lo=min(ed.trace); hi=max(ed.trace);
if ~(isfinite(lo)&&isfinite(hi)&&hi>lo), lo=0; hi=1; end
yyaxis(ax,'left')
plot(ax,ed.time,ed.trace,'-','Color',[0.25 0.35 0.55])
hold(ax,'on')
yb=lo-0.03*(hi-lo);
for i=1:nEp
    plot(ax,[ed.epochStart(i) ed.epochStart(i)],[lo hi],'--k')
    if ed.epochTrust(i), c=[0.15 0.60 0.25]; else, c=[0.85 0.20 0.15]; end
    plot(ax,[ed.epochStart(i) ed.epochEnd(i)],[yb yb],'-','Color',c,'LineWidth',3)
end
plot(ax,[ed.epochEnd(end) ed.epochEnd(end)],[lo hi],'--k')
ylabel(ax,ed.traceLabel)
hold(ax,'off')
%the bar hangs BELOW the trace, so the left limits are set from both together rather
%than left to a tight fit that would clip whichever was drawn last
ylim(ax,[lo-0.06*(hi-lo) hi+0.04*(hi-lo)])
yyaxis(ax,'right')
plot(ax,ed.time,ed.timeLossSeries,'-')
ylabel(ax,'Time loss, s')
xlabel(ax,'Time, s')
xlim(ax,[ed.time(1) ed.time(end)])
end

% =====================================================================
function drawConfidencePage(rep,ed)
%drawConfidencePage  '_rep_nvcconfidence.jpg' - THE PAGE THE REDESIGN IS FOR.
%
%   The [nSeg x nEp] confidence map says which epochs were poor; the stacked bar beside
%   it says WHY, by counting for every epoch the share of segments each factor was the
%   weakest for.  That is the whole difference between "this epoch is bad" and "this
%   epoch is bad BECAUSE its T90 ran into the finale".
%
%   THE FACTORS ARE NOT STORED IN THE PRODUCT - nine to fourteen more [nSeg x nEp]
%   arrays to answer a question ConfMin plus this page already answer.  The core
%   returns them so that this page and the editor can draw them.
fh=[];
try
    fh=reportFigure(rep,'nvcconfidence');
    tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
    nEp=numel(ed.epochTrust);

    ax=nexttile(tl,1);
    imagesc(ax,double(ed.conf))
    clim(ax,[0 1]); colormap(ax,parula); colorbar(ax)
    xlabel(ax,'Epoch'); ylabel(ax,'Segment')
    title(ax,sprintf('Confidence  (%.1f %% of pairs trusted)', ...
        100*mean(ed.trust,'all')))

    ax=nexttile(tl,2);
    [share,names]=weakestShare(ed);
    %ONE COLOUR PER FACTOR, SET EXPLICITLY.  There are nine factors by default and
    %fourteen with the deviation group on, and the default colour order has seven - so
    %left alone this page gives fT90 and fSnrFn the same orange, on the one page whose
    %entire purpose is to say WHICH check failed.
    b=bar(ax,1:nEp,share,'stacked');
    cm=turbo(numel(names));
    for k=1:numel(b), b(k).FaceColor=cm(k,:); b(k).EdgeColor='none'; end
    xlim(ax,[0.5 nEp+0.5]); ylim(ax,[0 1])
    xlabel(ax,'Epoch'); ylabel(ax,'Share of segments')
    legend(ax,names,'Location','eastoutside','Box','off','Interpreter','none')
    title(ax,'Which check was the weakest')

    reportSave(rep,fh,'nvcconfidence');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function [share,names]=weakestShare(ed)
%weakestShare  Per epoch, the share of segments each factor was the minimum for.
names=ed.factorNames;
nF=numel(names);
[nUnit,nEp]=size(ed.conf);
F=zeros(nUnit,nEp,nF);
for k=1:nF, F(:,:,k)=double(ed.factors.(names{k})); end
[~,w]=min(F,[],3);
share=zeros(nEp,nF);
for k=1:nF, share(:,k)=sum(w==k,1)'/max(1,nUnit); end
end

% =====================================================================
function drawResponsePage(rep,results,a,ed,L)
%drawResponsePage  '_rep_nvcresponse.jpg' - the trusted epochs and their average.
%
%   THERE IS NO FIT ON THIS PAGE ANY MORE.  What it draws is the field-mean of every
%   trusted epoch, their average on top, the four windows the markers were measured in
%   and the three times the median segment delivered its response by.
%
%   A PAGE MAY AVERAGE; THE ANALYSIS MAY NOT.  This average is drawn from the cut
%   traces at the moment the page is made and is stored nowhere - every marker in the
%   product was measured on an individual epoch.  Said here because the obvious way to
%   "tidy" this page is to feed it an averaging pre-pass, which is exactly the thing
%   the redesign removed.
fh=[];
try
    k=find(ed.epochTrust);
    if isempty(k), k=(1:numel(ed.epochTrust))'; end
    Y=double(results.(a.name));
    E=zeros(size(a.cut,1),numel(k));
    for j=1:numel(k), E(:,j)=mean(Y(a.cut(:,k(j)),:),2,'omitnan'); end
    t=L.time(:);

    fh=reportFigure(rep,'nvcresponse');
    tl=tiledlayout(fh,1,3,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1,[1 2]);
    hold(ax,'on')
    xregion(ax,t(find(L.blIdx,1)),t(find(L.blIdx,1,'last')), ...
        'FaceColor',[0.45 0.45 0.45],'FaceAlpha',0.08);
    xregion(ax,L.tStim,L.tStim+L.D,'FaceColor',[0.30 0.55 0.85],'FaceAlpha',0.14);
    xregion(ax,L.tStim+L.D,t(find(L.pkIdx,1,'last')), ...
        'FaceColor',[0.30 0.55 0.85],'FaceAlpha',0.06);
    xregion(ax,t(find(L.finIdx,1)),t(find(L.finIdx,1,'last')), ...
        'FaceColor',[0.45 0.45 0.45],'FaceAlpha',0.08);
    plot(ax,t,E,'-','Color',[0.70 0.70 0.75],'LineWidth',0.5)
    plot(ax,t,mean(E,2,'omitnan'),'-','Color',[0.20 0.35 0.60],'LineWidth',2)
    [tv,tn]=medianTimes(a,ed);
    for j=1:numel(tv)
        xline(ax,L.tStim+tv(j),'--','Color',[0.85 0.20 0.15],'Label',tn{j}, ...
            'LabelOrientation','horizontal','Interpreter','none');
    end
    hold(ax,'off')
    xlabel(ax,'Time in the epoch, s'); ylabel(ax,ed.traceLabel)
    axis(ax,'tight'); grid(ax,'on')
    title(ax,sprintf('%d trusted epochs of %s, and their average', ...
        numel(k),a.name))

    ax=nexttile(tl,3);
    axis(ax,'off')
    text(ax,0.0,1.0,responseLines(a,ed),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName','Courier New','FontSize',10,'Interpreter','none')
    title(ax,'Median segment, trusted epochs')

    reportSave(rep,fh,'nvcresponse');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function [tv,tn]=medianTimes(a,ed)
%medianTimes  The median segment's timing markers over the trusted epochs, on the
%   STIMULUS clock - which is why the caller adds tStim before drawing them.
k=ed.epochTrust; if ~any(k), k=true(size(k)); end
tv=zeros(1,numel(a.L.pcts)); tn=cell(1,numel(a.L.pcts));
for j=1:numel(a.L.pcts)
    nm=sprintf('T%d',a.L.pcts(j));
    v=double(blockPlane(a.M,a.L.blockNames,nm,a.nSeg,a.nEp));
    tv(j)=median(v(:,k),'all','omitnan');
    tn{j}=nm;
end
ok=isfinite(tv); tv=tv(ok); tn=tn(ok);
end

% =====================================================================
function txt=responseLines(a,ed)
%responseLines  The median segment's markers, written for a reader rather than a
%   parser.  Every name comes from the layout, so a protocol with five timing levels
%   prints five lines without a word changing here.
k=ed.epochTrust; if ~any(k), k=true(size(k)); end
lines={sprintf('epochs      %d of %d',sum(k),numel(k))};
for nm=[{'Bl','BlStd','St','StRel','Pk','PkRel','Auc','AucRel','Fn'}, ...
        arrayfun(@(p)sprintf('T%d',p),a.L.pcts,'UniformOutput',false)]
    v=double(blockPlane(a.M,a.L.blockNames,nm{1},a.nSeg,a.nEp));
    lines=[lines,{sprintf('%-11s %.4g',nm{1},median(v(:,k),'all','omitnan'))}]; %#ok<AGROW>
end
lines=[lines,{''}, ...
    {sprintf('%-11s %.3f','Conf',median(double(a.C.conf(:,k)),'all','omitnan'))}, ...
    {sprintf('%-11s %.3f','ConfMin',median(double(a.C.confMin(:,k)),'all','omitnan'))}, ...
    {sprintf('%-11s %.0f %%','trusted',100*mean(a.C.trust(:,k),'all'))}];
txt=strjoin(lines,newline);
end

% =====================================================================
function drawTrialsPage(rep,a,ed)
%drawTrialsPage  '_rep_nvctrials.jpg' - the response epoch by epoch across the session.
%
%   Averaging twenty repetitions before measuring them made habituation, drift and a
%   single bad repetition look identical - like nothing.  Here each epoch is a point:
%   the median across segments with the interquartile band around it, so a response
%   that fades over the session is a slope a reader sees rather than a number nobody
%   computed.
%
%   IT PLOTS THE RELATIVE AMPLITUDES.  StRel and AucRel are the two markers with
%   validity 0.95, and AucRel is the only one that carries amplitude and persistence in
%   a single number.  It is drawn from the marker BLOCK rather than from the tree, so
%   the page is the same whatever s.segNvcReturn was set to.
fh=[];
try
    fh=reportFigure(rep,'nvctrials');
    tl=tiledlayout(fh,2,1,'TileSpacing','compact','Padding','compact');
    ep=1:a.nEp;
    want={'StRel','Response, fraction of baseline'; ...
          'AucRel','Extra perfusion, s of resting flow'};
    for kRow=1:2
        v=double(blockPlane(a.M,a.L.blockNames,want{kRow,1},a.nSeg,a.nEp));
        p=prctile(v,[25 50 75],1);
        ax=nexttile(tl,kRow);
        hold(ax,'on')
        fill(ax,[ep flip(ep)],[p(1,:) flip(p(3,:))],[0.30 0.55 0.85], ...
            'FaceAlpha',0.18,'EdgeColor','none');
        plot(ax,ep,p(2,:),'-o','Color',[0.20 0.35 0.60],'MarkerFaceColor','w', ...
            'LineWidth',1.2,'MarkerSize',4)
        %the legend names what was DRAWN: a session in which every epoch was trusted
        %has no untrusted marker, and a third entry for it would be a key to nothing
        names={'middle half of the segments','median segment'};
        if any(~ed.epochTrust)
            plot(ax,ep(~ed.epochTrust),p(2,~ed.epochTrust),'x', ...
                'Color',[0.85 0.20 0.15],'LineWidth',1.5,'MarkerSize',10)
            names=[names,{'not trusted'}]; %#ok<AGROW>
        end
        hold(ax,'off')
        ylabel(ax,want{kRow,2}); grid(ax,'on')
        xlim(ax,[0.5 a.nEp+0.5])
        %the band is the point of the panel, so it gets room rather than being fitted
        %flush to the box where its edges read as the axes rather than as the data
        b=[min(p(1,:)) max(p(3,:))];
        if all(isfinite(b)) && b(2)>b(1), ylim(ax,b+0.06*diff(b)*[-1 1]); end
        if kRow==1
            legend(ax,names,'Location','best','Box','off')
            title(ax,sprintf('Response epoch by epoch (%s)',a.name))
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
function p=nvcPrefix(sigName)
%nvcPrefix  The branch prefix of a signal: ns flow, nd diameter, ng guided contrast.
switch sigName
    case 'dvsDiameter', p='nd';
    case 'gsData',      p='ng';
    otherwise,          p='ns';
end
end

% =====================================================================
function s=nvcDefaults(s)
%nvcDefaults  Resolve every default HERE, so the saved settings record what ran.
%   segNvcReturn defaults when absent OR empty; ppxNvcReturn defaults ONLY when the
%   field is ABSENT - an explicit [] must stay empty (per-pixel off), exactly like
%   runPulsatility's s.ppxPulsReturn.
%   THE CORES DEFAULT THE SAME VALUES INTERNALLY.  They are repeated here rather than
%   left to them because the settings file is the record of what a recording was
%   processed with, and a default that only exists inside a core is not in it.
if ~isfield(s,'nvcPeakGraceSec')||isempty(s.nvcPeakGraceSec), s.nvcPeakGraceSec=2; end
if ~isfield(s,'nvcAreaPcts')    ||isempty(s.nvcAreaPcts),     s.nvcAreaPcts=[10 50 90]; end
if ~isfield(s,'nvcSignals')     ||isempty(s.nvcSignals)
    s.nvcSignals={'sData','gsData','dvsData','dvsDiameter'};
end
if ~isfield(s,'segNvcReturn')   ||isempty(s.segNvcReturn)
    s.segNvcReturn={'levels','amplitudes','times'};
end
if ~isfield(s,'ppxNvcReturn'), s.ppxNvcReturn=[]; end
%the confidence.  Two thresholds and two scales, which with rejectFirstEpoch and the
%two epoch-area settings is the user's whole control over the trust decision.
if ~isfield(s,'nvcConfThreshold')   ||isempty(s.nvcConfThreshold),    s.nvcConfThreshold=0.6; end
if ~isfield(s,'nvcConfMinThreshold')||isempty(s.nvcConfMinThreshold), s.nvcConfMinThreshold=0.2; end
if ~isfield(s,'nvcReturnScale')     ||isempty(s.nvcReturnScale),      s.nvcReturnScale=5; end
if ~isfield(s,'nvcDevRules')        ||isempty(s.nvcDevRules),         s.nvcDevRules=false; end
if ~isfield(s,'nvcDevScale')        ||isempty(s.nvcDevScale),         s.nvcDevScale=3; end
if ~isfield(s,'rejectFirstEpoch')   ||isempty(s.rejectFirstEpoch),    s.rejectFirstEpoch=1; end
%the epoch-level rule
if ~isfield(s,'nvcEpochAreaFrac') ||isempty(s.nvcEpochAreaFrac),  s.nvcEpochAreaFrac=0.10; end
if ~isfield(s,'nvcEpochConnected')||isempty(s.nvcEpochConnected), s.nvcEpochConnected=true; end
%OFF BY DEFAULT HERE, ON BY DEFAULT IN A PROTOCOL.  A wrapper that blocks on a figure
%because a field was omitted cannot be run from a script; the launcher and the workbench
%preset both switch the editor on, which is where a person is choosing.
if ~isfield(s,'nvcEditEpochs')   ||isempty(s.nvcEditEpochs),    s.nvcEditEpochs=false; end
if ~isfield(s,'nvcRepresentative')||isempty(s.nvcRepresentative),s.nvcRepresentative=false; end
s.nvcEditEpochs   =logical(s.nvcEditEpochs);
s.nvcRepresentative=logical(s.nvcRepresentative);
%PARALLELISM IS OPTIONAL.  Each switch is a BOUND on its parfor (Inf workers or 0),
%never a branch: parfor(...,0) runs the identical loop body serially IN THE CLIENT and
%starts no pool at all.
if ~isfield(s,'parforNvcSegments')||isempty(s.parforNvcSegments), s.parforNvcSegments=true; end
if ~isfield(s,'parforNvcPixels')  ||isempty(s.parforNvcPixels),   s.parforNvcPixels=true; end

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
%   The geometry is this step's own input - there is no producing step to inherit it
%   from - so a missing or mis-sized entry has to be a named message rather than a
%   subscript error inside the loop.
need={'epochsN','epochDurationSec','epochStimStartSec','stimDurationSec'};
for k=1:numel(need)
    if ~isfield(s,need{k})||isempty(s.(need{k}))||~isscalar(s.(need{k}))|| ...
            ~isfinite(s.(need{k}))||s.(need{k})<=0
        unit='seconds';
        if strcmp(need{k},'epochsN'), unit='epochs'; end
        error('runNVC:missingGeometry', ...
            's.%s is required - a positive number of %s.',need{k},unit);
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
