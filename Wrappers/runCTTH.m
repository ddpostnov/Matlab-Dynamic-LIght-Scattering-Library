%runCTTH  Model-free bolus transit: delay, mean transit time and the width of it
%
%   runCTTH(s,fNames) loads every '*_b_I_r.mat' bolus product in fNames and measures,
%   per segment and optionally per pixel, when the tracer arrived, how long it took to
%   pass and how widely spread those times were - each of them referenced to the
%   recording's OWN arterial input, which the step derives for itself.
%
%   THERE IS NO LANDMARK MACHINERY LEFT, AND NO FILTERING.  The step this replaces
%   hampel-despiked, median-filtered and Savitzky-Golay smoothed every trace, then read
%   five landmarks off it with findpeaks.  Measured (research-ctth.md section 10) the
%   whole filtering chain moves T50 by 0.0035 s against a Delay spread of 0.31 s - about
%   1 % - and the per-frame spatial median it also applied is what made the per-pixel
%   pass expensive.  The markers here are integrals and interpolated crossings, which
%   average noise rather than amplify it, so all of it is gone: s.medSpace, s.medTime,
%   s.sgOrder, s.sgFrame, s.hampWin, s.hampSig, s.slopeWin, s.promFrac, s.minStep and
%   s.baseWin retire with the machinery they served, and nothing reads them afterwards.
%
%   THE MARKER SET IS TWENTY-ONE SCALARS PER UNIT, and its contract is getBolusMetrics.
%   Three of them are transit times and they are the only ones comparable between
%   recordings: Delay (a difference of 50 % times), Mtt (a difference of first moments)
%   and Cth (the root of a difference of variances - the H in CTTH).
%
%   WHAT IS NOT HERE, BY THE AUTHOR'S DECISION OF 07-Aug-2026.  On this preparation the
%   recirculation OVERLAPS a ~4 s infusion, so the first pass never separates: the field
%   mean settles at 75 % of its own peak and the only trough after the peak is 0.9 %
%   prominent.  An area with no tail to cut is undefined, so there is NO first-pass area,
%   NO rCBV and NO recirculation onset - the old tReCircB was measuring the first noise
%   wiggle past the maximum, at a median time that WAS the peak.  What replaces them is
%   BvRel, the plateau step relative to the input's, which is available BECAUSE the first
%   pass failed to separate rather than in spite of it.
%
%   BvRel IS A RELATIVE MAP AND IT IS EXPOSED TO HAZE.  Widefield fluorescence floods
%   everything with out-of-plane and scattered light from the vessels, and measured on
%   the reference recording a vessel WALL reads 0.42 against a lumen's 0.45 where the
%   true blood-volume fraction differs by more than an order of magnitude.  So BvRel
%   orders units within one recording and is not a blood-volume map until background
%   removal has run.  Every TIME marker is unaffected: a haze constant in time cancels in
%   the Step-normalisation.
%
%   Cth IS NaN BELOW A FLOOR THE STEP CALIBRATES PER RECORDING, and that is the author's
%   decision R4-3 rather than a caveat in a document.  The width comes out of a
%   difference of two numbers that agree to one per cent, so a 1 % error in the plateau
%   level moves it by more than the answer - and no per-trace statistic can see that,
%   because it is systematic.  getBolusCthFloor pushes kernels of known spread through
%   this recording's own input on its own time base and reads back the smallest one that
%   survives; below it Cth is NaN, and results.ctth.cthFloor is stored so a Cth column
%   always has its error bar beside it.
%
%   THE INPUT IS DERIVED, FROM SEGMENTS, ONCE.  No ROI is drawn (author's decision R4-2):
%   getBolusInput takes the earliest-arriving units AFTER an amplitude gate, which is the
%   part that is not optional - ungated, "earliest arriving" selects the DIMMEST units,
%   whose 10 % threshold sits inside their own noise, and returns an input LATER than the
%   whole field.  It runs on the SEGMENT traces even when per-pixel maps were asked for:
%   the input is a property of the recording, and deriving it twice, once per unit size,
%   would be two origins for one measurement.
%
%   IT SAYS SO WHEN THE FRAMES HAVE BEEN CLEANED, AND IT DOES NOT REFUSE.  Background
%   removal offers to take the glow around the vessels off every frame as well as off the
%   picture, and doing that moves every time below by a knowable and large amount without
%   changing how the result looks.  Per the author (2026-08-07) that combination is the
%   user's to make, so this step runs it and warns once per file - see
%   warnIfCubeWasCleaned at the foot of this file for the numbers.  There is deliberately
%   no error and no registry incompatibility.
%
%   THE OPTIONAL MODEL IS NOT BUILT.  The author's brief asks for a gamma-variate fit as
%   a selectable second method beside the model-free one.  fitGammaVariatePerPixel, which
%   the launcher used to call, is not in this repository, and nothing in the research
%   session fitted a gamma - so there is no measurement to ship it against and no
%   s.ctthMethod here.  A one-valued setting would be the stub this branch has been
%   removing, so the method arrives with the fit or not at all.
%
%   THE RUN, PER FILE
%     load the triplet, and warn if the frames were cleaned
%     derive the arterial input from the SEGMENT traces
%     SETUP the layout on it - it measures the input with the same code path as a unit
%     calibrate this recording's Cth floor through that input
%     per signal: the [nUnit x nMetric] marker block, then ONE confidence call on it
%     NaN every Cth below the floor
%     the tree, the five table columns, the optional per-pixel maps, four pages, save
%
%   TWO LAYOUTS, ONE TIME BASE, AND THEY ARE NOT OPTIONAL.  getBolusConfidence needs
%   Step, BlStd, TPeak, Delay, the first level time and the three lowercase shape checks,
%   so a run with s.segCtthReturn {'times'} would have no confidence at all.  SETUP is
%   therefore called TWICE: once with the COMPLETE marker set, which is what every block
%   is measured with, and once with the user's own selection, which is what NAMES the
%   tree columns.  This file types no marker name of its own.
%
%   OUTPUT TREE  (RESULTS.ctth; every float SINGLE)
%       .time       [nT x 1]    the recording clock the input is quoted on
%       .aif        [nT x 1]    the derived input, baseline-referenced
%       .aifUnits   [nSeg x 1]  LOGICAL - which segments it came from
%       .cthFloor   scalar      the smallest spread this recording resolves, s
%       .cthUsable  scalar      LOGICAL - false when it resolves none at all
%       .inputM1    scalar      the input's OWN centroid, s.  STORED ON PURPOSE: it is
%                               how much of every absolute time was the infusion, and a
%                               reader asking why T50 is not comparable between
%                               recordings needs it
%       .esMetrics  <- THE PRODUCT.  One sub-struct per analysed signal, each field
%                      [nUnit x 1] and row-aligned to sMetrics / dvsMetrics:
%           .sData       .<bs*>    .dvsData  .<bs*>    .dvsDiameter .<bd*>
%                      the requested markers, then ALWAYS <prefix>Conf, <prefix>ConfMin
%                      and <prefix>Trust - never gated by s.segCtthReturn, because a
%                      level selection that removed them would leave the tree unreadable
%                      rather than smaller
%       .ppx        (only when s.ppxCtthReturn is non-empty) [Y x X] per requested
%                      marker, plus bsConf and bsConfMin.  Background pixels are NaN in
%                      ALL of them, confidence included: an unmeasured pixel has no
%                      transit time, and a zero confidence there would be a measurement
%                      that was never made rather than one that failed
%
%     A DIAMETER TRACE GETS THE bd PREFIX AND MEANS SOMETHING ELSE.  dvsDiameter is a
%     width in pixels, so its "arrival" is when the vessel changed calibre and not when
%     the tracer reached it.  It is measured because the machinery is agnostic and the
%     number is occasionally wanted; a bd* column must never be pooled with a bs* one.
%
%   METRICS-TABLE DUPLICATION (a small key set, same names as the tree, single)
%     results.sMetrics   (from sData)      : bsDelay bsMtt bsCth bsBvRel bsConf
%     results.dvsMetrics (from dvsData)    : the same five bs*
%                        (from dvsDiameter): the five bd*
%     bsCth is NaN below the floor, exactly as in the tree.  THE LEVELS GATE THE TREE,
%     NOT THE TABLE - s.segCtthReturn decides only what reaches the saved tree, and a
%     table column must not appear and disappear with a display setting.  Same rule as
%     runNVC, runPulsatility and runFitVasoreactivity.
%
%   INPUTS
%     s        parameter struct with fields
%              % how a curve is read
%                . ctthLevelPcts   the levels the times are read at, whole percentages
%                                  in (0,100)  (default [10 25 50 75 90])
%                . ctthPlateauSec  the end window the final level is averaged over, s
%                                  (default 1)
%                . ctthGuardSec    quiet time left between the pre-injection window and
%                                  the arrival, s (default 0.5)
%                . ctthSlopeSec    window the two slope markers are measured over, s
%                                  (default 0.1)
%                . ctthBaselineSec [t1 t2] pre-injection window; empty (the default)
%                                  derives it from the input's own onset
%              % where the arterial input comes from
%                . ctthInputAmpPct amplitude gate, percentile of the plateau steps
%                                  (default 75)
%                . ctthInputPct    share of the gated units taken, percentile of their
%                                  arrival (default 5)
%              % how much to trust a unit
%                . ctthConfThreshold    geometric-mean threshold (default 0.6)
%                . ctthConfMinThreshold weakest-factor threshold (default 0.2)
%                . ctthSettleFrac       how much the curve may still be rising at the
%                                       end, as a share of its own step (default 1)
%                . ctthInputScale       how far ahead of the artery a unit may fill, s
%                                       (default 0.5)
%              % the width floor
%                . ctthCthTol      relative error a probe may carry and still count as
%                                  resolved (default 0.25)
%                . ctthCthProbes   the probe spreads, s (default 0.1:0.1:3)
%              % what is kept
%                . segCtthReturn   cell subset of {'levels','amplitudes','times',
%                                  'transit','shape'} (default all five)
%                . ppxCtthReturn   per-pixel level cell / GATE ({} = off, the default)
%                . parforCTTHPixels logical, default true - a WORKER BOUND on the
%                                  per-pixel loop, never a branch
%     fNames   cell array of '*_b_I_r.mat' paths - the RESULTS member of each bolus
%              product.  The SOURCE cube and the SETTINGS are named from it.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUT FILES (side-effects; SOURCE is never written)
%       *_b_I_r.mat   RESULTS  - RESULTS.ctth.* and the bs*/bd* metrics columns
%       *_b_I_s.mat   SETTINGS - field settings.runCTTH added
%       *_rep_ctth-input.jpg      where the arterial curve came from, and how far it
%                                 sits from the whole field
%       *_rep_ctth-markers.jpg    the field curve with its times and its plateau
%       *_rep_ctth-confidence.jpg the seven checks, and which one was the weakest
%       *_rep_ctth-maps.jpg       Delay, Mtt, Cth, BvRel and Conf as pictures
%
%   EXAMPLE
%     s.libraryFolder=libraryFolder;
%     s.ctthLevelPcts=[10 25 50 75 90];
%     s.segCtthReturn={'levels','amplitudes','times','transit','shape'};
%     s.ppxCtthReturn={};              % non-empty also makes a picture of every marker
%     fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat');
%     runCTTH(s,fNames(:));
%
%   DEPENDS ON
%     Core/Intensity/getBolusInput, getBolusMetrics, getBolusMoments,
%     getBolusConfidence, getBolusCthFloor; core LSCI library utilities.
%
% See also: getBolusMetrics, getBolusConfidence, getBolusInput, getBolusCthFloor,
%           runIntensityBolus, runBackgroundRemoval, runSegmentation, getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - How a curve is read
% s.ctthLevelPcts=[10,25,50,75,90]; % the tracer has filled this far, in per cent, at
%                         % each time that gets reported
% s.ctthPlateauSec=1;     % how much of the end of the recording is averaged to read the
%                         % final level, seconds
% s.ctthGuardSec=0.5;     % quiet time left between the pre-injection period and the arrival
% s.ctthSlopeSec=0.1;     % how long a stretch the steepest rise and fall are measured over
% %ADJUSTED IF NECESSARY - Where the arterial curve comes from
% s.ctthInputAmpPct=75;   % how bright a vessel has to be to count towards it
% s.ctthInputPct=5;       % how many of the earliest-filling vessels make it up
% %ADJUSTED IF NECESSARY - How much to trust a vessel
% s.ctthConfThreshold=0.6;    % how good all the checks together have to be
% s.ctthConfMinThreshold=0.2; % how good the single worst check has to be
% s.ctthSettleFrac=1;     % how much the curve may still be rising at the end
% s.ctthInputScale=0.5;   % how far ahead of the artery a vessel may fill, seconds
% s.ctthCthTol=0.25;      % how close the spread check has to come before the recording
%                         % is judged long enough for it
% %ADJUSTED IF NECESSARY - Which levels to store
% s.segCtthReturn={'levels','amplitudes','times','transit','shape'};
% s.ppxCtthReturn={};     % {} = off. Non-empty measures every pixel and stores a picture
% s.parforCTTHPixels=true; % false measures one at a time and starts no parallel pool

%------------- BEGIN CODE --------------
function runCTTH(s,fNames)

if ~all( cellfun(@(f) isempty(f) || contains(f,'_b_I_r.mat'), fNames(:)) )
    error('runCTTH:badInput', ...
        ['One or more *non-empty* entries are not a bolus product.  This step reads ' ...
         'the RESULTS member of a "_b_I" triplet - list them with ' ...
         'getFileNamesList(resultsFolder,''*_b_I_r.mat'').']);
end

s=ctthDefaults(s);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Transit time',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end          % cooperative cancel between files
    if ~isempty(fNames{fidx})

        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        %SOURCE is opened for the per-pixel path and for nothing else, and that path is
        %off by default - so the usual run never pays the multi-gigabyte load.
        if ~isempty(s.ppxCtthReturn)
            load(getProductPath(s.fName,'d'),'source');
        end
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        warnIfCubeWasCleaned(settings,s.fName);
        requireBolus(results,s.fName);

        %runCTTH owns RESULTS.ctth entirely, and its own table columns the same way.
        %Both are dropped before anything is written, so a product processed by an
        %earlier version of the step cannot keep that version's columns beside the new
        %ones - which is one file disagreeing with itself.  The retired landmark set
        %(t0B, tUpslopeB, tPeakB, tReCircB, tBaselineB and the four v*B) goes with them.
        results=dropOwnBranch(results);
        results=dropOwnColumns(results,'sMetrics');
        results=dropOwnColumns(results,'dvsMetrics');

        time=double(results.time(:));

        % =============================================================
        % The arterial input, from the SEGMENT traces.  One origin per recording.
        % =============================================================
        inputName=inputSignal(results,s.fName);
        [aif,iInfo]=getBolusInput(time,double(results.(inputName)),s);

        %TWO LAYOUTS, ONE TIME BASE.  layoutCalc is the COMPLETE marker set and is what
        %every block is measured with, because the confidence needs all of it;
        %layoutTree is the user's own selection and its scalarNames are exactly the
        %columns that reach the saved tree.  See the header.
        sCalc=s; sCalc.segCtthReturn={};
        L    =getBolusMetrics(time,aif,sCalc);
        LTree=getBolusMetrics(time,aif,s);

        % =============================================================
        % This recording's own width floor, calibrated through its own input.
        % =============================================================
        F=getBolusCthFloor(time,aif,s);
        if ~F.usable
            %NOT AN ERROR.  Every time marker survives a record too short for a width,
            %so the step reports the times and says the width is not available.
            warning('runCTTH:cthBelowFloor', ...
                'For %s, %s.',shortName(s.fName),F.words);
        end

        results.ctth.time     =single(time);
        results.ctth.aif      =single(aif);
        results.ctth.aifUnits =logical(iInfo.units(:));
        results.ctth.cthFloor =single(F.floor);
        results.ctth.cthUsable=logical(F.usable);
        results.ctth.inputM1  =single(L.in.m1);

        % =============================================================
        % Per signal: the marker block, one confidence call, then the floor.
        % =============================================================
        A=cell(1,numel(s.ctthSignals));
        for kSig=1:numel(s.ctthSignals)
            a=measureSignal(results,s.ctthSignals{kSig},s,sCalc,L,LTree,F);
            if isempty(a), continue, end
            results.ctth.esMetrics.(a.name)=a.tree;
            A{kSig}=rmfield(a,'tree');
        end
        A=A(~cellfun(@isempty,A));
        if isempty(A)
            error('runCTTH:noTraces', ...
                ['%s carries none of %s, so there is nothing to measure.  Run the ' ...
                 'segmentation first.'],s.fName,strjoin(s.ctthSignals,' / '));
        end

        results=writeTables(results,A);

        % =============================================================
        % Per-pixel maps, gated by a non-empty s.ppxCtthReturn.
        % =============================================================
        if ~isempty(s.ppxCtthReturn)
            results.ctth.ppx=perPixel(source.data,results.cMask>0,L,F,s);
        end

        % =============================================================
        % The four pages.
        % =============================================================
        drawInputPage(rep,results,aif,iInfo,L,inputName);
        drawMarkersPage(rep,results,A,L);
        drawConfidencePage(rep,A);
        drawMapsPage(rep,results,A,F);

        settings.runCTTH=reportSettings(s);
        reportWriting(rep);
        save(s.fName,'results','-v7.3','-nocompression');
        save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
        clear source
    end
end
reportClose(rep);

end


% ====================================================================== %
function warnIfCubeWasCleaned(settings,fName)
%warnIfCubeWasCleaned  Say once, per file, that these landmarks came off cleaned frames.
%
%   IT WARNS AND IT DOES NOT BLOCK, and that is the author's ruling of 2026-08-07:
%   "If a user chooses apply to source - they can still run CTTH.  It is their choice.
%   What is correct is another question."  An earlier draft of the background-removal
%   work specified a hard error here and it was overruled, so there is no error, no
%   registry incompatibility and nothing greyed out.
%
%   WHAT WOULD OTHERWISE BE LOST IS THAT THE RESULT LOOKS COMPLETELY NORMAL.  R1
%   measured, on the same segments of the same recording, an arrival shifted by
%   +2171 ms - 217 frames - a segment-to-segment CV of 95% in the first-pass area, and
%   no calibre at which it becomes safe: 30-50 um vessels had 12.5% of segments within
%   one frame.  Nothing about the numbers says so afterwards.  So the choice is
%   recorded in the settings sidecar, which is what makes such a run auditable, and
%   this line is what makes it noticed.
%
%   IT IS AN ORDINARY MATLAB WARNING AND NOT A REPORT LINE.  The reporting module is
%   three lines per recording by design and this is not a stage; runIntensityBolus
%   already uses warning this way, so it is not a new idiom.  Reading a settings field a
%   step wrote is a gate on a live CHOICE, not a back-compatibility fallback - a product
%   either carries that settings group or it does not, and both states are current.
if ~isfield(settings,'runBackgroundRemoval'), return; end
b=settings.runBackgroundRemoval;
if ~isfield(b,'applyToSource') || ~isequal(b.applyToSource,true), return; end
[~,n,e]=fileparts(char(fName));
warning('runCTTH:cleanedCube', ...
    ['The frames of %s have had the glow around the vessels taken off them.  ' ...
     'Measured on a recording processed both ways, that moves the arrival of the ' ...
     'bolus about two seconds later and scatters the first-pass area from one ' ...
     'vessel to the next; the times below are still produced and still look ' ...
     'ordinary.'], [n e]);
end


% ====================================================================== %
function a=measureSignal(results,sigName,s,sCalc,L,LTree,F)
%measureSignal  One signal: the [nUnit x nMetric] block, then ONE confidence call.
%
%   THE CONFIDENCE IS COMPUTED ON THE WHOLE BLOCK AND NOT PER UNIT, which here is a
%   matter of one call rather than of necessity - unlike the NVC step there is no epoch
%   axis and therefore no factor that needs the other columns.  It is written this way
%   because the core takes a block and returning to a per-unit call would be a second
%   spelling of the same arithmetic.
a=[];
if ~isfield(results,sigName) || isempty(results.(sigName)), return, end

Y=double(results.(sigName));
names=L.blockNames;
nUnit=size(Y,2);

M=nan(nUnit,numel(names),'single');
for i=1:nUnit
    m=getBolusMetrics(Y(:,i),L,sCalc);
    for j=1:numel(names), M(i,j)=m.(names{j}); end
end

C=getBolusConfidence(M,names,L,s);

%R4-3, AND IT IS ONE LINE.  A width below what this recording can resolve is not a
%small width, it is no measurement - so it is NaN in the block, which is what the tree
%and the table are both written from, and the floor it failed is stored beside them.
M=gateCth(M,names,F);

%THE TREE IS WRITTEN FROM LTree.scalarNames AND THE BLOCK IS INDEXED BY NAME, so this
%file types no marker name of its own and a gated level cannot silently shift a column.
%Conf / ConfMin / Trust are ordinary members with the branch prefix and are NEVER gated,
%which is what gets them into the explorer and the export sheet with no special case.
prefix=ctthPrefix(sigName);
T=struct();
for k=1:numel(LTree.scalarNames)
    nm=LTree.scalarNames{k};
    T.([prefix nm])=blockCol(M,names,nm);
end
T.([prefix 'Conf'])   =C.conf;
T.([prefix 'ConfMin'])=C.confMin;
T.([prefix 'Trust'])  =C.trust;

a.name  =sigName;
a.prefix=prefix;
a.L     =L;
a.M     =M;
a.C     =C;
a.nUnit =nUnit;
a.tree  =T;
end


% ====================================================================== %
function M=gateCth(M,names,F)
%gateCth  NaN every Cth the recording cannot resolve.  See getBolusCthFloor.
k=find(strcmp(names,'Cth'),1);
if isempty(k), return, end
v=M(:,k);
v(~(v>=F.floor))=NaN;           % written positively, so a NaN Cth stays NaN
M(:,k)=v;
end


% ====================================================================== %
function results=writeTables(results,A)
%writeTables  The five key columns into sMetrics / dvsMetrics.
%   Row order matches the 1:nUnit loop, the same segment invariant runPulsatility and
%   runNVC rely on.  bsCth is NaN below the floor because the block it is read from
%   already is.
for k=1:numel(A)
    a=A{k};
    switch a.name
        case 'sData',                   tbl='sMetrics';
        case {'dvsData','dvsDiameter'}, tbl='dvsMetrics';
        otherwise, continue
    end
    if ~isfield(results,tbl), continue, end
    for key={'Delay','Mtt','Cth','BvRel'}
        results.(tbl).([a.prefix key{1}])=blockCol(a.M,a.L.blockNames,key{1});
    end
    results.(tbl).([a.prefix 'Conf'])=a.C.conf;
end
end


% ====================================================================== %
function results=dropOwnBranch(results)
%dropOwnBranch  Remove results.ctth and the retired per-pixel landmark maps.
%   The old step wrote imgT0B / imgTUpslopeB / imgTPeakB / imgTReCircB / imgTBaselineB
%   and their imgV*B twins straight onto results.  Nothing reads them now, and a product
%   carrying both them and results.ctth would be one file describing two different
%   analyses.
if isfield(results,'ctth'), results=rmfield(results,'ctth'); end
have=fieldnames(results)';
old=have(~cellfun(@isempty,regexp(have,'^img[TV](0|Upslope|Peak|ReCirc|Baseline)B$','once')));
if ~isempty(old), results=rmfield(results,old); end
end


% ====================================================================== %
function results=dropOwnColumns(results,tbl)
%dropOwnColumns  Remove every column this step owns from one metrics table.
%   TWO GENERATIONS OF THEM.  bs / bd are runCTTH's branch prefixes and nobody else's -
%   the naming grammar puts a capital straight after them, so the match is exact rather
%   than a prefix sweep that would also take a column called `bsl`.  The second pattern
%   is the RETIRED landmark set of the step this replaces; it is named explicitly
%   because those columns have no prefix grammar to match on, and leaving them would put
%   a tReCircB that was measuring noise next to a Delay that is not.
if ~isfield(results,tbl), return, end
T=results.(tbl);
if istable(T), have=T.Properties.VariableNames; else, have=fieldnames(T)'; end
mine=have(~cellfun(@isempty,regexp(have,'^b[sd][A-Z]','once')));
old =have(ismember(have,{'t0B','tUpslopeB','tPeakB','tReCircB','tBaselineB', ...
                         'v0B','vUpslopeB','vPeakB','vReCircB'}));
gone=unique([mine,old]);
if isempty(gone), return, end
if istable(T), results.(tbl)=removevars(T,gone); else, results.(tbl)=rmfield(T,gone); end
end


% ====================================================================== %
function ppx=perPixel(cube,mask,L,F,s)
%perPixel  A picture of every requested marker.  THE ONE PATH THAT READS SOURCE.
%
%   Every pixel goes through the core the segments did, on the same clock and with the
%   same complete layout, so a map and a scalar cannot mean two different things.  A
%   BACKGROUND PIXEL IS WRITTEN BY NO ITERATION AND STAYS NaN, confidence included.
%
%   IT IS AFFORDABLE WITHOUT THE POOL, and that is what deleting the per-frame spatial
%   median bought: 0.19 ms per unit is 5.8 min for 1.68 M pixels serially and 0.4 min on
%   sixteen workers.  s.parforCTTHPixels is therefore a worker BOUND on the same loop
%   body rather than the difference between feasible and not - parfor(...,0) runs the
%   identical body in the client and starts no pool.
[Y,X,nT]=size(cube);
D=reshape(cube,Y*X,nT);
maskLin=mask(:);
names=L.blockNames;
sPix=s; sPix.segCtthReturn=s.ppxCtthReturn;
LTree=getBolusMetrics(L.time,L.inputTrace,sPix);

Mp=nan(Y*X,numel(names),'single');
nw=0; if s.parforCTTHPixels, nw=Inf; end        % a worker BOUND, not a branch
sPar=reportSettings(s);
parfor (p=1:Y*X, nw)
    if ~maskLin(p), continue, end               % background stays NaN
    m=getBolusMetrics(double(D(p,:)).',L,sPar);
    row=nan(1,numel(names),'single');
    for j=1:numel(names), row(j)=m.(names{j}); end
    Mp(p,:)=row;
end

Cp=getBolusConfidence(Mp,names,L,s);
Mp=gateCth(Mp,names,F);
bg=~maskLin;
conf=Cp.conf;       conf(bg)=NaN;
confMin=Cp.confMin; confMin(bg)=NaN;

ppx=struct();
for k=1:numel(LTree.scalarNames)
    nm=LTree.scalarNames{k};
    ppx.(['bs' nm])=reshape(blockCol(Mp,names,nm),Y,X);
end
ppx.bsConf   =reshape(conf,Y,X);
ppx.bsConfMin=reshape(confMin,Y,X);
end


% ====================================================================== %
function v=blockCol(M,names,want)
%blockCol  One marker of the block as a [nUnit x 1] single, indexed BY NAME.
%   The block is always measured with the complete set, so a miss here is a programming
%   error rather than a level the user switched off - and it is named as one.
k=find(strcmp(names,want),1);
if isempty(k)
    error('runCTTH:missingMetric', ...
        ['The marker block has no ''%s''.  The block is measured with the COMPLETE ' ...
         'layout; only the tree is gated.'],want);
end
v=M(:,k);
end


% ====================================================================== %
function name=inputSignal(results,fName)
%inputSignal  The traces the arterial input is derived from.
%   sData, else dvsData.  results.sMap is the segmentation the recording carries and
%   sData is what it indexes, so this is a segmented-trace rule.  dvsDiameter is NEVER
%   it: a diameter is a width in pixels and its "earliest arrival" is when a vessel
%   dilated, which is not what an input function is.
for want={'sData','dvsData'}
    if isfield(results,want{1}) && ~isempty(results.(want{1}))
        name=want{1}; return
    end
end
error('runCTTH:noSegmentTraces', ...
    ['%s carries neither results.sData nor results.dvsData, so the arterial input ' ...
     'cannot be derived.  Run the segmentation first.'],fName);
end


% ====================================================================== %
function requireBolus(results,fName)
%requireBolus  The two things a product must carry before a transit time means anything.
if ~isfield(results,'time') || isempty(results.time)
    error('runCTTH:noTime', ...
        '%s carries no results.time, so there is no clock to measure a transit on.',fName);
end
if ~isfield(results,'cMask') || isempty(results.cMask)
    error('runCTTH:noCMask', ...
        ['%s carries no results.cMask, so nothing says which pixels are tissue.  Run ' ...
         'the segmentation first.'],fName);
end
end


% ====================================================================== %
function p=ctthPrefix(sigName)
%ctthPrefix  The branch prefix of a signal: bs bolus trace, bd bolus diameter.
%   The same grammar as the NVC step's ns / nd, and for the same reason: the physical
%   quantity is in the column name, and a diameter must never be pooled with an intensity.
switch sigName
    case 'dvsDiameter', p='bd';
    otherwise,          p='bs';
end
end


% ====================================================================== %
function s=ctthDefaults(s)
%ctthDefaults  Resolve every default HERE, so the saved settings record what ran.
%   THE CORES DEFAULT THE SAME VALUES INTERNALLY.  They are repeated here rather than
%   left to them because the settings file is the record of what a recording was
%   processed with, and a default that only exists inside a core is not in it.
%
%   THERE IS NO s.calcData.  The step this replaces gated the per-pixel pass on
%   "regions" | "all" AND this session's brief also specifies s.ppxCtthReturn as the
%   selector.  Two fields for one decision can disagree - "all" with an empty selector
%   computes nothing, an empty "regions" with a selector computes nothing either - so
%   the library's own idiom is used instead: ppxCtthReturn is the GATE and the selector
%   at once, exactly as ppxNvcReturn, ppxPulsReturn and ppxVsmReturn are in the three
%   steps that already work this way.  calcData retires with the old wrapper.
if ~isfield(s,'ctthLevelPcts')  ||isempty(s.ctthLevelPcts),  s.ctthLevelPcts=[10 25 50 75 90]; end
if ~isfield(s,'ctthPlateauSec') ||isempty(s.ctthPlateauSec), s.ctthPlateauSec=1;   end
if ~isfield(s,'ctthGuardSec')   ||isempty(s.ctthGuardSec),   s.ctthGuardSec=0.5;   end
if ~isfield(s,'ctthSlopeSec')   ||isempty(s.ctthSlopeSec),   s.ctthSlopeSec=0.1;   end
if ~isfield(s,'ctthBaselineSec'),                            s.ctthBaselineSec=[]; end
if ~isfield(s,'ctthInputAmpPct')||isempty(s.ctthInputAmpPct),s.ctthInputAmpPct=75; end
if ~isfield(s,'ctthInputPct')   ||isempty(s.ctthInputPct),   s.ctthInputPct=5;     end
if ~isfield(s,'ctthConfThreshold')   ||isempty(s.ctthConfThreshold),    s.ctthConfThreshold=0.6;    end
if ~isfield(s,'ctthConfMinThreshold')||isempty(s.ctthConfMinThreshold), s.ctthConfMinThreshold=0.2; end
if ~isfield(s,'ctthSettleFrac') ||isempty(s.ctthSettleFrac), s.ctthSettleFrac=1;   end
if ~isfield(s,'ctthInputScale') ||isempty(s.ctthInputScale), s.ctthInputScale=0.5; end
if ~isfield(s,'ctthCthTol')     ||isempty(s.ctthCthTol),     s.ctthCthTol=0.25;    end
if ~isfield(s,'ctthCthProbes'),                              s.ctthCthProbes=[];   end
if ~isfield(s,'segCtthReturn')  ||isempty(s.segCtthReturn)
    s.segCtthReturn={'levels','amplitudes','times','transit','shape'};
end
%ppxCtthReturn defaults ONLY when the field is ABSENT - an explicit {} must stay empty
%(per-pixel off), exactly like runPulsatility's s.ppxPulsReturn.
if ~isfield(s,'ppxCtthReturn'), s.ppxCtthReturn={}; end
if ~isfield(s,'parforCTTHPixels')||isempty(s.parforCTTHPixels), s.parforCTTHPixels=true; end
s.parforCTTHPixels=logical(s.parforCTTHPixels);

%THE SIGNALS ARE NOT A SETTING.  A bolus product carries whichever traces the
%segmentation wrote, and there is no reason to measure fewer of them: the marker cost is
%0.19 ms per unit and the dvs tables are two orders of magnitude smaller than sData.
s.ctthSignals={'sData','dvsData','dvsDiameter'};
end


% ====================================================================== %
function drawInputPage(rep,results,aif,iInfo,L,inputName)
%drawInputPage  '_rep_ctth-input.jpg' - WHERE THE ARTERIAL CURVE CAME FROM.
%
%   IT MUST SHOW BOTH CURVES, because the distance between them is the thing the reader
%   has to believe.  Measured on the reference recording the whole field is 0.674 s later
%   than the derived input - which is LARGER than the entire arterial-to-parenchyma delay
%   it would be subtracted from - so a reader who cannot see that gap has no way to tell
%   a derived input from a mis-derived one.
%
%   A PAGE IS A BY-PRODUCT AND MUST NEVER KILL A RUN.  The analysis is complete by the
%   time this is attempted, so a failed drawing is swallowed.
fh=[];
try
    fh=reportFigure(rep,'ctth-input');
    tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1);
    img=double(results.imgI);
    imagesc(ax,img)
    clim(ax,[prctile(img(:),1) prctile(img(:),99)])
    colormap(ax,gray)
    axis(ax,'image')
    hold(ax,'on')
    %the chosen units painted on the picture, not listed beside it
    sel=paintUnits(results.sMap,iInfo.units);
    hSel=imagesc(ax,cat(3,ones(size(sel)),0.25*ones(size(sel)),0.15*ones(size(sel))));
    set(hSel,'AlphaData',0.85*double(sel));
    hold(ax,'off')
    title(ax,sprintf('%d of %d vessels make up the arterial curve', ...
        iInfo.nUnits,numel(iInfo.units)))

    ax=nexttile(tl,2);
    t=L.time;
    fieldY=mean(double(results.(inputName)),2,'omitnan');
    fieldY=fieldY-mean(fieldY(L.blIdx));
    plot(ax,t,fieldY,'-','Color',[0.55 0.55 0.60],'LineWidth',1.2)
    hold(ax,'on')
    plot(ax,t,aif,'-','Color',[0.80 0.20 0.15],'LineWidth',1.8)
    oF=getBolusMoments(t,fieldY,L.plateauSec);
    xline(ax,L.in.m1,'-','Color',[0.80 0.20 0.15]);
    xline(ax,oF.m1,'-','Color',[0.55 0.55 0.60]);
    hold(ax,'off')
    grid(ax,'on')
    xlabel(ax,'Time, s'); ylabel(ax,'Intensity above the pre-injection level')
    legend(ax,{'the whole field','the arterial curve'},'Location','southeast','Box','off')
    title(ax,sprintf(['The field fills %.3f s after the artery; %.3f s of every ' ...
        'absolute time is the infusion'],oF.m1-L.in.m1,L.in.m1))

    reportSave(rep,fh,'ctth-input');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end


% ====================================================================== %
function drawMarkersPage(rep,results,A,L)
%drawMarkersPage  '_rep_ctth-markers.jpg' - the curve, its times, and its plateau.
%
%   THE UNRETURNED FRACTION IS ON THIS PAGE FOR THE SAME REASON IT IS ON THE ENTRY
%   STEP'S: it is the one number that says whether a first-pass area would have meant
%   anything, and every refusal in this step follows from it.
fh=[];
try
    a=A{1};
    y=mean(double(results.(a.name)),2,'omitnan');
    Bl=mean(y(L.blIdx)); Fn=mean(y(L.endIdx)); Pk=max(y);
    unret=(Fn-Bl)/max(Pk-Bl,eps);
    t=L.time;

    fh=reportFigure(rep,'ctth-markers');
    tl=tiledlayout(fh,1,3,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1,[1 2]);
    plot(ax,t,y-Bl,'-','Color',[0.20 0.35 0.60],'LineWidth',1.6)
    hold(ax,'on')
    xregion(ax,t(find(L.blIdx,1)),t(find(L.blIdx,1,'last')), ...
        'FaceColor',[0.45 0.45 0.45],'FaceAlpha',0.10);
    xregion(ax,t(find(L.endIdx,1)),t(find(L.endIdx,1,'last')), ...
        'FaceColor',[0.45 0.45 0.45],'FaceAlpha',0.10);
    yline(ax,Fn-Bl,'-','Label','plateau','Color',[0.15 0.45 0.75],'LineWidth',1.4, ...
        'LabelHorizontalAlignment','left');
    [tv,tn]=medianTimes(a);
    for j=1:numel(tv)
        xline(ax,tv(j),'--','Color',[0.85 0.20 0.15],'Label',tn{j}, ...
            'LabelOrientation','horizontal','Interpreter','none');
    end
    hold(ax,'off')
    grid(ax,'on')
    xlabel(ax,'Time, s'); ylabel(ax,'Intensity above the pre-injection level')
    title(ax,sprintf(['%s, the median vessel''s times - %.3f of the rise had not come ' ...
        'back by the end'],a.name,unret))

    ax=nexttile(tl,3);
    axis(ax,'off')
    text(ax,0.0,1.0,markerLines(A,L),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName','Courier New','FontSize',9,'Interpreter','none')
    title(ax,'Median vessel')

    reportSave(rep,fh,'ctth-markers');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end


% ====================================================================== %
function [tv,tn]=medianTimes(a)
%medianTimes  The median unit's T0, level times and TPeak, on the recording clock.
want=[{'T0'},arrayfun(@(p)sprintf('T%d',p),a.L.pcts,'UniformOutput',false),{'TPeak'}];
tv=nan(1,numel(want)); tn=want;
for j=1:numel(want)
    tv(j)=median(double(blockCol(a.M,a.L.blockNames,want{j})),'omitnan');
end
ok=isfinite(tv); tv=tv(ok); tn=tn(ok);
end


% ====================================================================== %
function txt=markerLines(A,L)
%markerLines  The median unit's markers, written for a reader rather than a parser.
%   Every name comes from the layout, so a protocol with three timing levels prints
%   three lines without a word changing here.
lines={};
for k=1:numel(A)
    a=A{k};
    if k>1, lines=[lines,{''}]; end %#ok<AGROW>
    lines=[lines,{sprintf('%s  (%d)',a.name,a.nUnit)}]; %#ok<AGROW>
    want=[{'Step','PkRel','BvRel','T0'}, ...
          arrayfun(@(p)sprintf('T%d',p),L.pcts,'UniformOutput',false), ...
          {'TPeak','Delay','Mtt','Cth','EdgeWid','RiseUp','FallDn'}];
    for j=1:numel(want)
        v=double(blockCol(a.M,a.L.blockNames,want{j}));
        lines=[lines,{sprintf('  %-9s %9.3f',want{j},median(v,'omitnan'))}]; %#ok<AGROW>
    end
    lines=[lines,{sprintf('  %-9s %9.3f','Conf',median(double(a.C.conf),'omitnan'))}, ...
        {sprintf('  %-9s %9.0f %%','trusted',100*mean(a.C.trust))}]; %#ok<AGROW>
end
txt=strjoin(lines,newline);
end


% ====================================================================== %
function drawConfidencePage(rep,A)
%drawConfidencePage  '_rep_ctth-confidence.jpg' - the seven checks, and which one bound.
%
%   The distributions say how each check scored; the bar beside them says which check was
%   the WEAKEST for what share of the units, which is the difference between "these
%   vessels are poor" and "these vessels are poor BECAUSE the curve had not settled".
%
%   NO FACTOR IDS ON THE PAGE.  getBolusConfidence carries its own phrasebook, beside the
%   loop that builds the factors, precisely so this page never holds a second list that
%   goes stale the first time a check is added.
fh=[];
try
    a=A{1};
    names=a.C.factorNames;
    words=a.C.factorWords;
    nF=numel(names);
    Fm=zeros(a.nUnit,nF);
    for k=1:nF, Fm(:,k)=double(a.C.factors.(names{k})); end

    fh=reportFigure(rep,'ctth-confidence');
    tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1);
    %ONE ROW PER CHECK, drawn as the middle half with the median on it: seven box-like
    %summaries read at a glance where seven histograms do not.
    hold(ax,'on')
    q=prctile(Fm,[5 25 50 75 95],1);
    for k=1:nF
        plot(ax,[q(1,k) q(5,k)],[k k],'-','Color',[0.65 0.65 0.70],'LineWidth',1);
        plot(ax,[q(2,k) q(4,k)],[k k],'-','Color',[0.30 0.55 0.85],'LineWidth',6);
        plot(ax,q(3,k),k,'o','MarkerFaceColor','w','MarkerEdgeColor',[0.20 0.35 0.60], ...
            'MarkerSize',6,'LineWidth',1.2);
    end
    hold(ax,'off')
    ylim(ax,[0.4 nF+0.6]); xlim(ax,[-0.02 1.02]); grid(ax,'on')
    yticks(ax,1:nF); yticklabels(ax,words);
    xlabel(ax,'Score, 0 to 1')
    title(ax,sprintf('Every check, over %d vessels  (%.1f %% trusted)', ...
        a.nUnit,100*mean(a.C.trust)))

    ax=nexttile(tl,2);
    [~,wk]=min(Fm,[],2);
    share=zeros(1,nF);
    for k=1:nF, share(k)=mean(wk==k); end
    b=barh(ax,1:nF,share);
    b.FaceColor=[0.30 0.55 0.85]; b.EdgeColor='none';
    ylim(ax,[0.4 nF+0.6]); xlim(ax,[0 1]); grid(ax,'on')
    %the rows line up with the panel beside it, so a second copy of the same seven
    %sentences would only take room away from the bars
    yticks(ax,1:nF); yticklabels(ax,[]);
    xlabel(ax,'Share of vessels it was the weakest for')
    title(ax,'Which check bound')

    reportSave(rep,fh,'ctth-confidence');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end


% ====================================================================== %
function drawMapsPage(rep,results,A,F)
%drawMapsPage  '_rep_ctth-maps.jpg' - the five markers as pictures.
%
%   THE Cth PANEL CARRIES ITS FLOOR IN ITS TITLE, per R4-3.  A Cth map without the floor
%   it had to clear is exactly the picture-without-its-null this package keeps refusing:
%   the number is plausible whether or not the recording could resolve it, and only the
%   floor says which.
%
%   IT IS PAINTED FROM THE SEGMENTS unless the per-pixel maps were asked for, in which
%   case those are drawn instead.  Either way the page exists on every run, because a
%   step whose picture is optional is a step nobody looks at.
fh=[];
%the maps are painted through results.sMap, which indexes sData - so sData is the signal
%that can be painted, and it is the one drawn whenever it is present
k=find(strcmp(cellfun(@(x)x.name,A,'UniformOutput',false),'sData'),1);
if isempty(k), k=1; end
a=A{k};
try
    want={'Delay','Mtt','Cth','BvRel'};
    ttl={'Delay, s','Mean transit time, s','Spread of transit times, s', ...
         'Relative blood volume'};
    havePpx=isfield(results,'ctth') && isfield(results.ctth,'ppx');

    %TWO ROWS OF THREE, not one row of five.  These panels are square images on a wide
    %page, so a single row leaves two thirds of the page white and shrinks every map to
    %a thumbnail - and a map nobody can read is not a report page.
    fh=reportFigure(rep,'ctth-maps');
    tl=tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');
    for k=1:numel(want)
        ax=nexttile(tl,k);
        if havePpx && isfield(results.ctth.ppx,['bs' want{k}])
            img=double(results.ctth.ppx.(['bs' want{k}]));
        else
            img=paintValues(results.sMap,double(blockCol(a.M,a.L.blockNames,want{k})));
        end
        drawMap(ax,img)
        t=ttl{k};
        if strcmp(want{k},'Cth')
            if F.usable
                t=sprintf('%s\n(floor %.1f s - below it, blank)',t,F.floor);
            else
                t=sprintf('%s\n(this recording resolves none)',t);
            end
        end
        title(ax,t)
    end
    ax=nexttile(tl,5);
    if havePpx && isfield(results.ctth.ppx,'bsConf')
        img=double(results.ctth.ppx.bsConf);
    else
        img=paintValues(results.sMap,double(a.C.conf));
    end
    drawMap(ax,img); clim(ax,[0 1])
    title(ax,'Confidence')

    reportSave(rep,fh,'ctth-maps');
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end


% ====================================================================== %
function drawMap(ax,img)
%drawMap  One marker picture, scaled to its own middle 98 % so one outlier cannot flatten
%   it, with the unmeasured pixels left white rather than painted as a value.
imagesc(ax,img,'AlphaData',~isnan(img))
v=img(isfinite(img));
if ~isempty(v) && max(v)>min(v)
    lo=prctile(v,1); hi=prctile(v,99);
    if hi>lo, clim(ax,[lo hi]); end
end
axis(ax,'image'); colormap(ax,turbo); colorbar(ax)
xticks(ax,[]); yticks(ax,[])
end


% ====================================================================== %
function img=paintValues(sMap,v)
%paintValues  One value per segment, painted back through the segmentation's label map.
%   Background is NaN, which drawMap leaves blank.
img=nan(size(sMap));
lab=double(sMap);
ok=lab>=1 & lab<=numel(v);
img(ok)=v(lab(ok));
end


% ====================================================================== %
function bw=paintUnits(sMap,units)
%paintUnits  The selected segments as a mask on the picture.
lab=double(sMap);
bw=false(size(lab));
ok=lab>=1 & lab<=numel(units);
bw(ok)=units(lab(ok));
end


% ====================================================================== %
function n=shortName(f)
%shortName  Name and extension - what the user sees on disk.
[~,b,e]=fileparts(char(f));
n=[b e];
end
%------------- END OF CODE --------------
