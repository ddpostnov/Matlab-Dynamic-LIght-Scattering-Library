%runVasomotion  Wavelet analysis of vasomotion amplitude and its dynamics
%
%   runVasomotion(s,fNames) loads every *_BFI_r.mat file in fNames (the
%   blood-flow-index datasets produced by runBFI) and characterises the
%   low-frequency oscillations of the signal with a continuous wavelet
%   transform (analytic Morlet, 'amor').  For each spatial ROI (results
%   fields sData, dvsData, dvsDiameter and, if present, the full-temporal-
%   resolution gsData from runGuidedContrast/runGuidedIntensity - analysed on
%   its own time base results.gsTime) the relative fluctuation
%
%        d(t) = (x(t) - c) / c ,   c = (moving) mean or median of x
%
%   (the centring statistic c is chosen by s.normalisation: global 'mean' or
%   'median', or moving 'mmean' or 'mmedian' with window s.normsize; median
%   by default) is transformed and its amplitude integrated over the vasomotion
%   band s.vFR (VB) and a comparison band s.cFR (CB) to yield instantaneous
%   amplitude envelopes.  From those the shared core getVasomotionMetrics derives
%   per-band envelope moments and percentiles, VB frequency-shape/multiplicity
%   metrics (centroid, spread, peak count), per-frequency amplitude-spectrum
%   moments, VB amplitude-percentile-resolved spectra, independent per-band Otsu
%   flare/silence clustering, a band-limited reconstruction and the decimated
%   amplitude/phase spectrum.  The results are assembled into the band-branched tree
%   RESULTS.vasomotion (see OUTPUTS), and a handful of summary columns are also
%   duplicated into sMetrics/dvsMetrics.
%
%   Which analysis levels are computed and stored is selected by s.segVsmReturn
%   (bands/moments/series/clustering/reconstruction/spectrum; see below) - by
%   default all six.  Which data-type signals are analysed is selected by
%   s.vsmSignals (subset of {sData,dvsData,dvsDiameter,gsData}; default all
%   present).  When s.ppxVsmReturn is non-empty a LEAN per-pixel analysis is also
%   run TRULY per pixel: every pixel is analysed on its own and stored in
%   RESULTS.vasomotion.ppx, with no spatial averaging.  The per-pixel path is
%   deliberately minimal - only the cheap band-amplitude scalars ([Y x X] maps) and
%   an optional decimated amp/phase spectrum ([Y x X x nF x nD] cube); it NEVER runs
%   the expensive VB peak count, the time-series, clustering, reconstruction or
%   per-frequency moments (those are the ~99% per-pixel cost or huge cubes).  The
%   levels are selected by s.ppxVsmReturn, restricted to the lean set {bands,spectrum}.
%   Both per-pixel paths print ~1% completion
%   progress.  A TEMPORARY per-segment averaging demo (s.ppxSegmentAveraging,
%   scaffolding to be removed after validation) can additionally combine the
%   per-pixel complex spectra within each segment - coherently and/or
%   incoherently - into RESULTS.vasomotion.ppxs.coherent / .ppxs.incoherent to
%   compare against the segment-first RESULTS.vasomotion.sData; it is implemented
%   as the self-contained, deletable subfunction runSegmentAveragingTEMP (which
%   does not call the shared core - it carries private frozen copies of the maths).
%
%   INPUTS
%     s        parameter struct with fields
%                • vFR              [lo hi] vasomotion band (VB), Hz
%                • cFR              [lo hi] comparison band (CB), Hz
%                • wFR              [lo hi] wavelet frequency limits, Hz
%                • wVPO             wavelet voices per octave
%                • normalisation    centring statistic for the fluctuation:
%                                   global 'mean'/'median' or moving 'mmean'/
%                                   'mmedian' (default 'median')
%                • normsize         window length for 'mmean'/'mmedian'
%                                   (>=3; inf or 0 -> global statistic)
%                • tgtFS            target sampling frequency for the kept
%                                   spectrum, Hz  ([] = no decimation)
%                • pcts             percentile edges (e.g. 0:10:100)
%                • otsuMaxN         max number of Otsu thresholds to test
%                • otsuElbow        elbow threshold for the threshold count
%                • nPeakProm        findpeaks MinPeakProminence for the VB peak
%                                   count N(t), as a fraction of the per-time VB
%                                   band max-min range (default 0.10 = 10%)
%                • vsmSignals       cell subset of {'sData','dvsData',
%                                   'dvsDiameter','gsData'} naming which data-type
%                                   signals to analyse; default (absent) all four
%                                   (those actually present in results are used)
%                • ppxVsmReturn     GATE + selector for the LEAN true per-pixel analysis -
%                                   there is no analysePerPixel flag.  NON-EMPTY runs
%                                   it, EMPTY ([]) skips it; the value is restricted to the
%                                   lean subset {'bands','spectrum'}: 'bands' -> the cheap
%                                   band-amplitude [Y x X] scalar maps (VB & CB ampMean/
%                                   ampStd always, + ampSkew/centroid-group/ampPct; NEVER
%                                   the VB peak count), 'spectrum' -> the decimated amp/
%                                   phase cube [Y x X x nF x nD].  series/clustering/
%                                   reconstruction/moments are rejected (they never run per
%                                   pixel).  Default (absent) {'bands'}.
%                • ppxSegmentAveraging  TEMPORARY scaffolding: cell subset of
%                                   {'coherent','incoherent'} (default [] = off)
%                                   enabling the per-segment averaging demo ->
%                                   RESULTS.vasomotion.ppxs.coherent / .incoherent
%                • segVsmReturn     cell subset of {'bands','moments','series',
%                                   'clustering','reconstruction','spectrum'} naming
%                                   which per-segment analysis levels to compute and
%                                   store; default (absent) is the COMPLETE set (all
%                                   six levels, incl. the amplitude/phase spectrum grid)
%                • parforVasomotionSegments   logical, default true - run the
%                                   PER-SEGMENT loop in parallel
%                • parforVasomotionPixels     logical, default true - run the true
%                                   PER-PIXEL loop in parallel
%                • parforVasomotionAveraging  logical, default true - run the
%                                   per-segment AVERAGING loop in parallel
%                                   Each is a WORKER BOUND on its own parfor, not a
%                                   branch: false runs the identical loop body
%                                   serially in the client and starts no pool.  One
%                                   field per loop, so the three can be set apart.
%     fNames   cell array of *_BFI_r.mat paths.
%                • Optional workbench hooks in s (no-op when absent):
%                  s.stageFcn(stage,detail), s.cancelFcn()->tf.  Cancel is checked
%                  between files, never inside the parfor.
%
%   OUTPUTS  (RESULTS.vasomotion; band-branched tree.  <sig> = sData, dvsData,
%            dvsDiameter or gsData - one row/slice per segment; nSeg segments,
%            nF = numel(.f).  VB = s.vFR band, CB = s.cFR band, held as struct
%            BRANCHES so field names carry no VB/CB suffix.)
%
%     ROOT axes (shared by the results.time-base signals + ppx/ppxs; gsData
%     carries its OWN copy on the results.gsTime base):
%       .f            [nF x 1]   frequency grid, Hz (descending)
%       .timeWT       [nT x 1]   cone-of-influence-trimmed analysed time, s
%       .timeDWT      [nD x 1]   decimated time base for spectrum.amp/phase, s
%       .pctCenters   [1 x nB]   amplitude-percentile bin centres (from s.pcts)
%
%     Per-signal sub-tree RESULTS.vasomotion.<sig> (present per s.vsmSignals),
%     each with four shape branches gated by s.segVsmReturn:
%       .scalars.VB / .scalars.CB   [nSeg x 1] temporal reductions of the band
%           envelope: ampMean/ampStd/ampSkew, ampPct [nSeg x numel(s.pcts)]
%           ('bands'); VB ONLY: fCentMean/fCentStd, fSprdMean/fSprdStd, shapePeak
%           (=fCentMean/fSprdMean), nPeakMean/nPeakStd ('bands'); durFlareMean/Std,
%           durSilenceMean/Std, ampFlareMean/Std, ampSilenceMean/Std ('clustering';
%           CB from its OWN independent Otsu mask)
%       .fVectors   [nSeg x nF] band-agnostic per-frequency moments ampMean/ampStd/
%           ampSkew ('moments'); .fVectors.VB.ampMeanPct [nSeg x nF x nB] ('moments');
%           .fVectors.VB/.CB .ampFlare/.ampSilence [nSeg x nF] ('clustering')
%       .timeVectors.VB   [nT x nSeg] amp/fCent/fSprd/nPeak ('series'), rData
%           ('reconstruction'), maskFlare ('clustering');
%       .timeVectors.CB   [nT x nSeg] amp ('series'), maskFlare ('clustering')
%       .spectrum.amp / .spectrum.phase  [nSeg x nF x nD] amplitude & phase of the
%           decimated CWT ('spectrum'); at full resolution (s.tgtFS empty) amp=|CWT| /
%           phase=angle(CWT), decimated they are averaged separately (amplitude by moving
%           mean, phase by circular mean) - a deliberate new baseline (not abs/angle of a
%           decimated complex grid)
%     gsData additionally carries its own .f/.timeWT/.timeDWT/.pctCenters.
%
%     LEAN per-pixel RESULTS.vasomotion.ppx (non-empty s.ppxVsmReturn), gated by the
%       lean s.ppxVsmReturn, with the leading segment dimension replaced by the [Y x X]
%       pixel grid.  ONLY the band-amplitude scalars ('bands': .scalars.VB/.CB ampMean/
%       ampStd, ampSkew, ampPct [Y x X x npc] and the VB centroid group - NO nPeak) and,
%       if requested, the decimated .spectrum.amp/.phase [Y x X x nF x nD] ('spectrum');
%       NaN outside segmented/finite pixels.  NO fVectors, NO timeVectors.  It REUSES the
%       root axes (no own f/timeWT/timeDWT/pctCenters).
%
%     TEMPORARY per-segment averaging RESULTS.vasomotion.ppxs.coherent /
%       .incoherent (s.ppxSegmentAveraging; scaffolding to be removed): per-segment
%       sub-trees on the root axes built from an averaged MAGNITUDE spectrum, so they
%       carry scalars + fVectors + timeVectors (amp/fCent/fSprd/nPeak/maskFlare) but
%       NO complex.WT and NO timeVectors.VB.rData.
%
%     Metrics tables (summary columns duplicated from the tree for quick access):
%       sData        -> results.sMetrics  : ampMeanVB, ampMeanCB, fCentMean,
%                                            shapePeak, nPeakMean
%       dvsData      -> results.dvsMetrics: ampMeanVB, ampMeanCB, fCentMean,
%                                            shapePeak, nPeakMean
%       dvsDiameter  -> results.dvsMetrics: the SAME set with a _diam suffix
%                                            (ampMeanVB_diam, ...) to avoid collision
%       gsData       -> NO table (there is no gMetrics).
%     Row order of every RESULTS.vasomotion.<sig> array MATCHES the metrics-table
%     segment order (the 1:nSeg parfor preserves it).
%
%   OUTPUT FILES (side-effects)
%       *_BFI_d.mat   SOURCE   - BFI cube (read only for per-pixel / averaging; not modified)
%       *_BFI_r.mat   RESULTS  - vasomotion results in RESULTS.vasomotion.* (see OUTPUTS);
%                                summary columns also appended to sMetrics/dvsMetrics
%       *_BFI_s.mat   SETTINGS - field settings.runVasomotion added
%
%   EXAMPLE
%     s.vFR=[0.05,0.25]; s.cFR=[0.4,0.6]; s.wFR=[0.01,1]; s.wVPO=10;
%     s.tgtFS=1; s.pcts=0:10:100; s.otsuMaxN=5; s.otsuElbow=0.05; s.nPeakProm=0.10;
%     s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'};
%     s.ppxVsmReturn={'bands'};          % LEAN per-pixel maps ON -> RESULTS.vasomotion.ppx ([] = off)
%     s.ppxSegmentAveraging=[];          % TEMPORARY per-segment averaging demo (off)
%     s.segVsmReturn={'bands','moments','series','clustering','reconstruction','spectrum'};
%     files = dir(fullfile(dataRoot,'*_BFI_r.mat'));
%     runVasomotion(s, fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     Core/Vasomotion/getVasomotionMetrics (shared wavelet vasomotion core),
%     MATLAB Wavelet Toolbox (cwtfilterbank/wt/icwt), Image Processing
%     Toolbox (multithresh/bwareaopen), Statistics/Signal Toolboxes
%     (skewness/prctile/findpeaks) and the Parallel Computing Toolbox (parfor over
%     segments, and over pixels for the per-pixel twin); core LSCI library utilities.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Frequency bands
% s.vFR=[0.05,0.25];      % vasomotion band VB [lo hi], Hz
% s.cFR=[0.4,0.6];        % comparison band CB [lo hi], Hz
% s.wFR=[0.01,1];         % wavelet frequency limits [lo hi], Hz
% s.wVPO=10;              % wavelet voices per octave
% %ADJUSTED IF NECESSARY - Normalisation, percentiles and spectrum decimation
% s.normalisation='median'; % global 'mean'/'median' or moving 'mmean'/'mmedian'
% s.normsize=101;         % window for 'mmean'/'mmedian' (>=3; inf or 0 -> global)
% s.tgtFS=1;              % target sampling frequency for the kept spectrum, Hz
% s.pcts=0:10:100;        % percentile edges for percentile-resolved spectra
% %ADJUSTED IF NECESSARY - Flare/silence segmentation (Otsu + elbow) and peak count
% s.otsuMaxN=5;           % max number of Otsu thresholds to test
% s.otsuElbow=0.05;       % elbow threshold for choosing the number of thresholds
% s.nPeakProm=0.10;       % VB peak-count prominence, fraction of per-time band max-min range
% %ADJUSTED IF NECESSARY - Which signals / outputs to retain
% s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'};  % data-type signals to analyse
% s.ppxVsmReturn={'bands'};  % NON-EMPTY = LEAN per-pixel maps ON (RESULTS.vasomotion.ppx),
%                         % [] = off (there is no analysePerPixel flag). Restricted to the
%                         % lean set: 'bands' (cheap band-amplitude scalar maps, NO peak
%                         % count) and/or 'spectrum' (decimated amp/phase cube [Y x X x nF x
%                         % nD]). series/clustering/reconstruction/moments never run per pixel.
% s.ppxSegmentAveraging=[]; % TEMPORARY scaffolding: {'coherent','incoherent'} or [] (off)
% %s.segVsmReturn selects which analysis levels to compute/store: 'bands' (scalars.VB/CB),
% %'moments' (fVectors ampMean/Std/Skew + VB.ampMeanPct), 'series' (timeVectors
% %amp/fCent/fSprd/nPeak), 'clustering' (flare/silence scalars+spectra+maskFlare),
% %'reconstruction' (timeVectors.VB.rData), 'spectrum' (spectrum.amp/phase). Default = all six.
% s.segVsmReturn={'bands','moments','series','clustering','reconstruction','spectrum'};
% %ADJUSTED IF NECESSARY - Parallel execution (all three default to true when absent)
% s.parforVasomotionSegments=true;   % per-segment loop
% s.parforVasomotionPixels=true;     % true per-pixel loop
% s.parforVasomotionAveraging=true;  % per-segment averaging loop (TEMPORARY block)
% %false runs that loop one item at a time in this MATLAB and starts no parallel pool -
% %slower, but it costs no worker processes, which is what a small machine needs.

function runVasomotion(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_BFI_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_BFI_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_t_BFI_r.mat'').']);
end

%centring statistic used to normalise each signal into the relative fluctuation
%(x-c)/c is selected by s.normalisation/normsize and applied PER SERIES inside
%getVasomotionMetrics (per-column, so bit-identical to the old whole-matrix
%normalisation).  The defaults are set here so they are recorded in the saved
%settings.  'mean'/'median' use a single global value per signal; 'mmean'/
%'mmedian' use a moving statistic with a centred window of length s.normsize
%(inf or 0 -> global; intermediate values must be >=3).
if ~isfield(s,'normalisation') || isempty(s.normalisation)
    s.normalisation='median';
end
if ~isfield(s,'normsize') || isempty(s.normsize)
    s.normsize=inf;
end

%s.vsmSignals selects which data-type signals are analysed (subset of the four
%signal keys).  Default (absent/empty) is the complete set; signals not actually
%present in the loaded results are silently skipped by the loop below.
allSignals={'sData','dvsData','dvsDiameter','gsData'};
if ~isfield(s,'vsmSignals') || isempty(s.vsmSignals)
    s.vsmSignals=allSignals;
end
if ischar(s.vsmSignals) || isstring(s.vsmSignals), s.vsmSignals=cellstr(s.vsmSignals); end
if ~iscellstr(s.vsmSignals) || ~all(ismember(s.vsmSignals,allSignals))
    error('s.vsmSignals must be a cell subset of {''sData'',''dvsData'',''dvsDiameter'',''gsData''}.');
end

%TEMPORARY SCAFFOLDING (to be removed after validation): s.ppxSegmentAveraging
%opts into the per-segment averaging demo - a cell subset of {'coherent',
%'incoherent'} (default [] = off).  'coherent' combines the per-pixel COMPLEX
%wavelet coefficients within a segment (via runSegmentation's per-segment sStat,
%recovered from the settings file) and stores the resulting segment metrics in
%RESULTS.vasomotion.ppxs.coherent (the historical transform-then-combine product);
%'incoherent' averages the per-pixel amplitudes (abs first) and stores the parallel
%sibling RESULTS.vasomotion.ppxs.incoherent.  Once this block is deleted, the
%per-pixel analysis (non-empty s.ppxVsmReturn) is PURELY per-pixel (the ppx tree
%below) with no averaging.
if ~isfield(s,'ppxSegmentAveraging'), s.ppxSegmentAveraging=[]; end
if ~isempty(s.ppxSegmentAveraging)
    if ischar(s.ppxSegmentAveraging) || isstring(s.ppxSegmentAveraging)
        s.ppxSegmentAveraging=cellstr(s.ppxSegmentAveraging);
    end
    if ~iscellstr(s.ppxSegmentAveraging) || ...
            ~all(ismember(s.ppxSegmentAveraging,{'coherent','incoherent'}))
        error('s.ppxSegmentAveraging must be a cell subset of {''coherent'',''incoherent''} (or [] = off).');
    end
end

%the LEAN per-pixel maps are gated by s.ppxVsmReturn: NON-EMPTY turns the true
%per-pixel analysis on (there is no separate analysePerPixel flag), EMPTY turns it off.
%s.ppxVsmReturn is kept SEPARATE from the per-segment s.segVsmReturn and is RESTRICTED to
%the lean set {'bands','spectrum'} (validated in the per-pixel block below): the per-pixel
%path never runs the expensive peak count or the series/clustering/reconstruction/moments
%levels.  Its default (absent) is {'bands'} = the cheap band-amplitude [Y x X] scalar maps.
%Absent-only default: an explicit [] must stay empty (off), so do NOT default on isempty.
if ~isfield(s,'ppxVsmReturn'), s.ppxVsmReturn={'bands'}; end

%PARALLELISM IS OPTIONAL, ONE FIELD PER LOOP.  This wrapper has three genuinely
%different parfors - over segments, over pixels, and over the pixels of one segment
%inside the averaging block - so it carries three switches rather than one.  Each is
%a BOUND on its own parfor (Inf workers or 0), never a branch: parfor(...,0) runs the
%identical loop body serially IN THE CLIENT and starts no pool at all, so not one
%line of the loop bodies is duplicated to support the serial case.  Default true, so
%an existing settings file that carries none of these fields behaves exactly as
%before; false is the escape hatch for a machine that cannot afford a 16-worker pool.
if ~isfield(s,'parforVasomotionSegments')  || isempty(s.parforVasomotionSegments)
    s.parforVasomotionSegments=true;
end
if ~isfield(s,'parforVasomotionPixels')    || isempty(s.parforVasomotionPixels)
    s.parforVasomotionPixels=true;
end
if ~isfield(s,'parforVasomotionAveraging') || isempty(s.parforVasomotionAveraging)
    s.parforVasomotionAveraging=true;
end

%which analysis levels to compute/store is selected by s.segVsmReturn (a cell subset
%of {'bands','moments','series','clustering','reconstruction','spectrum'}); resolved
%once per time base into layout.want by getVasomotionMetrics (s.segVsmReturn is the sole
%selector; absent -> the documented default, all six).

%vasomotion and comparison bands must lie within the wavelet frequency limits,
%otherwise the interpolation onto the wavelet grid would silently yield NaN
if s.vFR(1)<s.wFR(1) || s.vFR(2)>s.wFR(2) || s.cFR(1)<s.wFR(1) || s.cFR(2)>s.wFR(2)
    error('s.vFR and s.cFR must lie within s.wFR.');
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  Cancel is only
% checked between files - the per-pixel parfor below runs silently to the end.
rep=reportOpen(s,'Vasomotion',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        if ~isempty(s.ppxVsmReturn) || ~isempty(s.ppxSegmentAveraging)
            load(getProductPath(s.fName,'d'),'source')
        end
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        %shared wavelet setup + ROOT axes for the results.time base (built once and
        %reused for every results.time-base signal + the per-pixel/averaging paths).
        %gsData lives on results.gsTime and gets its own layout/axes inside the loop.
        layoutT=getVasomotionMetrics(results.time,s);
        axT=vsmAxes(layoutT,s);
        results.vasomotion.f=axT.f;
        results.vasomotion.timeWT=axT.timeWT;
        results.vasomotion.timeDWT=axT.timeDWT;
        results.vasomotion.pctCenters=axT.pctCenters;

        for k=1:numel(s.vsmSignals)
            sigName=s.vsmSignals{k};
            if ~isfield(results,sigName), continue; end
            %gsData is stored at the full frame rate (results.gsTime); the other traces
            %share the decimated results.time.  Pick the matching layout/axes so each is
            %analysed at its own sampling frequency.
            if strcmp(sigName,'gsData')
                if ~isfield(results,'gsTime')
                    warning('results.gsData present but results.gsTime is missing; skipping gsData.');
                    continue
                end
                layout=getVasomotionMetrics(results.gsTime,s);
                ax=vsmAxes(layout,s);
            else
                layout=layoutT; ax=axT;
            end
            want=layout.want;   %which analysis levels s.segVsmReturn selects
            if want.spectrum && ~isempty(s.tgtFS) && s.tgtFS>layout.fs
                error('Target sampling frequency cannot be higher than source sampling frequency')
            end
            f=layout.f; timeVSM=layout.timeVSM;
            nF=numel(f); nT=numel(timeVSM); nD=numel(ax.timeDWT);
            npc=numel(s.pcts); nB=npc-1;

            %pull the signal matrix into a plain local: indexing the struct field
            %results.(sigName)(:,i) inside parfor would broadcast the whole results struct.
            sig=results.(sigName);
            if numel(sig)==0, continue; end
            nSeg=size(sig,2);

            %per-segment accumulators (leading dim = segment; timeVectors [nT x nSeg]).
            %nSeg is small, so preallocate every level unconditionally and let the tree
            %assembler drop the unselected branches - the parfor writes are still gated.
            %Zero-default (invalid segments keep 0): band amp mean/std and the
            %flare/silence duration+amplitude clustering scalars.
            [vbAmpMean,vbAmpStd,cbAmpMean,cbAmpStd, ...
             vbDurFlareMean,vbDurFlareStd,vbDurSilenceMean,vbDurSilenceStd, ...
             vbAmpFlareMean,vbAmpFlareStd,vbAmpSilenceMean,vbAmpSilenceStd, ...
             cbDurFlareMean,cbDurFlareStd,cbDurSilenceMean,cbDurSilenceStd, ...
             cbAmpFlareMean,cbAmpFlareStd,cbAmpSilenceMean,cbAmpSilenceStd]=deal(zeros(nSeg,1,'single'));
            %NaN-default (undefined for a non-finite series): skewness, percentiles and
            %the VB frequency-shape/multiplicity scalars (matches the old marker layout).
            [vbAmpSkew,cbAmpSkew,fCentMean,fCentStd,fSprdMean,fSprdStd, ...
             shapePeak,nPeakMean,nPeakStd]=deal(nan(nSeg,1,'single'));
            vbAmpPct=nan(nSeg,npc,'single'); cbAmpPct=nan(nSeg,npc,'single');
            [spctMean,spctStd,spctSkew,vbSpctFlare,vbSpctSilence, ...
             cbSpctFlare,cbSpctSilence]=deal(zeros(nSeg,nF,'single'));
            vbSpctPct=zeros(nSeg,nF,nB,'single');
            [tsVB,tsCB,fCentTV,fSprdTV,nPeakTV,flareVB,flareCB]=deal(zeros(nT,nSeg,'single'));
            rData=zeros(nT,nSeg,class(sig));
            [sAmp,sPhase]=deal(zeros(nSeg,nF,nD,'single'));   %amplitude & phase grids (was one complex WT)

            %broadcast level flags for the parfor (collect only what is selected).  Band
            %scalars are gated per FINE token (bandsAmp/Skew/Shape/Pct/Peak) mirroring the
            %core, so each sub-group's accumulators only fill when that token is requested.
            wantBandsAmp=want.bandsAmp; wantBandsSkew=want.bandsSkew; wantBandsShape=want.bandsShape;
            wantBandsPct=want.bandsPct; wantBandsPeak=want.bandsPeak;
            wantMoments=want.moments; wantSeries=want.series;
            wantClustering=want.clustering; wantRecon=want.reconstruction; keepSpectrum=want.spectrum;

            nwSeg=0; if s.parforVasomotionSegments, nwSeg=Inf; end   %worker bound, not a branch
            parfor (i=1:nSeg, nwSeg)
                m=getVasomotionMetrics(sig(:,i),layout,s);
                if wantRecon
                    rData(:,i)=m.rData;   %already rescaled per series (NaN for a non-finite series, as before)
                end
                if m.valid
                    if wantMoments
                        spctMean(i,:)=m.spctMean; spctStd(i,:)=m.spctStd; spctSkew(i,:)=m.spctSkew;
                        vbSpctPct(i,:,:)=m.vbSpctPct;
                    end
                    if wantBandsAmp
                        vbAmpMean(i)=m.vbAmpMean; vbAmpStd(i)=m.vbAmpStd;
                        cbAmpMean(i)=m.cbAmpMean; cbAmpStd(i)=m.cbAmpStd;
                    end
                    if wantBandsSkew
                        vbAmpSkew(i)=m.vbAmpSkew; cbAmpSkew(i)=m.cbAmpSkew;
                    end
                    if wantBandsPct
                        vbAmpPct(i,:)=m.vbAmpPct; cbAmpPct(i,:)=m.cbAmpPct;
                    end
                    if wantBandsShape
                        fCentMean(i)=m.fCentMean; fCentStd(i)=m.fCentStd; fSprdMean(i)=m.fSprdMean; fSprdStd(i)=m.fSprdStd;
                        shapePeak(i)=m.shapePeak;
                    end
                    if wantBandsPeak
                        nPeakMean(i)=m.nPeakMean; nPeakStd(i)=m.nPeakStd;
                    end
                    if wantSeries
                        tsVB(:,i)=m.ts; tsCB(:,i)=m.tsc; fCentTV(:,i)=m.fCent; fSprdTV(:,i)=m.fSprd; nPeakTV(:,i)=m.nPeak;
                    end
                    if wantClustering
                        vbSpctFlare(i,:)=m.vbSpctFlare; vbSpctSilence(i,:)=m.vbSpctSilence;
                        cbSpctFlare(i,:)=m.cbSpctFlare; cbSpctSilence(i,:)=m.cbSpctSilence;
                        vbDurFlareMean(i)=m.vbDurFlareMean; vbDurFlareStd(i)=m.vbDurFlareStd;
                        vbDurSilenceMean(i)=m.vbDurSilenceMean; vbDurSilenceStd(i)=m.vbDurSilenceStd;
                        vbAmpFlareMean(i)=m.vbAmpFlareMean; vbAmpFlareStd(i)=m.vbAmpFlareStd;
                        vbAmpSilenceMean(i)=m.vbAmpSilenceMean; vbAmpSilenceStd(i)=m.vbAmpSilenceStd;
                        cbDurFlareMean(i)=m.cbDurFlareMean; cbDurFlareStd(i)=m.cbDurFlareStd;
                        cbDurSilenceMean(i)=m.cbDurSilenceMean; cbDurSilenceStd(i)=m.cbDurSilenceStd;
                        cbAmpFlareMean(i)=m.cbAmpFlareMean; cbAmpFlareStd(i)=m.cbAmpFlareStd;
                        cbAmpSilenceMean(i)=m.cbAmpSilenceMean; cbAmpSilenceStd(i)=m.cbAmpSilenceStd;
                        flareVB(:,i)=m.flare; flareCB(:,i)=m.cbFlare;
                    end
                    if keepSpectrum
                        sAmp(i,:,:)=m.amp; sPhase(i,:,:)=m.phase;
                    end
                end
            end

            %pack the accumulators and assemble the band-branched sub-tree.  vsmTree gates
            %every branch by `want`, so unselected levels are simply not stored.
            acc=struct('vbAmpMean',vbAmpMean,'vbAmpStd',vbAmpStd,'vbAmpSkew',vbAmpSkew,'vbAmpPct',vbAmpPct, ...
                'cbAmpMean',cbAmpMean,'cbAmpStd',cbAmpStd,'cbAmpSkew',cbAmpSkew,'cbAmpPct',cbAmpPct, ...
                'fCentMean',fCentMean,'fCentStd',fCentStd,'fSprdMean',fSprdMean,'fSprdStd',fSprdStd, ...
                'shapePeak',shapePeak,'nPeakMean',nPeakMean,'nPeakStd',nPeakStd, ...
                'vbDurFlareMean',vbDurFlareMean,'vbDurFlareStd',vbDurFlareStd, ...
                'vbDurSilenceMean',vbDurSilenceMean,'vbDurSilenceStd',vbDurSilenceStd, ...
                'vbAmpFlareMean',vbAmpFlareMean,'vbAmpFlareStd',vbAmpFlareStd, ...
                'vbAmpSilenceMean',vbAmpSilenceMean,'vbAmpSilenceStd',vbAmpSilenceStd, ...
                'cbDurFlareMean',cbDurFlareMean,'cbDurFlareStd',cbDurFlareStd, ...
                'cbDurSilenceMean',cbDurSilenceMean,'cbDurSilenceStd',cbDurSilenceStd, ...
                'cbAmpFlareMean',cbAmpFlareMean,'cbAmpFlareStd',cbAmpFlareStd, ...
                'cbAmpSilenceMean',cbAmpSilenceMean,'cbAmpSilenceStd',cbAmpSilenceStd, ...
                'spctMean',spctMean,'spctStd',spctStd,'spctSkew',spctSkew,'vbSpctPct',vbSpctPct, ...
                'vbSpctFlare',vbSpctFlare,'vbSpctSilence',vbSpctSilence, ...
                'cbSpctFlare',cbSpctFlare,'cbSpctSilence',cbSpctSilence, ...
                'tsVB',tsVB,'tsCB',tsCB,'fCentTV',fCentTV,'fSprdTV',fSprdTV,'nPeakTV',nPeakTV, ...
                'rData',rData,'flareVB',flareVB,'flareCB',flareCB,'amp',sAmp,'phase',sPhase);
            shp=identityShapes();
            T=assembleVasomotionTree(acc,want,shp);
            if strcmp(sigName,'gsData')   %gsData carries its OWN axes (results.gsTime base)
                T.f=ax.f; T.timeWT=ax.timeWT; T.timeDWT=ax.timeDWT; T.pctCenters=ax.pctCenters;
            end
            results.vasomotion.(sigName)=T;

            %summary columns duplicated into the metrics tables (DATA-MODEL section 11):
            %sData -> sMetrics; dvsData/dvsDiameter -> dvsMetrics (diameter set _diam-suffixed
            %to avoid a collision); gsData writes no table.  Requires a band scalar level; the
            %columns read the always-preallocated accumulators (NaN/0 where a fine token was
            %not requested), so the default all-bands run fills every column.
            if want.anyBand
                switch sigName
                    case 'sData'
                        results.sMetrics.ampMeanVB=vbAmpMean;  results.sMetrics.ampMeanCB=cbAmpMean;
                        results.sMetrics.fCentMean=fCentMean;  results.sMetrics.shapePeak=shapePeak;
                        results.sMetrics.nPeakMean=nPeakMean;
                    case 'dvsData'
                        results.dvsMetrics.ampMeanVB=vbAmpMean;  results.dvsMetrics.ampMeanCB=cbAmpMean;
                        results.dvsMetrics.fCentMean=fCentMean;  results.dvsMetrics.shapePeak=shapePeak;
                        results.dvsMetrics.nPeakMean=nPeakMean;
                    case 'dvsDiameter'
                        results.dvsMetrics.ampMeanVB_diam=vbAmpMean;  results.dvsMetrics.ampMeanCB_diam=cbAmpMean;
                        results.dvsMetrics.fCentMean_diam=fCentMean;  results.dvsMetrics.shapePeak_diam=shapePeak;
                        results.dvsMetrics.nPeakMean_diam=nPeakMean;
                    % gsData: no table (no gMetrics)
                end
            end
        end


        % % For future - delay estimate based on correlation of reconstructed data
        % figure
        % for k=1:1:size(rData,3)
        % delay=zeros(size(rData,2),1);
        % cc=zeros(size(rData,2),1);
        % for i=1:1:size(rData,2)
        % [xc, lags] = xcorr(rData(:,i,k),rData(:,14,k), 400, 'normalized');
        % [cc(i), max_idx] = max(xc);
        % delay(i) = lags(max_idx) / fs;
        % end
        % subplot(1,2,1)
        % imagesc(delay(results.sMap))
        % axis image
        % subplot(1,2,2)
        % imagesc(cc(results.sMap))
        % axis image
        % title(num2str(rFR(k)))
        % pause(0.1)
        % end

        if ~isempty(s.ppxVsmReturn) || ~isempty(s.ppxSegmentAveraging)
            if ~isfield(results,'sMap')
                error('results.sMap not found; runSegmentation must be run before per-pixel analysis.');
            end
            %per-pixel analysis works on source.data (decimated results.time). The shared
            %setup below (pixel timecourses D, segment labels sMapLin, pixel count npx) feeds
            %both the true per-pixel maps (a flat pixel parfor mirroring the per-segment loop)
            %and the TEMPORARY per-segment averaging (which iterates segment by segment,
            %holding only one segment's pixels at a time).
            sz=size(source.data);
            nSeg=double(max(results.sMap(:)));
            D=reshape(source.data,[],sz(3));   %[nPix x T] pixel timecourses
            sMapLin=results.sMap(:);           %[nPix x 1] segment labels (0=background)
            npx=sz(1)*sz(2);

            % =================================================================
            % TRUE per-pixel analysis: run the ANALYSIS per pixel and store the result as
            % [Y x X (x ...)] maps/cubes (NO spatial averaging) in RESULTS.vasomotion.ppx -
            % the per-pixel TWIN of the per-segment sub-tree (pixels in place of segments).
            % =================================================================
            if ~isempty(s.ppxVsmReturn)
                %The true per-pixel path is LEAN (R3): per pixel it computes only WT -> amp/phase
                %-> cheap band-amplitude scalars, NEVER the ~99%-cost VB peak count and never the
                %time-series / clustering / reconstruction / per-frequency-moment cubes.  So
                %s.ppxVsmReturn is restricted to {'bands','spectrum'}: 'bands' maps to the cheap
                %band scalars (VB & CB ampMean/ampStd ALWAYS + ampSkew + the frequency-shape
                %centroid group + ampPct - all trapz/prctile/skewness, measured cost-comparable to
                %the per-pixel wt), EXCLUDING the expensive bandsPeak; 'spectrum' -> the decimated
                %amp/phase cube [Y x X x nF x nD].  Rejected levels error clearly.  Retained
                %expressions are identical to the per-segment block; only the pixel iteration, the
                %[Y x X] reshape (via ppxShapes) and the NaN (not zero) preallocation differ, so
                %non-segmented / non-finite pixels read NaN.  ppx REUSES the root axes - it stores
                %no f/timeWT/timeDWT/pctCenters of its own.
                ppxReq=s.ppxVsmReturn;
                if ischar(ppxReq)||isstring(ppxReq), ppxReq=cellstr(ppxReq); end
                if ~iscellstr(ppxReq)
                    error('s.ppxVsmReturn must be a cell subset of {''bands'',''spectrum''} (the lean per-pixel set).');
                end
                bad=~ismember(ppxReq,{'bands','spectrum'});
                if any(bad)
                    error(['s.ppxVsmReturn is restricted to the LEAN set {''bands'',''spectrum''} for the true ' ...
                        'per-pixel path; invalid level(s): %s.  (series/clustering/reconstruction/moments and the ' ...
                        'VB peak count never run per pixel.)'],strjoin(ppxReq(bad),', '));
                end
                %map the lean tokens to the core's FINE band tokens (bandsPeak deliberately absent)
                bandTok={}; specTok={};
                if ismember('bands',ppxReq),    bandTok={'bandsAmp','bandsSkew','bandsShape','bandsPct'}; end
                if ismember('spectrum',ppxReq), specTok={'spectrum'}; end
                sPix=s; sPix.segVsmReturn=[bandTok,specTok];
                layoutPix=getVasomotionMetrics(results.time,sPix);
                wantPix=layoutPix.want;
                if wantPix.spectrum && ~isempty(s.tgtFS) && s.tgtFS>layoutPix.fs
                    error('Target sampling frequency cannot be higher than source sampling frequency')
                end
                Y=sz(1); X=sz(2); fP=layoutPix.f; nFreqP=numel(fP);
                nTdecP=numel(axT.timeDWT); npcP=numel(s.pcts);

                %lean per-pixel accumulators [npx x ...] preallocated as NaN (background/invalid
                %pixels read NaN), each gated by its FINE band token so an unselected sub-group
                %costs no memory.  Off-groups get an [npx x 1] real placeholder so the parfor
                %always has a first-dim-npx sliced variable to classify; its guarded write never
                %runs, and the tree assembler never reads it (gated by wantPix).
                off=nan(npx,1,'single');
                if wantPix.bandsAmp
                    [pxVbAmpMean,pxVbAmpStd,pxCbAmpMean,pxCbAmpStd]=deal(nan(npx,1,'single'));
                else
                    [pxVbAmpMean,pxVbAmpStd,pxCbAmpMean,pxCbAmpStd]=deal(off);
                end
                if wantPix.bandsSkew
                    [pxVbAmpSkew,pxCbAmpSkew]=deal(nan(npx,1,'single'));
                else
                    [pxVbAmpSkew,pxCbAmpSkew]=deal(off);
                end
                if wantPix.bandsShape
                    [pxFCentMean,pxFCentStd,pxFSprdMean,pxFSprdStd,pxShapePeak]=deal(nan(npx,1,'single'));
                else
                    [pxFCentMean,pxFCentStd,pxFSprdMean,pxFSprdStd,pxShapePeak]=deal(off);
                end
                if wantPix.bandsPct
                    pxVbAmpPct=nan(npx,npcP,'single'); pxCbAmpPct=nan(npx,npcP,'single');
                else
                    [pxVbAmpPct,pxCbAmpPct]=deal(off);
                end
                if wantPix.spectrum
                    pxAmp=nan(npx,nFreqP,nTdecP,'single'); pxPhase=nan(npx,nFreqP,nTdecP,'single');
                else
                    [pxAmp,pxPhase]=deal(off);
                end

                %broadcast level flags for the parfor (fine band tokens; bandsPeak never set)
                wBandsAmp=wantPix.bandsAmp; wBandsSkew=wantPix.bandsSkew; wBandsShape=wantPix.bandsShape;
                wBandsPct=wantPix.bandsPct; wSpectrum=wantPix.spectrum;

                nwPix=0; if s.parforVasomotionPixels, nwPix=Inf; end   %worker bound, not a branch
                parfor (p=1:npx, nwPix)
                    if sMapLin(p)==0, continue; end          %background: leave NaN
                    mp=getVasomotionMetrics(single(D(p,:))',layoutPix,sPix);
                    if mp.valid
                        if wBandsAmp
                            pxVbAmpMean(p)=mp.vbAmpMean; pxVbAmpStd(p)=mp.vbAmpStd;
                            pxCbAmpMean(p)=mp.cbAmpMean; pxCbAmpStd(p)=mp.cbAmpStd;
                        end
                        if wBandsSkew
                            pxVbAmpSkew(p)=mp.vbAmpSkew; pxCbAmpSkew(p)=mp.cbAmpSkew;
                        end
                        if wBandsShape
                            pxFCentMean(p)=mp.fCentMean; pxFCentStd(p)=mp.fCentStd; pxFSprdMean(p)=mp.fSprdMean; pxFSprdStd(p)=mp.fSprdStd;
                            pxShapePeak(p)=mp.shapePeak;
                        end
                        if wBandsPct
                            pxVbAmpPct(p,:)=mp.vbAmpPct; pxCbAmpPct(p,:)=mp.cbAmpPct;
                        end
                        if wSpectrum
                            pxAmp(p,:,:)=mp.amp; pxPhase(p,:,:)=mp.phase;
                        end
                    end
                end

                accP=struct('vbAmpMean',pxVbAmpMean,'vbAmpStd',pxVbAmpStd,'vbAmpSkew',pxVbAmpSkew,'vbAmpPct',pxVbAmpPct, ...
                    'cbAmpMean',pxCbAmpMean,'cbAmpStd',pxCbAmpStd,'cbAmpSkew',pxCbAmpSkew,'cbAmpPct',pxCbAmpPct, ...
                    'fCentMean',pxFCentMean,'fCentStd',pxFCentStd,'fSprdMean',pxFSprdMean,'fSprdStd',pxFSprdStd, ...
                    'shapePeak',pxShapePeak,'amp',pxAmp,'phase',pxPhase);
                shpP=ppxShapes(Y,X,nFreqP,npcP,nTdecP);
                results.vasomotion.ppx=assembleVasomotionTree(accP,wantPix,shpP);
            end

            % =================================================================
            % TEMPORARY per-segment averaging (SCAFFOLDING - remove after validation).
            % Reproduces the historical transform-then-combine product for comparison
            % against the segment-first RESULTS.vasomotion.sData; gated by
            % s.ppxSegmentAveraging.  Now fully self-contained in the deletable
            % subfunction runSegmentAveragingTEMP (it calls getVasomotionMetrics ZERO
            % times, using its own private tmp* copies of the wavelet setup + metric
            % core).  Deleting that subfunction, its tmp* helpers, this call and the
            % s.ppxSegmentAveraging plumbing removes the demo entirely - nothing else
            % moves.  It still builds per-segment sub-trees on the ROOT axes from an
            % averaged MAGNITUDE spectrum (scalars + fVectors + timeVectors, but NO
            % complex.WT and NO timeVectors.VB.rData) via the SHARED tree assembler.
            % =================================================================
            if ~isempty(s.ppxSegmentAveraging)
                results=runSegmentAveragingTEMP(s,settings,D,sMapLin,nSeg,results);
            end
        end

        settings.runVasomotion=reportSettings(s);
        reportWriting(rep);
        save(s.fName,'results','-v7.3','-nocompression');
        save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end

end
reportClose(rep);

end

% =====================================================================
function ax=vsmAxes(layout,s)
%vsmAxes  Shared axes for one time base: frequency grid, analysed time, decimated
%   time (for spectrum.amp/phase) and amplitude-percentile bin centres.  timeDWT uses the
%   SAME even-window decimation as the core, so the axis matches m.amp/m.phase.
ax.f=layout.f;
ax.timeWT=layout.timeVSM;
if ~isempty(s.tgtFS)
    avgN=floor(layout.fs./s.tgtFS/2)*2+1;
    tw=movmean(layout.timeVSM,avgN,'Endpoints','discard'); tw=tw(1:avgN:end);
else
    tw=layout.timeVSM;
end
ax.timeDWT=tw;
ax.pctCenters=(s.pcts(1:end-1)+s.pcts(2:end))./2;
end

% =====================================================================
function shp=identityShapes()
%identityShapes  Reshape handles for the per-segment / per-Y / averaging trees:
%   the accumulators are already in tree shape ([nSeg x ...], timeVectors [nT x nSeg]),
%   so every handle is the identity.
[shp.sc,shp.scPct,shp.fv,shp.pc,shp.tv,shp.spec]=deal(@(A)A);
end

% =====================================================================
function shp=ppxShapes(Y,X,nF,npc,nD)
%ppxShapes  Reshape handles mapping the [npx x ...] LEAN per-pixel accumulators onto the
%   [Y x X x ...] map/cube layout (npx = Y*X pixels are the leading dim, column-major, so
%   reshape places pixel p at (y,x) exactly as sMap does).  The lean ppx tree carries only
%   scalars (+ percentiles) and the optional decimated amp/phase spectrum, so just those
%   three shape classes are provided (no per-frequency / time-vector cubes).
shp.sc   =@(A)reshape(A,Y,X);
shp.scPct=@(A)reshape(A,Y,X,npc);
shp.spec =@(A)reshape(A,Y,X,nF,nD);
end


% #####################################################################
% ## TEMPORARY per-segment averaging cluster (SCAFFOLDING).          ##
% ## Deletes as ONE unit: runSegmentAveragingTEMP + every tmp* helper##
% ## below + comboName + the call site under `if ~isempty(...)` + the##
% ## s.ppxSegmentAveraging plumbing.  It calls getVasomotionMetrics  ##
% ## ZERO times - the tmp* helpers are FROZEN copies of that core's  ##
% ## vsmSetup / makeNrmStat / resolveNPeakProm / vsmSpectrumMetrics / ##
% ## clusterEnvelope / peakCountSeries; do NOT keep them in sync with ##
% ## the core, they die with the block.  It DOES reuse the shared,   ##
% ## non-temporary identityShapes + assembleVasomotionTree (kept).   ##
% #####################################################################
function results=runSegmentAveragingTEMP(s,settings,D,sMapLin,nSeg,results)
%runSegmentAveragingTEMP  TEMPORARY per-segment averaging demo (self-contained).
%   For each segment it transforms every member pixel to its UNdecimated complex CWT
%   (inlined copy of the old core returnComplex step), combines those across the
%   segment (coherent: runSegmentation's sStat then abs; incoherent: mean of |CWT|),
%   and runs the magnitude-spectrum metric core (tmpSpectrumMetrics) on the result,
%   writing RESULTS.vasomotion.ppxs.coherent / .incoherent on the ROOT axes.  The
%   maths is identical to the pre-decouple block, so coh_median stays bit-exact and
%   coh_mean stays within the existing soft rtol.
doCoh=ismember('coherent',s.ppxSegmentAveraging);
doInc=ismember('incoherent',s.ppxSegmentAveraging);
%the coherent per-pixel COMPLEX combine reuses the SAME per-segment reduction statistic
%runSegmentation used to build the original segment traces (results.sData), recovered
%from the settings file (settings.runSegmentation.sStat), so ppxs.coherent stays
%comparable to the segment-first sData instead of depending on a separate parameter.
%'mean' -> complex mean; otherwise (runSegmentation's median default) -> robust
%component-wise (real/imag) complex median.
if doCoh
    if isfield(settings,'runSegmentation') && isfield(settings.runSegmentation,'sStat') ...
            && strcmp(settings.runSegmentation.sStat,'mean')
        ppwtStat=@(x) mean(x,3);
    else
        ppwtStat=@(x) median(real(x),3)+1i*median(imag(x),3);
    end
end
%private wavelet setup + peak-prominence fraction (frozen copies of the core helpers).
%The centring-stat handle is rebuilt per worker inside the loop (as the old core did),
%so no function handle is broadcast into the parfor.  The averaging always computes the
%FULL metric set (bands/spectrum/series/clustering), so the ppxs sub-trees carry
%scalars+fVectors+timeVectors but no complex/rData (those need the per-series complex path).
layout=tmpVsmSetup(results.time,s);
nPeakProm=tmpResolveNPeakProm(s);
wantAvg=struct('bands',true,'spectrum',true,'series',true,'clustering',true, ...
    'reconstruction',false,'complex',false);
fb=layout.fb; fwt=layout.fwt; coi=layout.coi; f=layout.f;
idxsVFR=layout.idxsVFR; idxsCFR=layout.idxsCFR; fs=layout.fs;
keep=coi<0.05;
nT=sum(keep); nf=numel(f);
nwAvg=0; if s.parforVasomotionAveraging, nwAvg=Inf; end   %worker bound, not a branch
shpA=identityShapes();

for combo=1:2
    if combo==1 && ~doCoh, continue; end
    if combo==2 && ~doInc, continue; end
    %accumulators for this combine (metric levels only; no complex/rData).  Same
    %zero/NaN default split as the per-segment block (segments with no valid pixels
    %keep the default).
    [avbAmpMean,avbAmpStd,acbAmpMean,acbAmpStd, ...
     avbDurFlareMean,avbDurFlareStd,avbDurSilenceMean,avbDurSilenceStd, ...
     avbAmpFlareMean,avbAmpFlareStd,avbAmpSilenceMean,avbAmpSilenceStd, ...
     acbDurFlareMean,acbDurFlareStd,acbDurSilenceMean,acbDurSilenceStd, ...
     acbAmpFlareMean,acbAmpFlareStd,acbAmpSilenceMean,acbAmpSilenceStd]=deal(zeros(nSeg,1,'single'));
    [avbAmpSkew,acbAmpSkew,afCentMean,afCentStd,afSprdMean,afSprdStd, ...
     ashapePeak,anPeakMean,anPeakStd]=deal(nan(nSeg,1,'single'));
    avbAmpPct=nan(nSeg,numel(s.pcts),'single'); acbAmpPct=nan(nSeg,numel(s.pcts),'single');
    [aspctMean,aspctStd,aspctSkew,avbSpctFlare,avbSpctSilence, ...
     acbSpctFlare,acbSpctSilence]=deal(zeros(nSeg,nf,'single'));
    avbSpctPct=zeros(nSeg,nf,numel(s.pcts)-1,'single');
    [atsVB,atsCB,afCentTV,afSprdTV,anPeakTV,aflareVB,aflareCB]=deal(zeros(nT,nSeg,'single'));

    for kk=1:nSeg
        pix=find(sMapLin==kk);
        nPix=numel(pix);
        if nPix==0, continue; end
        segData=D(pix,:);              %[nPix x T] sliced input for parfor
        wttsPix=complex(nan(nf,nT,nPix,'single'));
        parfor (pp=1:nPix, nwAvg)
            %inlined copy of the old core returnComplex step: (x-c)/c -> wt -> coi-trim
            %-> interp1 onto the grid f (UNdecimated complex coefficients for the combine)
            nrmStat=tmpMakeNrmStat(s);
            X=single(segData(pp,:))';
            nrm=nrmStat(X,1);
            sub=single((X-nrm)./nrm);
            if all(isfinite(sub))
                wtts=wt(fb,sub);
                wtts=wtts(:,keep);
                wttsPix(:,:,pp)=interp1(fwt,wtts,f);
            end
        end
        validPix=squeeze(all(isfinite(wttsPix),[1 2]));
        if ~any(validPix), continue; end
        if combo==1
            %coherent: combine the COMPLEX coefficients (runSegmentation's sStat), then abs
            wttsSeg=abs(ppwtStat(wttsPix(:,:,validPix)));
        else
            %incoherent: mean of the per-pixel amplitudes (abs first); coherent statistic unused
            wttsSeg=mean(abs(wttsPix(:,:,validPix)),3);
        end
        m=tmpSpectrumMetrics(wttsSeg,f,idxsVFR,idxsCFR,s.vFR,s.cFR,s.pcts,fs,s.otsuMaxN,s.otsuElbow,nPeakProm,wantAvg);
        avbAmpMean(kk)=m.vbAmpMean; avbAmpStd(kk)=m.vbAmpStd; avbAmpSkew(kk)=m.vbAmpSkew; avbAmpPct(kk,:)=m.vbAmpPct;
        acbAmpMean(kk)=m.cbAmpMean; acbAmpStd(kk)=m.cbAmpStd; acbAmpSkew(kk)=m.cbAmpSkew; acbAmpPct(kk,:)=m.cbAmpPct;
        afCentMean(kk)=m.fCentMean; afCentStd(kk)=m.fCentStd; afSprdMean(kk)=m.fSprdMean; afSprdStd(kk)=m.fSprdStd;
        ashapePeak(kk)=m.shapePeak; anPeakMean(kk)=m.nPeakMean; anPeakStd(kk)=m.nPeakStd;
        aspctMean(kk,:)=m.spctMean; aspctStd(kk,:)=m.spctStd; aspctSkew(kk,:)=m.spctSkew; avbSpctPct(kk,:,:)=m.vbSpctPct;
        avbSpctFlare(kk,:)=m.vbSpctFlare; avbSpctSilence(kk,:)=m.vbSpctSilence;
        acbSpctFlare(kk,:)=m.cbSpctFlare; acbSpctSilence(kk,:)=m.cbSpctSilence;
        avbDurFlareMean(kk)=m.vbDurFlareMean; avbDurFlareStd(kk)=m.vbDurFlareStd;
        avbDurSilenceMean(kk)=m.vbDurSilenceMean; avbDurSilenceStd(kk)=m.vbDurSilenceStd;
        avbAmpFlareMean(kk)=m.vbAmpFlareMean; avbAmpFlareStd(kk)=m.vbAmpFlareStd;
        avbAmpSilenceMean(kk)=m.vbAmpSilenceMean; avbAmpSilenceStd(kk)=m.vbAmpSilenceStd;
        acbDurFlareMean(kk)=m.cbDurFlareMean; acbDurFlareStd(kk)=m.cbDurFlareStd;
        acbDurSilenceMean(kk)=m.cbDurSilenceMean; acbDurSilenceStd(kk)=m.cbDurSilenceStd;
        acbAmpFlareMean(kk)=m.cbAmpFlareMean; acbAmpFlareStd(kk)=m.cbAmpFlareStd;
        acbAmpSilenceMean(kk)=m.cbAmpSilenceMean; acbAmpSilenceStd(kk)=m.cbAmpSilenceStd;
        atsVB(:,kk)=m.ts; atsCB(:,kk)=m.tsc; afCentTV(:,kk)=m.fCent; afSprdTV(:,kk)=m.fSprd; anPeakTV(:,kk)=m.nPeak;
        aflareVB(:,kk)=m.flare; aflareCB(:,kk)=m.cbFlare;
    end

    accA=struct('vbAmpMean',avbAmpMean,'vbAmpStd',avbAmpStd,'vbAmpSkew',avbAmpSkew,'vbAmpPct',avbAmpPct, ...
        'cbAmpMean',acbAmpMean,'cbAmpStd',acbAmpStd,'cbAmpSkew',acbAmpSkew,'cbAmpPct',acbAmpPct, ...
        'fCentMean',afCentMean,'fCentStd',afCentStd,'fSprdMean',afSprdMean,'fSprdStd',afSprdStd, ...
        'shapePeak',ashapePeak,'nPeakMean',anPeakMean,'nPeakStd',anPeakStd, ...
        'vbDurFlareMean',avbDurFlareMean,'vbDurFlareStd',avbDurFlareStd, ...
        'vbDurSilenceMean',avbDurSilenceMean,'vbDurSilenceStd',avbDurSilenceStd, ...
        'vbAmpFlareMean',avbAmpFlareMean,'vbAmpFlareStd',avbAmpFlareStd, ...
        'vbAmpSilenceMean',avbAmpSilenceMean,'vbAmpSilenceStd',avbAmpSilenceStd, ...
        'cbDurFlareMean',acbDurFlareMean,'cbDurFlareStd',acbDurFlareStd, ...
        'cbDurSilenceMean',acbDurSilenceMean,'cbDurSilenceStd',acbDurSilenceStd, ...
        'cbAmpFlareMean',acbAmpFlareMean,'cbAmpFlareStd',acbAmpFlareStd, ...
        'cbAmpSilenceMean',acbAmpSilenceMean,'cbAmpSilenceStd',acbAmpSilenceStd, ...
        'spctMean',aspctMean,'spctStd',aspctStd,'spctSkew',aspctSkew,'vbSpctPct',avbSpctPct, ...
        'vbSpctFlare',avbSpctFlare,'vbSpctSilence',avbSpctSilence, ...
        'cbSpctFlare',acbSpctFlare,'cbSpctSilence',acbSpctSilence, ...
        'tsVB',atsVB,'tsCB',atsCB,'fCentTV',afCentTV,'fSprdTV',afSprdTV,'nPeakTV',anPeakTV, ...
        'flareVB',aflareVB,'flareCB',aflareCB);
    results.vasomotion.ppxs.(comboName(combo))=assembleVasomotionTree(accA,wantAvg,shpA);
end
end

% =====================================================================
function layout=tmpVsmSetup(tvec,s)
%tmpVsmSetup  TEMPORARY frozen copy of getVasomotionMetrics/vsmSetup (fb/fwt/coi/f grid).
fs=1./mean(diff(tvec));
fb=cwtfilterbank('SignalLength',numel(tvec), ...
    'SamplingFrequency',fs, ...
    'Wavelet','amor', ...
    'FrequencyLimits',s.wFR, ...
    'VoicesPerOctave',s.wVPO);
[~,fwt,coi]=wt(fb,single(tvec));
f=sort(unique(cat(1,fwt(:),s.vFR',s.cFR')),'descend');
idxsVFR=f>=s.vFR(1) & f<=s.vFR(2);
idxsCFR=f>=s.cFR(1) & f<=s.cFR(2);
timeVSM=tvec(coi<0.05);
layout.fb=fb; layout.fwt=fwt; layout.coi=coi; layout.f=f;
layout.idxsVFR=idxsVFR; layout.idxsCFR=idxsCFR; layout.fs=fs;
layout.timeVSM=timeVSM;
end

% =====================================================================
function nrmStat=tmpMakeNrmStat(s)
%tmpMakeNrmStat  TEMPORARY frozen copy of getVasomotionMetrics/makeNrmStat.
if ~isfield(s,'normalisation') || isempty(s.normalisation)
    s.normalisation='median';
end
if ~isfield(s,'normsize') || isempty(s.normsize)
    s.normsize=inf;
end
switch s.normalisation
    case 'mean'
        nrmStat=@(x,dim) mean(x,dim,'omitmissing');
    case 'median'
        nrmStat=@(x,dim) median(x,dim,'omitmissing');
    case 'mmean'
        if isinf(s.normsize) || s.normsize==0
            nrmStat=@(x,dim) mean(x,dim,'omitmissing');
        elseif s.normsize<3
            error('s.normsize must be >=3 (or inf/0 for a global statistic) when using ''mmean''/''mmedian''.');
        else
            nrmStat=@(x,dim) movmean(x,s.normsize,dim,'omitmissing');
        end
    case 'mmedian'
        if isinf(s.normsize) || s.normsize==0
            nrmStat=@(x,dim) median(x,dim,'omitmissing');
        elseif s.normsize<3
            error('s.normsize must be >=3 (or inf/0 for a global statistic) when using ''mmean''/''mmedian''.');
        else
            nrmStat=@(x,dim) movmedian(x,s.normsize,dim,'omitmissing');
        end
    otherwise
        error('s.normalisation must be ''mean'', ''median'', ''mmean'' or ''mmedian''.');
end
end

% =====================================================================
function p=tmpResolveNPeakProm(s)
%tmpResolveNPeakProm  TEMPORARY frozen copy of getVasomotionMetrics/resolveNPeakProm.
if ~isfield(s,'nPeakProm') || isempty(s.nPeakProm)
    p=0.10;
    return
end
p=s.nPeakProm;
if ~(isnumeric(p)&&isscalar(p)&&isfinite(p)&&p>=0)
    error('getVasomotionMetrics:nPeakProm', ...
        's.nPeakProm must be a finite scalar >= 0 (fraction of the per-time VB band max-min range).');
end
end

% =====================================================================
function m = tmpSpectrumMetrics(wtts,f,idxsVFR,idxsCFR,vFR,cFR,pcts,fs,maxN,otsuElbow,nPeakProm,want)
%tmpSpectrumMetrics  TEMPORARY frozen copy of getVasomotionMetrics/vsmSpectrumMetrics.
m=struct();

needVB    = want.bands || want.series || want.clustering || want.spectrum;  %VB envelope ts
needCB    = want.bands || want.series || want.clustering;                   %CB envelope tsc
needShape = want.bands || want.series;                                      %VB fc/sprd/nPeak

% ---- (spectrum) per-frequency temporal moments of |CWT| ----
if want.spectrum
    m.spctMean=mean(wtts,2)';
    m.spctStd =std(wtts,0,2)';
    m.spctSkew=skewness(wtts,0,2)';
end

% ---- Group-A base: band-integrated amplitude envelopes a_B(t) (trapz form) ----
if needVB
    AbV=wtts(idxsVFR,:); fV=f(idxsVFR);
    ts=-trapz(fV,AbV,1); ts=ts(:)./(vFR(2)-vFR(1));
    if want.bands || want.spectrum, vpcts=prctile(ts,pcts); end     %shared by vbAmpPct + vbSpctPct binning
end
if needCB
    AbC=wtts(idxsCFR,:); fC=f(idxsCFR);
    tsc=-trapz(fC,AbC,1); tsc=tsc(:)./(cFR(2)-cFR(1));
end

% ---- Group-A base: VB instantaneous frequency-shape series (VB band only) ----
if needShape
    den=trapz(fV,AbV,1);                                    %[1 x nT] (<=0 on the descending grid)
    fc=trapz(fV,fV.*AbV,1)./den;                            %centroid; sign cancels
    sprd=sqrt(trapz(fV,(fV-fc).^2.*AbV,1)./den);            %spread (Hz)
    bad=~(-den>0);                                          %band integral <=0 or non-finite -> degenerate
    fc(bad)=NaN; sprd(bad)=NaN;
    fc=fc(:); sprd=sprd(:);
    nPeak=tmpPeakCountSeries(AbV,nPeakProm);                %[nT x 1], NaN where no-peak/degenerate
end

% ---- (series) store the Group-A time series (axis = timeVSM) ----
if want.series
    m.ts=ts;        %VB envelope a_vFR(t)
    m.tsc=tsc;      %CB envelope a_cFR(t)
    m.fCent=fc;     %VB instantaneous centroid
    m.fSprd=sprd;   %VB instantaneous spread
    m.nPeak=nPeak;  %VB instantaneous peak count
end

% ---- (bands) scalar temporal reductions of the Group-A series ----
if want.bands
    m.vbAmpMean=mean(ts);       m.vbAmpStd=std(ts);         %== old adVFR
    m.vbAmpSkew=skewness(ts,0);                             %== old skewVFR
    m.vbAmpPct=vpcts;                                       %== old pctsVFR
    m.cbAmpMean=mean(tsc);      m.cbAmpStd=std(tsc);        %== old adCFR
    m.cbAmpSkew=skewness(tsc,0);                            %== old skewCFR
    m.cbAmpPct=prctile(tsc,pcts);                          %== old pctsCFR
    m.fCentMean=mean(fc,'omitnan');    m.fCentStd=std(fc,0,'omitnan');
    m.fSprdMean=mean(sprd,'omitnan');  m.fSprdStd=std(sprd,0,'omitnan');
    m.shapePeak=m.fCentMean./m.fSprdMean;                  %peak sharpness (how well-defined the rhythm)
    m.nPeakMean=mean(nPeak,'omitnan'); m.nPeakStd=std(nPeak,0,'omitnan');
end

% ---- (spectrum) VB amplitude-percentile-resolved MEAN spectrum ----
if want.spectrum
    nB=numel(pcts)-1;
    pctSpc=zeros(numel(f),nB);
    for ii=1:nB
        sel=ts>=vpcts(ii) & ts<(vpcts(ii+1)+eps);
        pctSpc(:,ii)=mean(wtts(:,sel),2);                  %== old spctPCTS(:,:,1) (mean page only)
    end
    m.vbSpctPct=pctSpc;
end

% ---- (clustering) independent per-band Otsu-elbow flare/silence segmentation ----
if want.clustering
    cV=tmpClusterEnvelope(ts, wtts,fs,vFR,maxN,otsuElbow); %VB (retained: identical to old VB path)
    cC=tmpClusterEnvelope(tsc,wtts,fs,vFR,maxN,otsuElbow); %CB (SAME procedure, independent mask)
    % VB (bit-identical to the pre-redesign VB markers)
    m.flare=cV.mask;
    m.vbSpctFlare=cV.spctFlare;          m.vbSpctSilence=cV.spctSilence;   %== old spctFLR/SLC(:,1)
    m.vbDurFlareMean=cV.durFlare(1);     m.vbDurFlareStd=cV.durFlare(2);   %== old statFLR(2:3)
    m.vbDurSilenceMean=cV.durSilence(1); m.vbDurSilenceStd=cV.durSilence(2); %== old statSLC(2:3)
    m.vbAmpFlareMean=cV.ampFlare(1);     m.vbAmpFlareStd=cV.ampFlare(2);   %== old adFVFR
    m.vbAmpSilenceMean=cV.ampSilence(1); m.vbAmpSilenceStd=cV.ampSilence(2); %== old adSVFR
    % CB (independent CB mask -> CB amplitude inside CB's own clusters; new semantics)
    m.cbFlare=cC.mask;
    m.cbSpctFlare=cC.spctFlare;          m.cbSpctSilence=cC.spctSilence;
    m.cbDurFlareMean=cC.durFlare(1);     m.cbDurFlareStd=cC.durFlare(2);
    m.cbDurSilenceMean=cC.durSilence(1); m.cbDurSilenceStd=cC.durSilence(2);
    m.cbAmpFlareMean=cC.ampFlare(1);     m.cbAmpFlareStd=cC.ampFlare(2);
    m.cbAmpSilenceMean=cC.ampSilence(1); m.cbAmpSilenceStd=cC.ampSilence(2);
end
end

% =====================================================================
function nPeak=tmpPeakCountSeries(Ab,prom)
%tmpPeakCountSeries  TEMPORARY frozen copy of getVasomotionMetrics/peakCountSeries.
nT=size(Ab,2);
nPeak=nan(nT,1);
for tt=1:nT
    col=Ab(:,tt);
    if ~all(isfinite(col)), continue; end         %non-finite column -> NaN
    rng=max(col)-min(col);
    try
        pk=findpeaks(col,'MinPeakProminence',prom*rng);
    catch
        continue                                   %degenerate/errored -> NaN
    end
    if ~isempty(pk), nPeak(tt)=numel(pk); end      %no peaks -> NaN (preset)
end
end

% =====================================================================
function c = tmpClusterEnvelope(env,wtts,fs,vFR,maxN,otsuElbow)
%tmpClusterEnvelope  TEMPORARY frozen copy of getVasomotionMetrics/clusterEnvelope.
nF=size(wtts,1);

% ---- elbow-selected Otsu threshold -> raw flare mask ----
optN=zeros(maxN,1);
for n=1:maxN, [~,optN(n)]=multithresh(env,n); end
optN=diff(optN);
optN=find(optN<otsuElbow,1,'first');
if isempty(optN)
    mask=false(size(env));
else
    minT=ceil(fs./min(vFR)); %minimum flare duration
    thr=multithresh(env,optN);
    mask=env>thr(1);
    tmp=padarray(mask,[minT,0],"replicate");
    tmp=bwareaopen(tmp,minT)>0;
    mask=tmp(minT+1:end-minT);
end

% ---- per-cluster MEAN spectra + duration/amplitude mean+std ----
spctFlare=zeros(nF,1); spctSilence=zeros(nF,1);
durFlare=[NaN,NaN];    durSilence=[NaN,NaN];      %[mean std] (s)
ampFlare=[0,0];        ampSilence=[0,0];          %[mean std]
if any(mask)
    spctFlare=mean(wtts(:,mask),2);
    ampFlare=[mean(env(mask)),std(env(mask))];
end
if any(~mask)
    spctSilence=mean(wtts(:,~mask),2);
    ampSilence=[mean(env(~mask)),std(env(~mask))];
end
if any(mask) && any(~mask)
    toggleIdx=find(diff(mask(:)'));
    internalLengths=diff(toggleIdx);
    internalStates=mask(toggleIdx(1:end-1)+1);
    durFLR=internalLengths(internalStates==1);
    durSLC=internalLengths(internalStates==0);
    durFlare=[mean(durFLR),std(durFLR)]/fs;
    durSilence=[mean(durSLC),std(durSLC)]/fs;
end

c.mask=mask;
c.spctFlare=spctFlare; c.spctSilence=spctSilence;
c.durFlare=durFlare;   c.durSilence=durSilence;
c.ampFlare=ampFlare;   c.ampSilence=ampSilence;
end

% =====================================================================
function nm=comboName(combo)
%comboName  Map the averaging combine index to its ppxs sub-tree field name.
%   TEMPORARY - part of the per-segment averaging cluster (dies with the block).
if combo==1, nm='coherent'; else, nm='incoherent'; end
end
