%runIntensityPulsatility - How much the plasma volume in each vessel pulses over one beat.
%
%   THE FLUORESCENCE PULSATILITY STEP, AND IT IS A SEPARATE WRAPPER RATHER THAN A
%   REGISTRY TWIN [D13].  runVasomotion names the physical quantity in a TREE NODE
%   (results.vasomotion.sData / .dvsData / .dvsDiameter), so one wrapper serves both
%   modalities and a trace is a trace.  runPulsatility names it in the COLUMN PREFIX -
%   its own header says so: "a field-name PREFIX encoding the physical quantity
%   (ps = pulsatile flow, pd = pulsatile diameter)".  Pointed at a '_c_I' product it
%   would write psPI meaning a pulsatile INTENSITY: one column name carrying two
%   physical quantities, in a table that is pooled across recordings and exported to
%   one spreadsheet.  That is the hazard this library already names for ns* against
%   nd*, so the branch gets its own wrapper and its own prefix.
%
%   THE SEPARATION IS AT THE WRAPPER, NOT AT THE SCIENCE.  Every number here comes from
%   the SHARED core Core/Pulsatility/getPulsatilityMetrics, which takes one averaged-
%   cycle waveform and is agnostic about what it measures.  Nothing of its maths is
%   copied and there is no fluorescence fork of it.
%
%   THE PREFIX IS pv, FOR PULSATILE VOLUME.  Fluorescence intensity from an
%   intravascular tracer is proportional to how much labelled plasma sits in the light
%   path, so what these columns describe is a plasma-volume waveform - not a flow
%   (ps*, the speckle branch's) and not a diameter (pd*).  It reads as one family with
%   those two and it says what was measured.  A prefix reaches every settings file,
%   every export column and every explorer axis, and this library does not carry a
%   second spelling of anything, so renaming it later costs a re-processing.
%
%   THE INPUT IS ALREADY ONE AVERAGED CYCLE.  runIntensityInternalCycle folds hundreds
%   of accepted beats onto one phase axis of s.nPhaseBins BIN CENTRES, in ABSOLUTE
%   intensity [D11] - so PI = (max-min)/mean means here exactly what it means on a
%   '_c_BFI' product and the two branches' numbers go in one table.  The bin centres
%   matter to this step: the period is time(end)+time(1) rather than time(end), which
%   getPulsatilityMetrics now owns for both conventions.
%
%   THE MODEL IS OFF BY DEFAULT, AND THAT IS A MEASUREMENT RATHER THAN A PREFERENCE.
%   At about 100 frames a second against a heart of about 10 there are roughly TEN
%   independent samples in a beat, which supports the fundamental and four harmonics
%   and nothing beyond (02-cardiac.md section 3).  A five-harmonic model has ELEVEN
%   free parameters, so on this product it would interpolate the beat rather than
%   describe it and its R2 would be about 1 whatever the data said.  The markers are
%   therefore taken on the averaged cycle itself - which is already a mean over
%   hundreds of real beats, so the smoothing a fit would have provided has been done
%   on measurements instead of on a model, and source.sem says by how much.  Ask for
%   s.segPulsReturn {'markers','model','reconstruction'} and the fit runs.
%
%   THERE IS NO pd* COLUMN ON THIS BRANCH, and it is the one column decision this step
%   takes on its own.  runDynamicSegmentation's per-frame diameter walks INTEGER
%   indices on a profile interpolated by s.pInterpF, so it cannot change by less than
%   0.25 px.  On the fluorescence reference recording the cardiac diameter change is
%   0.014 to 0.152 px zero-to-peak - 1.6 to 18 times SMALLER than one step of that
%   instrument - the per-frame noise is 4.6 to 45 times the modulation, and nothing in
%   any pulsatility wrapper has a matched control to say so.  A pdPI taken there is
%   (max-min)/mean of the quantisation wearing a physiological name.  Where a cardiac
%   diameter change IS measured on this branch is runMotionEnhancement: sub-pixel wall
%   positions, against the scrambled-timing copies stored beside the cycle, with the
%   translation share beside them (wmDilRel, wmDilShare).
%
%   SO A '_c_I' PRODUCT CARRIES THREE FAMILIES OF PER-SEGMENT COLUMN AND THEY MUST NOT
%   BE POOLED:
%       pv*   here            the pulsatile plasma-volume waveform
%       wm*   wall motion     wall displacement, sub-pixel, with a floor and a split
%       ps*   NOWHERE         it is the speckle branch's flow prefix and this step is
%                             the reason it never appears on an intensity product
%
%   IT BUILDS NO CINE AND DETECTS NO HEART RATE.  runMotionEnhancement owns the
%   magnification, the wall columns and the optional movie (master plan section 7,
%   Q1); this is a trace step and stays one - a trace step has nowhere to keep the
%   matched control every motion number needs.  The rate on the report page is READ off
%   the product (results.rateMedian), never detected here: there are three cardiac-gate
%   implementations in this library already and this is not a fourth.
%
% Syntax:
%    runIntensityPulsatility(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%             .nHarm          harmonics in the sine model, when it is asked for
%                             (default 5).
%             .segPulsReturn  per-segment levels, a cell subset of
%                             {'markers','model','reconstruction'}.  Default
%                             {'markers'} - see the header.
%             .ppxPulsReturn  per-pixel GATE: {'markers'} makes the [Y x X] marker
%                             maps, [] skips them.  Default {'markers'}, and it is
%                             restricted to that one level: a per-pixel harmonic fit
%                             on this product is the same over-parameterised model as
%                             above, run a million times.
%             .fitOptions     (optional) fitoptions override, passed to the core.
%    fNames - cell array of '*_c_I_r.mat' paths.
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.  Cancel is checked between files.
%
% Outputs:
%    (none) - writes IN PLACE into the product it was given, NON-DESTRUCTIVELY: the
%    source member is read and never re-saved.
%
%    results.pulsatility          the tree, one sub-tree per analysed signal
%      .time       [nT x 1]       the product's own phase axis, single
%      .harmonics  [nHarm x 1]    column index for scalars.hAmp / hPhase
%      .sData      .scalars pvMin pvMax pvTimeMin pvTimeMax pvMean pvStd pvPI pvRI
%                          pvSymRatio (+ hAmp hPhase pvR2 with 'model')
%                  .timeVectors.fData [nT x nSeg] with 'reconstruction'
%      .dvsData    the same, on runDynamicSegmentation's traces
%      .ppx        .scalars the same nine as [Y x X] maps
%
%    results.sMetrics   (from sData)   gains  pvPI pvRI pvMean pvTimeMin pvTimeMax
%                                             pvSymRatio
%    results.dvsMetrics (from dvsData) gains  the same six
%    <stem>_rep_pulsatility.jpg
%
%    THE TREE HOLDS EVERYTHING AND THE TABLE HOLDS A KEY SET, under the same names -
%    the pattern runPulsatility set and runMotionEnhancement copied.  The tree's rows
%    are in the metrics tables' segment order, so a join is by row.
%
% Example:
%    s.libraryFolder = libraryFolder;
%    s.segPulsReturn = {'markers'};
%    s.ppxPulsReturn = {'markers'};
%    runIntensityPulsatility(s, getFileNamesList(rootFolder,'*_c_I_r.mat'));
%
% Dependencies: Core/Pulsatility/getPulsatilityMetrics; Core/Reporting; MATLAB Curve
%               Fitting Toolbox ONLY when the model level is asked for.
% See also: getPulsatilityMetrics, runPulsatility, runIntensityInternalCycle,
%           runMotionEnhancement, runVasomotion, runSegmentation
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Harmonic model
% s.nHarm=5;              % harmonics in y=SUM a_n*sin(2*pi*n*x+b_n)+c, when the model
%                         % is asked for.  A fluorescence beat carries about ten
%                         % independent samples, so five harmonics is the most such a
%                         % model can even be identified from
% %ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
% s.segPulsReturn={'markers'};  % markers (the scalars) is the default and runs no fit.
%                         % Add 'model' (hAmp/hPhase/R2) and/or 'reconstruction'
%                         % (timeVectors.fData) to fit every segment's beat
% %ADJUSTED IF NECESSARY - Per-pixel maps (GATE)
% s.ppxPulsReturn={'markers'};  % NON-EMPTY = the [Y x X] marker maps ON
%                         % (RESULTS.pulsatility.ppx), [] = off.  'markers' is the only
%                         % level this branch offers per pixel

%------------- BEGIN CODE --------------
function runIntensityPulsatility(s,fNames)

% THE CARDIAC INTENSITY PRODUCT AND NOTHING ELSE.  '_a_I' is a frame series with no
% averaged beat in it and '_b_I' is an injection; folding either of them is
% runIntensityInternalCycle's job and this step does not do it a second time.
if ~all( cellfun(@(x) isempty(x) || contains(x,'_c_I_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_c_I_r.mat".  This step ' ...
        'reads the averaged cardiac cycle of a fluorescence recording - list them ' ...
        'with getFileNamesList(rootFolder,''*_c_I_r.mat'').']);
end

s = withDefaults(s);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  'Pulsatility' is the
% REGISTRY LABEL, not the function name, and it is deliberately the same word the
% speckle step uses - only one of the two survives the modality filter.
rep = reportOpen(s,'Pulsatility',fNames);

for fidx = 1:1:numel(fNames)
    if reportCancelled(rep), break; end          % cooperative cancel between files
    if ~isempty(fNames{fidx})

    s.fName = fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results source settings

    load(getProductPath(s.fName,'s'),'settings');
    load(s.fName,'results');

    if ~isfield(results,'sMetrics') || ~isfield(results,'sData')
        error('runIntensityPulsatility:notSegmented', ...
            ['%s has no vessels on it yet, so there are no traces to measure a beat ' ...
             'in.  Run Segmentation on it first.'], shortName(s.fName));
    end

    % THE SOURCE IS READ ONLY FOR THE PER-PIXEL MAPS, AND ONLY ITS FIRST MEMBER.  A
    % '_c_I_d' carries the cycle, its per-pixel standard error and s.nControls
    % scrambled-timing copies of the cycle - about five times the cycle on a default
    % run.  The controls are runMotionEnhancement's floor and the SEM describes the
    % cycle rather than being a second one, so neither is read here; they are dropped
    % as soon as the file is open so the peak is the cycle alone.
    if ~isempty(s.ppxPulsReturn)
        load(getProductPath(s.fName,'d'),'source');
        if isfield(source,'control'), source.control = []; end
        if isfield(source,'sem'),     source.sem     = []; end
    end

    % runIntensityPulsatility owns RESULTS.pulsatility on this branch entirely; rebuild
    % it from scratch so no stale sub-branch survives a re-run with different levels.
    if isfield(results,'pulsatility')
        results = rmfield(results,'pulsatility');
    end

    % SETUP once for the product's own phase axis.  All the per-segment signals and the
    % per-pixel maps share it, and the period it derives is the one place the bin-centre
    % convention is accounted for.
    layout = getPulsatilityMetrics(results.time,s);
    nHarm  = layout.nHarm;  want = layout.want;  nT = numel(results.time);

    % shared root axes (single; the tree's OWN copy - results.time stays double)
    results.pulsatility.time      = single(results.time(:));
    results.pulsatility.harmonics = single((1:nHarm)');

    % =============================================================
    % Per-segment analysis.  TWO SIGNALS, NOT THREE: dvsDiameter is deliberately absent
    % and the header says why in numbers.  Both of these are intensity traces, so both
    % carry the pv prefix and the sub-tree key is what tells them apart.
    % =============================================================
    sigNames = {'sData','dvsData'};
    for kSig = 1:numel(sigNames)
        sigName = sigNames{kSig};
        if ~isfield(results,sigName), continue; end
        sigMat = results.(sigName);
        if isempty(sigMat), continue; end
        nSeg = size(sigMat,2);

        % NaN-preallocate every accumulator; an invalid segment is written by no branch
        % and therefore stays NaN across every field (invalid -> NaN).
        [mMin,mMax,mTimeMin,mTimeMax,mMean,mStd,mPI,mRI,mSym,mR2] = ...
            deal(nan(nSeg,1,'single'));
        mHAmp   = nan(nSeg,nHarm,'single');
        mHPhase = nan(nSeg,nHarm,'single');
        mFData  = nan(nT,nSeg,'single');

        % SERIAL, and on this branch that is not a compromise.  The default level runs
        % no fit at all, so the loop is a handful of reductions per segment and a pool
        % would cost more to start than the whole step takes; when the model IS asked
        % for, the Trust-Region NLS is BLAS-threading sensitive and a worker converges
        % about one single ULP away from the client, which is why the speckle step
        % keeps its per-segment loop serial too.
        for i = 1:nSeg
            m = getPulsatilityMetrics(sigMat(:,i),layout,s);
            if m.valid
                mMin(i)=m.min; mMax(i)=m.max; mTimeMin(i)=m.timeMin; mTimeMax(i)=m.timeMax;
                mMean(i)=m.mean; mStd(i)=m.std; mPI(i)=m.PI; mRI(i)=m.RI; mSym(i)=m.symRatio;
                if want.model
                    mHAmp(i,:)=m.hAmp; mHPhase(i,:)=m.hPhase; mR2(i)=m.R2;
                end
                if want.reconstruction
                    mFData(:,i)=m.fData;
                end
            end
        end

        % assemble the sub-tree; each level is gated by `want`.
        T = struct();
        if want.markers
            T.scalars.pvMin      = mMin;
            T.scalars.pvMax      = mMax;
            T.scalars.pvTimeMin  = mTimeMin;
            T.scalars.pvTimeMax  = mTimeMax;
            T.scalars.pvMean     = mMean;
            T.scalars.pvStd      = mStd;
            T.scalars.pvPI       = mPI;
            T.scalars.pvRI       = mRI;
            T.scalars.pvSymRatio = mSym;
        end
        if want.model
            T.scalars.hAmp   = mHAmp;
            T.scalars.hPhase = mHPhase;
            T.scalars.pvR2   = mR2;
        end
        if want.reconstruction
            T.timeVectors.fData = mFData;
        end
        results.pulsatility.(sigName) = T;

        % THE KEY SET, DUPLICATED INTO THE METRICS TABLES UNDER THE SAME NAMES.  The
        % markers are always computed by the core, so this runs whatever level was
        % asked for.  Row order matches the 1:nSeg loop (the segment invariant).
        switch sigName
            case 'sData'
                results.sMetrics.pvPI=mPI;             results.sMetrics.pvRI=mRI;
                results.sMetrics.pvMean=mMean;         results.sMetrics.pvTimeMin=mTimeMin;
                results.sMetrics.pvTimeMax=mTimeMax;   results.sMetrics.pvSymRatio=mSym;
            case 'dvsData'
                results.dvsMetrics.pvPI=mPI;           results.dvsMetrics.pvRI=mRI;
                results.dvsMetrics.pvMean=mMean;       results.dvsMetrics.pvTimeMin=mTimeMin;
                results.dvsMetrics.pvTimeMax=mTimeMax; results.dvsMetrics.pvSymRatio=mSym;
        end
    end

    % =============================================================
    % Per-pixel marker maps (gated by a non-empty s.ppxPulsReturn).  VECTORISED over the
    % whole cube and therefore free; there is no per-pixel fit on this branch and no
    % parfor, because the model this step already declines to run per segment would be
    % the same eleven parameters against ten samples, a million times over.
    % =============================================================
    if ~isempty(s.ppxPulsReturn)
        if ~isfield(results,'sMap')
            error('runIntensityPulsatility:noSegmentMap', ...
                ['%s carries no segment map, so there is nothing to say which pixels ' ...
                 'are vessels.  Run Segmentation on it first.'], shortName(s.fName));
        end
        W   = source.data;
        Tc  = layout.T;                       % the period, bin-centre convention included
        inv = all(isnan(W),3);                % [Y x X] all-NaN cycles (invalid -> NaN)

        [maxMap,maxIdx] = max(W,[],3);
        [minMap,minIdx] = min(W,[],3);
        tMaxMap = source.time(maxIdx);
        tMinMap = source.time(minIdx);
        wrap    = tMinMap>tMaxMap;
        tMinMap(wrap) = tMinMap(wrap)-Tc;     % wrap-shift so tMin<=tMax
        meanMap = mean(W,3,'omitnan');
        stdMap  = std(W,0,3,'omitnan');
        piMap   = (maxMap-minMap)./meanMap;
        riMap   = (maxMap-minMap)./maxMap;
        ascMap  = tMaxMap-tMinMap;
        symMap  = (Tc-ascMap)./ascMap;

        % an invalid pixel gets time(1)/Inf out of the vectorised max/min - override
        % every marker to NaN, matching the per-segment invalid rule.
        minMap(inv)=NaN;  maxMap(inv)=NaN;  tMinMap(inv)=NaN;  tMaxMap(inv)=NaN;
        meanMap(inv)=NaN; stdMap(inv)=NaN;  piMap(inv)=NaN;    riMap(inv)=NaN;
        symMap(inv)=NaN;

        ppx = struct();
        ppx.scalars.pvMin      = single(minMap);
        ppx.scalars.pvMax      = single(maxMap);
        ppx.scalars.pvTimeMin  = single(tMinMap);
        ppx.scalars.pvTimeMax  = single(tMaxMap);
        ppx.scalars.pvMean     = single(meanMap);
        ppx.scalars.pvStd      = single(stdMap);
        ppx.scalars.pvPI       = single(piMap);
        ppx.scalars.pvRI       = single(riMap);
        ppx.scalars.pvSymRatio = single(symMap);
        results.pulsatility.ppx = ppx;
    end

    writePulsatilityPage(rep,results);

    settings.runIntensityPulsatility = reportSettings(s);
    reportWriting(rep);
    %NON-DESTRUCTIVE: SOURCE (_d) is never re-saved - only RESULTS and SETTINGS.
    save(s.fName,'results','-v7.3','-nocompression');
    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
    reportSaved(rep);
    clear source
    end
end
reportClose(rep);

end

% =====================================================================
function s = withDefaults(s)
%withDefaults  Everything this step has, defaulted where a hand-written s omits one.
%   THE TWO SELECTORS DEFAULT DIFFERENTLY AND DELIBERATELY.  segPulsReturn defaults
%   when absent OR empty, because an empty level list would compute nothing at all.
%   ppxPulsReturn defaults ONLY when the field is ABSENT - an explicit [] must stay
%   empty, since [] is how the per-pixel maps are turned off; that is exactly the rule
%   runPulsatility and runVasomotion use for their own per-pixel gates.
if ~isfield(s,'nHarm') || isempty(s.nHarm)
    s.nHarm = 5;
end
if ~isfield(s,'segPulsReturn') || isempty(s.segPulsReturn)
    s.segPulsReturn = {'markers'};
end
if ~isfield(s,'ppxPulsReturn')
    s.ppxPulsReturn = {'markers'};
end
% THE PER-PIXEL PATH OFFERS ONE LEVEL, and a level it does not offer is refused by name
% rather than ignored - runVasomotion's lean per-pixel set is the precedent.  A
% per-pixel harmonic fit here would be an eleven-parameter model against about ten
% independent samples, run once per pixel: the wrong answer, slowly.
if ~isempty(s.ppxPulsReturn)
    lv = s.ppxPulsReturn;
    if ischar(lv)||isstring(lv), lv = cellstr(lv); end
    if ~iscellstr(lv) || ~all(strcmp(lv,'markers'))
        error('runIntensityPulsatility:ppxLevel', ...
            ['The per-pixel maps on this branch are the model-free markers and ' ...
             'nothing else.  Use {''markers''} to make them or [] to skip them.']);
    end
    s.ppxPulsReturn = lv;
end
end

% =====================================================================
function writePulsatilityPage(rep,results)
%writePulsatilityPage  The pulsatility index on the field, the beats it came from, how
%   it varies with calibre, and its distribution.
%
%   IT TAKES NO SETTINGS, and that is worth one line: nothing this step is configured
%   with changes what the page shows.  There is no gate to mark, no calibration to
%   convert with, and the level selectors decide what is STORED rather than what is
%   drawn - the markers behind every panel are computed whatever level was asked for.
%
%   THE HEADLINE GOES IN THE SUBTITLE AND NOT IN A TITLE.  reportSave writes the
%   recording's name across every page with sgtitle, which attaches to whichever grid
%   the figure holds - so a title set here would be silently overwritten and the page
%   would export with the numbers simply missing.  And a text object does not wrap
%   itself on a fixed canvas, so the lines are broken here or they lose both ends.
fh = reportFigure(rep,'pulsatility');
tl = tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

M = results.sMetrics;
% A LEVEL SET THAT PRODUCED NO COLUMN STILL GETS A PAGE.  results.sData can be present
% and empty (a product whose segmentation found nothing), in which case the loop above
% wrote no column at all and the page must say that rather than throw on a missing name.
if ismember('pvPI',M.Properties.VariableNames)
    pi_ = double(M.pvPI);
else
    pi_ = nan(height(M),1);
end
ok  = isfinite(pi_);
dPx = double(M.diameter);

% ---- 1. the index painted on the vessels it was measured on --------------------
ax  = nexttile(tl,1);
val = pi_;  val(~ok) = NaN;
img = nan(size(results.sMap));
lab = double(results.sMap);
inS = lab>0 & lab<=numel(val);
img(inS) = val(lab(inS));
h = imagesc(ax,100*img); set(h,'AlphaData',~isnan(img));
axis(ax,'image'); colormap(ax,parula);
if any(ok), clim(ax,[0 max(prctile(100*pi_(ok),99),eps)]); end
cb = colorbar(ax); cb.Label.String = 'Pulsatility index, %';
title(ax,sprintf('%s carry a pulsatility index',plural(nnz(ok),'vessel','vessels')))

% ---- 2. the beats the index was taken from ------------------------------------
% Every segment as a fractional change from its own mean, so a bright vessel and a dim
% one are on one axis, with the median beat over the top of them.
ax = nexttile(tl,2);
Y  = double(results.sData);
mu = mean(Y,1,'omitnan');
Y  = 100*(Y - mu)./mu;
tt = double(results.time);
okS = false(1,size(Y,2));
nS  = min(numel(ok),size(Y,2));
okS(1:nS) = ok(1:nS);
keep = find(okS & all(isfinite(Y),1));
if isempty(keep)
    plot(ax,tt,nan(size(tt)),'-');
else
    plot(ax,tt,Y(:,keep),'-','Color',[0.6 0.7 0.85],'LineWidth',0.5); hold(ax,'on')
    plot(ax,tt,median(Y(:,keep),2),'-','Color',[0.15 0.2 0.5],'LineWidth',2); hold(ax,'off')
end
grid(ax,'on'); xlabel(ax,'Time through one beat, s'); ylabel(ax,'Change from the mean, %')
title(ax,'The averaged beat of every vessel')

% ---- 3. the index against calibre ---------------------------------------------
% IN PIXELS, because this step never learns how big a pixel is: it measures no length,
% so asking for a calibration it would only put on an axis label would be a required
% setting bought for nothing.
ax = nexttile(tl,3);
[xC,yC] = orNaN(dPx(ok),100*pi_(ok));
plot(ax,xC,yC,'o','MarkerSize',4,'MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','none');
grid(ax,'on'); xlabel(ax,'Vessel width, pixels'); ylabel(ax,'Pulsatility index, %')
title(ax,'How the index varies with vessel width')

% ---- 4. and its distribution ---------------------------------------------------
ax = nexttile(tl,4);
histogram(ax,100*pi_(ok),'FaceColor',[0.2 0.4 0.8],'EdgeColor','none');
grid(ax,'on'); xlabel(ax,'Pulsatility index, %'); ylabel(ax,'Vessels')
title(ax,'How many vessels pulse how much')

notes = wrapLines([{headline(results,pi_,ok)}, pageNotes()],230);
subtitle(tl,notes,'FontSize',8,'Interpreter','none');

reportSave(rep,fh,'pulsatility');
end

% =====================================================================
function t = headline(results,pi_,ok)
%headline  One line a reader can quote, with the beat it was taken from named.
if ~any(ok)
    t = 'No vessel on this recording carries a pulsatility index.';
else
    t = sprintf(['%s pulse by %.2f to %.2f %% of their own brightness over one beat, ' ...
        '%.2f %% in the median one.'], plural(nnz(ok),'vessel','vessels'), ...
        100*min(pi_(ok)), 100*max(pi_(ok)), 100*median(pi_(ok)));
end
% THE RATE IS READ OFF THE PRODUCT, NEVER DETECTED HERE.  The cardiac gate has three
% implementations in this library already and a fourth in a consumer is how two steps
% come to disagree about which beats they averaged.
if isfield(results,'rateMedian') && isfinite(results.rateMedian)
    t = sprintf('%s Averaged from %d beats at %.2f per second.', ...
        t, round(getfielddef(results,'nCycles',NaN)), results.rateMedian);
end
end

% =====================================================================
function c = pageNotes()
%pageNotes  The sentences that say what these numbers are and what they are NOT.
c = {['The index is how much the labelled plasma in a vessel changes over one beat, ' ...
      'as a fraction of its own average. It is a volume, not a speed and not a width ' ...
      '- a vessel can pulse strongly here and change its diameter hardly at all.'], ...
     ['How much each vessel''s WIDTH changed over the same beat is measured by the ' ...
      'wall motion step and reported in its own columns, against a copy of the ' ...
      'recording with the heartbeat''s timing scrambled away. There is deliberately ' ...
      'no width column here: the per-frame diameter this library measures cannot ' ...
      'change by less than a quarter of a pixel, and a cardiac width change on this ' ...
      'kind of recording is one and a half to eighteen times smaller than that.'], ...
     ['Every beat above is an average of hundreds of real ones, so what is drawn is ' ...
      'already far quieter than any single heartbeat was.']};
end

% =====================================================================
function [x,y] = orNaN(x,y)
%orNaN  One invisible point instead of nothing, so a series always has a handle.
if isempty(x), x = NaN; y = NaN; end
end

% =====================================================================
function c = wrapLines(c,n)
%wrapLines  Break each line on a word boundary so nothing runs off the page.
%   The canvas is a fixed 1400 px and a text object does not wrap itself: unwrapped, a
%   long line loses BOTH of its ends.
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
function t = plural(n,one,many)
%plural  "1 vessel" rather than "1 vessels", on a page a biologist reads.
if n==1, t = sprintf('%d %s',n,one); else, t = sprintf('%d %s',n,many); end
end

% =====================================================================
function v = getfielddef(st,f,d)
%getfielddef  One field of a struct, or a stated default.
if isfield(st,f) && ~isempty(st.(f)), v = st.(f); else, v = d; end
end

% =====================================================================
function n = shortName(p)
%shortName  A file's name without its folder, for a message a person reads.
[~,b,e] = fileparts(char(p));
n = [b e];
end
%------------- END OF CODE --------------
