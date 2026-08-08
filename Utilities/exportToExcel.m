%exportToExcel  Dump processed metrics and traces to an .xlsx workbook
%
%   exportToExcel(fNames) opens each processed dataset in fNames - the RESULTS
%   member of the triplet (*_BFI_r.mat or later) - converts numeric
%   category codes to readable text, builds several summary tables, and
%   writes them as separate sheets in a single Excel file whose name mirrors
%   the input file (e.g.  Foo_BFI_r.mat → Foo_BFI.xlsx):
%
%       Sheet             Content
%       ─────────────────────────────────────────────────────────────────
%       sMetrics          raw per-ROI metrics table
%       sData             time-series per ROI           (Time + ROI-cols)
%       sMetricsROI       area/length-summed, area-weighted ROI averages
%       sDataROI          area-weighted mean time-series per text label
%       sMetricsType      (opts.groupBy='type') the same aggregation grouped by
%                         VESSEL type (artery/vein/…, results.sMetrics.type)
%       sDataType         (opts.groupBy='type') mean time-series per vessel type
%       dvsMetrics        (if present) raw dynamic-vessel metrics
%       dvsData           (if present) flow traces per vessel ROI
%       dvsDiameter       (if present) diameter traces per vessel ROI
%       pulsatility       (if present) per-segment pulsatility tree scalars,
%                         aligned to sMetrics rows by idx
%       dvsPulsatility    (if present) per dynamic vessel: every pulsatility signal
%                         the file carries, side by side, aligned to dvsMetrics rows
%       wallMotion        (if present) per segment: the wall-motion measurement
%                         BEFORE its gates, the floor it was judged against, and
%                         which gate refused a row
%       topology          (if present) one row per ANALYSED AREA - the whole
%                         segmented field, then one per drawn region - with the
%                         vascular densities and the units they are in
%       nvcTrials         (if present) LONG: one row per segment per stimulus
%                         repetition - idx, the repetition, where it was cut, whether
%                         it was kept, and every response marker of the flow (ns*)
%                         and guided-contrast (ng*) traces, each with the two
%                         confidence numbers beside it
%       dvsNvcTrials      (if present) the same, one row per tracked vessel per
%                         repetition: flow (ns*) and diameter (nd*)
%
%   NOTE  runPulsatility is now non-destructive: the sData / dvsData / dvsDiameter
%   trace sheets carry the RAW averaged cycle (its harmonic-fit reconstruction
%   lives in results.pulsatility, not in these traces); sMetrics / dvsMetrics gain
%   the pulsatility columns.  The per-pixel maps results.pulsatility.ppx
%   are image data and are NOT exported here.
%
%   THE PULSATILITY PREFIX NAMES THE QUANTITY AND THREE OF THEM EXIST.  ps* is a
%   pulsatile FLOW, pd* a pulsatile DIAMETER, pv* a pulsatile PLASMA VOLUME - the
%   fluorescence branch's, where an intravascular tracer makes the intensity
%   proportional to the labelled plasma - and NO TWO OF THEM MAY BE POOLED, in a column
%   or in a figure.  The producers write the prefix into the marker names, so this file
%   names none of them: it reads the prefix off the tree and uses it for the fitted
%   model's hAmp / hPhase, which are the only two leaves the producers leave unprefixed
%   and therefore the only two that could otherwise be stacked across branches in one
%   merged column.  A wm* column, from the wall-motion step, is a fourth quantity again:
%   a wall DISPLACEMENT, peak to peak, per segment.
%
%   A MIXED WORKING SET IS AN ORDINARY CASE NOW, and two things make it work.  Files
%   legitimately carry disjoint column sets - a speckle recording has a BFI and a
%   pv*-free metrics table, an intensity one has neither BFI nor ps* - so every column
%   test here asks the table what it holds rather than naming a column, and merged mode
%   takes the UNION of the sheets it stacks (stackTables), a file missing a column
%   contributing empty cells there rather than dropping everybody else's.
%
%   THE NVC SHEETS ARE THE ONLY LONG ONES, and they have to be.  Every other metric
%   sheet is one row per segment, because every other analysis gives a segment ONE
%   answer; runNVC gives it one per stimulus repetition, and that second axis is the
%   whole point of the step - a response that fades over a session is a slope along
%   it, and an average of twenty repetitions hides it completely.  A wide sheet would
%   need a column per (metric, repetition) and could not be filtered or grouped, so
%   the repetition becomes a COLUMN and the rows multiply: 1940 segments x 20
%   repetitions is 38 800 rows, well inside Excel's 1 048 576.
%
%   A REPRESENTATIVE-REPETITION PRODUCT WRITES THE SAME SHEET, ONE ROW PER SEGMENT.
%   runNVC's optional collapse averages the trusted repetitions over each other and
%   replaces the recording with that average, so the product has no repetitions left
%   to be long in: the epoch columns are ABSENT rather than filled with blanks, and
%   the markers are the markers of the average.  It is also the only place those
%   markers appear at all - the metrics table carries five of them and this sheet
%   carries every one.
%
%   The trusted-repetition mean of five markers lands in sMetrics / dvsMetrics (and
%   so in sMetricsROI / sMetricsType) on its own, because runNVC writes those into
%   the tables as ordinary columns; their names come from the analysis settings, so
%   this file names none of them.  The per-pixel maps results.nvc.ppx are not
%   exported: they are image data, exactly like results.pulsatility.ppx.
%
%   MYOGRAPH RECORDINGS TAKE A SECOND SHEET SET, chosen FROM THE DATA and never
%   from the file name: a recording whose RESULTS carry analysed WINDOWS is a
%   myograph recording (pressure or wire), it has no per-segment table for the
%   sheets above to be about, and it writes these instead.  A pressure myograph
%   keeps them in a flat results.intervals; a wire myograph splits them by CHANNEL
%   (results.channel(i).intervals(k)) because one LabChart file is several
%   chambers, and myographIntervals flattens the two into one list of rows:
%
%       Sheet             Content
%       ─────────────────────────────────────────────────────────────────
%       settings          one row per parameter: which step wrote it, its name and
%                         its value, plus the recording's analysed window and code
%                         version - read from the triplet's SETTINGS member
%       comments          WIRE MYOGRAPH ONLY: the operator's LabChart log - time,
%                         text, channel and record - which is what the analysed
%                         windows were placed against.  A pressure myograph has no
%                         log and skips the sheet
%       intervals         one row per analysed window (per CHANNEL, on a wire
%                         myograph): its name, the channel it is about, what the
%                         operator assigned it to, start, end, duration, frame
%                         count, and the statistics of every diameter measured
%                         inside it
%       propagation       interval × diameter measure: speed and direction, the
%                         confidence interval, R², p, and the evidence behind them
%       vasomotion        interval × signal × unit: the band-suffixed vasomotion
%                         scalars (scalars.VB/CB.*), the same markers a speckle
%                         recording reports per segment
%       ampPct            interval × signal × unit × band × percentile: the band
%                         envelope amplitude at each percentile level
%       spectra           interval × signal: amplitude spectrum mean and standard
%                         deviation vs frequency, averaged over units, decimated
%       ampPctSpectra     interval × signal: the vasomotion-band spectrum SPLIT BY
%                         envelope amplitude - one column per amplitude bin, so a
%                         change in the size of the oscillation can be told from a
%                         change in its shape.  Decimated in frequency
%       diameterTraces    interval × time: the line-averaged diameter of each
%                         measure, decimated
%
%   THE SHEET LIMITS ARE RESPECTED BY DECIMATING, and the factor is written into
%   the sheet, so a reader can always see how coarse a curve is: both spectra
%   sheets are decimated in frequency and the diameter traces in time (opts.fDecim /
%   opts.tDecim override the automatic choice), and no sheet exceeds opts.maxRows.
%
%   A SIGNAL is what was analysed: the diameter MEASURE ('outer' / 'mid' /
%   'inner') for a pressure myograph, the CHANNEL for a wire one.  A UNIT is one
%   result inside it - a single line-averaged trace, or one image row when the
%   vasomotion was run per line.
%
%   INPUT
%     fNames   cell array of *.mat file paths.  Each must contain SOURCE,
%              RESULTS, SETTINGS structures created by the LSCI pipeline.
%     opts     (optional) struct selecting WHAT, HOW AVERAGED and in WHICH format
%              to write.  Omitting it - or passing an empty / no-field struct -
%              reproduces the full legacy behaviour above exactly (the five
%              launcher callers do; guiExport is the interactive front-end that
%              fills the struct in).  EVERY field is additive and default-off.
%              Fields:
%                .sheets  cellstr / string of sheet names to write, a subset of
%                         the UNION of both sheet sets - {sMetrics, sData,
%                          sMetricsROI, sDataROI, sMetricsType, sDataType,
%                          dvsMetrics, dvsData, dvsDiameter, pulsatility,
%                          dvsPulsatility, wallMotion, topology, nvcTrials,
%                          dvsNvcTrials} for a segmented recording and
%                         {settings, comments, intervals, propagation, vasomotion,
%                          ampPct, spectra, ampPctSpectra, diameterTraces} for a
%                         myograph one.
%                         Absent / empty = the default set for .groupBy (below).
%                         A selected sheet whose data is absent from a given file
%                         is simply skipped, exactly as in the full path - which is
%                         also how a mixed list of both kinds writes each file the
%                         sheets that are about it and no others.
%                .format  output extension, '.xlsx' (default) or '.xls'.
%                .groupBy the ROW GRANULARITY of the metric sheets, and with it
%                         the DEFAULT sheet set:
%                           ''       (default) per-segment rows - all 11 legacy
%                                    sheets, including the label-averaged ROI pair,
%                                    plus the 9 myograph sheets, so a bare call
%                                    writes whichever of them the file is about;
%                           'label'  average over results.sMetrics.label - the
%                                    sMetricsROI / sDataROI pair alone (the SAME
%                                    aggregation the legacy path writes, not a
%                                    second implementation);
%                           'type'   average over results.sMetrics.type, the
%                                    VESSEL type written by setVesselTypes - the
%                                    sMetricsType / sDataType pair alone.
%                         Both averaging axes are per-SEGMENT axes, so a myograph
%                         recording contributes nothing to them; its own averaging
%                         axis is the interval, which it always has.
%                         .sheets, when given, still selects freely across all
%                         twenty-two names.
%                .weightByArea  logical, TRUE by default: aggregated rows are
%                         weighted by results.sMetrics.area, the weights the ROI
%                         sheets have always used.  FALSE gives a plain unweighted
%                         mean.  It applies to the .groupBy AXIS ONLY - every other
%                         aggregation in the same call keeps its historical area
%                         weighting, so a legacy sheet can never change silently.
%                .columns cellstr subset of the METRIC-SHEET variable names to
%                         export.  Absent = all.  Taken LITERALLY (include 'idx' /
%                         'label' / 'type' yourself if you want them) and applied to
%                         the metric sheets only - sMetrics / sMetricsROI /
%                         sMetricsType for a segmented recording, intervals /
%                         propagation / vasomotion for a myograph one.  The trace,
%                         dvs, spectra and percentile sheets are untouched.  A
%                         column absent from a given file is skipped for that file,
%                         exactly as an absent sheet already is; a sheet left with
%                         no column at all is skipped for that file.
%                .fDecim  MYOGRAPH ONLY.  Keep every Nth frequency point of the
%                         spectra sheet.  Absent = automatic (about 60 points per
%                         curve), which is what keeps a sheet readable.
%                .tDecim  MYOGRAPH ONLY.  Keep every Nth sample of the diameter
%                         traces.  Absent = automatic (about 5000 points).
%                .maxRows MYOGRAPH ONLY.  Row ceiling per sheet (default 1e6, under
%                         Excel's 1,048,576); a sheet that reaches it is written
%                         truncated rather than failing.
%                .merge   logical.  FALSE (default) = one workbook per input file,
%                         as always.  TRUE = ONE workbook for the whole call, its
%                         sheets accumulated across files and every row prefixed
%                         with 'file' plus - when .labels supplies them - 'animal',
%                         'type' (the RECORDING type) and 'group' (the experimental
%                         group), so the comparison labels survive into the sheet
%                         and the workbook is directly usable for statistics.
%                         A per-segment 'type' column (the VESSEL type) is renamed
%                         'vesselType' in merged mode so the two axes stay distinct.
%                         Sheets are unioned column-wise: a file missing a column
%                         contributes empty cells there rather than dropping out.
%                .labels  optional 1xN struct array PARALLEL to fNames carrying the
%                         per-file labels for merged mode: .animal, .type and
%                         .group (or .expGroup, the session's name for it).  The
%                         GUI fills it from the workbench session; a launcher
%                         calling exportToExcel bare supplies nothing.
%                .outFile destination of the merged workbook.  Absent = a
%                         'mergedExport' workbook beside the first input file.
%              The first sheet actually written to a file uses 'replacefile', in
%              BOTH modes, so a re-export with a narrower selection never leaves
%              stale sheets behind.
%
%   OUTPUT
%     None – one Excel workbook per input file (or one merged workbook) is
%     written to disk.
%
%   EXAMPLE
%     files = dir(fullfile(dataRoot,'*BFI_r.mat'));
%     fNames = fullfile({files.folder}',{files.name}');
%     exportToExcel(fNames);                              % full legacy workbook
%     exportToExcel(fNames, struct('sheets',{{'sMetrics','pulsatility'}}));
%                                                         % just those two sheets
%     exportToExcel(fNames, struct('groupBy','type','weightByArea',false));
%                                                         % unweighted per-vessel-type means
%     exportToExcel(fNames, struct('merge',true,'labels',L,'outFile',out));
%                                                         % one workbook for statistics
%     exportToExcel({'Mouse1_MYO_r.mat'});                 % the myograph workbook
%
%   DEPENDS ON
%     MATLAB R2019b+ (for writetable with 'Sheet' option) and data schema
%     used throughout the Dynamic Light Scattering Imaging toolbox.
%
% See also: guiExport, setVesselTypes, guiWorkbench, wbSession, myographProduct,
%           runNVC, runIntensityPulsatility, runMotionEnhancement,
%           runTopologyAnalysis, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

function exportToExcel(fNames, opts)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_r.mat".  Export reads ' ...
        'the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_BFI_r.mat'').']);
end

if nargin<2, opts = struct(); end
O = parseExportOpts(opts, fNames);                  % the legacy 9 sheets; '.xlsx' default
wantSheet = @(nm) any(strcmp(nm, O.sheets));
sink = newSink(O);                                  % per-file writer OR merged accumulator

for fidx=1:1:numel(fNames)
     if ~isempty(fNames{fidx})
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    fName=fNames{fidx};
    load(fName,'results');

    fName=strrep(fName,'_r.mat',O.ext);
    sink=openFile(sink,fName,fNames{fidx},fidx);        % target file / merged row prefix

    % WHICH SHEET SET, DECIDED BY THE DATA.  A myograph recording carries analysed
    % WINDOWS and no per-segment table; a segmented recording carries the table and
    % no windows.  Dispatching on the window fields rather than on the file name is
    % what lets one call take a mixed list and give each file the sheets that are
    % about it - and it is the only place either kind is recognised.  BOTH are
    % tested: a wire myograph stores its windows under results.channel and leaves
    % results.intervals empty, so a test on the flat field alone would hand a wire
    % recording to the segment exporter.
    if isfield(results,'intervals') || isfield(results,'channel')
        sink=emitMyographSheets(sink,results,O,wantSheet,fNames{fidx});
    else
        sink=emitSegmentSheets(sink,results,O,wantSheet);
    end
     end
end

finishSink(sink);                                   % merged mode: write the one workbook
end

function sink = emitSegmentSheets(sink, results, O, wantSheet)
%emitSegmentSheets  THE per-segment sheet set: the nine legacy sheets plus the two
%   opt-in vessel-type ones and the two the fluorescence branch added.  This is the
%   body a bare exportToExcel(fNames) has always run - the myograph branch is beside
%   it, never inside it.
%
%   A PRODUCT THAT HAS NOT BEEN SEGMENTED HAS NONE OF THESE SHEETS TO WRITE, and on the
%   fluorescence branch that is an ordinary file rather than a broken one:
%   Launcher_intensity's export cell globs every product of every branch, and an
%   angiogram straight from the entry step carries a picture and a clock and no table.
%   It is skipped exactly as an absent sheet is, which leaves no workbook for that file
%   rather than an empty one - and, more to the point, leaves the OTHER files of the
%   call exported instead of stopping the batch on a missing column.

    if ~isfield(results,'sMetrics') || ~istable(results.sMetrics) ...
            || height(results.sMetrics)==0
        return
    end

    catNames = ["background"  "parenchyma"  "unsegmented" ...
        "outerWall"   "innerWall"   "lumen"];          % 1×6

    results.sMetrics.category(isnan( results.sMetrics.category))=0;
    results.sMetrics.category = categorical( ...
        results.sMetrics.category, 0:5, catNames);          % numeric → cat


    T=results.sMetrics;
    T(end,:) = []; % crappy bug workaround
    T(end+1,:) = results.sMetrics(end,:);% crappy bug workaround

    sink=emitSheet(sink,pickColumns(T,O.columns),'sMetrics',wantSheet);


    idxAll    = results.sMetrics.idx(:)';
    goodMask  = measuredRows(results.sMetrics);
    if isfield(results,'sData') && isfield(results,'time')
        roiIdx    = idxAll(goodMask);
        roiData   = results.sData(:,goodMask);
        roiNames  = matlab.lang.makeValidName(compose("ROI %04d",roiIdx) );
        T = [ table(results.time(:),'VariableNames', {'Time'}), array2table(roiData,'VariableNames',roiNames) ];

        sink=emitSheet(sink,T,'sData',wantSheet);
    end

    % ---- averages over the LABEL axis (the historical ROI sheets) -----------
    % One aggregation routine serves both axes: this call is the legacy
    % area-weighted 'label' grouping, unchanged.
    if wantSheet('sMetricsROI') || wantSheet('sDataROI')
        [Tm,Td] = aggregateSegments(results,'label','type',axisWeight(O,'label'));
        if ~isempty(Tm)
            sink=emitSheet(sink,pickColumns(Tm,O.columns),'sMetricsROI',wantSheet);
            sink=emitSheet(sink,Td,'sDataROI',wantSheet);
        end
    end

    % ---- averages over the VESSEL-TYPE axis (opts.groupBy='type') -----------
    % The axis guiExplore already compares on (artery / vein / …), written by
    % setVesselTypes; same machinery, grouping variable swapped.
    if wantSheet('sMetricsType') || wantSheet('sDataType')
        [Tm,Td] = aggregateSegments(results,'type','label',axisWeight(O,'type'));
        if ~isempty(Tm)
            sink=emitSheet(sink,pickColumns(Tm,O.columns),'sMetricsType',wantSheet);
            sink=emitSheet(sink,Td,'sDataType',wantSheet);
        end
    end


    if isfield(results,"dvsMetrics")
        T=results.dvsMetrics;
        T(end,:) = []; % crappy bug workaround
        T(end+1,:) = results.dvsMetrics(end,:);% crappy bug workaround
        sink=emitSheet(sink,T,'dvsMetrics',wantSheet);

        idxAll   = results.dvsMetrics.idx(:)';       % row vector
        goodMask = measuredRows(results.dvsMetrics);
        roiIdx   = idxAll(goodMask);
        roiNames = matlab.lang.makeValidName(compose("ROI %04d", roiIdx));

        if isfield(results,'dvsData') && isfield(results,'time')
            T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
                array2table(results.dvsData(:,goodMask), ...
                'VariableNames',roiNames) ];

            sink=emitSheet(sink,T,'dvsData',wantSheet);
        end

        if isfield(results,'dvsDiameter') && isfield(results,'time')
            T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
                array2table(results.dvsDiameter(:,goodMask), ...
                'VariableNames',roiNames) ];

            sink=emitSheet(sink,T,'dvsDiameter',wantSheet);
        end
    end

    % ---- pulsatility summary sheets (per-segment results.pulsatility scalars) ---
    % Surface the tree's per-segment markers + harmonic coefficients as focused
    % sheets, row-aligned to sMetrics / dvsMetrics.  The per-pixel maps
    % results.pulsatility.ppx are image data and are intentionally NOT exported.
    if isfield(results,'pulsatility')
        if isfield(results.pulsatility,'sData') && isfield(results.pulsatility.sData,'scalars')
            Tp = pulsScalarTable(results.pulsatility.sData.scalars);
            if ~isempty(Tp)
                sink=emitSheet(sink,[idTable(results.sMetrics), Tp],'pulsatility',wantSheet);
            end
        end
        if isfield(results,'dvsMetrics') && isfield(results.pulsatility,'dvsData') ...
                && isfield(results.pulsatility.dvsData,'scalars')
            Tp = [idTable(results.dvsMetrics), ...
                  pulsScalarTable(results.pulsatility.dvsData.scalars)];
            if isfield(results.pulsatility,'dvsDiameter') ...
                    && isfield(results.pulsatility.dvsDiameter,'scalars')
                Tp = [Tp, pulsScalarTable(results.pulsatility.dvsDiameter.scalars)];
            end
            sink=emitSheet(sink,Tp,'dvsPulsatility',wantSheet);
        end
    end

    % ---- the wall motion, one row per segment (results.wallMotion) --------------
    % The metrics sheet already carries the eight wm* columns, which are the answer
    % AFTER three gates and are NaN in all eight wherever any of them refused.  This is
    % the evidence underneath: the same measurements before the gates, the floor each
    % was judged against, and which gate refused a row.  A per-row struct is a sheet.
    sink=emitLazy(sink,wantSheet,'wallMotion',@() wallMotionTable(results));

    % ---- the vascular density, one row per analysed area (results.topology) -----
    % NOT one row per segment, and that is the whole shape of it: a density is a
    % property of an AREA - the whole analysed field first, then one row per drawn
    % region - so it joins nothing and stands on its own.
    sink=emitLazy(sink,wantSheet,'topology',@() topologyTable(results));

    % ---- the per-repetition NVC sheets (results.nvc.esMetrics) ------------------
    % LONG, one row per (unit, repetition) - see the header for why this is the one
    % shape that can carry the epoch axis.  Built lazily, because a 38 800-row table
    % is not something a narrower selection should pay for.
    if isfield(results,'nvc') && isfield(results,'sMetrics')
        sink=emitLazy(sink,wantSheet,'nvcTrials', ...
            @() nvcTrialTable(results,results.sMetrics,{'sData','gsData'}));
    end
    if isfield(results,'nvc') && isfield(results,'dvsMetrics')
        sink=emitLazy(sink,wantSheet,'dvsNvcTrials', ...
            @() nvcTrialTable(results,results.dvsMetrics,{'dvsData','dvsDiameter'}));
    end
end

function T = nvcTrialTable(results, M, sigs)
%nvcTrialTable  ONE ROW PER UNIT PER STIMULUS REPETITION, for the signals that share
%   one metrics table.  sData and gsData are both per SEGMENT and go in the sMetrics
%   sheet; dvsData and dvsDiameter are both per tracked VESSEL and go in the
%   dvsMetrics one - the same pairing dvsPulsatility already uses, and for the same
%   reason: two signals measured on one set of rows are columns of one sheet, not two
%   sheets that have to be joined by hand afterwards.
%
%   THE COLUMNS COME FROM THE TREE AND NEVER FROM A LIST HERE.  runNVC emits 17
%   numbers per signal with the default levels - the fourteen markers, which
%   s.segNvcReturn gates, and the two confidence numbers and the trust flag, which it
%   never does - and three of the fourteen are named from s.nvcAreaPcts, so a list
%   typed in this file could not even be written down without reading the settings.
%   The names are already prefixed (ns / ng / nd) by the producer, so two signals
%   never collide in one sheet even where they measure the same quantity - which is
%   exactly what keeps a diameter response out of a flow column.
%
%   EVERY REPETITION IS A ROW, UNTRUSTED ONES INCLUDED, and epochTrust says which the
%   recording kept.  Dropping them here would hide the thing they were dropped for,
%   and runNVC deliberately measures them all for the same reason.  epochTrust is the
%   RECORDING'S decision, one flag per repetition; whether a repetition was any good
%   for THIS unit is the ns*Trust column beside the markers, and the two are
%   different questions.
%
%   A COLLAPSED PRODUCT TAKES THE ONE BRANCH IN THIS FUNCTION.  It has no
%   .epochStart, so there is nothing to number, nothing to place on the recording
%   clock and no per-repetition decision: the three epoch columns are left OUT rather
%   than written as blanks, and the sheet is one row per unit.
T = table();
N = results.nvc;
if ~isfield(N,'esMetrics'), return; end

nUnit = height(M);
if nUnit==0, return; end
nEp   = numel(fieldOrEmpty(N,'epochStart'));
collapsed = nEp==0;
if collapsed, nEp = 1; end

% The identity block, one row per (unit, repetition).  THE ROWS RUN UNIT-FASTEST,
% which is not a free choice: a [nUnit x nEp] metric flattens column-major, so v(:)
% already runs that way and the two orders have to agree or every number lands
% against the wrong repetition.
T = idTable(M);
T = T(repmat((1:nUnit)', nEp, 1), :);
if ~collapsed
    T.epoch      = reshape(repmat(1:nEp, nUnit, 1), [], 1);
    T.epochStart = reshape(repmat(double(N.epochStart(:))', nUnit, 1), [], 1);
    trust = [];
    if isfield(N,'epochTrust'), trust = reshape(logical(N.epochTrust),1,[]); end
    if numel(trust)==nEp, T.epochTrust = reshape(repmat(trust, nUnit, 1), [], 1); end
end

nMet = 0;
for k = 1:numel(sigs)
    if ~isfield(N.esMetrics, sigs{k}), continue; end
    S  = N.esMetrics.(sigs{k});
    fn = fieldnames(S);
    for j = 1:numel(fn)
        v = double(S.(fn{j}));
        % A metric of a different height belongs to a different set of units and is
        % skipped rather than reshaped into this sheet - the same "skip what does not
        % fit this file" rule every other sheet here follows.
        if size(v,1)~=nUnit || size(v,2)~=nEp, continue; end
        T.(fn{j}) = v(:);
        nMet = nMet + 1;
    end
end
% No metric of these signals fits this table, so the sheet would be an index and
% nothing else - skipped, like an absent sheet.
if nMet==0, T = table(); end
end

%% ===================== THE MYOGRAPH SHEET SET ========================== %%
function sink = emitMyographSheets(sink, results, O, wantSheet, rPath)
%emitMyographSheets  THE myograph sheet set: the windows a recording was analysed
%   in, what was measured inside each of them, and what it was all measured with.
%
%   A myograph recording has no segments, so there is no per-segment table for the
%   sheets above to be about; its rows are INTERVALS, and inside an interval a
%   SIGNAL (the diameter measure, or the channel) and a UNIT (the line-averaged
%   trace, or one image row when the vasomotion ran per line).  Everything else -
%   which sheets are written, the merged row prefix, the column narrowing, the
%   workbook - is the ordinary machinery, unchanged.
%
%   THE SETTINGS SHEET IS THE ONE THING THAT NEEDS A SECOND FILE.  Parameters live
%   in the triplet's own SETTINGS member, not in its RESULTS, so this branch opens
%   it - and only it, and only when that sheet is wanted.  A recording whose
%   settings file has gone missing simply contributes no settings sheet, which is
%   the same "skip what is absent" rule everything else here follows.
%
%   THE TABLES ARE BUILT ONLY WHEN THEY ARE ASKED FOR.  The spectra and the traces
%   are the expensive ones in this file and a narrower selection must not pay for
%   them, so each goes through emitLazy rather than being built and then dropped.
%
%   THE WINDOWS COME THROUGH myographIntervals, which is where the two shapes of the
%   tree are known about: a pressure myograph's flat results.intervals and a wire
%   myograph's results.channel(i).intervals(k).  What comes back is one flat list
%   with the channel written onto every element, so a sheet is one row per window
%   per channel and nothing here has to ask which myograph it is looking at.
iv = myographIntervals(results);
if isempty(iv), return; end

sink = emitLazy(sink,wantSheet,'settings',       @() settingsTable(results,rPath));
sink = emitLazy(sink,wantSheet,'comments',       @() commentsTable(results));
sink = emitLazy(sink,wantSheet,'intervals',      @() pickColumns(intervalsTable(iv),O.columns));
sink = emitLazy(sink,wantSheet,'propagation',    @() pickColumns(propagationTable(iv),O.columns));
sink = emitLazy(sink,wantSheet,'vasomotion',     @() pickColumns(vasomotionTable(iv),O.columns));
sink = emitLazy(sink,wantSheet,'ampPct',         @() ampPctTable(iv,O));
sink = emitLazy(sink,wantSheet,'spectra',        @() spectraTable(iv,O));
sink = emitLazy(sink,wantSheet,'ampPctSpectra',  @() ampPctSpectraTable(iv,O));
sink = emitLazy(sink,wantSheet,'diameterTraces', @() diameterTraceTable(iv,O));
end

function T = settingsTable(results, rPath)
%settingsTable  WHAT THE RECORDING WAS PROCESSED WITH: one row per parameter, so
%   the sheet stacks across recordings in a merged workbook and two protocols can be
%   read side by side.  (The old myograph exporter wrote this as a two-column block
%   per file with a '=== file N ===' separator, which cannot be merged or filtered.)
%
%   THE 'step' COLUMN IS THE FUNCTION THAT WROTE THE PARAMETERS, verbatim - that is
%   what settings.<field> is named after, and this sheet's whole use is tracing a
%   number back to the call that would reproduce it, or pasting it into a launcher.
%   Plus a leading 'recording' block for the provenance the settings do not carry:
%   which absolute window was analysed, and which version of the code did it.
T = table();
C = {};
% ---- the recording's own provenance ----
C = addSetting(C,'recording','timeCrop',      fieldOrEmpty(results,'timeCrop'));
if isfield(results,'meta') && isstruct(results.meta)
    for f = {'formatVersion','codeVersion','createdTimestamp'}
        C = addSetting(C,'recording',f{1},fieldOrEmpty(results.meta,f{1}));
    end
end
% ---- one block per step that wrote settings ----
S = loadSettings(rPath);
fn = fieldnames(S);
for i = 1:numel(fn)
    st = S.(fn{i});
    if ~isstruct(st) || ~isscalar(st), continue; end
    p = fieldnames(st);
    for j = 1:numel(p)
        C = addSetting(C,fn{i},p{j},st.(p{j}));
    end
end
if isempty(C), return; end
T = cell2table(C,'VariableNames',{'step','parameter','value'});
end

function T = commentsTable(results)
%commentsTable  THE OPERATOR'S LOG, as LabChart recorded it: what was done to the
%   preparation and when.  It is the record the analysed windows were placed
%   against - a drug added, a wash started, a stimulus given - and until now it
%   never left the interval editor, so a finished export could not say what any
%   window WAS.  Long rather than wide, so it stacks across recordings in a merged
%   workbook and can be filtered like every other sheet here.
%
%   Times are ABSOLUTE seconds from the start of the recording, the same base the
%   intervals sheet uses, so the two can be read against each other directly.  A
%   pressure myograph has no operator log and simply contributes no sheet, which is
%   the ordinary "skip what is absent" rule.
T = table();
if ~isfield(results,'comments'), return; end
cm = results.comments;
if isempty(cm) || ~isstruct(cm), return; end
n = numel(cm);
T = table('Size',[n 4],'VariableTypes',{'double','string','string','double'}, ...
    'VariableNames',{'time_s','text','channel','record'});
for k = 1:n
    T.time_s(k)  = fieldNum(cm(k),'time');
    T.text(k)    = string(charOrEmpty(fieldOrEmpty(cm(k),'text')));
    T.channel(k) = string(charOrEmpty(fieldOrEmpty(cm(k),'channel')));
    T.record(k)  = fieldNum(cm(k),'record');
end
end

function s = charOrEmpty(v)
%charOrEmpty  A comment's text or channel as char, whatever the SDK handed back.
s = '';
if ischar(v), s = v; elseif isstring(v) && ~isempty(v), s = char(v(1)); end
end

function S = loadSettings(rPath)
%loadSettings  The SETTINGS member of this recording's triplet, or an empty struct.
S = struct();
if isempty(rPath), return; end
sName = getProductPath(rPath,'s');
if ~isfile(sName), return; end
try
    L = load(sName,'settings');
    if isfield(L,'settings') && isstruct(L.settings), S = L.settings; end
catch
end
end

function C = addSetting(C, step, name, v)
%addSetting  One parameter row, with the value rendered as text a person can read
%   AND paste back into a script.  Anything that is not a plain value - a nested
%   struct, a function handle a caller failed to strip - is skipped rather than
%   printed as a class name, because a settings sheet that lies is worse than one
%   with a gap in it.
if isempty(v) || isa(v,'function_handle') || isstruct(v), return; end
if isnumeric(v) || islogical(v)
    txt = mat2str(double(v));
elseif ischar(v)
    txt = v;
elseif isstring(v)
    txt = strjoin(cellstr(v(:)'),', ');
elseif iscellstr(v)                  % the shape settings use for a list of names
    txt = strjoin(v(:)',', ');
else
    return
end
C(end+1,:) = {step, name, txt};
end

function v = fieldOrEmpty(st,f)
v = [];
if isstruct(st) && isfield(st,f), v = st.(f); end
end

function T = ampPctSpectraTable(iv,O)
%ampPctSpectraTable  The vasomotion-band spectrum SPLIT BY ENVELOPE AMPLITUDE:
%   fVectors.VB.ampMeanPct is the mean spectrum of the moments when the band
%   envelope sat in each amplitude bin, so it separates "the oscillation is always
%   this shape, just bigger sometimes" from "big and small episodes are different
%   oscillations" - which one number per interval cannot.
%
%   Averaged over the units and DECIMATED IN FREQUENCY, one column per bin labelled
%   by its bin CENTRE.  NOTE those centres are not the percentile LEVELS that index
%   the scalar ampPct sheet - they key different things and are deliberately kept in
%   different sheets.  The columns are fixed by the FIRST tree that has them and
%   later ones are padded or trimmed to match, so one recording's sheet has one
%   shape; a recording binned differently unions its own columns in a merged book.
T = table();
C = {}; labs = {}; nB = 0;
for k = 1:numel(iv)
    if atCap(C,O), break, end
    V = subStruct(iv(k),'vasomotion');
    sigs = fieldnames(V);
    for j = 1:numel(sigs)
        if atCap(C,O), break, end
        v = V.(sigs{j});
        if ~isstruct(v) || ~isfield(v,'fVectors') || ~isfield(v.fVectors,'VB') || ...
                ~isfield(v.fVectors.VB,'ampMeanPct') || isempty(v.fVectors.VB.ampMeanPct)
            continue
        end
        f = double(v.f(:)); nf = numel(f);
        if nf==0, continue; end
        M = squeeze(mean(double(v.fVectors.VB.ampMeanPct),1,'omitnan'));   % [nF x nBin]
        if isvector(M), M = M(:); end
        if isempty(labs)
            nB = size(M,2); labs = binLabels(v,nB);
        end
        dec = decimation(O.fDecim, nf, 60);
        for q = 1:dec:nf
            C(end+1,:) = [{intervalName(iv(k),k), signalName(v,sigs{j}), f(q)}, ...
                num2cell(padRow(M(q,:),nB)), {dec}]; %#ok<AGROW>
        end
    end
end
if isempty(C), return; end
T = cell2table(C,'VariableNames',[{'interval','signal','f_Hz'}, labs, {'fDecim'}]);
end

function labs = binLabels(v,nB)
%binLabels  A column name per amplitude bin, from its centre; a bin the tree gives
%   no centre for is named by its index rather than left nameless.
labs = cell(1,nB);
pc = [];
if isfield(v,'pctCenters'), pc = double(v.pctCenters(:)); end
for b = 1:nB
    if b<=numel(pc), labs{b} = sprintf('pct_%g',pc(b)); else, labs{b} = sprintf('bin%d',b); end
end
labs = matlab.lang.makeValidName(labs);      % a fractional centre is not a name
labs = matlab.lang.makeUniqueStrings(labs);
end

function y = padRow(x,n)
%padRow  The first n values of a row, NaN-padded - a fixed sheet width.
x = double(x(:)'); y = nan(1,n); m = min(n,numel(x)); y(1:m) = x(1:m);
end

function sink = emitLazy(sink, wantSheet, sheet, builder)
%emitLazy  Build one sheet only when the selection actually wants it.
if ~wantSheet(sheet), return; end
sink = emitSheet(sink, builder(), sheet, wantSheet);
end

function T = intervalsTable(iv)
%intervalsTable  One row per analysed window: WHAT THE WINDOW IS, and the level and
%   the swing of every diameter that was measured inside it.
%   The times are ABSOLUTE seconds from the start of the recording (spec §3), so a
%   window can always be found again in the original footage.  A wire myograph has
%   no diameter and simply leaves those columns empty, as an absent column already
%   is everywhere else here.
%
%   THE 'channel' COLUMN IS THE AXIS, not the operator's assignment: it says which
%   channel this row is about.  Two chambers whose baseline windows are both called
%   'baseline' are two rows that differ only in it, so a sheet without it could not
%   be read at all; 'scope' is the assignment the operator made, and says whether
%   that window belongs to the one channel or to the whole file.
n = numel(iv);
T = table('Size',[n 8], ...
    'VariableTypes',{'string','string','string','double','double','double','double','double'}, ...
    'VariableNames',{'interval','channel','scope','tStart_s','tEnd_s','duration_s','nFrames','nY'});
[meas,stats] = diameterAxes(iv);
for m = 1:numel(meas)
    for q = 1:numel(stats)
        T.([meas{m} '_' stats{q}]) = nan(n,1);
    end
end
for k = 1:n
    d = subStruct(iv(k),'diameter');
    T.interval(k)   = string(intervalName(iv(k),k));
    T.channel(k)    = string(charOrEmpty(fieldOrEmpty(iv(k),'channelName')));
    T.scope(k)      = string(intervalScope(iv(k)));
    T.tStart_s(k)   = fieldNum(iv(k),'tStart');
    T.tEnd_s(k)     = fieldNum(iv(k),'tEnd');
    T.duration_s(k) = T.tEnd_s(k) - T.tStart_s(k);
    T.nFrames(k)    = frameCount(iv(k),d);
    if isfield(d,'nY'), T.nY(k) = double(d.nY); end
    if ~isfield(d,'stats') || ~isstruct(d.stats), continue; end
    for m = 1:numel(meas)
        if ~isfield(d.stats,meas{m}), continue; end
        S = d.stats.(meas{m});
        for q = 1:numel(stats)
            T.([meas{m} '_' stats{q}])(k) = fieldNum(S,stats{q});
        end
    end
end
end

function T = propagationTable(iv)
%propagationTable  One row per interval × analysed measure: the answer (speed and
%   direction), the interval around it, and the evidence the confidence was built
%   from - which is what makes a reported speed something a reader can judge.
C = {};
for k = 1:numel(iv)
    p = subStruct(iv(k),'propagation');
    sigs = fieldnames(p);
    for j = 1:numel(sigs)
        q = p.(sigs{j});
        if ~isstruct(q) || isempty(q), continue; end
        m = subStruct(q,'metrics');
        C(end+1,:) = { intervalName(iv(k),k), sigs{j}, ...
            fieldNum(q,'speed'), fieldStr(q,'speedUnit'), fieldStr(q,'direction'), ...
            ciBound(q,1), ciBound(q,2), fieldNum(q,'R2'), fieldNum(q,'pValue'), ...
            fieldNum(q,'confidence'), fieldStr(q,'confidenceLevel'), ...
            fieldNum(m,'medianCorr'), fieldNum(m,'rowFraction'), ...
            fieldNum(m,'totalLagSamples'), fieldNum(q,'belowResolution'), ...
            fieldNum(q,'domFreq'), fieldNum(q,'nRows'), ...
            strjoin(fieldCellstr(q,'qualityFlags'),'; '), fieldStr(q,'confidenceText') }; %#ok<AGROW>
    end
end
T = cellRows(C,{'interval','measure','speed','speedUnit','direction', ...
    'speedCI_lo','speedCI_hi','R2','pValue','confidence','confidenceLevel', ...
    'medianCorr','rowFraction','totalLagSamples','belowResolution','domFreq', ...
    'nRows','qualityFlags','confidenceText'});
end

function T = vasomotionTable(iv)
%vasomotionTable  One row per interval × signal × unit, carrying the vasomotion
%   markers of both bands.  These are the SAME markers a speckle recording reports
%   per segment - the analysis is the shared core - so a myograph result and an
%   imaging result can be put side by side in one statistics package.  The band is
%   the name suffix: VB is the vasomotion band, CB the comparison band.
mk = vsmMarkerNames();
C  = {};
for k = 1:numel(iv)
    V = subStruct(iv(k),'vasomotion');
    sigs = fieldnames(V);
    for j = 1:numel(sigs)
        v = V.(sigs{j});
        if ~isstruct(v) || ~isfield(v,'scalars'), continue; end
        for u = 1:vsmUnitCount(v)
            row = [{intervalName(iv(k),k), signalName(v,sigs{j}), u}, ...
                   arrayfun(@(b) bandScalar(v,mk{b,2},mk{b,3},u), 1:size(mk,1), 'UniformOutput',false)];
            C(end+1,:) = row; %#ok<AGROW>
        end
    end
end
T = cellRows(C,[{'interval','signal','unit'}, mk(:,1)']);
end

function T = ampPctTable(iv,O)
%ampPctTable  The band-envelope amplitude at each PERCENTILE LEVEL, long: one row
%   per interval × signal × unit × band × level.  scalars.<band>.ampPct is
%   prctile(a(t),levels) of the band-mean envelope, so it says how much of the
%   recording sat at each amplitude rather than only what the average was.
%
%   NOTE the levels here are not the bin CENTRES that key the percentile-resolved
%   SPECTRA: a RESULTS file carries no settings struct, so the levels are recovered
%   from how many of them there are, exactly as guiExplore recovers them.
C = {};
for k = 1:numel(iv)
    if atCap(C,O), break, end
    V = subStruct(iv(k),'vasomotion');
    sigs = fieldnames(V);
    for j = 1:numel(sigs)
        if atCap(C,O), break, end
        v = V.(sigs{j});
        if ~isstruct(v) || ~isfield(v,'scalars'), continue; end
        for bc = {'VB','CB'}
            if atCap(C,O), break, end
            band = bc{1};
            if ~isfield(v.scalars,band) || ~isfield(v.scalars.(band),'ampPct'), continue; end
            P = double(v.scalars.(band).ampPct);          % [nUnit x nLevel]
            if isempty(P), continue; end
            lev = percentileLevels(size(P,2));
            for u = 1:size(P,1)
                if atCap(C,O), break, end
                for q = 1:numel(lev)
                    C(end+1,:) = {intervalName(iv(k),k), signalName(v,sigs{j}), u, ...
                        band, lev(q), P(u,q)}; %#ok<AGROW>
                end
            end
        end
    end
end
T = cellRows(C,{'interval','signal','unit','band','percentile','amp'});
end

function T = spectraTable(iv,O)
%spectraTable  The amplitude spectrum of every interval × signal, averaged over its
%   units and DECIMATED IN FREQUENCY so the sheet stays readable.  The factor is
%   written into the sheet rather than left implicit - a curve nobody can tell the
%   resolution of is not evidence.
C = {};
for k = 1:numel(iv)
    if atCap(C,O), break, end
    V = subStruct(iv(k),'vasomotion');
    sigs = fieldnames(V);
    for j = 1:numel(sigs)
        if atCap(C,O), break, end
        v = V.(sigs{j});
        if ~isstruct(v) || ~isfield(v,'fVectors') || ~isfield(v.fVectors,'ampMean'), continue; end
        f = double(v.f(:)); nf = numel(f);
        if nf==0, continue; end
        am = mean(double(v.fVectors.ampMean),1,'omitnan'); am = am(:);
        as = nan(nf,1);
        if isfield(v.fVectors,'ampStd')
            as = mean(double(v.fVectors.ampStd),1,'omitnan'); as = as(:);
        end
        dec = decimation(O.fDecim, nf, 60);
        for q = 1:dec:nf
            C(end+1,:) = {intervalName(iv(k),k), signalName(v,sigs{j}), f(q), ...
                pick(am,q), pick(as,q), dec}; %#ok<AGROW>
        end
    end
end
T = cellRows(C,{'interval','signal','f_Hz','ampMean','ampStd','fDecim'});
end

function T = diameterTraceTable(iv,O)
%diameterTraceTable  The line-averaged diameter of every measure, DECIMATED IN TIME:
%   time runs down the rows because each interval has its own stretch of the
%   recording and they do not share a time base, and the three measures run across
%   because they are sampled together.
[meas,~] = diameterAxes(iv);
if isempty(meas), T = table(); return; end
C = {};
for k = 1:numel(iv)
    if atCap(C,O), break, end
    d = subStruct(iv(k),'diameter');
    if ~isfield(d,'time') || isempty(d.time), continue; end
    t = double(d.time(:)); nT = numel(t);
    dec = decimation(O.tDecim, nT, 5000);
    for q = 1:dec:nT
        row = [{intervalName(iv(k),k), t(q)}, ...
               cellfun(@(m) traceValue(d,m,q), meas, 'UniformOutput',false), {dec}];
        C(end+1,:) = row; %#ok<AGROW>
    end
end
T = cellRows(C,[{'interval','t_s'}, meas, {'tDecim'}]);
end

% ---------------------------------------------------------------------------
function [meas,stats] = diameterAxes(iv)
%diameterAxes  The diameter MEASURES and the STATISTIC names this recording actually
%   carries, unioned over its intervals in the order the tree lists them.  READ OFF
%   d.stats AND NEVER OFF fieldnames(d), which is what keeps the diameter branch's
%   arrays - the walls, the per-line diameters, the quality flags - out of a sheet:
%   a sheet is what a protocol is compared on, and a [frames x lines] array is not
%   that.  Adding a statistic adds a column; adding an array adds nothing.
meas = {}; stats = {};
for k = 1:numel(iv)
    d = subStruct(iv(k),'diameter');
    if ~isfield(d,'stats') || ~isstruct(d.stats), continue; end
    fn = fieldnames(d.stats)';
    meas = [meas, fn(~ismember(fn,meas))]; %#ok<AGROW>
    for m = 1:numel(fn)
        sn = fieldnames(d.stats.(fn{m}))';
        stats = [stats, sn(~ismember(sn,stats))]; %#ok<AGROW>
    end
end
end

function n = vsmUnitCount(v)
%vsmUnitCount  How many results one signal's tree holds: one for a line-averaged
%   trace or a channel, one per image row when the vasomotion ran per line.
n = 0;
if isfield(v.scalars,'VB') && isfield(v.scalars.VB,'ampMean')
    n = numel(v.scalars.VB.ampMean);
elseif isfield(v.scalars,'CB') && isfield(v.scalars.CB,'ampMean')
    n = numel(v.scalars.CB.ampMean);
end
end

function y = bandScalar(v,band,name,u)
%bandScalar  One per-unit band marker, NaN when this tree does not carry it.
y = NaN;
if isfield(v.scalars,band) && isfield(v.scalars.(band),name)
    d = double(v.scalars.(band).(name));
    if numel(d)>=u, y = d(u); end
end
end

function mk = vsmMarkerNames()
%vsmMarkerNames  The vasomotion markers a flat sheet carries, and where each lives
%   in the band-branched tree: {column name, band, field}.  The frequency-shape and
%   multiplicity markers exist for the vasomotion band only; the flare / silence
%   clustering exists for both, each from its own threshold.
mk = { 'ampMeanVB','VB','ampMean';   'ampStdVB','VB','ampStd';   'ampSkewVB','VB','ampSkew'; ...
       'ampMeanCB','CB','ampMean';   'ampStdCB','CB','ampStd';   'ampSkewCB','CB','ampSkew'; ...
       'fCentMeanVB','VB','fCentMean';   'fCentStdVB','VB','fCentStd'; ...
       'fSprdMeanVB','VB','fSprdMean';   'fSprdStdVB','VB','fSprdStd'; ...
       'shapePeakVB','VB','shapePeak';   'nPeakMeanVB','VB','nPeakMean';   'nPeakStdVB','VB','nPeakStd'; ...
       'durFlareMeanVB','VB','durFlareMean';       'durFlareStdVB','VB','durFlareStd'; ...
       'durSilenceMeanVB','VB','durSilenceMean';   'durSilenceStdVB','VB','durSilenceStd'; ...
       'ampFlareMeanVB','VB','ampFlareMean';       'ampFlareStdVB','VB','ampFlareStd'; ...
       'ampSilenceMeanVB','VB','ampSilenceMean';   'ampSilenceStdVB','VB','ampSilenceStd'; ...
       'durFlareMeanCB','CB','durFlareMean';       'durFlareStdCB','CB','durFlareStd'; ...
       'durSilenceMeanCB','CB','durSilenceMean';   'durSilenceStdCB','CB','durSilenceStd'; ...
       'ampFlareMeanCB','CB','ampFlareMean';       'ampFlareStdCB','CB','ampFlareStd'; ...
       'ampSilenceMeanCB','CB','ampSilenceMean';   'ampSilenceStdCB','CB','ampSilenceStd' };
end

function lev = percentileLevels(nP)
%percentileLevels  The percentile LEVELS ampPct was evaluated at.  A RESULTS file
%   carries no settings struct, so they are recovered from their count - which is
%   right for the default 0:10:100 and for any other evenly spaced choice, and is
%   the same fallback guiExplore uses.
lev = linspace(0,100,max(nP,1));
end

function tf = atCap(C,O)
%atCap  Has this sheet reached its row ceiling?  Tested at the top of EVERY loop
%   level rather than only the innermost, so accumulation actually stops instead of
%   the outer loops carrying on past the cap.
tf = size(C,1) >= O.maxRows;
end

function d = decimation(requested, n, target)
%decimation  Keep every dth point: what the caller asked for, or enough to bring a
%   curve down to about `target` points.
if ~isempty(requested), d = max(1,round(double(requested(1)))); return; end
d = max(1,ceil(n/target));
end

function v = traceValue(d,measure,q)
%traceValue  One sample of one line-averaged diameter trace.
v = NaN;
if isfield(d,measure)
    x = double(d.(measure));
    if numel(x)>=q, v = x(q); end
end
end

function T = cellRows(C, names)
%cellRows  A table from accumulated rows, or an EMPTY table when nothing was
%   accumulated - which is how a sheet this file has no data for gets skipped,
%   exactly like an absent sheet.
if isempty(C), T = table(); return; end
T = cell2table(C,'VariableNames',names);
end

function nm = intervalName(ivk,k)
%intervalName  What the window is called, or what it is when nobody named it.
nm = '';
if isfield(ivk,'name'), nm = char(string(ivk.name)); end
if isempty(nm), nm = sprintf('interval %d',k); end
end

function s = intervalScope(ivk)
%intervalScope  What the operator assigned this window to.  A wire myograph window
%   left on every channel is repeated under each of them, and telling that apart
%   from four windows that were each assigned to one chamber is the difference
%   between one decision and four.
c = fieldCellstr(ivk,'channels');
if isempty(c), s = 'all channels'; else, s = strjoin(c,'; '); end
if isempty(charOrEmpty(fieldOrEmpty(ivk,'channelName'))), s = ''; end
end

function nm = signalName(v,fld)
%signalName  What was analysed.  A channel's real name is kept beside its tree
%   because 'Force 1 (mN)' does not survive being made into a field name; a diameter
%   measure is its own field name already.
nm = fld;
if isfield(v,'channelName') && ~isempty(v.channelName), nm = char(string(v.channelName)); end
end

function n = frameCount(ivk,d)
%frameCount  How many samples the window holds, asked of the three places that can
%   answer, in this order: its range into the ORIGINAL recording, the diameter
%   measured in it, and - for a wire recording - the samples it kept.
%
%   THE THIRD IS NOT REDUNDANT.  A wire window carries a frame range only when every
%   channel of the file shares one sampling rate: with two rates there is no single
%   sample index that would be true for both, and a wrong one would be worse than
%   none, so cutMyographIntervals leaves it empty.  Such a window has no diameter
%   branch either - a wire myograph measures none - so without this last fallback a
%   mixed-rate recording would report an empty count for every window, which is the
%   one shape of file where the count is least guessable by hand.
n = NaN;
fr = double(fieldOrEmpty(ivk,'frames'));
if numel(fr) == 2
    n = fr(2) - fr(1) + 1; return
end
if isfield(d,'time') && ~isempty(d.time)
    n = numel(d.time); return
end
sm = subStruct(ivk,'samples');
if isfield(sm,'time') && ~isempty(sm.time)
    n = numel(sm.time);
end
end

function S = subStruct(st,f)
%subStruct  One branch of the interval tree as a struct, empty when it was never
%   written - so a caller can ask for its fields without testing twice.
S = struct();
if isstruct(st) && isfield(st,f) && isstruct(st.(f)) && isscalar(st.(f)), S = st.(f); end
end

function v = fieldNum(st,f)
v = NaN;
if isstruct(st) && isfield(st,f) && ~isempty(st.(f)) && (isnumeric(st.(f))||islogical(st.(f)))
    x = double(st.(f)); v = x(1);
end
end

function v = fieldStr(st,f)
v = '';
if isstruct(st) && isfield(st,f) && (ischar(st.(f))||isstring(st.(f))), v = char(st.(f)); end
end

function v = fieldCellstr(st,f)
v = {};
if isstruct(st) && isfield(st,f) && iscell(st.(f)), v = cellfun(@(x)char(string(x)),st.(f)(:)','UniformOutput',false); end
end

function v = ciBound(q,i)
v = NaN;
if isfield(q,'speedCI') && numel(q.speedCI)>=i, v = double(q.speedCI(i)); end
end

function v = pick(x,q)
if numel(x)>=q, v = x(q); else, v = NaN; end
end

%% ===================== shared machinery ================================ %%
function Ti = idTable(M)
% Identifier columns for a pulsatility sheet: idx plus type/label when present.
% Row order is the segment order shared by results.sMetrics / dvsMetrics and the
% results.pulsatility.<sig> tree (DATA-MODEL segment-order invariant).
Ti = table(double(M.idx(:)),'VariableNames',{'idx'});
if ismember('type', M.Properties.VariableNames),  Ti.type  = string(M.type(:));  end
if ismember('label',M.Properties.VariableNames),  Ti.label = string(M.label(:)); end
end

function Tp = pulsScalarTable(S)
% Flatten one results.pulsatility.<sig>.scalars struct (per-segment) into a table:
% each [nSeg x 1] scalar becomes a column under its own name; the harmonic
% coefficients hAmp/hPhase [nSeg x nHarm] expand to <prefix>hAmp1..N /
% <prefix>hPhase1..N.  Non-[nSeg x 1] extras are skipped; returns an empty table when
% nothing usable is present.
%
% THE PREFIX IS NOT PASSED IN, IT IS READ OFF THE STRUCT, and that is the difference
% between this and the version that took it as an argument.  Every marker the producers
% write already carries it - psPI is a pulsatile FLOW, pdPI a pulsatile DIAMETER, pvPI a
% pulsatile PLASMA VOLUME, and no two of the three may be pooled - so a literal here is
% a second place for the same fact to be written down and a second way for it to be
% wrong.  It was wrong: the dvsPulsatility sheet passed 'ps' unconditionally, which on a
% fluorescence product put a column named for a pulsatile flow on a sheet whose every
% other column is a plasma volume.
%
% AND hAmp / hPhase ARE WHY THIS MATTERS ON THE OTHER SHEET TOO.  They are the only two
% leaves the producers leave unprefixed, so a MERGED workbook - one row per file - would
% otherwise stack a speckle recording's flow harmonics and a fluorescence recording's
% plasma-volume harmonics into one column called hAmp1.  Tagging them with the prefix
% their siblings carry is what keeps the two apart, and it is why the tag is not
% conditional on the sheet holding two signals.
Tp = table();
if ~isstruct(S), return; end
harmPrefix = harmonicPrefix(S);
fn = fieldnames(S); cols = {}; names = {};
for i = 1:numel(fn)
    v = S.(fn{i});
    if isempty(v), continue; end
    if any(strcmp(fn{i},{'hAmp','hPhase'}))
        for k = 1:size(v,2)
            cols{end+1}  = double(v(:,k));                        %#ok<AGROW>
            names{end+1} = sprintf('%s%s%d',harmPrefix,fn{i},k);  %#ok<AGROW>
        end
    elseif iscolumn(v)
        cols{end+1}  = double(v(:));                              %#ok<AGROW>
        names{end+1} = fn{i};                                     %#ok<AGROW>
    end
end
if ~isempty(cols)
    Tp = array2table([cols{:}],'VariableNames',names);
end
end

function pre = harmonicPrefix(S)
%harmonicPrefix  The two-letter quantity prefix this signal's OWN markers carry - 'ps'
%   for a pulsatile flow, 'pd' for a pulsatile diameter, 'pv' for a pulsatile plasma
%   volume.  '' when there is nothing to read it off, which is the only case in which
%   nothing is claimed: a struct holding harmonics and no markers at all, and one
%   holding two prefixes at once, which no producer writes.
pre = '';
t = regexp(fieldnames(S), '^(p[a-z])[A-Z]', 'tokens', 'once');
t = t(~cellfun(@isempty, t));
if isempty(t), return; end
u = unique(cellfun(@(c) c{1}, t, 'UniformOutput', false));
if isscalar(u), pre = u{1}; end
end

function m = measuredRows(T)
%measuredRows  WHICH ROWS OF A METRICS TABLE ARE A REAL, MEASURED SEGMENT, as a row
%   logical over the table's height.  A row whose idx is missing or zero is an empty
%   table slot; and on a recording that carries a flow index, a row whose index was
%   never computed is not a segment anybody measured.
%
%   THE FLOW INDEX IS ASKED FOR RATHER THAN NAMED, and that is what lets one export take
%   a mixed working set.  The fluorescence branch writes NO BFI column at all - there is
%   nothing to invert, an intravascular tracer already makes the intensity proportional
%   to the plasma volume - so naming it here refused every intensity product outright,
%   with an error that stopped the whole batch rather than skipping one file.  Asking
%   the table what it holds is the same statement this file already makes about
%   dvsMetrics and about every absent sheet.
idx = T.idx(:)';
m   = ~isnan(idx) & idx~=0;
if ismember('BFI', T.Properties.VariableNames)
    m = m & ~isnan(T.BFI(:)');
end
end

function T = wallMotionTable(results)
%wallMotionTable  ONE ROW PER SEGMENT: the wall-motion measurement BEFORE its gates, the
%   floor it was judged against, and which gate refused it.
%
%   THE ROW ORDER IS THE METRICS TABLE'S, and the producer asserts it (wallMotion.idx ==
%   sMetrics.idx), so the join is by row and no key column is needed beyond the idx that
%   every metric sheet already carries.
%
%   THE COLUMNS COME FROM THE TREE AND NEVER FROM A LIST HERE, exactly as the NVC sheets'
%   do: a field that is one value per row is a column, and a field that is not - the
%   pixel size, how many matched controls each vessel was compared with, the two band
%   edges of a continuous run - describes the whole measurement rather than a vessel.
%   The one place that rule is ambiguous is a recording with a single segment, where a
%   scalar and a column are the same shape; that is not a recording anyone analyses.
%
%   THE REASON COLUMN IS KEPT AND IT IS TEXT.  'why' is the producer's own sentence for a
%   segment that was never offered a cut at all, and a sheet that dropped it would leave
%   a reader with a row of NaNs and nothing to read.
T = table();
if ~isfield(results,'wallMotion') || ~isstruct(results.wallMotion), return; end
if ~isfield(results,'sMetrics')   || ~istable(results.sMetrics),    return; end
W = results.wallMotion;
n = height(results.sMetrics);
if n==0, return; end

T   = idTable(results.sMetrics);
n0  = width(T);
fn  = fieldnames(W);
for i = 1:numel(fn)
    if strcmp(fn{i},'idx'), continue; end            % idTable already carries it
    v = W.(fn{i});
    if (isnumeric(v) || islogical(v)) && isvector(v) && numel(v)==n
        T.(fn{i}) = double(v(:));
    elseif iscell(v) && numel(v)==n && all(cellfun(@(x) ischar(x)||isstring(x), v(:)))
        T.(fn{i}) = string(v(:));
    end
end
if width(T)<=n0, T = table(); end                    % nothing per-row: no sheet
end

function T = topologyTable(results)
%topologyTable  results.topology.metrics as the step wrote it - ONE ROW PER ANALYSED
%   AREA, the whole segmented field first and then one row per drawn region - with the
%   units the numbers are in written onto every row.
%
%   THE UNITS ARE PART OF THE SHEET AND NOT A NOTE SOMEWHERE.  The same column is
%   micrometres on a calibrated recording and pixels on one whose micrometres per pixel
%   nobody supplied, the choice is made per recording, and a merged workbook stacks both
%   kinds - so a sheet that did not carry them would be a column of numbers that cannot
%   be read.
T = table();
if ~isfield(results,'topology') || ~isstruct(results.topology), return; end
if ~isfield(results.topology,'metrics') || ~istable(results.topology.metrics), return; end
T = results.topology.metrics;
u = fieldOrEmpty(results.topology,'units');
if ~isstruct(u), return; end
for f = {'lengthUnit','areaUnit','densityUnit'}
    if isfield(u,f{1}) && (ischar(u.(f{1})) || isstring(u.(f{1})))
        T.(f{1}) = repmat(string(u.(f{1})), height(T), 1);
    end
end
end

function [Tm, Td] = aggregateSegments(results, keyVar, otherVar, weightByArea)
%aggregateSegments  THE per-segment aggregation, one axis at a time.
%   Groups the valid, named segments of results.sMetrics by keyVar (the LABEL for
%   the historical ROI sheets, the VESSEL TYPE for the new ones) and returns
%     Tm  the metric table  - area and length SUMMED, every other numeric column
%                             averaged with the given weights, plus the modal
%                             value of otherVar (the companion axis) per group;
%     Td  the trace table   - Time + one weighted-mean time series per group.
%   Weights are results.sMetrics.area normalised within the group (weightByArea
%   true, the library's convention and what the ROI sheets have always used) or a
%   plain 1/n (false).  Both tables come back EMPTY when the file carries no
%   keyVar column, so the caller simply skips the sheets - exactly how an absent
%   sheet is already handled.  Called with ('label','type',true) this is the
%   legacy sMetricsROI / sDataROI code path, byte for byte.
Tm = table(); Td = table();
V = results.sMetrics.Properties.VariableNames;
if ~ismember(keyVar, V), return; end

goodMask = measuredRows(results.sMetrics);

keep = goodMask' & strlength(results.sMetrics.(keyVar)) > 0  ...
    & ~ismember(results.sMetrics.category, ["outerWall","innerWall","unsegmented"]);

M        = results.sMetrics(keep ,:);
areaAll  = M.area;                         % weights   (unchanged)
lenAll   = M.length;                       % new       (to be summed)
if ismember('RegOverlap', M.Properties.VariableNames)
    overlap   = M.RegOverlap;                       % new       (to be summed)
end
[G,labelList] = findgroups(M.(keyVar));
nG = numel(labelList);

% numeric variables to weight-average   (skip idx, area, length)
isNum   = varfun(@isnumeric,M,'OutputFormat','uniform');
numVar  = M.Properties.VariableNames(isNum);
numVar(ismember(numVar,{'idx','area','length','RegOverlap','closest vessel idx','nearestVesIdx','RegID'})) = [];

hasOther = ismember(otherVar, M.Properties.VariableNames);
if hasOther, keyNames = {keyVar, otherVar};
else,        keyNames = {keyVar};
end

% pre-allocate output: additional 'length' column
T = table('Size',[nG numel(numVar)+2+numel(keyNames)], ... % +area +length +key(+other)
    'VariableTypes',[repmat("double",1,numel(numVar)+2) repmat("string",1,numel(keyNames))], ...
    'VariableNames',[numVar {'area','length'} keyNames]);

for g = 1:nG
    rows = (G==g);
    areaSum        = sum(areaAll(rows));          % -------- SUM area
    lenSum         = sum(lenAll(rows));           % -------- SUM length


    T.area(g)      = areaSum;
    T.length(g)    = lenSum;

    if ismember('RegOverlap', M.Properties.VariableNames)
        overlapSum=sum(overlap(rows));
        T.RegOverlap(g)    = overlapSum;
    end

    w              = groupWeights(areaAll(rows), areaSum, weightByArea);
    for v = 1:numel(numVar)                       % weighted means
        col        = numVar{v};
        T{g,col}   = sum( M.(col)(rows) .* w );
    end

    if hasOther
        T.(otherVar)(g) = string( mode( categorical(M.(otherVar)(rows)) ) );
    end
    T.(keyVar)(g) = labelList(g);
end


Tm = movevars(T,keyNames,'Before',1);              % final layout

if ~isfield(results,'sData'), return; end

timeVec  = results.time(:);
sDataUse = results.sData(:,keep);
sigAgg   = zeros(numel(timeVec), nG);

for g = 1:nG
    rows = (G==g);
    w    = groupWeights(areaAll(rows), sum(areaAll(rows)), weightByArea);  % same weights
    sigAgg(:,g) = sDataUse(:,rows) * w;                 % N×k × k×1 = N×1
end

sigNames = matlab.lang.makeValidName( compose("%s", labelList) );

Td = [ table(timeVec,'VariableNames',{'Time'}), ...
    array2table(sigAgg,'VariableNames',sigNames) ];
end

function w = groupWeights(areaRows, areaSum, weightByArea)
%groupWeights  Normalised within-group weights: area-proportional (the ROI-sheet
%   convention) or a plain unweighted mean.
if weightByArea
    w = areaRows ./ areaSum;                            % normalised weights
else
    w = ones(numel(areaRows),1) ./ numel(areaRows);
end
end

function w = axisWeight(O, axisName)
%axisWeight  opts.weightByArea applies to the .groupBy AXIS ONLY - any other
%   aggregation in the same call keeps its historical area weighting, so asking
%   for an unweighted vessel-type average can never silently change the ROI sheets.
w = true;
if strcmp(O.groupBy, axisName), w = O.weightByArea; end
end

function T = pickColumns(T, cols)
%pickColumns  Keep only the requested sMetrics columns, in the table's own order.
%   Absent request = keep everything.  A requested column the file does not have
%   is simply not there - the same "skip what is absent" rule the sheets follow.
if isempty(cols), return; end
v = T.Properties.VariableNames;
T = T(:, v(ismember(v, cols)));
end

function [O] = parseExportOpts(opts, fNames)
%parseExportOpts  Normalise the optional selection struct into ONE contract:
%   what to write, how to average it, which columns, and where it goes.  Absent /
%   empty opts reproduce today's full behaviour.
%     opts.sheets       names to write (cellstr / string / char).  [] or absent =
%                       the default set for .groupBy.
%     opts.format       output workbook extension.  '.xlsx' (default) or '.xls' -
%                       both multi-sheet spreadsheet formats writetable supports.
%     opts.groupBy      '' | 'label' | 'type'   (row granularity + default sheets)
%     opts.weightByArea logical, default true   (applies to the groupBy axis)
%     opts.columns      cellstr subset of the sMetrics variable names
%     opts.merge        logical, default false  (one workbook for the whole call)
%     opts.labels       1xN struct array parallel to fNames (.animal/.type/.group)
%     opts.outFile      destination of the merged workbook
%     opts.fDecim       myograph: keep every Nth frequency point ([] = automatic)
%     opts.tDecim       myograph: keep every Nth time sample    ([] = automatic)
%     opts.maxRows      myograph: row ceiling per sheet (default 1e6)
O = struct('sheets',{{}},'ext','.xlsx','groupBy','','weightByArea',true, ...
    'columns',{{}},'merge',false,'labels',[],'outFile','', ...
    'fDecim',[],'tDecim',[],'maxRows',1e6);
if isempty(opts) || ~isstruct(opts), O.sheets = defaultSheets(O.groupBy); return; end

if isfield(opts,'format') && ~isempty(opts.format)
    e = lower(char(string(opts.format)));
    if e(1) ~= '.', e = ['.' e]; end
    if ~any(strcmp(e,{'.xlsx','.xls'}))
        error('exportToExcel:format', ...
            'opts.format must be ''.xlsx'' or ''.xls'' (multi-sheet spreadsheet); got ''%s''.', e);
    end
    O.ext = e;
end

if isfield(opts,'groupBy') && ~isempty(opts.groupBy)
    g = char(string(opts.groupBy));
    if ~any(strcmp(g,{'label','type'}))
        error('exportToExcel:groupBy', ...
            'opts.groupBy must be '''', ''label'' or ''type''; got ''%s''.', g);
    end
    O.groupBy = g;
end

if isfield(opts,'weightByArea') && ~isempty(opts.weightByArea)
    O.weightByArea = logical(opts.weightByArea(1));
end

% the myograph decimation controls: [] keeps the automatic choice, which is what a
% caller that has never heard of them gets
if isfield(opts,'fDecim') && ~isempty(opts.fDecim), O.fDecim = double(opts.fDecim(1)); end
if isfield(opts,'tDecim') && ~isempty(opts.tDecim), O.tDecim = double(opts.tDecim(1)); end
if isfield(opts,'maxRows') && ~isempty(opts.maxRows)
    O.maxRows = max(1,double(opts.maxRows(1)));
end

if isfield(opts,'columns') && ~isempty(opts.columns)
    O.columns = cellstr(string(opts.columns));
    O.columns = O.columns(:)';
end

if isfield(opts,'sheets') && ~isempty(opts.sheets)
    O.sheets = cellstr(string(opts.sheets));
    O.sheets = O.sheets(:)';
else
    O.sheets = defaultSheets(O.groupBy);
end

if isfield(opts,'merge') && ~isempty(opts.merge)
    O.merge = logical(opts.merge(1));
end

if isfield(opts,'labels') && ~isempty(opts.labels)
    if ~isstruct(opts.labels)
        error('exportToExcel:labels','opts.labels must be a struct array parallel to fNames.');
    end
    if numel(opts.labels) ~= numel(fNames)
        error('exportToExcel:labels', ...
            'opts.labels has %d entries but fNames has %d - they must be parallel.', ...
            numel(opts.labels), numel(fNames));
    end
    O.labels = opts.labels(:)';
end

if isfield(opts,'outFile') && ~isempty(opts.outFile)
    O.outFile = char(string(opts.outFile));
    [~,~,oe] = fileparts(O.outFile);
    if isempty(oe)
        O.outFile = [O.outFile O.ext];
    elseif ~any(strcmpi(oe,{'.xlsx','.xls'}))
        error('exportToExcel:outFile', ...
            'opts.outFile must be an .xlsx / .xls workbook; got ''%s''.', oe);
    end
end
if O.merge && isempty(O.outFile)
    O.outFile = defaultMergedName(fNames, O.ext);
end
end

function names = defaultSheets(groupBy)
%defaultSheets  What .groupBy asks for when opts.sheets is not given: EVERYTHING a
%   file can be about for per-recording rows, and the matching aggregated PAIR for
%   an averaging axis.  The two vessel-type sheets are therefore OPT-IN.
%
%   THE DEFAULT IS THE UNION OF BOTH SHEET SETS, and it costs the segmented path
%   nothing: a segmented recording has no intervals, so the nine myograph names it
%   now also asks for are skipped for it exactly as an absent sheet always was, and
%   a bare exportToExcel(fNames) still writes the historical nine and nothing else.
%   What it buys is that the same bare call on a myograph recording writes the
%   myograph workbook, instead of an empty one - the selection is the union and the
%   DATA decides which half of it appears.
switch groupBy
    case 'label', names = {'sMetricsROI','sDataROI'};
    case 'type',  names = {'sMetricsType','sDataType'};
    otherwise,    names = [legacySheets() myographSheets()];
end
end

function names = legacySheets()
%legacySheets  The sheets a bare exportToExcel(fNames) writes for a segmented
%   recording, in write order.  THE definition of "everything" for one.  The two NVC
%   sheets are the newest and, like every other, they are simply absent from a file
%   the step was never run on.
names = {'sMetrics','sData','sMetricsROI','sDataROI','dvsMetrics','dvsData', ...
         'dvsDiameter','pulsatility','dvsPulsatility','wallMotion','topology', ...
         'nvcTrials','dvsNvcTrials'};
end

function names = myographSheets()
%myographSheets  The nine sheets a myograph recording writes, in write order.  Only a
%   WIRE myograph carries an operator log, so 'comments' is the one of them that a
%   pressure recording routinely skips.
names = {'settings','comments','intervals','propagation','vasomotion','ampPct', ...
         'spectra','ampPctSpectra','diameterTraces'};
end

function out = defaultMergedName(fNames, ext)
%defaultMergedName  Where a merged workbook lands when the caller names no file:
%   beside the first real input, so it is found next to the data it summarises.
out = ['mergedExport' ext];
for i = 1:numel(fNames)
    if ~isempty(fNames{i})
        out = fullfile(fileparts(fNames{i}), ['mergedExport' ext]);
        return
    end
end
end

function sink = newSink(O)
%newSink  The writer.  Unmerged it is a thin wrapper over writetable, one target
%   file at a time; merged it accumulates every sheet across the whole call and
%   writes once at the end (the only way row counts and column unions can be
%   right).  Either way the FIRST sheet written to a workbook uses 'replacefile'.
sink = struct('merge',O.merge,'outFile',O.outFile,'labels',O.labels, ...
    'target','','wroteAny',false,'prefix',table(),'order',{{}}, ...
    'data',containers.Map('KeyType','char','ValueType','any'));
end

function sink = openFile(sink, outName, inName, fidx)
%openFile  Start a new input file: point the writer at its workbook (unmerged) or
%   build the row prefix its rows will carry (merged).
if sink.merge
    sink.prefix = prefixRow(inName, sink.labels, fidx);
else
    sink.target   = outName;
    sink.wroteAny = false;                          % first written sheet -> replacefile
end
end

function P = prefixRow(inName, labels, fidx)
%prefixRow  The leading identity columns of a merged row: the recording name
%   always, plus the three comparison labels when the caller supplied them.  The
%   'type' here is the RECORDING type (spec §2); a per-segment VESSEL type column
%   is renamed 'vesselType' when the two meet (addPrefix).
[~,nm,ex] = fileparts(inName);
P = table(string(regexprep([nm ex],'_[drs]\.mat$','')),'VariableNames',{'file'});
if isempty(labels) || fidx > numel(labels), return; end
L = labels(fidx);
P.animal = string(labelField(L,{'animal'}));
P.type   = string(labelField(L,{'type'}));
P.group  = string(labelField(L,{'group','expGroup'}));
end

function v = labelField(L, names)
%labelField  First present field of a per-file label struct, as char ('' if none).
v = '';
for i = 1:numel(names)
    if isfield(L,names{i}) && ~isempty(L.(names{i}))
        v = char(string(L.(names{i}))); return
    end
end
end

function T = addPrefix(T, P)
%addPrefix  Stamp the merged row prefix onto one sheet, renaming any data column
%   that would collide with it - a per-segment 'type' becomes 'vesselType', so a
%   merged workbook never confuses the vessel type with the recording type.
if isempty(P) || width(P)==0, return; end
v  = T.Properties.VariableNames;
pv = P.Properties.VariableNames;
for i = 1:numel(v)
    if ~any(strcmp(v{i}, pv)), continue; end
    if strcmp(v{i},'type') && ~any(strcmp('vesselType',v))
        v{i} = 'vesselType';
    else
        v{i} = matlab.lang.makeUniqueStrings([v{i} '_seg'], [pv v]);
    end
end
T.Properties.VariableNames = v;
T = [repmat(P, height(T), 1), T];
end

function sink = emitSheet(sink, T, sheet, wantSheet)
%emitSheet  Write (or accumulate) one SELECTED sheet.  The first sheet actually
%   written to a workbook uses 'replacefile' (a fresh file); later sheets append.
%   A selection that drops sMetrics therefore still starts from a clean file, and
%   the default (all-sheets) path is byte-identical to the historical "sMetrics
%   with 'replacefile' first, the rest appended" behaviour.  A table left with no
%   columns (opts.columns matched nothing in this file) is skipped, like an absent
%   sheet.  Merged, the sheet is stamped with the file's prefix and stacked onto
%   what earlier files contributed.
if ~wantSheet(sheet), return; end
if width(T)==0, return; end
if sink.merge
    T = addPrefix(T, sink.prefix);
    if isKey(sink.data, sheet)
        T = stackTables(sink.data(sheet), T);
    else
        sink.order{end+1} = sheet;
    end
    sink.data(sheet) = T;
    return
end
if ~sink.wroteAny
    writetable(T, sink.target, 'Sheet', sheet, 'WriteMode','replacefile');
else
    writetable(T, sink.target, 'Sheet', sheet);
end
sink.wroteAny = true;
end

function finishSink(sink)
%finishSink  Merged mode: write the single accumulated workbook, sheets in the
%   order they first appeared (the same order an unmerged export writes them).
if ~sink.merge || isempty(sink.order), return; end
first = true;
for i = 1:numel(sink.order)
    T = sink.data(sink.order{i});
    if width(T)==0, continue; end
    if first
        writetable(T, sink.outFile, 'Sheet', sink.order{i}, 'WriteMode','replacefile');
        first = false;
    else
        writetable(T, sink.outFile, 'Sheet', sink.order{i});
    end
end
end

function T = stackTables(A, B)
%stackTables  Vertically concatenate two sheets that need not share columns: the
%   union is taken and a file missing a column contributes empty cells there.
%   Without it one file lacking psPI would drop every other file's psPI, which is
%   exactly the "skip what is absent, per file" rule the sheets already follow.
va = A.Properties.VariableNames;
vb = B.Properties.VariableNames;
addToB = setdiff(va, vb, 'stable');
addToA = setdiff(vb, va, 'stable');
for i = 1:numel(addToB), B.(addToB{i}) = fillLike(A.(addToB{i}), height(B)); end
for i = 1:numel(addToA), A.(addToA{i}) = fillLike(B.(addToA{i}), height(A)); end
order = [va, addToA];
T = [A(:,order); B(:,order)];
end

function col = fillLike(proto, n)
%fillLike  An n-row "absent" column of the same kind as proto.
w = size(proto,2);
if isnumeric(proto),        col = NaN(n,w);
elseif islogical(proto),    col = false(n,w);
elseif isstring(proto),     col = strings(n,w);
elseif iscategorical(proto),col = repmat(categorical(missing),n,w);
elseif iscell(proto),       col = repmat({''},n,w);
else,                       col = NaN(n,w);
end
end
