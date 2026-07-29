%exportToExcel  Dump processed LSCI metrics and traces to an .xlsx workbook
%
%   exportToExcel(fNames) opens each processed *.mat dataset in fNames
%   (contrast-derived *_BFI_d / _BFI_r / _BFI_s or later), converts numeric
%   category codes to readable text, builds several summary tables, and
%   writes them as separate sheets in a single Excel file whose name mirrors
%   the input file (e.g.  Foo_BFI_d.mat → Foo_BFI.xlsx):
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
%       pulsatility       (if present) per-segment pulsatility tree scalars for
%                         the flow signal (markers psMin..psSymRatio/psR2 +
%                         harmonic hAmp/hPhase), aligned to sMetrics rows by idx
%       dvsPulsatility    (if present) per dynamic vessel: flow (ps*) + diameter
%                         (pd*) pulsatility scalars, aligned to dvsMetrics rows
%
%   NOTE  runPulsatility is now non-destructive: the sData / dvsData / dvsDiameter
%   trace sheets carry the RAW averaged cycle (its harmonic-fit reconstruction
%   lives in results.pulsatility, not in these traces); sMetrics / dvsMetrics gain
%   the ps*/pd* pulsatility columns.  The per-pixel maps results.pulsatility.ppx
%   are image data and are NOT exported here.
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
%                         {sMetrics, sData, sMetricsROI, sDataROI, sMetricsType,
%                          sDataType, dvsMetrics, dvsData, dvsDiameter,
%                          pulsatility, dvsPulsatility}.
%                         Absent / empty = the default set for .groupBy (below).
%                         A selected sheet whose data is absent from a given file
%                         is simply skipped, exactly as in the full path.
%                .format  output extension, '.xlsx' (default) or '.xls'.
%                .groupBy the ROW GRANULARITY of the metric sheets, and with it
%                         the DEFAULT sheet set:
%                           ''       (default) per-segment rows - all 9 legacy
%                                    sheets, including the label-averaged ROI pair;
%                           'label'  average over results.sMetrics.label - the
%                                    sMetricsROI / sDataROI pair alone (the SAME
%                                    aggregation the legacy path writes, not a
%                                    second implementation);
%                           'type'   average over results.sMetrics.type, the
%                                    VESSEL type written by setVesselTypes - the
%                                    sMetricsType / sDataType pair alone.
%                         .sheets, when given, still selects freely across all
%                         eleven names.
%                .weightByArea  logical, TRUE by default: aggregated rows are
%                         weighted by results.sMetrics.area, the weights the ROI
%                         sheets have always used.  FALSE gives a plain unweighted
%                         mean.  It applies to the .groupBy AXIS ONLY - every other
%                         aggregation in the same call keeps its historical area
%                         weighting, so a legacy sheet can never change silently.
%                .columns cellstr subset of the sMetrics VARIABLE NAMES to export.
%                         Absent = all.  Taken LITERALLY (include 'idx' / 'label' /
%                         'type' yourself if you want them) and applied to the
%                         sMetrics / sMetricsROI / sMetricsType sheets only - the
%                         trace and dvs sheets are untouched.  A column absent from
%                         a given file is skipped for that file, exactly as an
%                         absent sheet already is; a sheet left with no column at
%                         all is skipped for that file.
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
%     files = dir(fullfile(dataRoot,'*BFI_d.mat'));
%     fNames = fullfile({files.folder}',{files.name}');
%     exportToExcel(fNames);                              % full legacy workbook
%     exportToExcel(fNames, struct('sheets',{{'sMetrics','pulsatility'}}));
%                                                         % just those two sheets
%     exportToExcel(fNames, struct('groupBy','type','weightByArea',false));
%                                                         % unweighted per-vessel-type means
%     exportToExcel(fNames, struct('merge',true,'labels',L,'outFile',out));
%                                                         % one workbook for statistics
%
%   DEPENDS ON
%     MATLAB R2019b+ (for writetable with 'Sheet' option) and data schema
%     used throughout the Dynamic Light Scattering Imaging toolbox.
%
% See also: guiExport, setVesselTypes, guiWorkbench, wbSession
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

function exportToExcel(fNames, opts)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".mat".');
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
    load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

    fName=strrep(fName,'_d.mat',O.ext);
    sink=openFile(sink,fName,fNames{fidx},fidx);        % target file / merged row prefix

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
    goodMask  = ~isnan(idxAll) & idxAll~=0 & ~isnan(results.sMetrics.BFI(:)');
    roiIdx    = idxAll(goodMask);
    roiData   = results.sData(:,goodMask);
    roiNames  = matlab.lang.makeValidName(compose("ROI %04d",roiIdx) );
    T = [ table(results.time(:),'VariableNames', {'Time'}), array2table(roiData,'VariableNames',roiNames) ];

    sink=emitSheet(sink,T,'sData',wantSheet);

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
        goodMask = ~isnan(idxAll) & idxAll~=0 & ~isnan(results.dvsMetrics.BFI(:)');
        roiIdx   = idxAll(goodMask);
        roiNames = matlab.lang.makeValidName(compose("ROI %04d", roiIdx));

        T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
            array2table(results.dvsData(:,goodMask), ...
            'VariableNames',roiNames) ];

        sink=emitSheet(sink,T,'dvsData',wantSheet);

        T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
            array2table(results.dvsDiameter(:,goodMask), ...
            'VariableNames',roiNames) ];

        sink=emitSheet(sink,T,'dvsDiameter',wantSheet);
    end

    % ---- pulsatility summary sheets (per-segment results.pulsatility scalars) ---
    % Surface the tree's per-segment markers + harmonic coefficients as focused
    % sheets, row-aligned to sMetrics / dvsMetrics.  The per-pixel maps
    % results.pulsatility.ppx are image data and are intentionally NOT exported.
    if isfield(results,'pulsatility')
        if isfield(results.pulsatility,'sData') && isfield(results.pulsatility.sData,'scalars')
            Tp = pulsScalarTable(results.pulsatility.sData.scalars,'');
            if ~isempty(Tp)
                sink=emitSheet(sink,[idTable(results.sMetrics), Tp],'pulsatility',wantSheet);
            end
        end
        if isfield(results,'dvsMetrics') && isfield(results.pulsatility,'dvsData') ...
                && isfield(results.pulsatility.dvsData,'scalars')
            Tp = [idTable(results.dvsMetrics), ...
                  pulsScalarTable(results.pulsatility.dvsData.scalars,'ps')];
            if isfield(results.pulsatility,'dvsDiameter') ...
                    && isfield(results.pulsatility.dvsDiameter,'scalars')
                Tp = [Tp, pulsScalarTable(results.pulsatility.dvsDiameter.scalars,'pd')];
            end
            sink=emitSheet(sink,Tp,'dvsPulsatility',wantSheet);
        end
    end
     end
end

finishSink(sink);                                   % merged mode: write the one workbook
end

function Ti = idTable(M)
% Identifier columns for a pulsatility sheet: idx plus type/label when present.
% Row order is the segment order shared by results.sMetrics / dvsMetrics and the
% results.pulsatility.<sig> tree (DATA-MODEL segment-order invariant).
Ti = table(double(M.idx(:)),'VariableNames',{'idx'});
if ismember('type', M.Properties.VariableNames),  Ti.type  = string(M.type(:));  end
if ismember('label',M.Properties.VariableNames),  Ti.label = string(M.label(:)); end
end

function Tp = pulsScalarTable(S,harmPrefix)
% Flatten one results.pulsatility.<sig>.scalars struct (per-segment) into a table:
% each [nSeg x 1] scalar becomes a column under its own name; the harmonic
% coefficients hAmp/hPhase [nSeg x nHarm] expand to <harmPrefix>hAmp1..N /
% <harmPrefix>hPhase1..N (harmPrefix disambiguates the bare flow vs diameter
% coefficients when a sheet carries both).  Non-[nSeg x 1] extras are skipped;
% returns an empty table when nothing usable is present.
Tp = table();
if ~isstruct(S), return; end
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

idxAll   = results.sMetrics.idx(:)';
goodMask = ~isnan(idxAll) & idxAll~=0 & ~isnan(results.sMetrics.BFI(:)');

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
O = struct('sheets',{{}},'ext','.xlsx','groupBy','','weightByArea',true, ...
    'columns',{{}},'merge',false,'labels',[],'outFile','');
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
%defaultSheets  What .groupBy asks for when opts.sheets is not given: the whole
%   legacy workbook for per-segment rows, and the matching aggregated PAIR for an
%   averaging axis.  The two vessel-type sheets are therefore OPT-IN - a bare
%   exportToExcel(fNames) writes the historical nine and nothing else.
switch groupBy
    case 'label', names = {'sMetricsROI','sDataROI'};
    case 'type',  names = {'sMetricsType','sDataType'};
    otherwise,    names = legacySheets();
end
end

function names = legacySheets()
%legacySheets  The nine sheets a bare exportToExcel(fNames) has always written,
%   in write order.  THE definition of "everything" for the default path.
names = {'sMetrics','sData','sMetricsROI','sDataROI','dvsMetrics','dvsData', ...
         'dvsDiameter','pulsatility','dvsPulsatility'};
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
P = table(string(regexprep([nm ex],'_d\.mat$','')),'VariableNames',{'file'});
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
