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
%     opts     (optional) struct selecting WHAT and in WHICH format to write.
%              Omitting it - or passing an empty / no-field struct - reproduces the
%              full legacy behaviour above exactly (the five launcher callers do).
%              Fields:
%                .sheets  cellstr / string of sheet names to write, a subset of
%                         {sMetrics, sData, sMetricsROI, sDataROI, dvsMetrics,
%                          dvsData, dvsDiameter, pulsatility, dvsPulsatility}.
%                         Absent / empty = all sheets.  A selected sheet whose data
%                         is absent from a given file is simply skipped, exactly as
%                         in the full path.
%                .format  output extension, '.xlsx' (default) or '.xls'.
%              The first sheet actually written to a file uses 'replacefile', so a
%              re-export with a narrower selection never leaves stale sheets behind.
%
%   OUTPUT
%     None – one Excel workbook per input file is written to disk.
%
%   EXAMPLE
%     files = dir(fullfile(dataRoot,'*BFI_d.mat'));
%     fNames = fullfile({files.folder}',{files.name}');
%     exportToExcel(fNames);                              % full legacy workbook
%     exportToExcel(fNames, struct('sheets',{{'sMetrics','pulsatility'}}));
%                                                         % just those two sheets
%
%   DEPENDS ON
%     MATLAB R2019b+ (for writetable with 'Sheet' option) and data schema
%     used throughout the Dynamic Light Scattering Imaging toolbox.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

function exportToExcel(fNames, opts)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".mat".');
end

if nargin<2, opts = struct(); end
[sheetSel, outExt] = parseExportOpts(opts);         % {} = all sheets; '.xlsx' default
wantSheet = @(nm) isempty(sheetSel) || any(strcmp(nm, sheetSel));

for fidx=1:1:numel(fNames)
     if ~isempty(fNames{fidx})
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    fName=fNames{fidx};
    load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

    fName=strrep(fName,'_d.mat',outExt);
    wroteAny=false;                                     % first written sheet -> replacefile

    catNames = ["background"  "parenchyma"  "unsegmented" ...
        "outerWall"   "innerWall"   "lumen"];          % 1×6

    results.sMetrics.category(isnan( results.sMetrics.category))=0;
    results.sMetrics.category = categorical( ...
        results.sMetrics.category, 0:5, catNames);          % numeric → cat


    T=results.sMetrics;
    T(end,:) = []; % crappy bug workaround
    T(end+1,:) = results.sMetrics(end,:);% crappy bug workaround
    
    wroteAny=emitSheet(fName,T,'sMetrics',wantSheet,wroteAny);


    idxAll    = results.sMetrics.idx(:)';
    goodMask  = ~isnan(idxAll) & idxAll~=0 & ~isnan(results.sMetrics.BFI(:)');
    roiIdx    = idxAll(goodMask);
    roiData   = results.sData(:,goodMask);
    roiNames  = matlab.lang.makeValidName(compose("ROI %04d",roiIdx) );
    T = [ table(results.time(:),'VariableNames', {'Time'}), array2table(roiData,'VariableNames',roiNames) ];
    
    wroteAny=emitSheet(fName,T,'sData',wantSheet,wroteAny);

    if ismember('label', results.sMetrics.Properties.VariableNames)
        

        keep = goodMask' & strlength(results.sMetrics.label) > 0  ...
            & ~ismember(results.sMetrics.category, ["outerWall","innerWall","unsegmented"]);

        M        = results.sMetrics(keep ,:);
        areaAll  = M.area;                         % weights   (unchanged)
        lenAll   = M.length;                       % new       (to be summed)
        if ismember('RegOverlap', M.Properties.VariableNames)
            overlap   = M.RegOverlap;                       % new       (to be summed)
        end
        [G,labelList] = findgroups(M.label);
        nG = numel(labelList);

        % numeric variables to weight-average   (skip idx, area, length)
        isNum   = varfun(@isnumeric,M,'OutputFormat','uniform');
        numVar  = M.Properties.VariableNames(isNum);
        numVar(ismember(numVar,{'idx','area','length','RegOverlap','closest vessel idx','nearestVesIdx','RegID'})) = [];

        % pre-allocate output: additional 'length' column
        T = table('Size',[nG numel(numVar)+4], ...        % +area +length +label +type
            'VariableTypes',[repmat("double",1,numel(numVar)+2) "string" "string"], ...
            'VariableNames',[numVar "area" "length" "label" "type"]);

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

            w              = areaAll(rows) ./ areaSum;    % normalised weights
            for v = 1:numel(numVar)                       % weighted means
                col        = numVar{v};
                T{g,col}   = sum( M.(col)(rows) .* w );
            end

            T.type(g)  = string( mode( categorical(M.type(rows)) ) );
            T.label(g) = labelList(g);
        end

        
        T = movevars(T,{'label','type'},'Before',1);      % final layout
        wroteAny=emitSheet(fName,T,'sMetricsROI',wantSheet,wroteAny);



        timeVec  = results.time(:);
        sDataUse = results.sData(:,keep);
        sigAgg   = zeros(numel(timeVec), nG);

        for g = 1:nG
            rows = (G==g);
            w    = areaAll(rows) ./ sum(areaAll(rows));         % same weights
            sigAgg(:,g) = sDataUse(:,rows) * w;                 % N×k × k×1 = N×1
        end

        sigNames = matlab.lang.makeValidName( compose("%s", labelList) );

        T = [ table(timeVec,'VariableNames',{'Time'}), ...
            array2table(sigAgg,'VariableNames',sigNames) ];
        
        wroteAny=emitSheet(fName,T,'sDataROI',wantSheet,wroteAny);
    end
    


    if isfield(results,"dvsMetrics")
        T=results.dvsMetrics;
        T(end,:) = []; % crappy bug workaround
        T(end+1,:) = results.dvsMetrics(end,:);% crappy bug workaround
        wroteAny=emitSheet(fName,T,'dvsMetrics',wantSheet,wroteAny);

        idxAll   = results.dvsMetrics.idx(:)';       % row vector
        goodMask = ~isnan(idxAll) & idxAll~=0 & ~isnan(results.dvsMetrics.BFI(:)');
        roiIdx   = idxAll(goodMask);
        roiNames = matlab.lang.makeValidName(compose("ROI %04d", roiIdx));

        T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
            array2table(results.dvsData(:,goodMask), ...
            'VariableNames',roiNames) ];
        
        wroteAny=emitSheet(fName,T,'dvsData',wantSheet,wroteAny);

        T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
            array2table(results.dvsDiameter(:,goodMask), ...
            'VariableNames',roiNames) ];
        
        wroteAny=emitSheet(fName,T,'dvsDiameter',wantSheet,wroteAny);
    end

    % ---- pulsatility summary sheets (per-segment results.pulsatility scalars) ---
    % Surface the tree's per-segment markers + harmonic coefficients as focused
    % sheets, row-aligned to sMetrics / dvsMetrics.  The per-pixel maps
    % results.pulsatility.ppx are image data and are intentionally NOT exported.
    if isfield(results,'pulsatility')
        if isfield(results.pulsatility,'sData') && isfield(results.pulsatility.sData,'scalars')
            Tp = pulsScalarTable(results.pulsatility.sData.scalars,'');
            if ~isempty(Tp)
                wroteAny=emitSheet(fName,[idTable(results.sMetrics), Tp],'pulsatility',wantSheet,wroteAny);
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
            emitSheet(fName,Tp,'dvsPulsatility',wantSheet,wroteAny);   % last sheet: return unneeded
        end
    end
     end
end
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

function [sheetSel, outExt] = parseExportOpts(opts)
%parseExportOpts  Normalise the optional selection struct into a sheet filter and
%   an output extension.  Absent / empty opts reproduce today's full behaviour.
%     opts.sheets : names to write (cellstr / string / char).  [] or absent = ALL.
%     opts.format : output workbook extension.  '.xlsx' (default) or '.xls' - both
%                   multi-sheet spreadsheet formats writetable's 'Sheet' supports.
sheetSel = {};
outExt   = '.xlsx';
if isempty(opts) || ~isstruct(opts), return; end
if isfield(opts,'sheets') && ~isempty(opts.sheets)
    sheetSel = cellstr(string(opts.sheets));
    sheetSel = sheetSel(:)';
end
if isfield(opts,'format') && ~isempty(opts.format)
    e = lower(char(string(opts.format)));
    if e(1) ~= '.', e = ['.' e]; end
    if ~any(strcmp(e,{'.xlsx','.xls'}))
        error('exportToExcel:format', ...
            'opts.format must be ''.xlsx'' or ''.xls'' (multi-sheet spreadsheet); got ''%s''.', e);
    end
    outExt = e;
end
end

function wroteAny = emitSheet(fName, T, sheet, wantSheet, wroteAny)
%emitSheet  Write one SELECTED sheet.  The first sheet actually written to a file
%   uses 'replacefile' (a fresh workbook); later sheets append.  A selection that
%   drops sMetrics therefore still starts from a clean file, and the default
%   (all-sheets) path is byte-identical to the historical "sMetrics with
%   'replacefile' first, the rest appended" behaviour.
if ~wantSheet(sheet), return; end
if ~wroteAny
    writetable(T, fName, 'Sheet', sheet, 'WriteMode','replacefile');
else
    writetable(T, fName, 'Sheet', sheet);
end
wroteAny = true;
end