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
%
%   INPUT
%     fNames   cell array of *.mat file paths.  Each must contain SOURCE,
%              RESULTS, SETTINGS structures created by the LSCI pipeline.
%
%   OUTPUT
%     None – one Excel workbook per input file is written to disk.
%
%   EXAMPLE
%     files = dir(fullfile(dataRoot,'*BFI_d.mat'));
%     exportToExcel(fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     MATLAB R2019b+ (for writetable with 'Sheet' option) and data schema
%     used throughout the Dynamic Light Scattering Imaging toolbox.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: This header was generated with ChatGPT and may contain minor
%   inconsistencies—please verify before distribution.
%   ----------------------------------------------------------------------

function exportToExcel(fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".mat".');
end

for fidx=1:1:numel(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    fName=fNames{fidx};
    clearvars results source settings
    load(fNames{fidx})
    load(strrep(fNames{fidx},'_d.mat','_s.mat'));
    load(strrep(fNames{fidx},'_d.mat','_r.mat'));

    fName=strrep(fName,'_d.mat','.xlsx');

    catNames = ["background"  "parenchyma"  "unsegmented" ...
        "outerWall"   "innerWall"   "lumen"];          % 1×6

    valid   = ~isnan(results.sMetrics.category);                % ignore NaNs
    results.sMetrics.category = categorical( ...
        results.sMetrics.category, 0:5, catNames);          % numeric → cat


    T=results.sMetrics;
    T(end,:) = []; % crappy bug workaround
    T(end+1,:) = results.sMetrics(end,:);% crappy bug workaround
    writetable(T,fName,'Sheet','sMetrics','WriteMode','replacefile');


    idxAll    = results.sMetrics.idx(:)';
    goodMask  = ~isnan(idxAll);
    roiIdx    = idxAll(goodMask);
    roiData   = results.sData(:,goodMask);
    roiNames  = matlab.lang.makeValidName(compose("ROI %04d",roiIdx) );
    T = [ table(results.time(:),'VariableNames', {'Time'}), array2table(roiData,'VariableNames',roiNames) ];
    writetable(T,fName,'Sheet','sData');

    if ismember('label', results.sMetrics)

        keep =  strlength(results.sMetrics.label) > 0  ...
            & ~ismember(results.sMetrics.category, ["outerWall","innerWall"]);

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
        writetable(T,fName,'Sheet','sMetricsROI');



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
        writetable(T,fName,'Sheet','sDataROI');
    end


    if isfield(results,"dvsMetrics")
        T=results.dvsMetrics;
        T(end,:) = []; % crappy bug workaround
        T(end+1,:) = results.dvsMetrics(end,:);% crappy bug workaround
        writetable(T,fName,'Sheet','dvsMetrics');

        idxAll   = results.dvsMetrics.idx(:)';       % row vector
        goodMask = ~isnan(idxAll);                   % ignore NaN indices
        roiIdx   = idxAll(goodMask);
        roiNames = matlab.lang.makeValidName(compose("ROI %04d", roiIdx));

        T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
            array2table(results.dvsData(:,goodMask), ...
            'VariableNames',roiNames) ];
        writetable(T,fName,'Sheet','dvsData');

        T = [ table(results.time(:), 'VariableNames',{'Time'}), ...
            array2table(results.dvsDiameter(:,goodMask), ...
            'VariableNames',roiNames) ];
        writetable(T,fName,'Sheet','dvsDiameter');
    end

end
end