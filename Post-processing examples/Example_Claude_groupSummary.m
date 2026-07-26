% Example_Claude_groupSummary - Compare key vessel parameters across groups of recordings.
%
% WHAT THIS SCRIPT DOES
%   Collects many results files (*_BFI_r.mat) from a folder, splits them into
%   groups with getFileNamesList, and plots the main vessel parameters for each
%   group in separate panels. It ALWAYS shows BFI and diameter, and ADDS a
%   pulsatility panel and a vasomotion panel only if those parameters exist in
%   the files. Each panel shows the per-group average with error bars and every
%   recording as a dot.
%
% WHEN TO USE IT
%   You processed several recordings (e.g. different animals or conditions) and
%   want one figure that compares groups (e.g. control vs treated).
%
% HOW TO USE IT
%   1. Set rootFolder to the folder that contains your *_BFI_r.mat files.
%   2. Set groupPattern to a short regular expression that picks the GROUP name
%      out of each file name (e.g. 'WT|KO' or 'Ctrl|Stroke'). Every unique match
%      becomes one group. Leave it '' to treat all files as a single group.
%   3. Run the script.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar script)
%   "Write a simple, well-commented MATLAB script for non-expert users. Using a
%    helper called getFileNamesList(rootFolder, glob, groupPattern) that returns
%    a cell array with one row of file paths per group, load every results .mat
%    file (each has a struct 'results' with a table results.sMetrics that has one
%    row per vessel/tissue segment; column 'category' is 5 for vessels; columns
%    'BFI' and 'diameter' always exist; columns 'psPI' and 'ampMeanVB' may exist).
%    For each recording, average each parameter over its vessel segments. Then
%    draw one figure with a panel per parameter family (BFI and diameter always;
%    pulsatility and vasomotion only if present), showing the per-group mean with
%    error bars and the individual recordings as dots. Keep the code short and
%    readable and skip missing parameters gracefully."
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.

%% 0. Point MATLAB to the library -----------------------------------------
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));

%% 1. Settings - EDIT THESE -----------------------------------------------
rootFolder   = 'C:\path\to\your\data';   % folder searched recursively
fileGlob     = '*_BFI_r.mat';            % which results files to use
groupPattern = 'WT|KO';                  % regexp picking the group label ('' = one group)

%% 2. Collect files and split into groups ---------------------------------
if isempty(groupPattern)
    flat        = getFileNamesList(rootFolder,fileGlob);            % N-by-1
    groupFiles  = {flat(:)};
    groupLabels = {'all'};
else
    fNames      = getFileNamesList(rootFolder,fileGlob,groupPattern); % nG-by-cols
    nG          = size(fNames,1);
    groupFiles  = cell(nG,1);
    groupLabels = cell(nG,1);
    for g = 1:nG
        r              = fNames(g,:);
        r              = r(~cellfun(@isempty,r));                   % drop empty cells
        groupFiles{g}  = r;
        [~,nm]         = fileparts(r{1});
        groupLabels{g} = regexp(nm,groupPattern,'match','once');   % label for the x-axis
    end
end
nGroups = numel(groupFiles);
if sum(cellfun(@numel,groupFiles))==0
    error('No files matched "%s" under "%s". Check rootFolder and fileGlob.',fileGlob,rootFolder);
end

%% 3. Read one summary value per recording --------------------------------
% Long lists (one entry per recording) plus a matching group index G.
G = []; BFI = []; DIAM = []; PI = []; VFR = [];
for g = 1:nGroups
    for f = 1:numel(groupFiles{g})
        S   = load(groupFiles{g}{f},'results');
        T   = S.results.sMetrics;
        ves = T.category==5;                      % vessel segments only
        G(end+1,1)    = g;                         %#ok<*SAGROW>
        BFI(end+1,1)  = colMean(T,'BFI',ves);
        DIAM(end+1,1) = colMean(T,'diameter',ves);
        PI(end+1,1)   = colMean(T,'psPI',ves);
        VFR(end+1,1)  = colMean(T,'ampMeanVB',ves);
    end
end

%% 4. Decide which panels to show -----------------------------------------
% {title, values, y-label}.  BFI and diameter are always shown.
panels = {'BFI',BFI,'Mean BFI (a.u.)'; 'Diameter',DIAM,'Mean diameter (px)'};
if any(isfinite(PI)),  panels(end+1,:) = {'Pulsatility psPI',PI,'Mean psPI'};      end
if any(isfinite(VFR)), panels(end+1,:) = {'Vasomotion ampMeanVB',VFR,'Mean ampMeanVB (a.u.)'};   end

%% 5. Plot -----------------------------------------------------------------
nP = size(panels,1);
figure('Color','w','Position',[80 80 340*nP 430]);
tiledlayout(1,nP,'TileSpacing','compact','Padding','compact');
for p = 1:nP
    nexttile;
    plotGroupMetric(G,panels{p,2},nGroups,groupLabels,panels{p,1},panels{p,3});
end

%% ---- local helpers ------------------------------------------------------
function v = colMean(T,name,rows)
% Mean of a table column over selected rows; NaN if the column is absent.
    if ismember(name,T.Properties.VariableNames)
        v = mean(T.(name)(rows),'omitnan');
    else
        v = NaN;
    end
end

function plotGroupMetric(G,vals,nGroups,labels,ttl,ylab)
% Bar of per-group mean + SEM, with the individual recordings as dots.
    mu = nan(nGroups,1); se = nan(nGroups,1);
    for g = 1:nGroups
        x = vals(G==g); x = x(isfinite(x));
        if ~isempty(x)
            mu(g) = mean(x);
            se(g) = std(x)/sqrt(numel(x));
        end
    end
    bar(1:nGroups,mu,'FaceColor',[0.80 0.85 0.95]); hold on
    errorbar(1:nGroups,mu,se,'k','LineStyle','none','LineWidth',1);
    for g = 1:nGroups
        x = vals(G==g); x = x(isfinite(x));
        scatter(g+0.12*(rand(numel(x),1)-0.5),x,18,'k','filled','MarkerFaceAlpha',0.6);
    end
    set(gca,'XTick',1:nGroups,'XTickLabel',labels);
    ylabel(ylab); title(ttl); grid on; axis padded
end
