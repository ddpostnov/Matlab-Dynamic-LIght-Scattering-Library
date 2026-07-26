% Example_Claude_vascularTreeStats - Summarise & visualise the vascular tree.
%
% WHAT THIS SCRIPT DOES
%   Loads one results file (*_r.mat) that has been through setVascularTree and
%   summarises the derived parent->daughter hierarchy: how many trees
%   (territories), roots (arterial inlets) and outlets (venous outlets) there
%   are, the branch-generation depth and the Horton-Strahler order
%   distribution. It then shows the hierarchy map and checks the hemodynamic
%   gradients the hierarchy is built on - pulsatility psPI and pulse-arrival
%   time versus branch generation.
%
% WHEN TO USE IT
%   You have run setVascularTree on a pulsatility ("_c") recording, so
%   results.hierarchy exists and results.sMetrics has the columns generation,
%   strahlerOrder, treeID, parentIdx, daughterIdx, hierarchyConfidence.
%
% HOW TO USE IT
%   1. Set fName to your *_BFI_r.mat results file.
%   2. Run the script.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar script)
%   "Write a simple, well-commented MATLAB script for non-expert users. It loads
%    a results .mat file with a struct 'results'. results.hierarchy holds a
%    vascular parent-daughter graph (fields nodeIds, edges [parent child w],
%    roots, outlets, nodeType, generation, strahler, treeID). results.sMetrics
%    is a table with one row per segment and columns category (5=vessel lumen,
%    3=wall, 1=parenchyma), type ('Artery'/'Vein'), BFI, 'psPI', psTimeMin,
%    generation, strahlerOrder, treeID, hierarchyConfidence, isRoot, isOutlet.
%    Print a summary (number of territories, roots, outlets, node counts by
%    type, generation range, Strahler histogram, median confidence). Then show
%    results.mapTree as an image and draw box plots of psPI and psTimeMin
%    versus generation for the vessel lumen rows. If results.hierarchy is
%    missing, stop with a clear message telling the user to run setVascularTree."
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.

%% 0. Point MATLAB to the library -----------------------------------------
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));

%% 1. Load one results file -----------------------------------------------
% EDIT THIS: full path to a *_BFI_r.mat results file processed by setVascularTree.
if ~exist('fName','var') || isempty(fName)
    fName = 'C:\path\to\your\data\recording_c_BFI_r.mat';
end
load(fName,'results');

if ~isfield(results,'hierarchy')
    error('No results.hierarchy found. Run setVascularTree on this recording first.');
end
H = results.hierarchy;
T = results.sMetrics;

%% 2. Printed summary ------------------------------------------------------
nt = H.nodeType;
fprintf('\nVascular tree summary for:\n  %s\n', fName);
fprintf('  nodes: %d  (Artery %d, Vein %d, Parenchyma %d, Uncertain %d)\n', ...
    numel(H.nodeIds), sum(nt=="Artery"), sum(nt=="Vein"), sum(nt=="Parench"), sum(nt=="Uncertain"));
fprintf('  edges: %d   territories (treeID): %d\n', size(H.edges,1), numel(unique(H.treeID)));
fprintf('  roots (arterial inlets): %d   outlets (venous): %d\n', numel(H.roots), numel(H.outlets));
fprintf('  generation depth: %d..%d (median %d)\n', min(H.generation), max(H.generation), round(median(H.generation)));
fprintf('  Strahler order histogram:');
for k = 1:max(H.strahler), fprintf('  [%d]=%d', k, sum(H.strahler==k)); end
fprintf('\n  median relation confidence: %.2f\n\n', median(H.confidence,'omitnan'));

%% 3. Hierarchy map --------------------------------------------------------
figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
if isfield(results,'mapTree')
    imagesc(results.mapTree,'AlphaData',~isnan(results.mapTree)); axis image off
    title('Horton-Strahler order (high = trunk)'); colorbar
    colormap(gca,parula)
end

%% 4. Do the hemodynamic gradients hold along the tree? -------------------
% Use vessel lumen rows (category 5) with a valid generation.
isLum = T.category==5 & isfinite(T.generation);
g     = T.generation(isLum);
haveP = ismember('psPI',T.Properties.VariableNames);
nexttile;
if haveP
    boxchart(categorical(g), T.('psPI')(isLum));
    xlabel('generation'); ylabel('psPI');
    title('Pulsatility vs branch generation'); grid on
end

% Pulse-arrival (foot) vs generation, split arterial/venous, in a 2nd figure.
if ismember('psTimeMin',T.Properties.VariableNames)
    figure('Color','w','Position',[120 120 700 420]);
    tp = string(T.type);
    isA = isLum & tp=="Artery"; isV = isLum & tp=="Vein";
    hold on
    if any(isA), plot(T.generation(isA), 1e3*T.psTimeMin(isA), 'o','Color',[.85 .2 .2],'DisplayName','Artery'); end
    if any(isV), plot(T.generation(isV), 1e3*T.psTimeMin(isV), 's','Color',[.2 .4 .85],'DisplayName','Vein'); end
    hold off; grid on; legend('Location','best')
    xlabel('generation'); ylabel('pulse foot arrival (ms)');
    title('Pulse-arrival delay along the tree');
end
