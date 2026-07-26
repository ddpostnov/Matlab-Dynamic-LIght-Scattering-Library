% Example_Claude_arteryVeinStats - Compare arteries vs veins in one recording.
%
% WHAT THIS SCRIPT DOES
%   Loads one results file (*_r.mat) and compares arteries against veins using
%   the vessel type assigned by setVesselTypes. It prints a small summary table
%   (per type: number of vessels and median BFI / diameter / pulsatility) and
%   draws box plots of each parameter for arteries and veins side by side.
%
% WHEN TO USE IT
%   You have run setVesselTypes on a recording (so results.sMetrics has a 'type'
%   column) and want to see how arterial and venous vessels differ.
%
% HOW TO USE IT
%   1. Set fName to your *_r.mat results file.
%   2. Run the script.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar script)
%   "Write a simple, well-commented MATLAB script for non-expert users. It loads
%    a results .mat file with a struct 'results' whose table results.sMetrics has
%    one row per vessel/tissue segment and includes a text column 'type' with
%    values like 'Artery','Vein','Uncertain', a numeric 'BFI' column, a
%    'diameter' column and sometimes a 'psPI' column. Keep only arteries and
%    veins, print a summary table (count and median BFI / diameter / PI per type)
%    and draw grouped box plots comparing arteries and veins for each available
%    parameter. If the 'type' column is missing, stop with a clear message
%    telling the user to run the vessel-typing step first."
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.

%% 0. Point MATLAB to the library -----------------------------------------
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));

%% 1. Load one results file -----------------------------------------------
% EDIT THIS: full path to a *_r.mat results file that already has vessel types.
fName = 'C:\path\to\your\data\recording_t_BFI_r.mat';
load(fName,'results');
T = results.sMetrics;

if ~ismember('type',T.Properties.VariableNames)
    error('No "type" column found. Run setVesselTypes on this recording first.');
end

%% 2. Keep arteries and veins only ----------------------------------------
tp   = string(T.type);
keep = tp=="Artery" | tp=="Vein";
T    = T(keep,:);
tp   = tp(keep);
if isempty(T)
    error('No segments labelled Artery or Vein were found.');
end

%% 3. Which parameters can we show? ---------------------------------------
% Column name  ->  axis label.  BFI is always present; the rest are optional.
params = {'BFI','BFI (a.u.)'};
if ismember('diameter',T.Properties.VariableNames), params(end+1,:) = {'diameter','Diameter (px)'}; end %#ok<*AGROW>
if ismember('psPI', T.Properties.VariableNames), params(end+1,:) = {'psPI','psPI'};        end

%% 4. Printed summary (medians) -------------------------------------------
types    = ["Artery";"Vein"];
nSeg     = arrayfun(@(x) sum(tp==x), types);
summaryT = table(types,nSeg,'VariableNames',{'type','nVessels'});
for p = 1:size(params,1)
    med = arrayfun(@(x) median(T.(params{p,1})(tp==x),'omitnan'), types);
    summaryT.(['median_' matlab.lang.makeValidName(params{p,1})]) = med;
end
disp('Artery vs Vein summary (medians):'); disp(summaryT)

%% 5. Box plots ------------------------------------------------------------
figure('Color','w','Position',[100 100 340*size(params,1) 400]);
tiledlayout(1,size(params,1),'TileSpacing','compact','Padding','compact');
for p = 1:size(params,1)
    nexttile;
    boxchart(categorical(tp),T.(params{p,1}));
    ylabel(params{p,2}); title(params{p,1}); grid on
end
