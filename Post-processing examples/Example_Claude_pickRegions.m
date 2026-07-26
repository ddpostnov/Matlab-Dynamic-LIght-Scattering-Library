% Example_Claude_pickRegions - Click labelled regions in a results file and read their metrics.
%
% WHAT THIS SCRIPT DOES
%   Loads ONLY a results file (*_r.mat) - it does NOT open the large 3-D data
%   file - shows the segment map, and lets you click any labelled vessel or
%   tissue segment. For every click it prints that segment's metrics from
%   results.sMetrics (index, type/label if available, diameter, area, mean BFI,
%   and pulsatility / vasomotion values when those exist).
%
% WHEN TO USE IT
%   You have a processed & segmented recording and just want to look up the
%   numbers for a particular vessel or region, quickly, without waiting for the
%   full dataset to load.
%
% HOW TO USE IT
%   1. Set fName below to your *_r.mat results file.
%   2. Run the script, then click segments on the image. Press Enter (with no
%      click) to finish.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar script)
%   "Write a simple, well-commented MATLAB script for users who are not good
%    with MATLAB. It loads only a results .mat file containing a struct
%    'results' with: results.sMap (a 2-D label map where each vessel/tissue
%    segment has an integer id and the background is NaN), results.sMetrics (a
%    table with one row per segment, always including 'idx','diameter','area'
%    and often 'BFI','psPI','type','label'), and optionally results.imgBFI
%    (a 2-D image). Display imgBFI (or the label map if it is missing) with
%    robust colour limits, overlay the segment outlines, then let the user click
%    segments one by one. For each click, work out which segment was clicked and
%    print its metrics in a readable way, skipping columns that do not exist.
%    Stop when the user presses Enter."
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.

%% 0. Point MATLAB to the library -----------------------------------------
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));

%% 1. Load ONLY the results file ------------------------------------------
% EDIT THIS: full path to a *_r.mat results file.
fName = 'C:\path\to\your\data\recording_t_BFI_r.mat';
load(fName,'results');

sMap = results.sMap;                          % label map (NaN background)
L    = sMap; L(isnan(L)) = 0; L = round(L);   % integer labels, 0 = background

%% 2. Background image -----------------------------------------------------
if isfield(results,'imgBFI')
    img = results.imgBFI;                     % nicer backdrop if available
else
    img = L;                                  % otherwise show the label map
end
m     = isfinite(img) & img~=0;
clim0 = prctile(img(m(:)),[5 99]);

figure('Color','w');
imagesc(img); axis image; colormap gray
if clim0(2)>clim0(1), clim(clim0); end
hold on; visboundaries(L>0,'Color',[1 0.4 0],'LineWidth',0.5);
title('Click a segment to read its metrics (press Enter to finish)')

%% 3. Click loop -----------------------------------------------------------
vn = results.sMetrics.Properties.VariableNames;
fprintf('\nClick segments on the figure. Press Enter (no click) to stop.\n');
while true
    [x,y] = ginput(1);
    if isempty(x), break; end                 % Enter pressed -> done
    row = round(y); col = round(x);
    if row<1 || col<1 || row>size(L,1) || col>size(L,2), continue; end

    lbl = L(row,col);
    if lbl==0
        fprintf('  (background - no segment here)\n'); continue;
    end
    r = find(results.sMetrics.idx==lbl,1);
    if isempty(r)
        fprintf('  segment %d is not in the metrics table\n',lbl); continue;
    end
    plot(col,row,'w+','MarkerSize',10,'LineWidth',1.5);

    fprintf('\n--- segment %d ---\n',lbl);
    printField(results.sMetrics,r,vn,'type',    '   type      : %s\n');
    printField(results.sMetrics,r,vn,'label',   '   label     : %s\n');
    printField(results.sMetrics,r,vn,'diameter','   diameter  : %.2f px\n');
    printField(results.sMetrics,r,vn,'area',    '   area      : %.0f px\n');
    printField(results.sMetrics,r,vn,'BFI',     '   mean BFI  : %.3g\n');
    printField(results.sMetrics,r,vn,'psPI',    '   psPI      : %.3f\n');
    printField(results.sMetrics,r,vn,'ampMeanVB','   vaso ampMeanVB: %.3g\n');
end
disp('Done.');

%% ---- local helper -------------------------------------------------------
function printField(T,row,vn,name,fmt)
% Print one metric, but only if that column exists in the table.
    if ~ismember(name,vn), return; end
    v = T.(name)(row);
    if isnumeric(v)
        fprintf(fmt,v);
    else
        fprintf(fmt,string(v));               % text columns (type / label)
    end
end
