% Example_Claude_plotROIs - Draw ROIs on a processed recording and plot / report their signals.
%
% WHAT THIS SCRIPT DOES
%   Loads one processed dataset (a *_BFI_d.mat "data" file and, if present, its
%   *_r.mat "results" file), shows the time-averaged blood-flow-index (BFI)
%   image with robust colour limits, lets you draw one or more regions of
%   interest (ROIs) with the mouse, and then:
%     - plots each ROI's mean BFI signal over time (raw and normalised), and
%     - prints a small table with each ROI's mean, std, min and max BFI.
%
% WHEN TO USE IT
%   You have processed a recording up to a *_BFI_d.mat file and want a quick,
%   manual look at how flow behaves in a few places you pick by hand (e.g. a
%   vessel and the tissue around it).
%
% HOW TO USE IT
%   1. Set fName below to your *_BFI_d.mat file (the *_r.mat results file next
%      to it is loaded automatically if it exists).
%   2. Run the script and draw each ROI with the mouse (click the corners,
%      then double-click to close it). A dialog then lets you add another
%      ROI, remove the last one, or finish.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar script)
%   "Write a simple, well-commented MATLAB script for users who are not
%    confident with MATLAB. It loads a processed laser-speckle dataset from a
%    .mat file that contains a struct 'source' with a 3-D array source.data
%    (Y x X x Time, blood-flow index) and a time vector source.time (seconds),
%    plus an optional struct 'results' that may contain results.cMask (a tissue
%    mask). Show the time-averaged image using robust colour limits (ignore
%    NaN/Inf and the top/bottom few percent). Let the user draw several ROIs
%    with the mouse using MATLAB's built-in drawpolygon/createMask (keep the
%    whole ROI-drawing loop inside the script - do not call a custom helper),
%    compute each ROI's mean signal over time, and plot the raw
%    and the mean-normalised signals in separate panels. Finally print a table
%    of each ROI's mean, std, min and max value."
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.

%% 0. Point MATLAB to the library -----------------------------------------
% Auto-detected from this file's location. If you COPY this script somewhere
% else, set libraryFolder to the library folder by hand instead.
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));

%% 1. Load one processed recording ----------------------------------------
% EDIT THIS: full path to a processed *_BFI_d.mat file.
fName = 'C:\path\to\your\data\recording_t_BFI_d.mat';

load(fName,'source');                                    % source.data, source.time
resultsFile = strrep(fName,'_d.mat','_r.mat');
if isfile(resultsFile)
    load(resultsFile,'results');
else
    results = struct();                                  % results are optional here
end

%% 2. Time-averaged BFI image with robust colour limits -------------------
img = mean(source.data,3,'omitnan');
if isfield(results,'cMask')
    mask = results.cMask>0;                              % tissue only
else
    mask = isfinite(img) & img>0;
end
clim0 = prctile(img(mask(:)),[5 99]);                    % robust display range

figH = figure('Color','w');
imagesc(img); axis image; clim(clim0); colorbar
title('Time-averaged BFI - draw your ROIs')

%% 3. Draw ROIs interactively (all control lives in this script) ----------
% Click the corners of a region, then double-click (or press Enter) to close
% it. After each one a dialog lets you add another ROI, remove the last one,
% or finish. This is a plain drawpolygon loop - no getAndShowROIS helper.
figure(figH);
colors = lines(999);
hPoly  = {};                                             % polygons drawn so far
hold on
keepGoing = true;
while keepGoing
    title(sprintf('Draw ROI %d: click corners, double-click to close', ...
                  numel(hPoly)+1))
    p = drawpolygon(gca,'Color',colors(numel(hPoly)+1,:));
    if size(p.Position,1) >= 3                           % a usable polygon
        hPoly{end+1} = p;                                           %#ok<SAGROW>
    else
        delete(p);                                       % discard empty draw
    end
    switch questdlg('Add another ROI, remove the last one, or finish?', ...
                    'ROI selection','Add another','Remove last','Finish','Add another')
        case 'Remove last'
            if ~isempty(hPoly), delete(hPoly{end}); hPoly(end) = []; end
        case 'Finish'
            keepGoing = false;
    end
end
hold off

% Turn the drawn polygons into a [Y x X x nROI] logical stack.
nROI = numel(hPoly);
if nROI==0
    disp('No ROIs were drawn - nothing to plot.'); return
end
ROIs = false(size(img,1),size(img,2),nROI);
for k = 1:nROI
    ROIs(:,:,k) = createMask(hPoly{k},img);              % polygon -> logical mask
end

%% 4. Extract each ROI's mean signal over time ----------------------------
t    = source.time(:);
flat = reshape(source.data,[],numel(t));                 % (Y*X) x Time
sig  = zeros(numel(t),nROI);
for k = 1:nROI
    pix      = squeeze(ROIs(:,:,k))==1;
    sig(:,k) = mean(flat(pix(:),:),1,'omitnan');
end

%% 5. Plot ROI signals -----------------------------------------------------
figure('Color','w','Position',[100 100 1100 500]);
tiledlayout(2,4,'TileSpacing','tight','Padding','tight');

nexttile([2 2]);                                         % image + ROI overlay
imagesc(img); axis image; clim(clim0);
hold on
for k = 1:nROI
    visboundaries(ROIs(:,:,k),'Color',colors(k,:));      % outline each ROI
end
hold off
title('ROIs')

nexttile([1 2]);                                         % raw signals
plot(t,sig,'LineWidth',1.1); axis tight; grid on
xlabel('Time (s)'); ylabel('BFI'); title('Raw ROI signals')
legend(compose('ROI %d',1:nROI),'Location','best')

nexttile([1 2]);                                         % mean-normalised signals
plot(t,sig./mean(sig,1,'omitnan'),'LineWidth',1.1); axis tight; grid on
xlabel('Time (s)'); ylabel('BFI / mean'); title('Mean-normalised ROI signals')

%% 6. Report ROI values ----------------------------------------------------
roiTable = table((1:nROI)', mean(sig,1,'omitnan')', std(sig,0,1,'omitnan')', ...
                 min(sig,[],1)', max(sig,[],1)', ...
                 'VariableNames',{'ROI','meanBFI','stdBFI','minBFI','maxBFI'});
disp('ROI summary:'); disp(roiTable)
