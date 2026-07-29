% Launcher_speckleSimulation - Generate dynamic laser speckle datasets.
%
%   Five scenarios of dynamic speckle generation built on getDynamicSpeckles:
%     STEP 1  Exposure-integrated speckle (finite exposure time, LSCI-style).
%     STEP 2  Instantaneous frames, slowly decorrelating (Brownian, DLSI-style).
%     STEP 3  Instantaneous frames, fully decorrelated between frames.
%     STEP 4  Fully decorrelated, swept over pixels-per-speckle {0.5,1,2,4,8}.
%     STEP 5  Perfect speckle of 5 grain sizes written as .rls recordings.
%
%   Camera noise and camera effects are OFF by default - see the commented
%   'addNoise' lines in each step to switch them on. Run STEP 0 once per MATLAB
%   session, then run the step cells (%%) you need. STEPs 1-4 save a .mat
%   dataset to outputFolder and show a validation figure; STEP 5 writes .rls
%   files. The data also stays in the workspace. STEPs 1-4 are modest (~1 min
%   each on a GPU); STEP 5 is heavy at its default size (see its note).
%
% See also: getDynamicSpeckles, writeRLS, readRLS, getNormalizedG2, getK
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN RUN THE STEP YOU NEED
% Library path is auto-detected (this launcher lives in <library>/Simulation).
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
% Where the generated .mat datasets are written (edit if desired):
outputFolder = fullfile(libraryFolder,'Simulation','SimulatedData');
if ~exist(outputFolder,'dir'); mkdir(outputFolder); end


%% STEP 1 - EXPOSURE-INTEGRATED SPECKLE (finite exposure time, LSCI-style)
close all; clearvars -except libraryFolder outputFolder
% Slowly decorrelating Brownian dynamics integrated over a finite exposure:
% each output frame is the sum of exposureTime/dT instantaneous sub-frames, so
% the speckle blurs and the spatial contrast drops below 1 as T/tauC grows.
tauC             = 15;   % decorrelation time, us
dT               = 1;    % instantaneous sub-frame step, us
exposureTime     = 20;   % camera exposure per output frame, us
pixelsPerSpeckle = 2;    % speckle-to-pixel ratio (one speckle spans 2 pixels)
pixelsN          = 128;  % output pixels per side
nFrames          = 40;   % number of exposed frames
particlesN       = 1000; % number of scattering particles

[I, info] = getDynamicSpeckles('motionType','brownian','tauC',tauC,'dT',dT, ...
    'exposureTime',exposureTime,'pixelsPerSpeckle',pixelsPerSpeckle, ...
    'pixelsN',pixelsN,'nFrames',nFrames,'particlesN',particlesN,'verbose',true);
% To add camera noise + 8-bit digitisation, append (OFF by default):
%     'addNoise',true, 'cameraParams',struct('photonsPerPixel',500)

save(fullfile(outputFolder,'speckles_integrated.mat'),'I','info','-v7.3');

% Validation: per-frame spatial contrast K = std/mean (< 1 because of exposure).
K = squeeze(std(single(I),0,[1 2])./mean(single(I),[1 2]));
figure;
subplot(1,2,1); imagesc(I(:,:,1)); axis image; colormap gray; colorbar;
xlabel('pixels'); ylabel('pixels'); title('Exposure-integrated frame');
subplot(1,2,2); plot(K,'.-'); xlabel('frame'); ylabel('spatial contrast K'); ylim([0 1]);
title(sprintf('mean K=%.3f   (T/\\tau_c=%.1f)', mean(K), exposureTime/tauC));
set(gcf,'color','w');


%% STEP 2 - INSTANTANEOUS FRAMES, SLOWLY DECORRELATING (Brownian, DLSI-style)
close all; clearvars -except libraryFolder outputFolder
% No integration: every frame is a single instantaneous speckle pattern, and
% Brownian motion makes the pattern decorrelate slowly over ~tauC.
tauC             = 15;   % decorrelation time, us
dT               = 1;    % time between frames, us
pixelsPerSpeckle = 2;
pixelsN          = 128;
nFrames          = 800;  % many frames so g2 is well sampled
particlesN       = 1000;

[I, info] = getDynamicSpeckles('motionType','brownian','tauC',tauC,'dT',dT, ...
    'exposureTime',[],'pixelsPerSpeckle',pixelsPerSpeckle,'pixelsN',pixelsN, ...
    'nFrames',nFrames,'particlesN',particlesN,'verbose',true);
% Camera noise is OFF by default; append 'addNoise',true to switch it on.

save(fullfile(outputFolder,'speckles_raw_brownian.mat'),'I','info','-v7.3');

% Validation: g2(tau) decay and the recovered decorrelation time.
lagMax = round(4*tauC/dT);
g2 = squeeze(mean(getNormalizedG2(I,lagMax,size(I,3),[2 2]),[1 2]));
t  = (0:lagMax)'*dT;
ft = fittype('1+beta*exp(-2*(x/tau))+c', 'options', ...
    fitoptions('Method','NonlinearLeastSquares','Lower',[0 -0.5 0], ...
    'Upper',[1 0.5 100],'StartPoint',[0.5 0 tauC]));
fo = fit(t, double(g2), ft);
figure; plot(t, g2, 'o'); hold on; plot(t, 1+fo.beta*exp(-2*(t/fo.tau))+fo.c, '-','LineWidth',1.5);
xlabel('lag \tau, \mus'); ylabel('g_2(\tau)'); legend('simulated','fit'); set(gcf,'color','w');
title(sprintf('Brownian: input \\tau_c=%.1f, recovered \\tau_c=%.2f \\mus', tauC, fo.tau));


%% STEP 3 - INSTANTANEOUS FRAMES, FULLY DECORRELATED BETWEEN FRAMES
close all; clearvars -except libraryFolder outputFolder
% Every frame is an independent speckle realisation (particles fully
% re-randomised between frames), so g2 drops to 1 already at lag 1.
pixelsPerSpeckle = 2;
pixelsN          = 128;
nFrames          = 200;
particlesN       = 1000;

[I, info] = getDynamicSpeckles('motionType','decorrelated','exposureTime',[], ...
    'pixelsPerSpeckle',pixelsPerSpeckle,'pixelsN',pixelsN,'nFrames',nFrames, ...
    'particlesN',particlesN,'verbose',true);
% Camera noise is OFF by default; append 'addNoise',true to switch it on.

save(fullfile(outputFolder,'speckles_decorrelated.mat'),'I','info','-v7.3');

% Validation: g2(0) shows the full speckle contrast, g2(>=1) ~ 1.
g2 = squeeze(mean(getNormalizedG2(I,10,size(I,3),[2 2]),[1 2]));
figure; plot(0:10, g2, 'o-'); xlabel('lag, frames'); ylabel('g_2'); ylim([0.9 max(g2)*1.05]);
set(gcf,'color','w'); title(sprintf('Fully decorrelated: g_2(0)=%.3f, g_2(1)=%.3f', g2(1), g2(2)));


%% STEP 4 - FULLY DECORRELATED, SWEEP PIXELS-PER-SPECKLE {0.5, 1, 2, 4, 8}
close all; clearvars -except libraryFolder outputFolder
% A single fine, large-FOV field is generated ONCE at the finest sampling and
% block-averaged (intensity) to each ratio. Because every dataset comes from
% the SAME field, the illumination and speckle statistics are identical across
% ratios - only the pixel/speckle sampling changes. This is the clean way to
% study the speckle-to-pixel ratio without distorting the illumination.
sampPerSpeckle = 16;   % fine-grid samples per speckle (covers the finest bin)
pixelsN        = 160;  % pixels of the finest (pps=8) set; must be a multiple of 16
nFrames        = 150;  % independent frames (per-pixel contrast statistics)
particlesN     = 1500;

[base, baseInfo] = getDynamicSpeckles('motionType','decorrelated', ...
    'pixelsPerSpeckle',8,'sampPerSpeckle',sampPerSpeckle,'pixelsN',pixelsN, ...
    'nFrames',nFrames,'particlesN',particlesN,'verbose',true);
% Camera noise is OFF by default; append 'addNoise',true to switch it on.

ppsList = [8 4 2 1 0.5];         % pixels per speckle for each dataset
binList = [1 2 4 8 16];          % block-average factor applied to the fine set
K   = zeros(size(ppsList));
Kth = ppsList./sqrt(1+ppsList.^2);   % theory K = (s/p)/sqrt(1+(s/p)^2)
for i = 1:numel(ppsList)
    I = blockBinStack(base, binList(i));
    pixelsPerSpeckle = ppsList(i);
    info = baseInfo; info.pixelsPerSpeckle = ppsList(i); info.pixelSize = baseInfo.speckleSize/ppsList(i);
    K(i) = mean(std(single(I),0,3)./mean(single(I),3), 'all');
    save(fullfile(outputFolder,sprintf('speckles_pps%02d.mat', round(ppsList(i)*10))), ...
        'I','pixelsPerSpeckle','info','-v7.3');
end

figure; plot(ppsList, K, 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerSize',8); hold on;
plot(sort(ppsList), sort(ppsList)./sqrt(1+sort(ppsList).^2), 'r-', 'LineWidth',1.5);
set(gca,'XScale','log'); xlabel('pixels per speckle (s/p)'); ylabel('spatial contrast K');
legend('simulated','theory (s/p)/\surd(1+(s/p)^2)','Location','southeast');
title('Contrast vs pixels-per-speckle (identical illumination)'); set(gcf,'color','w');


%% STEP 5 - PERFECT SPECKLE OF 5 GRAIN SIZES, SAVED AS .rls RECORDINGS
close all; clearvars -except libraryFolder outputFolder
% Generate a large fully-decorrelated stack at 8 pixels/speckle, derive the
% physically-correct coarser versions (4, 2, 1, 0.5 px/speckle) by block-
% averaging the intensity of the SAME field, centre-crop to the requested size,
% scale each to fill 0-255 and write it as a .rls recording. You get one file
% per grain size (all baseN/16 px), plus the 2 px/speckle set additionally at
% 1024 and 2048 px (2048 is its full, un-cropped block-averaged size). The
% illumination width is set wide so the envelope stays flat over the FOV.
% WARNING: at the 8192 default this renders one 8192x8192 frame at a time and
% keeps a ~1.7 GB 2048-px stack; expect several minutes (~5-10 min) on a strong
% GPU. Reduce baseN / nFrames / particlesN while experimenting (baseN=1024).
baseN      = 8192;               % finest set (8 px/speckle) is baseN x baseN
nFrames    = 100;                % frames per recording
particlesN = 1000;               % scattering particles
illumWidth = 8000;               % illuminated width (um); wide -> flat envelope

% Output recordings: name, block-average factor (grain = 8/factor px/speckle),
% and the square crop size (clamped to what baseN allows). The 2 px/speckle set
% (factor 4) is written at 512, 1024 and its full 2048 px.
specs = { 'simulated_8sp',       1,  512
          'simulated_4sp',       2,  512
          'simulated_2sp',       4,  512
          'simulated_2sp_1024',  4, 1024
          'simulated_2sp_2048',  4, 2048
          'simulated_1sp',       8,  512
          'simulated_05sp',     16,  512 };
binVals = [1 2 4 8 16];          % distinct block-average factors to compute once

% Fix the speckle size once (it depends only on the geometry) so every frame
% shares exactly the same speckle-to-pixel ratio:
[~, calInfo] = getDynamicSpeckles('pixelsN',256,'nFrames',1,'pixelsPerSpeckle',8, ...
    'sampPerSpeckle',8,'motionType','decorrelated','particlesN',particlesN, ...
    'sizeX',illumWidth,'sizeY',illumWidth,'rngSeed',100);
speckleSize = calInfo.speckleSize;

nSpec = size(specs,1);
outN  = zeros(1,nSpec);
sets  = cell(1,nSpec);
for s = 1:nSpec
    outN(s) = min(specs{s,3}, baseN/specs{s,2});     % clamp crop to available size
    sets{s} = zeros(outN(s),outN(s),nFrames,'single');
end

rng(100);                        % reproducible, independent frames
for f = 1:nFrames
    frame = getDynamicSpeckles('pixelsN',baseN,'nFrames',1,'pixelsPerSpeckle',8, ...
        'sampPerSpeckle',8,'speckleSize',speckleSize,'calibrateSpeckle',false, ...
        'motionType','decorrelated','particlesN',particlesN, ...
        'sizeX',illumWidth,'sizeY',illumWidth,'rngSeed',[]);
    ds = cell(1,numel(binVals));
    for b = 1:numel(binVals), ds{b} = blockBinStack(frame, binVals(b)); end
    for s = 1:nSpec
        bIdx = find(binVals==specs{s,2}, 1);
        sets{s}(:,:,f) = centreCrop(ds{bIdx}, outN(s));
    end
    if mod(f,10)==0 || f==nFrames, fprintf('STEP5 frame %d/%d\n', f, nFrames); end
end

% Scale each set to fill 0-255 (uint8) and write it as a .rls recording:
for s = 1:nSpec
    mn = min(sets{s}(:)); mx = max(sets{s}(:));
    sets{s} = uint8(round((sets{s}-mn)/(mx-mn)*255));
    writeRLS(fullfile(outputFolder,[specs{s,1},'.rls']), sets{s}, 100);
end

% Quality control (the five grain sizes at 512 px): single frame (top row) and
% mean over frames (bottom row).
qcIdx   = [1 2 3 6 7];
qcNames = {'8sp','4sp','2sp','1sp','05sp'};
figure('color','w');
for j = 1:5
    subplot(2,5,j);   imagesc(sets{qcIdx(j)}(:,:,1));          axis image off; title([qcNames{j} ' frame 1']);
    subplot(2,5,5+j); imagesc(mean(single(sets{qcIdx(j)}),3)); axis image off; title([qcNames{j} ' mean']);
end
colormap gray;


%% LOCAL FUNCTIONS
function S = blockBinStack(I, b)
% Block-average b-by-b pixels per frame (area integration of the intensity).
if b == 1, S = I; return; end
[y, x, t] = size(I);
S = reshape(mean(reshape(I, b, y/b, x, t), 1), y/b, x, t);
S = permute(S, [2 1 3]);
S = reshape(mean(reshape(S, b, x/b, y/b, t), 1), x/b, y/b, t);
S = permute(S, [2 1 3]);
end

function C = centreCrop(A, n)
% Central n-by-n crop of each frame of A (n clamped to the frame size).
[h, w, ~] = size(A);
n  = min([n, h, w]);
r0 = floor((h-n)/2) + 1;
c0 = floor((w-n)/2) + 1;
C  = A(r0:r0+n-1, c0:c0+n-1, :);
end
