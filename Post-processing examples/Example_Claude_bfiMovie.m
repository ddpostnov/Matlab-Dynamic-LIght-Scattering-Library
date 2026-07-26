% Example_Claude_bfiMovie - Play back (and optionally save) a BFI recording as a movie.
%
% WHAT THIS SCRIPT DOES
%   Loads one processed dataset (*_BFI_d.mat), fixes the colour limits once so
%   the scale does not flicker, and plays the blood-flow-index (BFI) images as a
%   movie with the current time shown in the title. It can optionally save the
%   movie to an MP4 file next to your data.
%
% WHEN TO USE IT
%   The best first step with any new recording: watch it to get a feel for the
%   data, spot movement or artefacts, and see where and when flow changes.
%
% HOW TO USE IT
%   1. Set fName to your *_BFI_d.mat file.
%   2. Optionally set frameStride (skip frames for speed) and saveVideo.
%   3. Run the script.
%
% GENERATIVE-AI PROMPT  (paste into Claude / ChatGPT to regenerate a similar script)
%   "Write a simple, well-commented MATLAB script for users who are not good
%    with MATLAB. It loads a .mat file with a struct 'source' containing a 3-D
%    array source.data (Y x X x Time, blood-flow index) and a time vector
%    source.time in seconds, and optionally a struct 'results' with
%    results.cMask (a tissue mask). Compute robust colour limits once (ignore
%    NaN/Inf and the top/bottom 1% of tissue pixels) so the scale is fixed for
%    the whole movie. Play the frames with the current time in the title, let
%    the user skip frames for speed, and give an option to save the movie as an
%    MP4 file."
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.

%% 0. Point MATLAB to the library -----------------------------------------
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));

%% 1. Settings - EDIT THESE -----------------------------------------------
fName       = 'C:\path\to\your\data\recording_t_BFI_d.mat';
frameStride = 1;      % show every Nth frame (use 2 or 5 for long recordings)
saveVideo   = false;  % true -> also write an MP4 next to the data file
fps         = 20;     % playback / video frames per second

%% 2. Load data ------------------------------------------------------------
load(fName,'source');
resultsFile = strrep(fName,'_d.mat','_r.mat');
if isfile(resultsFile), load(resultsFile,'results'); else, results = struct(); end

t  = source.time(:);
nT = numel(t);

%% 3. Robust colour limits (computed once from a few sample frames) -------
sampleIdx = round(linspace(1,nT,min(nT,20)));
sample    = source.data(:,:,sampleIdx);
if isfield(results,'cMask')
    m = repmat(results.cMask>0,1,1,numel(sampleIdx));
else
    m = isfinite(sample) & sample>0;
end
clim0 = prctile(sample(m & isfinite(sample)),[1 99]);
if ~(clim0(2)>clim0(1))                      % guard against a flat range
    clim0 = [min(sample(:)) max(sample(:))];
end

%% 4. Play (and optionally record) ----------------------------------------
if saveVideo
    vidFile = strrep(fName,'_d.mat','_movie.mp4');       %#ok<UNRCH>
    vw = VideoWriter(vidFile,'MPEG-4'); vw.FrameRate = fps; open(vw);
end

figH = figure('Color','w');
axH  = axes(figH);
imH  = imagesc(axH,source.data(:,:,1)); axis(axH,'image'); clim(axH,clim0);
colormap(axH,parula); colorbar(axH)
for k = 1:frameStride:nT
    set(imH,'CData',source.data(:,:,k));
    title(axH,sprintf('BFI    t = %.2f s    (frame %d / %d)',t(k),k,nT));
    drawnow;
    if saveVideo
        writeVideo(vw,getframe(figH));                   %#ok<UNRCH>

    else
        pause(1/fps);                        % play at roughly fps
    end
end

if saveVideo, close(vw); fprintf('Saved movie to %s\n',vidFile); end  %#ok<UNRCH>
