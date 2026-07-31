% Launcher_guided - Guided full-resolution per-segment traces (demo on test data).
%
%   Demonstrates the "guided" step of the pipeline on the bundled mouse-brain
%   test recording.  The standard pipeline turns the raw recording into
%   temporally decimated contrast and segments it (results.sMap / results.sData
%   at the decimated frame rate).  The guided step then goes back to the raw
%   recording and, using results.sMap as a guide, extracts one trace per segment
%   at the FULL frame rate (results.gsData / results.gsTime).
%
%   Run STEP 0 once per MATLAB session, then run the step cells (%%) in order.
%   STEP 1-3 reproduce the usual contrast -> regions -> segmentation chain,
%   STEP 4 is the guided step, and STEP 5-7 show the results (overview, a vessel
%   trace at full vs decimated resolution, and a dynamic perfusion movie
%   reconstructed from gsData + sMap).  Everything runs on the test data with no
%   editing required.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED THROUGH THE STEPS
% LIBRARY PATH - auto-detected from this file's location (this launcher lives in
% <library>\Launchers).  Replace with your own path if you prefer.
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library

% The demo runs on the bundled test recording.  Point rawName at your own .rls
% (or .cxd) file to use the guided step on your data.
rootFolder = fullfile(libraryFolder,'Test data');
rawName    = fullfile(rootFolder,'mouseBrain_crop512x512x1000.rls');

% Use the GPU for contrast if one is available, otherwise the CPU.
procType = 'cpu'; if gpuDeviceCount>0, procType = 'gpu'; end

%% STEP 1 Process the raw recording into temporally decimated contrast
clearvars -except libraryFolder rootFolder rawName procType
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.contrastType='temporal'; %'temporal' or 'spatial'
s.contrastKernel=25; %typical values: 25 for 'temporal', 5 or 7 for 'spatial'
s.decimFactor=25; %decimates the contrast. Output framerate = original framerate / s.decimFactor
s.decimMethod='sharp'; %'sharp' (temporal only, decimFactor multiple of kernel) or 'leaking'

%ADJUSTED IF NECESSARY - PERFORMANCE
s.procType=procType; %'gpu' if a high-end GPU is available, 'cpu' otherwise

%ADJUSTED IF NECESSARY - INITIAL MASKING
s.trustLimitsK=[0.001,0.8]; %min (fastest flows) and max (slowest flows) expected contrast
s.trustLimitsI=[3,254]; %min and max expected intensity
s.minTrust=[0.7,0.7]; %per-pixel fraction of non-zero / non-saturated frames required
s.manualMask=0; %allows manual subselection of the area to mask

runContrastFromRLS(s,{rawName}); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 1b (OPTIONAL) Collect the report pages of STEP 1 into one PDF
close all
clearvars -except libraryFolder rootFolder rawName procType

%A wrapper writes one <recording>_rep_<page>.jpg beside each recording and stops
%there - the document is assembled HERE, between the steps, so a 60-file step is one
%thing to page through instead of 60 files.  COPY THIS CELL after any step whose
%pages you want together and change the tail and the output name: after setRegions
%it is '_rep_regions.jpg', after runSegmentation '_rep_categories.jpg' and
%'_rep_segments.jpg', after runRegistration '_rep_registration.jpg', after
%setVesselTypes '_rep_vesseltypes.jpg'.
tail='_rep_contrast.jpg'; %the page STEP 1 writes, one per recording
pdfFolder=fileparts(rawName);
D=dir(fullfile(pdfFolder,['*',tail]));
makeReportPdf(fullfile({D.folder}',{D.name}'),fullfile(pdfFolder,'report_contrast.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 2 (OPTIONAL) Define segmentation regions - this demo uses the WHOLE window
close all
clearvars -except libraryFolder rootFolder rawName procType
s.libraryFolder=libraryFolder;

% Segmentation runs on the whole field of view by default, so the demo does nothing
% here.  To restrict it to one or more sub-regions, run setRegions: it opens an
% interactive editor (Add / Delete / Reset ROI + polygon/rectangle/square/ellipse/
% circle shape selector) and writes results.regionsMask.  The number of regions is
% however many you draw - draw nothing, or skip this step entirely (as here), to keep
% the whole window.  Uncomment to try it:
%
% dName=strrep(rawName,'.rls','_t_K_d.mat');
% setRegions(s,{dName}); %LAUNCHES THE INTERACTIVE ROI EDITOR

%% STEP 3 Perform segmentation (categories + labels; builds results.sMap / sMetrics / sData)
close all
clearvars -except libraryFolder rootFolder rawName procType
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CATEGORIZATION
s.trustLimitsK=[0.001,0.99];
s.lSizeN=61; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=7;  % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.3;  % Segmentation sensitivity - increase if missing vessels
s.sSizeScale=1; % scaler for small objects assignment to background or unrecognised regions
s.deSens=1;  % reduce sensitivity to small objects
s.lThinN=2;  % Large vessels thinning
s.imOpen=0;  % Small vessels thinning
s.iEdge=2;   % internal edges for segmented vessels
s.eEdge=2;   % external edges for segmented vessels
s.iniSizeN=7;% Odd number equal or larger than the spatial contrast kernel

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'};

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
s.sStat='median';    % 'median' (default) or 'mean' statistic for per-segment traces
s.sMinL=15;          % Minimum length for segments
s.prchNSize=50;      % Parenchymal pixels neighbourhood
s.correctNodes=true; % Branching correction
s.simR=0.3;          % minimal similarity ratio between branches to be the same vessel
s.difR=0.4;          % minimal difference ratio to be considered different vessels

dName=strrep(rawName,'.rls','_t_K_d.mat');
runSegmentation(s,{dName}); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 4 GUIDED STEP - full-resolution per-segment contrast from the raw file
close all
clearvars -except libraryFolder rootFolder rawName procType
s.libraryFolder=libraryFolder;

dName=strrep(rawName,'.rls','_t_K_d.mat');

% Pass the raw file explicitly (it is also auto-detected if omitted).  Use
% runGuidedIntensity instead for per-segment mean intensity (e.g. fluorescence).
runGuidedContrast(s,{dName},{rawName}); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 5 See the results - segmentation overview
close all
clearvars -except libraryFolder rootFolder rawName procType
dName=strrep(rawName,'.rls','_t_K_d.mat');
load(strrep(dName,'_d.mat','_r.mat'),'results');

sMap=single(results.sMap); sMap(sMap==0)=NaN;
imgK=results.imgK;                      % mean contrast (dark = fast flow = vessels)

% Golden-angle colour map so neighbouring segments look different
n=double(max(results.sMap(:)));
phi=(sqrt(5)-1)/2;
cmap=hsv2rgb([mod((0:n-1)'*phi,1) 0.8*ones(n,1) 0.9-0.2*mod((0:n-1)',2)]);

f=figure('Color','w','Name','Guided - segmentation overview');
f.WindowState='maximized';
t=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax1=nexttile(t);
imagesc(ax1,imcomplement(mat2gray(imgK,double(prctile(imgK(:),[1 99])))));
axis(ax1,'image'); colormap(ax1,'gray'); title(ax1,'Mean contrast (vessels bright)');
hold(ax1,'on'); visboundaries(ax1,results.sMap>0 & results.cMask==5,'Color','c','LineWidth',0.5); hold(ax1,'off');
ax2=nexttile(t);
h=imagesc(ax2,sMap); set(h,'AlphaData',~isnan(sMap));
axis(ax2,'image'); colormap(ax2,cmap); set(ax2,'Color',[1 1 1]);
title(ax2,sprintf('%d segments (%d vessels)',n,nnz(results.sMetrics.category==5)));
drawnow
print(f,strrep(dName,'_d.mat','_guided_overview.jpg'),'-djpeg','-r200');

%% STEP 6 See the results - one vessel at full vs decimated resolution
clearvars -except libraryFolder rootFolder rawName procType
dName=strrep(rawName,'.rls','_t_K_d.mat');
load(strrep(dName,'_d.mat','_r.mat'),'results');

% Pick the largest lumen (vessel) segment
isLumen=results.sMetrics.category==5;
area=results.sMetrics.area; area(~isLumen)=-Inf;
[~,li]=max(area);

% results.gsType tells us whether gsData holds contrast or intensity, so the
% plot always shows the right quantity (contrast -> BFI = 1/K^2).
isContrast = ~isfield(results,'gsType') || strcmp(string(results.gsType),"contrast");
tHD=results.gsTime;                                       % full temporal resolution
if isContrast
    yHD=1./results.gsData(:,li).^2; yLab='BFI (1/K^2)';
    yLD=1./results.sData(:,li).^2;  tLD=results.time;     % decimated (segmentation) rate
    guidedName=sprintf('Guided full-res, single-frame spatial contrast (%.0f fps)',1/median(diff(tHD)));
else
    yHD=results.gsData(:,li);        yLab='Mean intensity';
    yLD=[]; tLD=[];
    guidedName=sprintf('Guided full-res mean intensity (%.0f fps)',1/median(diff(tHD)));
end

% Simple amplitude spectrum of the full-resolution trace (reveals the heart rate)
x=yHD-mean(yHD,'omitnan'); x(~isfinite(x))=0;
fs=1/median(diff(tHD)); N=numel(x); fAx=(0:floor(N/2))*fs/N;
P=abs(fft(x)); P=P(1:numel(fAx));

f=figure('Color','w','Name','Guided - vessel trace'); f.WindowState='maximized';
t=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
ax1=nexttile(t);
plot(ax1,tHD,yHD,'-','Color',[0 0.45 0.74],'LineWidth',1); hold(ax1,'on');
if ~isempty(yLD)
    plot(ax1,tLD,yLD,'-o','Color',[0.85 0.33 0.10],'LineWidth',1.5,'MarkerSize',4);
    legend(ax1,{guidedName,sprintf('Decimated sData, temporal contrast (%.1f fps)',1/median(diff(tLD)))},'Location','best');
    subtitle(ax1,['The two use different contrast estimators, so absolute BFI differs; ',...
        'the envelopes agree while only the guided trace captures the heartbeat.']);
else
    legend(ax1,{guidedName},'Location','best');
end
hold(ax1,'off'); grid(ax1,'on'); xlim(ax1,[tHD(1) tHD(end)]);
xlabel(ax1,'Time, s'); ylabel(ax1,yLab);
title(ax1,sprintf('Segment %d - the guided trace resolves pulsatility',li));
ax2=nexttile(t);
plot(ax2,fAx,P,'k','LineWidth',1); grid(ax2,'on'); xlim(ax2,[0 min(20,fs/2)]);
xlabel(ax2,'Frequency, Hz'); ylabel(ax2,'Amplitude'); title(ax2,'Spectrum of the full-resolution trace');
drawnow
print(f,strrep(dName,'_d.mat','_guided_trace.jpg'),'-djpeg','-r200');

%% STEP 7 See the results - dynamic perfusion movie (from gsData + sMap)
clearvars -except libraryFolder rootFolder rawName procType
dName=strrep(rawName,'.rls','_t_K_d.mat');
load(strrep(dName,'_d.mat','_r.mat'),'results');

%EDIT THESE IF YOU LIKE
frameStride=2;   % show every Nth frame (2 keeps the movie short)
fps=25;          % playback / video frames per second
smoothWin=5;     % cosmetic temporal smoothing (frames); 1 = none
saveVideo=true;  % also write an MP4 next to the data

% Reconstruct the per-segment field and paint it back onto the segments.
% results.gsType decides what is shown: contrast -> perfusion (BFI = 1/K^2),
% intensity -> mean intensity.
isContrast = ~isfield(results,'gsType') || strcmp(string(results.gsType),"contrast");
if isContrast
    field=1./results.gsData.^2; cbLab='BFI (1/K^2)'; titlePrefix='Dynamic perfusion';
else
    field=results.gsData;       cbLab='Mean intensity'; titlePrefix='Dynamic intensity';
end
if smoothWin>1, field=movmean(field,smoothWin,1,'omitnan'); end
valid=results.sMap>0 & results.mask;
lut=zeros(size(results.sMap)); lut(valid)=double(results.sMap(valid)); % 0 = background
clim0=prctile(field(isfinite(field)),[2 98]);
if ~(clim0(2)>clim0(1)), clim0=[min(field(:)) max(field(:))]; end

frames=1:frameStride:size(field,1);
if saveVideo
    vidFile=strrep(dName,'_d.mat','_guided_perfusion.mp4');
    try
        vw=VideoWriter(vidFile,'MPEG-4'); vw.FrameRate=fps; open(vw);
    catch ME   % e.g. the file is open elsewhere or locked by a sync client
        warning('Could not open %s (%s). Playing the animation without saving.',vidFile,ME.message);
        saveVideo=false;
    end
end
% A fixed figure size keeps every captured frame identical (required by
% VideoWriter); maximised/resizable windows produce varying frame sizes.
f=figure('Color','k','Name','Guided - dynamic perfusion',...
    'MenuBar','none','ToolBar','none','Position',[100 100 1000 900]);
ax=axes(f);
img0=reshape([NaN,field(frames(1),:)],1,[]); img0=reshape(img0(lut+1),size(lut));
imH=imagesc(ax,img0); set(imH,'AlphaData',~isnan(img0));
axis(ax,'image'); axis(ax,'off'); colormap(ax,'turbo'); clim(ax,clim0);
set(ax,'Color','k'); ax.Toolbar.Visible='off';
cb=colorbar(ax); cb.Color='w'; cb.Label.String=cbLab;
ttl=title(ax,'','Color','w');
for k=frames
    row=[NaN,field(k,:)];
    img=reshape(row(lut+1),size(lut));
    set(imH,'CData',img,'AlphaData',~isnan(img));
    ttl.String=sprintf('%s   t = %.2f s   (frame %d / %d)',...
        titlePrefix,results.gsTime(k),k,size(field,1));
    drawnow;
    if saveVideo, writeVideo(vw,getframe(f)); else, pause(1/fps); end
end
if saveVideo, close(vw); fprintf('Saved perfusion movie to %s\n',vidFile); end

% %% STEP 8 (OPTIONAL) Guided mean INTENSITY instead of contrast
% % This replaces results.gsData with intensity and sets results.gsType="intensity".
% % STEP 6-7 read gsType, so re-running them now shows intensity (not 1/K^2) correctly.
% % Re-run STEP 4 to restore the contrast/BFI perfusion movie.
% close all
% clearvars -except libraryFolder rootFolder rawName procType
% s.libraryFolder=libraryFolder;
% dName=strrep(rawName,'.rls','_t_K_d.mat');
% runGuidedIntensity(s,{dName},{rawName}); %LAUNCHES THE PROCESSING ROUTINE
