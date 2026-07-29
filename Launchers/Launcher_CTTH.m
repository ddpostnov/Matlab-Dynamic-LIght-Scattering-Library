% Launcher_CTTH - Bolus-injection (CTTH) analysis from wide-field fluorescence.
%
%   Example launcher for bolus-injection analysis of wide-field fluorescence
%   data.  Run STEP 0 once per MATLAB session, then run the step cells (%%) in
%   order.
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
rootFolder = 'C:\Data\mia'; %root folder for the files lookup


%% STEP 1 Process .rls files to get the temoiral contrast for segmentation and vasomotion analysis
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

s.fBolus=[];%[301,1500]; %expected frame indexes for bolus injection, leave empty if unknown
s.fAngio=[];%[1600,2000]; % expected frame indexes for angiogram, leave empty if unknown

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*BB0.cxd'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runBolus(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Define segmentation regions (interactive ROI editor; optional - whole window if skipped)
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%REGION SELECTION - setRegions is fully interactive: it opens an ROI editor per file
%(Add ROI / Delete ROI / Reset ROIs + polygon/rectangle/square/ellipse/circle shape
%selector; select an ROI and press Delete to remove it).  The number of regions is
%however many you draw - there is no count to set - and nothing advances until you
%press Done.  Draw nothing, or skip this step, to keep the whole window (no region
%mask is written).

%SET FILE NAMES HERE - passed as one ROW (transpose) so a drawn ROI carries (editable)
%across all bolus files as a single group (the intensity _I path has no grouping id)
fNames=getFileNamesList(rootFolder,'*_b_I_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

setRegions(s,fNames(:)'); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 3 Perform segmentation (categories + labels + per-segment traces; fully automatic)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - CATEGORIZATION ADJUSTEMNTS (intensity _I path: no trust limits)
s.lSizeN=121; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=5; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.2; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
s.sSizeScale=1; % scaler for small objects assignment to background or to unregognized regions
s.deSens=1; %can be used to reduce sensitivity to small objects
s.lThinN=2; % Large vessels thinning
s.imOpen=0; % Small vessels thinning
s.iEdge=1; %setting internal edges for segmented vessels
s.eEdge=1; %setting external edges for segmented vessels

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}; %CATEGORIES

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%SET FILE NAMES HERE - FLAT (order-independent; grouping was setRegions' job in STEP 2)
fNames=getFileNamesList(rootFolder,'*_b_I_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 4 (OPTIONAL. Only use if 1 or more regions were defined in STEP 2) Split the regions.
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_b_I_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%USES THE SAME FILE NAMES AS ABOVE as STEP 2 (crops each file by its own regionsMask -> RoiN_ files)
splitRegions(s,fNames(:)); %LAUNCHES THE UTILITY ROUTINE

%% STEP 5 (OPTIONAL) Dynamic segmentation - per-frame vessel diameter / flow (heavy)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 3)
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
s.sMaxLBI=(1/5)./s.sMinL; %Max local bending (0 to pi per pixel)
s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
s.sMaxKK=0.3; %Max accepted std/mean for the initial contrast estimation
s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate

%ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
s.gSizeN=3;
s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
s.pInterpF=4; % leave as is

%SET FILE NAMES HERE (after STEP 4 this pattern also matches the RoiN_ crops)
fNames=getFileNamesList(rootFolder,'*_b_I_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runDynamicSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud.

files      = dir(fullfile(rootFolder,'**','*t_K_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_K_d.mat$', '_c_K_d.mat');
fNames=cat(1,fNamesPLS,fNamesVSM);
runRegistration(s,fNames); %LAUNCHES THE UTILITY ROUTINE

%%
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;


    s.calcData="all";                 % per-region + per-pixel
    s.medSpace=3; s.medTime=3;         % despeckle / despike
    s.sgOrder=2;  s.sgFrame=11;        % smoothing
    s.hampWin=3;  s.hampSig=3;         % motion-spike rejection
    s.slopeWin=9; s.promFrac=0.2;      % robust upslope / peak
    s.minStep=1;                       % uint16 resolution

fNames=getFileNamesList(rootFolder,'*_b_I_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.


runCTTH(s,fNames(:))

%% STEP 8 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*_K_d.mat'));
fNames     = fullfile({files.folder}', {files.name}');

runBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Perform vasomotion analysis
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
s.vFR=[0.05,0.25];
s.cFR=[0.4,0.6];
s.wFR=[0.01,1];
s.wVPO=10;
s.normalisation='median'; %'mean'/'median' (global) or 'mmean'/'mmedian' (moving). Change only if needed
s.normsize=101; %window for 'mmean'/'mmedian'; inf or 0 -> global. Change only if needed
s.tgtFS=1; %Hz
s.pcts=0:10:100;

s.otsuMaxN=5;
s.otsuElbow= 0.05;
s.nPeakProm=0.10;   %VB peak-count prominence = fraction of per-time band max-min range

s.ppxVsmReturn = [];         %per-pixel analysis OFF ([] = off); set non-empty (e.g. {'bands'}) to turn it ON -> RESULTS.vasomotion.ppx. No analysePerPixel flag.
s.ppxSegmentAveraging=[]; %TEMPORARY scaffolding (to be removed): per-segment averaging demo, subset of {'coherent','incoherent'}; [] = off. Change only if needed
s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'};  %which data-type signals get vasomotion analysis
%s.segVsmReturn selects which levels to store in results.vasomotion: 'bands' (scalars.VB/CB),
%'moments' (fVectors ampMean/Std/Skew + VB.ampMeanPct), 'series' (timeVectors amp/fCent/fSprd/nPeak),
%'clustering' (flare/silence scalars+spectra+maskFlare), 'reconstruction' (timeVectors.VB.rData),
%'spectrum' (spectrum.amp/.phase grid). Change only if needed.
s.segVsmReturn={'bands','moments','series','clustering','reconstruction'};

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*t_BFI_d.mat'));
fNames     = fullfile({files.folder}', {files.name}');

runVasomotion(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 10 Perform pulsatility analysis (strictly after conversion to BFI)
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED (OR VERIFIED) PER PROTOCOL - Harmonic model
s.nHarm=5;              % number of harmonics in y=SUM a_n*sin(2*pi*n*x+b_n)+c
%ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
s.segPulsReturn={'markers','model','reconstruction'};  % markers = model-free scalars
                        % (PI/RI/mean/extrema/timing/symmetry); model = harmonic
                        % hAmp/hPhase/R2 (runs the fit); reconstruction = the
                        % timeVectors.fData model cube (runs the fit). Default = all three.
%ADJUSTED IF NECESSARY - Per-pixel maps (GATE + selector; [] = per-pixel analysis off)
s.ppxPulsReturn={'markers'};  % NON-EMPTY = per-pixel marker maps ON
                        % (results.pulsatility.ppx); add 'model'/'reconstruction' to
                        % also fit every masked pixel (large full-resolution cubes).

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*c_BFI_d.mat')); %<---ALWAYS REFER TO "c_BFI_d.mat" files, but you may use regexp to define specific "c_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

runPulsatility(s, fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 11 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
s.useReference=true; 
s.libraryFolder=libraryFolder;

files      = dir(fullfile(rootFolder,'**','*t_BFI_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_BFI_d.mat$', '_c_BFI_d.mat');
fNames=cat(1,fNamesPLS,fNamesVSM);
s.refFName=fNames{1};
setVesselTypes(s,fNames);

%% STEP 11 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

%%

figure
plot(squeeze(max(source.data,[],[1,2])))
hold on
plot(squeeze(mean(source.data,[1,2])))
hold off

%%

dataBL=source.data(:,:,1:575);
dataRP=source.data(:,:,576:end-100);
dataFN=source.data(:,:,end-99:end);

figure
cval=[mean(dataRP(:,:,1),'all'),prctile(max(dataRP,[],3),99,'all')];
for i=1:1:size(dataRP,3)
imagesc(dataRP(:,:,i))
clim(cval)
axis image
colormap gray
title(num2str(i/100))
pause(0.05)
end


dataRPN=single(dataRP)-mean(single(dataBL),3);%
dataRPN=mat2gray(dataRPN);
fSize=151;
parfor i=1:1:size(dataRPN,3)
i
dataRPN(:,:,i)=dataRPN(:,:,i)-imopen(medfilt2(dataRPN(:,:,i),[15,15],'symmetric'),strel('disk',fSize));
end

figure
cval=[mean(dataRPN(:,:,1),'all'),prctile(max(dataRPN,[],3),99,'all')];
for i=1:1:size(dataRPN,3)
imagesc(dataRPN(:,:,i))
clim(cval)
axis image
colormap gray
title(num2str(i/100))
pause(0.5)
end
%%
img=img-imopen(medfilt2(img,[15,15],'symmetric'),strel('disk',fSize));
data=applyDirectionalFilter(single(dataRP), img);
[fitP, fitQ, fitG,bMin,bMax] = fitGammaVariatePerPixel(data, time);
%%
figure
imagesc(mean(data,3))
mask=roipoly;


%%
figure
img=squeeze(MTT);
imagesc(img)
clim(prctile(img(mask(:)),[0.1,99.9]))
colorbar
axis image
colormap jet
title('MTT')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);

%%
figure
img=squeeze(fitP(:,:,3));
imagesc(img)
clim([0,2.5])
colorbar
axis image
colormap jet
title('Arrival time')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);

figure
img=squeeze(fitP(:,:,5));
imagesc(img)
clim(prctile(img(mask(:)),[0.1,99.9]))
colorbar
axis image
colormap jet
title('Peak time')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);

figure
img=squeeze(fitP(:,:,4));
imagesc(img)
clim([1,4])
colorbar
axis image
colormap jet
title('Alpha')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);


%%
files      = dir(fullfile(rootFolder,'**','*t_BFI_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_BFI_d.mat$', '_c_BFI_d.mat');

clearvars vsm
for i=1:1:numel(fNamesVSM)
load(strrep(fNamesVSM{i},'_d.mat','_r.mat'),'results');
vsm(i)=results;
clearvars results
end

clearvars pls
for i=1:1:numel(fNamesPLS)
load(strrep(fNamesPLS{i},'_d.mat','_r.mat'),'results');
pls(i)=results;
clearvars results
end

f=vsm(1).vasomotion.f;
pctCenters=vsm(1).vasomotion.pctCenters;
types={'Arteries','Parenchyma','Veins'};
conds={'AWK','ISO','K/X'};

pctSpc=zeros(3,3,numel(f),numel(pctCenters));

ampVB=zeros(numel(vsm),3);
ampCB=zeros(numel(vsm),3);


for i=1:1:numel(vsm)
    ampVB(i,1)= mean(vsm(i).sMetrics.ampMeanVB(strcmp(vsm(i).sMetrics.type,"Artery")));
    ampVB(i,2)=mean(vsm(i).sMetrics.ampMeanVB(vsm(i).sMetrics.category==1));
    ampVB(i,3)= mean(vsm(i).sMetrics.ampMeanVB(strcmp(vsm(i).sMetrics.type,"Vein")));

    ampCB(i,1)= mean(vsm(i).sMetrics.ampMeanCB(strcmp(vsm(i).sMetrics.type,"Artery")));
    ampCB(i,2)=mean(vsm(i).sMetrics.ampMeanCB(vsm(i).sMetrics.category==1));
    ampCB(i,3)= mean(vsm(i).sMetrics.ampMeanCB(strcmp(vsm(i).sMetrics.type,"Vein")));

    pctSpc(i,1,:,:)= mean(vsm(i).vasomotion.sData.fVectors.VB.ampMeanPct(strcmp(vsm(i).sMetrics.type,"Artery"),:,:),1);
    pctSpc(i,2,:,:)= mean(vsm(i).vasomotion.sData.fVectors.VB.ampMeanPct(vsm(i).sMetrics.category==1,:,:),1);
    pctSpc(i,3,:,:)= mean(vsm(i).vasomotion.sData.fVectors.VB.ampMeanPct(strcmp(vsm(i).sMetrics.type,"Vein"),:,:),1);
end

figure
for i=1:1:3
for j=1:1:3
subplot(3,3,(i-1)*3+j)
plot(f,squeeze(pctSpc(i,j,:,:)))
title([conds{i},', ',types{j}])
ylim([0,0.2])
xlabel('Frequency, Hz')
ylabel('WT Amplitude')
end
end
set(gcf,'Color','w');

figure
for i=1:1:3
subplot(1,3,i)
plot(squeeze(ampVB(i,:)))
hold on
plot(squeeze(ampCB(i,:)))
hold off
xlim([0,4])
ylabel('Amplitude density')
legend({'VSM: 0.05-0.25Hz','CTR: 0.4-0.6Hz'})
xticks([1,2,3]);
xticklabels(types)
title(conds{i})
ylim([0.01,0.07])
end
set(gcf,'Color','w');
