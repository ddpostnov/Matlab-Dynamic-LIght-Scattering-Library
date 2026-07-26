% Launcher_vasoreactivity - Vasoreactivity / vasomotion analysis pipeline.
%
%   Example launcher for vasoreactivity or vasomotion analysis; for vasomotion,
%   >=25 fps raw data is advised.  Run STEP 0 once per MATLAB session, then run
%   the step cells (%%) in order.
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
rootFolder = 'C:\Dropbox\Work\Data'; %root folder for the files lookup

%% STEP 1 Process .rls files to get the contrast
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.contrastType='temporal'; %'temporal' or 'spatial'
s.contrastKernel=25; %typical values: 25 for 'temporal', 5 or 7 for 'spatial'
s.decimFactor=25; %decimates the contrast. Output framerate = original framerate / s.decimation
s.decimMethod='sharp'; %or  s.decimationMethod='leaking'; 'sharp' is only for temporal analysis and and s.decimation being a multiple integer of s.contrastKernel

%ADJUSTED IF NECESSARY - PERFORMANCE ADJUSTEMNTS
s.procType='gpu'; %use 'gpu' for spatial contrast type if high-end GPU is availible, 'cpu' otherwise

%ADJUSTED IF NECESSARY - INITIAL MASKING PARAMETERS
s.trustLimitsK=[0.001,0.5]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke
s.trustLimitsI=[5,250]; %minimum (first value) and maximum (second value) of expected intensity.
s.minTrust=[0.99,0.99]; %per-pixel trust limits in relation to the portion of frames with minimum (0) or maximum (usually 255) intensity.
s.manualMask=0; %allows manual subselection of the area to mask

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*BV.rls'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runContrast(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Define pixel categories (based on temporal contrast data)
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.trustLimitsK=[0.001,0.5]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke

%ADJUSTED IF NECESSARY - SEGMENTATION ADJUSTEMNTS
s.regionsN=1; %Numer of regions for manual selection. 0 if using entire window.
s.lSizeN=121; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=7; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.3; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
s.sSizeScale=1; % scaler for small objects assignment to background or to unregognized regions
s.deSens=1; %can be used to reduce sensitivity to small objects
s.lThinN=2; % Large vessels thinning 
s.imOpen=0; % Small vessels thinning 
s.iEdge=3; %setting internal edges for segmented vessels
s.eEdge=3; %setting external edges for segmented vessels

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}; %CATEGORIES

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat','[A-Z]+\d+'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    runCategories(s,fNames(i,:)');
end


%% STEP 3 (OPTIONAL. Only use if 1 or more regions are defined in step 2) Split the regions. 
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
splitRegions(s,fNames(:));

%% STEP 4 Perform segmentation
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - BASIC PARAMETERS
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.attmemptDS=false; %attempt to perform automated dynamic segmentation or not
s.sMinL=10; % Minimum length for segments
s.prchNSize=30; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
s.sMaxLBI=(1/5)./s.sMinL; %Max local bending (0 to pi per pixel)
s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
s.sMaxDK=0.2; %Max accepted std/mean for the initial diameter estimation
s.sMaxKK=0.3; %Max accepted std/mean for the initial contrast estimation
s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate


%ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
s.gSizeN=3;
s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
s.pInterpF=4; % leave as is

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
runSegmentation(s, fNames(:));

%% STEP 5 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud - same as in the segmentation step

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    runRegistration(s,fNames(i,:)');
end


%% STEP 6 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runBFI(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 (Optional) Perform vasomotion analysis
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

s.ppxVsmReturn = {'bands'};  %per-pixel analysis: NON-EMPTY = ON -> RESULTS.vasomotion.ppx (levels like segVsmReturn), [] = OFF. No analysePerPixel flag.
s.ppxSegmentAveraging=[]; %TEMPORARY scaffolding (to be removed): per-segment averaging demo, subset of {'coherent','incoherent'}; [] = off. Change only if needed
s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'};  %which data-type signals get vasomotion analysis
%s.segVsmReturn selects which levels to store in results.vasomotion: 'bands' (scalars.VB/CB),
%'moments' (fVectors ampMean/Std/Skew + VB.ampMeanPct), 'series' (timeVectors amp/fCent/fSprd/nPeak),
%'clustering' (flare/silence scalars+spectra+maskFlare), 'reconstruction' (timeVectors.VB.rData),
%'spectrum' (spectrum.amp/.phase grid). Change only if needed.
s.segVsmReturn={'bands','moments','series','clustering','reconstruction'};

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runVasomotion(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 8 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%%IF DOING IT FILE BY FILE
% s.useReference=false;
% s.refFName=''; %use '' instead of " "
%fNames=getFileNamesList(rootFolder,'*_c_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
%setVesselTypes(s,fNames(:));

%%IF using a reference
s.useReference=true; %Assumes PRE-registered files
%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_BFI_d.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 8 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

