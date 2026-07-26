% Launcher_Sonam_LSCI - Combined pulsatility and vasomotion pipeline (project-specific).
%
%   Project-specific launcher for combined pulsatility and vasomotion analysis
%   (assumes the corresponding acquisition requirements were met).  Run STEP 0
%   once per MATLAB session, then run the step cells (%%) in order.
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
rootFolder = 'O:\HE_BFI-DATA\Sonam\B3_Black batch\LSCI\T1\FM1'; %root folder for the files lookup


%% STEP 1 Process .rls files to get the temporal contrast for segmentation and vasomotion analysis
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
s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke
s.trustLimitsI=[1,255]; %minimum (first value) and maximum (second value) of expected intensity.
s.minTrust=[0.99,0.99]; %per-pixel trust limits in relation to the portion of frames with minimum (0) or maximum (usually 255) intensity.
s.manualMask=0; %allows manual subselection of the area to mask

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*.rls'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

getContrast(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

% STEP 2 Process .rls files to get the internal cycle data
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke
s.trustLimitsI=[1,255]; %minimum (first value) and maximum (second value) of expected intensity.
s.contrastKernelS=5; %contrast kernel for spatial (sLSCI) processing method
s.maxFrqIni=20; % initial max frequency of the activity of interest, Hz
s.minFrqIni=1; % initial min frequency of the activity of interest, Hz

%ADJUSTED IF NECESSARY - EXCLUSION CRITERIA
s.excludeFirstNCycles=0; %reject given number of cycles
s.coeffsSTD=[3,2,2,2,2,3,3,2,2]; %pulses rejection coefficients relative to the feature standard deviation
s.coeffsRel=[0.5,0.1]; %pulses rejection coefficients relative to the feature value
s.coeffsAbs=2; %pulses rejection coefficients relative to the absolute feature value

%ADJUSTED IF NECESSARY - CYCLE CALCULATION
s.method='sLSCIMM';%,'tLSCIMM','ltLSCIMM' %Typically 'sLSCIMM' is recommended. For high quality data 'ltLSCIMM' will produce better results. Other options are 'tLSCIMM' and 'sLSCIMMM'.
% method refers to spatial, temporal or lossless contrast calculation,
% while the MM or MMM refers to minimum to minimum stretching or minimum to
% maximum + maximum to minimum stretching.
s.decimationSpace=4; %spatial decimation used to conserve memory in the pre-processing steps
s.framesToAverage=1; %allows averaging multiple raw frames to artificially increase expsoure time
s.contrastKernelT=25; %contrast kernel for temporal (tLSCI) and lossless (ltLSCI) processing methods
s.contrastKernelPreproc=s.contrastKernelS; %contrast kernel used in preprocessing (spatial)
s.rangeFrq=1;%1/2; % relative frq range around the central frequency, Hz
s.interpFactor=10; %Sets the number of points that will replace two consequitive values during the interpolation sequence.
s.smoothCoef1=1/3; %in respect to minimum points per cycle value
s.minPromCoef=1/20;%1/2; % in respect to the std of the signal

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*BP.rls'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
getInternalCycle(s,fNames(:));
%% STEP 3 Define pixel categories (optionally can be done on temporal contrast data only)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke

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
    getCategories(s,fNames(i,:)');
end

% STEP 4 Assign same categories to the internal cycle (pulsatility) data 
% (OPTIONAL: used in combination with categories extration based on temporal contrast analysis only)
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
fNames=getFileNamesList(rootFolder,'*c_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
fNamesRef=regexprep(fNames, '\_c_K_d.mat$', '_t_K_d.mat');

assignCategories(fNames,fNamesRef); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 5 Perform segmentation
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - BASIC PARAMETERS
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.attmemptDS=false; %attempt to perform automated dynamic segmentation or not
s.sMinL=15; % Minimum length for segments
s.prchNSize=50; % Parenchymal pixels neighbourhoud.
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

%SET FILE NAMES HERE -
fNames=getFileNamesList(rootFolder,'*_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

getSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=50; % Parenchymal pixels neighbourhoud. Same as in segmentation

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_K_d.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    registerLSCItoLSCI(s,fNames(i,:)');
end

% STEP 8 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

getBFI(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

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

getVasomotion(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

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
fNames=getFileNamesList(rootFolder,'*_c_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

getPulsatility(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 11 Assign vessel types and regions of interest
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
fNames=getFileNamesList(rootFolder,'*_BFI_d.mat','[A-Z]+\d+','1BP_c_BFI_d\.mat');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 11 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

