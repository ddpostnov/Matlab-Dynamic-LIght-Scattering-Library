% Launcher_NVC - Externally-triggered neurovascular coupling (NVC) response analysis.
%
%   Example launcher for analysing externally-triggered NVC responses.  Edit the
%   file names and settings for your project.  Run STEP 0 once per MATLAB
%   session, then run the step cells (%%) in order.
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

getContrast(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

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
    getCategories(s,fNames(i,:)');
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
%% STEP 4 Get external cycle
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - STIM PARAMETERS
s.enablelRejectionModification=1;
%Type of stimulation start information. Use either 'manual' for a list of
%starting times to be used in the 'HH:mm:ss.SSS' format,or 'offset' for
%stimulation that starts at a fixed time from the recording start
%Note: stim start time corresponds to start of the first epoch, not the
%stimulation itself
s.stimStartType='offset';
%Set the stim offset (for offset stim start mode)
s.stimOffset=0; %seconds
%Set the list of stimulation start timestamps (for manual stim start mode)
stimStart{1}='09:23:31.346'; %'HH:mm:ss.SSS'

%define epochs (repeated stimulations) parameters
s.epochsN=20;
s.epochDurationSec=30;
%duration of single epoch, seconds
%time from start of the epoch considered to be baseline. Example [0,5]
%means that baseline starts with the epoch start (thus 0) and ends in 5
%seconds.
s.epochBaselineSec=[0,10];
s.epochStimStartSec=10; %time when stimulation actually starts

%time from the end of the epoch when flow is expected to return to baseline
%Example: [-5,0] means that finale starts 5 seconds before the end of the
%epoch and ends when the epoch ends.
s.epochFinaleSec=[-5,0];

s.maskType='cMask'; %'basic','cMask','selection';

%ADJUSTED IF NECESSARY - QUALITY CHECK
s.rejectBlCoef=1; %use Inf to disable rejection by this parameter
s.rejectEpochCoef=1; %use Inf to disable rejection by this parameter
s.rejectFinCoef=1; %use Inf to disable rejection by this parameter
s.rejectPeakCoef=1; %use Inf to disable rejection by this parameter
s.rejectBlSimCoef=1; %use Inf to disable rejection by this parameter
s.rejectSimCoef=1; %use Inf to disable rejection by this parameter
s.rejectTimeLoss=0.5; %allowed time loss due to grabbing faluere in seconds per epoch
s.rejectFirstEpoch=1; %always reject the first epoch

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

getExternalCycle(s,fNames(:));


%% STEP 5 Perform segmentation
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
fNames=getFileNamesList(rootFolder,'*_e_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
getSegmentation(s, fNames(:));

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
fNames=getFileNamesList(rootFolder,'*_e_K_d.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    registerLSCItoLSCI(s,fNames(i,:)');
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
fNames=getFileNamesList(rootFolder,'*_e_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

getBFI(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

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
fNames=getFileNamesList(rootFolder,'*_e_BFI_d.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 8 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_e_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

