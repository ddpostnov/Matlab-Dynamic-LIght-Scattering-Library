% Can be used for getting the feeling of your data without complex pre-processing. 
% A minimal Launcher that consists of contrast calculation and BFI conversion. 
% Won’t fit the requirements for most of the research projects


%% STEP 1 Process .rls files to get the contrast
%LIBRARY PATH - add YOUR path manualy here:
libraryFolder = 'C:\Dropbox\Work\GitHub\DDPLab-private\Dynamic LIght Scattering Library v2.0';
addpath(genpath(libraryFolder));
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.contrastType='temporal'; %'temporal' or 'spatial'
s.contrastKernel=25; %typical values: 25 for 'temporal', 5 or 7 for 'spatial'
s.decimation=25; %decimates the contrast. Output framerate = original framerate / s.decimation

%ADJUSTED IF NECESSARY - PERFORMANCE ADJUSTEMNTS
s.procType='fastcpu'; %use 'fastgpu' for spatial contrast type if high-end GPU is availible, 'fastcpu' otherwise
s.rawBatchSize=1000; %only affects processing speed, depends on availible memory (GPU and RAM)

%ADJUSTED IF NECESSARY - INITIAL MASKING PARAMETERS
s.minK=0.001; %expected s.minImum contrast, unless exposure time is too long or setup is malfunctioning values below 0.001 are not expected.
s.maxK=0.99; %anything above 1 is an artifact, for most of the applications 0.4-0.7 is expected s.maxImum.
s.minI=10; %expected s.minImum contrast, unless exposure time is too long or setup is malfunctioning values below 0.001 are not expected.
s.maxI=255; %anything above 1 is an artifact, for most of the applications 0.4-0.7 is expected s.maxImum.
s.minTrust=[0.68,0.68,0.68]; %in relation of uncertain frames to the total number
s.manualMask=0; %allows manual subselection of the area to mask

fNames{1}='C:\Dropbox\Work\Data\Mia\LSCI_20250317_1WTTM02BN.rls';

getContrast(s,fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
fNames{1}='C:\Dropbox\Work\Data\Mia\LSCI_20250317_1WTTM02BN_t_K_d.mat';

getBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

