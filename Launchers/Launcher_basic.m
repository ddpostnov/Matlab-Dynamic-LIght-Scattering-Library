% Can be used for getting the feeling of your data without complex pre-processing. 
% A minimal Launcher that consists of contrast calculation and BFI conversion. 
% Won’t fit the requirements for most of the research projects


%% STEP 1 Process .rls files to get the contrast
%LIBRARY PATH - add YOUR path manualy here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.contrastType='temporal'; %'temporal' or 'spatial'
s.contrastKernel=25; %typical values: 25 for 'temporal', 5 or 7 for 'spatial'
s.decimFactor=25; %decimates the contrast. Output framerate = original framerate / s.decimation
s.decimMethod='sharp'; %or  s.decimationMethod='leaking'; 'sharp' is only for temporal analysis and and s.decimation being a multiple integer of s.contrastKernel

%ADJUSTED IF NECESSARY - PERFORMANCE ADJUSTEMNTS
s.procType='gpu'; %use 'gpu' for spatial contrast type if high-end GPU is availible, 'cpu' otherwise

%ADJUSTED IF NECESSARY - INITIAL MASKING PARAMETERS
s.minK=0.001; %expected minimum contrast, unless exposure time is too long or setup is malfunctioning values below 0.001 are not expected.
s.maxK=0.99; %anything above 1 is an artifact, for most of the applications 0.4 would be expected, but can be up to 0.8 in stroke
s.minI=10; %expected minimum average intensity, defined by the amount of light and the exposure time during the recording. Can not be below 0, usually above 10 is expected
s.maxI=250; %expected maximum average intensity, define by the amount of light, the exposure time and the bit depth of the recording. Usually below 250 is expected.
s.minTrust=[0.99,0.99]; %per-pixel trust limits in relation to the portion of frames with minimum (0) or maximum (usually 255) intensity.
s.manualMask=0; %allows manual subselection of the area to mask

fNames{1}='C:\Dropbox\Work\Data\20230116_PSY01_a1_crop.rls';

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
fNames{1}='C:\Dropbox\Work\Data\20230116_PSY01_a1_crop_t_K_d.mat';

getBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

