% Either used to compare steady-state perfusion or responses to intervention.
% Use Launcher_slowDynamics.m. Ensure to edit the file names and settings according to the project needs.


%% STEP 1 Process .rls files to get the contrast
%LIBRARY PATH - add YOUR path manualy here:
libraryFolder = 'C:\Dropbox\Work\GitHub\DDPLab-private\Dynamic LIght Scattering Library v2.0';
addpath(genpath(libraryFolder));
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\DDPLab-private\Dynamic LIght Scattering Library v2.0';
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

%SET FILE NAMES HERE
%OPTION 1 - AUTOMATIC LOOKUP USING REGULAR EXPRESSIONS
s.reCalculate=false; %rewrites the files if existing
rootFolder = 'O:\HE_BFI-DATA\Rasmus\Pilot - Vasoactive - Diamox'; %root folder for the files lookup
files      = dir(fullfile(rootFolder,'**','*.rls')); %<--- use regexp to define the files of interest
fNames     = fullfile({files.folder}', {files.name}');
if ~s.reCalculate
    outNames  = regexprep(fNames, '\.rls$', '_t_K_d.mat');
    hasOut    = cellfun(@(f) isfile(f), outNames);
    fNames    = fNames(~hasOut);
end

%OPTION 2 - SET PATH MANUALY, REQUIRES COMMENTING OUT OPTION 1
% fNames{1}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";
% fNames{2}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_baseline.rls";
% fNames{3}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241218_ID251_24hfollowup.rls";

getContrast(s,fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Define pixel categories
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.maxK=0.3; % Maximum valid contrast - helps with initial masking
s.minK=0.0001; % Minimum valid contrast
s.regionsN=0; %Numer of regions for manual selection. 0 if using entire window.
s.lSizeN=151; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=21; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.2; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise

%ADJUSTED IF NECESSARY - SEGMENTATION ADJUSTEMNTS
s.lThinN=2; % Large vessels thinning (appears as internal edges)
s.sThinN=2; % Small vessels thinning (appears as internal edges)
s.iniSizeN=7; % Odd number equal or larger than the spatial contrast kernel

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','externalWalls','internalWalls','lumen'}; %CATEGORIES

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*t_K_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

getCategories(s,fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 3 (OPTIONAL. Only use if 1 or more regions are defined in step 2) Split the regions. 
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%USES THE SAME FILE NAMES AS ABOVE as STEP 2
splitRegions(s,fNames); %LAUNCHES THE UTILITY ROUTINE

%% STEP 4 Perform segmentation
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - BASIC PARAMETERS
s.attmemptDS=true; %attempt to perform automated dynamic segmentation or not
s.sMinL=15; % Minimum length for segments
s.prchNSize=30; % Parenchymal pixels neighbourhoud.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
s.sMaxLBI=(1/5)./s.sMinL; %Max local bending (0 to pi per pixel)
s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
s.sMaxDK=0.2; %Max accepted std/mean for the initial diameter estimation
s.sMaxKK=0.3; %Max accepted std/mean for the initial contrast estimation
s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate

%ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
s.pInterpF=10; % leave as is

%SET FILE NAMES HERE - IF REGIONS SPLIT WAS PERFORMED, OTHERWISE KEEP THE
%SAME NAMES
files      = dir(fullfile(rootFolder,'**','*t_K_d.mat')); 
fNames     = fullfile({files.folder}', {files.name}');

getSegmentation(s, fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 5 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud.

%SET FILE NAMES HERE.
%Registration has to be done ROI by ROI if split ROIS are used. It also has
%to be done ANIMAL by ANIMAL if multiple animals are compared. It can be
%convienintly done in a loop as below, but REQUIRES proper file naming.
for idxA=261 %animals index list, can include all of the animal indexes, e.g.  [1,2,3,4,5] etc
    for idxR=1 %ROIs index list - corresponds to the number of rois used in STEP 3
        files      = dir(fullfile(rootFolder,'**',sprintf('*_30min_*t_K_d.mat', idxR, idxA))); %Assuming the naming format as advised, e.g. LSCI_202304158_09WT_TN_MLH004BP and ROI splitting
        fNames     = fullfile({files.folder}', {files.name}');
        registerLSCItoLSCI(s,fNames); %LAUNCHES THE UTILITY ROUTINE
    end
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
files      = dir(fullfile(rootFolder,'**','*_30min_*t_K_d.mat'));
fNames     = fullfile({files.folder}', {files.name}');

getBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%%IF DOING IT FILE BY FILE
%  s.useReference=false; 
%  s.refFName=''; %use '' instead of " "
% files      = dir(fullfile(rootFolder,'**','Roi*c_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
% fNames     = fullfile({files.folder}', {files.name}');
% setTypes(s,fNames);

% %%IF using a reference
s.useReference=true; %Assumes PRE-registered files
%Settings tabs with a reference has to be done ROI by ROI if split ROIS are used. It also has
%to be done ANIMAL by ANIMAL if multiple animals are compared. It can be
%convienintly done in a loop as below, but REQUIRES proper file naming.
for idxA=261 %animals index list, can include all of the animal indexes, e.g.  [1,2,3,4,5] etc
    for idxR=1 %ROIs index list - corresponds to the number of rois used in STEP 3
        files      = dir(fullfile(rootFolder,'**',sprintf('*_30min_*t_BFI_d.mat', idxR, idxA))); %Assuming the naming format as advised, e.g. LSCI_202304158_09WT_TN_MLH004BP and ROI splitting
        fNames     = fullfile({files.folder}', {files.name}');
        s.refFName=fNames{1};
        setTypes(s,fNames); %LAUNCHES THE PROCESSING ROUTINE
    end
end

%% STEP 8 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*c_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

