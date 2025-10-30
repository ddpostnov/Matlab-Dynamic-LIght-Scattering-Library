% Pulsatility-enabled data
% Pulsatility analysis provides steady state pulsatility and blood flow information. Currently it is designed only for data with rich vascular features and high quality images.
% Are your recordings captured accordingly? 194 frames per second, cranial window etc?
% If you record multiple conditions – are they all in individual files?
% Use Launcher_pulsatility.m. Ensure to edit the file names and settings according to the project needs.


%% STEP 1 Process .rls files to get the internal cycle data
%LIBRARY PATH - add YOUR path manualy here:
libraryFolder = 'C:\Dropbox\Work\GitHub\DDPLab-private\Dynamic LIght Scattering Library v2.0';
addpath(genpath(libraryFolder));
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\DDPLab-private\Dynamic LIght Scattering Library v2.0';
addpath(genpath(libraryFolder));
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.trustLimitsK=[0.01,0.3]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke
s.trustLimitsI=[10,150];
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
s.minPromCoef=1/4;%1/2; % in respect to the std of the signal

%SET FILE NAMES HERE
%OPTION 1 - AUTOMATIC LOOKUP USING REGULAR EXPRESSIONS
s.reCalculate=false; %rewrites the files if existing
rootFolder = 'O:\HE_BFI-DATA\Sonam\B1_Red batch\LSCI\Female mouse 1'; %root folder for the files lookup
files      = dir(fullfile(rootFolder,'**','*BP.rls')); %<--- use regexp to define the files of interest
fNames     = fullfile({files.folder}', {files.name}');
if ~s.reCalculate
    outNames  = regexprep(fNames, '\.rls$', '_c_K_d.mat');
    hasOut    = cellfun(@(f) isfile(f), outNames);
    fNames    = fNames(~hasOut);
end

%OPTION 2 - SET PATH MANUALY, REQUIRES COMMENTING OUT OPTION 1
% fNames{1}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";
% fNames{2}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_baseline.rls";
% fNames{3}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241218_ID251_24hfollowup.rls";

getInternalCycle(s,fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Define pixel categories
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.maxK=0.4; % Maximum valid contrast - helps with initial masking
s.minK=0.0001; % Minimum valid contrast
s.regionsN=1; %Numer of regions for manual selection. 0 if using entire window.
s.lSizeN=61; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=11; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.3; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise

%ADJUSTED IF NECESSARY - SEGMENTATION ADJUSTEMNTS
s.lThinN=2; % Large vessels thinning (appears as internal edges)
s.sThinN=2; % Small vessels thinning (appears as internal edges)
s.iniSizeN=7; % Odd number equal or larger than the spatial contrast kernel

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','externalWalls','internalWalls','lumen'}; %CATEGORIES

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*c_K_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
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
files      = dir(fullfile(rootFolder,'**','*c_K_d.mat')); 
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
for idxA=24 %animals index list, can include all of the animal indexes, e.g.  [1,2,3,4,5] etc
    for idxR=1:1 %ROIs index list - corresponds to the number of rois used in STEP 3
        %files      = dir(fullfile(rootFolder,'**',sprintf('Roi%d*c_K_d.mat', idxR))); %Assuming the naming format as advised, e.g. LSCI_202304158_09WT_TN_MLH004BP and ROI splitting
        files      = dir(fullfile(rootFolder,'**',sprintf('Roi%d*LH%03d*c_K_d.mat', idxR, idxA))); %Assuming the naming format as advised, e.g. LSCI_202304158_09WT_TN_MLH004BP and ROI splitting
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
files      = dir(fullfile(rootFolder,'**','*c_K_d.mat'));
fNames     = fullfile({files.folder}', {files.name}');

getBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 Perform pulsatility analysis (strictly after conversion to BFI)
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED (OR VERIFIED) PER PROTOCOL - Waveform fitting
s.fitData="regions"; %"regions" - for all segmented data, "all" - for segmented data AND pixel by pixel fits, or "none"
%ADJUSTED IF NECESSARY - Waveform fitting
fitSettings={'Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off'};%,'Robust','Bisquare','MaxFunEvals',1200,'MaxIter',1000};
fitLimits={'Lower',[0,0,0,0,0,-pi,-pi,-pi,-pi,-pi,0],...
    'Upper',[Inf,Inf,Inf,Inf,Inf,pi,pi,pi,pi,pi,Inf],...
    'StartPoint',[0.9,0.8,0.7,0.6,0.5,0.1,0.2,0.3,0.4,0.5,100]};
s.fitOptions = fitoptions(fitSettings{:},fitLimits{:});
s.fitEquation = fittype('a1*sin(2*pi*x+b1)+a2*sin(4*pi*x+b2)+a3*sin(6*pi*x+b3)+a4*sin(8*pi*x+b4)+a5*sin(10*pi*x+b5)+c','options',s.fitOptions);
s.coefNames={'a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c','PR2'};

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*c_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

getPulsatility(s, fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 8 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%%IF DOING IT FILE BY FILE
% s.useReference=false; 
% s.refFName=''; %use '' instead of " "
%files      = dir(fullfile(rootFolder,'**','Roi*c_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
%fNames     = fullfile({files.folder}', {files.name}');
%setTypes(s,fNames);

%%IF using a reference
s.useReference=true; %Assumes PRE-registered files
%Settings tabs with a reference has to be done ROI by ROI if split ROIS are used. It also has
%to be done ANIMAL by ANIMAL if multiple animals are compared. It can be
%convienintly done in a loop as below, but REQUIRES proper file naming.
for idxA=24  %animals index list, can include all of the animal indexes, e.g.  [1,2,3,4,5] etc
    for idxR=1:1 %ROIs index list - corresponds to the number of rois used in STEP 3
        files      = dir(fullfile(rootFolder,'**','*c_BFI_d.mat')); %Assuming the naming format as advised, e.g. LSCI_202304158_09WT_TN_MLH004BP and ROI splitting
        %files      = dir(fullfile(rootFolder,'**',sprintf('Roi%d*LH%03d*c_BFI_d.mat', idxR, idxA)));  %Assuming the naming format as advised, e.g. LSCI_202304158_09WT_TN_MLH004BP and ROI splitting
        fNames     = fullfile({files.folder}', {files.name}');
        s.refFName=fNames{1};
        setTypes(s,fNames); %LAUNCHES THE PROCESSING ROUTINE
    end
end

%% STEP 9 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','Roi1*c_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE
