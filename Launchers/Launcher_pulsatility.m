% Launcher_pulsatility - Cardiac pulsatility analysis pipeline.
%
%   Example launcher for pulsatility analysis; a cranial window and >=194 fps raw
%   data are advised.  Run STEP 0 once per MATLAB session, then run the step
%   cells (%%) in order.
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
rootFolder = 'C:\Dropbox\Work\Data'; %root folder for the files lookup



%% STEP 1 Process .rls files to get the internal cycle data
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.contrastKernelS=5; %contrast kernel for spatial (sLSCI) processing method
s.maxFrqIni=20; % initial max frequency of the activity of interest, Hz
s.minFrqIni=1; % initial min frequency of the activity of interest, Hz

%ADJUSTED (OR VERIFIED) PER PROTOCOL - PULSE DETECTION
s.maskLimitsK=[0.01,0.3]; %contrast range of the pixels averaged to detect the pulse
s.maskLimitsI=[5,250]; %intensity range of the pixels averaged to detect the pulse

%ADJUSTED (OR VERIFIED) PER PROTOCOL - OUTPUT MASKING
s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast
s.trustLimitsI=[1,254]; %minimum (first value) and maximum (second value) of expected intensity.
s.minTrust=[0.99,0.99]; %per-pixel trust limits in relation to the portion of frames with minimum (0) or maximum (usually 255) intensity.

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
fNames=getFileNamesList(rootFolder,'*BP.rls'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
runInternalCycle(s,fNames(:));

%% STEP 1b (OPTIONAL) Collect the report pages of STEP 1 into one PDF
close all
clearvars -except fNames libraryFolder rootFolder

%A wrapper writes one <recording>_rep_<page>.jpg beside each recording and stops
%there - the document is assembled HERE, between the steps, so a 60-file step is one
%thing to page through instead of 60 files.  COPY THIS CELL after any step whose
%pages you want together and change the tail and the output name: after setRegions
%it is '_rep_regions.jpg', after runSegmentation '_rep_categories.jpg' and
%'_rep_segments.jpg', after runRegistration '_rep_registration.jpg', after
%setVesselTypes '_rep_vesseltypes.jpg'.
%runInternalCycle also writes '_rep_cycle-average.jpg' - copy the cell for it too.
tail='_rep_cycle-detect.jpg'; %the page STEP 1 writes, one per recording
pdfFolder=fileparts(char(fNames{find(~cellfun(@isempty,fNames(:)),1)}));
D=dir(fullfile(pdfFolder,['*',tail]));
makeReportPdf(fullfile({D.folder}',{D.name}'),fullfile(pdfFolder,'report_cycle-detect.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 2 Define segmentation regions (interactive ROI editor; optional - whole window if skipped)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%REGION SELECTION - setRegions is fully interactive: it opens an ROI editor per file
%(Add ROI / Delete ROI / Reset ROIs + polygon/rectangle/square/ellipse/circle shape
%selector; select an ROI and press Delete to remove it).  The number of regions is
%however many you draw - there is no count to set - and nothing advances until you
%press Done.  ROIs drawn on the first file of a group carry (editable) to the rest of
%the group and reset at the next group.  Draw nothing, or skip this step, to keep the
%whole window (no region mask is written).

%SET FILE NAMES HERE - GROUPED (rows = animal/FOV) so ROIs can carry within a group
fNames=getFileNamesList(rootFolder,'*_c_K_d.mat','[A-Z]+\d+'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE (setRegions iterates the groups itself - no for-loop)
setRegions(s,fNames);

%% STEP 3 Perform segmentation (categories + labels + per-segment traces; fully automatic)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.trustLimitsK=[0.001,0.5]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke

%ADJUSTED IF NECESSARY - CATEGORIZATION ADJUSTEMNTS
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

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.sMinL=15; % Minimum length for segments
s.prchNSize=30; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%SET FILE NAMES HERE - FLAT (order-independent; grouping was setRegions' job in STEP 2)
fNames=getFileNamesList(rootFolder,'*_c_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
runSegmentation(s, fNames(:));

%% STEP 4 (OPTIONAL. Only use if 1 or more regions were defined in STEP 2) Split the regions.
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_c_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE (crops each file by its own regionsMask -> RoiN_ files)
splitRegions(s,fNames(:));

%% STEP 5 (OPTIONAL) Dynamic segmentation - per-frame vessel diameter / flow (heavy)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 3)
s.sMinL=15; % Minimum length for segments
s.prchNSize=30; % Parenchymal pixels neighbourhoud.
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
fNames=getFileNamesList(rootFolder,'*_c_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
runDynamicSegmentation(s, fNames(:));

%% STEP 6 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud.

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_c_K_d.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    runRegistration(s,fNames(i,:)');
end


%% STEP 7 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_c_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runBFI(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 8 Perform pulsatility analysis (strictly after conversion to BFI)
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

runPulsatility(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Assign vessel types and regions of interest
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
fNames=getFileNamesList(rootFolder,'*_c_BFI_d.mat','[A-Z]+\d+','1BP_c_BFI_d\.mat');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 10 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_c_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames(:)); %LAUNCHES THE UTILITY ROUTINE
