% Launcher_pulsatility_vasomotion - Combined pulsatility and vasomotion pipeline.
%
%   Example launcher for combined pulsatility and vasomotion analysis (assumes
%   the corresponding acquisition requirements were met).  Run STEP 0 once per
%   MATLAB session, then run the step cells (%%) in order.
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
 libraryFolder = 'C:\Users\postn\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';%'C:\Users\AU707705\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
 addpath(genpath(libraryFolder));
 setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
 rootFolder = 'G:\BP'; %root folder for the files lookup
parallel.gpu.enableCUDAForwardCompatibility(true) 

%% STEP 1 Process .rls files to get the temporal contrast for segmentation and vasomotion analysis
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.contrastType='temporal'; %'temporal' or 'spatial'
s.contrastKernel=25; %typical values: 25 for 'temporal', 5 or 7 for 'spatial'
s.decimFactor=25; %decimates the contrast. Output framerate = original framerate / s.decimation
s.decimMethod='leaking'; %or  s.decimationMethod='leaking'; 'sharp' is only for temporal analysis and and s.decimation being a multiple integer of s.contrastKernel

%ADJUSTED IF NECESSARY - PERFORMANCE ADJUSTEMNTS
s.procType='gpu'; %use 'gpu' for spatial contrast type if high-end GPU is availible, 'cpu' otherwise

%ADJUSTED IF NECESSARY - INITIAL MASKING PARAMETERS
s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke
s.trustLimitsI=[1,254]; %minimum (first value) and maximum (second value) of expected intensity.
s.minTrust=[0.99,0.99]; %per-pixel trust limits in relation to the portion of frames with minimum (0) or maximum (usually 255) intensity.
s.manualMask=0; %allows manual subselection of the area to mask

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*.rls'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runContrastFromRLS',false);


runContrastFromRLS(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE
%
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
tail='_rep_contrast.jpg'; %the page STEP 1 writes, one per recording
pdfFolder=fileparts(char(fNames{find(~cellfun(@isempty,fNames(:)),1)}));
D=dir(fullfile(pdfFolder,['*',tail]));
makeReportPdf(fullfile({D.folder}',{D.name}'),fullfile(pdfFolder,'report_contrast.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 2 Process .rls files to get the internal cycle data
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
fNames=getFileNamesList(rootFolder,'*.rls'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runInternalCycle',false);

%RUN THE PROCESSING ROUTINE
runInternalCycle(s,fNames(:));
%% STEP 3 Define segmentation regions (interactive ROI editor; optional - whole window if skipped)
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

%SET FILE NAMES HERE - the TEMPORAL CONTRAST files only, GROUPED (rows = animal/FOV)
%so ROIs carry within a group.  The internal-cycle "_c" files are the SAME field of
%view, so they inherit the mask instead of being drawn again (see s.fNamesCopyTo below)
fNames=getFileNamesList(rootFolder,'*_t_K_d.mat','[A-Z]+\d+'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%Copy the regions drawn on each _t onto its paired _c (same idiom as STEP 5's
%segmentation copy) - one editor per recording instead of one per file.  Comment this
%line out, and list '*_K_d.mat' above, to draw on every file separately instead.
s.fNamesCopyTo=regexprep(fNames,'_t_K_d.mat$','_c_K_d.mat');

%RUN THE PROCESSING ROUTINE (setRegions iterates the groups itself - no for-loop)
setRegions(s,fNames);

%% STEP 4 (OPTIONAL. Only use if 1 or more regions were defined in STEP 3) Split the regions.
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%SET FILE NAMES HERE - both _t and _c now carry regionsMask (STEP 3 copied it onto _c),
%so every branch of a recording is cropped with the same regions
fNames=getFileNamesList(rootFolder,'*_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','splitRegions',false);

%RUN THE PROCESSING ROUTINE
splitRegions(s,fNames(:));

%% STEP 5 Perform segmentation on the temporal contrast, copy it onto the paired pulsatility files
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke

%ADJUSTED IF NECESSARY - CATEGORIZATION ADJUSTEMNTS
s.lSizeN=141; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=9; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.1; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
s.sSizeScale=1; % scaler for small objects assignment to background or to unregognized regions
s.deSens=1; %can be used to reduce sensitivity to small objects
s.lThinN=2; % Large vessels thinning
s.imOpen=0; % Small vessels thinning
s.iEdge=2; %setting internal edges for segmented vessels
s.eEdge=2; %setting external edges for segmented vessels

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}; %CATEGORIES

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.sMinL=15; % Minimum length for segments
s.prchNSize=50; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%SET FILE NAMES HERE - segment the temporal-contrast _t files (FLAT, order-independent)
fNames=getFileNamesList(rootFolder,'Roi*_t_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
fNames=fNames(:);

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runSegmentation',false);

% Assign the SAME segmentation to the co-registered pulsatility _c files (replaces the
% old assignCategories step): each _c inherits the _t segmentation (cMask / sMap /
% sMetrics) and gets its OWN sData re-extracted from its own cube. Requires _t and _c
% to be co-registered / same FOV. Computed AFTER any removeProcessedFiles so the source
% and target lists stay aligned.
s.fNamesCopyTo=regexprep(fNames,'_t_K_d.mat$','_c_K_d.mat');

runSegmentation(s, fNames); %LAUNCHES THE PROCESSING ROUTINE


%% STEP 6 (OPTIONAL) Dynamic segmentation - per-frame vessel diameter / flow (heavy)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 4)
s.sMinL=15; % Minimum length for segments
s.prchNSize=50; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
s.sMaxLBI=(1/7)./s.sMinL; %Max local bending (0 to pi per pixel)
s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
s.sMaxKK=0.3; %Max accepted std/mean for the initial contrast estimation
s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate

%ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
s.gSizeN=3;
s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
s.pInterpF=4; % leave as is

%SET FILE NAMES HERE (after STEP 5 this pattern also matches the RoiN_ crops)
fNames=getFileNamesList(rootFolder,'Roi*_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runDynamicSegmentation',false);

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
s.prchNSize=50; % Parenchymal pixels neighbourhoud. Same as in segmentation
s.silent=true; %run registration silently and generate report images
s.forceMethod='correlation'; % 'intensity' or 'correlation' to force in silent mode
s.rotationLimit=45; % degrees; reject registrations rotating > 45 ([] = none)


%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'Roi*_K_d.mat','[A-Z]+\d+');%,'*BV_t_K_d\.mat$');

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runRegistration',false);

%OPTIONAL - backup files
m=~cellfun(@isempty,fNames); cellfun(@copyfile,fNames(m),strrep(fNames(m),'.mat','_bckp.mat'));

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    runRegistration(s,fNames(i,:)');
end

%% OPTIONAL - fix registration issues (MANUAL) - follows right after registration uses the same settings - therefore has to be launched sequentially or the settings have to be re-assigned

%Example of manual fixing of badly registered files
fNames=strrep(fNames,'.mat','_bckp.mat');
fileIndexes=[1,17,18]; % should be [1,N1,N2,N3,N4] where 1 is always present as a reference, N1, N2 etc are badly registered files that require manual registration

s.silent=false; %run registration silently and generate report images
runRegistration(s,fNames(1,[1,17,18])');

%Delete backups
fNames=getFileNamesList(rootFolder,'*bckp.mat');
cellfun(@delete,fNames);

%% STEP 8 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=false; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'Roi*_K_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runBFI',false);

runBFI(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Perform vasomotion analysis
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - BASIC PARAMETERS
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
s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'};  %which data-type signals get vasomotion analysis
%s.segVsmReturn selects which levels to store in results.vasomotion: 'bands' (scalars.VB/CB),
%'moments' (fVectors ampMean/Std/Skew + VB.ampMeanPct), 'series' (timeVectors amp/fCent/fSprd/nPeak),
%'clustering' (flare/silence scalars+spectra+maskFlare), 'reconstruction' (timeVectors.VB.rData),
%'spectrum' (spectrum.amp/.phase grid). Change only if needed.
s.segVsmReturn={'bands','moments','series','clustering','spectrum'};
s.ppxVsmReturn =[];%{'bands'};% [];% {'bands'};         %per-pixel analysis OFF ([] = off); set non-empty (e.g. {'bands'}) to turn it ON -> RESULTS.vasomotion.ppx. No analysePerPixel flag.
s.ppxSegmentAveraging=[]; %TEMPORARY scaffolding (to be removed): per-segment averaging demo, subset of {'coherent','incoherent'}; [] = off. Change only if needed

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'Roi*_t_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

% Uncomment below if you need to avoid re-processing processed files
%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runVasomotion',false);

runVasomotion(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

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
s.ppxPulsReturn={'markers'};%,'model'};  % NON-EMPTY = per-pixel marker maps ON
                        % (results.pulsatility.ppx); add 'model'/'reconstruction' to
                        % also fit every masked pixel (large full-resolution cubes).

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'Roi*_c_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','runPulsatility',false);

runPulsatility(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

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
fNames=getFileNamesList(rootFolder,'Roi*_BFI_d.mat','[A-Z]+\d+','1BP_c_BFI_d\.mat');

%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','setVesselTypes',false);

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 12 Derive & edit the vascular parent-daughter hierarchy (after vessel types)
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - HIERARCHY DERIVATION
s.autoOnly=false;            % true = derive & save without opening the GUI
s.phiWeights=[1 1 1];        % relative weight of [foot(psTimeMin) peak(psTimeMax) -psPI] in the flow potential
s.useHarmonicPhase=false;    % include the fundamental phase b1 in the flow potential if present
s.propagatePartners={'t','s'}; % after a "_c" file is derived, auto-copy the hierarchy to its _t/_s partners if they exist

%SET FILE NAMES HERE - use the pulsatility "_c" files (they carry the pulse timing)
fNames=getFileNamesList(rootFolder,'Roi*_c_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%fNames=removeProcessedFiles(fNames,'_d.mat','_s.mat','setVascularTree',false);

setVascularTree(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 13 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_BFI_d.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

