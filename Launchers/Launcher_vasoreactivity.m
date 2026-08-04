% Launcher_vasoreactivity - Vasoreactivity / vasomotion analysis pipeline.
%
%   Example launcher for vasoreactivity or vasomotion analysis; for vasomotion,
%   >=25 fps raw data is advised.  Run STEP 0 once per MATLAB session, then run
%   the step cells (%%) in order.
%
%   WHERE THE RESULTS GO.  STEP 0 names two folders: rootFolder, where the recordings
%   are, and resultsFolder, where everything this pipeline writes goes - the _d/_r/_s
%   triplets, the report pages and the PDFs - under the same subfolders the recording
%   had.  They start out as one folder, which writes beside the recordings exactly as
%   before; point resultsFolder elsewhere to leave the raw data untouched.  Only the
%   step that reads the recording has to be told; every step after it is handed a
%   product that is already there.
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
rootFolder    = 'C:\Dropbox\Work\Data'; %where the RECORDINGS are
resultsFolder = rootFolder;             %where the RESULTS go. Point it elsewhere to keep the raw data untouched

%% STEP 1 Process .rls files to get the contrast
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;
s.rootFolder=rootFolder; s.resultsFolder=resultsFolder; %only the step that reads the recording needs these - every step after it is already in the results tree

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

runContrastFromRLS(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 1b (OPTIONAL) Collect the report pages of STEP 1 into one PDF
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

%A wrapper writes one <recording>_rep_<page>.jpg where that recording's results go and
%stops there - the document is assembled HERE, between the steps, so a 60-file step is
%one thing to page through instead of 60 files.  COPY THIS CELL after any step whose
%pages you want together and change the tail and the output name: after setRegions
%it is '_rep_regions.jpg', after runSegmentation '_rep_categories.jpg' and
%'_rep_segments.jpg', after runRegistration '_rep_registration.jpg', after
%setVesselTypes '_rep_vesseltypes.jpg'.
tail='_rep_contrast.jpg'; %the page STEP 1 writes, one per recording
D=getFileNamesList(resultsFolder,['*',tail]); %EVERY page under the results folder, subfolders included
makeReportPdf(D,fullfile(resultsFolder,'report_contrast.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 2 Define segmentation regions (interactive ROI editor; optional - whole window if skipped)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;

%REGION SELECTION - setRegions is fully interactive: it opens an ROI editor per file
%(Add ROI / Delete ROI / Reset ROIs + polygon/rectangle/square/ellipse/circle shape
%selector; select an ROI and press Delete to remove it).  The number of regions is
%however many you draw - there is no count to set - and nothing advances until you
%press Done.  ROIs drawn on the first file of a group carry (editable) to the rest of
%the group and reset at the next group.  Draw nothing, or skip this step, to keep the
%whole window (no region mask is written).

%SET FILE NAMES HERE - GROUPED (rows = animal/FOV) so ROIs can carry within a group
fNames=getFileNamesList(resultsFolder,'*_t_K_r.mat','[A-Z]+\d+'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE (setRegions iterates the groups itself - no for-loop)
setRegions(s,fNames);

%% STEP 3 Perform segmentation (categories + labels + per-segment traces; fully automatic)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
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
s.sMinL=10; % Minimum length for segments
s.prchNSize=30; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%SET FILE NAMES HERE - FLAT (order-independent; grouping was setRegions' job in STEP 2)
fNames=getFileNamesList(resultsFolder,'*_t_K_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE
runSegmentation(s, fNames(:));

%% STEP 4 (OPTIONAL. Only use if 1 or more regions were defined in STEP 2) Split the regions.
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_t_K_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE (crops each file by its own regionsMask -> RoiN_ files)
splitRegions(s,fNames(:));

% Optional dynamic segmentation (per-frame vessel diameter / flow) is available as
% a further step - runDynamicSegmentation(s,fNames(:)) on these files (or the RoiN_
% crops above). It is heavy and not run in this launcher; see Launcher_pulsatility
% for a fully configured dynamic-segmentation cell.

%% STEP 5 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud - same as in the segmentation step

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_t_K_r.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    runRegistration(s,fNames(i,:)');
end


%% STEP 6 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_t_K_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runBFI(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 (Optional) Perform vasomotion analysis
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

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
fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runVasomotion(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 8 (Optional) Fit the response to a drug injection or a gas challenge
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;

%ADJUSTED PER PROTOCOL - WHEN THE DRUG WAS GIVEN.  One number per recording, in the
%same order as the file list below; a single number is used for all of them.  Nothing
%earlier in this pipeline records an injection time, so it is set here and nowhere else.
s.injectionSec=600; %seconds from the start of the recording
s.baselineSec=[60,540]; %the quiet stretch before the injection the response is measured against. It HAS TO END BEFORE the injection

%ADJUSTED (OR VERIFIED) PER PROTOCOL - RESPONSE MODEL
s.vrcModel='bateman'; %a rising and a falling phase on a drifting baseline. The reported amplitude is the peak change in the trace's own units
s.tgtFS=1; %Hz. The traces are averaged down to this rate before fitting - the response takes minutes, so 1 Hz is plenty and 45 minutes of raw frames would not fit in memory

%ADJUSTED IF NECESSARY - WHAT COUNTS AS A RESPONSE
s.vrcAucSec=2700; %how long after the injection the area under the response is taken over, seconds
s.vrcStealK=2; %how far below baseline, in baseline noise levels, counts as a dip
s.vrcStealFrac=0.15; %or this fraction of the response, whichever is the larger drop
s.vrcStealSec=60; %and how long it has to stay there before the recording is flagged
s.vrcWashFrac=0.1; %washout time = when the response is back to this fraction of its peak

%ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
s.segVrcReturn={'markers','model','reconstruction'}; %markers (model-free scalars), model (fitted parameters - runs the fit), reconstruction (the fitted curve). Default = all three
s.ppxVrcReturn=[]; %[] = off. NON-EMPTY fits every pixel of the field over the whole recording - expect many minutes

%ADJUSTED IF NECESSARY - Fitting
s.vrcStarts=16; %how many start points the fit is tried from. Fewer is faster and more likely to settle on a wrong answer

%ADJUSTED IF NECESSARY - Parallel execution
s.parforVrcSegments=true;
s.parforVrcPixels=true; %false fits one at a time in this MATLAB and starts no parallel pool

%SET FILE NAMES HERE - THE ORDER MATTERS: s.injectionSec is read position by position
fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runFitVasoreactivity(s,fNames(:));  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%%IF DOING IT FILE BY FILE
% s.useReference=false;
% s.refFName=''; %use '' instead of " "
%fNames=getFileNamesList(resultsFolder,'*_c_K_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
%setVesselTypes(s,fNames(:));

%%IF using a reference
s.useReference=true; %Assumes PRE-registered files
%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat','[A-Z]+\d+');

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 10 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

