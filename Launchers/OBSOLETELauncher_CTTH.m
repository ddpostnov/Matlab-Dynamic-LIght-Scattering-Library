% Launcher_CTTH - Bolus-injection (CTTH) analysis from wide-field fluorescence.
%
%   Example launcher for bolus-injection analysis of wide-field fluorescence
%   data.  Run STEP 0 once per MATLAB session, then run the step cells (%%) in
%   order.
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
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
rootFolder    = 'C:\Data\mia'; %where the RECORDINGS are
resultsFolder = rootFolder;    %where the RESULTS go. Point it elsewhere to keep the raw data untouched


%% STEP 1 Split each bolus recording into a full-speed bolus span and an angiogram
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;
s.rootFolder=rootFolder; s.resultsFolder=resultsFolder; %only the step that reads the recording needs these - every step after it is already in the results tree

s.fBolus=[];%[301,1500]; %frames kept at the full frame rate. Leave empty to mark them on
                         %the recording. KEEP 25-30 SECONDS OF THEM and start at least
                         %1.5 s BEFORE the dye arrives: a shorter span still gives the
                         %transit time but not the spread of transit times, and a span
                         %that starts after the injection is refused by STEP 6
s.fAngio=[];%[1600,2000]; %frames averaged into the picture the vessels are found on.
                         %Leave empty to use everything after the bolus span
s.pixelSize=2.5;         %micrometres per pixel on this microscope. It is recorded on the
                         %product, and the steps that work in micrometres refuse without it

%SET FILE NAMES HERE - 'BB' in the recording code names a bolus recording
fNames=getFileNamesList(rootFolder,'*BB*.cxd'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runIntensityBolus(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

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
tail='_rep_bolus.jpg'; %the page STEP 1 writes, one per recording
D=getFileNamesList(resultsFolder,['*',tail]); %EVERY page under the results folder, subfolders included
makeReportPdf(D,fullfile(resultsFolder,'report_bolus.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 2 Define segmentation regions (interactive ROI editor; optional - whole window if skipped)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;

%REGION SELECTION - setRegions is fully interactive: it opens an ROI editor per file
%(Add ROI / Delete ROI / Reset ROIs + polygon/rectangle/square/ellipse/circle shape
%selector; select an ROI and press Delete to remove it).  The number of regions is
%however many you draw - there is no count to set - and nothing advances until you
%press Done.  Draw nothing, or skip this step, to keep the whole window (no region
%mask is written).

%SET FILE NAMES HERE - passed as one ROW (transpose) so a drawn ROI carries (editable)
%across all bolus files as a single group (the intensity _I path has no grouping id)
fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

setRegions(s,fNames(:)'); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 3 Perform segmentation (categories + labels + per-segment traces; fully automatic)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - CATEGORIZATION ADJUSTEMNTS (intensity _I path: no trust limits)
s.lSizeN=121; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=5; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.2; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
s.sSizeScale=1; % scaler for small objects assignment to background or to unregognized regions
s.deSens=1; %can be used to reduce sensitivity to small objects
s.lThinN=2; % Large vessels thinning
s.imOpen=0; % Small vessels thinning
s.iEdge=1; %setting internal edges for segmented vessels
s.eEdge=1; %setting external edges for segmented vessels

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}; %CATEGORIES

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%SET FILE NAMES HERE - FLAT (order-independent; grouping was setRegions' job in STEP 2)
fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 4 (OPTIONAL. Only use if 1 or more regions were defined in STEP 2) Split the regions.
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%USES THE SAME FILE NAMES AS ABOVE as STEP 2 (crops each file by its own regionsMask -> RoiN_ files)
splitRegions(s,fNames(:)); %LAUNCHES THE UTILITY ROUTINE

%% STEP 5 (OPTIONAL) Dynamic segmentation - per-frame vessel diameter / flow (heavy)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 3)
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
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
fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runDynamicSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 6 Measure the bolus transit: arrival, mean transit time and the spread of it
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%EVERY TIME IS MEASURED AGAINST THE RECORDING'S OWN ARTERIAL CURVE, which the step works
%out for itself from the vessels that fill first - there is no region to draw.  That is
%what makes the delay, the mean transit time and the spread comparable between
%recordings; the plain times (T0, T10...T90, the time of the peak) are not, because the
%injection lasts several seconds and adds itself to every one of them.

%ADJUSTED IF NECESSARY - HOW A CURVE IS READ
s.ctthLevelPcts=[10,25,50,75,90]; % how far the tracer has filled, in per cent, at each
                        % time that gets reported
s.ctthPlateauSec=1;     % how much of the end of the recording is averaged to read the
                        % final level, seconds
s.ctthGuardSec=0.5;     % quiet time left between the pre-injection period and the arrival
s.ctthSlopeSec=0.1;     % how long a stretch the steepest rise and fall are measured over

%ADJUSTED IF NECESSARY - WHERE THE ARTERIAL CURVE COMES FROM
s.ctthInputAmpPct=75;   % how bright a vessel has to be to count towards it. Do not lower
                        % this far: without the brightness gate the "earliest-filling"
                        % vessels are the dimmest ones, which cross their own threshold
                        % on noise, and every transit time would be referenced to that
s.ctthInputPct=5;       % how many of the earliest-filling vessels make it up

%ADJUSTED IF NECESSARY - HOW MUCH TO TRUST A VESSEL
s.ctthConfThreshold=0.6;    % how good all the checks together have to be
s.ctthConfMinThreshold=0.2; % how good the single worst check has to be
s.ctthSettleFrac=1;     % how much the curve may still be rising at the end and still be used
s.ctthInputScale=0.5;   % how far ahead of the artery a vessel may fill, seconds
s.ctthCthTol=0.25;      % how close the spread check has to come before the recording is
                        % judged long enough for it. Below that the spread is left blank
                        % rather than reported small

%ADJUSTED IF NECESSARY - WHICH NUMBERS TO KEEP
s.segCtthReturn={'levels','amplitudes','times','transit','shape'};
s.ppxCtthReturn={};     % {} = the vessels alone, which is what a table is read from.
                        % Non-empty ALSO makes a picture of every listed marker, one
                        % value per pixel - about six minutes on a full field
s.parforCTTHPixels=true; % false measures one pixel at a time and starts no parallel pool

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runCTTH(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud.

files      = dir(fullfile(resultsFolder,'**','*t_K_r.mat')); %<---ALWAYS REFER TO "_K_r.mat" files, but you may use regexp to define specific "_K_r.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_K_r.mat$', '_c_K_r.mat');
fNames=cat(1,fNamesPLS,fNamesVSM);
runRegistration(s,fNames); %LAUNCHES THE UTILITY ROUTINE

%% STEP 8 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
files      = dir(fullfile(resultsFolder,'**','*_K_r.mat'));
fNames     = fullfile({files.folder}', {files.name}');

runBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Perform vasomotion analysis
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

s.ppxVsmReturn = [];         %per-pixel analysis OFF ([] = off); set non-empty (e.g. {'bands'}) to turn it ON -> RESULTS.vasomotion.ppx. No analysePerPixel flag.
s.ppxSegmentAveraging=[]; %TEMPORARY scaffolding (to be removed): per-segment averaging demo, subset of {'coherent','incoherent'}; [] = off. Change only if needed
s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'};  %which data-type signals get vasomotion analysis
%s.segVsmReturn selects which levels to store in results.vasomotion: 'bands' (scalars.VB/CB),
%'moments' (fVectors ampMean/Std/Skew + VB.ampMeanPct), 'series' (timeVectors amp/fCent/fSprd/nPeak),
%'clustering' (flare/silence scalars+spectra+maskFlare), 'reconstruction' (timeVectors.VB.rData),
%'spectrum' (spectrum.amp/.phase grid). Change only if needed.
s.segVsmReturn={'bands','moments','series','clustering','reconstruction'};

%SET FILE NAMES HERE
files      = dir(fullfile(resultsFolder,'**','*t_BFI_r.mat'));
fNames     = fullfile({files.folder}', {files.name}');

runVasomotion(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 10 Perform pulsatility analysis (strictly after conversion to BFI)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

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
files      = dir(fullfile(resultsFolder,'**','*c_BFI_r.mat')); %<---ALWAYS REFER TO "c_BFI_r.mat" files, but you may use regexp to define specific "c_BFI_r.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

runPulsatility(s, fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 11 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

s.libraryFolder=libraryFolder;
s.useReference=true; 
s.libraryFolder=libraryFolder;

files      = dir(fullfile(resultsFolder,'**','*t_BFI_r.mat')); %<---ALWAYS REFER TO "_K_r.mat" files, but you may use regexp to define specific "_K_r.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_BFI_r.mat$', '_c_BFI_r.mat');
fNames=cat(1,fNamesPLS,fNamesVSM);
s.refFName=fNames{1};
setVesselTypes(s,fNames);

%% STEP 11 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
files      = dir(fullfile(resultsFolder,'**','*_BFI_r.mat')); %<---ALWAYS REFER TO "_BFI_r.mat" files, but you may use regexp to define specific "_BFI_r.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

