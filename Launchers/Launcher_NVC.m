% Launcher_NVC - Externally-triggered neurovascular coupling (NVC) response analysis.
%
%   Example launcher for analysing externally-triggered NVC responses.  Edit the
%   file names and settings for your project.  Run STEP 0 once per MATLAB
%   session, then run the step cells (%%) in order.
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


%% STEP 3 Segment the temporal contrast (categories + labels; also builds the cMask STEP 5 needs)
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

%RUN THE PROCESSING ROUTINE (this writes the per-segment traces STEP 8 measures the responses in)
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

%SET FILE NAMES HERE - the contrast products, which is what registration has always
%taken. There is no epoch-averaged product any more, so there is nothing else to register.
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

%% STEP 7 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%%IF DOING IT FILE BY FILE
% s.useReference=false;
% s.refFName=''; %use '' instead of " "
%fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
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

%% STEP 8 Measure the neurovascular coupling response, every stimulus repetition on its own
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder
s.libraryFolder=libraryFolder;

%THIS IS THE WHOLE NVC ANALYSIS, and it runs on the recording itself. Nothing is
%averaged before it is measured: every repetition of the stimulus is cut out of the
%recording and characterised separately, so a response that fades over the session,
%or one bad repetition, is something you can see instead of something that quietly
%moved the average.
%
%THE RESPONSE IS MEASURED, NOT FITTED. Fourteen numbers describe each repetition of
%each segment - its resting level and how noisy that was, how far flow rose during the
%stimulus and at its peak, how much extra perfusion the whole repetition delivered, and
%when 10, 50 and 90 % of that had arrived. Two more numbers say how much to trust the
%repetition for that segment, and one more says whether the recording as a whole trusts
%it, from how much of the imaged area responded.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - WHERE THE REPETITIONS ARE
%Type of stimulation start information. Use either 'manual' for a list of
%starting times to be used in the 'HH:mm:ss.SSS' format, or 'offset' for
%stimulation that starts at a fixed time from the recording start.
%Note: this is the start of the first EPOCH, not of the stimulation inside it.
s.stimStartType='offset';
s.stimOffset=0; %seconds (for the offset mode)
%s.stimStartAll{1}='09:23:31.346'; %'HH:mm:ss.SSS', one per file (for the manual mode)

s.epochsN=20;           %how many repetitions
s.epochDurationSec=30;  %how long one repetition is, seconds
%time from the start of the epoch considered to be baseline. Example [0,10] means that
%the baseline starts with the epoch (thus 0) and ends 10 seconds in.
s.epochBaselineSec=[0,10];
s.epochStimStartSec=10; %when the stimulation actually starts within the epoch
s.stimDurationSec=5;    %how long the stimulation lasts, seconds. THIS STEP USES IT.
%time from the end of the epoch when flow is expected to be back to baseline.
%Example: [-5,0] means the finale starts 5 seconds before the end of the epoch.
s.epochFinaleSec=[-5,0];

%ADJUSTED IF NECESSARY - HOW A RESPONSE IS MEASURED
s.nvcPeakGraceSec=2;   %how long after the stimulation stops the peak may still arrive,
%seconds. On the reference recording the peak is inside two seconds for 95 % of segments.
s.nvcAreaPcts=[10,50,90]; %the response is timed by when these percentages of the
%accumulated flow increase had been delivered. Whole percentages; each one is stored as
%a measurement of its own, named after its level.

%ADJUSTED IF NECESSARY - WHAT TO ANALYSE AND WHAT TO STORE
s.nvcSignals={'sData','gsData','dvsData','dvsDiameter'};
s.segNvcReturn={'levels','amplitudes','times'}; %levels (resting and end-of-repetition
%flow and their noise), amplitudes (stimulus, peak and area) and times. All of them are
%always measured - this only decides what is kept. Default = all three.
s.ppxNvcReturn=[];     %per-pixel maps, one per repetition. Leave it empty to skip them:
%this is the slowest and largest thing the step can be asked for.

%ADJUSTED IF NECESSARY - HOW MUCH TO TRUST ONE REPETITION OF ONE SEGMENT
%Every repetition of every segment gets two numbers between 0 and 1. The first is the
%overall confidence, the second is the score of the single weakest check - so a
%repetition can be discarded for one clear reason instead of for general mediocrity.
%Nothing is deleted by a low score: every repetition is measured and stored either way.
s.nvcConfThreshold=0.6;    %the overall confidence a repetition needs to be used
s.nvcConfMinThreshold=0.2; %the lowest any single check may score
s.nvcReturnScale=5;    %how far flow may still be from its own baseline at the end of a
%repetition, in noise levels, before the score starts to fall. A short rest between
%stimuli lowers this for every segment alike.
s.nvcDevRules=false;   %true also compares each repetition against the other repetitions
%of the same segment. Off by default: the checks that are always on judge a repetition
%without reference to any other, so they behave the same on four repetitions and forty.
s.nvcDevScale=3;       %how far a repetition may then sit from the others, in robust
%noise levels. With only a handful of repetitions that spread is itself noisy: at five
%repetitions about one ordinary repetition in five already scores below a half.
s.rejectFirstEpoch=1;  %never trust the first repetition

%ADJUSTED IF NECESSARY - WHICH REPETITIONS THE RECORDING TRUSTS
%A repetition is judged by how much of the imaged AREA responded, never by an average
%over the segments: only a part of the field is ever a responder, so a tenth of it
%responding coherently can still be a perfectly good repetition.
s.nvcEpochAreaFrac=0.10;  %how much of the imaged area has to respond
s.nvcEpochConnected=true; %count only the largest connected responding region, so
%scattered segments that passed the checks by chance do not add up to a response
s.nvcEditEpochs=1;     %1 = look at every repetition and change the decisions yourself

%ADJUSTED IF NECESSARY - THE REPRESENTATIVE REPETITION
s.nvcRepresentative=false; %true averages the trusted repetitions and REPLACES THE
%RECORDING WITH THAT AVERAGE - the traces, the flow images and the clock. THIS CANNOT BE
%UNDONE: the individual repetitions are not kept, the recording can never be re-cut, and
%the only record of which repetitions went in is the *_rep_nvccuts.jpg page.

%SET FILE NAMES HERE - the whole recording, IN THE ORDER s.stimStartAll is written in
fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runNVC(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 8b (OPTIONAL) Collect the report pages of STEP 8 into one PDF
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder

%STEP 8 writes four pages per recording: where the repetitions were cut and which the
%recording trusted ('_rep_nvccuts.jpg'), how confident every segment was in every
%repetition and which check was the weakest ('_rep_nvcconfidence.jpg'), the trusted
%repetitions and their average ('_rep_nvcresponse.jpg'), and the response repetition by
%repetition across the session ('_rep_nvctrials.jpg'). Change the tail to collect a
%different one.
tail='_rep_nvctrials.jpg';
D=getFileNamesList(resultsFolder,['*',tail]); %EVERY page under the results folder, subfolders included
makeReportPdf(D,fullfile(resultsFolder,'report_nvctrials.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 9 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_t_BFI_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

