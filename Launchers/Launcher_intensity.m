% Launcher_intensity - Wide-field fluorescence (.cxd) analysis, the whole branch.
%
%   Example launcher for intravascular-tracer fluorescence recordings.  Run STEP 0
%   once per MATLAB session, then run the step cells (%%) in order - or only the
%   cells the recording is for.  THE VESSELS MUST BE BRIGHTER THAN THE TISSUE:
%   every step here reads the picture that way, and a recording the other way round
%   produces ordinary-looking numbers that are wrong.
%
%   THREE RECORDINGS' WORTH OF PRODUCTS COME OUT OF ONE .cxd, and which one you make
%   is the protocol's answer rather than a setting:
%
%     STEP 1  the ANGIOGRAM        "_a_I"  the vessels, and their slow oscillations
%     STEP 2  the CARDIAC CYCLE    "_c_I"  one averaged heartbeat, and the wall motion
%                                          and pulsing volume measured on it
%     STEP 3  the BOLUS            "_b_I"  an injection: where the tracer arrives,
%                                          when, and which vessel feeds which
%
%   They are co-registered - one recording through one objective - so the regions are
%   drawn once and inherited, and the vessel hierarchy is worked out on the injection
%   and copied onto the other two.  Everything from STEP 4 to STEP 8 is offered on all
%   three; the measurements after that each belong to one of them and say so.
%
%   WHERE THE RESULTS GO.  STEP 0 names two folders: rootFolder, where the recordings
%   are, and resultsFolder, where everything this pipeline writes goes - the _d/_r/_s
%   triplets, the report pages and the PDFs - under the same subfolders the recording
%   had.  They start out as one folder, which writes beside the recordings exactly as
%   before; point resultsFolder elsewhere to leave the raw data untouched.  Only the
%   steps that read the recording have to be told; every step after them is handed a
%   product that is already there.
%
%   ONE MICROSCOPE, ONE PIXEL SIZE, ASKED ONCE.  STEP 0 also sets how many micrometres
%   across a pixel is, and the four cells that measure a real distance take it from
%   there.  There is deliberately no default for it anywhere in the library: every
%   length in those steps scales with it, so a number carried over from another
%   microscope is a silent wrong answer.
%
%   READING THE RECORDING IS THE EXPENSIVE PART.  Opening a 60 GB .cxd takes about
%   seventeen minutes before the first frame arrives, so STEPS 1 to 3 each read it
%   exactly once and stream it a frame at a time.  Nothing after STEP 3 opens the
%   recording again.
%
% See also: runIntensity, runIntensityInternalCycle, runIntensityBolus, setRegions,
%           runBackgroundRemoval, runSegmentation, runDynamicSegmentation,
%           runTopologyAnalysis, runMotionEnhancement, runVasomotion,
%           runIntensityPulsatility, runCTTH, setVesselTypes, setVascularTree
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
pixelSize     = 2.5;                    %micrometres per pixel ON THIS MICROSCOPE. Measured, never
                                        %guessed - there is no default for it anywhere in the
                                        %library, because every length in the steps that use it
                                        %scales with it and a number from another objective is a
                                        %wrong answer with nothing looking wrong. STEP 3 records
                                        %it on the product; STEPS 5 and 9 REFUSE without it;
                                        %STEP 8 reports in pixels instead if it is left empty



%% STEP 1 Read the recording into an angiogram (the vessels, and a trace per frame)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;
s.rootFolder=rootFolder; s.resultsFolder=resultsFolder; %only the steps that read the recording need these - every step after them is already in the results tree

%ADJUSTED (OR VERIFIED) PER PROTOCOL
s.decimFactor=1; %frames averaged into each saved frame. Output framerate = original / s.decimFactor.
                 %Raise it to fit a long recording: the vessels do not move between neighbouring
                 %frames, and the picture does not change at all while the count divides evenly
s.dataTypeOut='single'; %class the frames are stored as; '' keeps the recording's own
s.saveSource=true; %false keeps ONLY the picture and the traces, at one frame of memory however
                   %long the recording is - which is how a recording too long to hold gives an
                   %angiogram at all. It is a DIFFERENT product, not a smaller one: everything
                   %that reads frames (STEPS 7, 9, 10, and cleaning the recording itself in
                   %STEP 5) has nothing to read afterwards, and reading the recording again is
                   %the only way back

%SET FILE NAMES HERE - narrow the pattern to the recordings this step is for
fNames=getFileNamesList(rootFolder,'*.cxd'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runIntensity(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Read the recording into one averaged heartbeat
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;
s.rootFolder=rootFolder; s.resultsFolder=resultsFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - PULSE DETECTION
s.maskPctI=85; %percentile of the picture above which a pixel is averaged to detect the pulse.
               %85 keeps the brightest 15 %, which is where the dye is
s.maxFrqIni=20; % initial max frequency of the activity of interest, Hz
s.minFrqIni=1; % initial min frequency of the activity of interest, Hz
s.rangeFrq=0.3; % relative frq range around the central frequency

%ADJUSTED IF NECESSARY - EXCLUSION CRITERIA
s.excludeFirstNCycles=0; %reject given number of cycles
s.coeffsSTD=[3,2,2,2,2,3,3,2,2]; %pulses rejection coefficients relative to the feature standard deviation
s.coeffsRel=[0.5,0.1]; %pulses rejection coefficients relative to the feature value
s.coeffsAbs=2; %pulses rejection coefficients relative to the absolute feature value
s.gateDetrend=3; %the dye fades over a recording; that fading is removed as a polynomial of this
                 %order before the heart rate is looked for
s.minPromCoef=1/4; %in respect to the std of the band-passed rhythm

%ADJUSTED IF NECESSARY - CYCLE CALCULATION
s.nPhaseBins=25; %points the averaged beat is described by. About ten of them are measured at
                 %these frame rates and the rest place those measurements on one common axis
s.nControls=3; %AVERAGED BEATS BUILT FROM THE SAME HEARTBEATS WITH THEIR TIMING SCRAMBLED, stored
               %beside the real one. They are what says how much of the averaged beat is the heart
               %and how much is the recording, and STEP 9 cannot run without them - a recording
               %made with 0 here is refused there by name. Each one adds as much again to the file
s.controlSeed=20260807; %fixes which scrambling is used, so re-running gives the same comparison
s.decimationSpace=4; %pixels grouped together WHILE THE HEARTBEAT IS FOUND. It costs memory only -
                     %the averaged beat is always at full resolution

%SET FILE NAMES HERE - narrow the pattern to the recordings this step is for
fNames=getFileNamesList(rootFolder,'*.cxd'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runIntensityInternalCycle(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 3 Split an injection recording into a full-speed bolus span and a picture
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;
s.rootFolder=rootFolder; s.resultsFolder=resultsFolder;

%THIS STEP IS INTERACTIVE UNLESS BOTH SPANS ARE FILLED IN, and they ship empty - the injection is
%not at the same frame in two recordings, so it is marked on the recording's own trace.
s.fBolus=[];%[301,1500]; %frames kept at the full frame rate. Leave empty to mark them on
                         %the recording. KEEP 25-30 SECONDS OF THEM and start at least
                         %1.5 s BEFORE the dye arrives: a shorter span still gives the
                         %transit time but not the spread of transit times, and a span
                         %that starts after the injection is refused by STEP 12
s.fAngio=[];%[1600,2000]; %frames averaged into the picture the vessels are found on.
                         %Leave empty to use everything after the bolus span
s.pixelSize=pixelSize;   %recorded on the product, so the steps that work in micrometres can be
                         %told what a pixel is on this microscope

%SET FILE NAMES HERE - 'BB' in the recording code names a bolus recording
fNames=getFileNamesList(rootFolder,'*BB*.cxd'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runIntensityBolus(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 3b (OPTIONAL) Collect the report pages of one step into one PDF
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize

%A wrapper writes one <recording>_rep_<page>.jpg where that recording's results go and
%stops there - the document is assembled HERE, between the steps, so a 60-file step is
%one thing to page through instead of 60 files.  COPY THIS CELL after any step whose
%pages you want together and change the tail and the output name.  On this branch the
%tails are: STEP 1 '_rep_intensity.jpg'; STEP 2 '_rep_cardiac-detect.jpg' and
%'_rep_cardiac-average.jpg'; STEP 3 '_rep_bolus.jpg'; STEP 4 '_rep_regions.jpg';
%STEP 5 '_rep_background.jpg'; STEP 6 '_rep_categories.jpg' and '_rep_segments.jpg';
%STEP 8 '_rep_topology.jpg'; STEP 9 '_rep_wall-motion.jpg'; STEP 11
%'_rep_pulsatility.jpg'; STEP 12 '_rep_ctth-input.jpg', '_rep_ctth-markers.jpg',
%'_rep_ctth-confidence.jpg' and '_rep_ctth-maps.jpg'; STEP 13 '_rep_vesseltypes.jpg'.
tail='_rep_cardiac-average.jpg'; %the page you want collected, one per recording
D=getFileNamesList(resultsFolder,['*',tail]); %EVERY page under the results folder, subfolders included
makeReportPdf(D,fullfile(resultsFolder,'report_cardiac-average.pdf')); %ASSEMBLES THE DOCUMENT

%% STEP 4 Define segmentation regions (interactive ROI editor; optional - whole window if skipped)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%REGION SELECTION - setRegions is fully interactive: it opens an ROI editor per file
%(Add ROI / Delete ROI / Reset ROIs + polygon/rectangle/square/ellipse/circle shape
%selector; select an ROI and press Delete to remove it).  The number of regions is
%however many you draw - there is no count to set - and nothing advances until you
%press Done.  ROIs drawn on the first file of a group carry (editable) to the rest of
%the group and reset at the next group.  Draw nothing, or skip this step, to keep the
%whole window (no region mask is written).

%SET FILE NAMES HERE - GROUPED (rows = animal) so the ROIs drawn on one recording's first
%product carry to its other two.  A FLAT list would carry one animal's regions onto the
%next animal's field of view, which is why the grouping id is not optional here.
fNames=getFileNamesList(resultsFolder,'*_I_r.mat','[A-Z]+\d+'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE (setRegions iterates the groups itself - no for-loop)
setRegions(s,fNames);

%TO DRAW ONCE PER ANIMAL INSTEAD OF ONCE PER PRODUCT: open the editor on the angiogram
%alone and let the other two branches inherit the mask verbatim.
% fNames=getFileNamesList(resultsFolder,'*_a_I_r.mat');
% s.fNamesCopyTo=cellfun(@(f) {regexprep(f,'_a_I_r\.mat$','_c_I_r.mat'), ...
%                              regexprep(f,'_a_I_r\.mat$','_b_I_r.mat')}, fNames, 'UniformOutput',false);
% setRegions(s,fNames(:));

%% STEP 5 (OPTIONAL) Take the glow off the picture the vessels are found on
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%WHAT THIS IS FOR: a clearer picture for the steps that work on vessel shape - the
%segmentation, the vascular density and the wall dynamics.  The haze it removes is the
%dye's own light, spread sideways by the tissue.  It cleans the PICTURE and stores the
%glow beside it, so it is exactly reversible and running it twice gives the second
%radius rather than both.

%ADJUSTED (OR VERIFIED) PER PROTOCOL
s.pixelSize=pixelSize; %every length below is in micrometres; there is no default
s.radiusUm=200;    % the width of the glow around the vessels.  Measured: 200 um leaves
                   % 0.2 % of the halo where 377.5 um leaves 10.9 %.  Set it above the
                   % widest vessel you want to keep
s.medianUm=37.5;   % smooths speckle before the glow is measured; keep it below the
                   % narrowest vessel you care about
s.applyToSource=false; % clean every FRAME too, not just the picture.  Needed for a
                   % per-frame diameter (STEP 7) and for wall motion (STEP 9), and it
                   % costs a transit time (STEP 12) about two seconds of arrival - so
                   % leave it off if you want timings from this recording.  Which was
                   % chosen is recorded on the product either way, and a second pass
                   % over the frames is refused

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runBackgroundRemoval(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 6 Perform segmentation (categories + labels + per-segment traces; fully automatic)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - CATEGORIZATION ADJUSTEMNTS (fluorescence path: no trust limits)
s.lSizeN=121; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=5; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.2; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
s.sSizeScale=1; % scaler for small objects assignment to background or to unregognized regions
s.deSens=1; %can be used to reduce sensitivity to small objects
s.lThinN=2; % Large vessels thinning
s.imOpen=0; % Small vessels thinning
s.iEdge=1; %setting internal edges for segmented vessels
s.eEdge=1; %setting external edges for segmented vessels
s.diffusionSchedule='single'; %'single' on this branch and 'multiscale' on the speckle one, which
                              %is what each was already doing when the file name chose it. The
                              %nine extra coarse-to-fine passes were written for speckle
                              %contrast's noise; a fluorescence picture does not have it. LEAVE
                              %IT SET - an absent field means 'multiscale'

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}; %CATEGORIES

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
%sMinL, prchNSize, correctNodes, simR AND difR ARE REPEATED IN STEPS 7, 8 AND 9, and all
%four steps grow the same per-segment labels from them - so numbers found with different
%values cannot be compared with each other.  STEP 8 says so on every page when it meets a
%batch that was not processed the same way.
s.sStat='median'; % Statistics used for calculation of traces per segment. 'median' or 'mean'. Median is used by default.
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction (e.g. when a vessel is suspected to be crossed by another vessel, rather than to branch)
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%SET FILE NAMES HERE - FLAT (order-independent; grouping was setRegions' job in STEP 4)
fNames=getFileNamesList(resultsFolder,'*_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 (OPTIONAL) Dynamic segmentation - per-frame vessel diameter (heavy)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%WHAT IT IS FOR ON THIS BRANCH: the slow diameter changes STEP 10 reads, which move whole
%pixels over tens of seconds.  It is NOT the instrument for a heartbeat - the trace it
%produces cannot change by less than a quarter of a pixel and a cardiac width change here
%is smaller than that, which is why STEP 9 exists and STEP 11 writes no diameter column.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 6)
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; % Enable/disable branching correction
s.simR=0.3; % minimal similarity ratio between branches to be considered the same vessel
s.difR=0.4; % minimal difference ratio to be considered different vessels

%ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
s.sMaxLBI=(1/5)./s.sMinL; %Max local bending (0 to pi per pixel)
s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
s.sMaxKK=0.3; %Max accepted std/mean for the initial estimation
s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate

%ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
s.gSizeN=3;
s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
s.pInterpF=4; % leave as is

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runDynamicSegmentation(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 8 Measure the vascular density of the field
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%WHAT THESE NUMBERS ARE NOT.  Junction density counts where centre-lines meet, and two
%vessels crossing at different depths look the same as one vessel branching - so compare
%fields, not animals.  The tissue-to-vessel distance counts vessels below the surface as
%if they were in the plane and does not resolve the smallest ones at all, so it is not an
%oxygen supply distance.  Length density falls short where the vasculature is densest,
%because overlapping vessels merge into one line.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 6)
s.sMinL=15; % Minimum length for segments
s.prchNSize=90; % Parenchymal pixels neighbourhoud.
s.correctNodes=true; s.simR=0.3; s.difR=0.4;  % branch-node correction

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CALIBRATION AND REPORTING
s.pixelSize=pixelSize; % leave EMPTY for an uncalibrated recording - the densities and
                       % distances are then reported in pixels, and the page says so
s.evdBeyond=25;     % the share of tissue further than this from the nearest vessel is
                    % reported.  25 um is provisional: at 50 um the fraction is 5e-5 on
                    % real fields, i.e. identically zero
s.calibreEdges=0:5:100; % edges of the vessel-width histogram, um
s.evdEdges=0:2:60;      % edges of the tissue-distance histogram, um
s.keepMaps=false;   % also store the picture of the distance to the nearest vessel

%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runTopologyAnalysis(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Measure how far each vessel's walls moved over the heartbeat
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%EVERY VESSEL IS COMPARED ONLY WITH ITSELF - the same measurement on the same vessel with
%the heartbeat's timing scrambled away, which is what STEP 2 stored beside the averaged
%beat.  A vessel that does not clear that reports nothing, not a small number, and on the
%reference recording that is most of them.  How far the whole vessel SLID is reported
%beside how far its walls moved, because the two are the same size here and a width change
%quoted without it can be the preparation moving.

%ADJUSTED (OR VERIFIED) PER RIG - THERE IS NO DEFAULT
s.pixelSize=pixelSize; %micrometres per pixel

%ADJUSTED IF NECESSARY - WHAT IS ALLOWED TO CARRY A NUMBER
s.minDiameterUm=37; %vessels narrower than this report nothing. 37 um is where the smallest
                    %vessel above its own comparison was found on the reference recording -
                    %a limit of that microscope and that preparation, and a finer pixel does
                    %not lower it
s.confMin=2; %a vessel must double its own scrambled-timing floor
s.minCohere=0; %refuse a vessel whose measured places disagree more than this
s.minCuts=3; %a vessel measured across fewer places than this is one place wearing a vessel's name

%ADJUSTED IF NECESSARY - WHERE THE MEASUREMENTS GO
s.cutSpacing=3; %arc length between measured places along a vessel, px
s.endMarginRadii=1.5; %no measurement within this many vessel radii of a vessel's end
s.spanRadii=3; %each measurement reaches this many radii either side of the centre line...
s.spanPadUm=30; %...plus this much tissue beyond it, micrometres
s.smoothRadii=1; %how much the vessel's path is smoothed, in radii
s.smoothMin=3; %...floored here, px
s.interpF=4; %points per pixel along the resampled centre line
s.maxCLR=Inf; %refuse a vessel more curved than this. Inf lets every shape through
s.cutStepPx=0.25; %sample spacing across a vessel, px
s.widthTol=[0.5,2]; %fitted diameter over the outlines', accepted range
s.edgeTol=1.0; %how blurred a wall may be, as a share of the vessel's own width

%ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH THE SEGMENTATION STEP 6)
s.sMinL=15; s.prchNSize=90; s.correctNodes=true; s.simR=0.3; s.difR=0.4;

%CONTINUOUS RECORDINGS ONLY - THERE IS NO DEFAULT
s.bandHz=[]; %the heartbeat's band, Hz. It is not needed on the averaged beat below, and it
             %is refused by name on an angiogram or a bolus - see the note at the end of
             %this cell before running it on one

%OPTIONAL - THE ILLUSTRATION, WHICH IS NEVER THE DELIVERABLE
s.writeVideo=false; %on: the magnified beat AND the magnified scrambled-timing comparison,
                    %written as a pair because a magnified movie looks convincing whether or
                    %not anything moved. A few hundred megabytes each
s.videoFolder=''; %empty writes them in a folder BESIDE the results folder, never inside it
s.alpha=[53.1,106.2,212.5]; %how much the movie exaggerates the movement, per detail scale.
                            %One setting suits one vessel width, and the width it is right
                            %for is written on every frame
s.levels=3:5; %which detail scales the movie exaggerates

%SET FILE NAMES HERE - THE AVERAGED BEAT, which is the only product of this branch on which
%this measurement has ever cleared its own comparison
fNames=getFileNamesList(resultsFolder,'*_c_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runMotionEnhancement(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%ON AN ANGIOGRAM OR A BOLUS the step measures the same vessels inside a band you name,
%against an empty band of the same width - and on the reference recording that came out at
%1.0 times its own comparison, i.e. nothing, while looking entirely convincing.  Run it
%that way only to check that, and take the band from the rate STEP 2 reports for the same
%animal (results.rateMedian, results.band).
% s.bandHz=[9.5,12.5];
% fNames=getFileNamesList(resultsFolder,'*_a_I_r.mat');
% runMotionEnhancement(s,fNames(:));

%% STEP 10 Perform vasomotion analysis (the angiogram - the only product with a clock)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - Frequency bands
s.vFR=[0.05,0.25];      % vasomotion band VB [lo hi], Hz
s.cFR=[0.4,0.6];        % comparison band CB [lo hi], Hz
s.wFR=[0.01,1];         % wavelet frequency limits [lo hi], Hz
s.wVPO=10;              % wavelet voices per octave

%ADJUSTED IF NECESSARY - Normalisation, percentiles and spectrum decimation
s.normalisation='median'; % global 'mean'/'median' or moving 'mmean'/'mmedian'
s.normsize=101;         % window for 'mmean'/'mmedian' (>=3; inf or 0 -> global)
s.tgtFS=1;              % target sampling frequency for the kept spectrum, Hz
s.pcts=0:10:100;        % percentile edges for percentile-resolved spectra

%ADJUSTED IF NECESSARY - Flare/silence segmentation (Otsu + elbow) and peak count
s.otsuMaxN=5;           % max number of Otsu thresholds to test
s.otsuElbow=0.05;       % elbow threshold for choosing the number of thresholds
s.nPeakProm=0.10;       % VB peak-count prominence, fraction of per-time band max-min range

%ADJUSTED IF NECESSARY - Which signals / outputs to retain
s.vsmSignals={'sData','dvsData','dvsDiameter','gsData'}; %the last two exist only where STEP 7
                        % has run; a signal a product does not carry is skipped
s.ppxVsmReturn=[];      % per-pixel analysis OFF ([] = off); non-empty (e.g. {'bands'}) turns it ON
s.ppxSegmentAveraging=[]; % TEMPORARY scaffolding: {'coherent','incoherent'} or [] (off)
s.segVsmReturn={'bands','moments','series','clustering','reconstruction'};

%SET FILE NAMES HERE - THE ANGIOGRAM ONLY.  The averaged beat is a tenth of a second on a
%phase axis, so there is no clock for a wavelet to run on, and a bolus is one to seven
%cycles of the slowest band edge - both are refused by name.
fNames=getFileNamesList(resultsFolder,'*_a_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runVasomotion(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 11 Measure how much the dye in each vessel rises and falls over one heartbeat
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%THIS IS A VOLUME, NOT A SPEED AND NOT A WIDTH.  Fluorescence from a tracer in the plasma
%is proportional to how much labelled plasma is in the light path, so the columns
%(pvPI, pvRI, pvMean, pvTimeMin, pvTimeMax, pvSymRatio) are a pulsatile plasma volume.
%How much each vessel's WIDTH changed over the same beat is STEP 9's answer, in its own
%columns - the two are different instruments and their numbers do not belong in one column.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - Harmonic model
s.nHarm=5;              % harmonics in the optional model. About ten points of a
                        % fluorescence beat are actually measured, so five is the most a
                        % model can be told apart from the beat itself

%ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
s.segPulsReturn={'markers'};  % the plain measurements, which need no fit and are the
                        % default. Add 'model' and/or 'reconstruction' to fit every beat -
                        % on this kind of recording that follows the beat rather than
                        % describing it

%ADJUSTED IF NECESSARY - Per-pixel maps (GATE; [] = per-pixel analysis off)
s.ppxPulsReturn={'markers'};  % NON-EMPTY = the per-pixel marker maps ON. 'markers' is the
                        % only level this branch offers per pixel

%SET FILE NAMES HERE - THE AVERAGED BEAT.  This is NOT runPulsatility: that wrapper writes
%ps* columns meaning a pulsatile FLOW, and it refuses a fluorescence product by name.
fNames=getFileNamesList(resultsFolder,'*_c_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

runIntensityPulsatility(s, fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 12 Measure the bolus transit: arrival, mean transit time and the spread of it
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
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

%% STEP 13 Assign vessel types (arteries and veins), from the tracer's arrival
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%THE GUESS IS DRIVEN BY WHEN THE TRACER REACHED EACH VESSEL, which STEP 12 measures - so
%the reference is the INJECTION product and the types painted on it are carried onto the
%animal's other products.  Run this after STEP 12.  On a recording with no injection the
%editor still opens and the types can be painted by hand, and the step says so rather than
%reporting every vessel as undecided.
s.useReference=true; %ASSUMES PRE-REGISTERED FILES - one recording's three products always are

%SET FILE NAMES HERE - GROUPED (rows = animal), with the injection product pinned FIRST of
%each row, which is what makes it the reference the row inherits from
fNames=getFileNamesList(resultsFolder,'*_I_r.mat','[A-Z]+\d+','_b_I_r\.mat$'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

%RUN THE PROCESSING ROUTINE (one animal at a time - the editor is per reference)
for i=1:1:size(fNames,1)
    s.refFName=fNames{i,1};
    setVesselTypes(s,fNames(i,:)');
end

%% STEP 14 Work out which vessel feeds which (the flow hierarchy)
close all
clearvars -except fNames libraryFolder rootFolder resultsFolder pixelSize
s.libraryFolder=libraryFolder;

%IT IS DERIVED ON THE INJECTION RECORDING AND COPIED ONTO THE OTHER TWO.  The direction
%comes from when the tracer reached each vessel, which spreads over seconds across the
%field; a heartbeat's phase is confined inside a tenth of a second and wraps, so two
%vessels a beat apart in true arrival cannot be told apart on it.  The step STOPS rather
%than guessing if the transit times are not there.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - HIERARCHY DERIVATION
s.autoOnly=false; % true to derive & save without opening the editor
s.connectivity=8; % 4 or 8 pixel adjacency between segments
s.minBorder=1; % minimum shared-border pixels to treat two segments as touching

%FOV BRIDGING - reconnect same-type vessels split by a crossing vessel in the
%2D projection.  Search this many pixels from the vessel tips / walls:
s.bridgeTipRadius=50; % px gap searched from the vessel ends (0 = off)
s.bridgeWallRadius=0; % px gap searched from the vessel sides (usually 0)

%LEAVE s.flowParams UNSET and each file is ordered by the columns its own product carries -
%on an injection product the tracer's arrival and its mean transit time.  Setting it by
%hand is used verbatim on every file, which is the one way to order a hierarchy by the
%heartbeat's timing instead:
% s.flowParams=struct('name',{'pvTimeMax','pvTimeMin','pvPI'},'label',{'peak','foot','PI'},'role',{'arrival','arrival','pulsatility'},'scope',{'','',''},'weight',{2,1,0.5},'enabled',{true,true,true});

%SET FILE NAMES HERE - THE INJECTION PRODUCT.  The angiogram and the averaged beat inherit
%the hierarchy afterwards, which the step does by itself.
fNames=getFileNamesList(resultsFolder,'*_b_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.

setVascularTree(s,fNames(:)); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 15 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
fNames=getFileNamesList(resultsFolder,'*_I_r.mat'); %if structured file names were used then the getFileNamesList function can be used to populate them correctly. Otherwise you can generate fNames list manually.
%INTERACTIVE ALTERNATIVE: run guiExport - the standalone export tool in front of this
%same routine (pick files / a folder / a workbench session, choose the parameters and
%the averaging, write one workbook per recording or one merged workbook for statistics)
exportToExcel(fNames(:)); %LAUNCHES THE UTILITY ROUTINE
