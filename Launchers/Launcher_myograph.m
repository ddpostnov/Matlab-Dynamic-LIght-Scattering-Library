% Launcher_myograph - Pressure and wire myograph analysis pipeline.
%
%   Example launcher for a myograph experiment: a folder of pressure-myograph
%   recordings (.avi and anything else VideoReader opens) goes in, and diameter,
%   propagation and vasomotion come out.  The WIRE myograph variant - LabChart
%   .adicht recordings - is the commented block at the end of each step.
%
%   Run STEP 0 once per MATLAB session, then run the step cells (%%) in order.
%
%   WHAT A MYOGRAPH RUN LEAVES BEHIND.  One triplet per recording -
%   <name>_MYO_d.mat (SOURCE), _MYO_r.mat (RESULTS), _MYO_s.mat (SETTINGS) - and
%   three lines of text per step in the command window.  It writes NO report images:
%   the figures come from guiExplore and the numbers from exportToExcel, both at the
%   end of this file.  A user coming from the speckle pipeline should not go looking
%   for pages that were never written.
%
%   STEPS 2, 3 AND 5 OPEN A WINDOW and wait for you, one recording at a time.  All
%   three are optional and they are the only interactive steps here; everything else
%   is headless and can be left running.  The same analysis is available with no code
%   at all through guiWorkbench, which is the ordinary way in - this file is for
%   scripting it.
%
%   THE ANALYSED SPAN IS A PROPERTY OF THE RECORDING, not of s.  A crop (step 2) or a
%   set of pre-set windows (step 3) is stored in that recording's results, and the
%   diameter step reads it from there.  So the same s can be applied to a whole
%   folder whose recordings were each cropped differently.
%
% See also: guiWorkbench, guiExplore, guiExport, exportToExcel, runMyographVideo,
%           runMyographDiameter, runMyographPropagation, runMyographVasomotion,
%           setMyographCrop, setMyographPresetIntervals, setMyographIntervals,
%           runLabChart
%
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
setLibraryPath(libraryFolder); %keeps .claude/ tooling copies OFF the path - they SHADOW the library
rootFolder = 'C:\Dropbox\Work\Data'; %root folder for the files lookup



%% STEP 1 Read the recordings and create one _MYO triplet per recording
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CALIBRATION AND GEOMETRY
s.pixelSize=0.62; %um per px. Leave empty or 0 for an uncalibrated recording - results are then reported in px
s.rowRange=[1,Inf]; %[lo hi] image rows the vessel occupies. Rows outside it are never measured

%SET FILE NAMES HERE
fNamesRaw=getFileNamesList(rootFolder,'*.avi'); %the RAW recordings

%RUN THE PROCESSING ROUTINE
runMyographVideo(s,fNamesRaw(:));

%from here on every step takes the PRODUCT, not the recording
fNames=getFileNamesList(rootFolder,'*_MYO_d.mat');

%WIRE MYOGRAph VARIANT - a LabChart recording instead of a video:
% s.records=[]; %which LabChart records to read; [] reads all of them
% s.channels={}; %names of the channels to keep; {} keeps every channel with samples
% fNamesRaw=getFileNamesList(rootFolder,'*.adicht');
% runLabChart(s,fNamesRaw(:));
% fNames=getFileNamesList(rootFolder,'*_MYO_d.mat');


%% STEP 2 (OPTIONAL, INTERACTIVE) Crop the recording in time
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%One window per recording, showing a brightness profile of the whole file. Drag a
%band over the part worth analysing and press Done. NOTHING IS MEASURED YET, so the
%point of the window is the video preview beside the curve.
%ALTERNATIVE TO STEP 3 - a recording has one set of windows, chosen once. Use this
%when you want ONE span and will cut it up later; use step 3 when you already know
%the windows (baseline / drug / washout) and want them measured as themselves.

s.profileSamples=1200; %points in the brightness profile. More is finer and slower to open

%RUN THE PROCESSING ROUTINE
setMyographCrop(s,fNames(:),fNamesRaw(:));


%% STEP 3 (OPTIONAL, INTERACTIVE) Pre-set the analysis windows
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%Same window as step 2, but every band you draw becomes an analysed window with its
%own name. The diameter step then measures INSIDE them and nowhere else, so ten
%minutes of a two-hour recording costs ten minutes of memory.
%MUTUALLY EXCLUSIVE WITH STEPS 2 AND 5 - the windows are chosen once.

s.profileSamples=1200;

%RUN THE PROCESSING ROUTINE
setMyographPresetIntervals(s,fNames(:),fNamesRaw(:));


%% STEP 4 Measure the diameter
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%MANDATORY for a pressure myograph, and the step everything else reads. It measures
%all THREE diameters - outer, wall-centre and luminal - within whichever span steps
%2 or 3 left behind, or over the whole file when neither ran.
%NOT FOR A WIRE MYOGRAPH, which records channels and measures no diameter - skip
%straight to step 6.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - WALL DETECTION
s.rowRange=[1,Inf]; %[lo hi] image rows the vessel occupies
s.edgeMode='mid'; %which measure is the DEFAULT analysed one: 'outer', 'mid' (wall centre) or 'inner'. All three are always measured
s.wallContrast=0.05; %minimum contrast for a row to count as a wall
s.smoothSigma=1.2; %Gaussian pre-smoothing, px
s.dustRadius=8; %size of the dark specks removed before detection, px

%ADJUSTED IF NECESSARY - WHAT THE DIAMETER IS ALLOWED TO DO
s.tSmoothHz=1; %how fast the diameter is allowed to change, Hz
s.minWallGap=3; %smallest diameter that can be real, px
s.tSpan=25; %temporal outlier-rejection window, frames
s.ySpan=31; %along-vessel outlier-rejection window, rows
s.outlierK=3; %outlier rejection threshold, in standard deviations

%RUN THE PROCESSING ROUTINE
runMyographDiameter(s,fNames(:),fNamesRaw(:));


%% STEP 5 (OPTIONAL, INTERACTIVE) Define the analysis windows on the measured diameter
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%One window per recording, on the diameter trace this time, with the Y-vs-time map
%below it and a video preview with the detected walls overlaid. This is the ordinary
%way to choose windows: you can see what the vessel did before deciding where
%baseline ended.
%CANNOT BE COMBINED WITH STEP 3. Moving a window drops the propagation and the
%vasomotion computed for the old one, so run steps 6 and 7 after this, not before.
%A WIRE MYOGRAPH USES THE SAME WINDOW - with its channels and the LabChart operator
%comments in place of the diameter and the video - and adds a CHANNEL column to the
%table: leave a window on 'all channels' or give it one chamber of its own. Its
%results are then keyed by channel: results.channel(i).intervals(k).

s.edgeMode='mid'; %which diameter the trace shows: 'outer', 'mid' (wall centre) or 'inner'
s.drawPoints=2000; %points drawn per channel while choosing windows. Drawing only - lower it if the window feels slow, it changes nothing that comes out

%RUN THE PROCESSING ROUTINE
setMyographIntervals(s,fNames(:),fNamesRaw(:));


%% STEP 6 Propagation along the vessel
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%How fast, and in which direction, the oscillation travels along the vessel -
%estimated from the lag at maximum correlation between locations, reported with R2,
%a p-value against a row-shuffled null, and a confidence sentence.
%NOT FOR A WIRE MYOGRAPH - it has one measurement point and nothing to propagate along.
%Increasing row index = downward in the image.
%INDEPENDENT OF STEP 7 - both read the diameter and write their own branch, so the
%order is the registry's rather than a dependency.

%ADJUSTED (OR VERIFIED) PER PROTOCOL
s.diameterMeasures={'mid'}; %which diameters to analyse
s.vFR=[0.05,0.25]; %vasomotion band - it sizes the lag search around the dominant oscillation
s.detrendSec=30; %window removing slow drift before correlating, s

%ADJUSTED IF NECESSARY - HOW MUCH EVIDENCE AN ESTIMATE NEEDS
s.propMinCoh=0.3; %how well a location must correlate with the reference to be used
s.propMinRows=20; %fewest locations an estimate may rest on
s.propNShuffle=200; %surrogate repeats behind the p-value
s.propMinMeasured=0.5; %least of a location's trace that must be a real measurement
s.propMaxLagFrac=0.5; %largest lag searched, as a fraction of the dominant period

%RUN THE PROCESSING ROUTINE
runMyographPropagation(s,fNames(:));


%% STEP 7 Vasomotion
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%THE SAME ANALYSIS THE SPECKLE PIPELINE RUNS, on the diameter instead of the flow -
%one <VSM> tree per analysed window per signal, so guiExplore and exportToExcel treat
%a myograph recording exactly like a segmented one.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - WHAT IS ANALYSED
s.diameterMeasures={'mid'}; %which diameters to analyse. Ignored by a wire myograph, which analyses its channels
s.perLine=false; %false: one trace per window. true: one per measured image row - slower, and how you see a wave along the vessel

%ADJUSTED (OR VERIFIED) PER PROTOCOL - THE BANDS
s.vFR=[0.05,0.25]; %vasomotion band, Hz
s.cFR=[0.4,0.6]; %comparison/control band, Hz
s.wFR=[0.01,1]; %full analysed range, Hz
s.wVPO=10; %wavelet voices per octave

%ADJUSTED IF NECESSARY - NORMALISATION AND STORAGE
s.normalisation='median'; %how the trace is normalised before the transform
s.normsize=inf; %window for that normalisation, samples
s.tgtFS=1; %how finely the spectrum is stored, Hz
s.pcts=0:10:100; %percentile bins of the band envelope
s.segVsmReturn={'bands','moments','series','clustering','spectrum'}; %which analysis levels are computed and stored
%'bands' is every band scalar at once. To skip the expensive peak count, name the
%groups instead: {'bandsAmp','bandsSkew','bandsShape','bandsPct'} - dropping 'bandsPeak'.
%Drop 'spectrum' if you do not need the time-frequency grid; it is the largest thing stored.
s.parforMyographLines=false; %true runs the per-line loop on several workers. Only worth it with perLine=true

%RUN THE PROCESSING ROUTINE
runMyographVasomotion(s,fNames(:));


%% STEP 8 Export the numbers to Excel
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%A bare call writes every sheet the recording is about, chosen from the data:
%settings, comments (wire myograph only - the LabChart operator log), intervals,
%propagation, vasomotion, ampPct, spectra, ampPctSpectra and diameterTraces.
%One workbook per recording.

exportToExcel(fNames(:));

%ONE MERGED WORKBOOK FOR STATISTICS instead, with each row labelled by animal, type
%and experimental group:
% opts.merge=true;
% opts.outFile=fullfile(rootFolder,'myograph_summary.xlsx');
% exportToExcel(fNames(:),opts);

%OR do it interactively, which is the same exportToExcel behind a window:
% guiExport;


%% STEP 9 Look at the results
close all
clearvars -except fNames fNamesRaw libraryFolder rootFolder

%THE ONLY PLACE A MYOGRAPH FIGURE COMES FROM - no step above writes an image.
%Two views are worth opening on every recording before you believe any number:
%  Diameter map    - the diameter along the vessel over time. It says whether
%                    detection held: a lost wall shows as a band of missing values.
%  Detected walls  - one frame with the detected edges over it. It says whether the
%                    detector found the right edges at all.
%Then: diameter traces (in px or % of baseline), the individual line traces,
%diameter statistics, propagation scalars and the per-line lag behind them, and every
%vasomotion view - compared across experimental group, animal, recording type or
%analysed WINDOW, which is usually the one a myograph protocol wants on the x-axis.

guiExplore;
