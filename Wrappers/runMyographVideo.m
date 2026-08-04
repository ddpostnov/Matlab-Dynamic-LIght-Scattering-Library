%runMyographVideo  Register a pressure-myograph recording and create its data product
%
%   runMyographVideo(s,fNames) opens every pressure-myograph video in fNames, reads
%   its header, and CREATES the recording's _MYO pair:
%
%       Mouse1.avi  ->  Mouse1_MYO_r.mat   Mouse1_MYO_s.mat
%
%   It is the ENTRY STEP of the pressure-myograph pipeline and the only function
%   that creates the pair.  Everything after it - the time crop, the intervals,
%   the diameter, the propagation and the vasomotion - opens what is here and writes
%   it back (see myographProduct).  fNames is a cell array of raw recording paths;
%   the wrapper iterates them itself, so there is no launcher for-loop.
%
%   IT IS CHEAP ON PURPOSE.  Only the container header is read, so registering a
%   folder of two-hour recordings costs seconds.  What it writes is
%   results.recording, THE RECORDING'S IDENTITY CARD - what the file is, not what is
%   in it - and it holds no arrays at all.  The diameter step is what measures, and
%   what it measures lives per window in results.intervals(k).diameter.
%
%   NOTHING COPIES THE RECORDING.  A myograph video is not changed by the analysis
%   and re-reading it is fast, so the product describes it and keeps the
%   MEASUREMENT, and copies neither.  There is no _MYO_d.mat.
%
%   THERE IS NO WHOLE-RECORDING TIME BASE, on purpose.  What a later step needs to
%   address the recording is .nFrames and .frameRate, which is what setMyographCrop
%   and setMyographPresetIntervals use; a window's own times are its
%   .diameter.time, and the stretches between windows are not stored at all.
%
%   TIMES ARE ABSOLUTE SECONDS from the start of the recording, never re-based to
%   zero, so an interval can always be located back in the original footage.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.runMyographVideo.
%       .pixelSize  (optional, default []) µm per px.  Empty or 0 means the
%                   recording is uncalibrated and results are reported in px.  It is
%                   fixed ONCE, here, because it is a property of the microscope and
%                   not of an analysis; the propagation step reads it back from
%                   results.recording rather than asking for it again.
%       .rowRange   (optional, default [1 Inf]) [lo hi] image rows that hold the
%                   vessel.  Rows outside it are not measured.
%     fNames   cell array of raw pressure-myograph recording paths (any container
%              VideoReader can open - .avi, .mp4, .mov, .mkv).  Empty cells are
%              skipped.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   RESULTS  results.recording (the file name, the modality,
%                        the frame rate, the image size, the row range, the pixel
%                        size, the frame count and the duration), results.timeCrop
%                        = [], an empty results.intervals and results.meta
%     <name>_MYO_s.mat   SETTINGS settings.runMyographVideo = s
%
%   EXAMPLE
%     s.pixelSize = 0.62;          % µm per px (omit for an uncalibrated recording)
%     s.rowRange  = [40 440];      % the rows the vessel occupies
%     D = dir(fullfile(rootFolder,'*.avi'));
%     runMyographVideo(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     myographProduct (Core/Myograph), Core/Reporting, and MATLAB's base
%     VideoReader.
%
% See also: runMyographDiameter, runMyographPropagation, runMyographVasomotion,
%           myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Calibration
% s.pixelSize=[];         % µm per px; empty or 0 -> results are reported in px
% s.rowRange=[1 Inf];     % [first last] image row holding the vessel

function runMyographVideo(s,fNames)

if ~all( cellfun(@(x) isempty(x) || ~endsWith(lower(char(x)),'.mat'), fNames(:)) )
    error('runMyographVideo reads the RAW recording; one or more entries are .mat files.');
end
if ~isfield(s,'pixelSize'), s.pixelSize=[]; end
if ~isfield(s,'rowRange') || isempty(s.rowRange), s.rowRange=[1 Inf]; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Video',fNames);
codeVer=codeVersion();                     % once per call - it shells out to git

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    s.fName=char(fNames{fidx});
    reportFile(rep,fidx,s.fName);
    clearvars results settings

    % ---- the container header: the frame count, the frame rate and the size ----
    v=VideoReader(s.fName);
    nFrames=frameCount(v);
    frameRate=v.FrameRate;

    % ---- RESULTS: the recording's identity card, no crop, no intervals yet ----
    results=struct();
    results.recording=struct( ...
        'fName',s.fName, ...
        'modality','PMYO', ...
        'frameRate',frameRate, ...
        'fs',frameRate, ...
        'size',[v.Height v.Width], ...
        'rowRange',s.rowRange, ...
        'pixelSize',s.pixelSize, ...
        'measures',{{}}, ...                % filled by the diameter step
        'nFrames',nFrames, ...
        'duration',nFrames/max(frameRate,eps));
    results.timeCrop=[];
    results.intervals=myographProduct('intervals');
    results.meta=struct('formatVersion',4,'codeVersion',codeVer, ...
        'createdTimestamp',timestamp());

    settings=struct('runMyographVideo',reportSettings(s));

    reportWriting(rep);
    % This is the ENTRY step, so this is where the recording's path becomes the pair's
    % path - and the only place that has to know a project may keep its results apart
    % from its recordings.  ONLY THE NAME OF THE PRODUCT MOVES: results.recording.fName
    % above still holds the raw path, which is what myographRecordingPath looks the
    % video up by and the only pointer back to it.  With no results folder set the name
    % comes back verbatim.
    myographProduct('save',getResultsPath(s.fName,s),results,settings);
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function n=frameCount(v)
%frameCount  The number of frames, from the index when the container has one and
%   from the duration when it does not (a variable-frame-rate container may not).
n=[];
try
    n=v.NumFrames;
catch
end
if isempty(n) || ~isfinite(n), n=floor(v.Duration*v.FrameRate); end
n=max(double(n),0);
end

% =====================================================================
function v=codeVersion()
%codeVersion  Short git hash of the library, or 'unknown' - the version stamp that
%   lets a result be traced back to the code that made it.
v='unknown';
try
    here=fileparts(mfilename('fullpath'));
    [st,o]=system(sprintf('git -C "%s" rev-parse --short HEAD',here));
    if st==0 && ~isempty(strtrim(o)), v=strtrim(o); end
catch
end
end

% =====================================================================
function t=timestamp()
t=char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
end
