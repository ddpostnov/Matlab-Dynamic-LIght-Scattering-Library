%runLabChart  Read a wire-myograph LabChart recording and create its data product
%
%   runLabChart(s,fNames) reads every LabChart '.adicht' recording in fNames and
%   CREATES that recording's _MYO pair:
%
%       Rat3.adicht  ->  Rat3_MYO_r.mat   Rat3_MYO_s.mat
%
%   It is the ENTRY STEP of the wire-myograph pipeline and the only function that
%   creates the pair for it - the pressure myograph's twin is runMyographVideo.
%   Everything after it (the intervals and the vasomotion) opens what is here and
%   writes it back, see myographProduct.  fNames is a cell array of raw recording
%   paths; the wrapper iterates them itself, so there is no launcher for-loop.
%
%   THE WHOLE RECORDING IS READ, unlike the pressure myograph.  A LabChart file is
%   a set of sampled traces, not a video: there is nothing to detect and nothing to
%   measure, so reading it IS the step, and s.records is how a long protocol file is
%   narrowed to the runs that matter.  Everything the later steps need is in the
%   product afterwards, and the '.adicht' is not opened again.
%
%   THE CHANNELS ARE THE MEASUREMENT, so they live in the RESULTS.  A myograph
%   product has no SOURCE member: results.channel(i) carries .name .units .fs .data
%   .time, and the interval step adds the .intervals analysed on it.  What the
%   recording IS - its name, its modality, its shared rate, its length and the
%   channels it holds by name - is results.recording, and that block holds no
%   arrays.
%
%   CHANNELS KEEP THEIR OWN SAMPLING RATES.  results.channel is the truth: one
%   element per channel, each with its own .fs, .data and .time.
%   results.recording.fs - the one rate the generic consumers may address the
%   recording by - is filled ONLY when every channel shares one, which is the common
%   case; when they differ it is left EMPTY rather than resampled, because a rate
%   that was chosen per channel is a decision about the measurement and not ours to
%   undo.
%
%   TIMES ARE ABSOLUTE SECONDS from the start of the first record, never re-based,
%   so a window can always be located back in the original recording and a comment
%   sits at the sample it describes.  Between LabChart records the base has real
%   gaps - the minutes the operator was stopped - and they are kept.
%
%   THE RECORDING CAN BE DECIMATED AS IT IS READ, and this is the only step that
%   offers it.  A wire myograph sampled at several kHz for four hours is a product
%   nothing downstream needs at that density; s.decimate reduces every channel to
%   about s.decimateFS hertz by replacing each block of samples with their mean or
%   their median, and the channels, the time base and the comments all come out of
%   readLabChart already reduced together.  Nothing after this step knows or asks:
%   it simply reads a recording that was sampled less densely.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.runLabChart.
%       .records  (optional, default []) which LabChart records to read; [] reads
%                 all of them.  Times stay absolute either way.
%       .channels (optional, default {}) names of the channels to keep; {} keeps
%                 every channel that carries samples.  A name that matches nothing
%                 is reported and ignored.
%       .decimate (optional, default false) reduce the sampling rate while reading.
%       .decimateFS (optional, default 100) the rate to reduce to, Hz.  APPROXIMATE:
%                 a whole number of samples is averaged into each kept point, so
%                 10 kHz asked down to 10 Hz gives exactly 10 Hz and asked down to
%                 3 Hz gives 3.0003 Hz.  A channel already at or below it is kept
%                 as recorded, and says so.
%       .decimateFilter (optional, default 'mean') 'mean' or 'median' - how the
%                 samples of one block are combined into the point replacing them.
%     fNames   cell array of raw '.adicht' paths.  Empty cells are skipped.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   RESULTS  results.recording (.modality='WMYO', the shared
%                        rate, the length and the channel names), results.channel
%                        with one element per channel carrying its samples,
%                        results.comments and results.blocks from the recording, an
%                        empty results.intervals, results.timeCrop = [] and
%                        results.meta
%     <name>_MYO_s.mat   SETTINGS settings.runLabChart = s
%     NO REPORT IMAGE.  Like every other myograph step, this one narrates the three
%     ordinary lines per recording and writes no page (author, 2026-08-02).  What
%     was recorded is looked at in the interval editor, which opens on the channels
%     themselves with the comments marked on them.
%
%   EXAMPLE
%     s.records  = [];               % every record in the file
%     s.channels = {};               % every channel
%     s.decimate = true;             % 10 kHz is more than anything downstream reads
%     s.decimateFS = 10;
%     s.decimateFilter = 'median';
%     D = dir(fullfile(rootFolder,'*.adicht'));
%     runLabChart(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     readLabChart (Core/Read files), myographProduct (Core/Myograph) and
%     Core/Reporting.  Through readLabChart it needs the ADInstruments SDK in
%     '3rd party/', which is Windows and 64-bit MATLAB only.
%
% See also: readLabChart, setMyographIntervals, runMyographVasomotion,
%           runMyographVideo, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - What is read
% s.records=[];           % which LabChart records to read; [] = all of them
% s.channels={};          % which channels to keep by name; {} = all of them
% s.decimate=false;       % reduce the sampling rate while reading
% s.decimateFS=100;       % the rate to reduce to, Hz (approximate)
% s.decimateFilter='mean';% 'mean' or 'median' filter before decimating

function runLabChart(s,fNames)

if ~all( cellfun(@(x) isempty(x) || ~endsWith(lower(char(x)),'.mat'), fNames(:)) )
    error('runLabChart reads the RAW recording; one or more entries are .mat files.');
end
if ~isfield(s,'records'), s.records=[]; end
if ~isfield(s,'channels') || isempty(s.channels), s.channels={}; end
if ~isfield(s,'decimate'), s.decimate=false; end
if ~isfield(s,'decimateFS'), s.decimateFS=100; end
if ~isfield(s,'decimateFilter') || isempty(s.decimateFilter), s.decimateFilter='mean'; end

% Built ONCE, before the loop: the decimation is a property of the protocol, not of
% a file.  [] is how readLabChart is told to keep every sample, so the checkbox is
% the only thing that decides - a rate left in the box does nothing on its own.
decim=[];
if ~isempty(s.decimate) && logical(s.decimate(1))
    decim=struct('fs',s.decimateFS,'filter',s.decimateFilter);
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'LabChart',fNames);
codeVer=codeVersion();                     % once per call - it shells out to git

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    s.fName=char(fNames{fidx});
    reportFile(rep,fidx,s.fName);
    clearvars results settings

    [channels,comments,blocks]=readLabChart(s.fName,s.records,decim);
    channels=keepChannels(channels,s.channels);
    if isempty(channels)
        error('runLabChart:noChannels', ...
            'No channel of %s carries any samples - nothing to analyse.',s.fName);
    end

    % ---- RESULTS: the identity card, the channels themselves, and no windows ----
    [flatT,flatFS]=flatBase(channels);
    results=struct();
    results.recording=struct( ...
        'fName',s.fName, ...
        'modality','WMYO', ...
        'frameRate',flatFS, ...             % a wire recording has no frames: it is
        'fs',flatFS, ...                    %   the shared rate, empty when they differ
        'measures',{{}}, ...                % a wire myograph measures no diameter
        'nFrames',numel(flatT), ...
        'duration',spanOf(channels), ...
        'channels',{reshape({channels.name},1,[])});
    results.timeCrop=[];
    results.comments=comments;
    results.blocks=blocks;                  % the LabChart record boundaries
    results.intervals=myographProduct('intervals');
    % THE CHANNELS ARE THE MEASUREMENT, and with no source member this is where they
    % live: one element per channel, each keeping its own rate, samples and times.
    % The axis exists from the start - a wire myograph stores its windows under it
    % too (see setMyographIntervals), and a field that appears only once a step has
    % run is a field every consumer has to guard against.
    results.channel=channelAxis(channels);
    results.meta=struct('formatVersion',4,'codeVersion',codeVer, ...
        'createdTimestamp',timestamp());

    settings=struct('runLabChart',reportSettings(s));

    reportWriting(rep);
    % This is the ENTRY step, so this is where the recording's path becomes the pair's
    % path - and the only place that has to know a project may keep its results apart
    % from its recordings.  ONLY THE NAME OF THE PRODUCT MOVES: results.recording.fName
    % above still holds the raw path, which is what myographRecordingPath looks the
    % recording up by and the only pointer back to it.  With no results folder set the
    % name comes back verbatim.
    myographProduct('save',getResultsPath(s.fName,s),results,settings);
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function channels=keepChannels(channels,wanted)
%keepChannels  The channel-name filter, applied AFTER the file was read.  A
%   LabChart file often carries a dozen traces of which a protocol uses two, and
%   naming the two here is what keeps the product that size.  A name that matches
%   nothing is said out loud rather than silently dropped: it is nearly always a
%   typo, and a silent one would look like a channel that was not recorded.
if isempty(wanted) || isempty(channels), return; end
wanted=reshape(cellstr(wanted),1,[]);
have={channels.name};
miss=wanted(~ismember(lower(strtrim(wanted)),lower(strtrim(have))));
if ~isempty(miss)
    warning('runLabChart:noSuchChannel', ...
        'This recording has no channel called %s; it holds %s.', ...
        strjoin(miss,', '),strjoin(have,', '));
end
channels=channels(ismember(lower(strtrim(have)),lower(strtrim(wanted))));
end

% =====================================================================
function chs=channelAxis(channels)
%channelAxis  results.channel, built from the prototype so every element carries
%   every field the rest of the library expects - including the .intervals the
%   interval step fills in later.  The samples are copied across as they were read:
%   a channel keeps its own rate, and nothing here resamples anything.
chs=myographProduct('channels');
for i=1:1:numel(channels)
    chs(i).name=char(channels(i).name);
    chs(i).units=char(fieldOr(channels(i),'units',''));
    chs(i).fs=double(fieldOr(channels(i),'fs',NaN));
    chs(i).data=channels(i).data(:);
    chs(i).time=channels(i).time(:);
    chs(i).intervals=myographProduct('intervals');
end
end

% =====================================================================
function dur=spanOf(channels)
%spanOf  How long the recording is: the UNION over the channels, because they may
%   differ in length when they differ in rate, and the shortest of them is not the
%   recording.
lo=inf; hi=-inf;
for i=1:1:numel(channels)
    t=double(channels(i).time);
    if isempty(t), continue; end
    lo=min(lo,t(1)); hi=max(hi,t(end));
end
if isfinite(lo) && isfinite(hi), dur=hi-lo; else, dur=0; end
end

% =====================================================================
function [time,fs]=flatBase(channels)
%flatBase  The ONE time base every channel shares, or nothing.  It exists so a
%   consumer that knows nothing about LabChart can address a wire myograph by
%   sample index, and it is filled ONLY when that reading would be true: one rate,
%   one length, one time base.  Otherwise both come back empty and the channels are
%   the only answer - see the header.
[time,fs]=deal([]);
if isempty(channels), return; end
rates=[channels.fs];
n=arrayfun(@(c) numel(c.data),channels);
if any(~isfinite(rates)) || max(rates)-min(rates)>1e-9*max(max(rates),1), return; end
if any(n~=n(1)) || n(1)==0, return; end
t=channels(1).time(:);
for i=2:1:numel(channels)
    if max(abs(channels(i).time(:)-t))>0.5/max(rates(1),eps), return; end
end
time=t; fs=double(rates(1));
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
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
