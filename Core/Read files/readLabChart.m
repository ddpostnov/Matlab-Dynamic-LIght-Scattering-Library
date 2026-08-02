%readLabChart - Read a LabChart .adicht recording via the ADInstruments SDK.
%
%   Reads an ADInstruments LabChart file: every channel on its own time base,
%   the operator's comments with their timing, and the record boundaries.  It is
%   the ONLY place in this library that touches the '+adi' package, so a change
%   in that SDK is a change to one file.
%
%   CHANNELS ARE NOT RESAMPLED.  LabChart routinely records channels at
%   different rates - a pressure trace at 1 kHz beside a temperature at 4 Hz is
%   ordinary - so each channel comes back with its own .fs, .data and .time and
%   the caller decides whether a common base exists.  Inventing one here would
%   either throw away samples or invent them.
%
%   TIME IS CONTINUOUS ACROSS RECORDS, in ABSOLUTE SECONDS FROM THE START OF THE
%   FIRST RECORD.  A LabChart file is a sequence of records (a new one whenever
%   the operator stopped and restarted, or changed a setting), each with its own
%   wall-clock start.  Every time this function returns - sample times, comment
%   times, block boundaries - is placed on that one base, so a comment can be
%   read against the samples it describes and a window can be located back in the
%   original recording.  THE BASE THEREFORE HAS GAPS: the minutes between a stop
%   and the next start are real and are not closed up.  Nothing may infer a
%   sampling rate from diff(time) across a record boundary - use channels(i).fs.
%
%   A CHANNEL THAT WAS NOT RECORDED IN A RECORD IS SKIPPED FOR IT, not zero-
%   filled: its samples simply do not cover that stretch, and a run of zeros
%   would be analysed as a signal that was flat rather than absent.
%
%   DECIMATION IS OPTIONAL, HAPPENS RECORD BY RECORD, AND INTERPOLATES NOTHING.
%   A wire myograph is routinely sampled at several kHz while the responses being
%   looked for last seconds, so the recording can be reduced to an approximate
%   rate HERE - the one place the full-rate samples exist, and the only place
%   where reducing them saves the memory the rest of the pipeline would carry.
%   Each kept point is the MEAN or the MEDIAN of the whole block of samples it
%   replaces; the blocks do not overlap, they cover the recording from end to end
%   (the last one is short when the samples do not divide evenly), and each point
%   is placed at the MIDDLE of the block it came from, so a decimated trace still
%   lies on top of the one it was made from.  The block is a whole number of
%   samples - round(recorded rate / requested rate) - so the rate that comes back
%   is CLOSE to the one asked for rather than equal to it, and a channel recorded
%   at or below the requested rate is kept as it was and said out loud.  EACH
%   RECORD IS REDUCED ON ITS OWN, so no kept point ever straddles the stop
%   between two records, and the comments are moved onto the samples that survive
%   (see snapComments).  blocks and meta are unchanged by it: they say when the
%   recording was made, not how densely it is now sampled.
%
%   ONE WINDOWS SETTING CAN STOP A BATCH RUN DEAD.  When the decimal symbol and the
%   list separator are the same character, ADInstruments' DLL shows a modal 'CHART
%   Locale' box the first time it opens a file - once per MATLAB session, from
%   inside a vendored binary, so nothing in MATLAB can dismiss it.  It concerns
%   LabChart's own arithmetic, not the samples returned here.  This function tests
%   for the condition BEFORE loading the DLL and says so once, with the setting that
%   removes it; see warnLocaleDialog below.
%
% Syntax:
%    [channels,comments,blocks,meta] = readLabChart(fName)
%    [channels,comments,blocks,meta] = readLabChart(fName, records)
%    [channels,comments,blocks,meta] = readLabChart(fName, records, decim)
%
% Inputs:
%    fName   - path to the .adicht file ('.adicht' is appended when missing).
%    records - (optional) which records to read ([] or omitted = all).  Times
%              stay absolute from the start of the FIRST record of the file, so
%              reading records 3-4 still reports them where they happened.
%    decim   - (optional) [] or omitted reads every sample as it was recorded.
%              A struct asks for the decimation described above:
%                 .fs      the rate to reduce every channel to, Hz.  Approximate:
%                          the nearest whole-number block size is used.
%                 .filter  'mean' | 'median' - how the samples of one block are
%                          combined into the single point that replaces them.
%
% Outputs:
%    channels - 1xN struct array, one element per channel that carries samples:
%                 .name    channel name as LabChart holds it
%                 .units   its unit, from the first record it appears in
%                 .fs      sampling rate, Hz - the rate the samples returned here
%                          actually have, i.e. AFTER decimation when it was asked for
%                 .data    [n x 1] the samples, records concatenated
%                 .time    [n x 1] their times, seconds on the base above
%                 .records which records contributed, in order
%    comments - 1xM struct array of the operator's comments:
%                 .time    seconds on the same base
%                 .text    the comment
%                 .channel the channel it was placed on, '' when file-wide
%                 .record  which record it belongs to
%    blocks   - [nRec x 2] the start and end second of each record read.
%    meta     - .fName .nRecordsInFile .records .recordStarts .recordDurations
%               .triggerDelays .sdkVersion .readTimestamp .decimation
%               (.decimation is [] when every sample was kept, otherwise the
%               requested rate, the filter and the block size used per channel)
%
% Example:
%    [ch,com,blk] = readLabChart('Rat3.adicht');
%    plot(ch(1).time, ch(1).data);  xline([com.time]);
%
%    % the same recording reduced from 10 kHz to about 10 Hz
%    ch = readLabChart('Rat3.adicht', [], struct('fs',10,'filter','median'));
%
% Dependencies: the ADInstruments SDK MATLAB package ('+adi', Jim Hokanson,
%    MIT), vendored in '3rd party/adinstruments_sdk_matlab'.  It is a MEX
%    wrapper over ADInstruments' own DLL and works on WINDOWS with 64-bit
%    MATLAB only.  Put the library on the path with setLibraryPath, which
%    resolves it; the '+adi' folder itself must NOT be added (it is a package,
%    reached through its parent).
% See also: readCXD, runLabChart, setMyographIntervals, setLibraryPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function [channels,comments,blocks,meta] = readLabChart(fName,records,decim)

if nargin<2, records=[]; end
if nargin<3, decim=[]; end
decim=checkDecimation(decim);
fName=char(fName);
if ~endsWith(lower(fName),'.adicht'), fName=[fName '.adicht']; end
if ~isfile(fName)
    error('readLabChart:noFile','No such LabChart file: %s',fName);
end
requireSDK();
warnLocaleDialog();

% remove_empty_channels is turned OFF so the channel numbering the file uses is
% the numbering that comes back - a comment names its channel by that number, and
% a filtered list would attribute comments to the wrong trace.  Channels that turn
% out to carry nothing are dropped below, by this function, after the samples have
% been counted.
f=adi.readFile(fName,'remove_empty_channels',false);

[offset,duration,trigger]=recordBase(f);
keep=recordsToRead(records,f.n_records);
blocks=[offset(keep)' offset(keep)'+duration(keep)'];

[channels,factors]=readChannels(f,keep,offset,decim);
comments=readComments(f,keep,offset);
comments=snapComments(comments,channels,factors);
meta=struct('fName',fName,'nRecordsInFile',f.n_records,'records',keep, ...
    'recordStarts',offset(keep),'recordDurations',duration(keep), ...
    'triggerDelays',trigger(keep),'sdkVersion',sdkVersion(), ...
    'readTimestamp',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')));
% Assigned rather than passed to struct(): .decimation is itself a struct, and a
% struct handed to struct() as a value is the one shape that does not mean "one
% field holding it".
meta.decimation=decimationRecord(decim,channels,factors);
end

% =====================================================================
function decim=checkDecimation(decim)
%checkDecimation  What "decimate to 10 Hz" has to say before a file is opened.
%   A rate that is missing, negative or typed as words, or a filter this reader
%   does not have, is a mistake worth stopping for: silently reading the file at
%   full rate would look like a decimation that had happened, and the difference
%   is a hundredfold in the size of every product downstream.
if isempty(decim), decim=[]; return; end
if ~isstruct(decim) || ~isscalar(decim) || ~isfield(decim,'fs') || ~isfield(decim,'filter')
    error('readLabChart:badDecimation', ...
        ['Decimation is asked for with a struct holding the rate to reduce to ' ...
         '(.fs, in Hz) and the filter to use (.filter, "mean" or "median").']);
end
if ~isnumeric(decim.fs) || ~isscalar(decim.fs) || ~isfinite(decim.fs) || decim.fs<=0
    error('readLabChart:badDecimation', ...
        'The frequency to decimate to must be a single positive number of hertz.');
end
filt=lower(strtrim(char(decim.filter)));
if ~any(strcmp(filt,{'mean','median'}))
    error('readLabChart:badDecimation', ...
        'The signal is filtered with "mean" or "median" before decimating, not "%s".',filt);
end
decim=struct('fs',double(decim.fs),'filter',filt);
end

% =====================================================================
function d=decimationRecord(decim,channels,factors)
%decimationRecord  What was actually done, for meta - [] when nothing was.
%   The requested rate is not the achieved one and the block size differs between
%   channels recorded at different rates, so the factors are reported per channel
%   beside the names they belong to.
d=[];
if isempty(decim) || isempty(channels) || ~any(factors>1), return; end
d=struct('fs',decim.fs,'filter',decim.filter,'factors',factors, ...
    'channels',{{channels.name}});
end

% =====================================================================
function requireSDK()
%requireSDK  One sentence, naming what to install and where it goes.  Not a
%   warning and not a fall-back: a LabChart file cannot be read any other way,
%   so continuing without the SDK could only produce an empty answer that looked
%   like a recording.
%   IT ASKS 'which', NOT 'exist'.  exist('adi.readFile','file') answers 0 for a
%   PACKAGE function that is perfectly reachable, so the obvious test would refuse
%   to read a file on a machine that has the SDK installed correctly.
if ~isempty(which('adi.readFile')), return; end
error('readLabChart:noSDK', ...
    ['Reading LabChart files needs the ADInstruments SDK for MATLAB ' ...
     '(https://github.com/JimHokanson/adinstruments_sdk_matlab), which belongs ' ...
     'in the library''s "3rd party" folder as "adinstruments_sdk_matlab"; it ' ...
     'runs on Windows with 64-bit MATLAB only.']);
end

% =====================================================================
function warnLocaleDialog()
%warnLocaleDialog  THE ONE WAY THIS READER CAN HANG, said out loud before it does.
%   ADInstruments' DLL checks Windows' regional settings the first time it opens a
%   file, and when the decimal symbol and the list separator are the SAME character
%   it puts up a modal 'CHART Locale' box.  That warning is about LabChart's own
%   arithmetic extension and means nothing for the samples this reader returns - but
%   it is MODAL, it comes from the DLL rather than from MATLAB, and a batch run with
%   nobody at the keyboard waits at it for ever.  It has already cost this project
%   one silent hour.
%
%   So the condition is tested HERE, in MATLAB, before the DLL is ever loaded, and
%   named once per session with the setting that removes it.  It cannot be suppressed
%   from inside MATLAB: the box lives in a vendored binary, and the cure is the
%   Windows setting itself.
persistent said
if ~isempty(said) || ~ispc, return; end
said=true;
try
    dec=winqueryreg('HKEY_CURRENT_USER','Control Panel\International','sDecimal');
    lst=winqueryreg('HKEY_CURRENT_USER','Control Panel\International','sList');
catch
    return                                  % no reading of the settings, no claim
end
if ~strcmp(char(dec),char(lst)), return; end
warning('readLabChart:localeDialog', ...
    ['Windows on this computer uses "%s" as both the decimal symbol and the list ' ...
     'separator, so the LabChart reader will show a "CHART Locale" warning box the ' ...
     'first time it opens a file. It is about LabChart''s own calculations and does ' ...
     'not affect the data read here, but it has to be clicked before reading ' ...
     'continues - so an unattended batch run will stop at it. To stop it appearing, ' ...
     'set the list separator to ";" in Windows: Settings > Time & language > ' ...
     'Language & region > Regional format > Additional settings > List separator.'], ...
     char(dec));
end

% =====================================================================
function [offset,duration,trigger]=recordBase(f)
%recordBase  WHERE EACH RECORD SITS on the file's own second base.
%   LabChart stamps every record with the wall-clock time it started, so the gap
%   between a stop and the next start is known and is kept.  When those stamps do
%   not make sense - missing, out of order, or placing a record before the
%   previous one had finished, which some converted files do - the records are
%   laid end to end instead, which loses the gaps but never overlaps two records.
n=f.n_records;
[offset,duration,trigger]=deal(zeros(1,n));
if n==0, return; end
starts=zeros(1,n);
for r=1:1:n
    duration(r)=double(f.records(r).duration);
    starts(r)=double(f.records(r).record_start);        % MATLAB datenum
    trigger(r)=double(fieldOr(f.records(r),'trigger_minus_rec_start',0));
end
offset=(starts-starts(1))*86400;                        % days -> seconds
ends=offset+duration;
ok=all(isfinite(offset)) && all(offset>=0) && all(diff(offset)>=0) && ...
   all(offset(2:end)>=ends(1:end-1)-1e-6);
if ~ok
    offset=[0 cumsum(duration(1:end-1))];
end
end

% =====================================================================
function keep=recordsToRead(records,n)
%recordsToRead  Which records were asked for, checked against what is there.
if isempty(records), keep=1:n; return; end
keep=round(double(records(:)'));
bad=keep<1 | keep>n;
if any(bad)
    error('readLabChart:badRecord', ...
        'This file has %d records; %s was asked for.',n,mat2str(unique(keep(bad))));
end
keep=unique(keep,'stable');
end

% =====================================================================
function [channels,factors]=readChannels(f,keep,offset,decim)
%readChannels  Every channel, its records concatenated onto the file's base.
%   THE TIME BASE IS BUILT FROM THE CHANNEL'S OWN SAMPLE PERIOD rather than read
%   back from the SDK, because the period is what makes the samples uniform: it
%   is exact inside a record, it is unambiguous about which record it belongs to,
%   and it cannot disagree with the .fs this function also reports.
%
%   DECIMATION IS APPLIED TO EACH RECORD AS IT ARRIVES, before it is appended, so
%   a long file never holds its full-rate concatenation - which is the point of
%   decimating here rather than afterwards.  factors comes back one per kept
%   channel: how many recorded samples each of its points now stands for.
channels=struct('name',{},'units',{},'fs',{},'data',{},'time',{},'records',{});
factors=[];
for c=1:1:f.n_channels
    cs=f.channel_specs(c);
    [d,t,recs,fs0,fs,un,fac]=deal([],[],[],[],[],'',1);
    for r=keep
        if ~f.isChannelInRecord(c,r), continue; end
        y=double(f.getChannelData(c,r,'return_object',false));
        y=y(:);
        if isempty(y), continue; end
        dt=perRecord(cs.dt,r);
        if ~isfinite(dt) || dt<=0, dt=1/max(perRecord(cs.fs,r),eps); end
        k=blockSize(dt,decim);
        [y,tr]=reduceRecord(y,dt,offset(r),k,decim);
        d=[d;y]; %#ok<AGROW>
        t=[t;tr]; %#ok<AGROW>
        recs(end+1)=r; %#ok<AGROW>
        if isempty(fs0)
            fs0=perRecord(cs.fs,r);              % as RECORDED - the rate to compare against
            fac=k;
            % The stated rate is the one the returned samples have.  With no
            % decimation that is the recorded rate unchanged; with it, the rate is
            % derived from the same dt the times were built from, so .fs and
            % diff(time) cannot disagree at the far end of a large block.
            fs=fs0; if k>1, fs=1/(k*dt); end
            un=char(unitOf(cs.units,r));
        elseif abs(perRecord(cs.fs,r)-fs0)>1e-9*max(fs0,1)
            warning('readLabChart:rateChanged', ...
                ['Channel "%s" was recorded at more than one sampling rate; ' ...
                 'its times are exact but its stated rate is the first one, %g Hz.'], ...
                char(cs.name),fs);
        end
    end
    if isempty(d), continue; end                 % nothing recorded: not a channel here
    if ~isempty(decim) && fac<=1
        % Said out loud, because a trace that stayed dense among decimated ones
        % looks like a fault otherwise.  A LabChart file routinely mixes a pressure
        % trace at 1 kHz with a temperature at 4 Hz, and the slow one has nothing to
        % give away.
        warning('readLabChart:notDecimated', ...
            ['Channel "%s" was recorded at %g Hz, which is already at or below the ' ...
             '%g Hz asked for, so it is kept as it was recorded.'], ...
            char(cs.name),double(fs),decim.fs);
    end
    channels(end+1)=struct('name',char(cs.name),'units',un,'fs',double(fs), ...
        'data',d,'time',t,'records',recs); %#ok<AGROW>
    factors(end+1)=fac; %#ok<AGROW>
end
end

% =====================================================================
function k=blockSize(dt,decim)
%blockSize  How many recorded samples one kept point stands for.
%   A WHOLE NUMBER, so the kept points stay uniformly spaced inside a record and
%   nothing has to be interpolated - which is why the rate that comes back is the
%   nearest achievable one rather than the one asked for.  1 means "keep this
%   record as it was": either no decimation was asked for, or the channel is
%   already at or below the requested rate and there is nothing to average.
k=1;
if isempty(decim) || ~isfinite(dt) || dt<=0, return; end
k=max(1,round(1/(decim.fs*dt)));
end

% =====================================================================
function [y,t]=reduceRecord(y,dt,t0,k,decim)
%reduceRecord  ONE record's samples and their times, decimated by k.
%   The blocks are consecutive and non-overlapping and the last one is short
%   rather than dropped, so the kept points cover the record from its first
%   sample to its last.  The short block is filled with NaN and the filter is
%   told to ignore it, which is also what makes a NaN already in the recording
%   cost its own sample rather than the whole block it sits in.
%
%   THE TIME OF A KEPT POINT IS THE MIDDLE OF ITS BLOCK - computed from the
%   indices rather than by averaging a time vector that would have to be built at
%   full rate first.  A block average timed at the block's START would slide the
%   whole trace half a block earlier, which at a factor of 1000 is visible.
n=numel(y);
if k<=1
    t=t0+(0:n-1)'*dt;
    return
end
nOut=ceil(n/k);
pad=nOut*k-n;
if pad>0, y=[y;nan(pad,1)]; end
y=reshape(y,k,nOut);
if strcmp(decim.filter,'median')
    y=median(y,1,'omitnan')';
else
    y=mean(y,1,'omitnan')';
end
first=(0:nOut-1)'*k;                     % first sample of each block, 0-based
last=min((1:nOut)'*k,n)-1;               % last sample of each block, 0-based
t=t0+dt*(first+last)/2;
end

% =====================================================================
function comments=snapComments(comments,channels,factors)
%snapComments  Move each comment onto a sample that still exists.
%   A comment is what a window is chosen against, and after decimation its
%   recorded time almost always falls BETWEEN two kept points - so it is moved to
%   the nearest of them, which is at most half a block away and is by definition
%   the resolution that was asked for.  A comment placed on a named channel is
%   snapped to THAT channel's points; a file-wide one goes to the finest grid
%   available, which is the smallest move that can be made.
%
%   ONLY WITHIN ONE BLOCK, though: a comment further than that from every kept
%   point is sitting in the stop between two records, and dragging it into a
%   neighbouring record would move it to a moment that is not the one it marks.
%   The order is restored afterwards, because two comments a fraction of a block
%   apart can otherwise swap when only one of them moves.
if isempty(comments) || isempty(channels) || ~any(factors>1), return; end
rates=[channels.fs];
rates(factors<=1)=-Inf;                  % only a decimated channel can be snapped to
[~,finest]=max(rates);
names=strtrim({channels.name});
for i=1:1:numel(comments)
    j=finest;
    k=find(strcmpi(strtrim(comments(i).channel),names),1);
    if ~isempty(k) && factors(k)>1, j=k; end
    tv=channels(j).time;
    [gap,idx]=min(abs(tv-comments(i).time));
    if gap<=1/channels(j).fs, comments(i).time=tv(idx); end
end
[~,ord]=sort([comments.time]);
comments=comments(ord);
end

% =====================================================================
function comments=readComments(f,keep,offset)
%readComments  The operator's marks, placed on the same base as the samples.
%   They are what a window is chosen against - 'drug in' is a comment, not a
%   feature of the trace - so a comment whose record was not read is dropped
%   rather than left at a time that no longer exists.
comments=struct('time',{},'text',{},'channel',{},'record',{});
if f.n_records==0, return; end
try
    T=f.getAllComments('output','table');
catch
    return                                      % a file with no comments at all
end
if isempty(T), return; end
names=f.channel_names;
for i=1:1:height(T)
    r=double(T.record(i));
    if ~any(keep==r), continue; end
    ch='';
    k=double(T.channel(i));                     % -1 = the comment is file-wide
    if k>=1 && k<=numel(names), ch=char(names{k}); end
    comments(end+1)=struct('time',offset(r)+double(T.time(i)), ...
        'text',char(string(T.text(i))),'channel',ch,'record',r); %#ok<AGROW>
end
if ~isempty(comments)
    [~,ord]=sort([comments.time]);
    comments=comments(ord);
end
end

% =====================================================================
function v=sdkVersion()
%sdkVersion  WHICH COPY OF THE SDK READ THIS FILE.  Upstream publishes no version
%   string, so the vendored folder is the answer; the commit it was taken at is
%   recorded in '3rd party/README.txt' beside it.
v='';
p=which('adi.readFile');
if ~isempty(p), v=fileparts(fileparts(p)); end
end

% =====================================================================
function v=perRecord(x,r)
%perRecord  A channel property that LabChart stores once per record.
if isempty(x), v=NaN; elseif isscalar(x), v=double(x); else, v=double(x(min(r,numel(x)))); end
end

% =====================================================================
function u=unitOf(units,r)
if isempty(units), u=''; return; end
if iscell(units), u=units{min(r,numel(units))}; else, u=units; end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
try
    v=st.(f);
catch
    v=d;
end
if isempty(v), v=d; end
end
%------------- END OF CODE --------------
