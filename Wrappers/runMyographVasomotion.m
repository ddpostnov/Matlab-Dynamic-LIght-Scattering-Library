%runMyographVasomotion  Wavelet vasomotion analysis of a myograph diameter trace
%
%   runMyographVasomotion(s,fNames) runs the vasomotion analysis on every interval
%   of every *_MYO_r.mat recording in fNames and writes the result into the
%   recording's own _MYO pair.  fNames is a cell array of *_MYO_r.mat paths; the
%   wrapper iterates them itself, so there is no launcher for-loop.
%
%   IT IS THE SAME ANALYSIS THE SPECKLE PIPELINE RUNS.  getMyographVasomotion calls
%   the shared core (getVasomotionMetrics) and the shared tree assembler, so a
%   myograph result and an LSCI result are the same tree of the same numbers, and
%   the defaults below are the speckle vasomotion step's defaults verbatim.  What a
%   myograph brings instead of vessel segments is the vessel's own diameter.
%
%   LINE-AVERAGED BY DEFAULT, PER LINE ON REQUEST.  By default the analysed signal
%   is the interval's line-averaged diameter - one oscillation for the vessel, which
%   is what a pressure myograph is set up to measure.  With s.perLine the analysis
%   runs on every measured image row instead, giving one unit per line: the exact
%   shape of a per-segment speckle result, and the same tree either way, so anything
%   that reads one reads the other.  BOTH ARE READ OUT OF THE WINDOW ITSELF -
%   results.intervals(k).diameter, the trace and the per-line array the diameter
%   step wrote - so a window is analysed on exactly the numbers stored under it.
%   Per line needs those arrays: a window written with s.keepArrays false has only
%   its trace, and this step says so rather than quietly answering with one unit.
%
%   WHICH DIAMETER IS ANALYSED.  The diameter step measures three (outer wall, wall
%   centre, lumen); s.diameterMeasures names the ones analysed here, wall centre by
%   default, and results nest per measure.  It is one protocol decision shared with
%   the propagation step.
%
%   A WIRE MYOGRAPH IS ANALYSED PER CHANNEL INSTEAD, and the CHANNEL IS THE OUTER
%   AXIS of its result tree - one LabChart file is several chambers, and a window may
%   belong to one of them:
%       results.channel(i).name
%       results.channel(i).intervals(k).vasomotion.<channel>
%   Each window is cut out of ITS CHANNEL'S OWN samples by time - the channels may be
%   sampled at different rates, so there is no shared frame index, and the samples
%   sit in results.channel(i) beside the windows they belong to.  The tree's own
%   field is the channel's name made into a legal one, and the real name is kept
%   beside it as .channelName so nothing is lost to the renaming; that keeps the
%   <VSM> shape identical to the speckle pipeline's results.vasomotion.<signal>,
%   which is what lets one set of plots and one exporter read both.  s.perLine does
%   not apply and is ignored.
%
%   INPUTS
%     s        parameter structure.  Carried through into
%              settings.runMyographVasomotion.  Defaults are filled when missing:
%                .diameterMeasures  which diameters to analyse, default {'mid'}
%                                   (pressure myograph only)
%                .perLine           false (one trace per interval) | true (one per
%                                   measured image row).  Ignored by a wire
%                                   myograph, which has one trace per channel
%                .vFR .cFR .wFR .wVPO      the frequency bands and the wavelet
%                .normalisation .normsize .tgtFS  how the trace is normalised and
%                                   how finely the spectrum is stored
%                .pcts .otsuMaxN .otsuElbow .nPeakProm   percentile bins and the
%                                   flare/silence clustering
%                .segVsmReturn      which analysis levels are computed and stored
%                .parforMyographLines  run the per-line loop on several workers
%              The speckle-only fields (vsmSignals, ppxVsmReturn) do not apply.
%     fNames   cell array of *_MYO_r.mat paths.  Empty cells are skipped.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   PRESSURE: results.intervals(k).vasomotion.<measure> - one
%                        vasomotion tree per analysed diameter measure (1 unit
%                        averaged, nY units per line).
%                        WIRE: results.channel(i).intervals(k).vasomotion.<channel>
%     <name>_MYO_s.mat   settings.runMyographVasomotion = s
%
%   EXAMPLE
%     s.diameterMeasures = {'mid'};
%     s.vFR = [0.05 0.25];   s.cFR = [0.4 0.6];
%     s.perLine = false;
%     D = dir(fullfile(rootFolder,'*_MYO_r.mat'));
%     runMyographVasomotion(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     getMyographVasomotion, myographMeasureIndex, myographChannelSamples,
%     myographProduct (Core/Myograph), the shared vasomotion core (Core/Vasomotion)
%     and Core/Reporting.  Needs the Wavelet Toolbox.
%
% See also: runMyographVideo, runLabChart, runMyographDiameter,
%           runMyographPropagation, getMyographVasomotion, runVasomotion,
%           myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - What is analysed
% s.diameterMeasures={'mid'};  % 'outer' | 'mid' | 'inner', one or several
% s.perLine=false;             % true = one result per measured image row
% s.vFR=[0.05 0.25];           % vasomotion band [lo hi], Hz
% s.cFR=[0.4 0.6];             % control band [lo hi], Hz
% %ADJUSTED IF NECESSARY - Analysis
% s.wFR=[0.01 1];              % full wavelet frequency range, Hz
% s.wVPO=10;                   % wavelet voices per octave
% s.tgtFS=1;                   % stored spectrum sampling rate, Hz ([] = full)

function runMyographVasomotion(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_MYO_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_MYO_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_MYO_r.mat'').']);
end
s=withDefaults(s);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Vasomotion',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    fName=char(fNames{fidx});
    reportFile(rep,fidx,fName);
    clearvars results settings

    [results,settings]=myographProduct('open',fName);
    sFile=s; sFile.fName=fName;

    if isWire(results)
        results=analyseChannels(sFile,results,fName);
    else
        if ~anyDiameter(results)
            error('runMyographVasomotion:noDiameter', ...
                'No diameter has been measured for %s yet - run the Diameter step first.',fName);
        end
        results=analyseDiameters(sFile,results,fName);
    end

    settings.runMyographVasomotion=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,results,settings);
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function results=analyseDiameters(s,results,fName)
%analyseDiameters  THE PRESSURE MYOGRAPH: one result per analysed diameter measure,
%   read out of the window's own diameter block.
rec=fieldOr(results,'recording',struct());
meas=cellstr(s.diameterMeasures);
for k=1:1:numel(results.intervals)
    d=fieldOr(results.intervals(k),'diameter',[]);
    if ~isstruct(d) || isempty(fieldOr(d,'time',[])), continue; end
    v=struct();
    for j=1:1:numel(meas)
        [~,nm]=myographMeasureIndex(meas{j},measuresOf(rec));
        sig=analysedSignal(d,nm,s.perLine,fName,results.intervals(k).name);
        if isempty(sig), continue; end
        v.(meas{j})=getMyographVasomotion(s,sig,double(d.time));
    end
    results.intervals(k).vasomotion=v;
end
end

% =====================================================================
function results=analyseChannels(s,results,fName)
%analyseChannels  THE WIRE MYOGRAPH: the CHANNEL is the outer axis, so the loop is
%   channel then window rather than window then channel.  Everything inside it is
%   the same analysis - the same core, the same bands, the same tree - which is what
%   makes a wire result, a pressure result and an LSCI result comparable.
%
%   EACH WINDOW IS CUT OUT OF ITS CHANNEL'S OWN SAMPLES, by TIME.  Channels may be
%   recorded at different rates, so there is no shared frame index to cut by; the
%   window's start and end are what the channel is asked for, and it answers with
%   the samples it actually has in that stretch.
%
%   THE TREE KEEPS ITS SIGNAL KEY even though the channel is now above it:
%   results.channel(i).intervals(k).vasomotion.<channel> is the same <VSM> shape as
%   results.vasomotion.<signal> in the speckle pipeline, and every plot and every
%   export in this library reads that shape.  Naming the axis twice costs a field;
%   breaking the shape would cost every consumer.
%
%   s.perLine DOES NOT APPLY HERE and is ignored: a channel is one trace, so there
%   are no lines to run it on.
%
%   THE WINDOW CARRIES ITS OWN SAMPLES.  The intervals step cut each channel to the
%   windows it was given and discarded the rest, so a window and the samples it is
%   analysed on are ONE element and there is no second list to match it against by
%   name.  A recording that has not been through that step yet still has the samples
%   whole on the channel, and intervalSamples reads either shape.
if ~isfield(results,'channel') || isempty(results.channel) || ...
        ~any(arrayfun(@(c) ~isempty(fieldOr(c,'intervals',[])),results.channel))
    warning('runMyographVasomotion:noWindows', ...
        ['No analysis windows are defined for %s, so nothing was analysed - ' ...
         'run the Intervals step on it first.'],fName);
    return
end
for i=1:1:numel(results.channel)
    nm=char(fieldOr(results.channel(i),'name',''));
    ivs=fieldOr(results.channel(i),'intervals',[]);
    [chData,~]=myographChannelSamples(results.channel(i));
    if isempty(chData)
        warning('runMyographVasomotion:noSamples', ...
            'Channel "%s" of %s carries no samples, so its windows were not analysed.', ...
            nm,fName);
        continue
    end
    fld=safeFields({nm});
    for k=1:1:numel(ivs)
        [sig,tt]=intervalSamples(results.channel(i),ivs(k));
        why=tooShortFor(sig,tt,s);
        if ~isempty(why)
            % A WINDOW CAN BE TOO SHORT TO MEAN ANYTHING, and a window on 'all
            % channels' reaches the slow channels too - a five-second window on a
            % 4 Hz channel is twenty samples.  It is said out loud and passed over
            % rather than pushed into the wavelet, which would either throw and take
            % the whole recording down with it or answer with a spectrum of one
            % point that nobody would question.
            warning('runMyographVasomotion:windowTooShort', ...
                'Window "%s" on channel "%s" of %s %s, so it was not analysed.', ...
                char(fieldOr(ivs(k),'name',sprintf('interval%d',k))),nm,fName,why);
            ivs(k).vasomotion=struct();
            continue
        end
        v=struct();
        v.(fld{1})=getMyographVasomotion(s,sig,tt);
        % The channel's own name is kept beside its tree: the field had to be made
        % into a legal MATLAB name, and 'Force 1 (mN)' does not survive that.
        v.(fld{1}).channelName=nm;
        ivs(k).vasomotion=v;
    end
    results.channel(i).intervals=ivs;
end
end

% =====================================================================
function fld=safeFields(names)
%safeFields  Channel names as struct fields.  makeValidName does the legalising and
%   makeUniqueStrings does what it cannot: 'Force 1' and 'Force-1' both legalise to
%   'Force_1', and silently overwriting one result with the other is exactly the
%   kind of loss the original names are kept to prevent.
fld=matlab.lang.makeValidName(reshape(cellstr(names),1,[]));
fld=matlab.lang.makeUniqueStrings(fld,{},namelengthmax);
end

% =====================================================================
function tf=isWire(results)
%isWire  A wire recording is one that has a CHANNEL AXIS.  The modality is written in
%   results.recording, but the test that decides what this wrapper can DO is what the
%   file holds, so it is asked of the data.  It is asked of the axis and not of the
%   samples: the intervals step moves the samples off the channel and into its
%   windows, and a recording that has been through it would otherwise be read as a
%   pressure myograph with nothing measured.
tf=false;
ch=fieldOr(results,'channel',[]);
if isempty(ch) || ~isstruct(ch), return; end
tf=true;
end

% =====================================================================
function tf=anyDiameter(results)
%anyDiameter  Has anything been measured on this pressure recording yet?
tf=false;
ivs=fieldOr(results,'intervals',[]);
if isempty(ivs) || ~isstruct(ivs), return; end
tf=any(arrayfun(@(iv) isstruct(fieldOr(iv,'diameter',[])) && ...
    ~isempty(fieldOr(iv.diameter,'time',[])),ivs));
end

% =====================================================================
function meas=measuresOf(rec)
%measuresOf  The measures, and their order.  results.recording is the one authority
%   on it, written by the diameter step in the same breath as the windows.
meas=reshape(cellstr(fieldOr(rec,'measures',{})),1,[]);
end

% =====================================================================
function why=tooShortFor(sig,tt,s)
%tooShortFor  Is there anything in this window to measure an oscillation in?  Stated
%   as what an ANSWER NEEDS rather than as what the wavelet happens to throw on:
%   a spectrum down to s.wFR(1) needs the window to hold at least one full cycle of
%   that frequency, and any spectrum at all needs enough samples to have one.  Below
%   either, what comes back is not a weak result, it is not a result.
%   Returns the reason as a phrase for the warning, or '' when the window is fine.
why='';
n=numel(sig);
if n<32
    why=sprintf('holds only %d samples',n); return
end
span=tt(end)-tt(1);
fLo=min(s.wFR);
if ~isfinite(fLo) || fLo<=0, return; end
if span*fLo<1
    why=sprintf(['is %.1f s long, less than one cycle of the %.3g Hz low end of ' ...
        'the analysed band'],span,fLo);
end
end

% =====================================================================
function [sig,tt]=intervalSamples(c,iv)
%intervalSamples  This channel's own samples inside this window.  THE WINDOW CARRIES
%   THEM once the intervals step has run: it cut the channel to the windows and threw
%   the rest away, so a window's samples are the window's own and nothing has to be
%   selected out of anything.  A recording that has not been through that step yet
%   still has them whole on the channel, and is cut here by TIME - the channels may
%   be sampled at different rates, so there is no shared index to cut by.
sig=[]; tt=[];
sm=fieldOr(iv,'samples',[]);
if isstruct(sm) && isscalar(sm) && ~isempty(fieldOr(sm,'data',[]))
    sig=double(sm.data(:));
    tt=double(sm.time(:));
    return
end
t=double(fieldOr(c,'time',[])); t=t(:);
t0=double(fieldOr(iv,'tStart',[])); t1=double(fieldOr(iv,'tEnd',[]));
if isempty(t) || isempty(t0) || isempty(t1), return; end
idx=t>=min(t0,t1) & t<=max(t0,t1);
sig=double(c.data(idx)); sig=sig(:);
tt=t(idx); tt=tt(:);
end

% =====================================================================
function s=withDefaults(s)
%withDefaults  The speckle vasomotion step's defaults, verbatim, plus the two
%   choices a myograph adds.  Filled only where the caller said nothing, so a
%   protocol that sets a band keeps it.
def.diameterMeasures={'mid'};   def.perLine=false;
def.vFR=[0.05 0.25];            def.cFR=[0.4 0.6];
def.wFR=[0.01 1];               def.wVPO=10;
def.normalisation='median';     def.normsize=101;      def.tgtFS=1;
def.pcts=0:10:100;              def.otsuMaxN=5;        def.otsuElbow=0.05;
def.nPeakProm=0.10;
def.segVsmReturn={'bands','moments','series','clustering','spectrum'};
def.parforMyographLines=true;
fnd=fieldnames(def);
for i=1:1:numel(fnd)
    if ~isfield(s,fnd{i}) || (isempty(s.(fnd{i})) && ~islogical(def.(fnd{i})))
        s.(fnd{i})=def.(fnd{i});
    end
end
end

% =====================================================================
function sig=analysedSignal(d,name,perLine,fName,ivName)
%analysedSignal  THE trace this interval is analysed on, out of its own diameter
%   block.
%   Line-averaged, it is the trace the diameter step already wrote - read back
%   rather than recomputed, so the spectrum belongs to exactly the curve the
%   diameter page plotted.  Per line, it is .lines.<measure>, which holds the
%   measured rows and only those: the rows outside the measured band carry an
%   interpolated fill and are not oscillations of anything.
%
%   PER LINE NEEDS THE ARRAYS, and a window written with s.keepArrays false has
%   none - there is no source copy to fall back on any more.  Said out loud and
%   passed over: silently answering with the pooled trace instead would report one
%   unit where the protocol asked for one per line, which is a different result
%   wearing the right shape.
sig=[];
if perLine
    if ~isfield(d,'lines') || ~isfield(d.lines,name)
        warning('runMyographVasomotion:noPerLineDiameter', ...
            ['Window "%s" of %s has no per-line diameter stored, so it could not be ' ...
             'analysed per line.  Run the Diameter step on it with ''Keep every ' ...
             'line''''s diameter'' switched on.'],char(ivName),fName);
        return
    end
    sig=double(d.lines.(name));
    return
end
if isfield(d,name), sig=double(d.(name)(:)); end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
