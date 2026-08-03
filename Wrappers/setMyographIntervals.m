%setMyographIntervals  Interactive, per-file definition of a myograph's analysis windows
%
%   setMyographIntervals(s,fNames) opens the interval editor on every *_MYO_r.mat
%   recording in fNames, on the diameter that has already been measured, and writes
%   the windows the operator defines back into that recording's own _MYO triplet.
%   ONE WINDOW PER RECORDING, and nothing advances until Done is pressed - the
%   wrapper loops the files itself, exactly as setRegions does, so there is no file
%   dropdown to get lost in and no launcher for-loop.  fNames is a cell array of
%   *_MYO_r.mat paths.
%
%   WHAT AN INTERVAL IS.  A named stretch of time that the later steps analyse on
%   its own: the propagation and the vasomotion run per interval, and the Excel
%   export and the explorer stratify by it.  It is how 'baseline', 'drug' and
%   'washout' become three answers from one recording instead of one average across
%   all three.
%
%   THE EDITED SET REPLACES WHAT WAS THERE, AND THE DIAMETER BLOCK IS RE-CUT WITH IT.
%   Every interval that comes back is re-sliced out of source: its .frames range into
%   source.time and its .diameter - the line-averaged traces and their statistics -
%   are rebuilt from the boundaries as they now stand, so the numbers under an
%   interval always describe the window that carries them.  ANY PROPAGATION OR
%   VASOMOTION AN INTERVAL CARRIED IS DROPPED for the same reason: its window moved,
%   so that analysis no longer describes it.  Re-run those steps.
%
%   THE RECORDING IS NOT RE-MEASURED AND NOT RE-CUT ON DISK.  source is read, never
%   written: the diameter was measured once, over the analysed span, and an interval
%   is a range into it.  A WINDOW CANNOT LEAVE THAT SPAN - after a time crop the trace
%   is only the cropped stretch, and the editor clips any band that would run past
%   either end - so an interval always names frames that exist.  Where the recording
%   jumps (disjoint pre-set intervals) the trace is drawn with a break, so the
%   stretches that exist are the ones on screen.
%
%   IT SERVES BOTH MYOGRAPHS, and the difference between them is ONE switch: what
%   the recording's own source.modality says it is.  A pressure myograph opens the
%   editor on the measured diameter with the video beside it; a wire myograph opens
%   the SAME editor on the recorded channels, with the LabChart comments beside them.
%   Everything else - the loop, the re-cut, the save, the settings - is one piece of
%   code, because an interval means the same thing in both.
%
%   A WIRE MYOGRAPH'S WINDOWS ARE STORED PER CHANNEL.  One LabChart file is several
%   chambers, and a window may belong to one of them - the drug went into chamber 2
%   at a different minute from chamber 3 - so the editor's table carries a CHANNEL
%   column ('all channels', or one named channel) and the result tree makes the
%   channel its outer axis:
%       results.channel(i).name        the channel
%       results.channel(i).intervals(k)   the windows analysed on it
%   A window left on 'all channels' appears under every channel, because that is
%   what it means; results.intervals stays empty for a wire recording, so there is
%   one place a window lives and not two.  A pressure myograph is unchanged: one
%   vessel, one flat results.intervals, and no results.channel.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.setMyographIntervals.
%       .edgeMode  (optional, default 'mid') which of the three measured diameters
%                  the trace shows - 'outer', 'mid' (the wall centre) or 'inner'.
%                  It is the ONE thing a protocol fixes about this step; everything
%                  else about it is the drawing, and the drawing lives in the
%                  results, the same split setRegions makes between s.nRegions and
%                  results.regionsMask.  A wire myograph ignores it: it has no
%                  diameter, and which channel a window is analysed on is chosen in
%                  the window's own table.
%       .drawPoints (optional, default 2000) HOW MUCH OF THE RECORDING IS DRAWN in
%                  the editor - the number of min/max buckets the visible span is
%                  cut into, per channel.  A LabChart file is tens of millions of
%                  samples and drawing them all is what made the window unusable;
%                  this is the ceiling on what reaches the renderer, re-applied
%                  every time the operator zooms, so detail comes back as they go
%                  in.  IT AFFECTS THE DRAWING ONLY.  The intervals are cut out of
%                  the raw samples and every analysis reads those, so lowering it on
%                  a slow machine changes what the operator sees and NOTHING about
%                  what comes out.
%     fNames   cell array of *_MYO_r.mat paths.  Empty cells are skipped.
%     fNamesRaw (optional) the matching raw recordings, for the preview.  Omitted,
%              each recording's own source.fName is used.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files, so Stop takes effect at the next
%     recording rather than in the middle of one).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   PRESSURE: results.intervals REPLACED with .name .tStart
%                        .tEnd .frames .diameter for every window defined.
%                        WIRE: results.channel REPLACED with .name .intervals, and
%                        results.intervals emptied.
%                        Either way .propagation and .vasomotion are cleared
%     <name>_MYO_s.mat   settings.setMyographIntervals = s
%
%   EXAMPLE
%     s.edgeMode = 'mid';
%     D = dir(fullfile(rootFolder,'*_MYO_r.mat'));
%     setMyographIntervals(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     editMyographIntervals, getMyographTrace, cutMyographSource,
%     myographIntervals, myographProduct (Core/Myograph) and Core/Reporting.
%
% See also: setMyographPresetIntervals, setMyographCrop, runMyographDiameter,
%           runMyographPropagation, runMyographVasomotion, editMyographIntervals,
%           runLabChart
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL
% s.edgeMode='mid';       % which diameter the trace shows: outer / mid / inner
% %ADJUSTED IF NECESSARY - Drawing only, never the data
% s.drawPoints=2000;      % points drawn per channel across the visible span

function setMyographIntervals(s,fNames,fNamesRaw)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_MYO_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_MYO_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_MYO_r.mat'').']);
end
if nargin<3, fNamesRaw={}; end
if ~isfield(s,'edgeMode') || isempty(s.edgeMode), s.edgeMode='mid'; end
if ~isfield(s,'drawPoints') || isempty(s.drawPoints), s.drawPoints=2000; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Intervals',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    fName=char(fNames{fidx});
    reportFile(rep,fidx,fName);
    clearvars source results settings

    [source,results,settings]=myographProduct('open',fName);
    % THE ONE PLACE THE TWO MYOGRAPHS DIFFER IN THIS WRAPPER.  Everything after it -
    % the editor, the re-cut, the save - is the same code for both, because an
    % interval is the same thing in both: a named stretch of time.
    switch upper(char(fieldOr(source,'modality','PMYO')))
        case 'WMYO'
            mode='WMYO';
            comments=fieldOr(results,'comments',[]);
        otherwise
            mode='PMYO';
            comments=[];
            if isempty(fieldOr(source,'data',[]))
                error('setMyographIntervals:noDiameter', ...
                    ['No diameter has been measured for %s yet - run the Diameter step ' ...
                     'first, or choose the windows before it with Pre-set intervals.'],fName);
            end
    end
    sFile=s; sFile.fName=fName;

    data=getMyographTrace(source,sFile,rawRecording(fNamesRaw,fidx),comments);
    [ivT,names,extra]=editMyographIntervals(data,editorOptions(results,fName,mode,s.drawPoints));

    % BOTH ARE ALWAYS WRITTEN, and exactly one of them is filled: a pressure myograph
    % has one vessel and a flat list of windows, a wire myograph has several chambers
    % and the CHANNEL is the outer axis.  Assigning the pair rather than branching is
    % what keeps a recording that changed modality (it cannot, but a file copied over
    % another can) from keeping the other one's windows.
    [results.intervals,results.channel]=cutMyographSource(source,ivT,names,{extra.channels});
    settings.setMyographIntervals=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,[],results,settings);   % SOURCE untouched
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function opts=editorOptions(results,fName,mode,drawPoints)
%editorOptions  How the editor opens on this recording: pre-loaded with the windows
%   it already carries, titled with the recording so a run of files never leaves any
%   doubt about which one is on screen.
[ivT,names,chans]=intervalsOf(results);
[~,stem]=fileparts(regexprep(fName,'_MYO_[drs]\.mat$',''));
opts=struct('mode',mode,'intervals',ivT,'names',{names},'channels',{chans}, ...
    'drawPoints',drawPoints,'title',['Intervals - ' stem]);
end

% =====================================================================
function [ivT,names,chans]=intervalsOf(results)
%intervalsOf  The windows already defined, as the editor pre-loads them.  A window
%   with no start or end is not a window yet and is not offered.
%
%   A WIRE MYOGRAPH'S WINDOWS ARE STORED PER CHANNEL, and a window left on 'all
%   channels' therefore sits under every one of them.  Offering the flat list back
%   would re-open the recording with four copies of it, so the copies are folded
%   here: a window whose own .channels is empty is the same window wherever it was
%   found, and is identified by its name and its boundaries.
ivT=zeros(0,2); names={}; chans={};
ivs=myographIntervals(results);
seen={};
for k=1:1:numel(ivs)
    t0=fieldOr(ivs(k),'tStart',[]); t1=fieldOr(ivs(k),'tEnd',[]);
    if isempty(t0) || isempty(t1), continue; end
    nm=char(fieldOr(ivs(k),'name','')); if isempty(nm), nm=sprintf('interval%d',k); end
    ch=fieldOr(ivs(k),'channels',{});
    if isempty(ch)
        key=sprintf('%s|%.9g|%.9g',nm,double(t0),double(t1));
        if any(strcmp(seen,key)), continue; end
        seen{end+1}=key; %#ok<AGROW>
    end
    ivT(end+1,:)=[double(t0) double(t1)]; %#ok<AGROW>
    names{end+1}=nm; %#ok<AGROW>
    chans{end+1}=ch; %#ok<AGROW>
end
end

% =====================================================================
function raw=rawRecording(fNamesRaw,fidx)
%rawRecording  The raw path the caller passed for this file, if any.  Empty falls
%   back to the one the entry step recorded (getMyographTrace resolves it).
raw='';
if numel(fNamesRaw)>=fidx && ~isempty(fNamesRaw{fidx}), raw=char(fNamesRaw{fidx}); end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
