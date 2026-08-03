%setMyographPresetIntervals  Choose a myograph's analysis windows BEFORE measuring it
%
%   setMyographPresetIntervals(s,fNames) opens the interval editor on every
%   *_MYO_r.mat recording in fNames while nothing has been measured yet, and writes
%   the windows the operator picks - names and times, nothing else - into that
%   recording's own _MYO triplet.  ONE WINDOW PER RECORDING, blocking, no file
%   dropdown: the wrapper loops the files, exactly as setRegions does.
%
%   WHY BEFORE.  The diameter step then measures ONLY inside these windows, and
%   source.time and source.data come back only that long (see runMyographDiameter's
%   order of precedence).  Measuring four 5-minute windows of a two-hour recording
%   costs twenty minutes of work and twenty minutes of memory rather than two hours
%   of both.  Everything outside them is never read.
%
%   WHAT THE OPERATOR IS LOOKING AT.  There is no diameter to plot yet, so the trace
%   is a per-frame BRIGHTNESS PROFILE of the recording - a coarse pre-pass, read
%   before the window opens, that shows where the protocol changed - and the video
%   preview carries no wall overlay, because no wall has been detected.  The video is
%   the point of this window: the windows are recognised in the recording and only
%   confirmed on the curve.  getMyographTrace documents what the pre-pass costs.
%
%   IT IS AN ALTERNATIVE TO THE TIME CROP, AND TO THE INTERVAL EDITOR.  One set of
%   windows per recording: choosing them here means not cropping (a crop is one
%   window chosen the same way) and not defining them again afterwards on the
%   measured diameter.  The Constructor enforces that - ticking this step unticks
%   both of the others.
%
%   INPUTS
%     s        parameter structure.  Carried through into
%              settings.setMyographPresetIntervals.
%       .profileSamples (optional, default 1200) how many points the brightness
%                  profile has.  More is a finer curve and a longer wait before the
%                  window opens; the wait grows with this number, not with the length
%                  of the recording.
%     fNames   cell array of *_MYO_r.mat paths.  Empty cells are skipped.
%     fNamesRaw (optional) the matching raw recordings.  Omitted, each recording's
%              own source.fName is used.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   results.intervals REPLACED with .name .tStart .tEnd only -
%                        no .frames and no .diameter, because nothing has been
%                        measured; results.timeCrop cleared, since these windows ARE
%                        the analysed span now
%     <name>_MYO_s.mat   settings.setMyographPresetIntervals = s
%
%   EXAMPLE
%     D = dir(fullfile(rootFolder,'*_MYO_r.mat'));
%     setMyographPresetIntervals(struct(), fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     editMyographIntervals, getMyographTrace, myographProduct (Core/Myograph) and
%     Core/Reporting.
%
% See also: setMyographCrop, setMyographIntervals, runMyographDiameter,
%           editMyographIntervals, getMyographTrace
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY
% s.profileSamples=1200;  % points in the brightness profile the windows are picked on

function setMyographPresetIntervals(s,fNames,fNamesRaw)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_MYO_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_MYO_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_MYO_r.mat'').']);
end
if nargin<3, fNamesRaw={}; end
if ~isfield(s,'profileSamples') || isempty(s.profileSamples), s.profileSamples=1200; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Pre-set intervals',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    fName=char(fNames{fidx});
    reportFile(rep,fidx,fName);
    clearvars source results settings

    [source,results,settings]=myographProduct('open',fName);
    sFile=s; sFile.fName=fName;

    data=getMyographTrace(source,sFile,rawRecording(fNamesRaw,fidx));
    [ivT,names]=editMyographIntervals(data,editorOptions(results,fName));

    results.intervals=windowsOnly(ivT,names);
    % These windows ARE the analysed span, so a crop left over from an earlier
    % decision would be a second, contradictory answer to the same question.
    results.timeCrop=[];
    settings.setMyographPresetIntervals=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,[],results,settings);   % SOURCE untouched
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function opts=editorOptions(results,fName)
%editorOptions  How the editor opens on this recording: pre-loaded with whatever
%   windows it already carries, titled with the recording so a run of files never
%   leaves any doubt about which one is on screen.
[ivT,names]=intervalsOf(results);
[~,stem]=fileparts(regexprep(fName,'_MYO_[drs]\.mat$',''));
opts=struct('mode','PMYO','intervals',ivT,'names',{names}, ...
    'title',['Pre-set intervals - ' stem]);
end

% =====================================================================
function [ivT,names]=intervalsOf(results)
%intervalsOf  The windows already defined, as the editor pre-loads them.
ivT=zeros(0,2); names={};
ivs=[];
if isfield(results,'intervals'), ivs=results.intervals; end
for k=1:1:numel(ivs)
    t0=fieldOr(ivs(k),'tStart',[]); t1=fieldOr(ivs(k),'tEnd',[]);
    if isempty(t0) || isempty(t1), continue; end
    ivT(end+1,:)=[double(t0) double(t1)]; %#ok<AGROW>
    nm=char(fieldOr(ivs(k),'name','')); if isempty(nm), nm=sprintf('interval%d',k); end
    names{end+1}=nm; %#ok<AGROW>
end
end

% =====================================================================
function out=windowsOnly(ivT,names)
%windowsOnly  The interval elements this step may write: a name and a window, and
%   nothing else.  .frames and .diameter belong to the measurement, which has not
%   happened; leaving them EMPTY rather than guessing them is what lets the diameter
%   step tell 'a window to measure' from 'a window already measured'.
out=myographProduct('intervals');
for k=1:1:size(ivT,1)
    out(k).name=char(names{k});
    out(k).tStart=ivT(k,1);
    out(k).tEnd=ivT(k,2);
    out(k).channels={};
    out(k).frames=[];
    out(k).diameter=[];
    out(k).propagation=[];
    out(k).vasomotion=[];
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
