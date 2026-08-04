%setMyographCrop  Choose the one stretch of a myograph recording that is analysed
%
%   setMyographCrop(s,fNames) opens the interval editor on every *_MYO_r.mat
%   recording in fNames restricted to a SINGLE window, and records that window as
%   results.timeCrop.  ONE WINDOW PER RECORDING, blocking, no file dropdown: the
%   wrapper loops the files, exactly as setRegions does.
%
%   IT IS A DECISION, NOT A CUT.  This step runs before anything has been measured,
%   so it shortens nothing: it records which window matters, and the diameter step
%   is what reads only those frames and measures only that long.
%
%   CHANGE THE CROP AFTER A DIAMETER EXISTS AND SOMETHING IS DISCARDED, because the
%   measurement lives in the windows (results.intervals(k).diameter) and there is no
%   source copy behind them.  Which thing is discarded depends on where the new crop
%   falls, and the difference is worth stating:
%
%     A CROP INSIDE WHAT WAS ALREADY MEASURED NARROWS IT.  The windows are cut to the
%     new crop and the frames outside it go; the measurement inside it is kept, and
%     nothing has to be measured again.  Tightening a crop after a first look is the
%     ordinary reason to come back to this step, and it costs the frames it drops
%     rather than the whole run.
%
%     A CROP THAT REACHES ANYWHERE ELSE DISCARDS THE MEASUREMENT.  The stretch it
%     asks for was never read, and a window that answered only the part of it that
%     happens to exist would under-report without saying so.  The windows are
%     cleared, and the diameter step measures the new one.
%
%   Either way it cannot be undone from the product: only the recording can bring
%   the frames back, so the recording has to be beside the results before this step
%   will discard anything, and a recording that travelled without it is left alone
%   and said so by name.  A crop confirmed unchanged discards nothing and is silent,
%   as is one chosen before anything has been measured - which is the ordinary case,
%   this step being the one that runs first.  Re-opened after the measurement was
%   cleared, the editor is back to the brightness profile, because that is honestly
%   all the product then holds.
%
%   THE CROP IS IN ABSOLUTE SECONDS from the start of the recording, like every other
%   time in this product, so it can always be located back in the original footage.
%
%   IT IS AN ALTERNATIVE TO PRE-SET INTERVALS.  A crop is one window; pre-set
%   intervals are several.  Choosing one means not choosing the other, and the
%   Constructor enforces it - ticking either unticks the other.
%
%   NO WINDOW DRAWN MEANS THE WHOLE RECORDING.  results.timeCrop is then emptied,
%   which is exactly what "no crop" means to the diameter step, in the same way that
%   drawing no ROI in setRegions means the whole field of view.
%
%   WHAT THE OPERATOR IS LOOKING AT.  Before the diameter exists the trace is a
%   per-frame brightness profile of the recording and the video preview carries no
%   wall overlay (see getMyographTrace).  Re-opened after a diameter has been
%   measured, it is the measured windows laid end to end, with the walls drawn on
%   the preview and a break in the line wherever the recording was not measured.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.setMyographCrop.
%       .profileSamples (optional, default 1200) how many points the brightness
%                  profile has, when there is no diameter yet to plot instead.
%       .edgeMode  (optional, default 'mid') which diameter is plotted when there IS
%                  one.
%     fNames   cell array of *_MYO_r.mat paths.  Empty cells are skipped.
%     fNamesRaw (optional) the matching raw recordings.  Omitted, each recording's
%              own results.recording.fName is used.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   results.timeCrop = [t0 t1] (or [] when nothing was drawn);
%                        results.intervals NARROWED to the new crop when it lies
%                        inside a measured window, and CLEARED when it does not
%     <name>_MYO_s.mat   settings.setMyographCrop = s
%
%   EXAMPLE
%     D = dir(fullfile(rootFolder,'*_MYO_r.mat'));
%     setMyographCrop(struct(), fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     editMyographIntervals, getMyographTrace, myographCropWindow,
%     cutMyographIntervals, myographRecordingPath, myographProduct (Core/Myograph)
%     and Core/Reporting.  The narrow-or-clear decision itself is myographCropWindow:
%     it is what decides whether a measurement is kept or thrown away, and this
%     wrapper blocks on the editor, so the decision lives where a test can reach it.
%
% See also: setMyographPresetIntervals, setMyographIntervals, runMyographDiameter,
%           editMyographIntervals, getMyographTrace, myographCropWindow
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY
% s.profileSamples=1200;  % points in the brightness profile the window is picked on

function setMyographCrop(s,fNames,fNamesRaw)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_MYO_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_MYO_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_MYO_r.mat'').']);
end
if nargin<3, fNamesRaw={}; end
if ~isfield(s,'profileSamples') || isempty(s.profileSamples), s.profileSamples=1200; end
if ~isfield(s,'edgeMode') || isempty(s.edgeMode), s.edgeMode='mid'; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Time crop',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    fName=char(fNames{fidx});
    reportFile(rep,fidx,fName);
    clearvars results settings

    [results,settings]=myographProduct('open',fName);
    sFile=s; sFile.fName=fName;
    raw=rawRecording(fNamesRaw,fidx);

    % THE ONE-WAY CHECK, BEFORE THE EDITOR AND ONLY WHEN THERE IS SOMETHING TO LOSE.
    % A crop chosen before anything has been measured discards nothing, which is what
    % this step usually is; a crop changed after a diameter exists cannot be taken
    % back without the recording, so the recording has to be there first.
    if measuredFrames(results)>0
        [p,wanted]=myographRecordingPath(results,fName,raw);
        if isempty(p)
            warning('setMyographCrop:recordingGone', ...
                ['Changing the analysed window throws away the measurements outside ' ...
                 'it.  Keep %s beside the results so they can be measured again, or ' ...
                 'leave the window as it is.'],wanted);
            continue
        end
    end

    data=getMyographTrace(results,sFile,raw);
    ivT=editMyographIntervals(data,editorOptions(results,fName));

    was=fieldOr(results,'timeCrop',[]);
    crop=cropOf(ivT);
    results.timeCrop=crop;
    if ~isequal(was,crop)
        results=applyCrop(results,crop,fName);
    end
    settings.setMyographCrop=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,results,settings);
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function results=applyCrop(results,crop,fName)
%applyCrop  What a moved crop does to the measurement already stored.  NARROW IT when
%   the new window lies inside one stretch that was measured - the frames outside go
%   and the rest is kept, which is what tightening a crop should cost.  CLEAR IT
%   otherwise: the new window reaches into footage that was never read, and a window
%   answering only the part of it that exists would under-report in silence.
%   Nothing at all measured is the ordinary case and is silent.
before=measuredFrames(results);
if before==0, results.intervals=myographProduct('intervals'); return; end
[~,b,e]=fileparts(char(fName)); shortName=[b e];

% THE DECISION ITSELF IS A CORE (myographCropWindow), not twelve lines here.  It is
% what decides whether a measurement is narrowed or thrown away, and this wrapper
% opens the interval editor - which blocks until Done - so a decision left inside it
% could only ever be reached by hand.
[ivT,names]=myographCropWindow(results,crop);
if ~isempty(ivT)
    [results.intervals,~,tally]=cutMyographIntervals(results,ivT,names);
    if tally.kept<tally.before
        warning('setMyographCrop:framesDiscarded', ...
            ['The analysed window of %s was narrowed, so its measurement now keeps ' ...
             '%d of the %d frames that were there; the rest was discarded and only ' ...
             'the recording can bring it back.'],shortName,tally.kept,tally.before);
    end
    return
end

% The analysed window moved outside what was measured, so every interval defined
% against the old one names frames that were never read.  They go; the diameter step
% then measures the new window as one interval, or the operator defines new ones.
results.intervals=myographProduct('intervals');
warning('setMyographCrop:intervalsCleared', ...
    ['The analysed window of %s moved outside what was measured, so the windows ' ...
     'defined on the old one were cleared, and the %d frames measured inside them ' ...
     'went with them.  Measure the diameter again before analysing it.'], ...
    shortName,before);
end

% =====================================================================
function n=measuredFrames(results)
%measuredFrames  How many frames of this recording are stored, over every window.
%   It is the number a discard is weighed against, and zero is what "nothing has been
%   measured yet" looks like from here.
n=0;
ivs=fieldOr(results,'intervals',[]);
if isempty(ivs) || ~isstruct(ivs), return; end
for k=1:1:numel(ivs)
    n=n+numel(double(fieldOr(fieldOr(ivs(k),'diameter',[]),'time',[])));
end
end

% =====================================================================
function opts=editorOptions(results,fName)
%editorOptions  ONE window, pre-loaded with the crop this recording already carries.
%   maxIntervals = 1 is what makes Add grey out once the window exists: a crop that
%   could be two windows would be a set of pre-set intervals, which is the other step.
crop=fieldOr(results,'timeCrop',[]);
ivT=zeros(0,2); names={};
if numel(crop)==2, ivT=double(crop(:))'; names={'analysed window'}; end
[~,stem]=fileparts(regexprep(fName,'_MYO_[drs]\.mat$',''));
opts=struct('mode','PMYO','intervals',ivT,'names',{names},'maxIntervals',1, ...
    'title',['Time crop - ' stem]);
end

% =====================================================================
function c=cropOf(ivT)
%cropOf  The single window as [t0 t1], or EMPTY when none was drawn - which is what
%   "the whole recording" means to every reader of results.timeCrop.
c=[];
if isempty(ivT), return; end
c=[min(ivT(1,:)) max(ivT(1,:))];
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
