%setMyographCrop  Choose the one stretch of a myograph recording that is analysed
%
%   setMyographCrop(s,fNames) opens the interval editor on every *_MYO_r.mat
%   recording in fNames restricted to a SINGLE window, and records that window as
%   results.timeCrop.  ONE WINDOW PER RECORDING, blocking, no file dropdown: the
%   wrapper loops the files, exactly as setRegions does.
%
%   IT IS A DECISION, NOT A CUT.  This step runs before anything has been measured,
%   so it does not shorten source: it records which window matters, and the diameter
%   step is what reads only those frames and writes a source.time / source.data that
%   long.  CHANGE THE CROP AFTER A DIAMETER EXISTS AND THE DIAMETER MUST BE MEASURED
%   AGAIN - the arrays on disk cover the old window and nothing can re-slice them
%   onto a window they never covered.  Because of that, changing the crop CLEARS the
%   intervals that were defined against the old one: they name frames of a recording
%   that is about to become a different length.  A crop confirmed unchanged clears
%   nothing.
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
%   measured, it is the measured diameter, with the walls drawn on the preview.
%
%   INPUTS
%     s        parameter structure.  Carried through into settings.setMyographCrop.
%       .profileSamples (optional, default 1200) how many points the brightness
%                  profile has, when there is no diameter yet to plot instead.
%       .edgeMode  (optional, default 'mid') which diameter is plotted when there IS
%                  one.
%     fNames   cell array of *_MYO_r.mat paths.  Empty cells are skipped.
%     fNamesRaw (optional) the matching raw recordings.  Omitted, each recording's
%              own source.fName is used.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   results.timeCrop = [t0 t1] (or [] when nothing was drawn);
%                        results.intervals CLEARED when the crop changed
%     <name>_MYO_s.mat   settings.setMyographCrop = s
%
%   EXAMPLE
%     D = dir(fullfile(rootFolder,'*_MYO_r.mat'));
%     setMyographCrop(struct(), fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     editMyographIntervals, getMyographTrace, myographProduct (Core/Myograph) and
%     Core/Reporting.
%
% See also: setMyographPresetIntervals, setMyographIntervals, runMyographDiameter,
%           editMyographIntervals, getMyographTrace
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

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
    clearvars source results settings

    [source,results,settings]=myographProduct('open',fName);
    sFile=s; sFile.fName=fName;

    data=getMyographTrace(source,sFile,rawRecording(fNamesRaw,fidx));
    ivT=editMyographIntervals(data,editorOptions(results,fName));

    was=fieldOr(results,'timeCrop',[]);
    crop=cropOf(ivT);
    had=~isempty(fieldOr(results,'intervals',[]));
    results.timeCrop=crop;
    if ~isequal(was,crop)
        % The analysed window moved, so every interval defined against the old one
        % names frames of a recording that no longer exists.  They go; the diameter
        % step then measures the new window as one interval, or the operator defines
        % new ones on it.  Nothing to clear is the ordinary case - a crop chosen
        % before anything was measured - and it is silent.
        results.intervals=myographProduct('intervals');
        if had
            warning('setMyographCrop:intervalsCleared', ...
                ['The analysed window of %s changed, so the intervals defined on the ' ...
                 'old one were cleared.  Measure the diameter again before ' ...
                 'analysing it.'],fName);
        end
    end
    settings.setMyographCrop=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,[],results,settings);   % SOURCE untouched
    reportSaved(rep);
end
reportClose(rep);
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
