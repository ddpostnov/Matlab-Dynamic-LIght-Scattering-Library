%runMyographDiameter  Measure the vessel diameter of a pressure-myograph recording
%
%   runMyographDiameter(s,fNames) measures the outer, wall-centre and luminal
%   diameter of every *_MYO_d.mat recording in fNames by reading its video, and
%   writes the measured arrays into the recording's own _MYO triplet.  fNames is a
%   cell array of *_MYO_d.mat paths; the wrapper iterates them itself, so there is
%   no launcher for-loop.  The video is found beside the product, or handed over as
%   the optional third argument.
%
%   WHAT IS MEASURED - the analysed span, in this order of precedence:
%     1. PRE-SET INTERVALS.  Windows chosen on the video before any diameter existed
%        (results.intervals, carrying tStart/tEnd): only those are read, and they are
%        measured as separate intervals.
%     2. A TIME CROP (results.timeCrop): only that window is read.
%     3. THE WHOLE FILE, when neither was chosen.  It is then measured as a single
%        interval called 'whole recording' - which is exactly what a recording with
%        no crop and no pre-allocation means.
%   The span is read from the RESULTS, never from a per-type setting: interval times
%   belong to one recording, and no two recordings of a protocol share them.
%
%   A CROP OR A SET OF INTERVALS MAKES THE RECORDING SHORTER, EVERYWHERE.  The
%   diameter is measured in ONE pass over the analysed span ("detect once, then
%   slice"), and source.time is REPLACED by that span: the crop, the pre-set
%   intervals concatenated, or the whole file.  Frames outside it are never read,
%   never stored and never returned - measuring 200 s of a two-hour recording costs
%   200 s of memory.  Two consequences:
%     * TIMES STAY ABSOLUTE seconds from the start of the recording.  They are not
%       re-based to zero, per file or per interval, so an interval can always be
%       located back in the original footage.
%     * WITH DISJOINT INTERVALS source.time HAS GAPS.  It is the intervals
%       concatenated, so the sampling is uniform WITHIN an interval and jumps
%       between them.  Every analysis reads one interval's contiguous .frames slice,
%       so nothing ever runs across a jump - but source.fs is the authority on the
%       sampling rate, and diff(source.time) is not.
%
%   THREE DIAMETERS, ALWAYS.  Each wall is a dark band rather than a line, so
%   getMyographDiameter returns all three of its edges - outer, wall centre and
%   lumen - in the third dimension of source.data, ordered like source.measures.
%   s.edgeMode names the one this step PLOTS and the later steps analyse by default;
%   it no longer selects what is measured.
%
%   LINE AVERAGES ARE THE INTERVAL'S TRACE.  results.intervals(k).diameter holds the
%   line-AVERAGED trace of each measure, averaged over the measured rows only
%   (s.rowRange): the rows outside it were never measured and carry an interpolated
%   fill, and averaging them in would drag the trace towards the edge of the vessel.
%   The full [frames x lines x 3] arrays live ONCE, in source; an interval carries
%   its frame range into them.
%
%   INPUTS
%     s        parameter structure - everything getMyographDiameter takes.  Carried
%              through into settings.runMyographDiameter.  The handful a protocol is
%              actually written around:
%                .rowRange     [lo hi] image rows the vessel occupies
%                .wallContrast minimum contrast for a row to count as a wall
%                .smoothSigma  Gaussian pre-smoothing, px
%                .dustRadius   size of the dark specks removed before detection, px
%                .tSmoothHz    how fast the diameter is allowed to change, Hz
%                .minWallGap   smallest diameter that can be real, px
%              s.intervals / s.intervalNames / s.timeCrop are NOT read: the analysed
%              span comes from the recording's own results (see above).
%     fNames   cell array of *_MYO_d.mat paths.  Empty cells are skipped.
%     fNamesRaw (optional) the matching raw recordings.  Omitted, each recording's
%              own source.fName is used, and then the video sitting beside the
%              product with the same stem.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files and inside the frame loop).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_d.mat   source.time replaced by the analysed span; source.data,
%                        .wallL, .wallR (all [T x nY x 3]), .mask [T x nY],
%                        .valid [T x 1], .measures and .rowRange written
%     <name>_MYO_r.mat   results.intervals(k).frames and .diameter (.time, .outer,
%                        .mid, .inner, .stats.<measure>, .nY)
%     <name>_MYO_s.mat   settings.runMyographDiameter = s
%
%   EXAMPLE
%     s.rowRange     = [40 440];
%     s.wallContrast = 0.05;
%     D = dir(fullfile(rootFolder,'*_MYO_d.mat'));
%     runMyographDiameter(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     getMyographDiameter, myographMeasureIndex, myographProduct (Core/Myograph),
%     Core/Reporting, and MATLAB's base VideoReader.
%
% See also: runMyographVideo, runMyographPropagation, runMyographVasomotion,
%           getMyographDiameter, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Wall detection
% s.rowRange=[1 Inf];     % [first last] image row holding the vessel
% s.wallContrast=0.05;    % minimum contrast for a row to count as a wall
% s.minWallGap=3;         % smallest diameter that can be real, px
% %ADJUSTED IF NECESSARY - Image cleaning
% s.smoothSigma=1.2;      % Gaussian pre-smoothing, px
% s.dustRadius=8;         % dark specks up to this size are removed (0 = off)
% s.tSmoothHz=1;          % the diameter cannot change faster than this, Hz

function runMyographDiameter(s,fNames,fNamesRaw)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_MYO_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_MYO_d.mat".');
end
if nargin<3, fNamesRaw={}; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Diameter',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    fName=char(fNames{fidx});
    reportFile(rep,fidx,fName);
    clearvars source results settings

    [source,results,settings]=myographProduct('open',fName);
    rawName=rawRecording(source,fNamesRaw,fidx,fName);

    % ---- the settings this recording is measured with ----
    sFile=s; sFile.fName=fName;
    if ~isfield(sFile,'rowRange') || isempty(sFile.rowRange)
        sFile.rowRange=fieldOr(source,'rowRange',[1 Inf]);
    end

    % ---- the analysed span: pre-set intervals, else the crop, else the whole file --
    sCore=sFile;
    [sCore.intervals,sCore.intervalNames,sCore.timeCrop,asked]=analysedSpan(results);
    sCore.progressFcn=@(varargin)[];        % cores stay silent: three lines per file

    ivs=getMyographDiameter(sCore,rawName);

    % ---- SOURCE: the analysed span, and the arrays measured over it ---------------
    % getMyographDiameter hands back one struct per interval, each exactly as long as
    % that interval has frames, so the span is a CONCATENATION of what it returned -
    % never a full-length array with holes in it.
    source.time=vertcat(ivs.time);
    source.data=cat(1,ivs.diameter);
    source.wallL=cat(1,ivs.idxL);
    source.wallR=cat(1,ivs.idxR);
    source.mask=logical(cat(1,ivs.mask));
    source.valid=logical(vertcat(ivs.valid));
    source.measures=ivs(1).measures;
    source.rowRange=sCore.rowRange;

    % ---- RESULTS: each interval's frame range into that span, and its traces ------
    rows=measuredRows(sCore.rowRange,size(source.data,2));
    results.intervals=fillIntervals(asked,ivs,rows);

    settings.runMyographDiameter=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,source,results,settings);
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function [ivT,ivNames,tCrop,asked]=analysedSpan(results)
%analysedSpan  WHAT IS MEASURED, in the author's order of precedence: pre-set
%   intervals, else the time crop, else the whole file.  Read from the RESULTS,
%   because interval times are a property of one recording and a crop is a decision
%   made on one recording - neither can live in a per-type setting.
%
%   AN INTERVAL THAT ALREADY CARRIES A DIAMETER IS STILL A WINDOW TO MEASURE.  A
%   re-run with different wall parameters has to measure the SAME windows, or the
%   second run of a protocol would silently analyse something else than the first;
%   the single 'whole recording' interval this step invents spans the whole file, so
%   re-measuring it is the identical pass.  The one case this leaves to its owner is
%   a crop CHANGED after a diameter run: the crop step redefines what the recording
%   is, and clearing the intervals it invalidates is part of setting a crop.
%
%   'asked' IS THE WINDOWS THAT WERE ASKED FOR, in the order they were asked - the
%   same array the times came out of, handed back so the caller cannot re-filter it
%   differently and pair a measured interval with another one's name.
ivT=[]; ivNames={}; tCrop=[]; asked=myographProduct('intervals');
ivs=[];
if isfield(results,'intervals'), ivs=results.intervals; end
if ~isempty(ivs)
    keep=false(1,numel(ivs));
    for k=1:1:numel(ivs)
        keep(k)=~isempty(fieldOr(ivs(k),'tStart')) && ~isempty(fieldOr(ivs(k),'tEnd'));
    end
    if any(keep)
        asked=ivs(keep);
        ivT=[cellfun(@double,{asked.tStart}'), cellfun(@double,{asked.tEnd}')];
        ivNames=cellfun(@char,{asked.name},'UniformOutput',false);
        return
    end
end
if isfield(results,'timeCrop') && numel(results.timeCrop)==2
    tCrop=results.timeCrop;
end
end

% =====================================================================
function out=fillIntervals(asked,ivs,rows)
%fillIntervals  One results.intervals element per measured interval: its frame range
%   into the span source.time now holds, and its line-averaged traces + statistics.
%   WITH NO INTERVALS DEFINED the whole span is a single interval, named for what it
%   is.  A window that WAS asked for keeps its name and its tStart/tEnd - those are
%   the user's definition of it, and this step measures the window rather than
%   redefines it; .frames and .diameter.time carry what was actually measured.
%   'asked' comes straight from analysedSpan, element for element with ivs.
out=myographProduct('intervals');
i0=1;
for k=1:1:numel(ivs)
    n=numel(ivs(k).time);
    if k<=numel(asked)
        out(k).name=asked(k).name;
        out(k).tStart=asked(k).tStart;
        out(k).tEnd=asked(k).tEnd;
        out(k).channels=fieldOr(asked(k),'channels',{});
    else
        nm=char(ivs(k).name);
        if isempty(nm) || strcmp(nm,'default'), nm='whole recording'; end
        out(k).name=nm;
        out(k).channels={};
        if n>0, out(k).tStart=ivs(k).time(1); out(k).tEnd=ivs(k).time(end); end
    end
    out(k).frames=[i0 i0+n-1];
    out(k).diameter=diameterBranch(ivs(k),rows);
    % the two analyses downstream are re-derived from the diameter that has just
    % been measured, so anything an earlier run left behind is now out of date
    out(k).propagation=[];
    out(k).vasomotion=[];
    i0=i0+n;
end
end

% =====================================================================
function d=diameterBranch(iv,rows)
%diameterBranch  The interval's own trace: each measure LINE-AVERAGED over the rows
%   that were really measured, plus the statistics a protocol is compared on.
%   The full per-line arrays are not copied here - they live once, in source.
d=struct('time',[],'nY',0,'stats',struct());
meas=iv.measures;
r=rows(rows<=size(iv.diameter,2));
if isempty(r), r=1:size(iv.diameter,2); end
d.time=iv.time(:);
d.nY=numel(r);
valid=logical(iv.valid(:));
for m=1:1:numel(meas)
    trace=mean(double(iv.diameter(:,r,m)),2,'omitnan');
    d.(meas{m})=trace;
    d.stats.(meas{m})=traceStats(trace,valid,iv.mask(:,r));
end
end

% =====================================================================
function st=traceStats(trace,valid,mask)
%traceStats  What a protocol compares: the level, the swing and how much of it was
%   actually seen.  OFF-FOV FRAMES ARE EXCLUDED from the level and the swing - a
%   wall that dilated out of view leaves an edge-clamped lower bound behind, and
%   averaging it in would report a constriction that never happened.  How many such
%   frames there were is not hidden: it is validFraction.
if isempty(valid) || numel(valid)~=numel(trace), valid=true(size(trace)); end
v=trace(valid); v=v(isfinite(v));
st=struct('mean',NaN,'std',NaN,'min',NaN,'max',NaN,'amplitude',NaN, ...
    'measuredFraction',NaN,'validFraction',NaN);
if ~isempty(v)
    st.mean=mean(v);  st.std=std(v);
    st.min=min(v);    st.max=max(v);
    st.amplitude=st.max-st.min;
end
if ~isempty(mask), st.measuredFraction=mean(double(mask(:))); end
if ~isempty(trace), st.validFraction=mean(double(valid)); end
end

% =====================================================================
function rows=measuredRows(rowRange,nY)
%measuredRows  The image rows the vessel was measured on.  Rows outside them carry
%   an interpolated fill, so nothing that averages along the vessel may include them.
if numel(rowRange)~=2, rowRange=[1 Inf]; end
r0=max(1,round(rowRange(1)));
r1=min(nY,round(rowRange(2)));
if ~(r1>=r0), rows=1:nY; else, rows=r0:r1; end
end

% =====================================================================
function raw=rawRecording(source,fNamesRaw,fidx,fName)
%rawRecording  THE VIDEO this product was made from: the list the caller passed, the
%   path the entry step recorded, or the video sitting beside the product.  The
%   recorded path is tried before the neighbours because it is the only one that
%   knows the original container's extension.
raw='';
if numel(fNamesRaw)>=fidx && ~isempty(fNamesRaw{fidx}), raw=char(fNamesRaw{fidx}); end
if isempty(raw) || ~isfile(raw), raw=char(fieldOr(source,'fName','')); end
if isempty(raw) || ~isfile(raw)
    [fPath,stem]=fileparts(regexprep(fName,'_MYO_d\.mat$',''));
    for ext={'.avi','.mp4','.mov','.mkv'}
        cand=fullfile(fPath,[stem ext{1}]);
        if isfile(cand), raw=cand; break; end
    end
end
if isempty(raw) || ~isfile(raw)
    error('runMyographDiameter:noRecording', ...
        ['The video this data came from was not found next to %s.  ' ...
         'Keep the recording beside its results, or pass the paths as the third argument.'],fName);
end
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
