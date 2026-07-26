%loadMyographResults  Load a *_myograph.mat result and migrate legacy files on the fly
%
%   [intervals,s,meta] = loadMyographResults(fName) loads a myograph result file
%   written by the launcher or runMyographFile and returns its per-interval
%   struct array, the parameter struct s and a meta struct.  Older files (written
%   before the meta/versioning was introduced) hold only intervals and s; they
%   still load - a meta struct is synthesised (frame rate recovered from the
%   interval time base, pixel size left empty/uncalibrated) and every interval is
%   given the full field set (missing vasomotion / prop left empty), so downstream code
%   (Explore, Excel export) can treat old and new files uniformly.
%
%   INPUT
%     fName   path to a *_myograph.mat file.
%   OUTPUT
%     intervals  struct array (name,time,idxL,idxR,diameter,mask,vasomotion,prop)
%     s          parameter struct that produced it ({} if absent)
%     meta       struct: formatVersion, fName, frameRate, pixelSize (µm/px; []
%                or 0 = uncalibrated), timeCrop, rowRange, codeVersion,
%                createdTimestamp (migratedFrom set when an old file was upgraded)
%
%   DEPENDS ON  base MATLAB only.  Self-contained.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

function [intervals,s,meta] = loadMyographResults(fName)

CURRENT_VERSION=2;                                   % current results-format version

% A result written by an older build may contain a parameter struct s whose GUI
% callbacks (progressFcn/cancelFcn/stageFcn) captured the workbench uifigure in their
% closures; load() then deserialises that figure as a stray "ghost" window (with
% "cannot load ... listener" warnings).  Snapshot the open figures, silence those
% warnings during load, then delete any window the load spawned.  New files are
% written handle-free (see stripFcnHandles) so this only ever fires on legacy files.
figsBefore=findall(groot,'Type','figure');
wState=warning('off','all'); wRestore=onCleanup(@()warning(wState)); % silence legacy load-listener warnings
S=load(char(fName));
clear wRestore;                                     % restore warnings immediately after the load
ghosts=setdiff(findall(groot,'Type','figure'),figsBefore);
if ~isempty(ghosts), delete(ghosts); end

if ~isfield(S,'intervals')
    error('loadMyographResults:noIntervals','%s has no ''intervals'' variable.',fName);
end
intervals=S.intervals;
if isfield(S,'s'),    s=stripFcnHandles(S.s); else, s=struct(); end
if isfield(S,'meta'), meta=S.meta; else, meta=struct(); end

meta=migrateMeta(meta,intervals,s,fName,CURRENT_VERSION);
intervals=ensureIntervalFields(intervals);
end

% =====================================================================
function meta=migrateMeta(meta,intervals,s,fName,current)
%MIGRATEMETA  fill / upgrade the meta struct for a legacy or partial file
if ~isfield(meta,'formatVersion') || isempty(meta.formatVersion)
    meta.formatVersion=1;                            % pre-versioning file
end
from=meta.formatVersion;
if ~isfield(meta,'fName')    || isempty(meta.fName),    meta.fName=char(fName); end
if ~isfield(meta,'frameRate')|| isempty(meta.frameRate)
    meta.frameRate=frameRateFromIntervals(intervals);
end
if ~isfield(meta,'pixelSize'),   meta.pixelSize=[]; end          % []/0 -> uncalibrated (px)
if ~isfield(meta,'timeCrop'),    meta.timeCrop=[]; end
if ~isfield(meta,'rowRange')
    if isfield(s,'rowRange') && ~isempty(s.rowRange), meta.rowRange=s.rowRange; else, meta.rowRange=[1 Inf]; end
end
if ~isfield(meta,'codeVersion'),      meta.codeVersion='legacy'; end
if ~isfield(meta,'createdTimestamp'), meta.createdTimestamp=''; end
if from<current
    meta.migratedFrom=from;
    meta.formatVersion=current;
end
end

% =====================================================================
function fr=frameRateFromIntervals(intervals)
fr=NaN;
if ~isempty(intervals) && isfield(intervals,'time') && numel(intervals(1).time)>1
    fr=1/median(diff(double(intervals(1).time)),'omitnan');
end
end

% =====================================================================
function intervals=ensureIntervalFields(intervals)
%ENSUREINTERVALFIELDS  guarantee the full field set so old files behave like new ones
need={'name','time','idxL','idxR','diameter','mask','valid','vasomotion','prop'};
for k=1:numel(need)
    if ~isfield(intervals,need{k})
        [intervals.(need{k})]=deal([]);
    end
end
% legacy files predate the off-FOV 'valid' flag: treat every measured frame as valid
for j=1:numel(intervals)
    if isempty(intervals(j).valid) && ~isempty(intervals(j).diameter)
        intervals(j).valid=true(size(intervals(j).diameter,1),1);
    end
end
end
