%analyzeMyographIntervals  Vasomotion + propagation on already-measured intervals
%
%   [intervals,tim] = analyzeMyographIntervals(s,intervals) runs the vasomotion
%   (getMyographVasomotion) and propagation (getMyographPropagation) analysis on a
%   struct array of intervals whose diameter has ALREADY been measured - no video
%   is read and no diameter is re-detected.  It is the analysis stage the GUI runs
%   after the intervals have been defined and corrected, and the shared core that
%   runMyographFile calls for its vasomotion/propagation stages.
%
%   s.runVasomotion / s.runPropagation (default true) gate the two steps.  An
%   optional stageFcn(stage,detail) receives per-interval progress messages.
%
%   INPUT
%     intervals  struct array with at least time and diameter (mask/idxL/idxR kept).
%   OUTPUT
%     intervals  same array with .vasomotion and/or .prop appended per interval.
%     tim        struct with vasomotion and propagation wall-clock seconds.
%
%   DEPENDS ON  getMyographVasomotion, getMyographPropagation (Core/Myograph/).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

function [intervals,tim] = analyzeMyographIntervals(s,intervals,stageFcn)

if nargin<3, stageFcn=[]; end
doV = ~isfield(s,'runVasomotion')  || isempty(s.runVasomotion)  || s.runVasomotion;
doP = ~isfield(s,'runPropagation') || isempty(s.runPropagation) || s.runPropagation;
tim=struct('vasomotion',0,'propagation',0);
nIv=numel(intervals);

if doV
    t0=tic;
    for iv=1:nIv
        chkCancel(s);
        stage(stageFcn,'vasomotion',sprintf('interval %d/%d (%s)',iv,nIv,intervals(iv).name));
        intervals(iv).vasomotion=getMyographVasomotion(s,intervals(iv).diameter,intervals(iv).time);
    end
    tim.vasomotion=toc(t0);
end
if doP
    t0=tic;
    for iv=1:nIv
        chkCancel(s);
        stage(stageFcn,'propagation',sprintf('interval %d/%d (%s)',iv,nIv,intervals(iv).name));
        intervals(iv).prop=getMyographPropagation(s,intervals(iv));
    end
    tim.propagation=toc(t0);
end
end

% =====================================================================
function chkCancel(s)
if isfield(s,'cancelFcn') && ~isempty(s.cancelFcn) && s.cancelFcn()
    error('analyzeMyographIntervals:cancelled','Analysis stopped by user.');
end
end

% =====================================================================
function stage(fcn,st,msg)
if ~isempty(fcn), try, fcn(st,msg); catch, end, else, fprintf('[%s] %s\n',st,msg); end
end
