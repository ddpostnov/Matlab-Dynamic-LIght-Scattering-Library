%analyzeMyographIntervals  Vasomotion + propagation on already-measured intervals
%
%   [intervals,tim] = analyzeMyographIntervals(s,intervals) runs the vasomotion
%   (getMyographVasomotion) and propagation (getMyographPropagation) analysis on a
%   struct array of intervals whose diameter has ALREADY been measured - no video
%   is read and no diameter is re-detected.  It is the analysis stage the GUI runs
%   after the intervals have been defined and corrected, and the shared core that
%   runMyographFile calls for its vasomotion/propagation stages.
%
%   getMyographDiameter measures three diameters (outer wall, wall centre, lumen)
%   and stores them in the 3rd dimension of interval.diameter.  ONE of them is
%   analysed: the one s.edgeMode names, wall centre by default (myographMeasureIndex
%   maps the setting onto the slice).  Intervals from before the three-measure
%   change carry a plain [nFrames x nY] diameter and are analysed as they are.
%
%   s.runVasomotion / s.runPropagation (default true) gate the two steps.  An
%   optional stageFcn(stage,detail) receives per-interval progress messages.
%
%   INPUT
%     intervals  struct array with at least time and diameter (mask/valid/idxL/idxR
%                kept, and mask/valid handed to the propagation estimator).
%   OUTPUT
%     intervals  same array with .vasomotion and/or .prop appended per interval.
%     tim        struct with vasomotion and propagation wall-clock seconds.
%
%   DEPENDS ON  getMyographVasomotion, getMyographPropagation,
%   myographMeasureIndex (Core/Myograph/).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

function [intervals,tim] = analyzeMyographIntervals(s,intervals,stageFcn)

if nargin<3, stageFcn=[]; end
doV = ~isfield(s,'runVasomotion')  || isempty(s.runVasomotion)  || s.runVasomotion;
doP = ~isfield(s,'runPropagation') || isempty(s.runPropagation) || s.runPropagation;
tim=struct('vasomotion',0,'propagation',0);
nIv=numel(intervals);
if nIv==0, return; end

% ---- which of the three diameters is analysed (same choice for every interval) ----
if isfield(s,'edgeMode'), em=s.edgeMode; else, em=''; end
if isfield(intervals,'measures'), meas=intervals(1).measures; else, meas={}; end
kM=myographMeasureIndex(em,meas);

if doV
    t0=tic;
    for iv=1:nIv
        chkCancel(s);
        stage(stageFcn,'vasomotion',sprintf('interval %d/%d (%s)',iv,nIv,intervals(iv).name));
        intervals(iv).vasomotion=getMyographVasomotion(s,pickMeasure(intervals(iv).diameter,kM),intervals(iv).time);
    end
    tim.vasomotion=toc(t0);
end
if doP
    t0=tic;
    for iv=1:nIv
        chkCancel(s);
        stage(stageFcn,'propagation',sprintf('interval %d/%d (%s)',iv,nIv,intervals(iv).name));
        intervals(iv).prop=getMyographPropagation(s,pickMeasure(intervals(iv).diameter,kM), ...
            intervals(iv).time,fieldOr(intervals(iv),'mask'),fieldOr(intervals(iv),'valid'));
    end
    tim.propagation=toc(t0);
end
end

% =====================================================================
function X=pickMeasure(A,k)
%PICKMEASURE  one measure out of a [frames x nY x 3] diameter (a legacy 2-D one is itself)
X=A(:,:,min(k,size(A,3)));
end

% =====================================================================
function v=fieldOr(st,f)
if isfield(st,f), v=st.(f); else, v=[]; end
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
