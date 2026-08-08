%getBolusInput - Derive the arterial input function from the units that fill first
%
%   Every transit marker in this branch is a DIFFERENCE against the arterial input, so
%   the input is the origin the whole measurement is quoted from.  The author's decision
%   of 07-Aug-2026 is that it is DERIVED AUTOMATICALLY from the earliest-arriving units
%   rather than drawn: measured on the reference recording, the derived input and a
%   hand-picked arterial blob agree to 0.035 s in their centroid and 0.078 s in their
%   50 % time, against an arterial-to-parenchyma delay of 0.5 s - so the automatic answer
%   costs nothing and removes a click from every recording.
%
%   THE AMPLITUDE GATE IS NOT OPTIONAL AND IT IS THE WHOLE FILE'S REASON TO EXIST.  An
%   ungated "earliest arriving" rule selects the DIMMEST units, because a dim unit's 10 %
%   threshold sits inside its own noise and is crossed at random early on.  MEASURED: the
%   ungated mask returned a centroid of 3.224 s - LATER than the whole field's 2.623 s and
%   1.3 s away from the true arterial 1.949 s, i.e. wrong by more than twice the signal.
%   With the gate it returns 1.914 s.  A step that derived its own input without this
%   would reference every marker to noise and look perfectly healthy doing it.
%
%   THE WHOLE FIELD IS NOT AN ACCEPTABLE INPUT and the number says why: its centroid is
%   0.674 s later than the artery's, which is larger than the entire arterial-to-
%   parenchyma delay it would be subtracted from.  Every transit marker would come out
%   negative.
%
%   IT IS TWO PASSES AND HAS TO BE.  The baseline window is derived from an onset, the
%   onset needs a curve, and the curve is what is being chosen - so a provisional pass on
%   the unit mean fixes the windows, and the real selection runs against those.
%
%   IT RUNS ON SEGMENTS, NOT PIXELS, EVEN WHEN THE CALLER WANTS PER-PIXEL MAPS.  The
%   input is a property of the RECORDING, so deriving it once per unit size would be two
%   different origins for one recording and the per-pixel maps would not be comparable
%   with the per-segment table.  runCTTH calls this once, on the segment traces.
%
% Syntax:
%    [aif,info] = getBolusInput(time,data,s)
%
% Inputs:
%    time - [nT x 1] recording clock, seconds from 0.
%    data - [nT x nUnit] the unit traces (segments).
%    s    - parameter struct.  Fields:
%             ctthInputAmpPct  amplitude gate, as a percentile of the units' plateau step
%                              (default 75 - the brighter half of the brighter half)
%             ctthInputPct     share of the GATED units taken, as a percentile of their
%                              arrival (default 5)
%             ctthPlateauSec   the plateau window, s (default 1)
%
% Outputs:
%    aif  - [nT x 1] the input function, BASELINE-REFERENCED.
%    info - struct: .units (logical [1 x nUnit]), .nUnits, .ampGate, .arrGate,
%           .stepRange, .arrRange.
%
% Example:
%    [aif,info] = getBolusInput(results.time, results.sData, s);
%    fprintf('input from %d segments\n', info.nUnits);
%
% See also: getBolusMetrics, getBolusMoments, getBolusCthFloor, getBolusConfidence, runCTTH
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function [aif,info] = getBolusInput(time,data,s)

t=double(time(:)); nT=numel(t); Y=double(data);
if size(Y,1)~=nT
    error('getBolusInput:size','data must be [nT x nUnit] on the given time base.');
end
ampPct=opt(s,'ctthInputAmpPct',75);
arrPct=opt(s,'ctthInputPct'   ,5);
wEnd  =opt(s,'ctthPlateauSec' ,1);

endIdx=t>=t(end)-wEnd;

% ---- pass 1: the unit mean fixes the pre-bolus window ----------------------
mu=mean(Y,2,'omitnan');
mu0=median(mu(t<=t(1)+0.2*(t(end)-t(1))));
muE=mean(mu(endIdx));
k=find(mu-mu0>=0.10*(muE-mu0),1,'first');
if isempty(k), k=max(2,round(0.1*nT)); end
blIdx=t<=t(k)-0.5;
if nnz(blIdx)<5, blIdx=t<=t(max(2,round(0.5*k))); end

% ---- pass 2: per-unit plateau step and 10 % arrival ------------------------
Bl  =mean(Y(blIdx,:),1,'omitnan');
Step=mean(Y(endIdx,:),1,'omitnan')-Bl;
D   =Y-Bl;
arr =nan(1,size(Y,2));
for j=1:size(Y,2)
    if ~(Step(j)>0), continue, end
    i=find(D(:,j)>=0.10*Step(j),1,'first');
    if ~isempty(i), arr(j)=t(i); end
end

% ---- the gate, then the selection -----------------------------------------
ok=isfinite(arr) & Step>0;
if ~any(ok)
    error('getBolusInput:noUnits', ...
        ['No segment of this recording rises above its own pre-injection level, so ' ...
         'there is no arterial curve to derive.  Either the bolus span misses the ' ...
         'injection or the segmentation found nothing vascular.']);
end
ampGate=prctile(Step(ok),ampPct);
gated=ok & Step>=ampGate;
if ~any(gated), gated=ok; end
arrGate=prctile(arr(gated),arrPct);
sel=gated & arr<=arrGate;
if ~any(sel)                                   % a percentile below one unit's worth
    [~,j]=min(arr+~gated*1e6); sel=false(size(ok)); sel(j)=true;
end

aif=mean(D(:,sel),2,'omitnan');
info.units=sel; info.nUnits=nnz(sel);
info.ampGate=ampGate; info.arrGate=arrGate;
info.stepRange=[min(Step(sel)) max(Step(sel))];
info.arrRange =[min(arr(sel))  max(arr(sel)) ];
end

% =====================================================================
function v=opt(s,name,dflt)
v=dflt;
if isfield(s,name)&&~isempty(s.(name)), v=double(s.(name)); end
end
