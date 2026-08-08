%getBolusCthFloor - The smallest transit-time spread THIS recording can resolve
%
%   THE DELIVERABLE IS NEVER THE MAP - IT IS THE MAP BESIDE THE FLOOR IT HAD TO CLEAR.
%   That is the rule MotionEnhancement was built on (a wall displacement is nothing until
%   it is quoted against its own null), and CTH needs it more than motion did, because
%   the way CTH fails is silent: on a short record the moment subtraction returns a
%   perfectly plausible number that is wrong by a factor of two, or a negative variance,
%   and no property of the measured curve says which.
%
%   IT IS A SYSTEMATIC ERROR, SO NO PER-TRACE STATISTIC CAN CATCH IT.  The variance comes
%   out of m2 = T^2 - 2*INT t*y dt / yT, which on the reference recording is
%   185.2 - 183.2 = 2.0: the answer is one per cent of the terms it is a difference of.
%   A 1 % error in the plateau level moves it by 1.8 s^2 - more than the answer.  The
%   curve can look perfectly settled and the number can still be meaningless, and it is:
%   MEASURED, the reference recording's field mean is flat to 0.45 % per second over its
%   last 3 s and still returns var = -4.9 s^2.
%
%   SO CALIBRATE INSTEAD.  Push kernels of KNOWN spread through the recording's OWN
%   measured input, on its OWN time base, and read them back with the same arithmetic the
%   step uses.  Whatever comes back wrong is not resolvable; the smallest spread that
%   comes back within tolerance is the floor.  It costs one convolution per probe and it
%   turns "CTH is 1.78 s" into "CTH is 1.78 s against a floor of 0.50 s", which is the
%   only form of that sentence anybody can act on.
%
%   THE INPUT IS HELD FLAT BEYOND THE RECORD, which is what a non-clearing tracer does,
%   so the test isolates the RECORD LENGTH and not the extrapolation.
%
%   THE FLOOR MOVES WITH THE TRANSIT DELAY, and the recording is what decides it.  At the
%   delays the reference recording shows (Mtt 0.65-0.83 s) it is 0.5-0.7 s; a preparation
%   with 3 s transits needs a longer span.  What settles Cth for good is 25-30 s of bolus
%   span, which is the entry step's business rather than this one's.
%
% Syntax:
%    F = getBolusCthFloor(t,aif,s)
%
% Inputs:
%    t   - [nT x 1] the recording clock, seconds from 0.
%    aif - [nT x 1] the measured input function, BASELINE-REFERENCED.
%    s   - parameter struct.  Fields:
%            ctthPlateauSec  the plateau window, s (default 1)
%            ctthCthTol      the relative error a probe may carry and still count as
%                            resolved (default 0.25, i.e. 25 %)
%            ctthCthProbes   the probe spreads to try, s (default 0.1:0.1:3)
%
% Outputs:
%    F - struct with .floor (smallest resolved spread, s; Inf when none resolved),
%        .usable (logical), .probeSd, .probeGot, .probeErr and .words - one sentence
%        saying what the recording can and cannot support.
%
% Example:
%    F = getBolusCthFloor(results.time,aif,s);
%    if ~F.usable, fprintf('%s\n',F.words); end
%
% See also: getBolusMoments, getBolusMetrics, getBolusInput, getBolusConfidence, runCTTH
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function F = getBolusCthFloor(t,aif,s)

t=double(t(:)); aif=double(aif(:));
dt=median(diff(t)); nT=numel(t);
wEnd=1; tol=0.25; sdP=0.1:0.1:3;
if isfield(s,'ctthPlateauSec')&&~isempty(s.ctthPlateauSec), wEnd=double(s.ctthPlateauSec); end
if isfield(s,'ctthCthTol')   &&~isempty(s.ctthCthTol),    tol =double(s.ctthCthTol);    end
if isfield(s,'ctthCthProbes')&&~isempty(s.ctthCthProbes), sdP =double(s.ctthCthProbes(:))'; end

% The probe kernels are gamma densities with a fixed mean and the requested spread; the
% mean is a delay and the moment subtraction removes it exactly, so it is the SPREAD the
% probe is about.  A long tail is built to convolve against and then cut back.
kMean=1.0;
nPad=ceil(60/dt);
aifL=[aif; repmat(aif(end),nPad,1)];             % held flat: a non-clearing tracer
tL=(0:dt:dt*(numel(aifL)-1))';
oIn=getBolusMoments(t,aif,wEnd);

F.probeSd=sdP; F.probeGot=nan(size(sdP)); F.probeErr=nan(size(sdP));
for i=1:numel(sdP)
    sh=(kMean/sdP(i))^2; th=kMean/sh;            % gamma with mean kMean, sd sdP(i)
    h=(tL.^(sh-1).*exp(-tL/th))/(gamma(sh)*th^sh);
    h(~isfinite(h))=0; h=h/(sum(h)*dt);
    c=conv(aifL,h)*dt; c=c(1:nT);
    o=getBolusMoments(t,c,wEnd);
    dv=o.var-oIn.var;
    got=sign(dv)*sqrt(abs(dv));
    F.probeGot(i)=got;
    F.probeErr(i)=(got-sdP(i))/sdP(i);
end

ok=abs(F.probeErr)<=tol;
% THE FLOOR IS THE START OF THE FIRST UNBROKEN RUN OF RESOLVED PROBES, not the first
% single one: a lone probe that lands inside tolerance while its neighbours do not is the
% error changing sign, and quoting it as a floor would be quoting a coincidence.
F.floor=Inf; run=0;
for i=1:numel(ok)
    if ok(i), run=run+1; else, run=0; end
    if run>=3, F.floor=sdP(i-2); break; end
end
F.usable=isfinite(F.floor);
if F.usable
    F.words=sprintf(['the recording resolves a transit-time spread down to %.1f s; ' ...
        'anything smaller is below its floor'],F.floor);
else
    F.words=sprintf(['the recording is too short to resolve any transit-time spread ' ...
        '- it holds %.1f s and the width needs the curve to settle first'],t(end)-t(1));
end
end
