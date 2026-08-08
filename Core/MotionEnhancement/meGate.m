%meGate - Recover the cardiac gate from the recording itself.
%
%   THE STEP THE WHOLE PACKAGE TURNS ON.  Amplified MRI does not magnify a raw time
%   series; it magnifies one averaged cardiac cycle reconstructed from hundreds of
%   beats, so the noise is already down by the root of that number before the
%   magnifier sees anything.  In MRI the gate is free - it is how a cine is acquired.
%   Here there is no physiological trace at all, so the gate has to come out of the
%   fluorescence intensity itself, and this is where it comes from.
%
%   Session 0 measured what it is worth on this recording, and it is not a nicety.
%   Band-passing a wall's position at the heart rate returns 1.08 to 1.25 times what
%   the identical measurement returns in an empty band - no detection, on any of five
%   vessels, and no better with a longer recording.  Folding the same walls onto
%   cardiac phase over 2 202 beats puts four of the five at 4.6 to 7.6 times their own
%   floor.  Same recording, same frames, same walls: the gate is the difference
%   between a measurement and nothing.
%
%   THE APPROACH AND THE VOCABULARY ARE runContrastInternalCycle's, DELIBERATELY.  This
%   library already solves this problem for LSCI acquisitions - central frequency from
%   the most prominent spectral peak, minima by findpeaks with a distance and a
%   prominence, a min-to-max-to-min pulsesList, nine features per cycle, and outliers
%   rejected against the OTHER cycles - and a second idiom for the same idea would be
%   a cost with no benefit.  What is different is the signal it is driven off
%   (fluorescence intensity rather than blood-flow index) and one rejection rule: rule
%   12 there guards a file read, and here it is the movement bout, which is the thing
%   that actually spoils a cycle in an awake animal.  That wrapper is NOT called -
%   it is .rls- and contrast-specific and reads a hundred gigabytes to get to this
%   point.
%
%   THE ONE STEP THAT COULD NOT BE PORTED IS THE CONDITIONING, AND IT WAS MEASURED
%   BEFORE IT WAS CHANGED.  runContrastInternalCycle smooths with a loess span of
%   floor(meanPointsPerCycle*s.smoothCoef1) and sets the peak prominence at
%   s.minPromCoef times the standard deviation of that smoothed trace.  Both fail
%   here, for reasons that are properties of this recording rather than of the rule:
%   at 99.96 Hz against an 11.32 Hz heart rate there are eight points per cycle, so
%   the span is two, which loess rounds down to one - no smoothing at all - and the
%   standard deviation of a fluorescence trace is dominated by bleaching, respiration
%   and vasomotion rather than by the pulse, so the prominence threshold comes out
%   several times the cardiac excursion.  Run that way on the reference recording it
%   found 1 167 of the 3 168 beats that are there, and 78 % of what it found had a
%   rate outside the accepted band because a "cycle" spanning three beats is three
%   times too long.
%
%   SO THE TRACE IS CONDITIONED WITH A ZERO-PHASE BAND-PASS AT THE RHYTHM'S OWN BAND,
%   PUT BACK ON ITS OWN LEVEL.  xOsc is the rhythm, filtfilt so there is no lag to
%   vary across the field; the level is a running mean over more than one cycle of
%   the slowest accepted rate, taken from the RAW trace so it is a real brightness.
%   Their sum is the trace the cycles and the features are read off, and it plays
%   exactly the part tsBFIF plays in runContrastInternalCycle: noise faster than the rhythm
%   gone, everything slower kept.  The prominence threshold then comes from the
%   rhythm's own amplitude rather than from the whole trace's variance, which is the
%   half of the rule that was doing the damage.  s.smoothCoef1 is NOT carried as a
%   setting nobody reads.
%
%   THE SLOW DECAY COMES OFF BEFORE THE SPECTRUM.  The dye bleaches, and a monotone
%   decay puts power straight into the low end of the search band.  A low-order
%   polynomial in time removes it and nothing near the cardiac or the vasomotion
%   band.  The raw trace is kept for the level and for feature 9 - which is about the
%   brightness the cycle actually sat at.
%
%   REPORT THE ACCEPTANCE RATE; DO NOT TUNE IT UP.  An awake animal's heart rate
%   wanders - 11.11 to 12.50 Hz interquartile on the reference recording, against the
%   steadier anaesthetised traces these coefficients were set on - so more cycles fail
%   here than there.  That is information about the preparation.  Loosening the
%   coefficients to make the number look better throws the information away and puts
%   the bad cycles into the average.
%
% Syntax:
%    gate = meGate(s, x, fs)
%    gate = meGate(s, x, fs, bad)
%
% Inputs:
%    s   - settings.  Fields read:
%          .minFrqIni .maxFrqIni  the band the rhythm's own frequency is looked for
%                                 in.  meEnhance sets these from s.band, so the same
%                                 machinery gates the vasomotion cycles.
%          .rangeFrq              the accepted band is that frequency plus and minus
%                                 this fraction of it.
%          .gateDetrend           order of the polynomial removed before the spectrum.
%          .minPromCoef           minimum peak prominence, as a fraction of the
%                                 rhythm's own standard deviation.
%          .coeffsSTD             one per feature: how many leave-one-out standard
%                                 deviations from the leave-one-out median a feature
%                                 may sit before the cycle is rejected.
%          .coeffsRel             (1) largest min-to-min level jump as a share of the
%                                 cycle's excursion; (2) largest departure of the
%                                 cycle's mean level from the recording's median.
%          .coeffsAbs             half-width, in cycles, of the run a rejection
%                                 spreads over.
%          .excludeFirstNCycles   rejected outright.
%    x   - [frames 1] the intensity trace the rhythm is detected in.
%    fs  - frames per second.
%    bad - optional [frames 1] logical from meEpochs.  A cycle that touches a
%          movement bout is rejected, with its own reason.  Absent means no bout is
%          known, which is not the same as none happening.
%
% Outputs:
%    gate - .f0             the rhythm's central frequency, Hz
%           .band           [min max] Hz a cycle's own rate has to fall in
%           .meanPointsPerCycle
%           .trace          the conditioned trace the cycles were found on
%           .rhythm         its band-passed part alone, which is what the prominence
%                           threshold was taken from
%           .cycles         [n 3] first minimum, maximum, second minimum, in frames
%           .cycleTime      [n 3] the same in seconds, with the two MINIMA refined
%                           to sub-frame precision - that is the axis to fold on
%           .cycleFrq       [n 1] each cycle's own rate, Hz
%           .features       [n 9]
%           .reject         [nRule n] logical
%           .reasons        {nRule 1} what each row means, for the report
%           .rejectCount    [nRule 1] how many cycles each rule caught, first
%                           rule to catch them not being exclusive
%           .accepted       [n 1] logical
%           .acceptedCycles [m 3] the survivors, in frames
%           .acceptedTime   [m 3] the survivors, in refined seconds
%           .acceptanceRate m/n
%           .beatTime       [n+1 1] the minima, refined seconds - the fiducials
%           .phase          [frames 1] cardiac phase in [0,1), NaN outside an
%                           accepted cycle
%
% Example:
%    tr   = load('rec_ME_trace.mat');
%    ep   = meEpochs(meSettings, probe.motion.displacement, fs);
%    gate = meGate(meSettings, tr.traceVes, fs, ep.bad);
%    fprintf('%.2f Hz, %d of %d cycles accepted\n', gate.f0, ...
%        size(gate.acceptedCycles,1), size(gate.cycles,1));
%
% Dependencies: getFFT; Signal Processing Toolbox (findpeaks, butter, filtfilt).
% See also: meEpochs, meCine, meEnhance, meProbe, runContrastInternalCycle,
%           runIntensityInternalCycle
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function gate = meGate(s, x, fs, bad)

x  = double(x(:));
nT = numel(x);
if nargin<4 || isempty(bad), bad = false(nT,1); end
bad = logical(bad(:));
if numel(bad) ~= nT
    error('meGate:bad','The bout flag is %d frames and the trace is %d.',numel(bad),nT);
end

time = (0:nT-1)'./fs;
xd   = x - polyBaseline(x, s.gateDetrend);

% ---- the rhythm's own frequency ---------------------------------------------
% The most prominent peak of a short-window spectrum inside the plausible band -
% runContrastInternalCycle's rule, with its window length, capped at what the trace can
% actually supply.
nShort = min(2.^nextpow2(fs/s.maxFrqIni*50), 2^floor(log2(nT)));
[pw,~,fShort] = getFFT(xd, fs, nShort, 'cpu');
lim = [find(fShort(:)>s.minFrqIni,1,'first'), find(fShort(:)>s.maxFrqIni,1,'first')-1];
if isempty(lim(1)) || isempty(lim(2)) || lim(2)<=lim(1)
    error('meGate:band','No spectral line fits between %g and %g Hz at %g Hz sampling.', ...
        s.minFrqIni, s.maxFrqIni, fs);
end
[~,iPk,~,prom] = findpeaks(pw(lim(1):lim(2)));
if isempty(iPk)
    error('meGate:noPeak','No peak at all between %g and %g Hz; the rhythm is not there.', ...
        s.minFrqIni, s.maxFrqIni);
end
[~,iBest] = max(prom);
f0     = fShort(iPk(iBest)+lim(1)-1);
minFrq = max(f0*(1-s.rangeFrq), s.minFrqIni);
maxFrq = min(f0*(1+s.rangeFrq), s.maxFrqIni);
mppc   = floor(fs/f0);

% ---- condition the trace ----------------------------------------------------
% The rhythm, zero-phase, plus the level it sits on.  See the header for what this
% replaces and for the measurement that justified replacing it.  filtfilt, never
% filter: a causal lag that varies across a field is the propagating wave this block
% must not invent, and a gate is not exempt from that just because it is a scalar.
xOsc = bandpass0(xd, [minFrq maxFrq], fs);
wLev = 2*floor(fs/minFrq)+1;                  % more than one cycle of the slowest rate
xf   = xOsc + movmean(x, wLev);

% ---- the minima, and the cycles between them --------------------------------
% THE MINIMA COME OFF THE RHYTHM ALONE, NOT OFF THE RHYTHM PLUS ITS LEVEL.  A minimum
% is found by its PROMINENCE, which is measured against the surrounding baseline - so
% a beat sitting on a rising level is less prominent than the same beat on a falling
% one, and the threshold then rejects beats for where they happen to sit rather than
% for what they are.  Measured on the reference recording, finding them on the summed
% trace lost 570 of 3 168 beats; finding them on xOsc loses none of consequence.  The
% FEATURES are still read off the summed trace, because a level is what several of
% them are about.
minProm = s.minPromCoef*std(xOsc);
[~,locsMin] = findpeaks(-xOsc,'MinPeakDistance',floor(fs/maxFrq),'MinPeakProminence',minProm);
nC = numel(locsMin)-1;
if nC < 3
    error('meGate:noCycles','Only %d cycles were found in %.1f s of trace.', max(nC,0), nT/fs);
end

cycles = zeros(nC,3);
for i = 1:nC
    seg = xf(locsMin(i):locsMin(i+1));
    [~,iMax] = max(seg);
    cycles(i,:) = [locsMin(i), locsMin(i)+iMax-1, locsMin(i+1)];
end

% ---- and where each minimum really is, between two frames -------------------
% THE FIDUCIAL IS REFINED TO SUB-FRAME PRECISION, AND THAT IS WHAT FILLS THE CINE'S
% PHASE AXIS.  Snap every beat's phase zero to the nearest frame and a frame's phase
% inside its beat is k/L for integer k and L - and with eight or nine frames to a beat
% there are only about sixteen values it can take, so the phase histogram comes out
% with holes in it and the cine's finer bins are filled by interpolating across them.
% A parabola through the three samples around the minimum puts the fiducial where the
% trace actually turns, which dithers the phases off the frame grid and is also
% simply a better answer: the heart does not beat on the frame clock.
tMin = (subSampleMinimum(xOsc, locsMin)-1)./fs;    % a 1-based index; time starts at 0

% ---- nine features per cycle, runContrastInternalCycle's nine ----------------
feat = zeros(nC,9);
for i = 1:nC
    a = cycles(i,1); b = cycles(i,2); c = cycles(i,3);
    three     = xf(cycles(i,:));
    feat(i,1) = mean(three);                                  % level
    feat(i,2) = std(three);
    feat(i,3) = max(three)-min(three);                        % excursion
    feat(i,4) = (tMin(i+1)-tMin(i))*fs;                       % length, frames
    feat(i,5) = abs(xf(a)-xf(c));                             % min-to-min jump
    feat(i,6) = sum(diff(xf(a:b))<0);                         % dips on the upstroke
    feat(i,7) = sum(diff(xf(b:c))>0);                         % rises on the decay
    feat(i,8) = mean(diff(xf(b:c)));                          % slope of the decay
    feat(i,9) = mean(x(a+1:c));                               % the brightness it sat at
end

% ---- one row per rejection rule, one column per cycle ------------------------
nRule  = 15;
reject = false(nRule,nC);
reject(1:9,:) = featureOutliers(feat, s.coeffsSTD);

reject(10,:) = (feat(:,5)./feat(:,3) > s.coeffsRel(1))';
cycleFrq     = fs./feat(:,4);
reject(11,:) = (cycleFrq < minFrq | cycleFrq > maxFrq)';

% Rule 12 is this block's own, and it is where the contrast cycle's file-read guard
% sits.  A cycle the animal moved through is not a cycle of the heart's shape: the
% field is displaced by many pixels, which is past the bound the magnifier is valid
% over, and averaging it in would put that displacement into the cine.
for i = 1:nC
    reject(12,i) = any(bad(cycles(i,1):cycles(i,3)));
end

reject(13,:) = (abs(1-feat(:,1)./median(feat(:,1))) > s.coeffsRel(2))';
reject(14,:) = any(~isfinite(feat),2)';

if s.excludeFirstNCycles > 0
    reject(1:14, 1:min(s.excludeFirstNCycles,nC)) = true;
end

% A single bad cycle in the middle of a good run is usually the detection slipping
% rather than the heart, and its neighbours inherit the doubt.
a = max(reject,[],1);
reject(15,:) = round(movmean(double(a),[s.coeffsAbs(1) s.coeffsAbs(1)])) > 0;

accepted = ~any(reject,1)';

% ---- per-frame cardiac phase -------------------------------------------------
% Only inside an ACCEPTED cycle, because a frame in a rejected one has no phase this
% block is willing to stand behind.  Against the REFINED boundaries, so a frame near
% the start of one beat and a frame near the start of the next are not given the same
% phase merely because both are the first frame after a minimum.
phase = nan(nT,1);
for i = find(accepted)'
    t0 = tMin(i); t1 = tMin(i+1);
    idx = find(time>=t0 & time<t1);
    phase(idx) = (time(idx)-t0)./(t1-t0);
end

gate = struct();
gate.f0                 = f0;
gate.band               = [minFrq maxFrq];
gate.meanPointsPerCycle = mppc;
gate.minProminence      = minProm;
gate.trace              = xf;
gate.rhythm             = xOsc;
gate.time               = time;
gate.cycles             = cycles;
gate.cycleTime          = [tMin(1:end-1), time(cycles(:,2)), tMin(2:end)];
gate.cycleFrq           = cycleFrq;
gate.features           = feat;
gate.reject             = reject;
gate.reasons            = ruleNames();
gate.rejectCount        = sum(reject,2);
gate.accepted           = accepted;
gate.acceptedCycles     = cycles(accepted,:);
gate.acceptedTime       = gate.cycleTime(accepted,:);
gate.acceptanceRate     = mean(accepted);
gate.beatTime           = tMin;
gate.phase              = phase;
gate.rateMedian         = median(cycleFrq(accepted));
gate.rateIQR            = prctile(cycleFrq(accepted),[25 75]);
gate.rateStd            = std(cycleFrq(accepted));
end

% =====================================================================
function p = subSampleMinimum(x, loc)
%subSampleMinimum  Where the trace really turns, in fractional samples.
%   A parabola through the sample at the minimum and its two neighbours.  The vertex
%   of that parabola is the answer for any smooth trace sampled finely enough that
%   three points describe the turn, which at eight samples to a beat is marginal and
%   still far better than snapping to the nearest one.  Clamped to half a sample
%   either way, because a vertex further out than that means the three points do not
%   describe a turn and the sample itself is the best available answer.
loc = loc(:);
p   = double(loc);
in  = loc>1 & loc<numel(x);
a = x(loc(in)-1);  b = x(loc(in));  c = x(loc(in)+1);
den = a - 2*b + c;
d = zeros(size(den));
ok = den ~= 0;
d(ok) = 0.5*(a(ok)-c(ok))./den(ok);
p(in) = p(in) + min(max(d,-0.5),0.5);
end

% =====================================================================
function y = bandpass0(x, band, fs)
%bandpass0  Zero-phase band-pass.  filtfilt, never filter.
w = band./(fs/2);
w = min(max(w,1e-6),0.999);
[bb,aa] = butter(2, w, 'bandpass');
y = filtfilt(bb, aa, double(x(:)));
end

% =====================================================================
function b = polyBaseline(x, order)
%polyBaseline  The slow decay of the dye, as a low-order polynomial in time.
n = numel(x);
t = ((1:n)'-1)./(n-1);
b = polyval(polyfit(t,x,order),t);
end

% =====================================================================
function reject = featureOutliers(feat, coeffsSTD)
%featureOutliers  Rules 1 to 9: a cycle whose feature stands out from the rest.
%   reject(k,i) is true when |f - median(f without i)| > coeffsSTD(k)*std(f without
%   i).  LEAVE-ONE-OUT, so a single wild cycle cannot widen the spread it is then
%   compared against and hide inside it - which is the whole point of the rule and
%   the reason it is not just a distance from the median of everything.
nF = size(feat,2);
nC = size(feat,1);
reject = false(nF,nC);
all1 = 1:nC;
for k = 1:nF
    v = feat(:,k);
    for i = 1:nC
        loo = v(all1~=i);
        m   = median(loo,'omitnan');
        sd  = std(loo,'omitnan');
        reject(k,i) = abs(v(i)-m) > coeffsSTD(k)*sd;   % false for a NaN, which rule 14 has
    end
end
end

% =====================================================================
function r = ruleNames()
%ruleNames  What each row of the rejection matrix means, in a reader's words.
r = {'level unlike the other cycles'
     'variability unlike the other cycles'
     'excursion unlike the other cycles'
     'length unlike the other cycles'
     'start and end levels unlike the other cycles'
     'dips on the upstroke unlike the other cycles'
     'rises on the decay unlike the other cycles'
     'decay slope unlike the other cycles'
     'brightness unlike the other cycles'
     'starts and ends at very different levels'
     'its own rate is outside the accepted band'
     'the animal moved during it'
     'sits at a brightness far from the recording''s'
     'a feature could not be computed'
     'next to other rejected cycles'};
end
%------------- END OF CODE --------------
