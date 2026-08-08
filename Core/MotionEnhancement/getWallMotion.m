%getWallMotion - Every cut of one segment, measured and collapsed into one number each.
%
%   THE PER-SEGMENT LAYER OVER meKymograph.  That function finds two walls on ONE
%   cut to a fraction of a pixel; this one runs it over the cuts getSegmentCuts
%   placed along a segment, refuses the cuts that do not cross a vessel, and
%   collapses what is left into a wall amplitude, a centre amplitude and the
%   agreement between the cuts.  It is called once with the product and once with
%   each matched control, and it MUST NOT KNOW WHICH IT HAS - same cuts, same span,
%   same wall fit, same rejections - because any bias meKymograph has then cancels
%   in the ratio between them, and that ratio is the only reason a control means
%   anything.
%
%   ONE CONVENTION, NAMED HERE, FOR EVERY COLUMN DOWNSTREAM.  Every amplitude this
%   returns is PEAK TO PEAK.  In 'cycle' mode that is 2*|c1| of the closed cycle's
%   fundamental, so |c1| alone is the zero-to-peak amplitude and is exactly the
%   quantity every number in Core/MotionEnhancement/README.md is quoted in - these
%   are TWICE those.  In 'band' mode it is twice the equivalent sinusoid's amplitude
%   inside the band, sqrt(2)*std of the band-passed trace, which is the same
%   convention meProbe uses.  A displacement is not a number until somebody says
%   which of the two it is.
%
%   THE FUNDAMENTAL, NOT A MAX MINUS MIN.  A cine's last bin is adjacent to its
%   first, so one Fourier coefficient carries the whole cycle; a peak-to-peak of the
%   samples lets one noisy bin set the answer on its own.
%
%   THE REDUCTION IS A MEAN OF THE CENTRED TRACES, AND IT CHANGES THE ANSWER BY A
%   FACTOR OF TEN.  Several cuts across one vessel see ONE pulse - a pressure wave
%   crosses a hundred micrometres in well under a millisecond against an
%   eighty-eight millisecond cycle - so cuts that agree reinforce and cuts that
%   disagree cancel, and the noise falls as the root of the number of cuts.
%   Averaging each cut's amplitude instead puts a floor under the noise that nothing
%   can lower: the mean of |c1| over pure noise is the noise amplitude, so a segment
%   with no motion would report a confident number, its control would report the
%   same number, and every ratio would be 1.0 for ever.  The mean of the traces is
%   what lets a control fall.  In 'cycle' mode it is EXACTLY the complex mean of the
%   per-cut fundamentals, because the transform is linear, so there is nothing to
%   choose between "average then measure" and "measure then average as vectors".
%
%   EACH CUT IS CENTRED ON ITS OWN MEAN FIRST, because a segment tapers: averaging
%   the raw traces would measure the taper and call it motion.  The resting widths
%   are kept separately - they are the diameter a fractional change is divided by.
%
%   AND THE AGREEMENT TEST COMES OUT FOR FREE.  The coherence is the reduced
%   amplitude over the mean of the per-cut amplitudes: 1 when every cut moves in
%   step, about 1/sqrt(nCut) when they are independent noise.  It is a per-segment
%   number computed from the same traces the amplitude comes from, not a second
%   pass.  It is NOT a formal statistic - the cuts of one vessel are correlated over
%   more than twelve pixels of arc, so its null value is not 1/sqrt(nCut) - and
%   nobody should read it as one.
%
%   NEAR THE FLOOR THE REDUCTION IS BIASED LOW, AND THE COHERENCE SAYS HOW MUCH.
%   Measured against a displacement known in closed form: at a coherence of 0.40 the
%   amplitude is 26 per cent low, at 0.97 it is 1.5 per cent low.  That is a
%   property of a vector mean - |mean(c)| < mean(|c|) whenever the phases scatter -
%   and what scatters them near the floor is noise.  The rows it biases are the rows
%   a confidence gate removes.
%
%   A CUT THAT DOES NOT LOOK LIKE A VESSEL IS REFUSED HERE, NOT AVERAGED AWAY.  A
%   mask, a skeleton and a distance transform will nominate a perpendicular that
%   crosses two vessels, or a place where the segmentation welded two things
%   together, and the wall fit there returns a number rather than an error.  Three
%   tests, all against what the SEGMENTATION already promised, so this takes the
%   library's word for what a vessel is instead of re-deciding:
%     - two walls found in every sample;
%     - the fitted diameter within a factor of the mask's own 2*radius;
%     - the wall's 10-to-90 edge no wider than the vessel it bounds.
%   The third is meProbe's looksLikeAVessel test and its threshold is NOT meProbe's
%   0.5: on this preparation the wall edge is 40 to 68 per cent of the calibre, so
%   0.5 rejects three of the five cuts the whole motion-enhancement package rests on
%   and cut 550 of 987 cuts across the reference field.  The test exists to catch a
%   cut that crosses no vessel at all, and an "edge" as wide as the vessel is that.
%
%   THE REJECTION IS DECIDED ON THE PRODUCT AND IMPOSED ON THE CONTROL.  A control
%   has its timing destroyed, not its picture, so the same cuts would pass on it -
%   but deciding separately would let the two runs average different populations,
%   and then the difference BETWEEN cuts appears in the ratio as detection.  The
%   caller passes the product's own keep vector in when it measures a control.
%
% Syntax:
%    m = getWallMotion(X,cuts,s)
%    m = getWallMotion(X,cuts,s,keep)
%
% Inputs:
%    X    - [rows columns nSample] DOUBLE.  The averaged cardiac cycle, one of its
%           matched controls, or a frame series.  Double on purpose: meKymograph
%           casts whatever it is handed, and a per-cut cast of a whole field would
%           be the entire cost of a sweep over a thousand segments.
%    cuts - getSegmentCuts' output for one segment, on X's own grid.
%    s    - parameter struct.  Reads:
%           .motionMode  'cycle' - X is a closed cycle, the amplitude is its
%                        fundamental.  'band' - X is a frame series, the amplitude
%                        is what sits inside s.bandHz and the floor is what sits in
%                        s.offBandHz.
%           .cutStepPx   sample spacing along a cut, px (0.25).
%           .widthTol    [lo hi] fitted diameter over the mask's, accepted range.
%           .edgeTol     10-to-90 wall edge, as a share of the diameter.
%           .bandHz .offBandHz .fps   'band' mode only.
%    keep - optional [nCut 1] logical to impose instead of deciding.
%
% Outputs:
%    m - .nCut         cuts that carried a measurement
%        .keep         [nCut 1] logical
%        .reject       [nCut 1] uint8 - 0 kept, 1 off the frame, 2 a sample with no
%                      wall, 3 not the width the mask promised, 4 an edge, not a wall
%        .wallAmp      peak-to-peak of the reduced half-width, px of X
%        .centreAmp    the same for the reduced centre, px of X
%        .wallNullAmp  'band' mode: the same measurement in the off band.  NaN in
%                      'cycle' mode, where the floor is a separate array
%        .centreNullAmp  likewise
%        .wallCohere .centreCohere   do the cuts agree, 0 to 1
%        .wallAmpCut .centreAmpCut   [nKept 1] per-cut peak-to-peak, for the spread
%        .wallTrace .centreTrace     [nSample 1] the reduced centred traces
%        .radiusPx     median resting half-width over the kept cuts, px of X
%        .radiusSdPx   its spread over the cuts - the segment's taper
%        .edgePx       median 10-to-90 wall edge width, px of X
%
% Example:
%    m  = getWallMotion(double(source.data), cuts, s);
%    mC = getWallMotion(double(source.control(:,:,:,1)), cuts, s, m.keep);
%    fprintf('%.4f px over %d cuts, x%.1f its control\n', ...
%        m.wallAmp, m.nCut, m.wallAmp/mC.wallAmp);
%
% Dependencies: meKymograph; Signal Processing Toolbox (butter, filtfilt) in 'band'
%               mode only.
% See also: getSegmentCuts, runMotionEnhancement, meKymograph, meCine, meProbe
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function m = getWallMotion(X,cuts,s,keep)

nCut    = numel(cuts);
nSample = size(X,3);
decide  = nargin<4 || isempty(keep);
if decide, keep = false(nCut,1); else, keep = logical(keep(:)); end

m = struct('nCut',0,'keep',keep,'reject',zeros(nCut,1,'uint8'), ...
    'wallAmp',NaN,'centreAmp',NaN,'wallNullAmp',NaN,'centreNullAmp',NaN, ...
    'wallCohere',NaN,'centreCohere',NaN,'wallAmpCut',[],'centreAmpCut',[], ...
    'wallTrace',[],'centreTrace',[],'radiusPx',NaN,'radiusSdPx',NaN,'edgePx',NaN);
if nCut==0, return; end

H  = nan(nSample,nCut);
C  = nan(nSample,nCut);
rP = nan(nCut,1);
eP = nan(nCut,1);
[nr,nc,~] = size(X);

for k = 1:nCut
    if ~decide && ~keep(k), continue; end
    c = cuts(k);
    c.step = s.cutStepPx;
    if c.center(1)-c.span<1 || c.center(1)+c.span>nr || ...
       c.center(2)-c.span<1 || c.center(2)+c.span>nc
        keep(k)=false;  m.reject(k)=1;  continue          % dropped, never clipped
    end

    kymo = meKymograph(X,c);
    hw   = kymo.halfWidth;
    ce   = kymo.centre;
    rP(k) = mean(hw,'omitnan');
    eP(k) = mean([kymo.widthL; kymo.widthR],'omitnan');

    if decide
        if ~all(isfinite(hw)) || ~all(isfinite(ce)) || ~isfinite(eP(k))
            m.reject(k)=2;  continue
        end
        rel = rP(k)/c.radius;
        if rel<s.widthTol(1) || rel>s.widthTol(2)
            m.reject(k)=3;  continue
        end
        if eP(k) > s.edgeTol*2*rP(k)
            m.reject(k)=4;  continue
        end
        keep(k) = true;
    end

    H(:,k) = hw;
    C(:,k) = ce;
end

m.keep = keep;
if ~any(keep), return; end

% ---- the reduction: every cut on its own mean, then one mean of the traces ----
Hk = H(:,keep) - mean(H(:,keep),1,'omitnan');
Ck = C(:,keep) - mean(C(:,keep),1,'omitnan');
m.nCut       = nnz(keep);
m.wallTrace  = mean(Hk,2,'omitnan');
m.centreTrace= mean(Ck,2,'omitnan');

% ---- and the amplitude convention, applied identically to both ---------------
switch lower(s.motionMode)
    case 'cycle'
        m.wallAmp     = cycleAmp(m.wallTrace);
        m.centreAmp   = cycleAmp(m.centreTrace);
        m.wallAmpCut  = arrayfun(@(j) cycleAmp(Hk(:,j)), (1:m.nCut).');
        m.centreAmpCut= arrayfun(@(j) cycleAmp(Ck(:,j)), (1:m.nCut).');
    case 'band'
        m.wallAmp      = bandAmp(m.wallTrace,  s.bandHz,   s.fps);
        m.centreAmp    = bandAmp(m.centreTrace,s.bandHz,   s.fps);
        m.wallNullAmp  = bandAmp(m.wallTrace,  s.offBandHz,s.fps);
        m.centreNullAmp= bandAmp(m.centreTrace,s.offBandHz,s.fps);
        m.wallAmpCut   = arrayfun(@(j) bandAmp(Hk(:,j),s.bandHz,s.fps), (1:m.nCut).');
        m.centreAmpCut = arrayfun(@(j) bandAmp(Ck(:,j),s.bandHz,s.fps), (1:m.nCut).');
    otherwise
        error('getWallMotion:motionMode', ...
            's.motionMode is ''cycle'' or ''band''; it is ''%s''.', char(s.motionMode));
end

% The reduced amplitude over the mean of the per-cut ones.  In 'cycle' mode this is
% identically the resultant length |sum(c1)|/sum(|c1|) of the per-cut fundamentals.
m.wallCohere   = m.wallAmp  /max(mean(m.wallAmpCut,  'omitnan'),eps);
m.centreCohere = m.centreAmp/max(mean(m.centreAmpCut,'omitnan'),eps);

m.radiusPx   = median(rP(keep),'omitnan');
m.radiusSdPx = std(rP(keep),'omitnan');
m.edgePx     = median(eP(keep),'omitnan');
end

% =====================================================================
function a = cycleAmp(v)
%cycleAmp  Peak to peak of a closed cycle's fundamental.
%   a = 2*|c1| with c1 = (2/n)*sum(v.*exp(-2i*pi*k/n)), which is meRun's own
%   meCycleAmp doubled.  A trace with a quarter of its samples missing is not a
%   cycle and returns NaN rather than a number nobody can trust.
v  = double(v(:));
ok = isfinite(v);
if nnz(ok) < 0.75*numel(v), a = NaN; return; end
v(~ok) = mean(v(ok));
n = numel(v);
a = 4*abs(sum(v.*exp(-2i*pi*(0:n-1).'./n))./n);
end

% =====================================================================
function a = bandAmp(x,band,fs)
%bandAmp  Peak to peak of the equivalent sinusoid inside a band, zero phase.
%   2*sqrt(2)*std of the band-passed trace, so a pure sinusoid of peak-to-peak A
%   reports A.  Gaps left by a sample whose walls could not be fitted are filled by
%   interpolation rather than dropped, because a band-pass has no idea what a gap is.
x  = double(x(:));
ok = isfinite(x);
if nnz(ok) < 0.5*numel(x), a = NaN; return; end
x(~ok) = interp1(find(ok),x(ok),find(~ok),'linear','extrap');
a = 2*sqrt(2)*std(bandPass0(x,band,fs));
end

% =====================================================================
function y = bandPass0(x,band,fs)
%bandPass0  Zero-phase band-pass.
%   filtfilt, never filter: a causal lag that varies across a field is a propagating
%   wave this measurement must not invent.
%
%   THIS IS THE FOURTH COPY OF TEN LINES, and that is worth a sentence rather than a
%   shrug.  meGate, meProbe and runIntensityInternalCycle each carry one as a
%   subfunction.  Extracting one shared core is the cheap half of the open question
%   about how many cardiac gates this library should have (Q16), and it was left
%   alone here for one reason: two of the three copies sit inside files that are
%   pinned BY VALUE, so moving them is a re-pin and a session of its own.
w = band./(fs/2);
w = min(max(w,1e-6),0.999);
[bb,aa] = butter(2,w,'bandpass');
y = filtfilt(bb,aa,double(x(:)));
end
%------------- END OF CODE --------------
