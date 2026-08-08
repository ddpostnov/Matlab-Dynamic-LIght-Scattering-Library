%meValidate - The calibration curve.  What comes out, for a known input and a known factor.
%
%   THE DELIVERABLE OF THE ENGINE, AND IT IS NOT THE MOVIE.  Motion magnification
%   always produces a moving picture: give it noise and a large enough factor and it
%   produces convincing wall motion from a recording that contains none.  So the
%   thing that makes this a library section rather than a demonstration is the
%   number beside the movie, and this is where the number is made - on a phantom
%   whose displacement is a constant in a formula rather than another estimate.
%
%   FOUR THINGS ARE MEASURED, AND THE LAST ONE IS THE ONE THAT MATTERS:
%
%   1. THE GAIN, against 1+alpha.  Wadhwa's Eq. 5 says the output has moved by
%      (1+alpha)*delta.  What actually comes out of a binomial Laplacian pyramid is
%      1 + eta*alpha, with eta below one and measured here per level, because the
%      pyramid's bands overlap and a band that is not amplified holds the edge back.
%      eta is a property of the representation, not a fitted correction: it is
%      constant in alpha, it is additive across levels, and it is the reason a
%      recommended alpha cannot be read off the theory alone.
%
%   2. THE SATURATION POINT, where alpha*delta approaches lambda/4 and the feature
%      breaks up.  On this recording it is nowhere near binding - Session 0 measured
%      wall displacements of 0.006 to 0.075 px, so the ceiling runs into the
%      hundreds - but it is where the curve ends and a future recording may live
%      there.
%
%   3. WHAT THE NOISE FLOOR IS, in the same pixels, and how it moves with the amount
%      of averaging that happened BEFORE the magnifier.  That is the axis the whole
%      package turns on: Session 0 found the wall motion invisible to any
%      single-beat measurement and visible only after folding 2 202 beats, so the
%      floor has to be reported against the averaging, not as one number.
%
%   4. WHETHER AMPLIFICATION IMPROVES DETECTION AT ALL, which is the question a
%      reader will ask and which the phantom can answer outright.  Both the signal
%      and the phase noise are multiplied by alpha, so the ratio between them is
%      expected to be flat in alpha - magnification makes a motion VISIBLE, it does
%      not make it DETECTABLE.  If that is what comes out, then averaging before the
%      magnifier is not one option among several, it is the only lever there is.
%
%   THE ESTIMATOR IS A PROJECTION ONTO THE KNOWN FREQUENCY, and its null is the same
%   projection at frequencies nothing occupies.  That is Session 0's method and it
%   is used here for the same reason: an amplitude read off a band-pass has a floor
%   of its own, and a measurement that does not carry its own null is not a
%   measurement.  A second null is run beside it - the identical pipeline on a
%   phantom whose amplitude is zero, control C1 - and the two should agree.
%
%   NOTHING HERE IS TUNED TO MAKE THE PHANTOM LOOK GOOD.  The curve IS the product.
%
% Syntax:
%    cal = meValidate(s, outFolder)
%
% Inputs:
%    s         - settings from meSettings.  The phantom, the magnifier and the
%                sweep fields are all read; see meSettings.
%    outFolder - where the report page and <stem>_ME_calibration.mat go.
%
% Outputs:
%    cal - .gain      the (amplitude x factor x level set) grid of measured gains
%          .eta       the share of the excess each level set delivers
%          .satAlpha  the factor at which the gain has fallen 10 per cent below its
%                     own linear trend, per amplitude - the saturation point
%          .satShift  the same as a commanded shift alpha*delta, which is what the
%                     lambda/4 bound is stated on
%          .gainBlur  gain against the phase-blur width, per level set
%          .signal .null .nullOff .detect
%                     the floor sweep, indexed (averaging, blur, level set, alpha).
%                     null is control C1 - the pipeline on a phantom that does not
%                     move - and nullOff is the signal run's own passband, which is
%                     Session 0's method.  Both are averaged over every bin of the
%                     passband, because one bin of a noise process scatters by its
%                     own size.
%          .floorPx .floorUm  the detection floor, same indexing
%          .riesz     the three-tap transform beside the exact one
%          .linear    the same numbers for meLinear, the noise-magnifying control
%          .settings  the settings this ran with
%
% Example:
%    s   = meSettings;
%    cal = meValidate(s, 'C:\work\motionEnhancement');
%
% Dependencies: mePhantom, meShift, meLinear, meKymograph, meSettings; the
%               reporting module (reportOpen, reportFile, reportFigure, reportSave,
%               reportClose).
% See also: mePhantom, meShift, meKymograph, meLinear, meSettings
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function cal = meValidate(s, outFolder)

if ~isfolder(outFolder), mkdir(outFolder); end
stem  = 'mePhantom';
fName = fullfile(outFolder,[stem '.mat']);

s.rootFolder    = outFolder;
s.resultsFolder = outFolder;
rep = reportOpen(s,'Motion enhancement calibration',{fName});
reportFile(rep,1,fName);

nMax   = maxLevels(s.phSize);
sets   = { 3, 4, 5, 3:5, 1:nMax };
setNm  = {'3','4','5','3:5','all'};
alphas = s.calAlpha(:)';
amps   = s.calAmp(:)';

cal            = struct();
cal.settings   = s;
cal.levelSets  = setNm;
cal.alpha      = alphas;
cal.amplitude  = amps;
cal.pixelSize  = s.pixelSize;

% ===== 1. the gain grid, with the photon noise turned down ===================
% The gain is a property of the representation, so it is measured where noise is
% not competing with it.  The floor is measured separately, at the real one.
fprintf('  gain grid: %d amplitudes x %d factors x %d level sets\n', ...
    numel(amps), numel(alphas), numel(sets));
gain = nan(numel(amps), numel(alphas), numel(sets));
raw  = nan(numel(amps),1);
for ia = 1:numel(amps)
    sp = s;  sp.phAmp = amps(ia);  sp.phPhotons = s.phPhotons*s.calCleanGain;
    [stack, truth] = mePhantom(sp);
    sp.passband = truth.passband;  sp.nPyrLevels = nMax;  sp.alpha = 0;  sp.levels = 1:nMax;
    [~, an] = meShift(stack, truth.fs, sp);
    raw(ia) = amplitudeOf(stack, truth);
    for is = 1:numel(sets)
        sp.levels = sets{is};
        for ip = 1:numel(alphas)
            sp.alpha = alphas(ip);
            gain(ia,ip,is) = amplitudeOf(meShift(an,sp), truth)/raw(ia);
        end
    end
    fprintf('    A = %.3f px: raw %.4f px (%.3f of truth), gain at alpha %d = %.2f\n', ...
        amps(ia), raw(ia), raw(ia)/amps(ia), alphas(end), gain(ia,end,4));
end
cal.gain    = gain;
cal.rawAmp  = raw;

% eta is read only where nothing can have saturated: the commanded shift alpha*A
% has to stay well inside a quarter of the shortest wavelength amplified, which is
% 16 px at level 3.  Two pixels is half of that, with room.
cal.etaPerAmp = nan(numel(amps), numel(sets));
for ia = 1:numel(amps)
    ok = alphas>0 & alphas.*amps(ia) <= 2;
    if ~any(ok), continue; end
    cal.etaPerAmp(ia,:) = squeeze(mean((gain(ia,ok,:)-1)./reshape(alphas(ok),1,[]),2,'omitnan'));
end
cal.eta = median(cal.etaPerAmp, 1, 'omitnan');

% Saturation: the factor at which the gain has fallen a tenth below the straight
% line eta predicts.  The trend uses the level set's own eta, taken across every
% amplitude, so a row with no un-saturated point of its own still has a reference.
cal.satAlpha = nan(numel(amps), numel(sets));
for ia = 1:numel(amps)
    for is = 1:numel(sets)
        pred = 1 + cal.eta(is).*alphas;
        bad  = find(gain(ia,:,is) < 0.9*pred & alphas>0, 1, 'first');
        if ~isempty(bad), cal.satAlpha(ia,is) = alphas(bad); end
    end
end
% The commanded shift at which it went, which is the quantity the bound is stated
% on.  If saturation is a property of the representation rather than of this
% phantom, this column is constant down each level set and can be read against
% lambda/4 directly.
cal.satShift = cal.satAlpha.*amps(:);

% ===== 2. what the blur and the level choice cost in gain ====================
% Both are chosen for noise, and both give up amplification to get it: a level that
% is not amplified holds the edge back, and a blur of the phase spreads the shift
% away from the wall that asked for it.  This panel and the floor below have to be
% read together, because neither is a decision on its own.
blurs  = s.calBlur(:)';
setsB  = {s.levels, 1:nMax};
setNmB = {sprintf('%d:%d',s.levels(1),s.levels(end)), sprintf('1:%d',nMax)};
fprintf('  blur: %d widths x %d level sets, on the clean phantom\n', numel(blurs), numel(setsB));
sp = s;  sp.phAmp = s.calFloorAmp;  sp.phPhotons = s.phPhotons*s.calCleanGain;
sp.nPyrLevels = nMax;
[stkC, truthC] = mePhantom(sp);
sp.passband = truthC.passband;
rawC = amplitudeOf(stkC, truthC);
gainBlur = nan(numel(blurs), numel(setsB));
for ib = 1:numel(blurs)
    for is = 1:numel(setsB)
        sp.blurSigma = blurs(ib);  sp.levels = setsB{is};  sp.alpha = 0;
        [~, anB] = meShift(stkC, truthC.fs, sp);
        sp.alpha = 20;
        gainBlur(ib,is) = amplitudeOf(meShift(anB,sp), truthC)/rawC;
    end
    fprintf('    blurSigma %.1f: gain at alpha 20 = %s\n', blurs(ib), ...
        mat2str(round(gainBlur(ib,:),2)));
end
cal.blur     = blurs;
cal.blurSets = setNmB;
cal.gainBlur = gainBlur;
cal.etaBlur  = (gainBlur-1)./20;

% ===== 3. the floor, against averaging, blur and levels ======================
% calSnrGain is a multiplier on the photon count, i.e. on the number of frames
% averaged into each frame the magnifier sees.  1 is one raw frame; the cine's own
% value comes from Session 0's beat count and is in meSettings.
nG = numel(s.calSnrGain);  nB = numel(blurs);  nS = numel(setsB);  nA = numel(alphas);
fprintf('  floor: %d averaging gains x %d blurs x %d level sets x %d factors\n', nG,nB,nS,nA);
sig = nan(nG,nB,nS,nA);  nul = sig;  nulOff = sig;
iB  = find(blurs==s.blurSigma,1);  if isempty(iB), iB = 1; end
for ig = 1:nG
    g  = s.calSnrGain(ig);
    sp = s;  sp.phPhotons = s.phPhotons*g;  sp.nPyrLevels = nMax;
    sp.phAmp = s.calFloorAmp;   [stk , truth] = mePhantom(sp);
    sp.phAmp = 0;                stk0         = mePhantom(sp);
    sp.passband  = truth.passband;
    if g == 1, stk1 = stk; stk1n = stk0; truth1 = truth; end
    for ib = 1:nB
        for is = 1:nS
            sp.blurSigma = blurs(ib);  sp.levels = setsB{is};  sp.alpha = 0;
            [~, anS] = meShift(stk , truth.fs, sp);
            [~, anN] = meShift(stk0, truth.fs, sp);
            for ip = 1:nA
                sp.alpha = alphas(ip);
                % Two independent floors: control C1 - the identical pipeline on a
                % phantom that does not move - and the signal run's own passband,
                % which is Session 0's method and needs no second recording.  They
                % measure the same thing by different routes and should agree.
                [a,aOff]            = amplitudeOf(meShift(anS,sp), truth);
                sig(ig,ib,is,ip)    = abs(a);
                nulOff(ig,ib,is,ip) = aOff;
                [~,c1]              = amplitudeOf(meShift(anN,sp), truth);
                nul(ig,ib,is,ip)    = c1;
            end
        end
    end
    fprintf('    x%-4g averaged (snr %5.0f/px): blur %.1f, levels %s, alpha 20 -> signal %.4f px, null %.4f px, ratio %.1f\n', ...
        g, sqrt(s.phPhotons*g), blurs(iB), setNmB{1}, ...
        sig(ig,iB,1,alphas==20), nul(ig,iB,1,alphas==20), ...
        sig(ig,iB,1,alphas==20)/nul(ig,iB,1,alphas==20));
end
cal.snrGain = s.calSnrGain(:)';
cal.signal  = sig;
cal.null    = nul;
cal.nullOff = nulOff;
cal.detect  = sig./nul;

% The floor is where the signal would equal twice its own null, and the response is
% linear in the displacement, so it scales straight off the point measured.
cal.floorPx = s.calFloorAmp.*2.*nul./sig;
cal.floorUm = cal.floorPx.*s.pixelSize;

% ===== 4. the three taps against the exact transform =========================
% The phantom does not depend on the transform, so the stacks already built serve.
fprintf('  Riesz: 3tap against exact\n');
cal.riesz = struct('mode',{'3tap','exact'},'gain',[],'floorPx',[],'seconds',[]);
for im = 1:2
    sp = s;  sp.riesz = cal.riesz(im).mode;  sp.nPyrLevels = nMax;
    sp.levels = s.levels;  sp.passband = truthC.passband;  sp.alpha = 0;
    tic; [~, anC] = meShift(stkC, truthC.fs, sp); cal.riesz(im).seconds = toc;
    sp.alpha = 20;
    cal.riesz(im).gain = amplitudeOf(meShift(anC,sp), truthC)/rawC;

    sp.alpha = 0;
    [~, anS] = meShift(stk1 , truth1.fs, sp);
    [~, anN] = meShift(stk1n, truth1.fs, sp);
    sp.alpha = 20;
    a      = abs(amplitudeOf(meShift(anS,sp), truth1));
    [~, n] = amplitudeOf(meShift(anN,sp), truth1);
    cal.riesz(im).floorPx = s.calFloorAmp*2*n/a;
    fprintf('    %-6s gain %.2f at alpha 20, floor %.4f px, analysis %.2f s\n', ...
        cal.riesz(im).mode, cal.riesz(im).gain, cal.riesz(im).floorPx, cal.riesz(im).seconds);
end

% ===== 5. the linear control =================================================
sp = s;  sp.nPyrLevels = nMax;  sp.levels = s.levels;  sp.passband = truth1.passband;
sp.alpha = 20;
aL      = abs(amplitudeOf(meLinear(stk1 , truth1.fs, sp), truth1));
[~, nL] = amplitudeOf(meLinear(stk1n, truth1.fs, sp), truth1);
cal.linear = struct('alpha',20,'gain',aL/abs(amplitudeOf(stk1,truth1)), ...
                    'floorPx',s.calFloorAmp*2*nL/aL);
fprintf('  meLinear at alpha 20: apparent gain %.2f, floor %.4f px (phase-based %.4f)\n', ...
    cal.linear.gain, cal.linear.floorPx, cal.floorPx(1,iB,1,alphas==20));

% ===== the page ==============================================================
drawCalibration(reportFigure(rep,'calibration'), cal, s, rep);
save(fullfile(outFolder,[stem '_ME_calibration.mat']),'cal','-v7');
reportClose(rep);
end

% =====================================================================
function [a, aFloor] = amplitudeOf(stack, truth)
%amplitudeOf  The displacement amplitude at the phantom's own frequency, signed,
%   and the same estimator's own floor.
%
%   THE FLOOR IS AVERAGED OVER EVERY BIN OF THE PASSBAND, not read off one.  A
%   single Fourier coefficient of a noise process is Rayleigh distributed: its
%   magnitude scatters by about its own size, so a floor taken from one bin carries
%   a factor of two and reports a detection limit that moves when nothing has
%   changed.  Averaging the thirty-odd bins the filter passes cuts that scatter by
%   root thirty and is what makes a floor comparable between two settings.
%
%   AND IT IS TAKEN INSIDE THE PASSBAND, because that is where the magnified phase
%   noise is.  Bins outside it carry only the raw edge-fit noise, which is a
%   different quantity and a smaller one.
k = meKymograph(stack, truth.cut);
x = k.(truth.measure);
a = projectAt(x, truth.freq, truth.fs);
if nargout>1
    df  = truth.fs/numel(x);
    kk  = round(truth.passband(1)/df):round(truth.passband(2)/df);
    kk  = kk(kk>=1 & abs(kk*df - truth.freq) > 1.5*df);
    v   = arrayfun(@(g) abs(projectAt(x,g,truth.fs)), kk.*df);
    aFloor = sqrt(mean(v.^2));
end
end

% =====================================================================
function a = projectAt(x, f, fs)
%projectAt  One Fourier coefficient, signed against sin(2*pi*f*t).  The phantom's
%   frequency sits on a bin exactly, so this leaks nothing.
x  = double(x(:));
ok = isfinite(x);
if ~all(ok), x(~ok) = interp1(find(ok),x(ok),find(~ok),'linear','extrap'); end
n = numel(x);  t = (0:n-1)'./fs;
a = -imag(2/n*sum((x-mean(x)).*exp(-1i*2*pi*f*t)));
end

% =====================================================================
function n = maxLevels(sz)
n = max(1, floor(log2(min(sz))) - 2);
end

% =====================================================================
function drawCalibration(fh, cal, s, rep)
%drawCalibration  One page: the curve, the saturation, the floor, and whether any
%   of it helps detection.
t = tiledlayout(fh,2,3,'Padding','compact','TileSpacing','compact');
al = cal.alpha;  amps = cal.amplitude;
iMid = find(amps>=0.05,1,'first');

nexttile(t);
loglog(1+al, 1+al,'k--','LineWidth',1.2); hold on
for is = 1:numel(cal.levelSets)
    loglog(1+al, squeeze(cal.gain(iMid,:,is)),'o-','LineWidth',1.3);
end
grid on; xlabel('1 + alpha'); ylabel('measured gain');
legend([{'1 + alpha'} cellfun(@(c) ['levels ' c], cal.levelSets,'UniformOutput',false)], ...
    'Location','northwest','Box','off');
title(sprintf('Gain at %.3f px, share %.2f on levels %s', ...
    amps(iMid), cal.eta(4), cal.levelSets{4}));

nexttile(t);
for ia = 1:numel(amps)
    loglog(1+al, squeeze(cal.gain(ia,:,4)),'o-','LineWidth',1.2); hold on
end
loglog(1+al, 1+al,'k--'); grid on
xlabel('1 + alpha'); ylabel('measured gain');
legend([arrayfun(@(a) sprintf('%.3f px',a), amps,'UniformOutput',false) {'1 + alpha'}], ...
    'Location','northwest','Box','off');
title(sprintf('Saturation, levels %s', cal.levelSets{4}));

iB = find(cal.blur==s.blurSigma,1);  if isempty(iB), iB = 1; end
iA = find(al==20,1);                 if isempty(iA), iA = numel(al); end

nexttile(t);
for ig = 1:numel(cal.snrGain)
    loglog(max(al,1), squeeze(cal.detect(ig,iB,1,:)),'o-','LineWidth',1.3); hold on
end
yline(2,'--r'); grid on
xlabel('alpha (0 plotted at 1)'); ylabel('signal over its own null');
legend(arrayfun(@(g) sprintf('x%g averaged',g), cal.snrGain,'UniformOutput',false), ...
    'Location','best','Box','off');
r0 = squeeze(cal.detect(1,iB,1,1));  rEnd = squeeze(cal.detect(1,iB,1,end));
title(sprintf('Detection at %.3f px: %.1fx better by alpha %d, then flat', ...
    s.calFloorAmp, rEnd/r0, al(end)));

nexttile(t);
for is = 1:numel(cal.blurSets)
    loglog(cal.snrGain, squeeze(cal.floorPx(:,iB,is,iA)),'o-','LineWidth',1.5); hold on
end
yline(0.0753,'--','Color',[0.1 0.5 0.2]);
yline(0.0058,':' ,'Color',[0.1 0.5 0.2]);
grid on; xlabel('frames averaged before the magnifier'); ylabel('detection floor, px');
legend([cellfun(@(c) ['levels ' c], cal.blurSets,'UniformOutput',false), ...
        {'largest vessel 0.075 px','smallest 0.006 px'}],'Location','best','Box','off');
title(sprintf('Floor over %d frames at snr %.0f/px', s.phFrames, sqrt(s.phPhotons)));

nexttile(t);
yyaxis left
plot(cal.blur, cal.gainBlur(:,1),'o-','LineWidth',1.5); hold on
plot(cal.blur, cal.gainBlur(:,2),'s--','LineWidth',1.2);
ylabel('gain at alpha 20');
yyaxis right
plot(cal.blur, squeeze(cal.floorPx(1,:,1,iA)),'o-','LineWidth',1.5); hold on
plot(cal.blur, squeeze(cal.floorPx(1,:,2,iA)),'s--','LineWidth',1.2);
set(gca,'YScale','log'); ylabel('floor at one frame, px');
grid on; xlabel('blurSigma, pixels of the level');
title(sprintf('The blur trade: %s solid, %s dashed', cal.blurSets{1}, cal.blurSets{2}));

nexttile(t);
b = [cal.riesz(1).gain cal.riesz(2).gain; cal.riesz(1).floorPx*1000 cal.riesz(2).floorPx*1000];
bar(b,'EdgeColor','none'); grid on
set(gca,'XTickLabel',{'gain at alpha 20','floor, px x 1000'});
legend({cal.riesz.mode},'Location','best','Box','off');
title(sprintf('Three taps against the exact transform, levels %s', cal.blurSets{1}));

reportSave(rep, fh, 'calibration');
end
%------------- END OF CODE --------------
