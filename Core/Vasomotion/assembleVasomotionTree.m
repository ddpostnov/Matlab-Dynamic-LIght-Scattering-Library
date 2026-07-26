%assembleVasomotionTree  Assemble one band-branched <VSM> tree from a flat metric bag
%
%   T = assembleVasomotionTree(acc,want,shp) builds a single vasomotion sub-tree
%   (the <VSM> unit: scalars / fVectors / timeVectors / complex, each split into the
%   VB and CB struct branches) from the flat accumulator bag `acc` produced by
%   getVasomotionMetrics.  It is the SINGLE source of truth for the tree field
%   layout, shared by every producer so their outputs are structurally identical:
%   runVasomotion's per-segment, per-pixel (ppx) and per-segment-averaging (ppxs)
%   paths, and getMyographVasomotion.  Only the leaf SHAPES differ between producers
%   (per-segment [nSeg x ...] vs per-pixel [Y x X x ...] maps/cubes); those are
%   supplied by the caller's `shp` reshape handles, so this function fixes the field
%   NAMES and branch structure once and for all.
%
%   Bands are held as struct BRANCHES (so field names carry no VB/CB suffix):
%   VB = vasomotion band (s.vFR), CB = comparison/control band (s.cFR).  Only VB
%   carries the frequency-shape / multiplicity fields (fCent*/fSprd*/shapePeak/
%   nPeak*).  Every branch is gated by the selected analysis level in `want`, so an
%   unrequested level is simply absent from T.
%
%   INPUTS
%     acc   flat struct of accumulated per-item quantities (fields named below).
%           Every field a selected level reads must be present; unselected levels'
%           fields may be omitted.  Leading dimension = items (segments / pixels / Y).
%             scalars   : vbAmpMean vbAmpStd vbAmpSkew vbAmpPct
%                         cbAmpMean cbAmpStd cbAmpSkew cbAmpPct
%                         fCentMean fCentStd fSprdMean fSprdStd shapePeak nPeakMean nPeakStd
%                         vbDurFlareMean/Std vbDurSilenceMean/Std vbAmpFlareMean/Std vbAmpSilenceMean/Std
%                         cbDurFlareMean/Std cbDurSilenceMean/Std cbAmpFlareMean/Std cbAmpSilenceMean/Std
%             fVectors  : spctMean spctStd spctSkew vbSpctPct
%                         vbSpctFlare vbSpctSilence cbSpctFlare cbSpctSilence
%             timeVectors: tsVB tsCB fCentTV fSprdTV nPeakTV rData flareVB flareCB
%             spectrum  : amp phase   (real amplitude & phase grids; was the complex WT)
%     want  struct of logical level flags from getVasomotionMetrics layout.want, in the
%           Session-2 FINE shape: the five band tokens .bandsAmp .bandsSkew .bandsShape
%           .bandsPct .bandsPeak (each gates its own scalar sub-group, so a lean caller
%           can request just .bandsAmp), plus .moments .series .clustering .reconstruction
%           .spectrum.  (The band scalars split lets the lean per-pixel path skip the
%           expensive nPeak/centroid work.)  For backward compatibility this function also
%           accepts the pre-redesign COARSE shape (.bands / .spectrum=moments /
%           .complex=amp-phase) still emitted by the frozen TEMPORARY per-segment averaging
%           block in runVasomotion; that shim is deleted together with the temp block.
%     shp   struct of reshape handles applied to each leaf by shape class:
%             .sc    scalar        [nItem x 1]      -> item layout (e.g. [Y x X])
%             .scPct percentile    [nItem x nPct]
%             .fv    per-frequency [nItem x nF]
%             .pc    pct spectrum  [nItem x nF x nBin]
%             .tv    time vector   (per-segment [nT x nItem] identity, or per-pixel
%                                   [nItem x nT] -> [Y x X x nT])
%             .spec  amp/phase grid [nItem x nF x nDec]
%           Use identity handles (@(A)A) when the accumulators are already in tree shape.
%
%   OUTPUT
%     T     the <VSM> sub-tree: T.scalars.{VB,CB}, T.fVectors[.{VB,CB}],
%           T.timeVectors.{VB,CB}, T.spectrum.amp/.phase (each present per `want`).  Shared
%           axes (f/timeWT/timeDWT/pctCenters) are the caller's responsibility - this
%           builds the four branches only.
%
%   DEPENDS ON
%     none (pure struct assembly).  Field set defined by
%     Docs/vasomotion-metrics-redesign/DATA-MODEL.md (git-excluded planning artifact).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

function T=assembleVasomotionTree(acc,want,shp)
T=struct();

% Normalise `want` onto the Session-2 fine tokens.  New callers (runVasomotion's
% per-segment and lean per-pixel paths, getMyographVasomotion) pass the fine shape
% directly.  The frozen TEMPORARY per-segment averaging block still passes the
% pre-redesign COARSE shape (.bands / .spectrum=moments / .complex=amp-phase); map it so
% that untouched caller keeps producing the same ppxs tree.  This legacy branch dies with
% the temp block - once it is gone, `want` is always the fine shape and this collapses to
% the else arm.
if isfield(want,'bands')                                   % legacy coarse shape (temp caller)
    [wBAmp,wBSkew,wBShape,wBPct,wBPeak]=deal(want.bands);
    wMoments =want.spectrum;                                % legacy 'spectrum' level == per-frequency moments
    wSpectrum=want.complex;                                 % legacy 'complex'  level == amp/phase grid
else                                                       % Session-2 fine shape
    wBAmp=want.bandsAmp; wBSkew=want.bandsSkew; wBShape=want.bandsShape;
    wBPct=want.bandsPct; wBPeak=want.bandsPeak;
    wMoments=want.moments; wSpectrum=want.spectrum;
end

% ---- scalars: band envelope reductions (finely gated) + VB frequency-shape + clustering ----
if wBAmp
    T.scalars.VB.ampMean=shp.sc(acc.vbAmpMean);   T.scalars.VB.ampStd=shp.sc(acc.vbAmpStd);
    T.scalars.CB.ampMean=shp.sc(acc.cbAmpMean);   T.scalars.CB.ampStd=shp.sc(acc.cbAmpStd);
end
if wBSkew
    T.scalars.VB.ampSkew=shp.sc(acc.vbAmpSkew);   T.scalars.CB.ampSkew=shp.sc(acc.cbAmpSkew);
end
if wBPct
    T.scalars.VB.ampPct=shp.scPct(acc.vbAmpPct);  T.scalars.CB.ampPct=shp.scPct(acc.cbAmpPct);
end
if wBShape
    T.scalars.VB.fCentMean=shp.sc(acc.fCentMean); T.scalars.VB.fCentStd=shp.sc(acc.fCentStd);
    T.scalars.VB.fSprdMean=shp.sc(acc.fSprdMean); T.scalars.VB.fSprdStd=shp.sc(acc.fSprdStd);
    T.scalars.VB.shapePeak=shp.sc(acc.shapePeak);
end
if wBPeak
    T.scalars.VB.nPeakMean=shp.sc(acc.nPeakMean); T.scalars.VB.nPeakStd=shp.sc(acc.nPeakStd);
end
if want.clustering
    T.scalars.VB.durFlareMean=shp.sc(acc.vbDurFlareMean);     T.scalars.VB.durFlareStd=shp.sc(acc.vbDurFlareStd);
    T.scalars.VB.durSilenceMean=shp.sc(acc.vbDurSilenceMean); T.scalars.VB.durSilenceStd=shp.sc(acc.vbDurSilenceStd);
    T.scalars.VB.ampFlareMean=shp.sc(acc.vbAmpFlareMean);     T.scalars.VB.ampFlareStd=shp.sc(acc.vbAmpFlareStd);
    T.scalars.VB.ampSilenceMean=shp.sc(acc.vbAmpSilenceMean); T.scalars.VB.ampSilenceStd=shp.sc(acc.vbAmpSilenceStd);
    T.scalars.CB.durFlareMean=shp.sc(acc.cbDurFlareMean);     T.scalars.CB.durFlareStd=shp.sc(acc.cbDurFlareStd);
    T.scalars.CB.durSilenceMean=shp.sc(acc.cbDurSilenceMean); T.scalars.CB.durSilenceStd=shp.sc(acc.cbDurSilenceStd);
    T.scalars.CB.ampFlareMean=shp.sc(acc.cbAmpFlareMean);     T.scalars.CB.ampFlareStd=shp.sc(acc.cbAmpFlareStd);
    T.scalars.CB.ampSilenceMean=shp.sc(acc.cbAmpSilenceMean); T.scalars.CB.ampSilenceStd=shp.sc(acc.cbAmpSilenceStd);
end

% ---- fVectors (moments: band-agnostic per-frequency moments + VB.ampMeanPct; clustering: flare/silence) ----
if wMoments
    T.fVectors.ampMean=shp.fv(acc.spctMean); T.fVectors.ampStd=shp.fv(acc.spctStd); T.fVectors.ampSkew=shp.fv(acc.spctSkew);
    T.fVectors.VB.ampMeanPct=shp.pc(acc.vbSpctPct);
end
if want.clustering
    T.fVectors.VB.ampFlare=shp.fv(acc.vbSpctFlare); T.fVectors.VB.ampSilence=shp.fv(acc.vbSpctSilence);
    T.fVectors.CB.ampFlare=shp.fv(acc.cbSpctFlare); T.fVectors.CB.ampSilence=shp.fv(acc.cbSpctSilence);
end

% ---- timeVectors (series: amp/fCent/fSprd/nPeak; reconstruction: rData; clustering: maskFlare) ----
if want.series
    T.timeVectors.VB.amp=shp.tv(acc.tsVB); T.timeVectors.VB.fCent=shp.tv(acc.fCentTV);
    T.timeVectors.VB.fSprd=shp.tv(acc.fSprdTV); T.timeVectors.VB.nPeak=shp.tv(acc.nPeakTV);
    T.timeVectors.CB.amp=shp.tv(acc.tsCB);
end
if want.reconstruction
    T.timeVectors.VB.rData=shp.tv(acc.rData);
end
if want.clustering
    T.timeVectors.VB.maskFlare=shp.tv(acc.flareVB);
    T.timeVectors.CB.maskFlare=shp.tv(acc.flareCB);
end

% ---- spectrum: amplitude & phase of the time-frequency grid (was the complex WT) ----
if wSpectrum
    T.spectrum.amp=shp.spec(acc.amp);
    T.spectrum.phase=shp.spec(acc.phase);
end
end
