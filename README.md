# Matlab Dynamic Light Scattering Library

A MATLAB toolbox for processing and analysing **laser speckle contrast imaging (LSCI)**
and related dynamic light scattering data — from raw recordings to blood-flow-index
maps, vascular segmentation, and pulsatility / vasomotion / neurovascular-coupling metrics.

> Headers, script formatting, and this README were prepared with the assistance of Claude Code.

---

## Overview

The library takes raw speckle recordings (`.rls`, `.mraw`, `.cxd`, `.dv`) and turns them into:

- **Contrast** and **blood flow index (BFI)** image sequences — `getK`, `getContrast`, `getContrastFromRLS`, `getBFI`;
- **Vascular segmentation** and per-segment traces — `getCategories`, `getSegmentation`, `setVesselTypes`;
- **Pulsatility**, **vasomotion**, **neurovascular-coupling** and **bolus / CTTH** metrics — `getPulsatility`, `getVasomotion`, `getInternalCycle`, `getExternalCycle`, `getCTTH`;
- optional **registration** across recordings (`registerLSCItoLSCI`, `registerRetinaLSCI`) and **Excel export** (`exportToExcel`).

A self-contained **myograph** toolset (`Myograph/`) measures vessel diameter from myograph videos and runs the same vasomotion analysis per user-defined interval.

---

## Requirements

- **MATLAB R2023b or newer** (some functions use recent language / graphics features).
- MathWorks toolboxes:
  - Image Processing Toolbox
  - Signal Processing Toolbox
  - Curve Fitting Toolbox
  - Statistics and Machine Learning Toolbox
  - Parallel Computing Toolbox — for the GPU paths (`procType='gpu'`) and `parfor`
  - Wavelet Toolbox — vasomotion analysis
- **Bio-Formats for MATLAB** (bundled in `3rd party/bfmatlab`) to read `.cxd` / `.dv` files.
- A CUDA-capable GPU is recommended for `procType='gpu'`; a CPU fallback is available.

---

## Repository structure

| Folder | Contents |
|---|---|
| `Basic functions/Laser Speckle Contrast Imaging` | Core contrast engine — `getK`, `getContrastFromRLS`, `getContrastFromMRAW`, `getEdgeSizeSLSCI` |
| `Basic functions/Read files` | Readers & file utilities — `readRLS`, `readMRAW`, `readCXD`, `readDV`, `cropRLS`, `fixMetaRLS`, `getPointerRLS`, `getFileNamesList`, `removeProcessedFiles` |
| `Basic functions/Dynamic Light Scattering Imaging` | g2 autocorrelation & DLSI/DCS model fitting — `getNormalizedG2`, `fitDLSI`, `fitNormalizedG2`, `getTauC` |
| `Basic functions/Post-processing` | FFT, speckle size, response features, overlap tables |
| `Basic functions` | Shared signal primitives — `getFFT`, `getVasomotionMetrics` (modular wavelet vasomotion core; `s.segVsmReturn` selects which levels — bands / moments / series / clustering / reconstruction / spectrum — are computed, default all six; `getVasomotion`'s per-pixel path is driven by the separate `s.ppxVsmReturn`; used by both `getVasomotion` and `getMyographVasomotion`), `assembleVasomotionTree` (builds the shared band-branched `results.vasomotion.<sig>` tree — scalars / fVectors / timeVectors / spectrum × VB/CB — from the core's flat metric bag), and `getPulsatilityMetrics` (modular harmonic pulsatility core mirroring `getVasomotionMetrics`: a two-mode SETUP/ANALYSIS that fits an `s.nHarm`-harmonic sine to one averaged cardiac cycle and returns model-free markers + harmonic coefficients; `s.segPulsReturn` / `s.ppxPulsReturn` select the levels — markers / model / reconstruction; used by `getPulsatility`) |
| `Basic functions/Registration` | Landmark / mask registration, retinal LSCI registration |
| `Basic functions/Visualisation and ROI selection` | Interactive ROI tools |
| `Processing routines` | High-level pipeline steps (contrast → categories → segmentation → BFI → cycles / pulsatility / vasomotion — `getVasomotion` writes the band-branched `results.vasomotion` tree per segment, and when `s.ppxVsmReturn` is non-empty also a LEAN per-pixel twin `results.vasomotion.ppx` — band-amplitude scalar `[Y×X]` maps plus an optional decimated `spectrum.amp`/`.phase` (`s.ppxVsmReturn` ∈ {`bands`,`spectrum`}); `getPulsatility` likewise writes the `results.pulsatility` tree per segment — `ps`/`pd`-prefixed markers + an `s.nHarm`-harmonic fit via the shared core `getPulsatilityMetrics` — and, when `s.ppxPulsReturn` is non-empty, the per-pixel twin `results.pulsatility.ppx`) |
| `Utility routines` | Cross-recording registration, region splitting, Excel export |
| `Launchers` | Ready-to-edit example pipelines — **start here** |
| `Post-processing examples` | Minimal, well-commented scripts for exploring processed results |
| `Myograph` | Myograph diameter + vasomotion suite (shares the wavelet core `getVasomotionMetrics` and the `assembleVasomotionTree` output tree; `getMyographVasomotion` returns one `<VSM>` tree stored as `intervals(iv).vasomotion`) |
| `Drafts` | Work-in-progress / superseded files — **not** part of the pipeline |
| `3rd party` | External libraries (Bio-Formats, superlets, …) — unmodified |

---

## Vasomotion output (`results.vasomotion`)

`getVasomotion` (LSCI / BFI) and `getMyographVasomotion` (diameter) share one wavelet
core (`getVasomotionMetrics`) and one output tree (`assembleVasomotionTree`). The LSCI
results file stores a **band-branched `results.vasomotion` tree**. The two frequency
bands — **VB** (vasomotion, `s.vFR`) and **CB** (comparison, `s.cFR`) — are struct
branches, so field names carry no VB/CB suffix.

**Root axes** (stored once, shared): `f` (frequency grid, Hz, descending), `timeWT`
(cone-of-influence-trimmed analysed time), `timeDWT` (decimated time base for
`spectrum.amp`/`.phase`), `pctCenters` (amplitude-percentile bin centres).

**Per-signal sub-trees** `results.vasomotion.<sig>`, one per analysed data-type signal
`<sig>` ∈ `{sData, dvsData, dvsDiameter, gsData}` (chosen by `s.vsmSignals`; `gsData`
carries its own axes on the `results.gsTime` base). Each holds four shape branches,
gated by `s.segVsmReturn`:

- **`scalars.VB` / `scalars.CB`** — temporal reductions of the band envelope:
  `ampMean` / `ampStd` / `ampSkew` / `ampPct`; VB-only frequency-shape & multiplicity
  `fCentMean/Std`, `fSprdMean/Std`, `shapePeak` (= `fCentMean/fSprdMean`),
  `nPeakMean/Std`; and clustering `durFlare*` / `durSilence*` / `ampFlare*` /
  `ampSilence*` (CB derived from its **own independent** Otsu flare/silence mask).
- **`fVectors`** (level `moments`) — per-frequency amplitude-spectrum moments `ampMean` / `ampStd` /
  `ampSkew`, plus `fVectors.VB.ampMeanPct` (percentile-resolved spectra) and
  `fVectors.{VB,CB}.ampFlare` / `.ampSilence`.
- **`timeVectors.VB`** — `amp` / `fCent` / `fSprd` / `nPeak` (`series`), `rData`
  (`reconstruction`), `maskFlare` (`clustering`); **`timeVectors.CB`** — `amp`,
  `maskFlare`.
- **`spectrum.amp`** / **`spectrum.phase`** (level `spectrum`) — the decimated wavelet
  **amplitude** and **phase** grids, each `[nSeg × nF × nD]`. Amplitude and phase are
  decimated SEPARATELY (recovered `|CWT|` and circular-mean phase), which preserves the
  fast-oscillation amplitude that coherent averaging of the complex coefficients would
  cancel — so `spectrum.amp` is a better magnitude estimate than `abs()` of a decimated
  complex grid.

**Per-pixel twin** `results.vasomotion.ppx` (when `s.ppxVsmReturn` is non-empty): a
**LEAN** per-pixel path — `s.ppxVsmReturn` is restricted to `{'bands','spectrum'}`. It
returns only band-amplitude `scalars` as `[Y×X]` maps (`ampMean` / `ampStd` always;
`ampSkew`, the `fCent*` / `fSprd*` / `shapePeak` centroid group and `ampPct` percentiles
are cost-gated) and, when `'spectrum'` is requested, the decimated `spectrum.amp` /
`spectrum.phase` as `[Y×X×nF×nD]` cubes. It has **no** `fVectors`, **no** `timeVectors`,
and never computes peaks (`nPeak`), series, clustering or reconstruction — those are
per-segment only. It reuses the root axes.

**Temporary per-pixel → segment averaging** `results.vasomotion.ppxs.coherent` /
`.incoherent` (when `s.ppxSegmentAveraging` is set; scaffolding to be removed): a
per-segment sub-tree built from an averaged **magnitude** spectrum, so it carries
`scalars` / `fVectors` / `timeVectors` but no `spectrum.amp`/`.phase` and no `timeVectors.VB.rData`.

**Summary columns** are duplicated from the tree into the per-segment metrics tables for
quick access: `sData → results.sMetrics` and `dvsData` + `dvsDiameter →
results.dvsMetrics` gain `ampMeanVB`, `ampMeanCB`, `fCentMean`, `shapePeak`,
`nPeakMean` (the `dvsDiameter` set is `_diam`-suffixed — `ampMeanVB_diam`, … — to share
`dvsMetrics` without collision); `gsData` writes **no** table. Row order matches the
metrics-table segment order.

Relative to the previous flat output, this redesign **drops** the kurtosis markers,
`bandRatio` / `flatness` / `peakRatio`, and the derivable
`occupancy` / `burstRate` / `medAmpFlare`, and **adds** the frequency-shape &
multiplicity metrics (`fCent*` / `fSprd*` / `shapePeak` / `nPeak*`), independent CB
flare/silence clustering, and the stored amplitude/phase wavelet grid (`spectrum.amp`/`.phase`). Myograph mirrors one `<sig>`
sub-tree as `intervals(iv).vasomotion` (single signal; no `ppx` / `ppxs`).

---

## Getting started

1. Clone the repository and keep **only one copy** on your machine — MATLAB caches paths and multiple versions cause conflicts.
2. Copy a launcher from `Launchers/` (or a script from `Post-processing examples/`) to your own working location; **leave the originals unchanged**.
3. In your copy, set `libraryFolder` to this repository, then run the `%% STEP` cells in order (run `STEP 0` once per MATLAB session).
4. Read the header (comments at the top) and inline comments of each function you use — do not run scripts blindly.

As a user you normally interact only with the **Launchers** and **Post-processing examples** folders; the processing steps you run are dictated by your protocol.

### File-naming convention

Each processing step appends flags to the original file name, e.g.:

- `_t_e_K_d.mat` — (t) temporal contrast, (e) estimated epoch, (K) contrast, (d) 3-D data;
- `_c_BFI_r.mat` — (c) cardiac cycle, (BFI) blood flow index, (r) results.

Most users only touch `_BFI_r.mat` (results) and `_BFI_d.mat` (3-D data) files. File selection is often done with wildcards / regular expressions — see `getFileNamesList`.

---

## Dependency graph

`A → B` means file *A* calls function *B*. The connected component is the main pipeline
(**Launchers → Processing routines → Basic functions**); isolated nodes are standalone
tools not hooked into that pipeline (the DLS-Imaging fitting functions, the file readers,
and the example scripts). `Drafts/` and `3rd party/` are excluded.

![File dependency graph](dependency_tree.png)

---

## Data acquisition protocols

Plan your file/folder naming in advance — it saves a lot of time later.

| Modality | Framerate | Duration | Key notes |
|---|---|---|---|
| **Calibration** | 20 Hz | 5 min | Paper sample (keep the same piece); exposure for ~30 % saturation (⟨I⟩ ≈ 90); large, evenly-lit FOV. |
| **Pulsatility** | 194 Hz | ≥ 2 min | Cranial window, awake optimal; FOV 1024×512; watch for dropped frames and physiological stability. |
| **Vasomotion** | ≥ 25 Hz (higher is better) | ≥ 10 min (30 optimal) | Cranial window, awake; a prior pulsatility recording helps map arteries vs veins. |
| **Neurovascular coupling (NVC)** | ≥ 25 Hz | 120 s baseline + 20 × (5 s @3 Hz stim + 25 s recovery) | Register stim timing to LSCI; longer recovery / more repeats are better; monitor animal state. |
| **Slow dynamics (e.g. vasoreactivity)** | ≥ 5 Hz | ≥ 30 min | Register stimulus timing; ensure position stability (else split into multiple files). |
| **Multiple conditions (e.g. stroke)** | — | — | Keep the system configuration constant; match FOV; prioritise tilt/focus matching over translation. |
| **Longitudinal / cross-group** | — | — | Calibrate at each timepoint (ideally each day); keep configuration constant; match FOV. |

---

## Launchers

Copy a launcher to your own location and edit it for your project.

| Launcher | Purpose |
|---|---|
| `Launcher_basic` | Minimal contrast → BFI, for a first look at your data. |
| `Launcher_pulsatility` | Cardiac pulsatility (cranial window, ≥ 194 fps advised). |
| `Launcher_vasoreactivity` | Vasoreactivity / vasomotion (≥ 25 fps for vasomotion). |
| `Launcher_pulsatility_vasomotion`, `Launcher_pulsatility_vasomotion_Line` | Combined pulsatility + vasomotion (line-scan variant). |
| `Launcher_NVC` | Externally-triggered neurovascular-coupling responses. |
| `Launcher_CTTH` | Bolus-injection (CTTH) analysis from wide-field fluorescence. |
| `Launcher_Sonam_LSCI` | Project-specific combined pulsatility + vasomotion. |
| `Myograph/Launcher_Myograph` | Myograph vessel-diameter + vasomotion (self-contained). |

## Post-processing examples

Small, heavily-commented scripts that consume processed `_BFI_d.mat` / `_BFI_r.mat` files.
Each script's header explains what it does, when to use it, and embeds a ready-to-reuse
generative-AI prompt so you can regenerate or adapt it for your own data.

| Example | Purpose |
|---|---|
| `Example_Claude_plotROIs` | Draw ROIs on a recording, plot their signals and report per-ROI values. |
| `Example_Claude_pickRegions` | Click labelled segments in a results file and print their metrics (no data file needed). |
| `Example_Claude_groupSummary` | Compare BFI / diameter / pulsatility / vasomotion across groups of recordings. |
| `Example_Claude_arteryVeinStats` | Compare arteries vs veins in one recording (needs `setVesselTypes`). |
| `Example_Claude_bfiMovie` | Play back (and optionally export) a recording as a movie — a good first look. |

Getting a "feel" for your data with basic MATLAB skills is recommended. Follow any pipeline
to a `_BFI_d.mat` file (3-D data `X × Y × Time` plus a time vector) and go from there —
ROI selection, filtering, plotting as image or video. Each `Example_Claude_*` script keeps a
generative-AI prompt in its header that you can reuse and adapt.

---

## License

This work is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0
International (CC BY-NC-SA 4.0)** license — see [`LICENSE.md`](LICENSE.md). You may share and
adapt the material for **non-commercial** purposes, with attribution and under the same
license; commercial use requires separate permission from the author.

© 2026 Dmitry D. Postnov, CFIN, Aarhus University.

## Contact & citation

**Dmitry D. Postnov** — CFIN, Aarhus University — <dpostnov@cfin.au.dk>

If you use this toolbox in your research, please cite the relevant LSCI / DLSI publications
from the Postnov lab.
