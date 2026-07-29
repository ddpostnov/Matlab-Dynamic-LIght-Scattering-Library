# Matlab Dynamic Light Scattering Library

A MATLAB toolbox for processing and analysing **laser speckle contrast imaging (LSCI)**
and related dynamic light scattering data — from raw recordings to blood-flow-index
maps, vascular segmentation, and pulsatility / vasomotion / neurovascular-coupling metrics.

> Headers, script formatting, and this README were prepared with the assistance of Claude Code.

---

## Overview

The library takes raw speckle recordings (`.rls`, `.mraw`, `.cxd`, `.dv`) and turns them into:

- **Contrast** and **blood flow index (BFI)** image sequences — `getK`, `runContrastFromRLS`, `getContrastFromRLS`, `runBFI`;
- **Vascular segmentation** and per-segment traces — `setRegions`, `runSegmentation`, `runDynamicSegmentation`, `setVesselTypes`;
- **Pulsatility**, **vasomotion**, **neurovascular-coupling** and **bolus / CTTH** metrics — `runPulsatility`, `runVasomotion`, `runInternalCycle`, `runExternalCycle`, `runCTTH`;
- optional **registration** across recordings (`runRegistration`, `registerRetinaLSCI`) and **Excel export** (`exportToExcel`).

You drive that pipeline in one of two ways: from the **Processing Workbench** (`GUIs/guiWorkbench`) — an interactive file × step matrix that ticks, parametrises and runs the steps for you, and is the recommended starting point — or from a **launcher script** (`Launchers/`), the scripted alternative for reproducible, batch, or headless work. Both call exactly the same `run…` / `set…` steps.

A **myograph** toolset — the headless `Core/Myograph/` functions plus the `GUIs/guiMyograph` app — measures vessel diameter from myograph videos and runs the same vasomotion analysis per user-defined interval.

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

The top level is organised by **layer**: `Core` (the math), `Wrappers` (the pipeline
steps / seam), `Launchers` (orchestration templates), and the consumers (`GUIs`,
`Utilities`, `Simulation`).

| Folder | Contents |
|---|---|
| `Core/Read files` | Readers & file utilities — `readRLS`, `readMRAW`, `readCXD`, `readDV`, `cropRLS`, `fixMetaRLS`, `getPointerRLS`, `getFileNamesList`, `removeProcessedFiles` |
| `Core/Laser Speckle Contrast Imaging` | Contrast engine — `getK`, `getContrastFromRLS`, `getContrastFromMRAW`, `getSpeckleSize`, `getEdgeSizeSLSCI`; segmentation cores — `enhanceForDisplay`, `getPixelCategories` (5-level category mask), `getSegmentationLabels` (indexed vessel/parenchyma maps), `showSegmentsPreview` |
| `Core/Dynamic Light Scattering Imaging` | g2 autocorrelation & DLSI/DCS model fitting — `getNormalizedG2`, `fitDLSI`, `getTauC` |
| `Core/Registration` | Landmark / mask registration — `registerToReference`, `enhanceForRegistration`, `registerRetinaLSCI`, `manualByPointRegistration` |
| `Core/Vasomotion` | `getVasomotionMetrics` (modular wavelet vasomotion core; `s.segVsmReturn` selects which levels — bands / moments / series / clustering / reconstruction / spectrum — are computed, default all six; `runVasomotion`'s per-pixel path is driven by the separate `s.ppxVsmReturn`; used by both `runVasomotion` and `getMyographVasomotion`) and `assembleVasomotionTree` (builds the shared band-branched `results.vasomotion.<sig>` tree — scalars / fVectors / timeVectors / spectrum × VB/CB — from the core's flat metric bag) |
| `Core/Pulsatility` | `getPulsatilityMetrics` (modular harmonic pulsatility core mirroring `getVasomotionMetrics`: a two-mode SETUP/ANALYSIS that fits an `s.nHarm`-harmonic sine to one averaged cardiac cycle and returns model-free markers + harmonic coefficients; `s.segPulsReturn` / `s.ppxPulsReturn` select the levels — markers / model / reconstruction; used by `runPulsatility`) |
| `Core/Vasculature` | Vascular hierarchy derivation — `getVascularTree` (pure staged flow-potential parent→daughter tree with FOV-bridging, Horton-Strahler order & generation; consumed by `setVascularTree`), plus the helpers `orderForest`, `getMetric` and `defaultFlowParams` |
| `Core/Shared` | Shared signal primitives — `getFFT` |
| `Core/Myograph` | Myograph diameter / vasomotion / propagation suite, headless (shares the wavelet core `getVasomotionMetrics` and the `assembleVasomotionTree` output tree; `getMyographVasomotion` returns one `<VSM>` tree stored as `intervals(iv).vasomotion`) |
| `Wrappers` | High-level pipeline steps — the `run…` / `set…` functions that read and write the `_d`/`_r`/`_s` file triplet (contrast → regions → segmentation → BFI → cycles / pulsatility / vasomotion; segmentation is `setRegions` (interactive multi-ROI editor → `results.regionsMask`) → `runSegmentation` (fully automatic categorize + label + per-segment traces; `s.fNamesCopyTo` copies the segmentation onto co-registered siblings, replacing the old `assignCategories`) → optional `runDynamicSegmentation` (per-frame vessel diameter / flow) — `runVasomotion` writes the band-branched `results.vasomotion` tree per segment, and when `s.ppxVsmReturn` is non-empty also a LEAN per-pixel twin `results.vasomotion.ppx` — band-amplitude scalar `[Y×X]` maps plus an optional decimated `spectrum.amp`/`.phase` (`s.ppxVsmReturn` ∈ {`bands`,`spectrum`}); `runPulsatility` likewise writes the `results.pulsatility` tree per segment — `ps`/`pd`-prefixed markers + an `s.nHarm`-harmonic fit via the shared core `getPulsatilityMetrics` — and, when `s.ppxPulsReturn` is non-empty, the per-pixel twin `results.pulsatility.ppx`), plus `runRegistration`, `splitRegions`, and the guided front-ends (`runGuidedContrast`, `runGuidedIntensity`) |
| `Launchers` | Ready-to-edit example pipelines — the scripted way to drive the same steps |
| `GUIs` | Interactive apps — `guiWorkbench` (**the Processing Workbench: start here**; the file × step matrix that runs the LSCI pipeline, with Export and Explore tabs), `guiExport` (the **standalone** export tool: pick files, scan a folder or load a workbench session, choose parameters and how to average them, write one workbook per recording or one merged workbook for statistics — the interactive route to the same `exportToExcel` the launchers call) and `guiMyograph` (myograph workbench for `.avi` diameter / vasomotion / propagation; double-click `launchMyographWorkbench.vbs` to open without starting MATLAB by hand) |
| `GUIs/workbench` | The workbench's own components — the headless brain (`wbStepRegistry` the step specs, `wbDiscoverFiles`, `wbFileModel`, `wbStateEngine`, `wbSettingsModel`, `wbInvalidate`, `wbExecutor`, `wbArtifacts`, `wbModalGuard`, `wbSession`) plus `guiExplore`, the results explorer hosted in the Explore tab (still openable standalone, or embeddable elsewhere with `guiExplore('Parent',container)`) |
| `Utilities` | Terminal consumers of finished results — `exportToExcel` (`exportToExcel(fNames)` writes the full workbook; the optional `exportToExcel(fNames,opts)` selects `opts.sheets` / `opts.format`, averages over labels or vessel type with `opts.groupBy` / `opts.weightByArea`, subsets `opts.columns`, and merges every file into one labelled workbook with `opts.merge` / `opts.labels` / `opts.outFile`. `GUIs/guiExport` is the interactive front-end for exactly these options) |
| `Simulation` | Synthetic dynamic-speckle generation (`getDynamicSpeckles`, `Launcher_speckleSimulation`) — self-contained |
| `3rd party` | External libraries (Bio-Formats, superlets, …) — unmodified |

---

## Vasomotion output (`results.vasomotion`)

`runVasomotion` (LSCI / BFI) and `getMyographVasomotion` (diameter) share one wavelet
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
2. Put the library on the path: `addpath(genpath('<path to this repository>'))`.
3. Run `guiWorkbench` — the [Processing Workbench](#processing-workbench). Load your recordings, tick the steps your protocol needs, press **Run**, then use its **Export** and **Explore** tabs.
4. Prefer a script? Copy a launcher from `Launchers/` to your own working location (**leave the originals unchanged**), set `libraryFolder` in your copy, and run the `%% STEP` cells in order (`STEP 0` once per MATLAB session). See [Launchers](#launchers).
5. Either way, read the header (comments at the top) and inline comments of each function you use — do not run steps blindly.

The processing steps you run are dictated by your protocol, not by the entry point: the workbench and the launchers call the same `run…` / `set…` wrappers and write the same `_d` / `_r` / `_s` files, so you can move between them freely. The `.avi` myograph workflow is not part of the workbench — use `guiMyograph` for it.

### File-naming convention

Each processing step appends flags to the original file name, e.g.:

- `_t_e_K_d.mat` — (t) temporal contrast, (e) estimated epoch, (K) contrast, (d) 3-D data;
- `_c_BFI_r.mat` — (c) cardiac cycle, (BFI) blood flow index, (r) results.

Most users only touch `_BFI_r.mat` (results) and `_BFI_d.mat` (3-D data) files. File selection is often done with wildcards / regular expressions — see `getFileNamesList`.

---

## Dependency graph

`A → B` means file *A* calls function *B*. The connected component is the main pipeline
(**Launchers → Wrappers → Core**); isolated nodes are standalone
tools not hooked into that pipeline (the DLS-Imaging fitting functions and the file
readers). `Drafts/` and `3rd party/` are excluded.

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

## Processing Workbench

`guiWorkbench` is the primary entry point: one window that turns the launcher workflow
into a spreadsheet. **Rows are recordings** (blocked by animal, reference first),
**columns are pipeline steps**, and each **cell** is the state of one step for one
recording — *ready* (tick it to queue), *checked*, *running*, *done*, *stale* (a setting
or an upstream re-run invalidated it) or *unavailable* (its prerequisites are not met).
State is read from the files on disk, so a session you closed last week comes back with
everything already done still marked done.

Every file carries three **independent labels**, each a regexp match over its name that
you can override by hand:

| label | what it is | who uses it |
|-------|-----------|-------------|
| **animal** | the subject — same brain, same FOV, co-registerable. Owns **one reference recording**, valid for all of its files. | registration, vessel typing, the reference |
| **type** | the recording's experimental role (`BV`/`BN`/`BP`, `ctrl`/`stim`/`washout`, `1`/`2`/`3` — whatever your names carry) | the processing configuration |
| **group** | the **experimental** group, a comparison label (`KO`/`WT`, `pre`/`post`) | Export and Explore only; processing ignores it |

They do not nest: a group may span several animals and one animal may span several
groups, so none is ever derived from another. **There is no built-in list of animals,
types or groups** — the values are simply whatever your regexps match or you type, and a
name that matches nothing lands in an ordinary `(untyped)` / `(ungrouped)` bucket rather
than being an error.

```matlab
addpath(genpath('<path to this repository>'))
guiWorkbench
```

1. **1 - Files** — build and curate the working set. Point **root** at a folder, put an
   extension or glob in **files** (`*.rls`, `*_t_K_d.mat`) and press **Scan** to recurse
   the whole tree; **Add files…** and **Add folder…** add more by hand without throwing
   the scan away. The five regexp boxes (**Animal**, **Type**, **Rec. index**,
   **Group**, **Reference**) label what the scan finds — hover each for worked examples.
   The table below is editable: click a column header to sort, and edit `animal`,
   `type`, `index`, `group` or `modality` in place. To label many files at once — the
   fast way — select the rows and use **Quick assign**: pick the field, pick or type the
   value, press **Apply to selected rows**. A hand-made label beats the regexp and
   survives a re-scan. `modality` is a dropdown constrained by the file extension (a
   `.rls` can be `LSCI`/`HSLSCI`/`DLSI`, a video can be `WMYO`/`EPFL`/`LSCI`, …), so an
   impossible pairing is simply refused. **Delete selected** drops rows from the working
   set — **nothing is ever deleted from disk**.
   Tick `ref` to pin an animal's **reference recording**: one per animal, valid for all
   of its files whatever their type or group, replacing whatever the Reference regexp
   picked. It is stored as a *recording identity*, so each step still resolves the
   branch (`_t` / `_c` / …) it needs.
   The banner under the table is the **gate**: every file needs an animal, a type, an
   index and a group, and file names must be unique across the scanned tree (the
   workbench and the pipeline both identify recordings by name). Until that holds, the
   other tabs stay out of reach. An animal with no reference only warns — registration
   and vessel typing are the steps that need one.
   **Save session… / Load session…** store the curated list, its labels, the references,
   the settings and the cell states in a versioned `.mat` sidecar, so a reload
   reproduces the table with no re-scan and long analyses survive a MATLAB restart.
2. **2 - Process** — tick the cells you want (**Check all** / **Clear checks** queue or
   clear every runnable cell at once). Pick a step in the **Step** dropdown on the right
   to open its **settings panel**; values shared between steps propagate automatically,
   and editing a setting marks that step and everything downstream *stale* so you cannot
   silently mix parameters. The preset dropdown seeds the whole parameter set, and
   **Save preset… / Load preset…** store your own — the modern replacement for keeping an
   edited copy of a launcher. **Preview order** lists
   exactly what would run, in dependency order, without calling anything; **Run**
   executes it, streaming per-cell progress and the command-window messages into the log
   pane. **Stop** cancels cooperatively between files. A failed step marks its cell as an
   error and the batch continues. Finished cells with a report image become clickable and
   open it in the reports panel. Interactive steps (region drawing, vessel typing) open
   their own editor and grey the workbench out until you finish.
3. **3 - Export** — pick which of the loaded recordings to export, tick the sheets you
   want (**All** / **None**), choose the format, and press **Export selected**. This is a
   selection UI over `exportToExcel`, which remains callable directly as
   `exportToExcel(fNames)` or `exportToExcel(fNames,opts)`. For the full export —
   choosing parameters, averaging over labels or vessel type, and merging every file into
   one workbook whose rows carry animal / type / experimental group — run the standalone
   **`guiExport`** (`GUIs/guiExport.m`), which needs no workbench: point it at files, at a
   folder, or at a saved workbench session.
4. **4 - Explore** — the results explorer (`GUIs/workbench/guiExplore`) hosted in-tab.
   Press **Load workbench files & animals** to seed it with the `_r.mat` results of the
   files you loaded, one explorer group per animal, then plot single recordings or group
   comparisons and export publication figures. See
   [Exploring processed results](#exploring-processed-results).

The workbench covers the **LSCI** pipeline. It never reimplements any processing — it
orders, gates, parametrises and calls the same wrappers a launcher would.

---

## Launchers

The launchers are the **scripted alternative** to the workbench — fully supported, and the
right choice when you want a reproducible record of a pipeline, a headless/batch run, or a
starting point to modify. Copy a launcher to your own location and edit it for your project.

| Launcher | Purpose |
|---|---|
| `Launcher_basic` | Minimal contrast → BFI, for a first look at your data. |
| `Launcher_pulsatility` | Cardiac pulsatility (cranial window, ≥ 194 fps advised). |
| `Launcher_vasoreactivity` | Vasoreactivity / vasomotion (≥ 25 fps for vasomotion). |
| `Launcher_pulsatility_vasomotion` | Combined pulsatility + vasomotion. |
| `Launcher_NVC` | Externally-triggered neurovascular-coupling responses. |
| `Launcher_CTTH` | Bolus-injection (CTTH) analysis from wide-field fluorescence. |
| `Launcher_DLSI_basic` | DLSI pipeline: raw `.mraw` recordings → per-pixel g2 / decorrelation-time fit. |
| `Launcher_guided` | Guided full-resolution per-segment trace extraction (demo on the bundled test data). |

The **myograph** workflow has neither a launcher script nor a workbench column — open
`GUIs/guiMyograph` in MATLAB, or double-click `GUIs/launchMyographWorkbench.vbs`.

---

## Exploring processed results

Follow any pipeline to a `_BFI_d.mat` file (3-D data `X × Y × Time` plus a time vector)
and explore it with basic MATLAB skills — ROI selection, filtering, plotting as image or
video.

For the finished `_r.mat` results (per-segment BFI / diameter / pulsatility / vasomotion,
single files or group comparisons), use the **Explore** tab of `guiWorkbench`: it is
seeded from the recordings you already loaded (one explorer group per animal, which you
can relabel there), matches the plot type to the
data (time series, box plots, spectra, spectrograms, maps), and exports the figure at a
chosen DPI. The same explorer can be opened on its own with `guiExplore`
(`GUIs/workbench/guiExplore.m`) if you want to browse results without loading a workbench
session, and embedded elsewhere with `guiExplore('Parent',container)`.

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
