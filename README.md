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

You drive that pipeline in one of two ways: from the **Processing Workbench** (`GUIs/guiWorkbench`) — curate your recordings, configure the pipeline once per recording **type**, press Run, and it orders, parametrises and calls the steps for you; the recommended starting point — or from a **launcher script** (`Launchers/`), the scripted alternative for reproducible, batch, or headless work. Both call exactly the same `run…` / `set…` steps.

A **myograph** toolset runs the same vasomotion analysis per user-defined interval on both kinds of myograph recording, **in the Processing Workbench** like any other recording, each writing one `_MYO` triplet. A **pressure** myograph `.avi` goes optionally Pre-set intervals *or* Time crop, then Video → Diameter → Intervals → Propagation → Vasomotion, measuring the vessel's diameter from the video. A **wire** myograph LabChart `.adicht` goes LabChart → Intervals → Vasomotion: reading the file is the measurement, and each interval is analysed on the channels chosen for it. The headless `Core/Myograph/` functions are behind both, and the interactive steps open one window per recording (`editMyographIntervals` — the same window with the video column replaced by the channel list and the recording's comments), exactly as `setRegions` does. There is **one myograph implementation and one file format**: the standalone `guiMyograph` app and the `*_myograph.mat` format it wrote were retired to `Drafts/` once the workbench route was complete.

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
- **ADInstruments SDK for MATLAB** (bundled in `3rd party/adinstruments_sdk_matlab`) to read LabChart `.adicht` files — **Windows and 64-bit MATLAB only**, since it is a MEX wrapper over ADInstruments' own DLL. Only the wire-myograph pipeline needs it.
- A CUDA-capable GPU is recommended for `procType='gpu'`; a CPU fallback is available.

---

## Repository structure

The top level is organised by **layer**: `Core` (the math), `Wrappers` (the pipeline
steps / seam), `Launchers` (orchestration templates), and the consumers (`GUIs`,
`Utilities`, `Simulation`).

| Folder | Contents |
|---|---|
| `Core/Read files` | Readers & file utilities — `readRLS`, `readMRAW`, `readCXD`, `readDV`, `readLabChart` (LabChart `.adicht`: channels on their own sampling rates, comments and record boundaries — the only file that touches the `+adi` package), `cropRLS`, `fixMetaRLS`, `getPointerRLS`, `getFileNamesList`, `removeProcessedFiles`, `getProductPath` (the `_d`/`_r`/`_s` member beside a product), `getResultsPath` (where a file goes when the results are kept apart from the recordings — the one place that rule is written) |
| `Core/Laser Speckle Contrast Imaging` | Contrast engine — `getK`, `getContrastFromRLS`, `getContrastFromMRAW`, `getSpeckleSize`, `getEdgeSizeSLSCI`; segmentation cores — `enhanceForDisplay`, `getPixelCategories` (5-level category mask), `getSegmentationLabels` (indexed vessel/parenchyma maps), `showSegmentsPreview` |
| `Core/Dynamic Light Scattering Imaging` | g2 autocorrelation & DLSI/DCS model fitting — `getNormalizedG2`, `fitDLSI`, `getTauC` |
| `Core/Registration` | Landmark / mask registration — `registerToReference`, `enhanceForRegistration`, `registerRetinaLSCI`, `manualByPointRegistration` |
| `Core/Vasomotion` | `getVasomotionMetrics` (modular wavelet vasomotion core; `s.segVsmReturn` selects which levels — bands / moments / series / clustering / reconstruction / spectrum — are computed, default all six; `runVasomotion`'s per-pixel path is driven by the separate `s.ppxVsmReturn`; used by both `runVasomotion` and `getMyographVasomotion`) and `assembleVasomotionTree` (builds the shared band-branched `results.vasomotion.<sig>` tree — scalars / fVectors / timeVectors / spectrum × VB/CB — from the core's flat metric bag) |
| `Core/Pulsatility` | `getPulsatilityMetrics` (modular harmonic pulsatility core mirroring `getVasomotionMetrics`: a two-mode SETUP/ANALYSIS that fits an `s.nHarm`-harmonic sine to one averaged cardiac cycle and returns model-free markers + harmonic coefficients; `s.segPulsReturn` / `s.ppxPulsReturn` select the levels — markers / model / reconstruction; used by `runPulsatility`) |
| `Core/Vasculature` | Vascular hierarchy derivation — `getVascularTree` (pure staged flow-potential parent→daughter tree with FOV-bridging, Horton-Strahler order & generation; consumed by `setVascularTree`), plus the helpers `orderForest`, `getMetric` and `defaultFlowParams` |
| `Core/Shared` | Shared signal primitives — `getFFT` |
| `Core/Myograph` | Myograph diameter / vasomotion / propagation suite, headless (shares the wavelet core `getVasomotionMetrics` and the `assembleVasomotionTree` output tree; `getMyographVasomotion` returns one `<VSM>` tree stored as `intervals(iv).vasomotion`). `myographProduct` is the one place the `_MYO` triplet is opened, saved and named — exactly one function creates it (the entry step), every other step is load-modify-save; `editMyographIntervals` is the shared interval editor (one blocking window per recording, in a pressure or a wire mode) and `getMyographTrace` builds what it draws — the measured diameter, a coarse brightness profile read from the video when nothing has been measured yet, or the recorded channels with their comments for a wire myograph |
| `Wrappers` | High-level pipeline steps — the `run…` / `set…` functions that read and write the `_d`/`_r`/`_s` file triplet (contrast → regions → segmentation → BFI → cycles / pulsatility / vasomotion; segmentation is `setRegions` (interactive multi-ROI editor → `results.regionsMask`) → `runSegmentation` (fully automatic categorize + label + per-segment traces; `s.fNamesCopyTo` copies the segmentation onto co-registered siblings, replacing the old `assignCategories`) → optional `runDynamicSegmentation` (per-frame vessel diameter / flow) — `runVasomotion` writes the band-branched `results.vasomotion` tree per segment, and when `s.ppxVsmReturn` is non-empty also a LEAN per-pixel twin `results.vasomotion.ppx` — band-amplitude scalar `[Y×X]` maps plus an optional decimated `spectrum.amp`/`.phase` (`s.ppxVsmReturn` ∈ {`bands`,`spectrum`}); `runPulsatility` likewise writes the `results.pulsatility` tree per segment — `ps`/`pd`-prefixed markers + an `s.nHarm`-harmonic fit via the shared core `getPulsatilityMetrics` — and, when `s.ppxPulsReturn` is non-empty, the per-pixel twin `results.pulsatility.ppx`), plus `runRegistration`, `splitRegions`, the guided front-ends (`runGuidedContrast`, `runGuidedIntensity`) and the two myograph chains — pressure (`runMyographVideo` → optional `setMyographPresetIntervals` *or* `setMyographCrop` → `runMyographDiameter` → `setMyographIntervals` → `runMyographPropagation` / `runMyographVasomotion`) and wire (`runLabChart` → `setMyographIntervals` → `runMyographVasomotion`), all writing one `_MYO` triplet per recording; `setMyographIntervals` and `runMyographVasomotion` are literally the same steps for both and branch once on `source.modality`, and the `setMyograph…` steps open the shared interval editor once per recording |
| `Launchers` | Ready-to-edit example pipelines — the scripted way to drive the same steps |
| `GUIs` | Interactive apps — **three programs that share one session file**: `guiWorkbench` (**the Processing Workbench: start here**; it runs the LSCI pipeline and writes the session), `guiExport` (the **standalone** export tool: pick files, scan a folder or load a workbench session, choose parameters and how to average them, write one workbook per recording or one merged workbook for statistics — the interactive route to the same `exportToExcel` the launchers call) and `guiExplore` (the **standalone** results explorer: pick files, scan a folder or load a session, then plot and compare by experimental group / recording index / animal / recording type and export publication figures — and **the only place a myograph figure comes from**, since no myograph step writes a report page) |
| `GUIs/workbench` | The workbench's own components — the headless brain: `wbStepRegistry` (the step specs — **the linchpin**; adding a step or a modality is a data edit here), `wbDiscoverFiles`, `wbFileModel`, `wbTypeModel` (the animal / type / group label axes), `wbTypeSelection` (which steps each (type, product) row runs), `wbTypePresets` (the standard protocols), `wbPrereqs`, `wbStateEngine`, `wbSettingsModel`, `wbInvalidate`, `wbRefBranch` (which branch of an animal's reference a step takes), `wbRunRange` (the From/To rule), `wbExecutor` (the run loop — and the only place a step's actual branch products are resolved), `wbArtifacts`, `wbModalGuard` — plus `wbSession`, the versioned session file that is the **only** coupling between the three programs above |
| `Utilities` | Terminal consumers of finished results — `exportToExcel` (`exportToExcel(fNames)` writes the full workbook; the optional `exportToExcel(fNames,opts)` selects `opts.sheets` / `opts.format`, averages over labels or vessel type with `opts.groupBy` / `opts.weightByArea`, subsets `opts.columns`, and merges every file into one labelled workbook with `opts.merge` / `opts.labels` / `opts.outFile`. `GUIs/guiExport` is the interactive front-end for exactly these options) |
| `Simulation` | Synthetic dynamic-speckle generation (`getDynamicSpeckles`, `Launcher_speckleSimulation`) — self-contained |
| `3rd party` | External libraries (Bio-Formats, the ADInstruments LabChart SDK, superlets, …) — unmodified; origin, licence and version of each are in `3rd party/README.txt` |

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
3. Run `guiWorkbench` — the [Processing Workbench](#processing-workbench). Scan and label your recordings (**Files**), tick the steps each recording **type** runs (**Constructor**), press **Run** (**Process**), then press **Export…** or **Explore…** to open `guiExport` / `guiExplore` on the session it just saved.
4. Prefer a script? Copy a launcher from `Launchers/` to your own working location (**leave the originals unchanged**), set `libraryFolder`, `rootFolder` and `resultsFolder` in your copy, and run the `%% STEP` cells in order (`STEP 0` once per MATLAB session). See [Launchers](#launchers) and [Where results are written](#where-results-are-written).
5. Either way, read the header (comments at the top) and inline comments of each function you use — do not run steps blindly.

The processing steps you run are dictated by your protocol, not by the entry point: the workbench and the launchers call the same `run…` / `set…` wrappers and write the same `_d` / `_r` / `_s` files, so you can move between them freely. A pressure-myograph `.avi` runs in the workbench too, from raw video to `results.intervals(k).vasomotion`, including the interactive steps — the intervals, the time crop and the pre-set intervals each open one window per recording and then move on to the next.

### File-naming convention

Each processing step appends flags to the original file name, e.g.:

- `_t_e_K_d.mat` — (t) temporal contrast, (e) estimated epoch, (K) contrast, (d) 3-D data;
- `_c_BFI_r.mat` — (c) cardiac cycle, (BFI) blood flow index, (r) results.

**A product is named by its `_r.mat`.** Every step reads and writes the RESULTS file, while
only some of them open the 3-D data at all — so the file list you hand a step (`fNames`) is
always an `_r.mat` glob, and the step names the `_d` and `_s` members from it with
`getProductPath`. Most users only touch `_BFI_r.mat` (results) and `_BFI_d.mat` (3-D data)
files. File selection is often done with wildcards / regular expressions — see
`getFileNamesList`.

### Where results are written

By default everything is written **beside the recording it came from** — the `_d`/`_r`/`_s`
triplets, the `_rep_*.jpg` report pages and the PDFs assembled from them — and that is what
a project which never says otherwise gets. You can instead keep the results in a folder of
their own: name a **root folder** (where the recordings are) and a **results folder** (where
the results go), and the subpath is mirrored whole, so `<root>\A\B\Mouse1.rls` produces
`<results>\A\B\Mouse1_t_K_r.mat`. The raw recordings are then never written to at all.

In the **workbench** it is the second box on the *Files* tab (**results go in**); it follows
the folder you search until you point it elsewhere. In a **launcher** it is `resultsFolder`
in `STEP 0`, which starts out equal to `rootFolder`. Either way **only the step that reads
the recording has to be told** — every step after it is handed a product that is already in
the results tree, so nothing downstream has to know which arrangement you chose.
`getResultsPath` is the one place the rule lives.

Two things worth knowing before you point it elsewhere. A file you list by hand with its
full path is left where it is and writes beside itself, whatever the two folders say — that
is deliberate, since inventing a subpath for it would collide two same-named recordings
from two different disks into one folder. And **renaming a recording in the workbench
renames the recording only**: results already computed from it keep their old name and are
no longer linked to it. The confirmation dialog says so before anything moves.

### File format

Every product is written with `save(…, '-v7.3', '-nocompression')`. `-v7.3` is required for
the arrays this library produces (a `_d` cube routinely exceeds the 2 GB limit of the older
format); `-nocompression` turns off HDF5 deflate, which is **single-threaded** and dominates
the write. Measured on a 400 MB `_d` cube:

| | save | load | on disk |
|---|---|---|---|
| `-v7.3` | 10.3 s | 2.4 s | 352 MB |
| `-v7.3 -nocompression` | **1.2 s** | **1.3 s** | 427 MB |

**8× faster to write, ~2× faster to read, for 21% more disk** — and the read matters too,
since every later step loads these files again. Disk is the cheap resource here and
wall-clock during a 200-file batch is not, so the trade is taken *everywhere*, including on
the small `_r`/`_s` files: one rule covers every `save` in the library rather than a per-file
judgement call. If you are archiving a finished project to slow or metered storage,
re-compress there (load, then `save(…, '-v7.3')`) rather than changing the pipeline.

### Time base

`source.time` and `results.time` are **column** vectors, `[T × 1]`, in seconds — frames run
*down the rows*, matching `results.sData` / `dvsData` `[nT × nSeg]` and the `[nT × 1]` axes of
the pulsatility and vasomotion trees. Every producer writes that orientation
(`runContrastFromRLS`, `runInternalCycle`, `runExternalCycle`, `runBolus`, `runIntensity`), and
a consumer that needs a definite shape should still say `time(:)` rather than assume one.

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

`guiWorkbench` is the primary entry point: one window, three tabs, three questions —
**scan and curate the files**, **configure each recording type**, then **run and watch**.

Nothing is ticked per file. Processing choices are a property of the recording **type**,
not of the individual recording: you configure `BV` once and every `BV` in the project is
processed alike. Pressing **Run** *expands* that configuration into the ordered list of
(file, step) work items a launcher would have run, so a protocol of 200 recordings takes
the same handful of clicks as one of 6. State is read from the files on disk — *done*,
*stale* (a setting or an upstream re-run invalidated it), *queued*, *running*, *error*,
*skipped* — so a project you closed last week comes back with everything already done
still marked done, and **Run** resumes rather than repeats.

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

1. **1 - Files** — build and curate the working set. Point **look in** at a folder, put an
   extension or glob in **for files named** (`*.rls`, `*_t_K_r.mat`) and press **Scan** to
   recurse the whole tree; **Add files…** and **Add folder…** add more by hand without
   throwing the scan away. **results go in** is where everything the run writes is put; it
   follows **look in** until you change it, so by default the products land beside the
   recordings — set it elsewhere and your raw data is never written to (see [Where results
   are written](#where-results-are-written)). The five regexp boxes (**Animal**, **Type**, **Rec. index**,
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
   the configuration and the completion record in a versioned `.mat` sidecar, so a reload
   reproduces the table with no re-scan and long analyses survive a MATLAB restart.
   Beside them, **Export…** and **Explore…** make that session current and open
   `guiExport` / `guiExplore` **on it** — see [Three programs, one
   session](#three-programs-one-session).
2. **2 - Constructor** — say what each **type** runs. The tab asks two questions,
   because one raw recording can drive **two independent result sets**: the steps that
   read the recording itself (contrast → `_t`/`_s`, internal cycle → `_c`) each write a
   *new* triplet, and everything later appends to one of them.
   * **1 - Processing the recording itself** (top) — a `type × step` grid. Each tick creates that
     type's **product row** below. A `BP` that ticks both producers drives both
     pipelines from one recording, with no ambiguity to resolve by hand.
   * **Done once per animal** (beside it) — **registration** and **vessel typing** span the
     *animal*, not the type: they run over all of that animal's files whatever their
     type or group, against its pinned reference recording. Ticking one pulls its
     prerequisites into every type, since every type must produce what it reads.
   * **2 - Processing each result** (bottom) — one row per **(type, result)**, labelled with
     the type token and the product's stage flag (`_t`, `_s`, `_c`, `_e`). Tick the steps
     each pipeline runs. A cell is greyed with its reason when the step does not belong
     on that branch (vasomotion is contrast-side, pulsatility cardiac-side), or shown as
     *inherited* when the wrapper covers the whole recording in one call — region
     drawing, segmentation, split regions and BFI are drawn once and apply to every
     product, so they can never be on for one and off for another.
   * **Settings** (right) — pick a **Type** and a **Step** to edit its parameters, split
     into *Basic* (what a protocol actually tunes) and *Advanced*. Settings are keyed by
     **(step, type)**, so both rows of a type share them and a divergent animal is simply
     a second type. Values shared between steps propagate automatically, and an edit
     marks that step and everything downstream *stale* so you cannot silently mix
     parameters. **Copy configuration from** seeds a type from a standard protocol
     (Pulsatility, Vasomotion, NVC, Vasoreactivity, Pulsatility+Vasomotion) or from
     another type you already configured — a starting point, never a locked mode.
   * **What will run** at the bottom is where status is read: how many files each
     type covers, which reference *file* each per-animal step resolved to, and the
     warnings.
3. **3 - Process** — read-only. Choose how much to run, press **Run**, watch.
   * **From / To** slice the run by pipeline step, so you can take a protocol to
     segmentation, tune it, and carry on without unticking anything. **`From = Last
     valid`** *resumes*: it starts at the first configured step that is not finished and
     skips whatever is already on disk. Naming a step instead **forces** it: everything
     from there to **To** runs again, done or not — that is how you re-run segmentation
     with new settings, or re-open the vessel-type painter on an animal you have already
     typed. **To** offers only the steps at or after **From**.
   * **Preview order** lists exactly what would run, in execution order, without calling
     anything. **Run** executes it **step-major, the way a launcher does** — every
     recording is contrasted, then every one segmented — so the steps that read across a
     file set (registration, vessel typing, split regions) always find the whole set at
     the same level. **Stop** cancels cooperatively between files; a failed step marks
     its cell as an error and the batch continues.
   * **Two state tables**, mirroring the Constructor: *Recordings* on top (what the raw
     file itself runs, which goes first) **beside the run log**, and *Products* below
     (one row per `(recording, product)` the configuration creates — a product row
     appears as soon as its raw producer is ticked, long before the file exists). Each
     cell reads `·` / queued / running NN% / done / error / skipped. Select a row to see
     its path, its labels and its last error.
   * **Wipe all** deletes every `_d.mat` / `_r.mat` / `_s.mat` belonging to the
     recordings listed on the Files tab, so a protocol processed with the wrong settings
     can be started over in one action. The raw recordings are never touched, nothing
     outside the listed recordings is touched, and the count is shown for confirmation
     first. It cannot be undone.
   * **Reports** (right) — the report images the run wrote, newest first; double-click to
     open one in the OS viewer, which stays zoomable while processing continues. Tick
     **Create PDF reports** and each step's images are appended into **one PDF** when
     that step finishes — a 60-file column becomes one document to page through instead
     of 60 links. The dropdown filters the list to images or PDFs.
   * Interactive steps (region drawing, vessel typing, tree editing) open their own
     editor and grey the workbench out until you finish.

**Next** and **Exit workbench** sit in the same bottom-right corner of every tab (the
last tab has no *Next*). **Exit is one action**: it requests the same cooperative stop
the *Stop* button does, saves the session, and closes as soon as the current step lets
go. The session is rewritten on **every** state change — a scan, a label, a Constructor
tick, a settings edit, each finished cell — always to **`workbench-sessions/lastSession.mat`**
beside the library (a gitignored folder), and additionally to wherever you last saved or
loaded one, so there is always something to resume from without anyone pressing *Save*.

The workbench covers the **LSCI** pipeline. It never reimplements any processing — it
orders, gates, parametrises and calls the same wrappers a launcher would. It also does
**not** export or plot: those are two separate programs.

### Three programs, one session

Export and Explore are no longer tabs. They are standalone tools, and the **session file
is the only thing that travels between them and the workbench**:

| program | what it does | how it gets its files |
|---|---|---|
| `guiWorkbench` | processes recordings; **writes** the session on every state change | scan a root, add files by hand, or load a session |
| `guiExport` (`GUIs/guiExport.m`) | choose parameters and averaging, write one workbook per recording or one merged workbook for statistics | pick files, scan a folder, **or load a session** |
| `guiExplore` (`GUIs/guiExplore.m`) | plot and compare finished results, export publication figures | pick a file or folder, **or load a session** |

Pressing **Export…** or **Explore…** in the workbench's Files tab saves the session and
opens the tool on it: the file list arrives with each recording's **animal**, **type**,
**experimental group** and **index** already resolved, and with the record of which steps
actually ran — so the explorer offers only the plots the data can support. The workbench
keeps no handle to either tool; closing one leaves it running, and both open perfectly
well on their own:

```matlab
guiExport                          % or guiExport('<path to>/workbench_session.mat')
guiExplore                         % or guiExplore('<path to>/workbench_session.mat')
```

---

## Launchers

The launchers are the **scripted alternative** to the workbench — fully supported, and the
right choice when you want a reproducible record of a pipeline, a headless/batch run, or a
starting point to modify. Copy a launcher to your own location and edit it for your project.
`STEP 0` is where the three folders are named — `libraryFolder`, `rootFolder` (your
recordings) and `resultsFolder` (where the results go, equal to `rootFolder` until you
change it; see [Where results are written](#where-results-are-written)) — and it is run once
per MATLAB session.

They are also the **reference implementation of every pipeline**: the workbench's step
specs (`GUIs/workbench/wbStepRegistry`) are transcribed from these files — the same
wrappers, the same `%Example-of-s` defaults, the same file lists and the same order — so
if the two ever disagree, the launcher is right. Reading the one that matches your
protocol is the fastest way to understand what the workbench is doing on your behalf.

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
| `Launcher_myograph` | Pressure myograph (`.avi`) end to end, with the wire-myograph (`.adicht`) variant commented under each step. |

A pressure-myograph `.avi` and a wire-myograph LabChart `.adicht` are also processed in
the Processing Workbench (`GUIs/guiWorkbench`) like any other recording — pick the
**Pressure myograph** or **Wire myograph** protocol in the Constructor to tick the steps
each one needs. The launcher and the workbench call the same wrappers and write the same
`_MYO` triplet.

### What a myograph run leaves behind

**Text and result files, not report images.** Every step narrates three lines — Starting,
Writing results, Finished — to the command window and the workbench log, and writes to the
recording's `_MYO` files, wherever [the results go](#where-results-are-written). It writes
**no `_rep_*.jpg` pages**, deliberately: a myograph
check is worth looking at properly rather than as a thumbnail. A user who has processed
LSCI first should not go looking for pages that were never written.

- **The figures come from `guiExplore`.** Including the two detection checks a report page
  would have carried — the **diameter map** (the diameter along the vessel over time, which
  says whether detection held) and the **detected walls** over a frame (which says whether
  it found the right edges) — plus diameter traces in pixels or % of baseline, the
  individual line traces, diameter statistics, propagation scalars and the per-line lag
  behind them, and every vasomotion view.
- **The numbers come from `guiExport` / `exportToExcel`.** A bare `exportToExcel(fNames)`
  writes a myograph recording up to nine sheets, chosen from the data: `settings`,
  `comments` (wire myograph only — the LabChart operator log), `intervals`, `propagation`,
  `vasomotion`, `ampPct`, `spectra`, `ampPctSpectra` and `diameterTraces`.

---

## Exploring processed results

Follow any pipeline to a `_BFI_d.mat` file (the 3-D data `X × Y × Time` plus a time vector,
written beside the `_BFI_r.mat` its steps were given)
and explore it with basic MATLAB skills — ROI selection, filtering, plotting as image or
video.

For the finished `_r.mat` results (per-segment BFI / diameter / pulsatility / vasomotion,
single files or group comparisons), use **`guiExplore`** (`GUIs/guiExplore.m`). It needs no
workbench: point it at a file, at a folder (labelled by three regexp boxes), or at a saved
workbench **session**, from which it takes the file list together with each recording's
experimental group, animal, recording type and index. Those four axes are independent — an
animal may span groups and a group may span animals — and any of them can be the x-axis or
the colour. The explorer matches the plot type to the data (time series, box plots,
spectra, spectrograms, maps) and exports the figure at a chosen DPI. Hand-made groups
(select files, name them, **Create group**) override both the regexp and the session, and
survive a re-scan.

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
