# MotionEnhancement

Make pulsatility- and vasomotion-driven **vascular wall motion visible** in fluorescence
vascular recordings. Output AVI.

**This is a library module.** It lives at `Core/MotionEnhancement/` and its wrapper is
`Wrappers/runMotionEnhancement.m`, registry step *Wall motion*, which is what a user runs.
The wrapper reads a `_c_I` or `_a_I` product, places perpendicular cuts along every segment
the library's own segmentation found, measures both walls on every cut with `meKymograph`,
and writes eight per-segment columns and their evidence; the magnified AVI below is its
optional illustration, off by default. The stand-alone entry points — `mePassA`, `meProbe`,
`meRun` — stay, because a 64 GB recording is still read once by hand and a session that
wants to look at four cuts should not have to run a step.

The dependency used to be one-way and asserted by a grep. That grep is gone, retired in
the session that promoted this folder: the block was stand-alone *"until tested"*, it is
tested, and the library now calls it on purpose. What is still asserted is that no `me*`
name collides with anything else on the path (`testMotionEnhancementPins` check 8).

The one file the stand-alone path writes beside the videos is `<recording>_ME.mat` — **not**
a `_r.mat`/`_s.mat` triplet, because those suffixes are swept by the workbench and the
explorer and it is not a library product. The wrapper writes into the product instead.

Method: **retrospective self-gating, then phase-based magnification on a Riesz pyramid**
(Wadhwa et al., ICCP 2014), written from the published pseudocode. Amplified MRI amplifies a
cardiac-gated cine, never a raw series — that is the part of the method that matters most
here, and it is stage one.

---

## What it produces

Four products, and **every one ships with its own matched control**: the same run with the
timing destroyed and nothing else changed — same frames, same texture, same photon noise,
same α. A magnified movie is persuasive whether or not it is true, so the deliverable is
never the movie. It is the movie beside the ratio between it and its own control.

| product | mode | α per level (λ = 16, 32, 64 recording px) | ÷ its own control |
|---|---|---|---|
| `<rec>_ME_puls_cine.avi` | cardiac, cycle-averaged, 38–75 µm calibre | `[53.1 106.2 212.5]` | **5.3–19.1** |
| `<rec>_ME_puls_cine_smallVessels.avi` | cardiac, cycle-averaged, ~38 µm calibre | `[373.8 747.7 1495.3]` | **3.6–28.2** |
| `<rec>_ME_puls_cont.avi` | cardiac, continuous | `[53.1 106.2 212.5]` | **1.0 — not detected** |
| `<rec>_ME_vaso_cont.avi` | vasomotion, continuous | `[11.2 22.5 44.9]` | **0.4–0.8 — not detected** |

**One α serves one calibre.** The output is proportional to the displacement, so the
cardiac cine is delivered twice and each says on its own frames which vessels it is valid
for. In the small-vessel run the large vessels are 6.7–7.8× past their λ/4 bound and their
wall edge comes out 30 % and 83 % wider — measured, visible, and written on the frame.

---

## How to run it

```matlab
s = meSettings;
s.outFolder = 'D:\work\motionEnhancement';   % ~15 GB per recording. Not a synced folder.
mePassA(s, 'D:\data\rec.cxd');               % ONCE per recording. ~23 min.
probe = meProbe(s, 'D:\data\rec.cxd');       % ~40 s. Finds the vessels and the rhythms.
edit meRun                                    % then Run Section, cell by cell. ~4 min.
```

`mePassA` reads the recording **once** and writes two flat memory-mappable uint8 working
copies (2×2-binned whole field, full-resolution crop) plus a sidecar. Everything after it
reads those. On the reference recording Bio-Formats spent **16 min 43 s indexing the CXD
before the first plane**, which is the entire reason the working copies exist. Do not use
`readCXD` on a file this size — it still preallocates the whole stack in RAM.
`runIntensity` was rewritten to stream and no longer does.

`meRun` is the entry point: a launcher-style cell script. Edit the recording, the folder and
the two calibres at the top; the rest follows.

---

## The artifact floor — the number every claim has to clear

Measured on the reference recording, 2026-08-07, as an **input-equivalent wall displacement
in pixels of the recording**: the input that would come out at twice the control's own
output. Directly comparable with the calibration phantom's numbers.

| what the floor is measured on | floor, px per wall |
|---|---|
| photon noise alone, at the cine's own averaging | **0.0007** |
| **the delivered cardiac cine** — its own matched control | **0.0008 to 0.031, median 0.0042** |
| the same cine folded on a period no physiology occupies | median **0.0053** |
| **the continuous mode**, no cycle averaging | **0.33 to 0.54** |

Two things follow, and they are the whole reason this block ships controls:

- **The cine's floor is six to eight times photon noise.** What limits this measurement is
  vessel texture, red cells crossing the wall and residual bulk motion — not photons. More
  averaging buys less than √N suggests.
- **The continuous mode's floor is 5 to 100 times the wall motion it is looking at.** Its
  product is exactly its own control, and it is the most convincing-looking of the four
  videos.

**The measured wall motion, for scale:** 0.003 to 0.083 px per wall (0.008–0.21 µm), i.e. a
pulsatility of 0.5–2 %, on vessels of 37–75 µm. At 2.5 µm/px.

**The smallest calibre with a signal above its own control is 36.9 µm** — the delivered
cardiac cine clears its control there by 7.7×, and by 15.8× when only λ = 16 recording px is
amplified. The **raw** wall measurement does not reach below **38 µm**; below that it is the
folding and the filtering that see the vessel, not the estimator.

---

## What it must not be used for

- **Not as a measurement.** The magnified movie is an illustration. The measurement is
  `meKymograph` on the raw stack, folded on the gate — and it is in `<rec>_ME.mat`.
- **Not without its matched control.** Every product has one, the same size, meant to be
  watched beside it. A reader who cannot tell them apart has been told something true.
- **Not for vasomotion on a recording like this one.** There is no vasomotion cine —
  0 of 49 cycles survive, because a 5.7 s cycle cannot avoid an animal that moves every
  2.3 s — and no vasomotion wall dilation was detected: 0.4 to 0.8 of its own control, and
  what in-band motion there is moves the vessel's *centre* as much as its half-width, i.e.
  it is the preparation translating, not the wall changing calibre. α for that product comes
  from an **upper bound**, not from a detection.
- **Not at a calibre other than the one burnt on the frame.** Past λ/4 the wall deforms
  rather than moves.
- **Not for pulse-wave propagation.** No row-wise phase ramp is resolved on the gated cine
  (the interval spans zero), and magnifying does not add one — but this block was not built
  to measure propagation and a row regression over a vascular network is not a velocity.
  `setVascularTree` is the library's tool for that.
- **Not without a calibration.** Every micrometre in this document is at 2.5 µm/px, measured
  once for one rig. `runMotionEnhancement` refuses to run until somebody has given it a
  pixel size, and it ships with no default, for the same reason.
- **Not as a fence.** The wall mask says where the shift is *applied*, not where it *lands*.
  Measured on the delivered cine, the off-mask magnification decays by **×140 over 64
  recording px** — exactly λ at the coarsest amplified level — so **amplification is confined
  to within about 160 µm of a wall** and no further in. Quote that radius, not a ratio: the
  ratio is ×13.7 immediately outside the mask contour and ×1 900 beyond 64 px, and both are
  the same mask. On a dense field it is ×3.1 immediately outside, because almost everything
  there is within 64 px of a vessel.
---

## The files

| | |
|---|---|
| `meSettings` | every default in one place |
| `mePassA` `meReadRaw` | the CXD, once, into two flat working copies; the one reader of them |
| `meProbe` | size, rate, bleach, spectrum, photon SNR, wall width and wall displacement, with its own null |
| `meEpochs` `meGate` `meCine` | when the animal was still; the beats and their rejection; one averaged cycle + its SEM |
| `meEnhance` | **the two modes** — what the magnifier is run on, and every derivation returned in `info.settings` |
| `meRieszPyramid` `meExpand` `meCollapse` `mePhase` | the representation |
| `meTemporal` `meGlobal` `meBlur` `meMask` | zero-phase band-pass; bulk-motion removal in the phase domain; amplitude-weighted blur; the wall mask |
| `meShift` | **the magnifier** — the pipeline order lives here, and its second syntax re-shifts a cached analysis so an α sweep is free |
| `meLinear` | linear EVM, kept as the noise-magnifying baseline. Not a fallback |
| `meKymograph` | a cut across a vessel, and its two walls to a fraction of a pixel — the estimator every claim rests on |
| `getWallMotion` | **the per-segment layer over it**: every cut of one segment measured, the three cut rejections, the reduction as a vector mean, and the amplitude in either of two modes. Named `get*` rather than `me*` because it is a library core with a wrapper above it. Its partner `getSegmentCuts` is in `Core/Vasculature/`, because placing perpendicular cuts along a centre line is geometry and has a second caller waiting |
| `mePhantom` `meValidate` | known displacement in closed form; the calibration curve. `meValidate` is ~10 min and is **not** a test |
| `meWriteVideo` `meRun` | the AVI, with α and the calibre burnt on every frame; the entry point |

## Tests

| | s | |
|---|---:|---|
| `testMotionEnhancementFast` | 20 | phantoms only, no data on disk. Thirteen checks. **T1** |
| `testMotionEnhancementPins` | 19 | the calibration's own numbers pinned by value, plus the name-collision guard. **T1** |
| `testMotionEnhancementReal` | 155 | every wrapper run once on the working copies; **fails rather than skips** without them. **T2** |

`meValidate` (~10 min) and `meRun` (~4 min) are not tests and must not enter the per-edit
path.

---

## Two things that will bite

- **`s.levels` names a physical scale, not an index.** A pyramid level carries a wavelength
  in pixels *of the array it is built on*, and the two working copies do not share one.
  `s.levels` is stated in levels **of the recording** and `meEnhance` maps it onto whichever
  copy it reads; `info.lambdaFull` records the wavelengths on every run. Taking `3:5`
  literally on the binned copy amplifies twice the intended scale and returns 2.1× less,
  with no error and no warning.
- **The gain is `1 + η·α`, not `1 + α`.** A binomial Laplacian pyramid's bands overlap and a
  level that is not amplified holds the edge back. η is a property of the **structure** as
  much as of the code: 0.36 on the calibration phantom, **0.13–0.20 on real walls** at the
  same single α. Plan with the real ones.

## Licensing

The methods are patented by MIT (US 9,338,331; US 2014/0072190) and MIT's released *code*
is non-commercial-research only. This is an independent implementation from the published
papers and pseudocode, so no MIT-licensed source is in this repository and no licence rides
along. The patents still nominally read on the method. Recorded here so the author decides
rather than discovers.
