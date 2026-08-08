# Working rules for this library

MATLAB R2025b. Vascular Dynamics Analysis: LSCI / DLSI contrast pipelines, the
Processing Workbench, and the Myograph subsystem. Layers are
`Core/` → `Wrappers/` → `Launchers/` + `GUIs/`, with tests in `claude-tests/`
(git-excluded — **Grep respects `.git/info/exclude`, so pass that path explicitly**).

---

## 1. Testing effort is budgeted

**No single test invocation may exceed 20 minutes. Target 2–10.**

A run heading past 20 minutes is **abnormal**. It is a defect in the *test* — an
oversized fixture, a pin being run as a per-edit check — not a fact about the library.
Split it, shrink its fixture, or ask. **Never sit through it.** Nothing in this library
requires it: a user processes a whole recording in under 20 minutes.

**Test depth is proportional to the change.**

| Tier | What it is | Cost |
|---|---|---|
| **T0 — static** | `checkcode` on the changed files. No MATLAB run. | seconds |
| **T1 — targeted** | The *one* test file covering the changed surface. | < 2 min |
| **T2 — module fast** | `test<Module>Fast` | 2–6 min |
| **T3 — exhaustive** | `test<Module>Exhaustive` — pins, goldens, pools. **Announce first.** | 15 min+ |

### Where the common cases land

| Change | Tier |
|---|---|
| GUI label, user-facing string, layout tweak | **T1** — the module's GUI test. *Not* a suite. |
| Comment, header, docs | **T0** |
| A new wrapper; a registry entry | **T2** |
| A core numerics change | **T1** (that core's golden) **+ T2** |
| A refactor guarded by a bit-identity pin | **T3**, announced |
| Anything touching `parfor` / `wbPool` | **T3** with `includeParfor` |

Escalating a tier costs one line saying why. Dropping a tier costs nothing.

**Running everything is never the default.** Ask first. And never re-run a suite that
was already green just to print a tidier summary — you already have the result. Three
sweeps where one was needed is the single most expensive mistake made in this repo.

---

## 2. Announce anything over ~2 minutes

Before a MATLAB run you expect to take more than ~2 minutes, say three things:

- **what** is running,
- **which tier, and why that tier**,
- **how long** you expect it to take.

Then surface progress *while it runs*. **An hour of silence is a defect in its own
right**, independent of whether the result was green. `runTestList` prints each test
before it starts and carries a running total for exactly this reason.

---

## 3. Test entry points

| Module | Fast (per-edit) | Exhaustive (on purpose) |
|---|---|---|
| Workbench | `testWorkbenchFast` | `testWorkbenchExhaustive([includeParfor])` |
| Myograph | `testMyographFast` | `testMyographExhaustive` |
| LSCI | `testLSCIFast` | `testLSCIExhaustive` |
| EPFL (fluorescence) | `testEPFLFast` | `testEPFLExhaustive([includeParfor])` |
| Reporting | `testReportingFast` | — whole module is ~2 min |
| Registration | run the two files directly (~12 s) | — |
| Explore | inside `testWorkbenchFast` | `testExploreExhaustive` (~8 min, **T3**) |
| Launchers | `testLaunchers` (~0.3 s, in the EPFL exhaustive list) | — |

**`testEPFLFast` is a RUNNER over eight files, and it stopped being one of them at S13.** The
fluorescence per-edit path grew one file per session from S1 to S10 and each one earned its
place, but the largest of the eight also wore the runner's name — so there was no way to say
"run the per-edit path" that did not also mean "read two gigabytes". Its 204 checks are
`testEPFLSpine` now.

**They are listed cheapest first, and that is the shape to keep.** The first three together are
3.6 s and cover eight cores and two wrappers end to end, because none of them opens a recording;
the last is 90–118 s because it reads the 2.1 GB `.cxd` twice. **A new claim belongs in a file
that opens nothing, by default.**

| file | s | what only it covers |
|---|---:|---|
| `testBolusMetrics` | 1 | the five bolus cores, synthetic |
| `testTopologyMetrics` | 1 | `getTopologyMetrics`, synthetic |
| `testVesselTypesTree` | **1.6** | `getVascularCues`/`defaultFlowParams` over all five products, the artery/vein guess **with its negative control**, and `setVascularTree` headless on an intensity product. Reads no recording, which is why it is the cheapest file in the module |
| `testBackgroundRemoval` | 16 | needs a **bolus** product and a crop small enough to clean every frame of — its own 59 MB fixture through `makeBackgroundFixture` |
| `testIntensityPulsatility` | 21 | `runIntensityPulsatility` against **closed-form** markers, the `pv*`/`wm*`/no-`ps*` prefix guard on one table, and `runVasomotion`'s widened gate on an `_a_I` |
| `testMotionEnhancementStep` | ~40 | `runMotionEnhancement` end to end on `mePhantom` products |
| `testRunCTTH` | 45 | the transit-time wrapper on a real bolus product |
| `testEPFLSpine` | 90–118 | the `_a_I`/`_c_I` name model, the EPFL registry block **by name and in order**, the protocols, `runIntensity` and `runIntensityInternalCycle` end to end on the small `.cxd`, their pre-read refusals, polarity, and the `intensityTwinOf` contract |

**`testEPFLExhaustive` adds `testLaunchers` and `testEPFLReal`** — the claims that had only ever
been made on phantoms: a real `_c_I` folded out of the module's own recording through both
cardiac metric steps, `runTopologyAnalysis` as a wrapper (gate, page, mixed-settings warning),
the fluorescence hierarchy derived from a **measured** arrival, `setVesselTypes`' three widgets
through the `uiwait` it blocks on, and `runCTTH`'s per-pixel path at 1.68 M pixels.
`testEPFLExhaustive(true)` appends `testParforToggles`, whose `ctth` axis is this branch's.

They write to `tempdir`, never into the repository. `claude-docs/test-cost-inventory.md` has the
measured cost and the reasoning behind each split.

The explorer has no fast runner of its own: `testExploreTool` plus the `('synthetic')`
mode of `testExploreIndex` and `testExplorePlan` live in `testWorkbenchFast`.
`testExploreExhaustive` adds their real-set parts and `testExploreCascade`, and **needs
the reference data at `C:\Dropbox\Work\Data`** — those tests fail rather than skip
without it, on purpose.

`claude-docs/test-cost-inventory.md` holds the measured wall time of every test.
Consult it before running something you have not run before.

---

## 4. Fixtures

- **Never render a fixture inside a test.** Synthetic myograph recordings go through
  `myographCachedVideo`, which keys the cache on the *parameters* — change a parameter
  and it re-renders; change nothing and it copies. A guard like `if ~exist(vp,'file')`
  keys on the file *name* and will silently serve a stale recording after a parameter
  edit.
- **Keep test data small enough to prove the point and no larger.** Ten minutes of
  recording is plenty; two hours is not a test, it is a bill. If an assertion needs a
  long recording to hold, that is worth knowing and worth writing down.
- Use `-nocompression` on `save -v7.3` for large test artifacts.

---

## 5. The test tree does not grow by default

Every file under `claude-tests/` earns its place, and **scaffolding is retired in the
same session as the refactor it served** — not left for someone to find later. A test
suite that only ever grows stops being a safety net and becomes a tax.

A test or its scaffolding is **obsolete** when any of these is true:

- **The refactor it guarded has landed and stayed green.** A before/after pin is
  evidence for one change. Once that change is settled the pin is spent.
- **It cannot fail.** If it compares today's code against a golden it regenerates from
  today's code, it asserts nothing.
- **Its fixture no longer exists**, so it skips or throws on every run.
- **Its default arguments no longer match its golden** — it reports a library fault that
  is really a stale argument.
- **It is a one-off `bench*` or `check*` script** from a session that has closed.

**Never give a frozen copy the real function's name.** `setLibraryPath` now demotes the
`claude-*` trees to the end of the path, so a library file always outranks a tooling
copy, and its second output names any collision:
`[dropped,shadowed] = setLibraryPath(root)`. That is a backstop, not a licence — a
harness that shadows a wrapper is still a bug waiting for the one code path that adds
`claude-tests` itself. Name frozen copies distinctly (`runInternalCycleBefore.m`), never
`runInternalCycle.m`.

**`claude-tests/` is git-excluded, so deletion is irreversible.** Copy anything you
remove to a scratch folder outside the repo first, and say in your report what was
removed and why.

Adding a test means adding its measured cost to `claude-docs/test-cost-inventory.md`.
A test with no recorded cost cannot be tiered, and anything untiered ends up in the
per-edit path by accident.

---

## 6. Gotchas that have already cost time

- **`matlab.exe` on Windows returns immediately** unless given `-wait`; and
  `Start-Process` with redirected stdout makes it exit without running anything. Launch
  batch MATLAB through `cmd`: `matlab -wait -batch "…" > log 2>&1`.
- **Check for a running `MATLAB.exe` before starting one.** The author's session may be
  busy; do not attach to a shared session, it will queue behind their work.
- **PowerShell `-match` is case-insensitive.** Parsing test logs for `FAIL`/`SKIP` with
  it matches lowercase words inside check *names* and invents failures. Use `-cmatch`.
- **MATLAB integer arithmetic saturates, it does not wrap.** A `uint64` hash pins at
  `intmax` after a few characters and every key collides. Hash in `double`.
- **`VideoWriter` appends the profile's extension.** An `outPath` of `x.avi.partial`
  produces `x.avi.partial.avi`.
- `addpath(genpath(root))` sweeps `.claude/worktrees/*` onto the path and **shadows the
  library**. Always use `setLibraryPath(root)`, which drops `.claude`, demotes the
  `claude-*` trees below every library folder, and reports name collisions.

---

## 7. House style

- Every file carries the library header block (H1 summary, syntax, inputs/outputs,
  example, `See also`, author/copyright) and the `%------------- BEGIN CODE --------------`
  marker.
- User-visible strings are written for biologists: no bracketed explanations, no
  wrapper or settings-field names, registry labels rather than ids.
- Retire a file by renaming it `OBSOLETE<name>.m`, never by deleting it.

- **THE CUT IS ALWAYS CLEAN. No code in this library exists to read a file it no longer
  writes.** When a shape changes, the already-processed files are migrated or re-processed
  — they are not a constraint on the design. So: no fallback that reads an old field name,
  no second spelling of a setting, no "a file written before X was added is read off Y".
  Every such fallback is a second way for a file to disagree with itself, it can never be
  exercised once the data has moved, and it quietly outlives the reason it was written.
  When a revision is planned, migrating the data is part of the revision.

  Two things this rule does **not** cover, because they are not back-compatibility:
  a **setting** that produces a smaller product on purpose (`s.keepArrays false`) is a
  live choice and the code that reads its shape stays; and a field that is legitimately
  in one of **two current states** (a wire channel's samples, whole before the intervals
  step and inside the windows after it) has one reader for both, by design.

- Data model: `_d` / `_r` / `_s` triplets — **except the myograph, which is the PAIR
  `_MYO_r` + `_MYO_s`.** A myograph recording is not changed by the analysis and
  re-reading it is seconds, so the product describes the recording (`results.recording`)
  and keeps the measurement (`results.intervals(k)`), and copies neither. The speckle
  branch still has three members and that is not changing — sweep on `_MYO_d`, never on
  `_d.mat`. Wrappers are `run<X>` / `set<X>`, cores are
  `get<X>` / `fit<X>`.
