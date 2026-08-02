# Third-Party and External Scripts

## Important License Information

The files, scripts, and libraries contained within this directory are **not** original to this project and are **not** governed by the main `LICENSE` file (CC BY-NC-SA 4.0) located in the root of this repository.

These components were developed by third-party authors and are distributed under their own, original software licenses. This `README.md` file is provided for transparency and to ensure compliance with all third-party license agreements.

Users of this project must adhere to the specific license terms of each individual script. Please refer to the original source repositories or the included license files for full details on permitted use, modification, and distribution.

## What is vendored here

Nothing in this folder is edited. A fix belongs upstream; a local change would be
silently lost the next time a package is refreshed, and would make the version
recorded below a lie. The library routes *around* these packages: exactly one
library function touches each of them, and that function is named below.

### `bfmatlab` — Bio-Formats for MATLAB

- Origin: Open Microscopy Environment, https://www.openmicroscopy.org/bio-formats/
- Licence: GPL v2 (see the Bio-Formats distribution).
- Version: not recorded when it was vendored; `bfCheckJavaPath` reports the
  version of the bundled `bioformats_package.jar` at runtime.
- Reached from: `Core/Read files/readCXD.m` (Hamamatsu `.cxd` recordings).

### `adinstruments_sdk_matlab` — LabChart `.adicht` reader

- Origin: Jim Hokanson, https://github.com/JimHokanson/adinstruments_sdk_matlab
- Licence: MIT (`adinstruments_sdk_matlab/LICENSE`, Copyright (c) 2014 Jim Hokanson).
- Version: commit `2666f38e242ad3d57f012576c1a26328f54ea43b`, 2024-07-27, vendored
  2026-08-02. Upstream publishes no tagged release, so the commit IS the version.
- Platform: **Windows and 64-bit MATLAB only.** The package is a MEX wrapper
  (`+adi/private/sdk_mex.mexw64`) over ADInstruments' own `ADIDatIOWin64.dll`, both
  of which ship inside it. There is no source-only fallback and none is wanted:
  `.adicht` is a closed format and this DLL is the only thing that reads it.
- Reached from: `Core/Read files/readLabChart.m` (wire-myograph LabChart recordings).
  That is the ONLY place in this library that names `adi.*`.
- On the path: the package folder `+adi` must **not** be added to the MATLAB path -
  its PARENT (`3rd party/adinstruments_sdk_matlab`) is what makes `adi.readFile`
  resolve. `Utilities/setLibraryPath` gets this right without doing anything special,
  because `genpath` skips `+package` folders; do not "fix" that by adding `+adi`.
- Also bundled by upstream and harmless: `files/blank_labchart_8_file.adicht`, a
  24 KB empty LabChart 8 file used by its own tests, and two format PDFs.
