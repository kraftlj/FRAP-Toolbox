# FRAP-Toolbox User Guide

This is the main GitHub user guide for FRAP-Toolbox. It modernizes the original
2014 supplemental user guide while preserving the legacy MATLAB workflow,
original test-data recipes, and supplemental tables used for parity testing.

Use this guide when you want to run FRAP analysis with the modern Python tools,
understand which model to choose, or reproduce the original MATLAB examples.
For implementation details, see [frap_toolbox_py/README.md](../../frap_toolbox_py/README.md).
For modernization parity notes, see
[docs/user-guide-parity-testing.md](../user-guide-parity-testing.md).

## Contents

- [Overview](#overview)
- [Citing FRAP-Toolbox](#citing-frap-toolbox)
- [Which Workflow Should I Use?](#which-workflow-should-i-use)
- [Modern Python Quick Start](#modern-python-quick-start)
- [Running Analyses](#running-analyses)
- [Supported Image Formats](#supported-image-formats)
- [Data Availability](#data-availability)
- [ROI Guidance](#roi-guidance)
- [Results And Exports](#results-and-exports)
- [Legacy MATLAB User Guide](#legacy-matlab-user-guide)
- [Original Test Data Workflows](#original-test-data-workflows)
- [Troubleshooting](#troubleshooting)
- [Supplemental Tables](#supplemental-tables)
- [Supplemental Figure Legends](#supplemental-figure-legends)
- [References](#references)

## Overview

FRAP-Toolbox analyzes Fluorescence Recovery After Photobleaching experiments.
The toolbox includes diffusion and reaction-model workflows for extracting
diffusion coefficients, mobile fractions, and kinetic fit parameters from raw
microscopy image sequences.

The repository now contains two user paths:

- A modern Python package with a Streamlit app, command-line interface, BioIO
  image loading, export bundles, ROI-mask support, and parity tests.
- The original MATLAB source archived under `legacy/matlab/`, preserved for
  reproducibility and historical comparison.

The original 2014 PDF was written for the MATLAB application. This Markdown
guide keeps those details but presents the Python workflow first for new
analysis work.

## Citing FRAP-Toolbox

**FRAP-Toolbox: Software for the analysis of Fluorescence Recovery After
Photobleaching**

Authors: Lewis J. Kraft, Jacob Dowler, Charles A. Day, Minchul Kang, and Anne K.
Kenworthy.

If FRAP-Toolbox supports your work, please cite:

```text
Kraft LJ, Dowler J, Day CA, Kang M, Kenworthy AK. (2014). FRAP-Toolbox: Software for the analysis of Fluorescence Recovery After Photobleaching. https://github.com/kraftlj/FRAP-Toolbox (accessed Month Day, Year).
```

Correspondence: Anne K. Kenworthy, Department of Molecular Physiology and
Biophysics, 718 Light Hall, Vanderbilt University School of Medicine, Nashville,
TN 37221, [Anne.kenworthy@vanderbilt.edu](mailto:Anne.kenworthy@vanderbilt.edu).

## Which Workflow Should I Use?

| Workflow | Use it when | Entry point |
| --- | --- | --- |
| Modern local app | You want an interactive browser UI for local datasets, ROI selection, fitting, and exports. | `frap-toolbox-app` |
| Modern CLI | You want repeatable scripted analysis, batch processing, or CI-friendly validation. | `frap-toolbox` |
| Legacy MATLAB app | You need to reproduce the original 2014 workflow or compare against MATLAB-era behavior. | `legacy/matlab/Main_GUI.m` |

For new work, start with the modern Python app or CLI. Use the legacy MATLAB
workflow when reproducing the original guide, investigating parity differences,
or validating historical exports.

## Modern Python Quick Start

The Python port supports diffusion, Reaction 1, and Reaction 2 fitting. Commands
below assume macOS or Linux with `bash`; on Windows, use the equivalent virtual
environment activation command.

### 1. Create An Environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

### 2. Install The Package

For test and development work:

```bash
python -m pip install -e ".[test]"
```

For the local browser app with common ND2 support:

```bash
python -m pip install -e ".[app-nd2,test]"
```

Optional extras:

- `app` installs the local Streamlit app.
- `nd2` installs Nikon ND2 support through BioIO.
- `bioformats` installs Bio-Formats-backed readers.
- `legacy-io` installs the older AICSImageIO fallback.
- `app-full` installs the app plus ND2, Bio-Formats, and legacy fallback I/O.

### 3. Run The App

```bash
frap-toolbox-app
```

The app opens a local browser interface for selecting datasets, choosing a
model, defining ROIs, running fits, reviewing plots, and writing local export
bundles.

### 4. Run The CLI

Example diffusion run with the tracked sample and original guide's Venus
cytoplasm ROI:

```bash
frap-toolbox sample-data/Diffusion/Venus_Cytoplasm_1.lsm \
  --roi 256 23 9 \
  --post-bleach-frame 21 \
  --pre-bleach-count 10 \
  --background 0 \
  --use-adjacent-roi
```

Batch the local guide diffusion stacks, when the ignored `test-data/` archive is
available:

```bash
frap-toolbox test-data/Diffusion/Venus_Cytoplasm_*.lsm \
  --roi 256 23 9 \
  --post-bleach-frame 21 \
  --pre-bleach-count 10 \
  --background 0 \
  --use-adjacent-roi
```

Useful CLI options:

- `--model`: choose `diffusion`, `reaction1`, or `reaction2`.
- `--fit-mode`: choose `global`, `individual`, `average_curve`,
  `simplified_kang`, or `simplified_kang_global`. Defaults are `global` for
  diffusion, `individual` for Reaction 1, and `average_curve` for Reaction 2.
- `--optimizer-mode`: defaults to `modern`; use `legacy_matlab` only for
  historical parity investigations.
- `--normalize-by-cell`: normalize by a whole-cell ROI or mask.
- `--bleach-mask` and `--cell-mask`: load saved ROI masks, usually `.npz`.
- `--output-dir`: write an export bundle with parameters, FRAP series, optional
  profiles, ROI masks, and run metadata.

## Running Analyses

Every workflow needs the same analysis inputs:

- Raw image files for one or more FRAP recovery sequences.
- Model choice: diffusion, Reaction 1, or Reaction 2.
- Bleach ROI geometry or a saved bleach mask.
- Post-bleach frame number.
- Background intensity, often measured from unlabeled controls.
- Pre-bleach frame count for normalization.
- Whole-cell ROI or mask if whole-cell normalization is enabled.

### Model Selection

| Model | Use it for | Required ROI style |
| --- | --- | --- |
| Diffusion | Single-component Brownian diffusion and mobile-fraction estimates. | Circular bleach ROI. |
| Reaction 1 | Single-exponential recovery from one dominant binding/reaction process. | Circular or user-defined bleach ROI. |
| Reaction 2 | Two-exponential recovery where two kinetic components are needed. | Circular or user-defined bleach ROI. |

The diffusion model reports bleach-depth/profile parameters, diffusion
coefficient `D`, mobile fraction `MF`, corrected mobile fraction when requested,
and residual sum of squares. Reaction models report fitted kinetic parameters
and residual sum of squares.

### Fit Modes

The modern diffusion path supports several fit modes:

- `global`: shared geometry and shared `D`/`MF` across pooled curves.
- `individual`: fit each stack and summarize per-stack results.
- `average_curve`: fit the averaged FRAP/profile curves.
- `simplified_kang`: per-curve Kang half-time estimator.
- `simplified_kang_global`: pooled fit of the simplified profile/recovery
  equations.

Use `global` for the current default Python workflow. Use legacy MATLAB parity
tests when the goal is to reproduce historical parameter tables exactly.

## Supported Image Formats

The legacy MATLAB application opens raw files through Bio-Formats. The original
guide verified Zeiss `.lsm` and Nikon `.nd2` files.

The modern Python package uses BioIO-backed readers plus optional reader extras.
Install the `nd2`, `bioformats`, or `legacy-io` extras when your microscope
format needs them. Reader provenance is recorded so parity and fallback behavior
can be audited.

Large microscopy fixtures are intentionally not tracked in git. Public checkouts
skip fixture-dependent tests automatically. Local parity work expects the
ignored `test-data/` archive described in [docs/developer-setup.md](../developer-setup.md).

## Data Availability

This repository tracks one raw microscopy sample for demos and smoke tests:
[`sample-data/Diffusion/Venus_Cytoplasm_1.lsm`](../../sample-data/Diffusion/Venus_Cytoplasm_1.lsm).

The full original user-guide fixture set is about 1.2 GB and remains outside
git. Download the published Zenodo archive
[`FRAP-Toolbox Legacy User Guide Test Data`](https://doi.org/10.5281/zenodo.20344310)
and unpack it at the repository root to restore `test-data/`. Verify the
download against the manifest and checksums in
[docs/data-availability.md](../data-availability.md).
For a scripted restore, run `python scripts/restore_test_data.py` from the
repository root.

## ROI Guidance

### Diffusion ROIs

For diffusion analysis:

1. Use a circular bleach ROI.
2. Record the center `(x, y)` and radius in pixels.
3. Record the first post-bleach frame.
4. Use an adjacent ROI when corrected mobile fraction is required.
5. Keep the full cell in frame if whole-cell normalization is planned.

The original guide's first diffusion example uses center `(256, 23)`, radius
`9`, post-bleach frame `21`, background `0`, and `10` pre-bleach images.

### Reaction ROIs

For Reaction 1 and Reaction 2 analysis:

1. Use a circular bleach ROI or a user-defined bleach polygon/mask.
2. Use a whole-cell ROI or mask when whole-cell normalization is enabled.
3. Record the post-bleach frame and pre-bleach frame count.
4. Save hand-drawn masks when you need reproducible reruns.

The modern app can author preview-backed ROIs and can load saved masks. The
canonical `.npz` mask container is documented in
[docs/roi-mask-format.md](../roi-mask-format.md).

## Results And Exports

The modern app and CLI can write analysis export bundles containing:

- `result-parameters.csv`
- `frap-series.csv`
- `post-bleach-profile.csv`, for diffusion workflows
- `roi-masks.npz`
- `run-metadata.json`

The legacy MATLAB application writes tab-delimited text files. Common legacy
suffixes are:

- `*_Diffusion_Fit_Parameters.txt`
- `*_Diffusion_FRAP_datasets.txt`
- `*_Diffusion_Postbleach_profiles.txt`
- `*_Reaction_Fit_Parameters.txt`
- `*_Reaction_FRAP_datasets.txt`
- `*_Reaction2_Fit_Parameters.txt`, as used by the archived Reaction 2 guide
  export in this repository

## Legacy MATLAB User Guide

This section preserves the installation and usage material from the original
supplemental PDF. The external download URLs are historical references and may
not represent current distribution channels.

### Legacy Installation Options

The original MATLAB-era guide described two ways to run FRAP-Toolbox:

1. Use the source files with a full installation of MATLAB.
2. Use the standalone application after installing the royalty-free MATLAB
   Compiler Runtime.

Standalone application and source-code information are now directed to the
GitHub repository: [https://github.com/kraftlj/FRAP-Toolbox](https://github.com/kraftlj/FRAP-Toolbox).

### Legacy System Requirements

FRAP-Toolbox was tested on:

- Windows XP, 32-bit.
- Windows 7, 64-bit.
- Mac OS X 10.9.

The original guide referenced MATLAB 2013 system requirements from MathWorks:
[http://www.mathworks.com/support/sysreq/current_release/](http://www.mathworks.com/support/sysreq/current_release/).

### Running The MATLAB Source

In MATLAB, navigate to `legacy/matlab/` and run the archived entry point.

From the repository root, you can also add the archived source directory to the
MATLAB path:

```matlab
addpath(fullfile(pwd, 'legacy', 'matlab'));
Main_GUI;
```

The current repository entry point for the archived MATLAB application is
`legacy/matlab/Main_GUI.m`.

### Legacy Standalone Application On Windows

1. Download the installation files from the historical download page.
2. Move the `FRAP-Toolbox` folder to a suitable location, such as
   `C:\FRAP-Toolbox`.
3. Install the MATLAB Compiler Runtime by double-clicking
   `MCR_R2013a_win32_installer.exe` and following the on-screen instructions.
4. Open `classpath.txt` for editing. By default, this file was located at
   `C:\Program Files\MATLAB\MATLAB Compiler Runtime\v81\toolbox\local\classpath.txt`.
5. If needed, grant edit permissions through **Properties > Security > Edit**.
6. Append this line to `classpath.txt`:

   ```text
   C:\FRAP-Toolbox\loci_tools.jar
   ```

7. Save `classpath.txt`.
8. Run FRAP-Toolbox by double-clicking `FrapToolbox.exe`.

The `classpath.txt` path must match the actual location of `loci_tools.jar`.

### Legacy Standalone Application On macOS

1. Download the installation files from the historical download page.
2. Right-click `FRAPToolbox_Installer` and choose **Open**.
3. The installer creates `/Applications/FRAPToolbox` and installs the MATLAB
   Compiler Runtime in `/Applications/MATLAB/`.
4. Move `loci_tools.jar` into `/Applications/FRAPToolbox/`.
5. Open the original `classpath.txt` location:

   ```text
   /Applications/MATLAB_R2013b/toolbox/local/classpath.txt
   ```

   You may need to right-click `MATLAB_R2013b` and choose **Show Package
   Contents**. Administrative privileges may be required.

6. Append this line:

   ```text
   /Applications/FRAPToolbox/loci_tools.jar
   ```

7. Save `classpath.txt`.
8. Run the `FrapToolbox` app from `/Applications/FRAPToolbox/application/`.

### Legacy MATLAB Analysis Flow

The original MATLAB main window asks for:

- The directory containing raw FRAP data.
- One or more selected files.
- The analysis model.
- Bleach ROI geometry.
- Post-bleach frame number.
- Constant background intensity.
- Whether to normalize by whole-cell intensity.
- Number of pre-bleach images used for normalization.

The **Image Preview** button loads the first selected dataset, provides a stack
scrubber, and overlays the bleaching ROI. The **Next** button loads all selected
datasets through Bio-Formats and opens the data analysis window.

In the analysis window, users set initial guesses, bounds, and fit ranges; run
the optimizer; inspect fits and residuals; and save optimized parameters and
processed data as text files.

## Original Test Data Workflows

This section preserves the original guide workflows. The same parameters are
also used by modern parity tests when the ignored `test-data/` archive is
available.

### Diffusion Model Workflow

Using a Zeiss LSM 510 confocal microscope, a small 1 um circular region in the
cytoplasm of HeLa cells expressing Venus or Venus-ATG5 was photobleached. The
recoveries are fit well by the FRAP-Toolbox diffusion model.

1. For a one-file smoke test, use
   `sample-data/Diffusion/Venus_Cytoplasm_1.lsm`. For the full original dataset
   set, restore the external `test-data/` archive described in
   [docs/data-availability.md](../data-availability.md).
2. Open FRAP-Toolbox.
3. Select **Diffusion** from the model dropdown.
4. Select **Circle** for the ROI.
5. Set the post-bleach image to `21`.
6. Set the mean background intensity to `0`.
7. Do not normalize by the whole cell.
8. Use `10` pre-bleach images for normalization.
9. Select all test datasets named `Venus_Cytoplasm_*.lsm`.
10. Click **Next**.
11. Enter `256 23` for the bleach ROI center and `9` for the radius.
12. Choose **Yes** to calculate a corrected mobile fraction.
13. Configure curve-fitting parameters as shown in
    [Table S1](#table-s1-curve-fitting-parameters-for-the-diffusion-test-data).
14. Select **No** when asked whether to fit averaged data.
15. Make sure all datasets are selected and press **Run**.
16. Inspect the result windows for post-bleach profiles, FRAP curves, fits, and
    residuals.
17. Compare optimized parameters with
    [Table S2](#table-s2-optimized-curve-fitting-parameters-for-the-diffusion-test-data).
18. Save with a tag such as `Venus_Cytoplasm`. The MATLAB app appends
    `*_Diffusion_Fit_Parameters.txt`, `*_Diffusion_FRAP_datasets.txt`, and
    `*_Diffusion_Postbleach_profiles.txt`.

The same procedure applies to the Venus-ATG5 diffusion datasets.

### Reaction 1 Model Workflow

Using a Nikon Eclipse Ti confocal microscope, Venus was photobleached in the
nuclear region of HeLa cells using a user-defined bleach ROI. The recoveries are
fit well by the Reaction 1 model.

1. Restore the external `test-data/` archive described in
   [docs/data-availability.md](../data-availability.md).
2. Open FRAP-Toolbox.
3. Select **Reaction 1** from the model dropdown.
4. Select **User Defined** for the ROI.
5. Set the post-bleach image to `6`.
6. Set the mean background intensity to `0`.
7. Normalize by the whole cell.
8. Use `5` pre-bleach images for normalization.
9. Select all datasets named `Venus_*.nd2`.
10. Click **Next**.
11. Draw the bleaching ROI in the nuclear region, then draw the whole-cell ROI.
    The original guide shows examples in Figure S5A.
12. Configure fitting parameters as shown in
    [Table S3](#table-s3-curve-fitting-parameters-for-the-reaction-1-test-data).
13. Select **No** when asked whether to fit averaged data.
14. Make sure all datasets are selected and press **Run**.
15. Inspect FRAP curves, fits, averages, and residuals.
16. Compare optimized parameters with
    [Table S4](#table-s4-optimized-curve-fitting-parameters-for-the-reaction-1-test-data).
17. Save with a tag such as `Venus_NCTransport`. The MATLAB app appends
    `*_Reaction_Fit_Parameters.txt` and `*_Reaction_FRAP_datasets.txt`.

Because reaction ROIs are drawn manually, reproduced parameters may not be
identical to the guide values, but they should be close when equivalent ROIs are
used.

### Reaction 2 Model Workflow

Using a Nikon Eclipse Ti confocal microscope, Venus-ATG5 was photobleached in
the nuclear region of HeLa cells using a user-defined bleach ROI. The recoveries
are fit well by the Reaction 2 model.

1. Restore the external `test-data/` archive described in
   [docs/data-availability.md](../data-availability.md).
2. Open FRAP-Toolbox.
3. Select **Reaction 2** from the model dropdown.
4. Select **User Defined** for the ROI.
5. Set the post-bleach image to `6`.
6. Set the mean background intensity to `0`.
7. Normalize by the whole cell.
8. Use `5` pre-bleach images for normalization.
9. Select all datasets named `Venus-Atg5_*.nd2`.
10. Click **Next**.
11. Draw the bleaching ROI in the nuclear region, then draw the whole-cell ROI.
    The original guide shows examples in Figure S5B.
12. Configure fitting parameters as shown in
    [Table S5](#table-s5-curve-fitting-parameters-for-the-reaction-2-test-data).
13. Select **Yes** when asked whether to fit averaged data.
14. Make sure all datasets are selected and press **Run**.
15. Inspect FRAP curves, fits, averages, and residuals.
16. Compare optimized parameters with
    [Table S6](#table-s6-optimized-curve-fitting-parameters-for-the-reaction-2-test-data).
17. Save with a tag such as `Venus-Atg5_NCTransport`. The archived guide
    parameter export uses `*_Reaction2_Fit_Parameters.txt`.

## Troubleshooting

- If image files do not load, install the relevant reader extra, such as `nd2`
  or `bioformats`, and confirm the microscope file format is supported.
- If pixel-size warnings appear, the file metadata may not include usable
  physical pixel sizes. Fits can still run, but absolute distances may need
  manual verification.
- If a fit saturates bounds or covariance cannot be estimated, re-check ROI
  placement, post-bleach frame, background value, and fit ranges.
- For batch processing, all datasets should have matching acquisition structure.
- For reaction workflows, save and reuse ROI masks when exact reruns matter.
- For reproducible bugs, collect the error text, FRAP-Toolbox version, operating
  system, input settings, and a small non-working dataset if possible.

## Supplemental Tables

### Table S1. Curve Fitting Parameters For The Diffusion Test Data

| Parameter | Initial guess | Lower bound | Upper bound | Fixed/adjustable |
| --- | ---: | ---: | ---: | --- |
| `K` | 1 | 0 | Inf | Adjustable |
| `re` | 3 | 0 | Inf | Adjustable |
| `D` | 10 | 0 | Inf | Adjustable |
| `Mf` | 1 | 0 | 2 | Adjustable |
| `kdecay` | 1.00E-03 | 0 | Inf | Adjustable |

| Fit range | First data point | Last data point |
| --- | ---: | ---: |
| Profile Fit | 1 | 5626 |
| FRAP Fit | 1 | 600 |
| Decay Fit | 500 | 600 |
| Corrected MF | 540 | 600 |

### Table S2. Optimized Curve Fitting Parameters For The Diffusion Test Data

| FileNames | k | re | D | MF | MF Corrected | SS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Venus_Cytoplasm_1.lsm` | 1.49 | 3.52 | 10.71 | 0.779 | 0.944 | 1.50E-06 |
| `Venus_Cytoplasm_10.lsm` | 1.51 | 3.92 | 31.26 | 0.798 | 1.005 | 6.80E-07 |
| `Venus_Cytoplasm_11.lsm` | 1.25 | 4.25 | 69.16 | 0.825 | 0.991 | 1.56E-06 |
| `Venus_Cytoplasm_12.lsm` | 1.63 | 3.56 | 28.16 | 0.908 | 0.952 | 3.61E-07 |
| `Venus_Cytoplasm_13.lsm` | 1.42 | 3.63 | 37.93 | 0.842 | 0.963 | 1.97E-06 |
| `Venus_Cytoplasm_14.lsm` | 1.43 | 3.98 | 38.57 | 0.835 | 0.998 | 1.27E-06 |
| `Venus_Cytoplasm_15.lsm` | 1.37 | 3.86 | 32.21 | 0.795 | 0.944 | 1.70E-06 |
| `Venus_Cytoplasm_2.lsm` | 1.41 | 3.89 | 34.05 | 0.790 | 1.010 | 8.89E-07 |
| `Venus_Cytoplasm_3.lsm` | 1.24 | 4.33 | 34.78 | 0.795 | 1.001 | 1.73E-06 |
| `Venus_Cytoplasm_4.lsm` | 1.14 | 4.50 | 49.07 | 0.868 | 1.041 | 4.70E-07 |
| `Venus_Cytoplasm_5.lsm` | 1.34 | 3.91 | 35.12 | 0.809 | 1.008 | 2.00E-06 |
| `Venus_Cytoplasm_6.lsm` | 1.26 | 4.15 | 24.05 | 0.963 | 1.028 | 1.99E-06 |
| `Venus_Cytoplasm_7.lsm` | 1.22 | 4.06 | 37.39 | 0.857 | 1.055 | 1.80E-06 |
| `Venus_Cytoplasm_8.lsm` | 1.19 | 4.21 | 47.75 | 0.870 | 0.964 | 9.31E-07 |
| `Venus_Cytoplasm_9.lsm` | 1.41 | 3.99 | 24.77 | 0.874 | 0.993 | 6.77E-07 |
| `Avg.` | 1.35 | 3.98 | 35.67 | 0.841 | 0.993 | 1.40E-07 |
| `Venus-Atg5_Cytoplasm_1.lsm` | 2.13 | 3.03 | 12.36 | 0.761 | 0.960 | 2.87E-06 |
| `Venus-Atg5_Cytoplasm_10.lsm` | 2.18 | 2.99 | 8.78 | 0.902 | 1.019 | 3.03E-06 |
| `Venus-Atg5_Cytoplasm_11.lsm` | 1.92 | 3.26 | 13.38 | 0.866 | 1.049 | 3.25E-06 |
| `Venus-Atg5_Cytoplasm_12.lsm` | 2.08 | 2.96 | 10.59 | 0.848 | 1.002 | 3.66E-06 |
| `Venus-Atg5_Cytoplasm_13.lsm` | 2.09 | 3.07 | 10.11 | 0.791 | 0.985 | 2.99E-06 |
| `Venus-Atg5_Cytoplasm_14.lsm` | 2.60 | 2.84 | 5.37 | 0.760 | 0.995 | 8.72E-07 |
| `Venus-Atg5_Cytoplasm_15.lsm` | 1.92 | 3.15 | 17.84 | 0.800 | 0.978 | 3.13E-06 |
| `Venus-Atg5_Cytoplasm_2.lsm` | 2.30 | 3.04 | 8.90 | 0.805 | 0.983 | 2.08E-06 |
| `Venus-Atg5_Cytoplasm_3.lsm` | 2.03 | 3.13 | 13.77 | 0.824 | 1.009 | 3.04E-06 |
| `Venus-Atg5_Cytoplasm_4.lsm` | 2.20 | 2.90 | 13.39 | 0.809 | 1.016 | 3.59E-06 |
| `Venus-Atg5_Cytoplasm_5.lsm` | 2.30 | 2.83 | 7.97 | 0.770 | 0.936 | 2.46E-06 |
| `Venus-Atg5_Cytoplasm_6.lsm` | 2.30 | 2.93 | 7.18 | 0.844 | 0.972 | 1.25E-06 |
| `Venus-Atg5_Cytoplasm_7.lsm` | 2.09 | 3.10 | 10.46 | 0.854 | 1.049 | 2.11E-06 |
| `Venus-Atg5_Cytoplasm_8.lsm` | 2.20 | 2.91 | 9.52 | 0.821 | 0.991 | 2.71E-06 |
| `Venus-Atg5_Cytoplasm_9.lsm` | 2.01 | 3.13 | 6.03 | 0.859 | 1.006 | 1.96E-06 |
| `Avg.` | 2.16 | 3.02 | 10.38 | 0.821 | 0.997 | 2.50E-07 |

### Table S3. Curve Fitting Parameters For The Reaction 1 Test Data

| Parameter | Initial guess | Lower bound | Upper bound | Fixed/adjustable |
| --- | ---: | ---: | ---: | --- |
| `a` | 1 | 0 | Inf | Adjustable |
| `b` | 1 | 0 | Inf | Adjustable |
| `c` | 1 | 0 | Inf | Adjustable |
| `kdecay` | 1.00E-03 | 0 | Inf | Adjustable |

| Fit range | First data point | Last data point |
| --- | ---: | ---: |
| Profile Fit | 1 | 130 |
| Decay Fit | 135 | 185 |

### Table S4. Optimized Curve Fitting Parameters For The Reaction 1 Test Data

| FileNames | a | b | c | SS |
| --- | ---: | ---: | ---: | ---: |
| `Venus_1001.nd2` | 0.873656 | 0.762292 | 0.00844 | 8.33E-08 |
| `Venus_1002.nd2` | 0.855588 | 0.753077 | 0.015899 | 6.70E-08 |
| `Avg.` | 0.864622 | 0.757685 | 0.01217 | 2.10E-07 |

### Table S5. Curve Fitting Parameters For The Reaction 2 Test Data

| Parameter | Initial guess | Lower bound | Upper bound | Fixed/adjustable |
| --- | ---: | ---: | ---: | --- |
| `a` | 1 | 0 | Inf | Adjustable |
| `b` | 0.5 | 0 | Inf | Adjustable |
| `c` | 0.05 | 0 | Inf | Adjustable |
| `d` | 0.5 | 0 | Inf | Adjustable |
| `f` | 5.00E-04 | 0 | Inf | Adjustable |
| `kdecay` | 1.00E-03 | 0 | Inf | Adjustable |

| Fit range | First data point | Last data point |
| --- | ---: | ---: |
| Profile Fit | 1 | 185 |
| Decay Fit | 135 | 185 |

### Table S6. Optimized Curve Fitting Parameters For The Reaction 2 Test Data

| FileNames | a | b | c | d | f | SS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Venus-Atg5_1002.nd2` |  |  |  |  |  |  |
| `Venus-Atg5_1003.nd2` |  |  |  |  |  |  |
| `Avg.` | 1.470 | 0.2285 | 0.02078 | 1.129 | 0.0004425 | 1.78E-07 |

## Supplemental Figure Legends

The original PDF figures are not tracked in this repository. These legends are
preserved so the historical workflow remains understandable.

### Figure S1

FRAP-Toolbox begins with a main window requesting the data location, selected
files, model, bleach ROI geometry, first post-bleach image, background
fluorescence intensity, and normalization options. Batch processing requires all
datasets to have the same structure.

### Figure S2

Previewing a FRAP dataset allows the user to verify that basic inputs are
correct.

### Figure S3

The data analysis and visualization window allows the user to set initial
guesses, lower and upper bounds, and fit ranges; choose individual or averaged
fitting; exclude datasets; view optimized parameters; and save tab-delimited
output.

### Figure S4

Diffusion-model pop-up windows show radial post-bleach profiles, normalized FRAP
curves, fits, averages, and residuals.

### Figure S5

User-defined ROIs for Reaction 1 and Reaction 2 test data. Panel A contains
`Venus_1001.nd2` and `Venus_1002.nd2`. Panel B contains
`Venus-Atg5_1002.nd2` and `Venus-Atg5_1003.nd2`. Bleach ROIs are shown as
dashed black lines in the original figure, and whole-cell ROIs are shown as
dashed white lines.

## References

1. Linkert M, Rueden CT, Allan C, Burel J-M, Moore W, et al. (2010) Metadata
   matters: access to image data in the real world. *Journal of Cell Biology*
   189: 777-782.
2. Kraft L, Kenworthy A (2012) Imaging protein complex fluorescence recovery
   after photobleaching. *Journal of Biomedical Optics* 17: 11008.
3. Kraft LJ, Nguyen TA, Vogel SS, Kenworthy AK (2014) Size, stoichiometry, and
   organization of soluble LC3-associated complexes. *Autophagy* 10.
4. Day C, Kraft L, Kang M, Kenworthy A (2012) Analysis of protein and lipid
   dynamics using confocal fluorescence recovery after photobleaching (FRAP).
   *Current Protocols in Cytometry*, Chapter 2.
