# Remote MATLAB Parity Probe

This guide is for a remote colleague who can run MATLAB against the public
FRAP-Toolbox repository plus the ignored `test-data/` fixture archive. Its goal
is to collect the remaining MATLAB-only evidence needed to close the legacy
MATLAB to Python fit-result parity issue.

Current Python baseline from this checkout:

```bash
.venv/bin/python -m pytest -q frap_toolbox_py/tests/test_user_guide_parity.py
```

Expected result: `21 passed, 1 xfailed`. The remaining expected failure is
strict diffusion live-loader parity, currently localized to MATLAB circular ROI
rasterization. Model equations, timestamps, profile extraction, photodecay, and
exported-vector fit behavior are already covered by passing tests.

## What To Return

Return the whole `scratch/matlab-parity-output/` folder after running the probe.
It should contain:

- MATLAB and toolbox version facts, plus `which -all` output for `lsqcurvefit`,
  `lsqnonlin`, `poly2mask`, `roipoly`, and `bfopen`.
- MATLAB `poly2mask` outputs for the diffusion guide circular bleach ROI and
  adjacent ROI, including pixel counts, 1-based coordinates, and row intervals.
- Diffusion D/MF optimizer traces for `Venus_Cytoplasm_3.lsm`,
  `Venus_Cytoplasm_5.lsm`, and `Venus_Cytoplasm_9.lsm`, using the exported
  guide vectors and default `optimset('lsqcurvefit')` options.
- A short generated README describing which probe sections succeeded or failed.
- If time allows, a fresh Reaction 2 `*_Reaction2_FRAP_datasets.txt` export or
  saved `roipoly` bleach/cell masks for Reaction 1 and Reaction 2.

## Setup

1. Clone the public repository:

   ```bash
   git clone https://github.com/kraftlj/FRAP-Toolbox.git
   cd FRAP-Toolbox
   ```

2. Download and place the large fixture folder at repo root as `test-data/`.
   The expected layout includes:

   ```text
   test-data/Userguide.pdf
   test-data/Diffusion/*.lsm
   test-data/Diffusion/*_Diffusion_Fit_Parameters.txt
   test-data/Diffusion/*_Diffusion_FRAP_datasets.txt
   test-data/Diffusion/*_Diffusion_Postbleach_profiles.txt
   test-data/Reaction 1/*.nd2
   test-data/Reaction 1/*_Reaction_Fit_Parameters.txt
   test-data/Reaction 1/*_Reaction_FRAP_datasets.txt
   test-data/Reaction 2/*.nd2
   test-data/Reaction 2/*_Reaction2_Fit_Parameters.txt
   ```

3. In MATLAB, set the current folder to the repo root.

4. Run the tracked probe template:

   ```matlab
   addpath(fullfile(pwd, 'scripts'));
   matlab_parity_probe;
   ```

   The script writes only under `scratch/matlab-parity-output/`.

## Probe Details

The probe recreates the exact diffusion ROI code path from
`ROIinitialization_Diffusion.m`:

```matlab
x0 = 256;
y0 = 23;
R0 = 9;
t = 0:pi/20:2*pi;
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
bleachroimask = poly2mask(xi, yi, 512, 512);
adjacentroimask = poly2mask(xi + R0*2.5, yi, 512, 512);
```

It also reconstructs the D/MF weighted diffusion fit from `DiffusionModel_2.m`
using MATLAB-exported vectors:

```matlab
wf = f ./ (t + sum(f));
frapfun = @(p,t) ...
    (KangFRAP(t,re,rn,p(1),k).*p(2) + (1-p(2))*f(1)) ./ (t + sum(f));
p = lsqcurvefit(frapfun, [10, 1], t, wf, [0, 0], [Inf, 2], options);
```

For each target file, the probe saves:

- the guide point and its weighted SSE;
- a normal default `lsqcurvefit` result;
- a traced default `lsqcurvefit` result with an `OutputFcn`;
- a `.mat` file containing all raw `optimValues` structs recorded by MATLAB;
- a CSV with common trace fields for quick inspection.

## Optional Reaction Evidence

If there is enough time, run the legacy GUI workflow from `Main_GUI.m` for the
Reaction 2 guide data:

- Model: Reaction 2
- Files: `test-data/Reaction 2/Venus-Atg5_*.nd2`
- ROI mode: User Defined
- Post-bleach image: `6`
- Background: `0`
- Whole-cell normalization: on
- Pre-bleach images: `5`
- Fit averaged data: yes
- Save with tag `Venus-Atg5_NCTransport`

Copy the generated `Venus-Atg5_NCTransport_Reaction2_FRAP_datasets.txt` into
`scratch/matlab-parity-output/`.

If you draw Reaction 1 or Reaction 2 ROIs interactively, also save the resulting
`bleachroimask` and `cellroimask` variables as `.mat` files with the source
filename and image shape. These hand-drawn masks are the missing artifact for
exact raw ND2 reaction parity.

## Codex Prompt For The Remote Helper

```text
You are helping diagnose MATLAB-to-Python parity for FRAP-Toolbox.

Repo: https://github.com/kraftlj/FRAP-Toolbox
After cloning, place the downloaded large fixture folder at repo root as `test-data/`.

Read these files first:
- docs/matlab-to-python-port-testing.md
- docs/user-guide-parity-testing.md
- docs/roi-rasterization-survey.md
- docs/guide-refit-comparison.md
- docs/remote-matlab-parity-probe.md

Run the tracked MATLAB probe:
1. Open MATLAB at the repo root.
2. Run:
   addpath(fullfile(pwd, 'scripts'));
   matlab_parity_probe;
3. Inspect `scratch/matlab-parity-output/README.txt`.
4. Return the entire `scratch/matlab-parity-output/` folder.

If the script fails, preserve the generated output folder and include the full
MATLAB error text. Do not edit source files unless you are explicitly adding a
temporary local probe under `scratch/`.
```

## Acceptance Criteria

The evidence bundle is enough to close the issue if it contains:

- MATLAB's exact `poly2mask` mask for the diffusion guide ROI, including pixel
  count and row intervals.
- At least one complete MATLAB `lsqcurvefit` trace showing where the default
  solver stops for a high-gap guide file such as `Venus_Cytoplasm_3.lsm`.
- Confirmation of MATLAB default optimizer options for the installed release.
- Reaction 2 FRAP vector export or saved reaction ROI masks, if available.

Report the exact MATLAB release because solver behavior may differ from the
original R2013-era application.
