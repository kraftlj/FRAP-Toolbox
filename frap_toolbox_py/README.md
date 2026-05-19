# Running the FRAP-Toolbox Python Port

This guide explains how to set up a Python environment and run the modern
FRAP-Toolbox app and CLI implemented in `frap_toolbox_py`. The steps assume
macOS or Linux with bash,
but only the activation command differs on Windows.

## 1. Prerequisites
- Python 3.10 or newer
- A C/Fortran toolchain for SciPy (Xcode CLT on macOS, build-essential on Linux)
- Java is only needed if you install the optional Bio-Formats reader extra

## 2. Create and activate a virtual environment
```bash
python3.13 -m venv .venv
source .venv/bin/activate
```

If you prefer `conda`, create and activate an equivalent Python 3.10+ environment instead.

## 3. Install the package and dependencies
Install the project in editable mode so CLI and app changes are picked up
automatically. For headless development and CI, use the core test stack:

```bash
pip install --upgrade pip
pip install -e ".[test]"
```

Install the local browser app stack when you want to run the GUI:

```bash
pip install -e ".[app,test]"
```

Optional extras:
- `pip install -e ".[nd2]"` for Nikon ND2 files through BioIO.
- `pip install -e ".[bioformats]"` for Bio-Formats-backed readers.
- `pip install -e ".[legacy-io]"` for the older AICSImageIO fallback.
- `pip install -e ".[qt]"` for the experimental Qt GUI dependencies.
- `pip install -e ".[app-nd2,test]"` for the most common local app setup with
  ND2 files.
- `pip install -e ".[app-full,test]"` for Streamlit plus ND2, Bio-Formats, and
  legacy I/O fallback.

For local guide parity work with all ignored microscopy fixtures available:

```bash
pip install -e ".[test,nd2,bioformats,legacy-io]"
```

## 4. Run the local app
```bash
frap-toolbox-app
```

The app opens a local Streamlit browser interface for selecting datasets, choosing
`diffusion`, `reaction1`, or `reaction2`, entering ROI parameters, running the fit, and
reviewing diagnostic plots.

Current ROI options:
- Circular numeric bleach ROI, matching the CLI `--roi X Y RADIUS` input.
- Optional circular whole-cell ROI when whole-cell normalization is enabled.
- Saved mask upload for reaction bleach and whole-cell ROIs. The app accepts
  `.npy`, `.npz`, `.csv`, `.tsv`, and whitespace-delimited `.txt` masks. The
  canonical `.npz` container is documented in `docs/roi-mask-format.md`;
  generic single-array masks are accepted for exploratory use.

Hand-drawn ROI support is expected to land as a small workflow on top of saved masks:
draw bleach and cell ROIs in an image preview, rasterize them once, save the masks, and
reuse those masks for repeatable fitting. Until that branch is merged, users can export
binary masks from another image tool and load them through the saved-mask upload controls.

## 5. Recommended test dataset parameters
The original MATLAB user guide documents parameters for the Zeiss LSM diffusion
datasets used during port validation. If you have the local `test-data/` folder
available, use:
- ROI center: `(256, 23)`
- ROI radius: `9`
- Post-bleach frame: `21`
- Pre-bleach frames for normalization: `10`
- Background intensity: `0`
- Enable corrected mobile fraction by using an adjacent ROI (2.5× radius offset)

## 6. Run the CLI
All commands assume the project root as the working directory and the virtual environment
is active.

Single stack:
```bash
frap-toolbox test-data/Diffusion/Venus_Cytoplasm_1.lsm \
  --roi 256 23 9 \
  --post-bleach-frame 21 \
  --pre-bleach-count 10 \
  --background 0 \
  --use-adjacent-roi
```

Batch the local Zeiss diffusion stacks, if you have the ignored `test-data/`
fixtures in your checkout:
```bash
frap-toolbox test-data/Diffusion/Venus_Cytoplasm_*.lsm \
  --roi 256 23 9 \
  --post-bleach-frame 21 \
  --pre-bleach-count 10 \
  --background 0 \
  --use-adjacent-roi
```

Useful flags:
- `--adjacent-offset` (default `2.5`): controls spacing of the adjacent ROI used when
  computing corrected mobile fractions.
- `--fit-mode` (default `global`): choose `global` for a true shared `K`/`re`/`D`/`MF`
  fit across pooled curves, `individual` for per-stack fits averaged afterward,
  `average_curve` to fit only the averaged FRAP/profile curves, `simplified_kang`
  for the Kang half-time estimator, or `simplified_kang_global` for a pooled global
  fit of the simplified profile/recovery equations.
- `--normalize-by-cell`: supply if you draw or compute a whole-cell mask and want to
  normalize by it.
- `--model`: choose `diffusion`, `reaction1`, or `reaction2`. Reaction workflows can
  use `--cell-roi X Y RADIUS` when whole-cell normalization is enabled.
- `--bleach-mask` and `--cell-mask`: load saved `.npz` ROI masks instead of circular
  ROI parameters. See `docs/roi-mask-format.md` for the container contract.

Diffusion output includes bleach depth, effective radius, diffusion coefficient, and
mobile fraction. Reaction output includes the fitted reaction parameters and residual
SS. The command returns a non-zero exit code if dataset loading or fitting fails,
making it suitable for scripted validation.

## 7. Run unit tests
```bash
python -m pytest
```

This executes the regression tests under `frap_toolbox_py/tests`. Integration
tests that need the large local `test-data/` fixtures skip automatically when
those files are not present, which is the expected behavior in CI and public
checkouts.

When `test-data/` is present locally, the same command also runs guide-derived
MATLAB parity tests. To focus on that layer:

```bash
python -m pytest -q frap_toolbox_py/tests/test_user_guide_parity.py
```

Add new tests before contributing model changes or additional readers.

## 8. Package checks
Before publishing or opening a packaging PR, run:

```bash
python -m pip check
python -m pip install build twine
python -m build --sdist --wheel
python -m twine check dist/*
```

## 9. Troubleshooting tips
- If you see repeated `Could not parse tiff pixel size` warnings, the TIFF metadata did not
  include a usable physical pixel size. The fits still run; specify voxel sizes manually if
  you need absolute distance units.
- A failure such as `Covariance of the parameters could not be estimated` typically means
  the fit saturated bounds or the ROI was mis-specified. Re-check the ROI center, radius,
  and post-bleach frame.
- The default reader path is BioIO plus `bioio-tifffile`. Add the `nd2` or `bioformats`
  extras when your microscope format needs another BioIO plugin.

## 10. Next steps
- Explore `frap_toolbox_py/data/loading.py` and `frap_toolbox_py/models/diffusion.py` to
  understand how stacks are loaded and fitted.
- Compare CLI output to the legacy MATLAB results stored in `test-data/Diffusion/*` to
  verify numerical agreement.
- Use the provided commands as a template for processing your own FRAP datasets.
