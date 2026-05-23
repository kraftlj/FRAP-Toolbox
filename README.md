# FRAP-Toolbox

Modern Python tools for analyzing Fluorescence Recovery After Photobleaching
(FRAP) experiments.

FRAP-Toolbox now centers on the Python analysis engine, local browser app, and
command-line interface. The original MATLAB implementation is preserved under
`legacy/matlab/` for historical reproduction and parity checks.

## Quick Start

Create an environment and install the core package:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pytest -q
```

Install the local app stack when you want the browser UI:

```bash
python -m pip install -e ".[app-nd2,test]"
frap-toolbox-app
```

## Main Workflows

### Local App

Use `frap-toolbox-app` for an interactive local workflow. The Streamlit app
lets you choose diffusion, Reaction 1, or Reaction 2 models; define circular,
polygon, or saved-mask ROIs where supported; run fits; inspect diagnostic plots;
and write export bundles.

### CLI

Use `frap-toolbox` for repeatable analysis, batch processing, and regression
checks. Example diffusion run with the tracked sample:

```bash
frap-toolbox sample-data/Diffusion/Venus_Cytoplasm_1.lsm \
  --roi 256 23 9 \
  --post-bleach-frame 21 \
  --pre-bleach-count 10 \
  --background 0 \
  --use-adjacent-roi
```

Add `--output-dir` to write result parameters, FRAP series, ROI masks, and run
metadata.

## Optional Installs

- `python -m pip install -e ".[app]"` for the Streamlit app.
- `python -m pip install -e ".[nd2]"` for Nikon ND2 files through BioIO.
- `python -m pip install -e ".[bioformats]"` for Bio-Formats-backed readers.
- `python -m pip install -e ".[legacy-io]"` for the AICSImageIO fallback path.
- `python -m pip install -e ".[app-full,test]"` for the app plus all optional
  reader paths.

## Documentation

- Python setup, CLI usage, app workflow, exports, and troubleshooting:
  [`frap_toolbox_py/README.md`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/frap_toolbox_py/README.md)
- Main user guide:
  [`docs/user-guide/`](https://github.com/kraftlj/FRAP-Toolbox/tree/master/docs/user-guide)
- Local app deployment notes:
  [`docs/app-deployment.md`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/docs/app-deployment.md)
- Developer setup and parity testing:
  [`docs/developer-setup.md`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/docs/developer-setup.md)
- ROI mask format:
  [`docs/roi-mask-format.md`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/docs/roi-mask-format.md)

## Data Availability

The repository tracks one small raw microscopy sample for demos and smoke tests:
[`sample-data/Diffusion/Venus_Cytoplasm_1.lsm`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/sample-data/Diffusion/Venus_Cytoplasm_1.lsm).

The full legacy user-guide fixture archive is about 1.2 GB and remains outside
git. Download the Zenodo archive
[`FRAP-Toolbox Legacy User Guide Test Data`](https://doi.org/10.5281/zenodo.20344310)
and unpack it at the repository root to restore `test-data/`. Checksum and
layout details are in
[`docs/data-availability.md`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/docs/data-availability.md).
The helper `python scripts/restore_test_data.py` can download, verify, and
unpack the archive for local parity work.

## Legacy MATLAB Archive

The original MATLAB application is retained under
[`legacy/matlab/`](https://github.com/kraftlj/FRAP-Toolbox/tree/master/legacy/matlab)
as a historical artifact. Its archived entry point is
[`legacy/matlab/Main_GUI.m`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/legacy/matlab/Main_GUI.m).
New analysis work should
start with the Python app or CLI; use MATLAB only for reproducing old workflows,
checking parity, or investigating historical exports.

The previous website-style README content is preserved at
[`legacy/matlab/legacy-website-content.md`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/legacy/matlab/legacy-website-content.md).

## Citing FRAP-Toolbox

```text
Kraft LJ, Dowler J, Day CA, Kang M, Kenworthy AK. (2014). FRAP-Toolbox:
Software for the analysis of Fluorescence Recovery After Photobleaching.
https://github.com/kraftlj/FRAP-Toolbox (accessed Month Day, Year).
```

FRAP-Toolbox is released under the GNU General Public License.
Machine-readable citation metadata is available in
[`CITATION.cff`](https://github.com/kraftlj/FRAP-Toolbox/blob/master/CITATION.cff).
