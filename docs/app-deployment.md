# Local App Deployment

This note tracks the practical user path for the modern Python app.

For the new cloud-first pilot architecture, see `docs/cloud-first-web-app.md`.
The Streamlit app below remains the local workstation/developer surface while the
Cloud Run + GCS + Batch workflow matures toward end-user parity.

## Recommended Local Setup

Use this setup for lab workstations that need to open Nikon ND2 files:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[app-nd2,test]"
frap-toolbox-app
```

Use `python -m pip install -e ".[app-full,test]"` when a workstation also needs
Bio-Formats or the legacy AICSImageIO fallback.

## Current App Workflow

1. Choose the model: `diffusion`, `reaction1`, or `reaction2`.
2. Choose a data directory and select one or more image files.
3. Define the bleach ROI as a circular numeric ROI. Reaction workflows can also
   upload a saved bleach mask.
4. Enable whole-cell normalization when needed and define the cell ROI the same way.
5. Set post-bleach frame, pre-bleach frame count, background, and fit mode.
6. Run analysis and inspect the parameter table, fit curve, residuals, and diffusion
   post-bleach profile when applicable.

The Streamlit app intentionally follows the CLI defaults:

- `post-bleach-frame`: 21
- `pre-bleach-count`: 10
- `background`: 0
- `normalize-by-cell`: off
- `fit-mode`: `global` for diffusion, `individual` for Reaction 1, `average_curve` for
  Reaction 2

## ROI Workflow

The stable path available now is circular numeric ROI entry plus saved-mask
upload. Diffusion bleach ROIs currently stay on the circular numeric path because
the diffusion model still needs nominal radius geometry; reaction bleach and
whole-cell ROIs can use saved masks.

Saved masks may be:

- `.npz`: the canonical FRAP-Toolbox ROI mask container described in
  `docs/roi-mask-format.md`
- `.npy`: a single 2-D NumPy array
- generic `.npz`: an array named for the ROI role, an array named `mask`, or
  exactly one array
- `.csv`, `.tsv`, or whitespace-delimited `.txt`: numeric 2-D masks

Nonzero values are treated as included pixels. Mask dimensions must match the image
height and width.

## Hand-Drawn ROI Design

The legacy MATLAB reaction workflows used hand-drawn bleach and whole-cell ROIs. The
Python app should keep that interaction but make it reproducible:

1. Show a preview frame with time scrubbing.
2. Let users draw polygon, freehand, and circle ROIs for bleach and cell regions.
3. Rasterize each ROI into a binary mask.
4. Save the masks and small metadata beside the analysis output.
5. Reuse saved masks in future runs without redrawing.

The current app already accepts the standardized FRAP-Toolbox `.npz` mask
container, so the next GUI step is drawing/saving controls that call the shared
ROI mask persistence helpers.

## Deployment Notes

The app launcher is the console script:

```bash
frap-toolbox-app
```

It runs Streamlit against `frap_toolbox_py/app.py`. No cloud service is required; data
stays on the local workstation.
