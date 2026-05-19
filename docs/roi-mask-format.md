# ROI Mask Format

The Python modernization stores hand-drawn ROIs as compressed NumPy `.npz`
containers. The format is intentionally small and local so the GUI, CLI, and
parity tests can share the same masks without depending on MATLAB workspace
files.

## Container

Each mask file contains:

- `metadata_json`: JSON text with format metadata.
- `mask__<name>` arrays: one boolean array per named ROI.

The current metadata schema is:

```json
{
  "format": "frap-toolbox-roi-masks",
  "version": 1,
  "image_shape": [512, 512],
  "roi_kind": "manual",
  "source_file": "Venus_1002.nd2",
  "notes": "drawn from legacy user-guide panel",
  "masks": [
    {"name": "bleach", "shape": [512, 512], "pixel_count": 252},
    {"name": "cell", "shape": [512, 512], "pixel_count": 18429}
  ]
}
```

`roi_kind` is usually `manual`, but can also describe generated masks such as
`circular`, `polygon`, or `inferred_legacy`. Mask names are limited to letters,
numbers, underscores, and hyphens. All saved mask arrays must have dtype `bool`
and must match `image_shape`.

## Python Helpers

```python
from frap_toolbox_py.roi import save_roi_masks, load_roi_mask

save_roi_masks(
    "Venus_1002_rois.npz",
    {"bleach": bleach_mask, "cell": cell_mask},
    source_file="Venus_1002.nd2",
    notes="drawn from the Reaction 1 guide image",
)

bleach_mask = load_roi_mask("Venus_1002_rois.npz", "bleach", expected_shape=(512, 512))
cell_mask = load_roi_mask("Venus_1002_rois.npz", "cell", expected_shape=(512, 512))
```

## CLI Usage

The existing circular ROI path remains supported:

```bash
frap-toolbox image.lsm --roi 256 23 9 --post-bleach-frame 21
```

For hand-drawn or reconstructed masks, provide saved mask files instead:

```bash
frap-toolbox Venus_1002.nd2 \
  --model reaction1 \
  --bleach-mask Venus_1002_rois.npz \
  --cell-mask Venus_1002_rois.npz \
  --normalize-by-cell \
  --post-bleach-frame 6 \
  --pre-bleach-count 5
```

`--bleach-mask` loads the named `bleach` mask. `--cell-mask` loads the named
`cell` mask. If a file contains only one mask, the CLI accepts that single mask
for the requested role.

## Modernization Impact

The legacy MATLAB reaction workflows used hand-drawn `roipoly` bleach and
whole-cell ROIs. Those analysis masks were not embedded as reusable artifacts in
the raw ND2 files, so raw Reaction 1 and Reaction 2 parity cannot be exact until
the masks are saved or reconstructed.

This format provides the missing artifact boundary:

1. The GUI can draw and save `bleach` and `cell` masks once.
2. The CLI and model loaders can consume the same masks deterministically.
3. Raw ND2 parity tests can compare generated FRAP vectors against MATLAB guide
   exports without redoing interactive ROI selection.
