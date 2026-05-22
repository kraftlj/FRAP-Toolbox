# Sample Data

This directory contains the one raw microscopy file tracked directly in GitHub
for demos and smoke tests:

- `Diffusion/Venus_Cytoplasm_1.lsm`

The file is the canonical first diffusion example from the original
FRAP-Toolbox user guide. It is small enough for the repository at about 30 MB,
but the full legacy fixture archive is about 1.2 GB and remains external. The
complete archive is published on Zenodo at
<https://doi.org/10.5281/zenodo.20344310>. See
[`docs/data-availability.md`](../docs/data-availability.md) for the full archive
manifest and restore instructions.

Use this sample for quick CLI checks:

```bash
frap-toolbox sample-data/Diffusion/Venus_Cytoplasm_1.lsm \
  --roi 256 23 9 \
  --post-bleach-frame 21 \
  --pre-bleach-count 10 \
  --background 0 \
  --use-adjacent-roi
```
