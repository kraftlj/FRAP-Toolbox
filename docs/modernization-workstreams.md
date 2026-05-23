# Modernization Status

This file records the current mainline status after the Python-first repository
refocus.

## Mainline Direction

- `master` centers on the Python analysis package, CLI, and local Streamlit app.
- The original MATLAB implementation is archived under `legacy/matlab/` for
  historical reproduction and parity investigations.
- The cloud web pilot is preserved outside mainline on `codex/cloud-web-pilot`.

## Current Technical Facts

- Diffusion, Reaction 1, and Reaction 2 fitting are available through the Python
  engine, CLI, and local app.
- Saved ROI masks use the canonical `.npz` format documented in
  `docs/roi-mask-format.md` and can be consumed by CLI/app reaction workflows.
- Optional full fixture restoration is manual through
  `scripts/restore_test_data.py`; installs, imports, tests, and app startup do
  not download microscopy data.
- Reaction 1 model parity is pinned to MATLAB-exported guide vectors.
- Reaction 2 has model support and raw ND2 smoke coverage, but the fixture
  archive still lacks a Reaction 2 FRAP vector export.
- Strict diffusion live-loader parity has one expected `xfail` while circular
  ROI rasterization and normalized-vector parity are matched to the archived
  MATLAB exports.
- Exact raw reaction curve parity still needs the hand-drawn MATLAB bleach and
  whole-cell ROI masks collected through `roipoly`.

## Acceptance Gates

Routine mainline checks should include:

- `git diff --check` or the staged equivalent before committing.
- `.venv/bin/python -m pytest -q`, with only documented skips or xfails.
- `.venv/bin/python -m pip check`.
- Package build and metadata validation with `build` and `twine`.
- Markdown link validation with `python scripts/check_markdown_links.py`.
- Wheel-content smoke checks so tests, generated exports, sample microscopy
  data, and legacy MATLAB artifacts are not shipped in the Python wheel.

## Near-Term Opportunities

1. Close strict diffusion loader parity and remove the current xfail.
2. Capture or reconstruct Reaction 1 and Reaction 2 MATLAB analysis ROI masks.
3. Add a Reaction 2 FRAP vector export to the external fixture archive.
4. Keep the Python package distribution focused on installable app/library
   code, while repository docs retain the full historical and parity context.
