# Developer Setup, CI, and Packaging

This project centers on the modern Python package, CLI, and local Streamlit app.
The original MATLAB application is archived under `legacy/matlab/` for parity
and historical reproduction. Keep public CI lightweight and reproducible. The
repository tracks one raw microscopy sample under `sample-data/`; the full
legacy fixture archive remains external and is documented in
`docs/data-availability.md`.

## Supported Python Versions

The Python package declares `requires-python >=3.10`. CI tests Python 3.10,
3.12, and 3.13. Local work should use one of those versions even if an older
environment happens to import the package.

## Install Paths

Core analysis library:

```bash
python -m pip install -e .
```

Unit tests:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

Contributor tooling:

```bash
python -m pip install -e ".[dev]"
```

Local Streamlit app:

```bash
python -m pip install -e ".[app,test]"
frap-toolbox-app
```

Optional microscope readers:

```bash
python -m pip install -e ".[nd2]"
python -m pip install -e ".[bioformats]"
python -m pip install -e ".[legacy-io]"
```

Install all parity-oriented reader paths for local validation:

```bash
python -m pip install -e ".[test,nd2,bioformats,legacy-io]"
```

The `bioformats` extra requires a Java runtime. Use it only when BioIO's native
format plugins are insufficient.

## Clean Checkout CI

The GitHub Actions workflow in `.github/workflows/python-tests.yml` runs on a
clean checkout where `test-data/` is absent. This is intentional:

- microscopy fixtures are large and ignored by git;
- the tracked `sample-data/` file is available for smoke coverage;
- tests that require them must call `pytest.skip` when files are absent;
- unit and synthetic regression tests should still run normally;
- packaging checks build the source distribution and wheel, then run
  `twine check`.

The public gate is:

```bash
git diff --check
python -m pip install -e ".[dev]"
python -m ruff check .
python scripts/check_markdown_links.py
python -m pytest -q
python -m pip check
```

## Local MATLAB Parity Tests

If you have restored the ignored `test-data/` archive beside the repository
root, pytest will automatically run the guide-derived parity tests. Download the
external archive from the Zenodo record
[FRAP-Toolbox Legacy User Guide Test Data](https://doi.org/10.5281/zenodo.20344310),
then verify it against `docs/data-availability.md`. You can restore and verify
the archive with:

```bash
python scripts/restore_test_data.py
```

It should contain at least:

- `test-data/Userguide.pdf`
- `test-data/Diffusion/*_Diffusion_Fit_Parameters.txt`
- `test-data/Diffusion/*_Diffusion_FRAP_datasets.txt`
- `test-data/Diffusion/*_Diffusion_Postbleach_profiles.txt`
- `test-data/Diffusion/*.lsm`
- `test-data/Reaction 1/*_Reaction_Fit_Parameters.txt`
- `test-data/Reaction 1/*_Reaction_FRAP_datasets.txt`
- `test-data/Reaction 1/*.nd2`
- `test-data/Reaction 2/*_Reaction2_Fit_Parameters.txt`
- `test-data/Reaction 2/*.nd2`

Run the parity layer directly with:

```bash
python -m pytest -q frap_toolbox_py/tests/test_user_guide_parity.py
```

Known current caveats:

- strict live diffusion curve parity has one intentional `xfail` while circular
  ROI rasterization is being matched to the MATLAB-era mask;
- Reaction 1 refits match the MATLAB-exported FRAP vectors, but raw ND2 curve
  parity still needs saved bleach and whole-cell ROI masks;
- Reaction 2 has a parameter-table fixture and raw ND2 files, but still needs a
  FRAP vector export or stored ROI masks for guide-level curve parity.

If a remote collaborator has MATLAB installed, send them
`docs/remote-matlab-parity-probe.md`. The tracked
`scripts/matlab_parity_probe.m` script writes a self-contained evidence bundle
under `scratch/matlab-parity-output/` without modifying source files. The script
adds `legacy/matlab/` to the MATLAB path automatically.

## Packaging Checks

Run these before publishing a release or opening a packaging-focused PR:

```bash
python -m pip check
python -m pip install -e ".[dev]"
rm -rf dist
python -m build --sdist --wheel
python -m twine check dist/*
```

Do not commit `dist/`, `.venv/`, `test-data/`, generated exports, or scratch
debugging files.

## Release Checklist

Before tagging a release:

1. Update `CHANGELOG.md` and confirm `CITATION.cff` still matches the release.
2. Run the public gate, package checks, and wheel install smoke test.
3. Confirm `run-metadata.json` remains backward compatible apart from additive
   fields such as `analysis_bundle_version`.
4. Confirm the known parity caveats are current in `docs/modernization-workstreams.md`.
5. Publish artifacts through GitHub releases or package infrastructure, not by
   committing generated archives.
