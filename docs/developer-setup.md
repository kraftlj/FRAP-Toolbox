# Developer Setup, CI, and Packaging

This project is in transition from the original MATLAB application to a modern
Python package. Keep public CI lightweight and reproducible, and keep large
microscopy fixtures local unless the project later adopts an explicit fixture
distribution mechanism.

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
- tests that require them must call `pytest.skip` when files are absent;
- unit and synthetic regression tests should still run normally;
- packaging checks build the source distribution and wheel, then run
  `twine check`.

The public gate is:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
python -m pip check
```

## Local MATLAB Parity Tests

If you have the ignored `test-data/` archive beside the repository root, pytest
will automatically run the guide-derived parity tests. The local archive should
contain at least:

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
under `scratch/matlab-parity-output/` without modifying source files.

## Packaging Checks

Run these before publishing a release or opening a packaging-focused PR:

```bash
python -m pip check
python -m pip install build twine
rm -rf dist
python -m build --sdist --wheel
python -m twine check dist/*
```

Do not commit `dist/`, `.venv/`, `test-data/`, generated exports, or scratch
debugging files.
