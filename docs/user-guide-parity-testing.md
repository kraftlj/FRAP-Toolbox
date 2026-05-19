# User Guide Parity Testing

This test design treats `test-data/Userguide.pdf` and the exported files beside
it as the compatibility contract for the Python modernization. The goal is to
separate three questions that are easy to conflate:

- Can the modern app open the same images with BioIO-backed readers?
- Does the preprocessing reproduce the MATLAB-exported curves, profiles, and
  tables?
- When model outputs differ, is the difference caused by data handling, model
  math, or optimizer stopping behavior?

The tests should be optional in public source checkouts because `test-data/`
contains large local fixtures and is intentionally ignored by git.

## Guide-Derived Workflows

### Diffusion

The guide's diffusion example uses:

- Model: Diffusion
- Files: `test-data/Diffusion/Venus_Cytoplasm_*.lsm` and
  `test-data/Diffusion/Venus-Atg5_Cytoplasm_*.lsm`
- ROI mode: Circle
- Bleach ROI for `Venus_Cytoplasm_1.lsm`: center `(256, 23)`, radius `9`
- Post-bleach frame: `21` in the MATLAB UI, represented as zero-based index
  `20` in Python
- Background: `0`
- Whole-cell normalization: off
- Pre-bleach frames: `10`
- Corrected mobile fraction: on, using the adjacent ROI
- Fit setup from Table S1:
  - `k`: initial `1`, lower `0`, upper `Inf`, adjustable
  - `re`: initial `3`, lower `0`, upper `Inf`, adjustable
  - `D`: initial `10`, lower `0`, upper `Inf`, adjustable
  - `MF`: initial `1`, lower `0`, upper `2`, adjustable
  - `kdecay`: initial `1e-3`, lower `0`, upper `Inf`, adjustable
  - Profile fit range `1:5626`, FRAP fit range `1:600`,
    decay fit range `500:600`, corrected mobile fraction range `540:600`
- Expected exports:
  - `*_Diffusion_Fit_Parameters.txt`
  - `*_Diffusion_FRAP_datasets.txt`
  - `*_Diffusion_Postbleach_profiles.txt`

### Reaction 1

The guide's first reaction example uses:

- Model: Reaction 1
- Files: `test-data/Reaction 1/Venus_*.nd2`
- ROI mode: User Defined
- Post-bleach frame: `6` in the MATLAB UI, zero-based index `5` in Python
- Background: `0`
- Whole-cell normalization: on
- Pre-bleach frames: `5`
- Fit setup from Table S3:
  - `a`, `b`, `c`: initial `1`, lower `0`, upper `Inf`, adjustable
  - `kdecay`: initial `1e-3`, lower `0`, upper `Inf`, adjustable
  - Profile fit range `1:130`, decay fit range `135:185`
- Expected exports:
  - `Venus_NCTransport_Reaction_Fit_Parameters.txt`
  - `Venus_NCTransport_Reaction_FRAP_datasets.txt`

### Reaction 2

The guide's second reaction example uses:

- Model: Reaction 2
- Files: `test-data/Reaction 2/Venus-Atg5_*.nd2`
- ROI mode: User Defined
- Post-bleach frame: `6` in the MATLAB UI, zero-based index `5` in Python
- Background: `0`
- Whole-cell normalization: on
- Pre-bleach frames: `5`
- Fit averaged data: yes
- Fit setup from Table S5:
  - `a`: initial `1`, lower `0`, upper `Inf`, adjustable
  - `b`: initial `0.5`, lower `0`, upper `Inf`, adjustable
  - `c`: initial `0.05`, lower `0`, upper `Inf`, adjustable
  - `d`: initial `0.5`, lower `0`, upper `Inf`, adjustable
  - `f`: initial `5e-4`, lower `0`, upper `Inf`, adjustable
  - `kdecay`: initial `1e-3`, lower `0`, upper `Inf`, adjustable
  - Profile fit range `1:185`, decay fit range `135:185`
- Expected export:
  - `Venus-Atg5_NCTransport_Reaction2_Fit_Parameters.txt`

The reaction guide depends on manually drawn bleach and whole-cell ROIs shown in
Figure S5. The Reaction 1 exported vectors are now enough to exercise the
Python reaction fitter directly. The raw ND2 files are present for Reaction 1
and Reaction 2, and tests now verify the readable ND2 timestamp metadata against
the MATLAB-exported Reaction 1 time vector for `Venus_1002.nd2`. Raw-image
reaction parity now has a durable `.npz` mask format for stored bleach and
whole-cell ROI masks; Reaction 2 guide refitting still needs either its FRAP
vector export or saved ROIs for the ND2 files. `Venus_1001.nd2` currently raises a legacy ND2
metadata parser error in both BioIO and AICSImageIO paths, so Reaction 1 raw
loader tests use the readable `Venus_1002.nd2` fixture.

## Test Layers

### 1. Reference Export Contract

These tests parse the MATLAB-exported text files and assert that key guide values
are available and unchanged. They do not exercise Python analysis code; they
protect the historical ground truth.

Acceptance examples:

- `Venus_Cytoplasm_1.lsm` has diffusion parameters `k=1.49227`,
  `re=3.52497`, `D=10.71`, `MF=0.778875`, corrected `MF=0.944335`.
- Diffusion average rows match the user guide for both Venus and Venus-Atg5.
- Reaction 1 rows match `Venus_1001.nd2`, `Venus_1002.nd2`, and `Avg.`.
- Reaction 2 has blank per-file rows and an averaged fit row matching Table S6.

### 2. Reader And Metadata Parity

These tests load the original files through the modern image loader and check
metadata-sensitive outputs that must match MATLAB before any modeling work can
be trusted.

Acceptance examples:

- LSM timestamps reproduce the exported `Time` vector after shifting to the
  post-bleach frame.
- Physical pixel size for Zeiss LSM diffusion data comes from LSM metadata.
- The Python loader records reader provenance so BioIO and fallback behavior can
  be audited.

### 3. Curve And Profile Parity

These tests compare Python preprocessing output against the MATLAB export
vectors:

- Raw FRAP
- Normalized FRAP
- Corrected FRAP
- Distance
- Post-bleach profile
- Corrected mobile fraction

This is the most important gate for "same operation as MATLAB." It should be
tight once parity is complete: timestamps within `1e-6`, curve/profile vectors
within export rounding, and corrected mobile fraction within `1e-4` or better.

Current modernization status: timestamps and post-bleach profile vectors are
now reproduced to export-rounding precision for `Venus_Cytoplasm_1.lsm`.
Photodecay correction is also reproduced from the MATLAB-exported normalized
FRAP curves across the Venus Cytoplasm guide set. The strict live-loader curve
test remains marked `xfail` because circular ROI rasterization still differs
from the 2014 MATLAB `poly2mask` behavior, which affects raw and normalized
FRAP curves.

For `Venus_Cytoplasm_1.lsm`, the MATLAB-exported Raw FRAP vector implies a
252-pixel bleach ROI. When that inferred mask is supplied to the Python loader,
the live image data reproduces the MATLAB Time, Raw FRAP, Normalized FRAP,
Distance, and Post-bleach Profile vectors to export precision. That proves the
remaining live-loader gap is localized to rasterizing the user-entered circle,
not LSM reading, timestamp extraction, profile construction, or normalization.
The inferred guide mask covers these 1-based row and inclusive column intervals:

| Row | Columns |
| --- | --- |
| 14 | 253-259 |
| 15 | 252-260 |
| 16 | 250-262 |
| 17 | 249-263 |
| 18 | 249-263 |
| 19 | 248-264 |
| 20 | 248-264 |
| 21 | 248-264 |
| 22 | 248-265 |
| 23 | 248-265 |
| 24 | 248-264 |
| 25 | 248-264 |
| 26 | 248-264 |
| 27 | 249-263 |
| 28 | 250-262 |
| 29 | 250-262 |
| 30 | 252-260 |
| 31 | 254-258 |

See `docs/roi-rasterization-survey.md` for the follow-up comparison against
OpenCV, scikit-image, Pillow, and matplotlib rasterizers.

See `docs/reporting-metric-sensitivity.md` for the downstream comparison of
ROI choices against publication-facing `D` and corrected mobile fraction values.

### 4. Model Math Parity

Once preprocessing vectors match, model tests should use the MATLAB-exported
curves and fit parameters directly:

- Evaluate the Python forward model using MATLAB parameters.
- Recompute FRAP fit, post-bleach profile fit, residuals, and SSE.
- Assert vector/SSE parity with the MATLAB export.

This isolates model equation fidelity from reader and optimizer differences.
The automated parity suite now includes a guide-level check that drives the
legacy diffusion fitter with MATLAB-exported vectors and reproduces the
`Venus_Cytoplasm_1.lsm` fit parameters from Table S2. The suite now also drives
both diffusion export tables through the legacy photodecay and optimizer path,
matching the published D/MF tables to guide-level tolerance while preserving a
direct landscape check showing that several guide D values are historical
early-stopping points, not the weighted least-squares basin floor.

For the modern app path, the diffusion fitter now distinguishes five fit modes:
`individual`, `average_curve`, `global`, `simplified_kang`, and
`simplified_kang_global`. The `global` mode is true global fitting: it estimates
one shared bleach/profile geometry from pooled post-bleach profiles, then
minimizes one concatenated residual vector across all FRAP curves with shared
`D` and `MF`. The averaged FRAP curve is still evaluated for plots and
summaries, but it is not the optimization target unless `average_curve` is
selected explicitly. The `simplified_kang` mode implements the Kang et al.
confocal half-time estimator `D = (rn^2 + re^2) / (8 * tau_1/2)` using each
corrected FRAP curve's linearly interpolated half-recovery time and an
effective bleach radius fit with the simplified confocal post-bleach profile
equation `1 - K * exp(-2r^2 / re^2)`. The `simplified_kang_global` mode uses
the same simplified profile/recovery equations, but fits shared `K`, `re`, `D`,
and `MF` against pooled profile and FRAP residuals.

### 5. Historical Optimizer Mode

The diffusion optimizer can produce different `D` values when SciPy continues
past MATLAB's historical stopping point. The Python port exposes a
`legacy_matlab` optimizer mode for the diffusion fitter. Tests for parameter
tables should remain split into two modes:

- `legacy_matlab`: emulate MATLAB-era `lsqcurvefit` stopping behavior and
  converge toward the guide's published parameter tables. The current hard
  passing checks cover `Venus_Cytoplasm_1.lsm` and both Venus/Venus-ATG5
  diffusion tables through the full exported-vector pipeline.
- `modern`: use stricter SciPy defaults and assert lower residuals or stable
  regression values documented for the new app. The current modern global-fit
  regression checks both Venus and Venus-ATG5 guide exports and asserts that the
  reported SSE is the sum of per-curve residual sums, not an averaged-curve SSE.

The user-facing app should make this distinction explicit if both modes are
offered.

### 6. Reaction Model Parity

The Python port now includes Reaction 1 and Reaction 2 backend fitters. Reaction
1 reproduces the guide's exported `Venus_NCTransport` vectors through the
MATLAB-compatible optimizer path:

- Single-component model: `F(t) = a - b * exp(-c * t)`.
- Optimization residual: `(F(t) - FRAP) / (t + sum(FRAP))`.
- Whole-cell normalization skips photodecay correction, matching
  `Figure_GUI_Reaction.m`.
- Averaged output follows the MATLAB script's index-wise curve average and
  first dataset time vector; it does not interpolate timestamp differences.

Reaction 2 backend support is also present:

- Two-component model: `F(t) = a - b * exp(-c * t) - d * exp(-f * t)`.
- Optimization residual follows `ReactionModel2.m`, which is unweighted even
  though the exported diagnostic residuals and SS are weighted.
- The current archive lacks a Reaction 2 FRAP vector export, so tests cover the
  equation/fitting behavior synthetically, verify raw ND2 timestamp loading for
  both Reaction 2 files, and keep the guide parameter-table contract intact.

### 7. Output Writer Parity

When exporters are implemented in Python, tests should compare the public file
contract:

- Expected file names for each workflow
- Row labels and column labels
- Blank cells in averaged-only Reaction 2 output
- Numeric formatting close enough for downstream spreadsheet use

Exact whitespace should not be the primary contract unless users depend on it;
semantic table equality is a better modernization target.

## Automation Rules

- Tests that require `test-data/` must call `pytest.skip` when fixtures are
  absent.
- Tests should use the original fixture paths but must never commit copied image
  data or generated exports.
- Slow full-folder parity tests should be marked separately from one-file smoke
  parity tests.
- Strict parity tests may start as `xfail` while a subsystem is being ported,
  but each `xfail` must state the current gap it represents.
- Once a subsystem reaches parity, remove the `xfail` and make the test part of
  the normal gate.

## Recommended First Gates

1. Keep the export-contract tests passing for diffusion, Reaction 1, and
   Reaction 2.
2. Make the diffusion `Venus_Cytoplasm_1.lsm` time and metadata parity test part
   of the normal suite.
3. Convert the strict diffusion curve parity test from `xfail` to passing after
   circular ROI preprocessing is matched to MATLAB's exported FRAP curves.
4. Save or reconstruct reaction ROI masks using the `.npz` contract in
   `docs/roi-mask-format.md`, then add Reaction 1 and Reaction 2 loader parity
   tests against the exported `*_FRAP_datasets.txt` files. Reaction 2 still
   needs a FRAP vector export before curve-level guide parity can be asserted.
5. Add model math parity tests before tightening optimizer parameter parity.
