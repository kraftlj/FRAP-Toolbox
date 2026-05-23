# MATLAB-to-Python Port Testing

## Overview

This note captures the investigation into why the Python diffusion fitter in
`frap_toolbox_py` diverges from the parameters reported by the original MATLAB
FRAP-Toolbox. The work focused on reproducing MATLAB’s fit for the dataset
`Venus_Cytoplasm_3.lsm` and then generalised to the entire diffusion folder.

For the broader guide-derived parity plan, including Reaction 1 and Reaction 2
coverage, see `docs/user-guide-parity-testing.md`.

## Key Findings

- **Model parity** – Using MATLAB’s own parameters (`k`, `re`, `D`, `MF`) inside
the Python forward model reproduces MATLAB’s exported FRAP fit and SSE to the
precision of the saved text files. The weighting function `frap/(t + Σfrap)` is
identical across MATLAB and Python.
- **Equation audit** – The guide-table parameters now reproduce every exported
  Venus Cytoplasm post-bleach profile fit, FRAP fit, weighted residual vector,
  and sum-squared residual to the precision of MATLAB's tab-delimited exports.
  This covers the Kang FRAP summation, profile equation, time slicing, nominal
  radius scaling, `f(1)` baseline term, and `(t + Σfrap)` weighting.
- **Photodecay audit** – MATLAB's `PhotoDecay.m` also used `lsqcurvefit`.
  Matching that stage required two details: the guide's one-based inclusive
  `500:600` range maps to Python `slice(499, 600)`, and MATLAB averaged FRAP
  curves by frame after a loose 10% time-vector check rather than interpolating
  time vectors. With the shared legacy trust-region path, Python regenerates the
  archived corrected FRAP curves within about `4e-6` for Venus Cytoplasm
  (`kdecay = 0.0051274`) and about `6e-6` for Venus-ATG5
  (`kdecay = 0.0055365`).
- **Objective landscape** – When scanning the sum-of-squares (SSE) over the
  diffusion coefficient, the global minimum lies near `D ≈ 74 µm²/s`,
  significantly lower in SSE than MATLAB's published `D = 34.78 µm²/s`.
  A bounded exhaustive profile, dense two-dimensional grid, differential
  evolution, and local minimization from the guide point all agree on the
  higher-D basin floor for `Venus_Cytoplasm_3.lsm`.
- **Legacy optimizer emulation** – The guide values are not a separate local
  minimum. They sit on the same basin shoulder, where the absolute SSE
  difference is small enough for MATLAB-era `lsqcurvefit` stopping behavior.
  The Python `legacy_matlab` mode now uses a shared trust-region-reflective
  emulation for photodecay and the D/MF fit: Coleman-Li scaling,
  finite-difference Jacobians, a MATLAB-style first-order optimality check, and
  an absolute residual-sum change stop. The path is calibrated against the
  archived guide vectors while keeping the modern lower-SSE optimizer separate.
- **Trimming & weighting** – MATLAB’s GUI optionally trims late-time points and
 weights early frames more heavily via the same `(t + Σfrap)⁻¹` factor. Matching
 trim ranges in Python nudges the fitted `D` slightly but never fully aligns with
 MATLAB’s numbers—the optimiser still slides toward the lower SSE. Therefore,
 trimming does not explain the discrepancy; the legacy optimizer path remains
 the active parity target.

## Reproduction Steps

1. Load MATLAB’s corrected FRAP data from
   `test-data/Diffusion/Venus_Cytoplasm_Diffusion_FRAP_datasets.txt`.
2. Fix `k` and `re` to the values in
   `Venus_Cytoplasm_Diffusion_Fit_Parameters.txt`.
3. Fit `D` and `MF` by minimising the weighted residual
   `(KangFRAP·MF + (1 − MF)·f0)/(t + Σfrap) − corrected/(t + Σfrap)`.
4. Profile the objective directly:
   - The guide point reproduces MATLAB's exported FRAP fit and SSE.
   - The bounded global D/MF minimum is at a higher D with lower SSE.
5. Vary solver tolerances, iteration limits, finite-difference steps, solver
   methods, and parameter scaling to characterize which early-stopping paths
   land closest to MATLAB's published table.

## Optimizer Path Notes

MathWorks documents `lsqcurvefit` as using `lsqnonlin`'s
trust-region-reflective algorithm for unconstrained or bound-constrained least
squares by default, with `TolFun`/`FunctionTolerance`, `TolX`/`StepTolerance`,
and first-order optimality defaults around `1e-6`. It also notes that the
default trust-region-reflective algorithm is a subspace method based on the
Coleman-Li interior-reflective Newton method and uses PCG internally:

- <https://www.mathworks.com/help/optim/ug/lsqcurvefit.html>
- <https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html>

The default SciPy trust-region path for `Venus_Cytoplasm_3.lsm`, with `k` and
`re` fixed to the MATLAB table values, moves through approximately:

| Step | D | MF | Weighted SSE |
| --- | ---: | ---: | ---: |
| start | 10.000 | 1.0000 | `3.36034e-05` |
| 1 | 14.781 | 0.8336 | `5.24058e-06` |
| 2 | 18.097 | 0.8126 | `3.91122e-06` |
| 3 | 24.730 | 0.8021 | `2.60923e-06` |
| 4 | 37.997 | 0.7918 | `1.59561e-06` |
| 5 | 63.380 | 0.7843 | `1.15978e-06` |
| 6 | 74.523 | 0.7870 | `1.11430e-06` |

The guide value for that dataset is `D=34.7806`, `MF=0.795314`, and
`SSE=1.73002e-06`, so it sits on the same descent shoulder rather than at a
separate local basin. Similar behavior is seen for the other largest D gaps,
especially `Venus_Cytoplasm_5.lsm` and `Venus_Cytoplasm_9.lsm`.

Changing SciPy tolerances or capping function evaluations can make the average
`D` resemble the guide average, but no single tested SciPy threshold/scale
setting reproduces the individual guide values reliably. This points to a
MATLAB R2013b `lsqcurvefit` trust-region path detail, such as internal scaling,
finite-difference handling, or the exact absolute stopping test, rather than a
remaining discrepancy in the FRAP equations or preprocessing.

The implemented legacy path follows the MATLAB-style early stopping behavior
more closely. Running the full MATLAB-export pipeline, including regenerated
photodecay correction from normalized FRAP, now gives:

| Dataset | MATLAB avg D | Legacy Python avg D | Avg D % diff | Max per-file D diff | MATLAB avg MF | Legacy Python avg MF |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Venus Cytoplasm | 35.6650 | 35.6442 | -0.058% | 0.4662 | 0.840538 | 0.840525 |
| Venus-ATG5 Cytoplasm | 10.3752 | 10.3781 | +0.028% | 0.0063 | 0.820940 | 0.820940 |

A remote MATLAB R2026a probe confirmed that current MATLAB default
`lsqcurvefit` reproduces the guide D/MF shoulder points for the largest Venus
Cytoplasm gaps when driven directly from the exported vectors. MATLAB stopped
with exit flag `3` after 2-3 iterations because the relative change in residual
sum of squares was below `FunctionTolerance = 1e-6`, even though the global
objective profile remains lower at higher D:

| Dataset | Guide D | R2026a D | Guide MF | R2026a MF | Exit flag | Iterations | R2026a SS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Venus_Cytoplasm_3.lsm` | 34.7806 | 34.7805966 | 0.795314 | 0.7953138 | 3 | 3 | `1.730022e-06` |
| `Venus_Cytoplasm_5.lsm` | 35.1196 | 35.1195862 | 0.808686 | 0.8086854 | 3 | 3 | `1.995470e-06` |
| `Venus_Cytoplasm_9.lsm` | 24.7671 | 24.7671242 | 0.873840 | 0.8738397 | 3 | 2 | `6.769451e-07` |

The same probe reported MATLAB R2026a Optimization Toolbox `26.1`, Image
Processing Toolbox `26.1`, default `Algorithm: trust-region-reflective`,
`TolFun: 1e-6`, `TolX: 1e-6`, `FinDiffType: forward`, and `TypicalX:
ones(numberofvariables,1)`.

## Practical Outcomes

- The Python implementation is faithful to MATLAB's model; the large Venus
  divergence was due to optimizer termination criteria, not FRAP equation
  differences.
- The default modern optimizer should remain available because it finds a
  lower-SSE solution. The `legacy_matlab` optimizer mode is for reproducing
  archived/public FRAP-Toolbox outputs and intentionally returns historical
  shoulder points.
- The modern app defaults to true global diffusion fitting. It fits shared `K`
  and `re` from pooled post-bleach profiles, then shared `D` and `MF` by
  minimizing the concatenated residual vectors from each input curve. The
  averaged curve is used only for diagnostics unless the user explicitly selects
  `average_curve`.
- A `simplified_kang` mode is also available for the Traffic 2012 confocal FRAP
  half-time estimator. This path is intentionally not a nonlinear curve fit: it
  estimates `tau_1/2` by linearly interpolating between corrected recovery
  samples, using the final included recovery value as `F_infinity`, and computes
  `D` from `(rn^2 + re^2) / (8 * tau_1/2)`. Its `re` is fit from the
  simplified confocal post-bleach profile `1 - K * exp(-2r^2 / re^2)`, not the
  legacy FRAP-Toolbox nested-exponential profile used by the Kang-series
  nonlinear fitter.
- A `simplified_kang_global` mode fits the simplified profile and recovery
  equations globally, using pooled profile residuals for shared `K`/`re` and
  pooled FRAP residuals for shared `D`/`MF`.
- To obtain the lower SSE solution inside MATLAB, tighten `lsqcurvefit`
  tolerances (`FunctionTolerance`, `StepTolerance`, `OptimalityTolerance`) or
  inspect the solver's iterative display/output function on the original MATLAB
  version.
- Reaction 1 and Reaction 2 backend fitters are now ported. Reaction 1 is
  pinned to the exported `Venus_NCTransport` guide vectors; Reaction 2 is pinned
  with synthetic equation/fitting tests and raw ND2 stack/timestamp smoke tests
  because this checkout includes its parameter table and raw ND2 files, but not
  its FRAP vector export.
- Reaction 1 raw diagnostics now show that the embedded ND2 stimulation polygon
  does not reproduce the exported `Venus_1002.nd2` Raw FRAP curve. Simple
  threshold-derived bleach/cell candidates are close but not exact, reinforcing
  that the missing artifact is the hand-drawn MATLAB analysis ROI, not the
  reaction model.

## Next Steps

- Keep `modern` as the default optimizer mode and reserve `legacy_matlab` for
  explicit CLI/app parity runs against historical MATLAB outputs.
- Document the tolerance sensitivity in user-facing docs so analysts understand
  why diffusion estimates can differ between MATLAB and Python runs.
- Continue improving strict raw-loader parity, especially circular ROI
  rasterization and normalized-vector parity from the original LSM files.
- Use `docs/roi-mask-format.md` to save or reconstruct user-defined bleach/cell
  ROI masks for Reaction 1 and Reaction 2 so the raw ND2 workflows can be
  tested end to end; the raw files are now available, but the manually drawn
  MATLAB analysis masks are not.
