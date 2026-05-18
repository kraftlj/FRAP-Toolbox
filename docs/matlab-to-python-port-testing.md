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
the Python forward model reproduces MATLAB’s corrected FRAP curve and SSE to
machine precision. The weighting function `frap/(t + Σfrap)` is identical across
MATLAB and Python.
- **Objective landscape** – When scanning the sum-of-squares (SSE) over the
  diffusion coefficient, the global minimum lies near `D ≈ 74 µm²/s`,
  significantly lower in SSE than MATLAB's published `D = 34.78 µm²/s`.
  A bounded exhaustive profile, dense two-dimensional grid, differential
  evolution, and local minimization from the guide point all agree on the
  higher-D basin floor for `Venus_Cytoplasm_3.lsm`.
- **Impact of solver tolerances** – The guide values are not a separate local
  minimum. They sit on the same basin shoulder, where the absolute SSE
  difference is small enough to be plausible for MATLAB-era `lsqcurvefit`
  termination behavior. Matching historical output therefore requires emulating
  the old solver's stopping path, not merely running a stronger optimizer.
- **Trimming & weighting** – MATLAB’s GUI optionally trims late-time points and
 weights early frames more heavily via the same `(t + Σfrap)⁻¹` factor. Matching
 trim ranges in Python nudges the fitted `D` slightly but never fully aligns with
 MATLAB’s numbers—the optimiser still slides toward the lower SSE. Therefore,
 trimming does not explain the discrepancy; stopping tolerances do.

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
5. Vary solver tolerances and iteration limits to characterize which
   early-stopping paths land closest to MATLAB's published table.

## Practical Outcomes

- The Python implementation is faithful to MATLAB's model; the divergence is
  due to optimizer termination criteria, not mathematical differences.
- To reproduce MATLAB numbers in Python, the legacy mode needs a deliberate
  MATLAB-era stopping-path emulation. To obtain the lower SSE solution inside
  MATLAB, tighten `lsqcurvefit` tolerances (`FunctionTolerance`,
  `StepTolerance`, `OptimalityTolerance`) or inspect the solver's iterative
  display/output function on the original MATLAB version.

## Next Steps

- Decide which set of tolerances should be the default for the Python CLI/GUI.
- Document the tolerance sensitivity in user-facing docs so analysts understand
  why diffusion estimates can differ between MATLAB and Python runs.
- If alignment with historical MATLAB outputs is essential, keep a switch that
  emulates the MATLAB stopping criteria and document that it intentionally
  returns the historical shoulder point rather than the lower-SSE modern
  optimum.
