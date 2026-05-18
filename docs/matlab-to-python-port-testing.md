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
	significantly lower in SSE than MATLAB’s published `D = 34.78 µm²/s`. Gradient
	estimates at MATLAB’s point stay around `−5.7×10⁻⁸`, which is beneath
	MATLAB’s default optimality tolerance (1e−6). SciPy’s tighter default tolerances
	keep iterating and reach the global minimum; MATLAB stops when the gradient
	drops below its default threshold.
- **Impact of solver tolerances** – With MATLAB’s defaults replicated in Python
 (`least_squares` with `gtol ≈ 1e−6`, limited iterations, or larger
 `diff_step`), Python also halts near `D ≈ 34` despite the lower SSE elsewhere.
 Conversely, tightening MATLAB’s `FunctionTolerance` / `OptimalityTolerance`
 allows it to find the `D ≈ 74` minimum.
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
4. Vary solver tolerances:
   - `gtol ≈ 1e−6` (MATLAB-style) → convergence near MATLAB’s published `D`.
   - Stricter tolerances → convergence near `D ≈ 74` with lower SSE.

## Practical Outcomes

- The Python implementation is faithful to MATLAB’s model; the divergence is due
  to optimiser termination criteria, not mathematical differences.
- To reproduce MATLAB numbers in Python, relax SciPy tolerances or cap the
  iteration budget. To obtain the lower SSE solution inside MATLAB, tighten
  `lsqcurvefit` tolerances (`FunctionTolerance`, `StepTolerance`,
  `OptimalityTolerance`) or scale parameters to keep gradients above the default
  stopping threshold.

## Next Steps

- Decide which set of tolerances should be the default for the Python CLI/GUI.
- Document the tolerance sensitivity in user-facing docs so analysts understand
  why diffusion estimates can differ between MATLAB and Python runs.
- If alignment with historical MATLAB outputs is essential, expose a switch that
  emulates the MATLAB stopping criteria.
