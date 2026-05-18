# Reporting Metric Sensitivity

This note checks whether the modernization preserves the values that are most
often reported from diffusion FRAP experiments: diffusion coefficient `D` and
corrected mobile fraction `MF Corrected`.

The comparison below uses the Venus Cytoplasm diffusion guide set with the guide
settings:

- Bleach ROI center `(256, 23)`, radius `9`
- Post-bleach frame `21` in the MATLAB UI, zero-based index `20` in Python
- Pre-bleach frames `10`
- Background `0`
- Adjacent ROI correction on
- Decay fit range `500:600` in MATLAB indexing
- Corrected mobile fraction range `540:600` in MATLAB indexing
- Individual dataset fitting, then averaging reported parameters

## ROI Options Checked

| ROI option | `Venus_Cytoplasm_1` D | Guide D | Avg D | Guide Avg D | Avg corrected MF | Guide Avg corrected MF |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Current `CircularROI` / `skimage.draw.disk` | 10.6162 | 10.7100 | 42.6070 | 35.6650 | 0.991987 | 0.993123 |
| 41-point circle polygon, zero-based vertices | 10.6378 | 10.7100 | 42.9706 | 35.6650 | 0.992123 | 0.993123 |
| Tuned polygon-containment diagnostic | 10.6751 | 10.7100 | 42.8142 | 35.6650 | 0.992928 | 0.993123 |

The corrected mobile fraction is robust to the practical ROI choice: the
current implementation's full-set average is within `0.0012` absolute units of
the MATLAB guide average, and most per-file differences are only a few
thousandths.

For `Venus_Cytoplasm_1.lsm`, `D` is also close under practical masks. The
current implementation gives `D=10.6162` versus the guide's `D=10.7100`, an
absolute difference of `0.0938`.

## Remaining D Gap

The full-set average `D` is not yet close enough for strict historical parity.
The important diagnostic is that the same gap remains when using the
MATLAB-exported FRAP and post-bleach profile vectors directly:

| Input to Python fitter | Avg D | Guide Avg D | Avg corrected MF | Guide Avg corrected MF |
| --- | ---: | ---: | ---: | ---: |
| MATLAB-exported vectors | 42.8221 | 35.6650 | 0.993123 | 0.993123 |

This means the remaining average-`D` discrepancy is not primarily caused by
BioIO loading, ROI rasterization, normalization, photodecay correction, or
post-bleach profile extraction. It is now localized to the historical fitting
procedure: MATLAB-era `lsqcurvefit` behavior, objective weighting, stopping
criteria, or another detail in how the old scripts converged for some datasets.

The largest exported-vector D differences are:

| Dataset | Python D from exported vectors | Guide D | Difference |
| --- | ---: | ---: | ---: |
| `Venus_Cytoplasm_3.lsm` | 74.5390 | 34.7806 | +39.7584 |
| `Venus_Cytoplasm_5.lsm` | 59.0905 | 35.1196 | +23.9709 |
| `Venus_Cytoplasm_9.lsm` | 31.7092 | 24.7671 | +6.9421 |
| `Venus_Cytoplasm_14.lsm` | 43.2788 | 38.5733 | +4.7055 |
| `Venus_Cytoplasm_13.lsm` | 42.5016 | 37.9302 | +4.5714 |

## Practical Recommendation

Keep the current `CircularROI` implementation for the modern app while clearly
documenting that it is not pixel-identical to MATLAB `poly2mask`. It has no
extra dependency, gives a close `D` for the guide's first worked example, and
keeps corrected mobile fractions very close to the guide table.

For publication-grade comparisons to historical MATLAB output, treat corrected
mobile fraction as currently reliable to guide-level tolerance, but keep
diffusion coefficient parity open until the legacy fitting behavior is matched
across the full guide set.
