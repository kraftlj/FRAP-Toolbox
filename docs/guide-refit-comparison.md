# Guide Refit Comparison

Generated from the MATLAB-exported diffusion vectors in
`test-data/Diffusion`.

Reaction 1 is now refit against the MATLAB-exported guide vectors. Reaction 2
model fitting is implemented, but the current `test-data` archive only contains
the Reaction 2 parameter table, not the FRAP vector export needed for parity
refitting.

## Methods Compared

- **Guide**: MATLAB-exported parameter table values.
- **Legacy parity**: Python refit with the MATLAB-style individual-fit path and
  legacy optimizer emulation.
- **Modern individual**: Python nonlinear Kang-series fit per curve with the
  modern optimizer, then parameter averaging.
- **Modern global**: true global Python nonlinear fit with one shared `K`,
  `re`, `D`, and `MF` across all curves.
- **Modern average curve**: nonlinear fit to the averaged FRAP/profile curves.
- **Simplified Kang**: Traffic 2012 confocal half-time estimator
  `D = (rn^2 + re^2) / (8 * tau_1/2)`, with `tau_1/2` linearly interpolated
  from each corrected FRAP curve and `re` fit from the simplified confocal
  post-bleach profile equation `1 - K * exp(-2r^2 / re^2)`.
- **Simplified Kang global**: true global fit of the simplified profile and
  recovery equations with one shared `K`, `re`, `D`, and `MF`.

`D CV%` and `MF CV%` are sample coefficients of variation across per-file
estimates. For true global and average-curve fits, there is only one group-level
estimate, so CV is either zero or not applicable.

For per-curve `simplified_kang`, `SS` is a diagnostic evaluation of the simplified
confocal recovery curve at the half-time-derived `D`; it is not the minimized
objective, because the simplified estimator is not a nonlinear curve fit.

## Group Summaries

### Venus Cytoplasm

| Mode | re | D | D %diff vs guide | D CV% | MF | MF %diff vs guide | MF CV% | Corrected MF | tau_1/2 | SS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Guide | 3.9839 | 35.6650 | +0.00% | 36.96 | 0.840538 | +0.00% | 6.04 | 0.993123 | NA | 1.40123e-07 |
| Legacy parity | 3.9845 | 35.6441 | -0.06% | 36.93 | 0.840525 | -0.00% | 6.04 | 0.993123 | NA | 1.39809e-07 |
| Modern individual | 3.9843 | 42.9794 | +20.51% | 37.94 | 0.840308 | -0.03% | 6.22 | 0.993123 | NA | 1.35695e-07 |
| Modern global | 3.9719 | 39.7211 | +11.37% | 0.00 | 0.833946 | -0.78% | 0.00 | 0.993123 | NA | 6.34719e-05 |
| Modern average curve | 3.9719 | 39.6361 | +11.13% | NA | 0.838851 | -0.20% | NA | 0.993123 | NA | 1.23068e-07 |
| Simplified Kang | 4.9183 | 40.1459 | +12.56% | 32.17 | 0.830142 | -1.24% | 10.44 | 0.993123 | 0.08506 | 2.50514e-07 |
| Simplified Kang global | 4.9014 | 40.5138 | +13.60% | 0.00 | 0.834755 | -0.69% | 0.00 | 0.993123 | 0.07714 | 6.38256e-05 |

### Venus-ATG5 Cytoplasm

| Mode | re | D | D %diff vs guide | D CV% | MF | MF %diff vs guide | MF CV% | Corrected MF | tau_1/2 | SS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Guide | 3.0176 | 10.3752 | +0.00% | 32.04 | 0.820940 | +0.00% | 5.06 | 0.996604 | NA | 2.49983e-07 |
| Legacy parity | 3.0180 | 10.3781 | +0.03% | 32.04 | 0.820942 | +0.00% | 5.06 | 0.996605 | NA | 2.48679e-07 |
| Modern individual | 3.0178 | 10.4666 | +0.88% | 32.80 | 0.821000 | +0.01% | 5.04 | 0.996605 | NA | 2.58861e-07 |
| Modern global | 3.0170 | 9.8325 | -5.23% | 0.00 | 0.814929 | -0.73% | 0.00 | 0.996605 | NA | 9.36128e-05 |
| Modern average curve | 3.0169 | 9.8691 | -4.88% | NA | 0.819305 | -0.20% | NA | 0.996605 | NA | 2.24960e-07 |
| Simplified Kang | 3.9287 | 10.7484 | +3.60% | 28.52 | 0.807160 | -1.68% | 6.70 | 0.996605 | 0.20372 | 4.94200e-07 |
| Simplified Kang global | 3.9211 | 10.0225 | -3.40% | 0.00 | 0.817062 | -0.47% | 0.00 | 0.996605 | 0.20395 | 9.38110e-05 |

## Per-File Diffusion Coefficients

### Venus Cytoplasm

| File | Guide D | Legacy D | Legacy %diff | Modern individual D | Modern indiv. %diff | Simplified Kang D | Simplified %diff | Simplified re | Simplified tau_1/2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Venus_Cytoplasm_1.lsm` | 10.7100 | 10.7135 | +0.03% | 12.9080 | +20.52% | 17.6545 | +64.84% | 4.4917 | 0.14977 |
| `Venus_Cytoplasm_2.lsm` | 34.0451 | 33.5788 | -1.37% | 37.9232 | +11.39% | 31.3142 | -8.02% | 4.7755 | 0.09494 |
| `Venus_Cytoplasm_3.lsm` | 34.7806 | 34.9065 | +0.36% | 74.2913 | +113.60% | 65.0303 | +86.97% | 5.2754 | 0.05537 |
| `Venus_Cytoplasm_4.lsm` | 49.0692 | 49.0824 | +0.03% | 52.1108 | +6.20% | 48.0686 | -2.04% | 5.4993 | 0.08119 |
| `Venus_Cytoplasm_5.lsm` | 35.1196 | 35.2454 | +0.36% | 58.6534 | +67.01% | 49.2573 | +40.26% | 4.8330 | 0.06176 |
| `Venus_Cytoplasm_6.lsm` | 24.0467 | 24.0511 | +0.02% | 27.9189 | +16.10% | 34.1528 | +42.03% | 4.9629 | 0.09373 |
| `Venus_Cytoplasm_7.lsm` | 37.3884 | 37.4301 | +0.11% | 40.5191 | +8.37% | 29.1530 | -22.03% | 4.8568 | 0.10533 |
| `Venus_Cytoplasm_8.lsm` | 47.7528 | 47.7599 | +0.01% | 50.3619 | +5.46% | 50.0749 | +4.86% | 5.0724 | 0.06667 |
| `Venus_Cytoplasm_9.lsm` | 24.7671 | 24.7758 | +0.04% | 31.7733 | +28.29% | 35.6151 | +43.80% | 5.2603 | 0.10055 |
| `Venus_Cytoplasm_10.lsm` | 31.2611 | 31.2737 | +0.04% | 33.7196 | +7.86% | 30.6954 | -1.81% | 4.8474 | 0.09967 |
| `Venus_Cytoplasm_11.lsm` | 69.1622 | 69.0125 | -0.22% | 71.3889 | +3.22% | 61.9615 | -10.41% | 5.2596 | 0.05778 |
| `Venus_Cytoplasm_12.lsm` | 28.1604 | 28.1630 | +0.01% | 32.3602 | +14.91% | 28.4608 | +1.07% | 4.4687 | 0.09200 |
| `Venus_Cytoplasm_13.lsm` | 37.9302 | 37.9800 | +0.13% | 42.5788 | +12.26% | 37.8997 | -0.08% | 4.4944 | 0.06985 |
| `Venus_Cytoplasm_14.lsm` | 38.5733 | 38.6292 | +0.14% | 43.2895 | +12.23% | 41.3127 | +7.10% | 4.9615 | 0.07744 |
| `Venus_Cytoplasm_15.lsm` | 32.2078 | 32.0599 | -0.46% | 34.8941 | +8.34% | 41.5370 | +28.97% | 4.7155 | 0.06986 |

### Venus-ATG5 Cytoplasm

| File | Guide D | Legacy D | Legacy %diff | Modern individual D | Modern indiv. %diff | Simplified Kang D | Simplified %diff | Simplified re | Simplified tau_1/2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Venus-Atg5_Cytoplasm_1.lsm` | 12.3554 | 12.3593 | +0.03% | 12.5121 | +1.27% | 13.7025 | +10.90% | 4.0997 | 0.16224 |
| `Venus-Atg5_Cytoplasm_2.lsm` | 8.8966 | 8.8983 | +0.02% | 8.9045 | +0.09% | 7.8205 | -12.10% | 3.9290 | 0.26237 |
| `Venus-Atg5_Cytoplasm_3.lsm` | 13.7669 | 13.7690 | +0.02% | 13.7966 | +0.22% | 12.8070 | -6.97% | 3.9750 | 0.16376 |
| `Venus-Atg5_Cytoplasm_4.lsm` | 13.3891 | 13.3951 | +0.04% | 13.4102 | +0.16% | 9.8667 | -26.31% | 3.7739 | 0.19281 |
| `Venus-Atg5_Cytoplasm_5.lsm` | 7.9681 | 7.9743 | +0.08% | 7.9804 | +0.15% | 10.3400 | +29.77% | 3.7556 | 0.18233 |
| `Venus-Atg5_Cytoplasm_6.lsm` | 7.1834 | 7.1847 | +0.02% | 6.8185 | -5.08% | 8.1639 | +13.65% | 3.8142 | 0.23772 |
| `Venus-Atg5_Cytoplasm_7.lsm` | 10.4623 | 10.4660 | +0.04% | 10.8137 | +3.36% | 10.0452 | -3.99% | 3.9709 | 0.20838 |
| `Venus-Atg5_Cytoplasm_8.lsm` | 9.5222 | 9.5231 | +0.01% | 9.5273 | +0.05% | 8.7093 | -8.54% | 3.7077 | 0.21133 |
| `Venus-Atg5_Cytoplasm_9.lsm` | 6.0272 | 6.0295 | +0.04% | 5.9783 | -0.81% | 7.5391 | +25.08% | 3.9448 | 0.27422 |
| `Venus-Atg5_Cytoplasm_10.lsm` | 8.7795 | 8.7822 | +0.03% | 8.7524 | -0.31% | 11.3631 | +29.43% | 3.9995 | 0.18671 |
| `Venus-Atg5_Cytoplasm_11.lsm` | 13.3750 | 13.3780 | +0.02% | 14.2495 | +6.54% | 17.3602 | +29.80% | 4.2813 | 0.13902 |
| `Venus-Atg5_Cytoplasm_12.lsm` | 10.5893 | 10.5924 | +0.03% | 11.0118 | +3.99% | 10.2525 | -3.18% | 3.9560 | 0.20273 |
| `Venus-Atg5_Cytoplasm_13.lsm` | 10.1128 | 10.1150 | +0.02% | 10.1736 | +0.60% | 13.0002 | +28.55% | 4.0051 | 0.16364 |
| `Venus-Atg5_Cytoplasm_14.lsm` | 5.3656 | 5.3668 | +0.02% | 5.2887 | -1.43% | 5.8107 | +8.29% | 3.7526 | 0.32396 |
| `Venus-Atg5_Cytoplasm_15.lsm` | 17.8351 | 17.8381 | +0.02% | 17.7810 | -0.30% | 14.4458 | -19.00% | 3.9657 | 0.14454 |

## Reaction 1 Refit

Reaction 1 uses the MATLAB single-component recovery equation
`F(t) = a - b * exp(-c * t)`. The guide workflow has whole-cell normalization
enabled, so photodecay correction is skipped and the exported normalized FRAP
curves are fit directly. The Python parity path uses the same weighted residual
as `ReactionModel1.m`: `(F(t) - FRAP) / (t + sum(FRAP))`.

| File | Guide a | Python a | Guide b | Python b | Guide c | Python c | Guide SS | Python SS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Venus_1001.nd2` | 0.873656 | 0.873656 | 0.762292 | 0.762292 | 0.0084402 | 0.0084402 | 8.33101e-08 | 8.33105e-08 |
| `Venus_1002.nd2` | 0.855588 | 0.855589 | 0.753077 | 0.753077 | 0.0158993 | 0.0158992 | 6.70117e-08 | 6.70108e-08 |
| `Avg.` | 0.864622 | 0.864622 | 0.757685 | 0.757685 | 0.0121697 | 0.0121697 | 2.09934e-07 | 2.09931e-07 |

## Reaction 1 Raw ROI Diagnostics

The readable `Venus_1002.nd2` raw stack matches the MATLAB-exported time vector,
but the stored ND2 stimulation ROI does not reproduce the exported analysis
ROI. The table below compares raw ND2-derived candidate curves against the
MATLAB-exported `Venus_1002.nd2` vectors.

| Candidate | Target | Pixels | RMSE | MAE | Max abs. error | First post-bleach mean | Target first post-bleach mean | Corr. |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Embedded ND2 stimulation polygon, best common mapping | Raw FRAP | 17,253 | 183.8 | 157.6 | 646.8 | 761.1 | 114.3 | 0.9443 |
| Threshold bleach aid: `pre > 1000` and `post/pre < 0.1` | Raw FRAP | 6,653 | 11.69 | 10.12 | 33.5 | 134.8 | 114.3 | 0.9998 |
| Threshold cell aid: `300 < pre < 2200` and `post/pre < 0.5` | Cell | 5,337 | 26.98 | 22.30 | 76.89 | 446.0 | 518.6 | 0.9858 |

The threshold masks are useful reconstruction diagnostics, not canonical masks:
they are close enough to show that the raw image data are in the right
neighborhood, but they are not export-rounding parity. Exact raw Reaction 1
parity still needs the original hand-drawn MATLAB analysis masks or ROI
vertices. `Venus_1001.nd2` still fails during ND2 metadata parsing with
`Unknown data type in metadata header: 0`, so that file remains a reader issue
rather than a model issue.

## Reaction Guide Status

| Workflow | Guide file(s) | Python refit status |
| --- | --- | --- |
| Reaction 1 | `test-data/Reaction 1/Venus_NCTransport_Reaction_Fit_Parameters.txt`, `*_FRAP_datasets.txt`, raw `Venus_*.nd2` files | Refit-ready from exported vectors; `Venus_1002.nd2` raw timestamps match the MATLAB time export; stimulation ROI and threshold diagnostics are now quantified; exact raw curve parity still needs stored user-defined bleach and cell ROI masks. |
| Reaction 2 | `test-data/Reaction 2/Venus-Atg5_NCTransport_Reaction2_Fit_Parameters.txt`, raw `Venus-Atg5_*.nd2` files | Backend model is ported; both raw ND2 files load with expected stack shape, integer pixels, pixel size, and bleach-relative timestamps; guide curve refit awaits Reaction 2 FRAP vector export or stored ROI masks for the ND2 files. |
