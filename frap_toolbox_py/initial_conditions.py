from __future__ import annotations

import numpy as np


def radial_profile(
    image: np.ndarray,
    pre_bleach_stack: np.ndarray,
    center: tuple[float, float],
    clamp: float | None = 2.0,
    invalid_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the normalized radial post-bleach profile.

    Mirrors the MATLAB implementation by averaging measurements at identical
    radial distances rather than binning into arbitrary width buckets. This
    preserves the nanometre-scale structure needed to reproduce the original
    fit parameters.
    """

    mean_pre_bleach = pre_bleach_stack.mean(axis=0)
    normalized = image / np.maximum(mean_pre_bleach, 1e-9)
    if clamp is not None:
        normalized = np.minimum(normalized, clamp)

    # MATLAB treats pixel centers as 1-based coordinates, so mirror that grid here
    # MATLAB's meshgrid(1:width, 1:height) uses xy-indexing (default)
    y_coords = np.arange(1, image.shape[0] + 1, dtype=float)
    x_coords = np.arange(1, image.shape[1] + 1, dtype=float)
    x_indices, y_indices = np.meshgrid(x_coords, y_coords)  # xy-indexing (default)
    radii = np.sqrt((x_indices - center[0]) ** 2 + (y_indices - center[1]) ** 2)

    flat_r = radii.ravel()
    flat_i = normalized.ravel()

    valid_mask = np.isfinite(flat_r) & np.isfinite(flat_i)
    if invalid_mask is not None:
        valid_mask &= ~invalid_mask.ravel()
    flat_r = flat_r[valid_mask]
    flat_i = flat_i[valid_mask]

    unique_radii, inverse = np.unique(flat_r, return_inverse=True)
    sums = np.bincount(inverse, weights=flat_i)
    counts = np.bincount(inverse)
    with np.errstate(invalid="ignore"):
        radial_means = sums / np.maximum(counts, 1)

    return unique_radii, radial_means
