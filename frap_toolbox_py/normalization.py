from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def normalize_frap(
    frap: np.ndarray,
    cell: Optional[np.ndarray],
    post_bleach_frame: int,
    pre_bleach_count: int,
    use_cell_normalization: bool,
    adjacent: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, Optional[float], Optional[np.ndarray]]:
    """Replicates the MATLAB normalization logic for FRAP datasets."""

    start = post_bleach_frame - pre_bleach_count
    end = post_bleach_frame
    if start < 0:
        raise ValueError("Pre-bleach frame window exceeds available frames.")

    if use_cell_normalization:
        if cell is None:
            raise ValueError("Cell intensities required for whole-cell normalization.")
        norm = frap / cell
        baseline = norm[start:end].mean()
        norm_frap = norm / baseline
        return norm_frap, None, None

    baseline = frap[start:end].mean()
    norm_frap = frap / baseline

    if adjacent is None:
        return norm_frap, None, None

    norm_adjacent = adjacent / adjacent[start:end].mean()
    # MATLAB used end - 3 * pre_bleach_count : end, which is inclusive at
    # both ends and therefore maps to one extra sample in Python slicing.
    tail_start = max(post_bleach_frame, len(norm_adjacent) - 3 * pre_bleach_count - 1)
    corrected_mf = 1 - (
        norm_adjacent[tail_start:].mean() - norm_frap[tail_start:].mean()
    )
    return norm_frap, corrected_mf, norm_adjacent
