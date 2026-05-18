from __future__ import annotations

import numpy as np
import pytest

from frap_toolbox_py.normalization import normalize_frap


def test_adjacent_mobile_fraction_window_matches_matlab_inclusive_tail():
    frap = np.linspace(10.0, 20.0, 100)
    adjacent = np.linspace(12.0, 25.0, 100)

    norm_frap, corrected_mf, norm_adjacent = normalize_frap(
        frap=frap,
        cell=None,
        post_bleach_frame=20,
        pre_bleach_count=10,
        use_cell_normalization=False,
        adjacent=adjacent,
    )

    matlab_tail_start = 100 - 3 * 10 - 1
    expected = 1 - (
        norm_adjacent[matlab_tail_start:].mean() - norm_frap[matlab_tail_start:].mean()
    )

    assert corrected_mf == pytest.approx(expected)
    assert len(norm_frap[matlab_tail_start:]) == 31
