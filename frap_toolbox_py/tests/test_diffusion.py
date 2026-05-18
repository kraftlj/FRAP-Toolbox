from __future__ import annotations

import numpy as np

from frap_toolbox_py.models.diffusion import kang_frap


def test_kang_frap_monotonic_recovery():
    t = np.linspace(0, 10, 100)
    curve = kang_frap(t, r_effective=3.0, r_nominal=3.0, diffusion=5.0, bleach_depth=1.0)
    assert curve.shape == t.shape
    assert curve[0] < curve[-1]
