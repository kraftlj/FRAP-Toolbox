from __future__ import annotations

import numpy as np
import pytest

from frap_toolbox_py.models.diffusion import (
    estimate_recovery_half_time,
    fit_diffusion_model,
    kang_frap,
    simplified_kang_diffusion_coefficient,
    simplified_kang_recovery,
)
from frap_toolbox_py.types import BasicInputs, DiffusionFitConfig, FitBounds, FRAPDataset


def test_kang_frap_monotonic_recovery():
    t = np.linspace(0, 10, 100)
    curve = kang_frap(t, r_effective=3.0, r_nominal=3.0, diffusion=5.0, bleach_depth=1.0)
    assert curve.shape == t.shape
    assert curve[0] < curve[-1]


def test_simplified_kang_diffusion_coefficient_uses_half_time_equation():
    diffusion = simplified_kang_diffusion_coefficient(
        r_nominal=1.0,
        r_effective=np.sqrt(7.0),
        half_time=1.0,
    )
    assert diffusion == pytest.approx(1.0)


def test_simplified_kang_recovery_hits_half_recovery_at_half_time():
    curve = simplified_kang_recovery(
        np.asarray([0.0, 1.0]),
        r_effective=np.sqrt(7.0),
        r_nominal=1.0,
        diffusion=1.0,
        mobile_fraction=1.0,
        initial=0.5,
    )
    np.testing.assert_allclose(curve, [0.5, 0.75], rtol=0, atol=1e-12)


def test_estimate_recovery_half_time_interpolates_midpoint():
    time = np.asarray([0.0, 1.0, 2.0])
    recovery = np.asarray([0.2, 0.5, 1.0])
    assert estimate_recovery_half_time(time, recovery, steady_state=1.0) == pytest.approx(1.2)


def test_simplified_kang_fit_mode_estimates_diffusion_without_curve_fitting():
    time = np.asarray([0.0, 1.0, 2.0])
    recovery = np.asarray([0.2, 0.6, 1.0])
    r_effective = np.sqrt(7.0)
    dataset = FRAPDataset(
        name="synthetic",
        time=time,
        frap=recovery,
        norm_frap=recovery,
        corrected_frap=recovery.copy(),
        radius=np.asarray([0.0]),
        post_bleach_profile=np.asarray([np.exp(-1.0)]),
        voxel_size_x=1.0,
        voxel_size_y=1.0,
    )
    inputs = BasicInputs(
        file_paths=[],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=0.0,
        post_bleach_frame=0,
        roi_definition=(0.0, 0.0, 1.0),
        pre_bleach_frame_count=1,
        use_adjacent_roi=False,
    )
    config = DiffusionFitConfig(
        k=FitBounds(1.0, 1.0, 1.0, "Fixed"),
        r_effective=FitBounds(r_effective, r_effective, r_effective, "Fixed"),
        diffusion_coefficient=FitBounds(10.0, 0.0, np.inf, "Adjustable"),
        mobile_fraction=FitBounds(1.0, 0.0, 2.0, "Adjustable"),
        decay_rate=FitBounds(0.0, 0.0, 0.0, "Fixed"),
        profile_range=(0, 1),
        frap_range=(0, 3),
        decay_fit_range=(0, 3),
        mobile_fraction_range=(2, 3),
        fit_averaged_data=False,
        fit_mode="simplified_kang",
    )

    result = fit_diffusion_model([dataset], inputs, config)

    assert result.half_time == pytest.approx(1.0)
    assert result.diffusion_coefficient == pytest.approx(1.0)
    assert result.mobile_fraction == pytest.approx(1.0)
    assert result.datasets[0].half_time == pytest.approx(1.0)
