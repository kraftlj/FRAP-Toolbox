from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from frap_toolbox_py.models.reaction import (
    default_reaction1_config,
    default_reaction2_config,
    fit_reaction_model,
    reaction1_curve,
    reaction2_curve,
)
from frap_toolbox_py.types import BasicInputs, FRAPDataset


def _reaction_inputs() -> BasicInputs:
    return BasicInputs(
        file_paths=[Path("synthetic")],
        model_index=2,
        roi_mode=2,
        normalize_by_cell=True,
        background_intensity=0.0,
        post_bleach_frame=0,
        roi_definition=(),
        pre_bleach_frame_count=1,
        use_adjacent_roi=False,
    )


def test_reaction1_curve_uses_single_exponential_recovery():
    t = np.asarray([0.0, 1.0])
    curve = reaction1_curve(t, a=1.0, b=0.5, c=np.log(2.0))
    np.testing.assert_allclose(curve, [0.5, 0.75], rtol=0, atol=1e-12)


def test_reaction2_curve_uses_two_exponential_recovery_components():
    t = np.asarray([0.0, 1.0])
    curve = reaction2_curve(t, a=1.0, b=0.25, c=np.log(2.0), d=0.25, f=np.log(4.0))
    np.testing.assert_allclose(curve, [0.5, 0.8125], rtol=0, atol=1e-12)


def test_reaction1_fit_recovers_synthetic_single_component_curve():
    time = np.linspace(0.0, 20.0, 30)
    expected = {"a": 0.9, "b": 0.7, "c": 0.2}
    curve = reaction1_curve(time, **expected)
    dataset = FRAPDataset(
        name="synthetic-r1",
        time=time,
        frap=curve,
        norm_frap=curve,
        corrected_frap=curve.copy(),
    )
    config = default_reaction1_config(len(time), post_bleach_frame=0)
    config.optimizer_mode = "modern"

    result = fit_reaction_model([dataset], _reaction_inputs(), config, model_order=1)

    for name, value in expected.items():
        assert result.parameters[name] == pytest.approx(value, rel=1e-5)
    assert result.sum_squared_residuals < 1e-20


def test_reaction2_average_fit_recovers_synthetic_two_component_curve():
    time = np.linspace(0.0, 250.0, 80)
    expected = {"a": 1.2, "b": 0.35, "c": 0.04, "d": 0.45, "f": 0.004}
    curve = reaction2_curve(time, **expected)
    datasets = [
        FRAPDataset(
            name=f"synthetic-r2-{idx}",
            time=time,
            frap=curve,
            norm_frap=curve,
            corrected_frap=curve.copy(),
        )
        for idx in range(2)
    ]
    config = default_reaction2_config(len(time), post_bleach_frame=0)
    config.optimizer_mode = "modern"
    config.fit_averaged_data = True

    result = fit_reaction_model(datasets, _reaction_inputs(), config, model_order=2)

    for name, value in expected.items():
        assert result.parameters[name] == pytest.approx(value, rel=1e-4)
    assert result.sum_squared_residuals < 1e-20
