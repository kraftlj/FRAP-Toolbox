from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pytest
from scipy.optimize import OptimizeWarning

from frap_toolbox_py.data.loading import load_diffusion_datasets
from frap_toolbox_py.models.diffusion import default_diffusion_config, fit_diffusion_model
from frap_toolbox_py.roi import CircularROI, adjacent_circle
from frap_toolbox_py.types import BasicInputs


SAMPLE_PATH = Path("test-data/Diffusion/Venus_Cytoplasm_1.lsm")


def test_lsm_metadata_flow():
    if not SAMPLE_PATH.exists():
        pytest.skip("LSM integration fixture is not present in this checkout.")

    sample_path = SAMPLE_PATH.resolve()
    bleach_roi = CircularROI(256.0, 23.0, 9.0)
    adjacent_roi = adjacent_circle(bleach_roi, offset_factor=2.5)

    inputs = BasicInputs(
        file_paths=[sample_path],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=0.0,
        post_bleach_frame=20,
        roi_definition=(bleach_roi.center_x, bleach_roi.center_y, bleach_roi.radius),
        pre_bleach_frame_count=10,
        use_adjacent_roi=True,
    )

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Could not parse tiff pixel size",
            module="(aicsimageio|bioio_tifffile)",
        )
        datasets = load_diffusion_datasets(
            inputs,
            bleach_roi=bleach_roi.to_mask,
            adjacent_roi=adjacent_roi.to_mask,
        )
    assert datasets, "Expected loader to return at least one dataset."

    dataset = datasets[0]

    # Metadata provenance should record that LSM metadata drove the scaling.
    assert dataset.metadata.get("time_source") == "lsm_metadata"
    assert dataset.metadata.get("voxel_source") == "lsm_metadata"

    # Physical pixel size is reported in microns; Zeiss metadata stores it in meters.
    assert dataset.voxel_size_x == pytest.approx(0.109863, rel=1e-3)
    assert dataset.voxel_size_y == pytest.approx(0.109863, rel=1e-3)

    # Time vector should reflect the native LSM timestamps (≈24 ms spacing).
    assert len(dataset.time) == dataset.frap.shape[0]
    frame_spacing = np.diff(dataset.time[:8]).mean()
    assert frame_spacing == pytest.approx(0.0240, rel=5e-3)

    # Radial profile radii should be expressed in physical units (microns).
    assert dataset.radius[0] == pytest.approx(0.05493, rel=5e-2)

    # End-to-end regression: fitting should return finite parameters.
    config = default_diffusion_config(
        len(dataset.norm_frap),
        len(dataset.radius),
        inputs.post_bleach_frame,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=OptimizeWarning)
        result = fit_diffusion_model(datasets, inputs, config)
    assert np.isfinite(result.k)
    assert np.isfinite(result.r_effective)
    assert np.isfinite(result.diffusion_coefficient)
    assert np.isfinite(result.mobile_fraction)
    assert result.sum_squared_residuals >= 0.0
