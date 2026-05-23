from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from frap_toolbox_py.exports import (
    frap_series_table,
    post_bleach_profile_table,
    result_parameters_table,
    write_analysis_bundle,
)
from frap_toolbox_py.roi import load_roi_mask
from frap_toolbox_py.types import DiffusionFitResult, FRAPDataset


def _diffusion_result() -> DiffusionFitResult:
    dataset = FRAPDataset(
        name="synthetic.lsm",
        time=np.asarray([0.0, 1.0, 2.0]),
        frap=np.asarray([10.0, 6.0, 8.0]),
        norm_frap=np.asarray([1.0, 0.6, 0.8]),
        corrected_frap=np.asarray([1.0, 0.6, 0.8]),
        cell=np.asarray([100.0, 95.0, 90.0]),
        adjacent=np.asarray([0.9, 0.8, 0.85]),
        radius=np.asarray([0.0, 1.0]),
        post_bleach_profile=np.asarray([0.4, 0.8]),
        fit_time=np.asarray([0.0, 1.0]),
        frap_fit=np.asarray([0.6, 0.75]),
        frap_residuals=np.asarray([0.0, 0.01]),
        post_bleach_fit=np.asarray([0.45, 0.75]),
        post_bleach_residuals=np.asarray([-0.05, 0.05]),
    )
    return DiffusionFitResult(
        datasets=[dataset],
        k=0.7,
        r_effective=1.5,
        half_time=1.2,
        diffusion_coefficient=3.4,
        mobile_fraction=0.8,
        corrected_mobile_fraction=0.82,
        sum_squared_residuals=0.01,
        decay_rate=0.001,
        averaged_profile_radius=np.asarray([0.0, 1.0]),
        averaged_profile=np.asarray([0.4, 0.8]),
        averaged_profile_fit=np.asarray([0.45, 0.75]),
        averaged_profile_residuals=np.asarray([-0.05, 0.05]),
        averaged_time=np.asarray([0.0, 1.0]),
        averaged_frap=np.asarray([0.6, 0.8]),
        averaged_frap_fit=np.asarray([0.6, 0.75]),
        averaged_frap_residuals=np.asarray([0.0, 0.01]),
    )


def test_export_tables_include_expected_semantic_columns():
    result = _diffusion_result()

    parameters = result_parameters_table(result)
    frap = frap_series_table(result)
    profile = post_bleach_profile_table(result)

    assert {"Parameter", "Value"} <= set(parameters.columns)
    assert "Diffusion coefficient (D)" in parameters["Parameter"].tolist()
    assert {"Dataset", "Time", "Raw FRAP", "Corrected FRAP", "FRAP Fit"} <= set(frap.columns)
    assert {"Dataset", "Radius", "Post-bleach Profile", "Profile Fit"} <= set(profile.columns)
    assert "Average" in frap["Dataset"].tolist()


def test_write_analysis_bundle_writes_analysis_export_contract(tmp_path: Path):
    result = _diffusion_result()
    mask = np.zeros((5, 5), dtype=bool)
    mask[1:3, 1:3] = True

    paths = write_analysis_bundle(
        tmp_path,
        result,
        model="diffusion",
        files=[Path("synthetic.lsm")],
        settings={"post_bleach_frame": 21, "pre_bleach_count": 10},
        roi_sources={"bleach": "Circular numeric ROI"},
        fit_mode="global",
        roi_definitions={
            "bleach": {
                "type": "circle",
                "center_x": 3.0,
                "center_y": 3.0,
                "radius": 2.0,
            }
        },
        roi_masks={"bleach": mask},
        roi_extra_metadata={
            "frame_index": 20,
            "roi_source": {"bleach": "Circular numeric ROI"},
            "model": "diffusion",
            "created_by": "frap-toolbox-app",
        },
    )

    expected_files = {
        "result-parameters.csv",
        "frap-series.csv",
        "post-bleach-profile.csv",
        "roi-masks.npz",
        "run-metadata.json",
    }
    assert expected_files <= {path.name for path in paths.values()}
    pd.read_csv(tmp_path / "result-parameters.csv")
    pd.read_csv(tmp_path / "frap-series.csv")
    pd.read_csv(tmp_path / "post-bleach-profile.csv")

    metadata = json.loads((tmp_path / "run-metadata.json").read_text(encoding="utf-8"))
    assert metadata["analysis_bundle_version"] == 1
    assert metadata["model"] == "diffusion"
    assert metadata["fit_mode"] == "global"
    assert metadata["created_by"] == "frap-toolbox-app"
    assert {
        "created_at",
        "created_by",
        "exports",
        "files",
        "fit_mode",
        "model",
        "package_version",
        "roi_sources",
        "settings",
    } <= set(metadata)
    assert metadata["roi_definitions"]["bleach"]["radius"] == 2.0
    assert metadata["exports"]["parameters"] == "result-parameters.csv"
    assert np.array_equal(load_roi_mask(tmp_path / "roi-masks.npz", "bleach"), mask)
