from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from frap_toolbox_py import cli
from frap_toolbox_py.roi import load_roi_mask
from frap_toolbox_py.types import DiffusionFitResult, FRAPDataset, ReactionFitResult


def _dataset(name: str = "synthetic.lsm") -> FRAPDataset:
    return FRAPDataset(
        name=name,
        time=np.asarray([0.0, 1.0, 2.0]),
        frap=np.asarray([10.0, 6.0, 8.0]),
        norm_frap=np.asarray([1.0, 0.6, 0.8]),
        corrected_frap=np.asarray([1.0, 0.6, 0.8]),
        radius=np.asarray([0.0, 1.0]),
        post_bleach_profile=np.asarray([0.4, 0.8]),
        fit_time=np.asarray([0.0, 1.0]),
        frap_fit=np.asarray([0.6, 0.75]),
        frap_residuals=np.asarray([0.0, 0.01]),
        post_bleach_fit=np.asarray([0.45, 0.75]),
        post_bleach_residuals=np.asarray([-0.05, 0.05]),
        metadata={"image_shape": (6, 6)},
    )


def test_cli_output_dir_writes_diffusion_export_bundle(monkeypatch, tmp_path: Path):
    dataset = _dataset()
    result = DiffusionFitResult(
        datasets=[dataset],
        k=0.7,
        r_effective=1.5,
        half_time=1.2,
        diffusion_coefficient=3.4,
        mobile_fraction=0.8,
        corrected_mobile_fraction=None,
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

    monkeypatch.setattr(cli, "load_diffusion_datasets", lambda *args, **kwargs: [dataset])
    monkeypatch.setattr(cli, "fit_diffusion_model", lambda *args, **kwargs: result)

    output_dir = tmp_path / "exports"
    cli.main(
        [
            "fake.lsm",
            "--roi",
            "3",
            "3",
            "2",
            "--post-bleach-frame",
            "1",
            "--pre-bleach-count",
            "1",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert (output_dir / "result-parameters.csv").exists()
    assert (output_dir / "frap-series.csv").exists()
    assert (output_dir / "post-bleach-profile.csv").exists()
    assert load_roi_mask(output_dir / "roi-masks.npz", "bleach").shape == (6, 6)
    metadata = json.loads((output_dir / "run-metadata.json").read_text(encoding="utf-8"))
    assert metadata["created_by"] == "frap-toolbox-cli"
    assert metadata["model"] == "diffusion"


def test_cli_output_dir_writes_reaction_export_bundle(monkeypatch, tmp_path: Path):
    dataset = _dataset("synthetic.nd2")
    result = ReactionFitResult(
        model_order=1,
        datasets=[dataset],
        parameters={"a": 0.9, "b": 0.4, "c": 0.2},
        decay_rate=0.0,
        sum_squared_residuals=0.02,
        averaged_time=np.asarray([0.0, 1.0]),
        averaged_frap=np.asarray([0.6, 0.8]),
        averaged_frap_fit=np.asarray([0.62, 0.78]),
        averaged_frap_residuals=np.asarray([0.01, -0.01]),
    )

    monkeypatch.setattr(cli, "load_reaction_datasets", lambda *args, **kwargs: [dataset])
    monkeypatch.setattr(cli, "fit_reaction_model", lambda *args, **kwargs: result)

    output_dir = tmp_path / "reaction-exports"
    cli.main(
        [
            "fake.nd2",
            "--model",
            "reaction1",
            "--roi",
            "3",
            "3",
            "2",
            "--post-bleach-frame",
            "1",
            "--pre-bleach-count",
            "1",
            "--fit-mode",
            "individual",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert (output_dir / "result-parameters.csv").exists()
    assert (output_dir / "frap-series.csv").exists()
    assert load_roi_mask(output_dir / "roi-masks.npz", "bleach").shape == (6, 6)
    metadata = json.loads((output_dir / "run-metadata.json").read_text(encoding="utf-8"))
    assert metadata["model"] == "reaction1"
    assert metadata["fit_mode"] == "individual"
