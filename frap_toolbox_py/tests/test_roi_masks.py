from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from frap_toolbox_py import cli
from frap_toolbox_py.roi import load_roi_mask, load_roi_masks, save_roi_mask, save_roi_masks
from frap_toolbox_py.types import FRAPDataset


def test_roi_mask_npz_roundtrip_preserves_masks_and_metadata(tmp_path):
    bleach = np.zeros((4, 5), dtype=bool)
    bleach[1:3, 2:4] = True
    cell = np.ones((4, 5), dtype=bool)
    path = tmp_path / "venus_1002_rois.npz"

    save_roi_masks(
        path,
        {"bleach": bleach, "cell": cell},
        roi_kind="manual",
        source_file="Venus_1002.nd2",
        notes="drawn from legacy user-guide panel",
        extra_metadata={"post_bleach_frame": 6},
    )

    loaded = load_roi_masks(path)

    assert loaded.image_shape == (4, 5)
    assert loaded.metadata["format"] == "frap-toolbox-roi-masks"
    assert loaded.metadata["version"] == 1
    assert loaded.metadata["roi_kind"] == "manual"
    assert loaded.metadata["source_file"] == "Venus_1002.nd2"
    assert loaded.metadata["notes"] == "drawn from legacy user-guide panel"
    assert loaded.metadata["extra"] == {"post_bleach_frame": 6}
    assert {entry["name"] for entry in loaded.metadata["masks"]} == {"bleach", "cell"}
    assert np.array_equal(loaded.require_mask("bleach"), bleach)
    assert np.array_equal(loaded.require_mask("cell"), cell)


def test_roi_mask_helpers_reject_non_bool_masks(tmp_path):
    mask = np.ones((3, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="dtype bool"):
        save_roi_mask(tmp_path / "bad.npz", "bleach", mask)


def test_roi_mask_helpers_reject_shape_mismatches(tmp_path):
    mask = np.ones((3, 3), dtype=bool)
    path = tmp_path / "mask.npz"

    with pytest.raises(ValueError, match="expected"):
        save_roi_mask(path, "bleach", mask, image_shape=(4, 4))

    save_roi_mask(path, "bleach", mask)
    with pytest.raises(ValueError, match="declares image shape"):
        load_roi_mask(path, "bleach", expected_shape=(4, 4))


def test_cli_loads_saved_masks_before_calling_reaction_backend(monkeypatch, tmp_path, capsys):
    bleach = np.zeros((4, 4), dtype=bool)
    bleach[1:3, 1:3] = True
    cell = np.ones((4, 4), dtype=bool)
    mask_path = tmp_path / "reaction_rois.npz"
    save_roi_masks(mask_path, {"bleach": bleach, "cell": cell}, source_file="fake.nd2")

    dataset = FRAPDataset(
        name="fake.nd2",
        time=np.asarray([0.0, 1.0, 2.0]),
        frap=np.asarray([10.0, 6.0, 8.0]),
        norm_frap=np.asarray([1.0, 0.6, 0.8]),
        corrected_frap=np.asarray([1.0, 0.6, 0.8]),
        cell=np.asarray([100.0, 95.0, 90.0]),
    )
    captured = {}

    def fake_load_reaction_datasets(inputs, bleach_roi, cell_roi):
        captured["inputs"] = inputs
        captured["bleach_roi"] = bleach_roi
        captured["cell_roi"] = cell_roi
        return [dataset]

    def fake_fit_reaction_model(datasets, inputs, config, model_order):
        captured["fit_model_order"] = model_order
        return SimpleNamespace(
            parameters={"a": 0.8, "b": 0.2, "c": 0.1},
            decay_rate=0.0,
            sum_squared_residuals=0.0,
        )

    monkeypatch.setattr(cli, "load_reaction_datasets", fake_load_reaction_datasets)
    monkeypatch.setattr(cli, "fit_reaction_model", fake_fit_reaction_model)

    cli.main(
        [
            "fake.nd2",
            "--model",
            "reaction1",
            "--bleach-mask",
            str(mask_path),
            "--cell-mask",
            str(mask_path),
            "--normalize-by-cell",
            "--post-bleach-frame",
            "1",
            "--pre-bleach-count",
            "1",
        ]
    )

    assert captured["inputs"].roi_mode == 2
    assert tuple(captured["inputs"].roi_definition) == ()
    assert np.array_equal(captured["bleach_roi"], bleach)
    assert np.array_equal(captured["cell_roi"], cell)
    assert captured["fit_model_order"] == 1
    assert "Reaction 1 model fit results" in capsys.readouterr().out


def test_cli_rejects_mask_only_diffusion_workflows(tmp_path):
    bleach = np.ones((4, 4), dtype=bool)
    mask_path = tmp_path / "diffusion_bleach.npz"
    save_roi_mask(mask_path, "bleach", bleach)

    with pytest.raises(SystemExit, match="diffusion requires --roi geometry"):
        cli.main(
            [
                "fake.lsm",
                "--model",
                "diffusion",
                "--bleach-mask",
                str(mask_path),
            ]
        )
