from __future__ import annotations

from io import BytesIO
from pathlib import Path

import numpy as np
import pytest

from frap_toolbox_py.app import (
    _build_inputs,
    _default_fit_mode,
    _load_mask_array,
    _resolve_roi,
)
from frap_toolbox_py.roi import CircularROI


class NamedBytesIO(BytesIO):
    def __init__(self, name: str):
        super().__init__()
        self.name = name


def test_default_fit_modes_match_cli_model_defaults():
    assert _default_fit_mode("diffusion") == "global"
    assert _default_fit_mode("reaction1") == "individual"
    assert _default_fit_mode("reaction2") == "average_curve"


def test_build_inputs_uses_zero_based_post_bleach_frame(tmp_path: Path):
    image_path = tmp_path / "example.nd2"
    inputs = _build_inputs(
        [image_path],
        "reaction2",
        roi_mode=2,
        roi_definition=(),
        post_bleach_frame=21,
        pre_bleach_count=10,
        background=3.0,
        normalize_by_cell=True,
        use_adjacent_roi=False,
    )

    assert inputs.model_index == 3
    assert inputs.post_bleach_frame == 20
    assert inputs.background_intensity == pytest.approx(3.0)
    assert inputs.normalize_by_cell is True
    assert inputs.file_paths == [image_path.resolve()]


def test_load_mask_array_accepts_npz_mask_key():
    expected = np.asarray([[0, 1], [2, 0]], dtype=np.uint8)
    buffer = NamedBytesIO("roi.npz")
    np.savez(buffer, mask=expected)

    loaded = _load_mask_array(buffer, label="bleach ROI")

    np.testing.assert_array_equal(loaded, expected > 0)


def test_load_mask_array_accepts_csv_masks():
    buffer = NamedBytesIO("roi.csv")
    buffer.write(b"0,1,0\n1,1,0\n")

    loaded = _load_mask_array(buffer)

    np.testing.assert_array_equal(
        loaded,
        np.asarray([[False, True, False], [True, True, False]]),
    )


def test_load_mask_array_rejects_multi_array_npz_without_mask_key():
    buffer = NamedBytesIO("roi.npz")
    np.savez(buffer, bleach=np.ones((2, 2)), cell=np.ones((2, 2)))

    with pytest.raises(ValueError, match="'mask' array"):
        _load_mask_array(buffer)


def test_resolve_circular_roi_returns_mask_factory_and_cli_definition():
    roi_source = "Circular numeric ROI"
    circular = CircularROI(center_x=5.0, center_y=6.0, radius=3.0)

    factory, roi_mode, roi_definition = _resolve_roi(roi_source, circular, None, "bleach ROI")

    assert roi_mode == 1
    assert roi_definition == (5.0, 6.0, 3.0)
    assert factory((12, 12)).any()
