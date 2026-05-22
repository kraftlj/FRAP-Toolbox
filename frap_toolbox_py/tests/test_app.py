from __future__ import annotations

from io import BytesIO
from pathlib import Path

import numpy as np
import pytest

from frap_toolbox_py.app import (
    SOFTWARE_AUTHORS,
    SOFTWARE_CITATION,
    ROI_SOURCE_DRAW_POLYGON,
    _add_polygon_point,
    _build_inputs,
    _click_to_circular_roi,
    _clear_polygon_points,
    _contrast_preview_frame,
    _default_fit_mode,
    _load_mask_array,
    _polygon_points_to_mask,
    _resolve_roi,
    _undo_polygon_point,
)
from frap_toolbox_py.roi import CircularROI, save_roi_masks


class NamedBytesIO(BytesIO):
    def __init__(self, name: str):
        super().__init__()
        self.name = name


def test_default_fit_modes_match_cli_model_defaults():
    assert _default_fit_mode("diffusion") == "global"
    assert _default_fit_mode("reaction1") == "individual"
    assert _default_fit_mode("reaction2") == "average_curve"


def test_app_attribution_uses_five_author_software_order():
    assert SOFTWARE_AUTHORS == (
        "Lewis J. Kraft, Jacob Dowler, Charles A. Day, Minchul Kang, and Anne K. Kenworthy"
    )
    assert SOFTWARE_CITATION == (
        "Kraft LJ, Dowler J, Day CA, Kang M, Kenworthy AK. (2014). "
        "FRAP-Toolbox: Software for the analysis of Fluorescence Recovery After Photobleaching. "
        "https://github.com/kraftlj/FRAP-Toolbox (accessed Month Day, Year)."
    )


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


def test_load_mask_array_accepts_canonical_roi_mask_npz(tmp_path: Path):
    expected = np.asarray([[False, True], [True, False]], dtype=bool)
    path = tmp_path / "roi.npz"
    save_roi_masks(path, {"bleach": expected}, source_file="fake.nd2")
    buffer = NamedBytesIO("roi.npz")
    buffer.write(path.read_bytes())

    loaded = _load_mask_array(buffer, mask_name="bleach", label="bleach ROI")

    np.testing.assert_array_equal(loaded, expected)


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

    factory, roi_mode, roi_definition = _resolve_roi(roi_source, circular, None, "bleach", "bleach ROI")

    assert roi_mode == 1
    assert roi_definition == (5.0, 6.0, 3.0)
    assert factory((12, 12)).any()


def test_preview_frame_contrast_returns_rgb_uint8_image():
    frame = np.asarray([[0.0, 5.0], [10.0, 20.0]])

    preview = _contrast_preview_frame(frame)

    assert preview.shape == (2, 2, 3)
    assert preview.dtype == np.uint8
    assert preview.max() > preview.min()


def test_click_to_circular_roi_uses_one_based_center_coordinates():
    roi = _click_to_circular_roi({"x": 0, "y": 3}, radius=4.0, image_shape=(10, 10))

    assert roi.center_x == pytest.approx(1.0)
    assert roi.center_y == pytest.approx(4.0)
    assert roi.radius == pytest.approx(4.0)


def test_polygon_point_helpers_add_undo_and_clear_points():
    points = _add_polygon_point([], {"x": 1, "y": 2}, (8, 8))
    points = _add_polygon_point(points, {"x": 10, "y": -2}, (8, 8))

    assert points == [(1.0, 2.0), (7.0, 0.0)]
    assert _undo_polygon_point(points) == [(1.0, 2.0)]
    assert _clear_polygon_points() == []


def test_polygon_points_to_mask_validates_minimum_points_and_shape():
    with pytest.raises(ValueError, match="at least three"):
        _polygon_points_to_mask([(1.0, 1.0), (2.0, 2.0)], (6, 6))

    mask = _polygon_points_to_mask([(1.0, 1.0), (4.0, 1.0), (1.0, 4.0)], (6, 6))

    assert mask.shape == (6, 6)
    assert mask.dtype == np.bool_
    assert mask.any()


def test_resolve_polygon_roi_returns_user_defined_mask():
    mask, roi_mode, roi_definition = _resolve_roi(
        ROI_SOURCE_DRAW_POLYGON,
        None,
        None,
        "bleach",
        "bleach ROI",
        polygon_points=[(1.0, 1.0), (4.0, 1.0), (1.0, 4.0)],
        image_shape=(6, 6),
    )

    assert roi_mode == 2
    assert roi_definition == ()
    assert mask.shape == (6, 6)
    assert mask.any()


def test_resolve_polygon_roi_requires_explicit_close():
    with pytest.raises(ValueError, match="Close the bleach ROI polygon"):
        _resolve_roi(
            ROI_SOURCE_DRAW_POLYGON,
            None,
            None,
            "bleach",
            "bleach ROI",
            polygon_points=[(1.0, 1.0), (4.0, 1.0), (1.0, 4.0)],
            image_shape=(6, 6),
            polygon_closed=False,
        )
