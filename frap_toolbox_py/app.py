from __future__ import annotations

from pathlib import Path
from typing import BinaryIO, Union
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning

try:
    import streamlit as st
except ImportError:  # pragma: no cover - app extra is optional for library use
    st = None

try:
    from streamlit_image_coordinates import streamlit_image_coordinates
except ImportError:  # pragma: no cover - app extra is optional for library use
    streamlit_image_coordinates = None

from frap_toolbox_py.data.loading import load_diffusion_datasets, load_reaction_datasets
from frap_toolbox_py.data.loading import _open_image
from frap_toolbox_py.exports import result_parameters_table, write_analysis_bundle
from frap_toolbox_py.models.diffusion import default_diffusion_config, fit_diffusion_model
from frap_toolbox_py.models.reaction import (
    default_reaction1_config,
    default_reaction2_config,
    fit_reaction_model,
)
from frap_toolbox_py.roi import CircularROI, PolygonROI, adjacent_circle, load_roi_mask
from frap_toolbox_py.types import BasicInputs


SUPPORTED_EXTENSIONS = [".lsm", ".tif", ".tiff", ".ome.tif", ".ome.tiff", ".nd2"]
MASK_EXTENSIONS = [".npy", ".npz", ".csv", ".txt", ".tsv"]
MODEL_LABELS = {
    "diffusion": "Diffusion",
    "reaction1": "Reaction 1",
    "reaction2": "Reaction 2",
}
MODEL_INDEX = {
    "diffusion": 1,
    "reaction1": 2,
    "reaction2": 3,
}
DEFAULT_MODEL_DIRS = {
    "diffusion": ("test-data", "Diffusion"),
    "reaction1": ("test-data", "Reaction 1"),
    "reaction2": ("test-data", "Reaction 2"),
}
DIFFUSION_FIT_MODES = [
    "global",
    "individual",
    "average_curve",
    "simplified_kang",
    "simplified_kang_global",
]
REACTION_FIT_MODES = ["individual", "average_curve"]
ROI_SOURCE_CIRCULAR = "Circular numeric ROI"
ROI_SOURCE_DRAW_CIRCLE = "Draw circle on preview"
ROI_SOURCE_DRAW_POLYGON = "Draw polygon on preview"
ROI_SOURCE_SAVED_MASK = "Saved mask upload"
SOFTWARE_AUTHORS = (
    "Lewis J. Kraft, Jacob Dowler, Charles A. Day, Minchul Kang, and Anne K. Kenworthy"
)
SOFTWARE_CITATION = (
    "Kraft LJ, Dowler J, Day CA, Kang M, Kenworthy AK. (2014). "
    "FRAP-Toolbox: Software for the analysis of Fluorescence Recovery After Photobleaching. "
    "http://www.fraptoolbox.com (accessed Month Day, Year)."
)

MaskSource = Union[BinaryIO, Path, str]


def _default_data_dir(model_key: str = "diffusion") -> str:
    parts = DEFAULT_MODEL_DIRS.get(model_key, DEFAULT_MODEL_DIRS["diffusion"])
    candidate = Path.cwd().joinpath(*parts)
    return str(candidate if candidate.exists() else Path.cwd())


def _find_files(directory: Path, extensions: list[str]) -> list[Path]:
    if not directory.exists() or not directory.is_dir():
        return []

    normalized = tuple(ext.lower() for ext in extensions)
    files = []
    for path in directory.iterdir():
        name = path.name.lower()
        if path.is_file() and any(name.endswith(ext) for ext in normalized):
            files.append(path)
    return sorted(files, key=lambda path: path.name.lower())


def _contrast_preview_frame(frame: np.ndarray) -> np.ndarray:
    image = np.asarray(frame, dtype=float)
    if image.ndim != 2:
        raise ValueError("Preview frames must be 2-D.")
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return np.zeros((*image.shape, 3), dtype=np.uint8)

    low, high = np.percentile(finite, [1, 99])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        low = float(np.nanmin(finite))
        high = float(np.nanmax(finite))
    if high <= low:
        scaled = np.zeros_like(image, dtype=np.uint8)
    else:
        scaled = np.clip((image - low) / (high - low), 0.0, 1.0)
        scaled = (scaled * 255).astype(np.uint8)
    return np.repeat(scaled[:, :, None], 3, axis=2)


def _preview_stack_shape(path: Path) -> tuple[int, int, int]:
    image, _ = _open_image(path)
    stack = np.asarray(image.get_image_data("TYX", C=0, Z=0))
    if stack.ndim != 3:
        raise ValueError("Expected stack with dimensions (T, Y, X).")
    return tuple(int(value) for value in stack.shape)  # type: ignore[return-value]


def _preview_frame(path: Path, frame_index: int) -> np.ndarray:
    image, _ = _open_image(path)
    stack = np.asarray(image.get_image_data("TYX", C=0, Z=0))
    if stack.ndim != 3:
        raise ValueError("Expected stack with dimensions (T, Y, X).")
    clipped = min(max(int(frame_index), 0), stack.shape[0] - 1)
    return stack[clipped]


def _clip_zero_based_point(x: float, y: float, image_shape: tuple[int, int]) -> tuple[float, float]:
    height, width = image_shape
    return (
        float(np.clip(x, 0, width - 1)),
        float(np.clip(y, 0, height - 1)),
    )


def _click_to_polygon_point(click: dict[str, object], image_shape: tuple[int, int]) -> tuple[float, float]:
    return _clip_zero_based_point(float(click["x"]), float(click["y"]), image_shape)


def _click_to_circular_roi(click: dict[str, object], radius: float, image_shape: tuple[int, int]) -> CircularROI:
    x, y = _clip_zero_based_point(float(click["x"]), float(click["y"]), image_shape)
    return CircularROI(center_x=x + 1.0, center_y=y + 1.0, radius=radius)


def _add_polygon_point(
    points: list[tuple[float, float]],
    click: dict[str, object],
    image_shape: tuple[int, int],
) -> list[tuple[float, float]]:
    return [*points, _click_to_polygon_point(click, image_shape)]


def _undo_polygon_point(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    return points[:-1]


def _clear_polygon_points() -> list[tuple[float, float]]:
    return []


def _polygon_points_to_mask(points: list[tuple[float, float]], image_shape: tuple[int, int]) -> np.ndarray:
    if len(points) < 3:
        raise ValueError("Draw at least three polygon points before running analysis.")
    mask = PolygonROI(points).to_mask(image_shape)
    if not mask.any():
        raise ValueError("Drawn polygon ROI is empty.")
    return mask


def _blend_mask(rgb: np.ndarray, mask: np.ndarray, color: tuple[int, int, int], alpha: float = 0.35) -> np.ndarray:
    if mask.shape != rgb.shape[:2]:
        return rgb
    output = rgb.astype(float, copy=True)
    color_array = np.asarray(color, dtype=float)
    output[mask] = output[mask] * (1.0 - alpha) + color_array * alpha
    return output.astype(np.uint8)


def _mark_point(rgb: np.ndarray, x: float, y: float, color: tuple[int, int, int], radius: int = 3) -> np.ndarray:
    output = rgb.copy()
    center_x = int(round(x))
    center_y = int(round(y))
    height, width = output.shape[:2]
    y0 = max(0, center_y - radius)
    y1 = min(height, center_y + radius + 1)
    x0 = max(0, center_x - radius)
    x1 = min(width, center_x + radius + 1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    disk_mask = (yy - center_y) ** 2 + (xx - center_x) ** 2 <= radius**2
    patch = output[y0:y1, x0:x1]
    patch[disk_mask] = color
    return output


def _draw_line(rgb: np.ndarray, start: tuple[float, float], end: tuple[float, float], color: tuple[int, int, int]) -> np.ndarray:
    output = rgb.copy()
    x0, y0 = start
    x1, y1 = end
    steps = max(int(abs(x1 - x0)), int(abs(y1 - y0)), 1)
    xs = np.linspace(x0, x1, steps + 1)
    ys = np.linspace(y0, y1, steps + 1)
    height, width = output.shape[:2]
    for x, y in zip(xs, ys):
        col = int(round(x))
        row = int(round(y))
        if 0 <= row < height and 0 <= col < width:
            output[row, col] = color
    return output


def _overlay_roi_preview(
    preview_rgb: np.ndarray,
    *,
    circular_rois: dict[str, CircularROI],
    polygon_points: dict[str, list[tuple[float, float]]],
    closed_polygons: dict[str, bool] | None = None,
    masks: dict[str, np.ndarray],
) -> np.ndarray:
    output = preview_rgb.copy()
    closed_polygons = closed_polygons or {}
    colors = {
        "bleach": (230, 88, 57),
        "cell": (49, 131, 196),
    }
    for role, mask in masks.items():
        output = _blend_mask(output, mask, colors.get(role, (80, 180, 120)))
    for role, roi in circular_rois.items():
        mask = roi.to_mask(output.shape[:2])
        output = _blend_mask(output, mask, colors.get(role, (80, 180, 120)))
        output = _mark_point(output, roi.center_x - 1.0, roi.center_y - 1.0, colors.get(role, (80, 180, 120)))
    for role, points in polygon_points.items():
        color = colors.get(role, (80, 180, 120))
        is_closed = closed_polygons.get(role, False)
        if is_closed and len(points) >= 3:
            try:
                output = _blend_mask(output, _polygon_points_to_mask(points, output.shape[:2]), color, alpha=0.25)
            except ValueError:
                pass
        for point in points:
            output = _mark_point(output, point[0], point[1], color)
        for start, end in zip(points, points[1:]):
            output = _draw_line(output, start, end, color)
        if is_closed and len(points) >= 3:
            output = _draw_line(output, points[-1], points[0], color)
    return output


def _default_fit_mode(model_key: str) -> str:
    if model_key == "reaction1":
        return "individual"
    if model_key == "reaction2":
        return "average_curve"
    return "global"


def _build_inputs(
    files: list[Path],
    model_key: str,
    roi_mode: int,
    roi_definition: tuple[float, ...],
    post_bleach_frame: int,
    pre_bleach_count: int,
    background: float,
    normalize_by_cell: bool,
    use_adjacent_roi: bool,
) -> BasicInputs:
    return BasicInputs(
        file_paths=[path.resolve() for path in files],
        model_index=MODEL_INDEX[model_key],
        roi_mode=roi_mode,
        normalize_by_cell=normalize_by_cell,
        background_intensity=background,
        post_bleach_frame=max(post_bleach_frame - 1, 0),
        roi_definition=roi_definition,
        pre_bleach_frame_count=pre_bleach_count,
        use_adjacent_roi=use_adjacent_roi,
    )


def _reset_mask_source(source: MaskSource) -> None:
    if hasattr(source, "seek"):
        source.seek(0)


def _mask_suffix(source: MaskSource) -> str:
    name = getattr(source, "name", None)
    if name is None:
        name = str(source)
    return Path(name).suffix.lower()


def _coerce_mask(array: np.ndarray, label: str = "mask") -> np.ndarray:
    mask = np.asarray(array)
    if mask.ndim != 2:
        raise ValueError(f"{label} must be a 2-D array.")
    if mask.dtype == bool:
        return mask
    if not np.all(np.isfinite(mask)):
        raise ValueError(f"{label} contains non-finite values.")
    return mask > 0


def _load_mask_array(source: MaskSource, mask_name: str = "mask", label: str = "mask") -> np.ndarray:
    suffix = _mask_suffix(source)
    _reset_mask_source(source)

    if suffix == ".npz":
        try:
            return load_roi_mask(source, mask_name, allow_single=True)
        except ValueError:
            _reset_mask_source(source)
        with np.load(source, allow_pickle=False) as loaded:
            keys = list(loaded.files)
            if mask_name in keys:
                array = loaded[mask_name]
            elif "mask" in keys:
                array = loaded["mask"]
            elif len(keys) == 1:
                array = loaded[keys[0]]
            else:
                raise ValueError(
                    f"{label} NPZ must be a FRAP-Toolbox ROI mask file, "
                    f"contain {mask_name!r}, contain a 'mask' array, or contain exactly one array."
                )
    elif suffix == ".npy":
        array = np.load(source, allow_pickle=False)
    elif suffix == ".csv":
        array = np.loadtxt(source, delimiter=",")
    elif suffix == ".tsv":
        array = np.loadtxt(source, delimiter="\t")
    elif suffix == ".txt":
        array = np.loadtxt(source)
    else:
        raise ValueError(f"Unsupported {label} format. Use one of: {', '.join(MASK_EXTENSIONS)}.")

    _reset_mask_source(source)
    return _coerce_mask(array, label)


def _circular_mask_factory(roi: CircularROI):
    def factory(shape):
        return roi.to_mask(shape)

    return factory


def _build_diffusion_result(
    files: list[Path],
    bleach_roi,
    roi_mode: int,
    roi_definition: tuple[float, ...],
    post_bleach_frame: int,
    pre_bleach_count: int,
    background: float,
    normalize_by_cell: bool,
    cell_roi,
    use_adjacent_roi: bool,
    adjacent_offset: float,
    max_profile_radius: float | None,
    fit_mode: str,
):
    if use_adjacent_roi and roi_mode != 1:
        raise ValueError("Adjacent ROI correction currently requires a circular numeric bleach ROI.")

    adjacent_factory = None
    if use_adjacent_roi:
        adjacent_factory = adjacent_circle(
            CircularROI(*roi_definition),
            offset_factor=adjacent_offset,
        ).to_mask

    inputs = _build_inputs(
        files,
        "diffusion",
        roi_mode,
        roi_definition,
        post_bleach_frame,
        pre_bleach_count,
        background,
        normalize_by_cell,
        use_adjacent_roi,
    )

    datasets = load_diffusion_datasets(
        inputs,
        bleach_roi=bleach_roi,
        cell_roi=cell_roi,
        adjacent_roi=adjacent_factory,
    )
    if not datasets:
        raise ValueError("No datasets were loaded.")

    config = default_diffusion_config(
        frap_length=len(datasets[0].norm_frap),
        profile_length=len(datasets[0].radius),
        post_bleach_frame=inputs.post_bleach_frame,
    )

    if max_profile_radius is not None and np.isfinite(max_profile_radius) and max_profile_radius > 0:
        profile_indices = np.where(datasets[0].radius <= max_profile_radius)[0]
        if profile_indices.size > 2:
            config.profile_range = (int(profile_indices[0]), int(profile_indices[-1]) + 1)
    config.fit_mode = fit_mode

    return fit_diffusion_model(datasets, inputs, config)


def _build_reaction_result(
    files: list[Path],
    model_key: str,
    bleach_roi,
    roi_mode: int,
    roi_definition: tuple[float, ...],
    post_bleach_frame: int,
    pre_bleach_count: int,
    background: float,
    normalize_by_cell: bool,
    cell_roi,
    fit_mode: str,
):
    inputs = _build_inputs(
        files,
        model_key,
        roi_mode,
        roi_definition,
        post_bleach_frame,
        pre_bleach_count,
        background,
        normalize_by_cell,
        False,
    )

    datasets = load_reaction_datasets(
        inputs,
        bleach_roi=bleach_roi,
        cell_roi=cell_roi,
    )
    if not datasets:
        raise ValueError("No datasets were loaded.")

    if model_key == "reaction1":
        config = default_reaction1_config(
            frap_length=len(datasets[0].norm_frap),
            post_bleach_frame=inputs.post_bleach_frame,
        )
        model_order = 1
    else:
        config = default_reaction2_config(
            frap_length=len(datasets[0].norm_frap),
            post_bleach_frame=inputs.post_bleach_frame,
        )
        model_order = 2

    if fit_mode == "individual":
        config.fit_averaged_data = False
    elif fit_mode == "average_curve":
        config.fit_averaged_data = True

    return fit_reaction_model(datasets, inputs, config, model_order=model_order)


def _plot_frap(result):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=[3, 1])
    ax, residual_ax = axes

    for dataset in result.datasets:
        ax.plot(dataset.time, dataset.corrected_frap, marker="o", linestyle="none", markersize=3, alpha=0.45)

    ax.plot(result.averaged_time, result.averaged_frap_fit, color="black", linewidth=2)
    ax.set_ylabel("Normalized intensity")
    ax.grid(True, alpha=0.25)

    residual_ax.axhline(0, color="black", linewidth=1, alpha=0.5)
    residual_ax.plot(result.averaged_time, result.averaged_frap_residuals, marker="o", linestyle="none", markersize=3)
    residual_ax.set_xlabel("Time after bleach (s)")
    residual_ax.set_ylabel("Residual")
    residual_ax.grid(True, alpha=0.25)
    fig.tight_layout()
    return fig


def _plot_profile(result):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=[3, 1])
    ax, residual_ax = axes

    ax.plot(result.averaged_profile_radius, result.averaged_profile, marker="o", linestyle="none", markersize=3)
    ax.plot(result.averaged_profile_radius, result.averaged_profile_fit, color="black", linewidth=2)
    ax.set_ylabel("Normalized intensity")
    ax.grid(True, alpha=0.25)

    residual_ax.axhline(0, color="black", linewidth=1, alpha=0.5)
    residual_ax.plot(result.averaged_profile_radius, result.averaged_profile_residuals, marker="o", linestyle="none", markersize=3)
    residual_ax.set_xlabel("Radial distance (um)")
    residual_ax.set_ylabel("Residual")
    residual_ax.grid(True, alpha=0.25)
    fig.tight_layout()
    return fig


def _result_table(result) -> pd.DataFrame:
    return result_parameters_table(result)


def _render_roi_controls(
    prefix: str,
    default_radius: float,
    *,
    allow_saved_mask: bool = True,
    allow_polygon: bool = True,
):
    key_prefix = prefix.lower()
    options = [ROI_SOURCE_CIRCULAR, ROI_SOURCE_DRAW_CIRCLE]
    if allow_polygon:
        options.append(ROI_SOURCE_DRAW_POLYGON)
    if allow_saved_mask:
        options.append(ROI_SOURCE_SAVED_MASK)

    roi_source = st.radio(
        f"{prefix} ROI source",
        options,
        horizontal=False,
        key=f"{key_prefix}_roi_source",
    )

    if roi_source == ROI_SOURCE_SAVED_MASK:
        uploaded = st.file_uploader(
            f"{prefix} mask",
            type=[ext.removeprefix(".") for ext in MASK_EXTENSIONS],
            key=f"{key_prefix}_mask_upload",
            help=(
                "Accepted formats are FRAP-Toolbox ROI mask NPZ files, NPY, CSV, "
                "TSV, or whitespace-delimited TXT. Generic NPZ files may use "
                "the requested ROI name, a 'mask' array, or one array."
            ),
        )
        return roi_source, None, uploaded

    if roi_source == ROI_SOURCE_DRAW_POLYGON:
        points = st.session_state.get(f"{key_prefix}_polygon_points", [])
        button_col_1, button_col_2, button_col_3 = st.columns(3)
        if button_col_1.button("Undo point", key=f"{key_prefix}_undo_point"):
            st.session_state[f"{key_prefix}_polygon_points"] = _undo_polygon_point(points)
            st.session_state[f"{key_prefix}_polygon_closed"] = False
        if button_col_2.button("Clear ROI", key=f"{key_prefix}_clear_points"):
            st.session_state[f"{key_prefix}_polygon_points"] = _clear_polygon_points()
            st.session_state[f"{key_prefix}_polygon_closed"] = False
        if button_col_3.button("Close ROI", key=f"{key_prefix}_close_points"):
            st.session_state[f"{key_prefix}_polygon_closed"] = len(points) >= 3
        point_count = len(st.session_state.get(f"{key_prefix}_polygon_points", []))
        closed_text = "closed" if st.session_state.get(f"{key_prefix}_polygon_closed", False) else "open"
        st.caption(f"{point_count} polygon points, {closed_text}")
        return roi_source, None, None

    if f"{key_prefix}_center_x" not in st.session_state:
        st.session_state[f"{key_prefix}_center_x"] = 256.0
    if f"{key_prefix}_center_y" not in st.session_state:
        st.session_state[f"{key_prefix}_center_y"] = 23.0 if key_prefix == "bleach" else 256.0
    if f"{key_prefix}_radius" not in st.session_state:
        st.session_state[f"{key_prefix}_radius"] = default_radius
    roi_col_1, roi_col_2 = st.columns(2)
    center_x = roi_col_1.number_input(
        f"{prefix} center X",
        min_value=0.0,
        step=1.0,
        key=f"{key_prefix}_center_x",
    )
    center_y = roi_col_2.number_input(
        f"{prefix} center Y",
        min_value=0.0,
        step=1.0,
        key=f"{key_prefix}_center_y",
    )
    radius = st.number_input(
        f"{prefix} radius",
        min_value=0.1,
        step=0.5,
        key=f"{key_prefix}_radius",
    )
    return roi_source, CircularROI(center_x=center_x, center_y=center_y, radius=radius), None


def _resolve_roi(
    roi_source: str,
    circular_roi: CircularROI | None,
    uploaded_mask,
    mask_name: str,
    label: str,
    *,
    polygon_points: list[tuple[float, float]] | None = None,
    image_shape: tuple[int, int] | None = None,
    polygon_closed: bool = True,
):
    if roi_source == ROI_SOURCE_SAVED_MASK:
        if uploaded_mask is None:
            raise ValueError(f"Upload a {label} mask before running analysis.")
        return _load_mask_array(uploaded_mask, mask_name=mask_name, label=label), 2, ()

    if roi_source == ROI_SOURCE_DRAW_POLYGON:
        if image_shape is None:
            raise ValueError(f"Preview an image before drawing a {label}.")
        if not polygon_closed:
            raise ValueError(f"Close the {label} polygon before running analysis.")
        return _polygon_points_to_mask(polygon_points or [], image_shape), 2, ()

    if circular_roi is None:
        raise ValueError(f"Define a circular {label} ROI before running analysis.")
    return _circular_mask_factory(circular_roi), 1, (circular_roi.center_x, circular_roi.center_y, circular_roi.radius)


def _role_key(role: str) -> str:
    return role.lower()


def _polygon_points_for_role(role: str) -> list[tuple[float, float]]:
    points = st.session_state.get(f"{_role_key(role)}_polygon_points", [])
    return [(float(x), float(y)) for x, y in points]


def _polygon_closed_for_role(role: str) -> bool:
    return bool(st.session_state.get(f"{_role_key(role)}_polygon_closed", False))


def _set_circle_from_click(role: str, click: dict[str, object], radius: float, image_shape: tuple[int, int]) -> None:
    roi = _click_to_circular_roi(click, radius, image_shape)
    key_prefix = _role_key(role)
    st.session_state[f"{key_prefix}_center_x"] = roi.center_x
    st.session_state[f"{key_prefix}_center_y"] = roi.center_y
    st.session_state[f"{key_prefix}_radius"] = roi.radius


def _append_polygon_click(role: str, click: dict[str, object], image_shape: tuple[int, int]) -> None:
    key_prefix = _role_key(role)
    st.session_state[f"{key_prefix}_polygon_points"] = _add_polygon_point(_polygon_points_for_role(role), click, image_shape)
    st.session_state[f"{key_prefix}_polygon_closed"] = False


def _render_preview_and_roi_drawing(
    files: list[Path],
    post_bleach_frame: int,
    normalize_by_cell: bool,
    roi_sources: dict[str, str],
) -> tuple[int | None, tuple[int, int] | None]:
    if not files:
        st.info("Select at least one dataset to preview images and draw ROIs.")
        return None, None

    try:
        frame_count, height, width = _preview_stack_shape(files[0])
    except Exception as exc:
        st.warning(f"Could not load preview image: {exc}")
        return None, None

    default_frame = min(max(int(post_bleach_frame), 1), frame_count)
    frame_number = st.slider(
        "Preview frame",
        min_value=1,
        max_value=frame_count,
        value=default_frame,
        step=1,
    )
    image_shape = (height, width)

    try:
        preview_rgb = _contrast_preview_frame(_preview_frame(files[0], frame_number - 1))
    except Exception as exc:
        st.warning(f"Could not render preview image: {exc}")
        return frame_number - 1, image_shape

    circular_rois: dict[str, CircularROI] = {}
    polygon_points: dict[str, list[tuple[float, float]]] = {}
    closed_polygons: dict[str, bool] = {}
    for role, source in roi_sources.items():
        if source in {ROI_SOURCE_CIRCULAR, ROI_SOURCE_DRAW_CIRCLE}:
            key_prefix = _role_key(role)
            circular_rois[role] = CircularROI(
                center_x=float(st.session_state.get(f"{key_prefix}_center_x", 256.0)),
                center_y=float(st.session_state.get(f"{key_prefix}_center_y", 23.0 if role == "bleach" else 256.0)),
                radius=float(st.session_state.get(f"{key_prefix}_radius", 9.0 if role == "bleach" else 150.0)),
            )
        elif source == ROI_SOURCE_DRAW_POLYGON:
            polygon_points[role] = _polygon_points_for_role(role)
            closed_polygons[role] = _polygon_closed_for_role(role)

    overlay = _overlay_roi_preview(
        preview_rgb,
        circular_rois=circular_rois,
        polygon_points=polygon_points,
        closed_polygons=closed_polygons,
        masks={},
    )

    roles = ["bleach"]
    if normalize_by_cell:
        roles.append("cell")
    active_role = st.radio("Active ROI", roles, horizontal=True, key="active_roi_role")
    active_source = roi_sources.get(active_role, ROI_SOURCE_CIRCULAR)

    click = None
    if streamlit_image_coordinates is None:
        st.image(overlay, caption=files[0].name)
        st.warning(
            "Install the app extra with streamlit-image-coordinates to draw ROIs on the preview."
        )
    else:
        click = streamlit_image_coordinates(overlay, key="roi_preview_coordinates")

    if click and active_source in {ROI_SOURCE_DRAW_CIRCLE, ROI_SOURCE_DRAW_POLYGON}:
        token = (active_role, click.get("x"), click.get("y"), click.get("time"))
        if st.session_state.get("last_roi_click_token") != token:
            if active_source == ROI_SOURCE_DRAW_CIRCLE:
                radius = float(st.session_state.get(f"{_role_key(active_role)}_radius", 9.0))
                _set_circle_from_click(active_role, click, radius, image_shape)
            else:
                _append_polygon_click(active_role, click, image_shape)
            st.session_state["last_roi_click_token"] = token
            st.rerun()

    return frame_number - 1, image_shape


def main() -> None:
    if st is None:
        raise SystemExit(
            "Streamlit is not installed. Install the app extra with "
            '`python -m pip install -e ".[app]"`.'
        )

    st.set_page_config(page_title="FRAP Toolbox", layout="wide")
    st.title("FRAP Toolbox")
    st.caption(f"Authors: {SOFTWARE_AUTHORS}")
    st.caption(f"Citation: {SOFTWARE_CITATION}")

    with st.sidebar:
        st.header("Dataset")
        model_key = st.selectbox(
            "Model",
            list(MODEL_LABELS.keys()),
            format_func=lambda key: MODEL_LABELS[key],
            index=0,
        )
        directory = Path(st.text_input("Directory", value=_default_data_dir(model_key))).expanduser()
        extensions = st.multiselect("Formats", SUPPORTED_EXTENSIONS, default=[".lsm", ".tif", ".tiff", ".nd2"])
        available_files = _find_files(directory, extensions)
        selected_names = st.multiselect(
            "Files",
            [path.name for path in available_files],
            default=[available_files[0].name] if available_files else [],
        )
        selected_files = [path for path in available_files if path.name in selected_names]

        st.header("ROI")
        bleach_source, bleach_circle, bleach_upload = _render_roi_controls(
            "Bleach",
            default_radius=9.0,
            allow_saved_mask=model_key != "diffusion",
            allow_polygon=model_key != "diffusion",
        )

        st.header("Analysis")
        post_bleach_frame = st.number_input("Post-bleach frame", min_value=1, value=21, step=1)
        pre_bleach_count = st.number_input("Pre-bleach frames", min_value=1, value=10, step=1)
        background = st.number_input("Background", min_value=0.0, value=0.0, step=1.0)
        normalize_by_cell = st.checkbox("Whole-cell normalization", value=False)

        cell_source = None
        cell_circle = None
        cell_upload = None
        if normalize_by_cell:
            cell_source, cell_circle, cell_upload = _render_roi_controls("Cell", default_radius=150.0)

        if model_key == "diffusion":
            circular_bleach = bleach_source in {ROI_SOURCE_CIRCULAR, ROI_SOURCE_DRAW_CIRCLE}
            use_adjacent_roi = st.checkbox("Adjacent ROI correction", value=False, disabled=not circular_bleach)
            adjacent_offset = st.slider("Adjacent ROI offset", min_value=1.0, max_value=5.0, value=2.5, step=0.1)
            max_profile_radius = st.number_input("Max profile radius (um)", min_value=0.0, value=0.0, step=0.5)
            fit_mode = st.selectbox(
                "Fit mode",
                DIFFUSION_FIT_MODES,
                index=DIFFUSION_FIT_MODES.index(_default_fit_mode(model_key)),
                help=(
                    "Global concatenates residuals across curves. Average curve reproduces MATLAB's averaged-data "
                    "option. Simplified Kang can run per curve or as a pooled global fit."
                ),
            )
        else:
            use_adjacent_roi = False
            adjacent_offset = 2.5
            max_profile_radius = 0.0
            fit_mode = st.selectbox(
                "Fit mode",
                REACTION_FIT_MODES,
                index=REACTION_FIT_MODES.index(_default_fit_mode(model_key)),
            )

        run_requested = st.button("Run analysis", type="primary", icon=":material/play_arrow:")

    roi_sources = {"bleach": bleach_source}
    if normalize_by_cell and cell_source is not None:
        roi_sources["cell"] = cell_source
    preview_frame_index, preview_image_shape = _render_preview_and_roi_drawing(
        selected_files,
        int(post_bleach_frame),
        normalize_by_cell,
        roi_sources,
    )

    if run_requested:
        if not selected_files:
            st.error("Select at least one dataset.")
            return

        try:
            bleach_roi, roi_mode, roi_definition = _resolve_roi(
                bleach_source,
                bleach_circle,
                bleach_upload,
                "bleach",
                "bleach ROI",
                polygon_points=_polygon_points_for_role("bleach"),
                image_shape=preview_image_shape,
                polygon_closed=_polygon_closed_for_role("bleach"),
            )
            cell_roi = None
            if normalize_by_cell:
                cell_roi, _, _ = _resolve_roi(
                    cell_source,
                    cell_circle,
                    cell_upload,
                    "cell",
                    "cell ROI",
                    polygon_points=_polygon_points_for_role("cell"),
                    image_shape=preview_image_shape,
                    polygon_closed=_polygon_closed_for_role("cell"),
                )

            profile_limit = max_profile_radius if max_profile_radius > 0 else None
            spinner_label = f"Running {MODEL_LABELS[model_key].lower()} fit..."
            with st.spinner(spinner_label):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.filterwarnings("always")
                    warnings.filterwarnings("ignore", category=OptimizeWarning)
                    if model_key == "diffusion":
                        result = _build_diffusion_result(
                            selected_files,
                            bleach_roi,
                            roi_mode,
                            roi_definition,
                            int(post_bleach_frame),
                            int(pre_bleach_count),
                            float(background),
                            normalize_by_cell,
                            cell_roi,
                            use_adjacent_roi,
                            float(adjacent_offset),
                            profile_limit,
                            fit_mode,
                        )
                    else:
                        result = _build_reaction_result(
                            selected_files,
                            model_key,
                            bleach_roi,
                            roi_mode,
                            roi_definition,
                            int(post_bleach_frame),
                            int(pre_bleach_count),
                            float(background),
                            normalize_by_cell,
                            cell_roi,
                            fit_mode,
                        )
        except Exception as exc:
            st.error(str(exc))
            return

        st.session_state["fit_result"] = result
        st.session_state["fit_model"] = model_key
        st.session_state["fit_warnings"] = [str(warning.message) for warning in caught]
        roi_masks = {}
        if preview_image_shape is not None:
            if isinstance(bleach_roi, np.ndarray):
                roi_masks["bleach"] = bleach_roi
            else:
                roi_masks["bleach"] = bleach_roi(preview_image_shape)
            if normalize_by_cell and cell_roi is not None:
                if isinstance(cell_roi, np.ndarray):
                    roi_masks["cell"] = cell_roi
                else:
                    roi_masks["cell"] = cell_roi(preview_image_shape)

        st.session_state["fit_context"] = {
            "model": model_key,
            "files": [path.resolve() for path in selected_files],
            "settings": {
                "post_bleach_frame": int(post_bleach_frame),
                "pre_bleach_count": int(pre_bleach_count),
                "background": float(background),
                "normalize_by_cell": bool(normalize_by_cell),
                "use_adjacent_roi": bool(use_adjacent_roi),
                "adjacent_offset": float(adjacent_offset),
                "max_profile_radius": None if max_profile_radius <= 0 else float(max_profile_radius),
            },
            "roi_sources": roi_sources,
            "fit_mode": fit_mode,
            "roi_masks": roi_masks,
            "roi_extra_metadata": {
                "frame_index": preview_frame_index,
                "roi_source": roi_sources,
                "model": model_key,
                "created_by": "frap-toolbox-app",
            },
        }

    result = st.session_state.get("fit_result")
    if result is None:
        st.info("No fit results yet.")
        return

    metrics = _result_table(result)
    st.dataframe(metrics, hide_index=True, width="stretch")

    if getattr(result, "averaged_profile", None) is None:
        st.subheader("FRAP Fit")
        st.pyplot(_plot_frap(result), clear_figure=True)
    else:
        plot_col_1, plot_col_2 = st.columns(2)
        with plot_col_1:
            st.subheader("FRAP Fit")
            st.pyplot(_plot_frap(result), clear_figure=True)
        with plot_col_2:
            st.subheader("Post-Bleach Profile")
            st.pyplot(_plot_profile(result), clear_figure=True)

    warnings_out = st.session_state.get("fit_warnings", [])
    if warnings_out:
        with st.expander("Warnings"):
            for warning in warnings_out:
                st.write(warning)

    with st.expander("Export"):
        default_output = Path.cwd() / "frap-toolbox-output"
        output_dir = Path(st.text_input("Output directory", value=str(default_output))).expanduser()
        if st.button("Write export bundle"):
            context = st.session_state.get("fit_context")
            if context is None:
                st.error("Run an analysis before exporting.")
            else:
                try:
                    paths = write_analysis_bundle(
                        output_dir,
                        result,
                        model=context["model"],
                        files=context["files"],
                        settings=context["settings"],
                        roi_sources=context["roi_sources"],
                        fit_mode=context["fit_mode"],
                        roi_masks=context["roi_masks"],
                        roi_extra_metadata=context["roi_extra_metadata"],
                    )
                except Exception as exc:
                    st.error(str(exc))
                else:
                    st.success(f"Wrote export bundle to {output_dir}")
                    st.write({name: str(path) for name, path in paths.items()})


if __name__ == "__main__":
    main()
