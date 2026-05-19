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

from frap_toolbox_py.data.loading import load_diffusion_datasets, load_reaction_datasets
from frap_toolbox_py.models.diffusion import default_diffusion_config, fit_diffusion_model
from frap_toolbox_py.models.reaction import (
    default_reaction1_config,
    default_reaction2_config,
    fit_reaction_model,
)
from frap_toolbox_py.roi import CircularROI, adjacent_circle
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


def _load_mask_array(source: MaskSource, label: str = "mask") -> np.ndarray:
    suffix = _mask_suffix(source)
    _reset_mask_source(source)

    if suffix == ".npz":
        with np.load(source, allow_pickle=False) as loaded:
            keys = list(loaded.files)
            if "mask" in keys:
                array = loaded["mask"]
            elif len(keys) == 1:
                array = loaded[keys[0]]
            else:
                raise ValueError(f"{label} NPZ must contain a 'mask' array or exactly one array.")
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
    if hasattr(result, "parameters"):
        rows = [
            {"Parameter": f"Reaction {result.model_order} parameter {name}", "Value": float(value)}
            for name, value in result.parameters.items()
        ]
        rows.extend(
            [
                {"Parameter": "Photodecay rate", "Value": float(result.decay_rate)},
                {"Parameter": "Sum squared residuals", "Value": float(result.sum_squared_residuals)},
            ]
        )
        return pd.DataFrame(rows)

    values = {
        "Bleach depth (k)": result.k,
        "Effective radius (re)": result.r_effective,
        "Half time (tau_1/2)": result.half_time,
        "Diffusion coefficient (D)": result.diffusion_coefficient,
        "Mobile fraction (MF)": result.mobile_fraction,
        "Corrected mobile fraction": result.corrected_mobile_fraction,
        "Photodecay rate": result.decay_rate,
        "Sum squared residuals": result.sum_squared_residuals,
    }
    return pd.DataFrame(
        [{"Parameter": label, "Value": None if value is None else float(value)} for label, value in values.items()]
    )


def _render_roi_controls(prefix: str, default_radius: float, allow_saved_mask: bool = True):
    roi_source = "Circular numeric ROI"
    if allow_saved_mask:
        roi_source = st.radio(
            f"{prefix} ROI source",
            ["Circular numeric ROI", "Saved mask upload"],
            horizontal=True,
            key=f"{prefix.lower()}_roi_source",
        )

    if roi_source == "Saved mask upload":
        uploaded = st.file_uploader(
            f"{prefix} mask",
            type=[ext.removeprefix(".") for ext in MASK_EXTENSIONS],
            key=f"{prefix.lower()}_mask_upload",
            help="Accepted formats are NPY, NPZ with key 'mask', CSV, TSV, or whitespace-delimited TXT.",
        )
        return roi_source, None, uploaded

    roi_col_1, roi_col_2 = st.columns(2)
    center_x = roi_col_1.number_input(f"{prefix} center X", min_value=0.0, value=256.0, step=1.0)
    center_y = roi_col_2.number_input(f"{prefix} center Y", min_value=0.0, value=23.0, step=1.0)
    radius = st.number_input(f"{prefix} radius", min_value=0.1, value=default_radius, step=0.5)
    return roi_source, CircularROI(center_x=center_x, center_y=center_y, radius=radius), None


def _resolve_roi(roi_source: str, circular_roi: CircularROI | None, uploaded_mask, label: str):
    if roi_source == "Saved mask upload":
        if uploaded_mask is None:
            raise ValueError(f"Upload a {label} mask before running analysis.")
        return _load_mask_array(uploaded_mask, label=label), 2, ()

    if circular_roi is None:
        raise ValueError(f"Define a circular {label} ROI before running analysis.")
    return _circular_mask_factory(circular_roi), 1, (circular_roi.center_x, circular_roi.center_y, circular_roi.radius)


def main() -> None:
    if st is None:
        raise SystemExit(
            "Streamlit is not installed. Install the app extra with "
            '`python -m pip install -e ".[app]"`.'
        )

    st.set_page_config(page_title="FRAP Toolbox", layout="wide")
    st.title("FRAP Toolbox")

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
        bleach_source, bleach_circle, bleach_upload = _render_roi_controls("Bleach", default_radius=9.0)

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
            circular_bleach = bleach_source == "Circular numeric ROI"
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

    if run_requested:
        if not selected_files:
            st.error("Select at least one dataset.")
            return

        try:
            bleach_roi, roi_mode, roi_definition = _resolve_roi(
                bleach_source,
                bleach_circle,
                bleach_upload,
                "bleach ROI",
            )
            cell_roi = None
            if normalize_by_cell:
                cell_roi, _, _ = _resolve_roi(
                    cell_source,
                    cell_circle,
                    cell_upload,
                    "cell ROI",
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


if __name__ == "__main__":
    main()
