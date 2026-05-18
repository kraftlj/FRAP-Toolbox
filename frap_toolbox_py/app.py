from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning
import streamlit as st

from frap_toolbox_py.data.loading import load_diffusion_datasets
from frap_toolbox_py.models.diffusion import default_diffusion_config, fit_diffusion_model
from frap_toolbox_py.roi import CircularROI, adjacent_circle
from frap_toolbox_py.types import BasicInputs


SUPPORTED_EXTENSIONS = [".lsm", ".tif", ".tiff", ".ome.tif", ".ome.tiff", ".nd2"]


def _default_data_dir() -> str:
    candidate = Path.cwd() / "test-data" / "Diffusion"
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


def _build_diffusion_result(
    files: list[Path],
    roi: CircularROI,
    post_bleach_frame: int,
    pre_bleach_count: int,
    background: float,
    use_adjacent_roi: bool,
    adjacent_offset: float,
    max_profile_radius: float | None,
):
    adjacent_factory = None
    if use_adjacent_roi:
        adjacent_factory = adjacent_circle(roi, offset_factor=adjacent_offset).to_mask

    inputs = BasicInputs(
        file_paths=[path.resolve() for path in files],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=background,
        post_bleach_frame=max(post_bleach_frame - 1, 0),
        roi_definition=(roi.center_x, roi.center_y, roi.radius),
        pre_bleach_frame_count=pre_bleach_count,
        use_adjacent_roi=use_adjacent_roi,
    )

    datasets = load_diffusion_datasets(
        inputs,
        bleach_roi=roi.to_mask,
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

    return fit_diffusion_model(datasets, inputs, config)


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
    values = {
        "Bleach depth (k)": result.k,
        "Effective radius (re)": result.r_effective,
        "Diffusion coefficient (D)": result.diffusion_coefficient,
        "Mobile fraction (MF)": result.mobile_fraction,
        "Corrected mobile fraction": result.corrected_mobile_fraction,
        "Photodecay rate": result.decay_rate,
        "Sum squared residuals": result.sum_squared_residuals,
    }
    return pd.DataFrame(
        [{"Parameter": label, "Value": None if value is None else float(value)} for label, value in values.items()]
    )


def main() -> None:
    st.set_page_config(page_title="FRAP Toolbox", layout="wide")
    st.title("FRAP Toolbox")

    with st.sidebar:
        st.header("Dataset")
        directory = Path(st.text_input("Directory", value=_default_data_dir())).expanduser()
        extensions = st.multiselect("Formats", SUPPORTED_EXTENSIONS, default=[".lsm", ".tif", ".tiff", ".nd2"])
        available_files = _find_files(directory, extensions)
        selected_names = st.multiselect(
            "Files",
            [path.name for path in available_files],
            default=[available_files[0].name] if available_files else [],
        )
        selected_files = [path for path in available_files if path.name in selected_names]

        st.header("ROI")
        roi_col_1, roi_col_2 = st.columns(2)
        center_x = roi_col_1.number_input("Center X", min_value=0.0, value=256.0, step=1.0)
        center_y = roi_col_2.number_input("Center Y", min_value=0.0, value=23.0, step=1.0)
        radius = st.number_input("Radius", min_value=0.1, value=9.0, step=0.5)

        st.header("Analysis")
        post_bleach_frame = st.number_input("Post-bleach frame", min_value=1, value=21, step=1)
        pre_bleach_count = st.number_input("Pre-bleach frames", min_value=1, value=10, step=1)
        background = st.number_input("Background", min_value=0.0, value=0.0, step=1.0)
        use_adjacent_roi = st.checkbox("Adjacent ROI correction", value=True)
        adjacent_offset = st.slider("Adjacent ROI offset", min_value=1.0, max_value=5.0, value=2.5, step=0.1)
        max_profile_radius = st.number_input("Max profile radius (um)", min_value=0.0, value=0.0, step=0.5)

        run_requested = st.button("Run analysis", type="primary", icon=":material/play_arrow:")

    if run_requested:
        if not selected_files:
            st.error("Select at least one dataset.")
            return

        roi = CircularROI(center_x=center_x, center_y=center_y, radius=radius)
        profile_limit = max_profile_radius if max_profile_radius > 0 else None
        with st.spinner("Running diffusion fit..."):
            with warnings.catch_warnings(record=True) as caught:
                warnings.filterwarnings("always")
                warnings.filterwarnings("ignore", category=OptimizeWarning)
                try:
                    result = _build_diffusion_result(
                        selected_files,
                        roi,
                        int(post_bleach_frame),
                        int(pre_bleach_count),
                        float(background),
                        use_adjacent_roi,
                        float(adjacent_offset),
                        profile_limit,
                    )
                except Exception as exc:
                    st.error(str(exc))
                    return

        st.session_state["fit_result"] = result
        st.session_state["fit_warnings"] = [str(warning.message) for warning in caught]

    result = st.session_state.get("fit_result")
    if result is None:
        st.info("No fit results yet.")
        return

    metrics = _result_table(result)
    st.dataframe(metrics, hide_index=True, width="stretch")

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
