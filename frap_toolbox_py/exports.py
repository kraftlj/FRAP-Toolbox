from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .roi import save_roi_masks


def package_version() -> str:
    try:
        return version("frap-toolbox-py")
    except PackageNotFoundError:  # pragma: no cover - editable edge case
        return "0.0.0+local"


def result_parameters_table(result: Any) -> pd.DataFrame:
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


def _array_value(values: Any, index: int) -> float | None:
    if values is None:
        return None
    array = np.asarray(values)
    if index >= array.size:
        return None
    value = array[index]
    if np.isscalar(value) and pd.isna(value):
        return None
    return float(value)


def _append_curve_rows(rows: list[dict[str, Any]], dataset: Any) -> None:
    lengths = [
        len(np.asarray(values))
        for values in [
            dataset.time,
            dataset.frap,
            dataset.norm_frap,
            dataset.corrected_frap,
            dataset.cell,
            dataset.adjacent,
            dataset.fit_time,
            dataset.frap_fit,
            dataset.frap_residuals,
        ]
        if values is not None
    ]
    for index in range(max(lengths or [0])):
        rows.append(
            {
                "Dataset": dataset.name,
                "Frame": index + 1,
                "Time": _array_value(dataset.time, index),
                "Raw FRAP": _array_value(dataset.frap, index),
                "Normalized FRAP": _array_value(dataset.norm_frap, index),
                "Corrected FRAP": _array_value(dataset.corrected_frap, index),
                "Cell": _array_value(dataset.cell, index),
                "Adjacent": _array_value(dataset.adjacent, index),
                "Fit Time": _array_value(dataset.fit_time, index),
                "FRAP Fit": _array_value(dataset.frap_fit, index),
                "FRAP Residual": _array_value(dataset.frap_residuals, index),
            }
        )


def frap_series_table(result: Any) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for dataset in result.datasets:
        _append_curve_rows(rows, dataset)

    averaged_time = getattr(result, "averaged_time", None)
    averaged_frap = getattr(result, "averaged_frap", None)
    averaged_fit = getattr(result, "averaged_frap_fit", None)
    averaged_residuals = getattr(result, "averaged_frap_residuals", None)
    if averaged_time is not None:
        for index in range(len(averaged_time)):
            rows.append(
                {
                    "Dataset": "Average",
                    "Frame": None,
                    "Time": _array_value(averaged_time, index),
                    "Raw FRAP": None,
                    "Normalized FRAP": None,
                    "Corrected FRAP": _array_value(averaged_frap, index),
                    "Cell": None,
                    "Adjacent": None,
                    "Fit Time": _array_value(averaged_time, index),
                    "FRAP Fit": _array_value(averaged_fit, index),
                    "FRAP Residual": _array_value(averaged_residuals, index),
                }
            )

    return pd.DataFrame(rows)


def _append_profile_rows(rows: list[dict[str, Any]], dataset: Any) -> None:
    if dataset.radius is None or dataset.post_bleach_profile is None:
        return
    lengths = [
        len(np.asarray(values))
        for values in [dataset.radius, dataset.post_bleach_profile, dataset.post_bleach_fit, dataset.post_bleach_residuals]
        if values is not None
    ]
    for index in range(max(lengths or [0])):
        rows.append(
            {
                "Dataset": dataset.name,
                "Index": index + 1,
                "Radius": _array_value(dataset.radius, index),
                "Post-bleach Profile": _array_value(dataset.post_bleach_profile, index),
                "Profile Fit": _array_value(dataset.post_bleach_fit, index),
                "Profile Residual": _array_value(dataset.post_bleach_residuals, index),
            }
        )


def post_bleach_profile_table(result: Any) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for dataset in result.datasets:
        _append_profile_rows(rows, dataset)

    averaged_radius = getattr(result, "averaged_profile_radius", None)
    averaged_profile = getattr(result, "averaged_profile", None)
    averaged_fit = getattr(result, "averaged_profile_fit", None)
    averaged_residuals = getattr(result, "averaged_profile_residuals", None)
    if averaged_radius is not None:
        for index in range(len(averaged_radius)):
            rows.append(
                {
                    "Dataset": "Average",
                    "Index": index + 1,
                    "Radius": _array_value(averaged_radius, index),
                    "Post-bleach Profile": _array_value(averaged_profile, index),
                    "Profile Fit": _array_value(averaged_fit, index),
                    "Profile Residual": _array_value(averaged_residuals, index),
                }
            )

    return pd.DataFrame(rows)


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def write_analysis_bundle(
    output_dir: Path | str,
    result: Any,
    *,
    model: str,
    files: Sequence[Path | str],
    settings: Mapping[str, Any],
    roi_sources: Mapping[str, str],
    fit_mode: str,
    roi_masks: Mapping[str, np.ndarray] | None = None,
    roi_extra_metadata: Mapping[str, Any] | None = None,
    created_by: str = "frap-toolbox-app",
) -> dict[str, Path]:
    """Write the local beta analysis export bundle and return created paths."""

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}
    parameters = result_parameters_table(result)
    parameters_path = destination / "result-parameters.csv"
    parameters.to_csv(parameters_path, index=False)
    paths["parameters"] = parameters_path

    frap_series = frap_series_table(result)
    frap_series_path = destination / "frap-series.csv"
    frap_series.to_csv(frap_series_path, index=False)
    paths["frap_series"] = frap_series_path

    profile = post_bleach_profile_table(result)
    if not profile.empty:
        profile_path = destination / "post-bleach-profile.csv"
        profile.to_csv(profile_path, index=False)
        paths["post_bleach_profile"] = profile_path

    roi_mask_path = None
    if roi_masks:
        first_mask = next(iter(roi_masks.values()))
        roi_mask_path = destination / "roi-masks.npz"
        save_roi_masks(
            roi_mask_path,
            roi_masks,
            image_shape=first_mask.shape,
            roi_kind="local-beta",
            source_file=Path(files[0]).name if files else None,
            extra_metadata={"created_by": created_by, **dict(roi_extra_metadata or {})},
        )
        paths["roi_masks"] = roi_mask_path

    metadata = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_by": created_by,
        "package_version": package_version(),
        "model": model,
        "fit_mode": fit_mode,
        "files": [str(Path(path)) for path in files],
        "settings": _json_safe(settings),
        "roi_sources": dict(roi_sources),
        "exports": {name: path.name for name, path in paths.items()},
    }
    if roi_mask_path is not None:
        metadata["roi_mask_file"] = roi_mask_path.name

    metadata_path = destination / "run-metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    paths["metadata"] = metadata_path
    return paths
