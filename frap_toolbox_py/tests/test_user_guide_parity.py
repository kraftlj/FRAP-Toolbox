from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pytest

from frap_toolbox_py.data.loading import _extract_time_vector, _open_image, load_diffusion_datasets
from frap_toolbox_py.models.diffusion import fit_diffusion_model, kang_frap
from frap_toolbox_py.models.reaction import fit_reaction_model
from frap_toolbox_py.roi import CircularROI, adjacent_circle
from frap_toolbox_py.types import BasicInputs, DiffusionFitConfig, FitBounds, FRAPDataset, ReactionFitConfig


FIXTURE_ROOT = Path("test-data")
DIFFUSION_ROOT = FIXTURE_ROOT / "Diffusion"
REACTION1_ROOT = FIXTURE_ROOT / "Reaction 1"
REACTION2_ROOT = FIXTURE_ROOT / "Reaction 2"
VENUS_CYTO_1 = "Venus_Cytoplasm_1.lsm"
POST_BLEACH_INDEX = 20


def _require_user_guide_fixtures(*paths: Path) -> None:
    required = (FIXTURE_ROOT / "Userguide.pdf", *paths)
    missing = [path for path in required if not path.exists()]
    if missing:
        pytest.skip("User guide parity fixtures are not present in this checkout.")


def _read_parameter_table(path: Path) -> dict[str, dict[str, float]]:
    rows = path.read_text().splitlines()
    headers = [cell.strip() for cell in rows[0].split("\t")]
    table: dict[str, dict[str, float]] = {}
    for row in rows[1:]:
        cells = [cell.strip() for cell in row.split("\t")]
        if not cells or not cells[0]:
            continue
        values: dict[str, float] = {}
        for header, cell in zip(headers[1:], cells[1:]):
            values[header] = float(cell) if cell else np.nan
        table[cells[0]] = values
    return table


def _read_vector_export(path: Path) -> dict[tuple[str, str], np.ndarray]:
    vectors: dict[tuple[str, str], np.ndarray] = {}
    for row in path.read_text().splitlines():
        cells = row.split("\t")
        if len(cells) < 3:
            continue
        filename = cells[0].strip()
        label = cells[1].strip()
        values = [float(cell) for cell in cells[2:] if cell.strip()]
        vectors[(filename, label)] = np.asarray(values, dtype=float)
    return vectors


def _open_image_or_skip(path: Path):
    try:
        return _open_image(path)
    except ImportError as exc:
        pytest.skip(f"Microscopy reader for {path.name} is not installed: {exc}")


def _curve_error_metrics(curve: np.ndarray, target: np.ndarray) -> dict[str, float]:
    diff = np.asarray(curve, dtype=float) - np.asarray(target, dtype=float)
    return {
        "rmse": float(np.sqrt(np.mean(diff * diff))),
        "mae": float(np.mean(np.abs(diff))),
        "max_abs": float(np.max(np.abs(diff))),
        "pre_mean": float(np.mean(curve[:5])),
        "post": float(curve[5]),
        "tail_mean": float(np.mean(curve[-10:])),
        "corr": float(np.corrcoef(curve, target)[0, 1]),
    }


def _reaction1_venus_1002_stack_and_exports():
    frap_export = REACTION1_ROOT / "Venus_NCTransport_Reaction_FRAP_datasets.txt"
    raw_path = REACTION1_ROOT / "Venus_1002.nd2"
    _require_user_guide_fixtures(frap_export, raw_path)

    frap_exports = _read_vector_export(frap_export)
    image, _ = _open_image_or_skip(raw_path)
    stack = np.asarray(image.get_image_data("TYX", C=0, Z=0), dtype=float)
    return raw_path, stack, frap_exports


def _nd2_stimulation_polygon_candidate_masks(path: Path, shape: tuple[int, int]) -> dict[str, np.ndarray]:
    nd2 = pytest.importorskip("nd2")
    from skimage.draw import polygon

    with nd2.ND2File(str(path)) as nd2_file:
        rois = list(nd2_file.rois.values())

    if not rois:
        pytest.skip(f"{path.name} has no ND2 ROI metadata.")

    points = rois[0].animParams[0].extrudedShape.basePoints
    x = np.asarray([point.x for point in points], dtype=float)
    y = np.asarray([point.y for point in points], dtype=float)
    height, width = shape

    coordinate_sets = {
        "xy_plus": ((x + 0.5) * (width - 1), (y + 0.5) * (height - 1)),
        "xy_flip": ((x + 0.5) * (width - 1), (0.5 - y) * (height - 1)),
        "xy_plus_512": ((x + 0.5) * width, (y + 0.5) * height),
        "xy_flip_512": ((x + 0.5) * width, (0.5 - y) * height),
        "swap_plus": ((y + 0.5) * (width - 1), (x + 0.5) * (height - 1)),
        "swap_flip": ((0.5 - y) * (width - 1), (x + 0.5) * (height - 1)),
    }

    masks = {}
    for label, (cols, rows) in coordinate_sets.items():
        rr, cc = polygon(rows, cols, shape=shape)
        mask = np.zeros(shape, dtype=bool)
        mask[rr, cc] = True
        masks[label] = mask
    return masks


def _threshold_bleach_candidate(stack: np.ndarray) -> np.ndarray:
    pre_bleach = stack[:5].mean(axis=0)
    post_bleach = stack[5]
    ratio = np.divide(
        post_bleach,
        pre_bleach,
        out=np.ones_like(post_bleach, dtype=float),
        where=pre_bleach > 0,
    )
    return (pre_bleach > 1000.0) & (ratio < 0.1)


def _threshold_cell_candidate(stack: np.ndarray) -> np.ndarray:
    pre_bleach = stack[:5].mean(axis=0)
    post_bleach = stack[5]
    ratio = np.divide(
        post_bleach,
        pre_bleach,
        out=np.ones_like(post_bleach, dtype=float),
        where=pre_bleach > 0,
    )
    return (pre_bleach > 300.0) & (pre_bleach < 2200.0) & (ratio < 0.5)


def _load_venus_cytoplasm_1():
    sample_path = (DIFFUSION_ROOT / VENUS_CYTO_1).resolve()
    bleach_roi = CircularROI(256.0, 23.0, 9.0)
    adjacent_roi = adjacent_circle(bleach_roi, offset_factor=2.5)
    inputs = BasicInputs(
        file_paths=[sample_path],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=0.0,
        post_bleach_frame=POST_BLEACH_INDEX,
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
        return load_diffusion_datasets(
            inputs,
            bleach_roi=bleach_roi.to_mask,
            adjacent_roi=adjacent_roi.to_mask,
        )[0]


def _inferred_legacy_venus_cytoplasm_1_bleach_mask(shape: tuple[int, int]) -> np.ndarray:
    """Mask inferred from the MATLAB Raw FRAP export for the guide fixture."""

    row_intervals = {
        14: (253, 259),
        15: (252, 260),
        16: (250, 262),
        17: (249, 263),
        18: (249, 263),
        19: (248, 264),
        20: (248, 264),
        21: (248, 264),
        22: (248, 265),
        23: (248, 265),
        24: (248, 264),
        25: (248, 264),
        26: (248, 264),
        27: (249, 263),
        28: (250, 262),
        29: (250, 262),
        30: (252, 260),
        31: (254, 258),
    }
    mask = np.zeros(shape, dtype=bool)
    for row, (start_col, end_col) in row_intervals.items():
        mask[row - 1, start_col - 1 : end_col] = True
    return mask


def _load_venus_cytoplasm_1_with_inferred_bleach_mask():
    sample_path = (DIFFUSION_ROOT / VENUS_CYTO_1).resolve()
    bleach_roi = CircularROI(256.0, 23.0, 9.0)
    adjacent_roi = adjacent_circle(bleach_roi, offset_factor=2.5)
    inputs = BasicInputs(
        file_paths=[sample_path],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=0.0,
        post_bleach_frame=POST_BLEACH_INDEX,
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
        return load_diffusion_datasets(
            inputs,
            bleach_roi=_inferred_legacy_venus_cytoplasm_1_bleach_mask,
            adjacent_roi=adjacent_roi.to_mask,
        )[0]


def _guide_diffusion_config(
    frap_exports: dict[tuple[str, str], np.ndarray],
    profile_exports: dict[tuple[str, str], np.ndarray],
    first_name: str = VENUS_CYTO_1,
) -> DiffusionFitConfig:
    return DiffusionFitConfig(
        k=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        r_effective=FitBounds(3.0, 0.0, np.inf, "Adjustable"),
        diffusion_coefficient=FitBounds(10.0, 0.0, np.inf, "Adjustable"),
        mobile_fraction=FitBounds(1.0, 0.0, 2.0, "Adjustable"),
        decay_rate=FitBounds(1e-3, 0.0, np.inf, "Adjustable"),
        profile_range=(0, len(profile_exports[(first_name, "Distance")])),
        frap_range=(POST_BLEACH_INDEX, len(frap_exports[(first_name, "Time")])),
        decay_fit_range=(499, 600),
        mobile_fraction_range=(539, 600),
        fit_averaged_data=False,
        optimizer_mode="legacy_matlab",
    )


def _build_diffusion_export_datasets(
    params: dict[str, dict[str, float]],
    frap_exports: dict[tuple[str, str], np.ndarray],
    profile_exports: dict[tuple[str, str], np.ndarray],
    *,
    use_corrected_as_normalized: bool = False,
) -> list[FRAPDataset]:
    datasets = []
    names = sorted(
        (name for name in params if name.endswith(".lsm")),
        key=lambda name: int(name.rsplit("_", 1)[1].split(".")[0]),
    )
    for name in names:
        time = frap_exports[(name, "Time")]
        norm_frap = (
            frap_exports[(name, "Corrected FRAP")].copy()
            if use_corrected_as_normalized
            else frap_exports[(name, "Normalized FRAP")]
        )
        datasets.append(
            FRAPDataset(
                name=name,
                time=time - time[0],
                frap=frap_exports[(name, "Raw FRAP")],
                norm_frap=norm_frap,
                corrected_frap=norm_frap.copy(),
                corrected_mobile_fraction=params[name]["MF Corrected"],
                voxel_size_x=0.10986328426900892,
                voxel_size_y=0.10986328426900892,
                radius=profile_exports[(name, "Distance")],
                post_bleach_profile=profile_exports[(name, "Post-bleach Profile")],
            )
        )
    return datasets


def _build_venus_cytoplasm_export_datasets(
    params: dict[str, dict[str, float]],
    frap_exports: dict[tuple[str, str], np.ndarray],
    profile_exports: dict[tuple[str, str], np.ndarray],
) -> list[FRAPDataset]:
    return _build_diffusion_export_datasets(params, frap_exports, profile_exports)


def _reference_diffusion_inputs() -> BasicInputs:
    return BasicInputs(
        file_paths=[Path("reference-export")],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=0.0,
        post_bleach_frame=POST_BLEACH_INDEX,
        roi_definition=(256.0, 23.0, 9.0),
        pre_bleach_frame_count=10,
        use_adjacent_roi=True,
    )


def _reference_reaction_inputs(model_index: int = 2) -> BasicInputs:
    return BasicInputs(
        file_paths=[Path("reference-export")],
        model_index=model_index,
        roi_mode=2,
        normalize_by_cell=True,
        background_intensity=0.0,
        post_bleach_frame=5,
        roi_definition=(),
        pre_bleach_frame_count=5,
        use_adjacent_roi=False,
    )


def _guide_reaction1_config() -> ReactionFitConfig:
    return ReactionFitConfig(
        a=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        b=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        c=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        decay_rate=FitBounds(1e-3, 0.0, np.inf, "Adjustable"),
        frap_range=(5, 130),
        decay_fit_range=(134, 185),
        fit_averaged_data=False,
        optimizer_mode="legacy_matlab",
    )


def _build_reaction_export_datasets(
    frap_exports: dict[tuple[str, str], np.ndarray],
    names: list[str],
) -> list[FRAPDataset]:
    return [
        FRAPDataset(
            name=name,
            time=frap_exports[(name, "Time")],
            frap=frap_exports[(name, "Raw FRAP")],
            cell=frap_exports[(name, "Cell")],
            norm_frap=frap_exports[(name, "Normalized FRAP")],
            corrected_frap=frap_exports[(name, "Corrected FRAP")].copy(),
        )
        for name in names
    ]


def test_user_guide_reference_parameter_tables_are_available():
    venus_params = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Fit_Parameters.txt"
    atg5_params = DIFFUSION_ROOT / "Venus-ATG5_Cytoplasm_Diffusion_Fit_Parameters.txt"
    reaction1_params = REACTION1_ROOT / "Venus_NCTransport_Reaction_Fit_Parameters.txt"
    reaction2_params = REACTION2_ROOT / "Venus-Atg5_NCTransport_Reaction2_Fit_Parameters.txt"
    _require_user_guide_fixtures(venus_params, atg5_params, reaction1_params, reaction2_params)

    venus = _read_parameter_table(venus_params)
    atg5 = _read_parameter_table(atg5_params)
    reaction1 = _read_parameter_table(reaction1_params)
    reaction2 = _read_parameter_table(reaction2_params)

    assert venus[VENUS_CYTO_1]["k"] == pytest.approx(1.49227)
    assert venus[VENUS_CYTO_1]["re"] == pytest.approx(3.52497)
    assert venus[VENUS_CYTO_1]["D"] == pytest.approx(10.71)
    assert venus[VENUS_CYTO_1]["MF"] == pytest.approx(0.778875)
    assert venus[VENUS_CYTO_1]["MF Corrected"] == pytest.approx(0.944335)
    assert venus["Avg."]["D"] == pytest.approx(35.665)

    assert atg5["Venus-Atg5_Cytoplasm_1.lsm"]["k"] == pytest.approx(2.13326)
    assert atg5["Avg."]["D"] == pytest.approx(10.3752)
    assert atg5["Avg."]["MF Corrected"] == pytest.approx(0.996604)

    assert reaction1["Venus_1001.nd2"]["a"] == pytest.approx(0.873656)
    assert reaction1["Venus_1002.nd2"]["c"] == pytest.approx(0.0158993)
    assert reaction1["Avg."]["SS"] == pytest.approx(2.09934e-07)

    assert np.isnan(reaction2["Venus-Atg5_1002.nd2"]["a"])
    assert np.isnan(reaction2["Venus-Atg5_1003.nd2"]["f"])
    assert reaction2["Avg."]["a"] == pytest.approx(1.4704)
    assert reaction2["Avg."]["f"] == pytest.approx(0.000442516)


def test_reaction1_nd2_time_metadata_matches_matlab_export_for_venus_1002():
    frap_export = REACTION1_ROOT / "Venus_NCTransport_Reaction_FRAP_datasets.txt"
    raw_path = REACTION1_ROOT / "Venus_1002.nd2"
    _require_user_guide_fixtures(frap_export, raw_path)

    frap_exports = _read_vector_export(frap_export)
    image, _ = _open_image_or_skip(raw_path)
    stack = np.asarray(image.get_image_data("TYX", C=0, Z=0), dtype=float)
    times, time_source = _extract_time_vector(image, stack.shape[0])

    assert stack.shape == (185, 512, 512)
    assert time_source == "ome_planes"
    np.testing.assert_allclose(
        times - times[5],
        frap_exports[("Venus_1002.nd2", "Time")],
        rtol=0,
        atol=1e-6,
    )


def test_reaction1_venus_1001_nd2_failure_is_reader_metadata_not_model_logic():
    raw_path = REACTION1_ROOT / "Venus_1001.nd2"
    _require_user_guide_fixtures(raw_path)

    try:
        image, _ = _open_image(raw_path)
    except ImportError as exc:
        if "Unknown data type in metadata header: 0" not in str(exc):
            pytest.skip(f"ND2 reader is not available for {raw_path.name}: {exc}")
        assert "Unknown data type in metadata header: 0" in str(exc)
        return

    try:
        stack = np.asarray(image.get_image_data("TYX", C=0, Z=0))
    except ValueError as exc:
        assert "Unknown data type in metadata header: 0" in str(exc)
        return

    assert stack.shape == (185, 512, 512)


def test_reaction1_venus_1002_embedded_stimulation_roi_does_not_reproduce_raw_export():
    raw_path, stack, frap_exports = _reaction1_venus_1002_stack_and_exports()
    target = frap_exports[("Venus_1002.nd2", "Raw FRAP")]
    masks = _nd2_stimulation_polygon_candidate_masks(raw_path, stack.shape[1:])

    candidates = []
    for label, mask in masks.items():
        curve = stack[:, mask].mean(axis=1)
        metrics = _curve_error_metrics(curve, target)
        candidates.append((metrics["rmse"], label, int(mask.sum()), metrics))

    best_rmse, best_label, best_pixels, best_metrics = min(candidates)

    assert best_label == "xy_plus"
    assert best_pixels == 17253
    assert best_rmse == pytest.approx(183.8, abs=0.5)
    assert best_metrics["post"] == pytest.approx(761.1, abs=0.5)
    assert abs(best_metrics["post"] - target[5]) > 600.0


def test_reaction1_venus_1002_threshold_inferred_bleach_mask_is_close_but_not_exact():
    _, stack, frap_exports = _reaction1_venus_1002_stack_and_exports()
    target = frap_exports[("Venus_1002.nd2", "Raw FRAP")]
    mask = _threshold_bleach_candidate(stack)
    curve = stack[:, mask].mean(axis=1)
    metrics = _curve_error_metrics(curve, target)

    assert int(mask.sum()) == 6653
    assert metrics["rmse"] == pytest.approx(11.69, abs=0.1)
    assert metrics["mae"] == pytest.approx(10.12, abs=0.1)
    assert metrics["max_abs"] == pytest.approx(33.5, abs=0.2)
    assert metrics["corr"] > 0.9998
    assert metrics["max_abs"] > 10.0


def test_reaction1_venus_1002_threshold_inferred_cell_mask_is_close_but_not_exact():
    _, stack, frap_exports = _reaction1_venus_1002_stack_and_exports()
    target = frap_exports[("Venus_1002.nd2", "Cell")]
    mask = _threshold_cell_candidate(stack)
    curve = stack[:, mask].mean(axis=1)
    metrics = _curve_error_metrics(curve, target)

    assert int(mask.sum()) == 5337
    assert metrics["rmse"] == pytest.approx(26.98, abs=0.1)
    assert metrics["mae"] == pytest.approx(22.30, abs=0.1)
    assert metrics["max_abs"] == pytest.approx(76.89, abs=0.2)
    assert metrics["corr"] > 0.985
    assert metrics["max_abs"] > 50.0


@pytest.mark.parametrize("raw_name", ["Venus-Atg5_1002.nd2", "Venus-Atg5_1003.nd2"])
def test_reaction2_nd2_time_metadata_uses_seconds_after_bleach(raw_name: str):
    raw_path = REACTION2_ROOT / raw_name
    _require_user_guide_fixtures(raw_path)

    image, _ = _open_image_or_skip(raw_path)
    raw_stack = image.get_image_data("TYX", C=0, Z=0)
    stack = np.asarray(raw_stack, dtype=float)
    times, time_source = _extract_time_vector(image, stack.shape[0])
    bleach_relative_time = times - times[5]

    pixel_sizes = image.physical_pixel_sizes
    assert np.issubdtype(raw_stack.dtype, np.integer)
    assert stack.shape == (185, 512, 512)
    assert time_source == "ome_planes"
    assert pixel_sizes.X == pytest.approx(0.15537014235055974)
    assert pixel_sizes.Y == pytest.approx(0.15537014235055974)
    assert stack[:5].mean() > stack[5].mean()
    assert bleach_relative_time[5] == pytest.approx(0.0, abs=1e-12)
    assert bleach_relative_time[6] == pytest.approx(10.0, abs=0.1)
    assert bleach_relative_time[-1] == pytest.approx(1790.0, abs=0.1)


def test_reaction1_legacy_fit_reproduces_user_guide_exported_vectors():
    params_path = REACTION1_ROOT / "Venus_NCTransport_Reaction_Fit_Parameters.txt"
    frap_export = REACTION1_ROOT / "Venus_NCTransport_Reaction_FRAP_datasets.txt"
    _require_user_guide_fixtures(params_path, frap_export)

    params = _read_parameter_table(params_path)
    frap_exports = _read_vector_export(frap_export)
    names = ["Venus_1001.nd2", "Venus_1002.nd2"]
    datasets = _build_reaction_export_datasets(frap_exports, names)

    result = fit_reaction_model(
        datasets,
        _reference_reaction_inputs(model_index=2),
        _guide_reaction1_config(),
        model_order=1,
    )

    for dataset in result.datasets:
        guide = params[dataset.name]
        assert dataset.reaction_parameters["a"] == pytest.approx(guide["a"], abs=1e-6)
        assert dataset.reaction_parameters["b"] == pytest.approx(guide["b"], abs=1e-6)
        assert dataset.reaction_parameters["c"] == pytest.approx(guide["c"], abs=1e-7)
        assert dataset.sum_squared_residuals == pytest.approx(guide["SS"], rel=2e-5)
        np.testing.assert_allclose(
            dataset.frap_fit,
            frap_exports[(dataset.name, "FRAP Fit")],
            rtol=0,
            atol=1e-6,
        )
        np.testing.assert_allclose(
            dataset.frap_residuals,
            frap_exports[(dataset.name, "Fit Residuals")],
            rtol=0,
            atol=1e-6,
        )

    guide_average = params["Avg."]
    assert result.parameters["a"] == pytest.approx(guide_average["a"], abs=1e-6)
    assert result.parameters["b"] == pytest.approx(guide_average["b"], abs=1e-6)
    assert result.parameters["c"] == pytest.approx(guide_average["c"], abs=1e-7)
    assert result.sum_squared_residuals == pytest.approx(guide_average["SS"], rel=5e-5)
    np.testing.assert_allclose(
        result.averaged_frap_fit,
        frap_exports[("Average", "FRAP Fit")],
        rtol=0,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        result.averaged_frap_residuals,
        frap_exports[("Average", "Fit Residuals")],
        rtol=0,
        atol=1e-6,
    )


def test_diffusion_lsm_time_metadata_matches_matlab_export():
    frap_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_FRAP_datasets.txt"
    _require_user_guide_fixtures(frap_export, DIFFUSION_ROOT / VENUS_CYTO_1)

    exports = _read_vector_export(frap_export)
    dataset = _load_venus_cytoplasm_1()

    matlab_time = exports[(VENUS_CYTO_1, "Time")]
    python_time = dataset.time - dataset.time[POST_BLEACH_INDEX]

    assert dataset.metadata.get("time_source") == "lsm_metadata"
    assert dataset.metadata.get("voxel_source") == "lsm_metadata"
    np.testing.assert_allclose(python_time, matlab_time, rtol=0, atol=1e-6)


def test_diffusion_corrected_mobile_fraction_stays_near_matlab_reference():
    params_path = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Fit_Parameters.txt"
    _require_user_guide_fixtures(params_path, DIFFUSION_ROOT / VENUS_CYTO_1)

    params = _read_parameter_table(params_path)
    dataset = _load_venus_cytoplasm_1()

    assert dataset.corrected_mobile_fraction == pytest.approx(
        params[VENUS_CYTO_1]["MF Corrected"],
        abs=1e-3,
    )


def test_diffusion_postbleach_profile_reproduces_matlab_export_for_venus_cytoplasm_1():
    profile_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Postbleach_profiles.txt"
    _require_user_guide_fixtures(profile_export, DIFFUSION_ROOT / VENUS_CYTO_1)

    profile_exports = _read_vector_export(profile_export)
    dataset = _load_venus_cytoplasm_1()

    np.testing.assert_allclose(
        dataset.radius,
        profile_exports[(VENUS_CYTO_1, "Distance")],
        rtol=0,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        dataset.post_bleach_profile,
        profile_exports[(VENUS_CYTO_1, "Post-bleach Profile")],
        rtol=0,
        atol=1e-6,
    )


def test_diffusion_loader_matches_exports_when_using_inferred_legacy_bleach_mask():
    frap_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_FRAP_datasets.txt"
    profile_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Postbleach_profiles.txt"
    params_path = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Fit_Parameters.txt"
    _require_user_guide_fixtures(frap_export, profile_export, params_path, DIFFUSION_ROOT / VENUS_CYTO_1)

    frap_exports = _read_vector_export(frap_export)
    profile_exports = _read_vector_export(profile_export)
    params = _read_parameter_table(params_path)
    dataset = _load_venus_cytoplasm_1_with_inferred_bleach_mask()
    python_time = dataset.time - dataset.time[POST_BLEACH_INDEX]
    inferred_mask = _inferred_legacy_venus_cytoplasm_1_bleach_mask((512, 512))

    assert inferred_mask.sum() == 252
    np.testing.assert_allclose(python_time, frap_exports[(VENUS_CYTO_1, "Time")], rtol=0, atol=1e-6)
    np.testing.assert_allclose(dataset.frap, frap_exports[(VENUS_CYTO_1, "Raw FRAP")], rtol=0, atol=1e-6)
    np.testing.assert_allclose(
        dataset.norm_frap,
        frap_exports[(VENUS_CYTO_1, "Normalized FRAP")],
        rtol=0,
        atol=1e-6,
    )
    np.testing.assert_allclose(dataset.radius, profile_exports[(VENUS_CYTO_1, "Distance")], rtol=0, atol=1e-6)
    np.testing.assert_allclose(
        dataset.post_bleach_profile,
        profile_exports[(VENUS_CYTO_1, "Post-bleach Profile")],
        rtol=0,
        atol=1e-6,
    )
    assert dataset.corrected_mobile_fraction == pytest.approx(
        params[VENUS_CYTO_1]["MF Corrected"],
        abs=1e-3,
    )


def test_diffusion_legacy_fit_reproduces_guide_venus_cytoplasm_1_from_matlab_exports():
    params_path = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Fit_Parameters.txt"
    frap_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_FRAP_datasets.txt"
    profile_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Postbleach_profiles.txt"
    _require_user_guide_fixtures(params_path, frap_export, profile_export)

    params = _read_parameter_table(params_path)
    frap_exports = _read_vector_export(frap_export)
    profile_exports = _read_vector_export(profile_export)
    datasets = _build_venus_cytoplasm_export_datasets(params, frap_exports, profile_exports)

    inputs = BasicInputs(
        file_paths=[Path("reference-export")],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=False,
        background_intensity=0.0,
        post_bleach_frame=POST_BLEACH_INDEX,
        roi_definition=(256.0, 23.0, 9.0),
        pre_bleach_frame_count=10,
        use_adjacent_roi=True,
    )
    config = _guide_diffusion_config(frap_exports, profile_exports)

    result = fit_diffusion_model(datasets, inputs, config)
    dataset = next(ds for ds in result.datasets if ds.name == VENUS_CYTO_1)
    expected = params[VENUS_CYTO_1]

    assert dataset.k == pytest.approx(expected["k"], abs=1e-4)
    assert dataset.r_effective == pytest.approx(expected["re"], abs=1e-3)
    assert dataset.diffusion_coefficient == pytest.approx(expected["D"], abs=1e-2)
    assert dataset.mobile_fraction == pytest.approx(expected["MF"], abs=1e-5)
    assert dataset.sum_squared_residuals == pytest.approx(expected["SS"], rel=1e-3)


def test_diffusion_guide_parameters_reproduce_exported_model_curves():
    params_path = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Fit_Parameters.txt"
    frap_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_FRAP_datasets.txt"
    profile_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Postbleach_profiles.txt"
    _require_user_guide_fixtures(params_path, frap_export, profile_export)

    params = _read_parameter_table(params_path)
    frap_exports = _read_vector_export(frap_export)
    profile_exports = _read_vector_export(profile_export)
    nominal_radius = 9.0 * 0.10986328426900892
    names = sorted(name for name in params if name.endswith(".lsm"))

    for name in names:
        expected = params[name]
        radius = profile_exports[(name, "Distance")]
        post_bleach_profile = profile_exports[(name, "Post-bleach Profile")]
        profile_fit = np.exp(
            -expected["k"] * np.exp(-2.0 * radius**2 / expected["re"] ** 2)
        )
        np.testing.assert_allclose(
            profile_fit,
            profile_exports[(name, "Profile Fit")],
            rtol=0,
            atol=2.5e-6,
        )
        np.testing.assert_allclose(
            post_bleach_profile - profile_fit,
            profile_exports[(name, "Fit Residuals")],
            rtol=0,
            atol=2.5e-6,
        )

        time = frap_exports[(name, "Time")]
        fit_time = time[POST_BLEACH_INDEX:] - time[POST_BLEACH_INDEX]
        corrected_frap = frap_exports[(name, "Corrected FRAP")][POST_BLEACH_INDEX:]
        weights = corrected_frap.sum()
        frap_initial = corrected_frap[0]
        frap_fit = (
            kang_frap(
                fit_time,
                expected["re"],
                nominal_radius,
                expected["D"],
                expected["k"],
            )
            * expected["MF"]
            + (1.0 - expected["MF"]) * frap_initial
        )
        weighted_residuals = corrected_frap / (fit_time + weights) - frap_fit / (
            fit_time + weights
        )

        np.testing.assert_allclose(
            fit_time,
            frap_exports[(name, "Fit Time")],
            rtol=0,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            frap_fit,
            frap_exports[(name, "FRAP Fit")],
            rtol=0,
            atol=2.5e-6,
        )
        np.testing.assert_allclose(
            weighted_residuals,
            frap_exports[(name, "Fit Residuals")],
            rtol=0,
            atol=7e-7,
        )
        assert np.dot(weighted_residuals, weighted_residuals) == pytest.approx(
            expected["SS"],
            rel=2e-3,
        )


def test_diffusion_d_mf_landscape_global_search_is_not_at_guide_d_for_venus_cytoplasm_3():
    params_path = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Fit_Parameters.txt"
    frap_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_FRAP_datasets.txt"
    _require_user_guide_fixtures(params_path, frap_export)

    params = _read_parameter_table(params_path)
    frap_exports = _read_vector_export(frap_export)
    name = "Venus_Cytoplasm_3.lsm"
    expected = params[name]
    time = frap_exports[(name, "Time")]
    corrected_frap = frap_exports[(name, "Corrected FRAP")]
    fit_time = time[POST_BLEACH_INDEX:] - time[POST_BLEACH_INDEX]
    fit_frap = corrected_frap[POST_BLEACH_INDEX:]
    weights = fit_frap.sum()
    target = fit_frap / (fit_time + weights)
    frap_initial = fit_frap[0]
    nominal_radius = 9.0 * 0.10986328426900892

    def residual_sum_squares(diffusion: float, mobile_fraction: float) -> float:
        prediction = (
            kang_frap(
                fit_time,
                expected["re"],
                nominal_radius,
                diffusion,
                expected["k"],
            )
            * mobile_fraction
            + (1.0 - mobile_fraction) * frap_initial
        ) / (fit_time + weights)
        residual = prediction - target
        return float(np.dot(residual, residual))

    guide_ss = residual_sum_squares(expected["D"], expected["MF"])

    # Profile the two-parameter objective over D by solving the linear least
    # squares problem for MF at each D. This is an exhaustive bounded check of
    # the D/MF basin, independent of the iterative optimizer path.
    diffusion_grid = np.linspace(20.0, 100.0, 801)
    profiled = []
    base = frap_initial / (fit_time + weights)
    centered_target = target - base
    for diffusion in diffusion_grid:
        slope = (
            kang_frap(
                fit_time,
                expected["re"],
                nominal_radius,
                diffusion,
                expected["k"],
            )
            - frap_initial
        ) / (fit_time + weights)
        mobile_fraction = np.dot(slope, centered_target) / np.dot(slope, slope)
        mobile_fraction = float(np.clip(mobile_fraction, 0.0, 2.0))
        profiled.append(
            (
                diffusion,
                mobile_fraction,
                residual_sum_squares(diffusion, mobile_fraction),
            )
        )

    best_diffusion, best_mobile_fraction, best_ss = min(profiled, key=lambda item: item[2])

    assert guide_ss == pytest.approx(expected["SS"], rel=1e-3)
    assert best_diffusion == pytest.approx(74.25, abs=0.2)
    assert best_mobile_fraction == pytest.approx(0.7873, abs=5e-4)
    assert best_ss < guide_ss * 0.75


def test_diffusion_legacy_optimizer_matches_guide_tables_from_matlab_exports():
    cases = [
        (
            "Venus_Cytoplasm",
            "Venus_Cytoplasm_1.lsm",
            0.50,
            1.5e-4,
            0.1,
        ),
        (
            "Venus-ATG5_Cytoplasm",
            "Venus-Atg5_Cytoplasm_1.lsm",
            0.01,
            1e-5,
            0.05,
        ),
    ]

    for prefix, first_name, max_d_error, max_mf_error, max_average_d_percent in cases:
        params_path = DIFFUSION_ROOT / f"{prefix}_Diffusion_Fit_Parameters.txt"
        frap_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_FRAP_datasets.txt"
        profile_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_Postbleach_profiles.txt"
        _require_user_guide_fixtures(params_path, frap_export, profile_export)

        params = _read_parameter_table(params_path)
        frap_exports = _read_vector_export(frap_export)
        profile_exports = _read_vector_export(profile_export)
        datasets = _build_diffusion_export_datasets(params, frap_exports, profile_exports)
        inputs = BasicInputs(
            file_paths=[Path("reference-export")],
            model_index=1,
            roi_mode=1,
            normalize_by_cell=False,
            background_intensity=0.0,
            post_bleach_frame=POST_BLEACH_INDEX,
            roi_definition=(256.0, 23.0, 9.0),
            pre_bleach_frame_count=10,
            use_adjacent_roi=True,
        )
        config = _guide_diffusion_config(frap_exports, profile_exports, first_name)

        result = fit_diffusion_model(datasets, inputs, config)
        guide_d = np.asarray([params[dataset.name]["D"] for dataset in result.datasets])
        fit_d = np.asarray([dataset.diffusion_coefficient for dataset in result.datasets])
        guide_mf = np.asarray([params[dataset.name]["MF"] for dataset in result.datasets])
        fit_mf = np.asarray([dataset.mobile_fraction for dataset in result.datasets])

        np.testing.assert_allclose(fit_d, guide_d, rtol=0, atol=max_d_error)
        np.testing.assert_allclose(fit_mf, guide_mf, rtol=0, atol=max_mf_error)
        assert abs((fit_d.mean() - guide_d.mean()) / guide_d.mean() * 100.0) <= max_average_d_percent


def test_diffusion_global_fit_concatenates_curves_instead_of_fitting_average_curve():
    cases = [
        (
            "Venus_Cytoplasm",
            "Venus_Cytoplasm_1.lsm",
            3.97188473,
            39.72106608,
            0.83394620,
            6.3471858e-05,
        ),
        (
            "Venus-ATG5_Cytoplasm",
            "Venus-Atg5_Cytoplasm_1.lsm",
            3.01699048,
            9.83250202,
            0.81492869,
            9.3612840e-05,
        ),
    ]

    for prefix, first_name, expected_re, expected_d, expected_mf, expected_concat_ss in cases:
        params_path = DIFFUSION_ROOT / f"{prefix}_Diffusion_Fit_Parameters.txt"
        frap_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_FRAP_datasets.txt"
        profile_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_Postbleach_profiles.txt"
        _require_user_guide_fixtures(params_path, frap_export, profile_export)

        params = _read_parameter_table(params_path)
        frap_exports = _read_vector_export(frap_export)
        profile_exports = _read_vector_export(profile_export)
        inputs = _reference_diffusion_inputs()

        global_datasets = _build_diffusion_export_datasets(params, frap_exports, profile_exports)
        global_config = _guide_diffusion_config(frap_exports, profile_exports, first_name)
        global_config.fit_mode = "global"
        global_config.optimizer_mode = "modern"
        global_config.fit_averaged_data = True
        global_result = fit_diffusion_model(global_datasets, inputs, global_config)

        average_datasets = _build_diffusion_export_datasets(params, frap_exports, profile_exports)
        average_config = _guide_diffusion_config(frap_exports, profile_exports, first_name)
        average_config.fit_mode = "average_curve"
        average_config.optimizer_mode = "modern"
        average_result = fit_diffusion_model(average_datasets, inputs, average_config)

        dataset_ss = sum(dataset.sum_squared_residuals for dataset in global_result.datasets)
        dataset_re = [dataset.r_effective for dataset in global_result.datasets]
        dataset_d = [dataset.diffusion_coefficient for dataset in global_result.datasets]
        dataset_mf = [dataset.mobile_fraction for dataset in global_result.datasets]

        assert global_result.r_effective == pytest.approx(expected_re, rel=1e-7)
        assert global_result.diffusion_coefficient == pytest.approx(expected_d, rel=1e-7)
        assert global_result.mobile_fraction == pytest.approx(expected_mf, rel=1e-7)
        assert global_result.sum_squared_residuals == pytest.approx(expected_concat_ss, rel=1e-7)
        assert global_result.sum_squared_residuals == pytest.approx(dataset_ss, rel=1e-12)
        np.testing.assert_allclose(dataset_re, global_result.r_effective, rtol=0, atol=1e-12)
        np.testing.assert_allclose(dataset_d, global_result.diffusion_coefficient, rtol=0, atol=1e-12)
        np.testing.assert_allclose(dataset_mf, global_result.mobile_fraction, rtol=0, atol=1e-12)
        assert global_result.sum_squared_residuals > average_result.sum_squared_residuals * 100.0


def test_diffusion_simplified_kang_uses_confocal_profile_equation():
    cases = [
        (
            "Venus_Cytoplasm",
            "Venus_Cytoplasm_1.lsm",
            4.9183,
            40.1459,
            0.08506,
        ),
        (
            "Venus-ATG5_Cytoplasm",
            "Venus-Atg5_Cytoplasm_1.lsm",
            3.9287,
            10.7484,
            0.20372,
        ),
    ]

    for prefix, first_name, expected_re, expected_d, expected_tau in cases:
        params_path = DIFFUSION_ROOT / f"{prefix}_Diffusion_Fit_Parameters.txt"
        frap_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_FRAP_datasets.txt"
        profile_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_Postbleach_profiles.txt"
        _require_user_guide_fixtures(params_path, frap_export, profile_export)

        params = _read_parameter_table(params_path)
        frap_exports = _read_vector_export(frap_export)
        profile_exports = _read_vector_export(profile_export)
        datasets = _build_diffusion_export_datasets(params, frap_exports, profile_exports)
        config = _guide_diffusion_config(frap_exports, profile_exports, first_name)
        config.fit_mode = "simplified_kang"
        config.optimizer_mode = "modern"

        result = fit_diffusion_model(datasets, _reference_diffusion_inputs(), config)

        assert result.r_effective == pytest.approx(expected_re, abs=5e-4)
        assert result.diffusion_coefficient == pytest.approx(expected_d, abs=5e-4)
        assert result.half_time == pytest.approx(expected_tau, abs=5e-5)


def test_diffusion_simplified_kang_global_pools_profile_and_recovery_curves():
    cases = [
        (
            "Venus_Cytoplasm",
            "Venus_Cytoplasm_1.lsm",
            4.90140590,
            40.51383439,
            0.83475502,
            0.07713860,
            6.3825569e-05,
        ),
        (
            "Venus-ATG5_Cytoplasm",
            "Venus-Atg5_Cytoplasm_1.lsm",
            3.92108297,
            10.02250597,
            0.81706244,
            0.20394796,
            9.3810997e-05,
        ),
    ]

    for prefix, first_name, expected_re, expected_d, expected_mf, expected_tau, expected_ss in cases:
        params_path = DIFFUSION_ROOT / f"{prefix}_Diffusion_Fit_Parameters.txt"
        frap_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_FRAP_datasets.txt"
        profile_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_Postbleach_profiles.txt"
        _require_user_guide_fixtures(params_path, frap_export, profile_export)

        params = _read_parameter_table(params_path)
        frap_exports = _read_vector_export(frap_export)
        profile_exports = _read_vector_export(profile_export)
        datasets = _build_diffusion_export_datasets(params, frap_exports, profile_exports)
        config = _guide_diffusion_config(frap_exports, profile_exports, first_name)
        config.fit_mode = "simplified_kang_global"
        config.optimizer_mode = "modern"
        config.fit_averaged_data = True

        result = fit_diffusion_model(datasets, _reference_diffusion_inputs(), config)
        dataset_ss = sum(dataset.sum_squared_residuals for dataset in result.datasets)

        assert result.r_effective == pytest.approx(expected_re, rel=1e-7)
        assert result.diffusion_coefficient == pytest.approx(expected_d, rel=1e-7)
        assert result.mobile_fraction == pytest.approx(expected_mf, rel=1e-7)
        assert result.half_time == pytest.approx(expected_tau, rel=1e-7)
        assert result.sum_squared_residuals == pytest.approx(expected_ss, rel=1e-7)
        assert result.sum_squared_residuals == pytest.approx(dataset_ss, rel=1e-12)
        np.testing.assert_allclose(
            [dataset.r_effective for dataset in result.datasets],
            result.r_effective,
            rtol=0,
            atol=1e-12,
        )


def test_diffusion_photodecay_stage_reproduces_corrected_frap_exports():
    cases = [
        ("Venus_Cytoplasm", VENUS_CYTO_1, 0.005127, 5e-6),
        ("Venus-ATG5_Cytoplasm", "Venus-Atg5_Cytoplasm_1.lsm", 0.005536, 8e-5),
    ]

    for prefix, first_name, expected_decay, corrected_atol in cases:
        params_path = DIFFUSION_ROOT / f"{prefix}_Diffusion_Fit_Parameters.txt"
        frap_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_FRAP_datasets.txt"
        profile_export = DIFFUSION_ROOT / f"{prefix}_Diffusion_Postbleach_profiles.txt"
        _require_user_guide_fixtures(params_path, frap_export, profile_export)

        params = _read_parameter_table(params_path)
        frap_exports = _read_vector_export(frap_export)
        profile_exports = _read_vector_export(profile_export)
        datasets = _build_diffusion_export_datasets(params, frap_exports, profile_exports)
        inputs = BasicInputs(
            file_paths=[Path("reference-export")],
            model_index=1,
            roi_mode=1,
            normalize_by_cell=False,
            background_intensity=0.0,
            post_bleach_frame=POST_BLEACH_INDEX,
            roi_definition=(256.0, 23.0, 9.0),
            pre_bleach_frame_count=10,
            use_adjacent_roi=True,
        )

        result = fit_diffusion_model(
            datasets,
            inputs,
            _guide_diffusion_config(frap_exports, profile_exports, first_name),
        )

        assert result.decay_rate == pytest.approx(expected_decay, rel=1e-4)
        for dataset in result.datasets:
            np.testing.assert_allclose(
                dataset.corrected_frap,
                frap_exports[(dataset.name, "Corrected FRAP")],
                rtol=0,
                atol=corrected_atol,
            )


@pytest.mark.xfail(
    reason=(
        "Strict MATLAB curve parity is the active preprocessing target; "
        "time/profile/photodecay-stage parity is in place, but circular ROI "
        "rasterization and normalized vector parity still differ from the "
        "2014 MATLAB exports."
    )
)
def test_diffusion_loader_strictly_reproduces_matlab_export_vectors_for_venus_cytoplasm_1():
    frap_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_FRAP_datasets.txt"
    profile_export = DIFFUSION_ROOT / "Venus_Cytoplasm_Diffusion_Postbleach_profiles.txt"
    _require_user_guide_fixtures(frap_export, profile_export, DIFFUSION_ROOT / VENUS_CYTO_1)

    frap_exports = _read_vector_export(frap_export)
    profile_exports = _read_vector_export(profile_export)
    dataset = _load_venus_cytoplasm_1()
    python_time = dataset.time - dataset.time[POST_BLEACH_INDEX]

    np.testing.assert_allclose(python_time, frap_exports[(VENUS_CYTO_1, "Time")], rtol=0, atol=1e-6)
    np.testing.assert_allclose(dataset.frap, frap_exports[(VENUS_CYTO_1, "Raw FRAP")], rtol=0, atol=1e-6)
    np.testing.assert_allclose(
        dataset.norm_frap,
        frap_exports[(VENUS_CYTO_1, "Normalized FRAP")],
        rtol=0,
        atol=1e-6,
    )
    np.testing.assert_allclose(dataset.radius, profile_exports[(VENUS_CYTO_1, "Distance")], rtol=0, atol=1e-6)
    np.testing.assert_allclose(
        dataset.post_bleach_profile,
        profile_exports[(VENUS_CYTO_1, "Post-bleach Profile")],
        rtol=0,
        atol=1e-6,
    )
