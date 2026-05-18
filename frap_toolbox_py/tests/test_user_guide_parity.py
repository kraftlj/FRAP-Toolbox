from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pytest

from frap_toolbox_py.data.loading import load_diffusion_datasets
from frap_toolbox_py.models.diffusion import fit_diffusion_model
from frap_toolbox_py.roi import CircularROI, adjacent_circle
from frap_toolbox_py.types import BasicInputs, DiffusionFitConfig, FitBounds, FRAPDataset


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


def _guide_diffusion_config(
    frap_exports: dict[tuple[str, str], np.ndarray],
    profile_exports: dict[tuple[str, str], np.ndarray],
) -> DiffusionFitConfig:
    return DiffusionFitConfig(
        k=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        r_effective=FitBounds(3.0, 0.0, np.inf, "Adjustable"),
        diffusion_coefficient=FitBounds(10.0, 0.0, np.inf, "Adjustable"),
        mobile_fraction=FitBounds(1.0, 0.0, 2.0, "Adjustable"),
        decay_rate=FitBounds(1e-3, 0.0, np.inf, "Adjustable"),
        profile_range=(0, len(profile_exports[(VENUS_CYTO_1, "Distance")])),
        frap_range=(POST_BLEACH_INDEX, len(frap_exports[(VENUS_CYTO_1, "Time")])),
        decay_fit_range=(499, 599),
        mobile_fraction_range=(539, 600),
        fit_averaged_data=False,
        optimizer_mode="legacy_matlab",
    )


def _build_venus_cytoplasm_export_datasets(
    params: dict[str, dict[str, float]],
    frap_exports: dict[tuple[str, str], np.ndarray],
    profile_exports: dict[tuple[str, str], np.ndarray],
) -> list[FRAPDataset]:
    datasets = []
    names = sorted(name for name in params if name.endswith(".lsm"))
    for name in names:
        time = frap_exports[(name, "Time")]
        datasets.append(
            FRAPDataset(
                name=name,
                time=time - time[0],
                frap=frap_exports[(name, "Raw FRAP")],
                norm_frap=frap_exports[(name, "Normalized FRAP")],
                corrected_frap=frap_exports[(name, "Normalized FRAP")].copy(),
                corrected_mobile_fraction=params[name]["MF Corrected"],
                voxel_size_x=0.10986328426900892,
                voxel_size_y=0.10986328426900892,
                radius=profile_exports[(name, "Distance")],
                post_bleach_profile=profile_exports[(name, "Post-bleach Profile")],
            )
        )
    return datasets


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


def test_diffusion_photodecay_stage_reproduces_corrected_frap_exports():
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

    result = fit_diffusion_model(datasets, inputs, _guide_diffusion_config(frap_exports, profile_exports))

    assert result.decay_rate == pytest.approx(0.005127, rel=1e-4)
    for dataset in result.datasets:
        np.testing.assert_allclose(
            dataset.corrected_frap,
            frap_exports[(dataset.name, "Corrected FRAP")],
            rtol=0,
            atol=5e-6,
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
