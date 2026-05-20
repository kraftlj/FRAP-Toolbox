from __future__ import annotations

import argparse
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import tempfile
from typing import Optional

import numpy as np
import pandas as pd

from frap_toolbox_py.app import (
    _build_diffusion_result,
    _build_reaction_result,
    _default_fit_mode,
    _plot_frap,
    _plot_profile,
    _result_table,
)
from frap_toolbox_py.cloud.manifest import OutputManifest, RoiGeometry, RunManifest, RunStatus, utc_now
from frap_toolbox_py.cloud.storage import GCSStorageBackend, LocalStorageBackend, StorageBackend
from frap_toolbox_py.roi import CircularROI, PolygonROI, load_roi_mask, save_roi_masks


def package_version() -> str:
    try:
        return version("frap-toolbox-py")
    except PackageNotFoundError:  # pragma: no cover - editable edge case
        return "0.0.0+local"


def join_uri(prefix: str, suffix: str) -> str:
    return f"{prefix.rstrip('/')}/{suffix.lstrip('/')}"


def storage_for_uri(uri: str, local_root: Optional[Path] = None) -> StorageBackend:
    if uri.startswith("gs://"):
        bucket = uri.removeprefix("gs://").split("/", 1)[0]
        return GCSStorageBackend(bucket)
    return LocalStorageBackend(local_root or Path(".frap-cloud-storage"))


def _load_manifest(storage: StorageBackend, uri: str) -> RunManifest:
    return RunManifest.from_dict(storage.read_json(uri))


def _materialize_inputs(manifest: RunManifest, storage: StorageBackend, work_dir: Path) -> list[Path]:
    input_dir = work_dir / "inputs"
    files: list[Path] = []
    for uploaded_file in manifest.files:
        destination = input_dir / uploaded_file.name
        files.append(storage.materialize(uploaded_file.object_uri, destination))
    return files


def _mask_factory_from_roi(roi: RoiGeometry, storage: StorageBackend, work_dir: Path):
    if roi.kind == "circle":
        circle = CircularROI(
            center_x=float(roi.center_x),
            center_y=float(roi.center_y),
            radius=float(roi.radius),
        )
        return circle.to_mask, 1, (circle.center_x, circle.center_y, circle.radius)

    if roi.kind in {"polygon", "freehand"}:
        polygon = PolygonROI(roi.points)
        return polygon.to_mask, 2, ()

    if roi.kind == "mask":
        if not roi.mask_uri:
            raise ValueError(f"{roi.role} mask ROI is missing mask_uri.")
        local_mask = storage.materialize(roi.mask_uri, work_dir / "rois" / Path(roi.mask_uri).name)
        mask = load_roi_mask(local_mask, roi.role, allow_single=True)
        return mask, 2, ()

    raise ValueError(f"Unsupported ROI kind {roi.kind!r}.")


def _roi_for_role(manifest: RunManifest, role: str) -> Optional[RoiGeometry]:
    for roi in manifest.rois:
        if roi.role == role:
            return roi
    return None


def _infer_image_shape(path: Path) -> tuple[int, int]:
    from frap_toolbox_py.data.loading import _open_image

    image, _ = _open_image(path)
    stack = np.asarray(image.get_image_data("TYX", C=0, Z=0))
    if stack.ndim != 3:
        raise ValueError("Expected stack with dimensions (T, Y, X).")
    return tuple(int(value) for value in stack.shape[1:])  # type: ignore[return-value]


def _save_generated_roi_masks(
    manifest: RunManifest,
    storage: StorageBackend,
    output_prefix: str,
    work_dir: Path,
    first_input: Path,
) -> dict[str, str]:
    shape = _infer_image_shape(first_input)
    masks: dict[str, np.ndarray] = {}
    for roi in manifest.rois:
        if roi.kind == "mask":
            if roi.mask_uri:
                masks[roi.role] = load_roi_mask(storage.materialize(roi.mask_uri, work_dir / "rois" / f"{roi.role}.npz"), roi.role, allow_single=True)
            continue
        factory, _, _ = _mask_factory_from_roi(roi, storage, work_dir)
        masks[roi.role] = factory(shape)

    if not masks:
        return {}

    local_path = work_dir / "outputs" / "roi-masks.npz"
    save_roi_masks(
        local_path,
        masks,
        image_shape=shape,
        roi_kind="cloud-browser",
        source_file=manifest.files[0].name,
        extra_metadata={"run_id": manifest.run_id},
    )
    output_uri = join_uri(output_prefix, "roi-masks.npz")
    storage.write_bytes(output_uri, local_path.read_bytes(), "application/octet-stream")
    return {name: output_uri for name in masks}


def _result_parameters(table: pd.DataFrame) -> dict[str, float | None]:
    parameters: dict[str, float | None] = {}
    for row in table.to_dict(orient="records"):
        value = row["Value"]
        parameters[row["Parameter"]] = None if pd.isna(value) else float(value)
    return parameters


def run_analysis(
    manifest_uri: str,
    *,
    output_prefix: Optional[str] = None,
    storage: Optional[StorageBackend] = None,
    local_storage_root: Optional[Path] = None,
    dry_run: bool = False,
) -> OutputManifest:
    storage = storage or storage_for_uri(manifest_uri, local_storage_root)
    manifest = _load_manifest(storage, manifest_uri)
    resolved_output_prefix = output_prefix or manifest.output_prefix

    with tempfile.TemporaryDirectory(prefix=f"frap-{manifest.run_id}-") as temp:
        work_dir = Path(temp)

        if dry_run:
            output = OutputManifest(
                run_id=manifest.run_id,
                status=RunStatus.SUCCEEDED,
                warnings=["Dry run completed without reading microscopy inputs."],
                container_version=package_version(),
                finished_at=utc_now(),
            )
            storage.write_json(join_uri(resolved_output_prefix, "output-manifest.json"), output.to_dict())
            return output

        input_files = _materialize_inputs(manifest, storage, work_dir)
        bleach_roi = _roi_for_role(manifest, "bleach")
        if bleach_roi is None:
            raise ValueError("Manifest is missing a bleach ROI.")
        bleach_factory, roi_mode, roi_definition = _mask_factory_from_roi(bleach_roi, storage, work_dir)

        cell_factory = None
        cell_roi = _roi_for_role(manifest, "cell")
        if cell_roi is not None:
            cell_factory, _, _ = _mask_factory_from_roi(cell_roi, storage, work_dir)

        settings = manifest.settings
        fit_mode = settings.fit_mode or _default_fit_mode(settings.model)
        if settings.model == "diffusion":
            result = _build_diffusion_result(
                input_files,
                bleach_factory,
                roi_mode,
                roi_definition,
                settings.post_bleach_frame,
                settings.pre_bleach_count,
                settings.background,
                settings.normalize_by_cell,
                cell_factory,
                settings.use_adjacent_roi,
                settings.adjacent_offset,
                settings.max_profile_radius,
                fit_mode,
            )
        else:
            result = _build_reaction_result(
                input_files,
                settings.model,
                bleach_factory,
                roi_mode,
                roi_definition,
                settings.post_bleach_frame,
                settings.pre_bleach_count,
                settings.background,
                settings.normalize_by_cell,
                cell_factory,
                fit_mode,
            )

        output_dir = work_dir / "outputs"
        output_dir.mkdir(parents=True, exist_ok=True)
        result_table = _result_table(result)
        result_table_uri = join_uri(resolved_output_prefix, "result-parameters.csv")
        storage.write_bytes(
            result_table_uri,
            result_table.to_csv(index=False).encode("utf-8"),
            "text/csv",
        )

        frap_plot = output_dir / "frap-fit.png"
        _plot_frap(result).savefig(frap_plot, dpi=160)
        frap_plot_uri = join_uri(resolved_output_prefix, "frap-fit.png")
        storage.write_bytes(frap_plot_uri, frap_plot.read_bytes(), "image/png")

        plots = {"frap_fit": frap_plot_uri}
        if getattr(result, "averaged_profile", None) is not None:
            profile_plot = output_dir / "post-bleach-profile.png"
            _plot_profile(result).savefig(profile_plot, dpi=160)
            profile_plot_uri = join_uri(resolved_output_prefix, "post-bleach-profile.png")
            storage.write_bytes(profile_plot_uri, profile_plot.read_bytes(), "image/png")
            plots["post_bleach_profile"] = profile_plot_uri

        roi_masks = _save_generated_roi_masks(manifest, storage, resolved_output_prefix, work_dir, input_files[0])
        output = OutputManifest(
            run_id=manifest.run_id,
            status=RunStatus.SUCCEEDED,
            parameters=_result_parameters(result_table),
            residual_summaries={"sum_squared_residuals": float(getattr(result, "sum_squared_residuals", 0.0))},
            result_tables={"parameters": result_table_uri},
            plots=plots,
            roi_masks=roi_masks,
            container_version=package_version(),
            input_hashes={uploaded_file.name: uploaded_file.checksum for uploaded_file in manifest.files if uploaded_file.checksum},
            finished_at=utc_now(),
        )
        storage.write_json(join_uri(resolved_output_prefix, "output-manifest.json"), output.to_dict())
        return output


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a cloud FRAP analysis manifest.")
    parser.add_argument("--manifest", required=True, help="Run manifest URI, usually gs://.../manifest.json")
    parser.add_argument("--output-prefix", help="Override the manifest output prefix.")
    parser.add_argument("--local-storage-root", type=Path, help="Root for local:// storage during tests/dev.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and write an output manifest without analysis.")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    try:
        run_analysis(
            args.manifest,
            output_prefix=args.output_prefix,
            local_storage_root=args.local_storage_root,
            dry_run=args.dry_run,
        )
    except Exception as exc:
        raise SystemExit(f"Cloud FRAP worker failed: {exc}") from exc


if __name__ == "__main__":
    main()
