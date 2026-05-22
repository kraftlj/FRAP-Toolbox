from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Iterable, List, Optional, Tuple, Union

import numpy as np
from ome_types import from_xml

try:
    from bioio import BioImage
except ImportError:  # pragma: no cover - exercised in legacy environments
    BioImage = None

try:
    from bioio_tifffile import Reader as BioIOTiffReader
except ImportError:  # pragma: no cover - plugin is optional at runtime
    BioIOTiffReader = None

try:
    from aicsimageio import AICSImage
    from aicsimageio.readers import TiffReader as AICSTiffReader
except ImportError:  # pragma: no cover - legacy fallback may not be installed
    AICSImage = None
    AICSTiffReader = None

try:
    import tifffile
except ImportError:  # pragma: no cover - optional dependency guards runtime envs
    tifffile = None

from ..initial_conditions import radial_profile
from ..normalization import normalize_frap
from ..types import BasicInputs, FRAPDataset

MaskFactory = Callable[[Tuple[int, int]], np.ndarray]
LOGGER = logging.getLogger(__name__)


def _make_mask(
    factory: Optional[Union[MaskFactory, np.ndarray]],
    shape: Tuple[int, int],
) -> Optional[np.ndarray]:
    if factory is None:
        return None
    if isinstance(factory, np.ndarray):
        if factory.shape != shape:
            raise ValueError("Provided mask array does not match image shape.")
        return factory.astype(bool)
    mask = factory(shape)
    if mask.shape != shape:
        raise ValueError("Mask factory returned mask with incorrect shape.")
    return mask.astype(bool)


def _open_image(path: Path) -> Tuple[Any, str]:
    errors: list[str] = []

    if BioImage is not None:
        try:
            return BioImage(path), "bioio"
        except Exception as exc:
            errors.append(f"BioIO default reader failed: {exc}")

        if BioIOTiffReader is not None:
            try:
                return BioImage(path, reader=BioIOTiffReader), "bioio_tifffile"
            except Exception as exc:
                errors.append(f"BioIO tifffile reader failed: {exc}")

    if AICSImage is not None:
        if AICSTiffReader is not None:
            try:
                return AICSImage(path, reader=AICSTiffReader), "aicsimageio_tifffile"
            except Exception as exc:
                errors.append(f"AICSImageIO tifffile reader failed: {exc}")

        try:
            return AICSImage(path), "aicsimageio"
        except Exception as exc:
            errors.append(f"AICSImageIO default reader failed: {exc}")

    details = "\n".join(errors) if errors else "No BioIO or legacy AICSImageIO reader is installed."
    raise ImportError(
        "Could not open the microscopy image. Install the BioIO reader plugin "
        "for this format, or install the legacy-io extra for AICSImageIO fallback.\n"
        f"{details}"
    )


def _seconds(value: Any, unit: Any = None) -> float:
    if hasattr(value, "total_seconds"):
        return float(value.total_seconds())

    numeric = float(value)
    if unit is None:
        return numeric

    unit_name = str(getattr(unit, "name", "") or "").lower()
    unit_value = str(getattr(unit, "value", unit) or "").lower()
    unit_text = f"{unit_name} {unit_value}"
    if "nanosecond" in unit_text or unit_value == "ns":
        return numeric * 1e-9
    if "microsecond" in unit_text or unit_value == "us":
        return numeric * 1e-6
    if "millisecond" in unit_text or unit_value == "ms":
        return numeric * 1e-3
    if "minute" in unit_text or unit_value == "min":
        return numeric * 60.0
    if "hour" in unit_text or unit_value == "h":
        return numeric * 3600.0
    return numeric


def _image_time_length(image: Any, frame_count: Optional[int]) -> int:
    if frame_count is not None:
        return frame_count
    return image.dims.T or 1


def _extract_time_vector(image: Any, frame_count: Optional[int] = None) -> Tuple[np.ndarray, str]:
    try:
        metadata = image.metadata
    except Exception as exc:
        LOGGER.debug("Could not parse image metadata for timestamps: %s", exc)
        return np.arange(_image_time_length(image, frame_count), dtype=float), "synthetic_linear"

    ome = None
    if isinstance(metadata, str):
        try:
            ome = from_xml(metadata)
        except Exception:
            ome = None
    elif hasattr(metadata, "images"):
        ome = metadata

    if ome is None:
        return np.arange(_image_time_length(image, frame_count), dtype=float), "synthetic_linear"

    pixels = ome.images[0].pixels
    planes = pixels.planes
    times = []
    by_t = {}
    for plane in planes:
        if plane.delta_t is None:
            continue
        key = plane.the_c, plane.the_z, plane.the_t
        if key[2] not in by_t:
            by_t[key[2]] = (plane.delta_t, getattr(plane, "delta_t_unit", None))
    if by_t:
        for t_index in range(_image_time_length(image, frame_count)):
            entry = by_t.get(t_index)
            if entry is None:
                times.append(times[-1] if times else 0.0)
            else:
                delta, unit = entry
                times.append(_seconds(delta, unit))
        return np.asarray(times, dtype=float), "ome_planes"

    exposure_time = pixels.time_increment
    if exposure_time is not None:
        count = _image_time_length(image, frame_count)
        return (
            np.linspace(
                0.0,
                _seconds(
                    exposure_time,
                    getattr(pixels, "time_increment_unit", None),
                )
                * (count - 1),
                count,
            ),
            "ome_time_increment",
        )

    return np.arange(_image_time_length(image, frame_count), dtype=float), "synthetic_linear"


def _ensure_microns(value: Optional[Union[float, int]]) -> Optional[float]:
    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(numeric) or numeric <= 0:
        return None
    if numeric < 1e-5:
        return numeric * 1e6
    return numeric


def _sanitize_lsm_times(raw_times: Iterable[float], frame_count: int) -> Optional[np.ndarray]:
    times = np.asarray(list(raw_times), dtype=float)
    if times.size == 0:
        return None
    times = times[:frame_count]
    if times.size < frame_count:
        return None
    times = times - times[0]
    if not np.all(np.isfinite(times)):
        return None
    if times.size > 1:
        diffs = np.diff(times)
        if np.nanmin(diffs) < 0:
            times = np.maximum.accumulate(times)
    return times


def _read_lsm_metadata(path: Path, frame_count: int) -> Tuple[Optional[float], Optional[float], Optional[np.ndarray]]:
    if tifffile is None:
        return None, None, None
    try:
        with tifffile.TiffFile(str(path)) as tf:
            if not getattr(tf, "is_lsm", False):
                return None, None, None
            metadata = tf.lsm_metadata
    except Exception as exc:  # pragma: no cover - handles corrupt/unexpected files
        LOGGER.debug("Failed to read LSM metadata from %s: %s", path, exc)
        return None, None, None

    if metadata is None:
        return None, None, None

    def _lsm_lookup(obj, *keys):
        for key in keys:
            if hasattr(obj, key):
                value = getattr(obj, key)
                if value is not None:
                    return value
            if isinstance(obj, dict) and key in obj:
                value = obj[key]
                if value is not None:
                    return value
        return None

    voxel_x = _ensure_microns(_lsm_lookup(metadata, "voxel_size_x", "VoxelSizeX"))
    voxel_y = _ensure_microns(_lsm_lookup(metadata, "voxel_size_y", "VoxelSizeY"))

    raw_times = _lsm_lookup(metadata, "time_stamps", "TimeStamps")
    times = None
    if raw_times is not None:
        times = _sanitize_lsm_times(raw_times, frame_count)

    return voxel_x, voxel_y, times


def load_diffusion_datasets(
    inputs: BasicInputs,
    bleach_roi: Union[MaskFactory, np.ndarray],
    cell_roi: Optional[Union[MaskFactory, np.ndarray]] = None,
    adjacent_roi: Optional[Union[MaskFactory, np.ndarray]] = None,
) -> List[FRAPDataset]:
    """Load FRAP datasets for the diffusion model."""

    datasets: List[FRAPDataset] = []
    for file_path in inputs.file_paths:
        path = Path(file_path)
        image, reader_backend = _open_image(path)
        LOGGER.debug("Loaded %s with %s.", path, reader_backend)
        raw_stack = image.get_image_data("TYX", C=0, Z=0)
        stack = np.asarray(raw_stack, dtype=float)
        if stack.ndim != 3:
            raise ValueError("Expected stack with dimensions (T, Y, X).")

        frame_count = stack.shape[0]
        lsm_voxel_x, lsm_voxel_y, lsm_times = _read_lsm_metadata(path, frame_count)

        times, time_source = _extract_time_vector(image, frame_count)
        times = np.asarray(times, dtype=float)

        if lsm_times is not None:
            times = lsm_times
            time_source = "lsm_metadata"

        if len(times) != frame_count:
            if len(times) == 1:
                times = np.linspace(0, times[0], frame_count)
                time_source = f"{time_source}_expanded"
            else:
                raise ValueError("Mismatch between time metadata and stack length.")

        image_shape = stack.shape[1:]
        bleach_mask = _make_mask(bleach_roi, image_shape)
        if bleach_mask is None:
            raise ValueError("Bleach ROI must be provided for diffusion loading.")

        cell_mask = _make_mask(cell_roi, image_shape)
        adjacent_mask = _make_mask(adjacent_roi, image_shape)

        frap_series = stack[:, bleach_mask].mean(axis=1)
        frap_series = frap_series - inputs.background_intensity

        cell_series: Optional[np.ndarray] = None
        if cell_mask is not None and inputs.normalize_by_cell:
            cell_series = stack[:, cell_mask].mean(axis=1)

        adjacent_series: Optional[np.ndarray] = None
        if adjacent_mask is not None and inputs.use_adjacent_roi:
            adjacent_series = stack[:, adjacent_mask].mean(axis=1) - inputs.background_intensity

        norm_frap, corrected_mf, norm_adjacent = normalize_frap(
            frap_series,
            cell_series,
            inputs.post_bleach_frame,
            inputs.pre_bleach_frame_count,
            inputs.normalize_by_cell,
            adjacent_series,
        )

        center_x: float
        center_y: float
        if inputs.roi_mode == 1:
            # The original MATLAB post-bleach profile uses a 1-based meshgrid
            # and, for the guide's Zeiss LSM examples, lands halfway between
            # image rows. Keeping this offset reproduces the exported Distance
            # and Post-bleach Profile tables.
            center_x = float(inputs.roi_definition[0])
            center_y = float(inputs.roi_definition[1]) - 0.5
        else:
            coords = np.argwhere(bleach_mask)
            if coords.size == 0:
                raise ValueError("Bleach ROI mask is empty; cannot determine centre.")
            center_y_zero, center_x_zero = coords.mean(axis=0)
            center_x = center_x_zero + 1.0
            center_y = center_y_zero + 1.0

        pre_start = inputs.post_bleach_frame - inputs.pre_bleach_frame_count
        pre_end = inputs.post_bleach_frame
        if pre_start < 0:
            raise ValueError("Pre-bleach frame range exceeds available frames.")
        profile_pre_start = pre_start - 1
        profile_pre_end = pre_end - 1
        if profile_pre_start < 0:
            raise ValueError("Post-bleach profile pre-bleach range exceeds available frames.")
        pre_bleach_stack = stack[profile_pre_start:profile_pre_end]
        post_bleach_image = stack[inputs.post_bleach_frame]

        raw_pre_stack = np.asarray(raw_stack[profile_pre_start:profile_pre_end])
        raw_post_image = np.asarray(raw_stack[inputs.post_bleach_frame])
        mean_pre_bleach_raw = raw_pre_stack.mean(axis=0)

        invalid_mask = np.zeros_like(post_bleach_image, dtype=bool)
        if np.issubdtype(raw_stack.dtype, np.integer):
            max_val = np.iinfo(raw_stack.dtype).max
            invalid_mask |= raw_post_image == max_val
            invalid_mask |= mean_pre_bleach_raw == max_val
        invalid_mask |= raw_post_image == 0
        invalid_mask |= mean_pre_bleach_raw == 0

        radial_r, radial_profile_values = radial_profile(
            post_bleach_image,
            pre_bleach_stack,
            center=(center_x, center_y),
            invalid_mask=invalid_mask,
        )

        pixel_sizes = image.physical_pixel_sizes
        ome_voxel_x = _ensure_microns(float(pixel_sizes.X)) if pixel_sizes.X is not None else None
        ome_voxel_y = _ensure_microns(float(pixel_sizes.Y)) if pixel_sizes.Y is not None else None

        voxel_x = lsm_voxel_x if lsm_voxel_x is not None else ome_voxel_x
        voxel_y = lsm_voxel_y if lsm_voxel_y is not None else ome_voxel_y

        if voxel_x is not None:
            radial_r = radial_r * voxel_x

        if lsm_voxel_x is not None or lsm_voxel_y is not None:
            voxel_source = "lsm_metadata"
        elif ome_voxel_x is not None or ome_voxel_y is not None:
            voxel_source = "ome_metadata"
        else:
            voxel_source = "unspecified"

        dataset = FRAPDataset(
            name=path.name,
            time=times,
            frap=frap_series,
            norm_frap=norm_frap,
            corrected_frap=norm_frap.copy(),
            cell=cell_series,
            adjacent=norm_adjacent,
            corrected_mobile_fraction=corrected_mf,
            voxel_size_x=voxel_x,
            voxel_size_y=voxel_y,
            radius=radial_r,
            post_bleach_profile=radial_profile_values,
            metadata={
                "source_path": str(path),
                "reader_backend": reader_backend,
                "time_source": time_source,
                "voxel_source": voxel_source,
                "image_shape": tuple(int(value) for value in image_shape),
            },
        )
        datasets.append(dataset)

    return datasets


def load_reaction_datasets(
    inputs: BasicInputs,
    bleach_roi: Union[MaskFactory, np.ndarray],
    cell_roi: Optional[Union[MaskFactory, np.ndarray]] = None,
) -> List[FRAPDataset]:
    """Load FRAP datasets for the reaction models."""

    datasets: List[FRAPDataset] = []
    for file_path in inputs.file_paths:
        path = Path(file_path)
        image, reader_backend = _open_image(path)
        LOGGER.debug("Loaded %s with %s.", path, reader_backend)
        raw_stack = image.get_image_data("TYX", C=0, Z=0)
        stack = np.asarray(raw_stack, dtype=float)
        if stack.ndim != 3:
            raise ValueError("Expected stack with dimensions (T, Y, X).")

        frame_count = stack.shape[0]
        lsm_voxel_x, lsm_voxel_y, lsm_times = _read_lsm_metadata(path, frame_count)
        times, time_source = _extract_time_vector(image, frame_count)
        times = np.asarray(times, dtype=float)

        if lsm_times is not None:
            times = lsm_times
            time_source = "lsm_metadata"

        if len(times) != frame_count:
            if len(times) == 1:
                times = np.linspace(0, times[0], frame_count)
                time_source = f"{time_source}_expanded"
            else:
                raise ValueError("Mismatch between time metadata and stack length.")

        image_shape = stack.shape[1:]
        bleach_mask = _make_mask(bleach_roi, image_shape)
        if bleach_mask is None:
            raise ValueError("Bleach ROI must be provided for reaction loading.")

        cell_mask = _make_mask(cell_roi, image_shape)
        frap_series = stack[:, bleach_mask].mean(axis=1) - inputs.background_intensity

        cell_series: Optional[np.ndarray] = None
        if inputs.normalize_by_cell:
            if cell_mask is None:
                raise ValueError("Cell ROI must be provided when whole-cell normalization is enabled.")
            cell_series = stack[:, cell_mask].mean(axis=1)

        norm_frap, _, _ = normalize_frap(
            frap_series,
            cell_series,
            inputs.post_bleach_frame,
            inputs.pre_bleach_frame_count,
            inputs.normalize_by_cell,
        )

        pixel_sizes = image.physical_pixel_sizes
        ome_voxel_x = _ensure_microns(float(pixel_sizes.X)) if pixel_sizes.X is not None else None
        ome_voxel_y = _ensure_microns(float(pixel_sizes.Y)) if pixel_sizes.Y is not None else None
        voxel_x = lsm_voxel_x if lsm_voxel_x is not None else ome_voxel_x
        voxel_y = lsm_voxel_y if lsm_voxel_y is not None else ome_voxel_y

        if lsm_voxel_x is not None or lsm_voxel_y is not None:
            voxel_source = "lsm_metadata"
        elif ome_voxel_x is not None or ome_voxel_y is not None:
            voxel_source = "ome_metadata"
        else:
            voxel_source = "unspecified"

        datasets.append(
            FRAPDataset(
                name=path.name,
                time=times,
                frap=frap_series,
                norm_frap=norm_frap,
                corrected_frap=norm_frap.copy(),
                cell=cell_series,
                voxel_size_x=voxel_x,
                voxel_size_y=voxel_y,
                metadata={
                    "source_path": str(path),
                    "reader_backend": reader_backend,
                    "time_source": time_source,
                    "voxel_source": voxel_source,
                    "image_shape": tuple(int(value) for value in image_shape),
                },
            )
        )

    return datasets
