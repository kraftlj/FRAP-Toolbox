from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Tuple

import numpy as np
from skimage.draw import disk, polygon

ROI_MASK_FORMAT = "frap-toolbox-roi-masks"
ROI_MASK_VERSION = 1
_MASK_KEY_PREFIX = "mask__"
_MASK_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")


@dataclass
class CircularROI:
    """Simple representation of a circular ROI in pixel coordinates."""

    center_x: float
    center_y: float
    radius: float

    def to_mask(self, shape: Tuple[int, int]) -> np.ndarray:
        # MATLAB specifies ROI centres in 1-based image coordinates where pixel
        # centres fall on integer locations. Convert to NumPy's 0-based index
        # space by subtracting one before rasterising the disk so the binary
        # mask matches the original toolbox geometry.
        row = self.center_y - 1.0
        col = self.center_x - 1.0
        rr, cc = disk((row, col), self.radius, shape=shape)
        mask = np.zeros(shape, dtype=bool)
        mask[rr, cc] = True
        return mask


@dataclass
class PolygonROI:
    """Representation of a user-defined polygon ROI."""

    vertices: Iterable[Tuple[float, float]]

    def to_mask(self, shape: Tuple[int, int]) -> np.ndarray:
        vertices = np.asarray(self.vertices, dtype=float)
        rr, cc = polygon(vertices[:, 1], vertices[:, 0], shape=shape)
        mask = np.zeros(shape, dtype=bool)
        mask[rr, cc] = True
        return mask


@dataclass(frozen=True)
class ROIMaskSet:
    """Saved ROI masks and JSON-compatible metadata."""

    masks: Mapping[str, np.ndarray]
    metadata: dict[str, Any]

    @property
    def image_shape(self) -> Tuple[int, int]:
        return _normalize_image_shape(self.metadata["image_shape"])

    def require_mask(
        self,
        name: str,
        *,
        expected_shape: Optional[Tuple[int, int]] = None,
        allow_single: bool = False,
    ) -> np.ndarray:
        """Return a named mask, optionally accepting a single-mask file."""

        if name in self.masks:
            mask = self.masks[name]
        elif allow_single and len(self.masks) == 1:
            mask = next(iter(self.masks.values()))
        else:
            available = ", ".join(sorted(self.masks)) or "<none>"
            raise ValueError(
                f"ROI mask file does not contain mask {name!r}. "
                f"Available masks: {available}."
            )

        if expected_shape is not None and tuple(mask.shape) != tuple(expected_shape):
            expected_shape_tuple = tuple(expected_shape)
            raise ValueError(
                f"ROI mask {name!r} has shape {tuple(mask.shape)}, "
                f"expected {expected_shape_tuple}."
            )
        return mask.copy()


def adjacent_circle(base_circle: CircularROI, offset_factor: float = 2.5) -> CircularROI:
    """Return an adjacent circular ROI offset along the x-axis."""

    return CircularROI(
        center_x=base_circle.center_x + base_circle.radius * offset_factor,
        center_y=base_circle.center_y,
        radius=base_circle.radius,
    )


def _normalize_image_shape(image_shape: Iterable[int]) -> Tuple[int, int]:
    shape = tuple(int(value) for value in image_shape)
    if len(shape) != 2 or any(value <= 0 for value in shape):
        raise ValueError("ROI mask image_shape must be a positive (Y, X) pair.")
    return shape  # type: ignore[return-value]


def _validate_mask_name(name: str) -> str:
    if not isinstance(name, str) or not name:
        raise ValueError("ROI mask names must be non-empty strings.")
    if not _MASK_NAME_PATTERN.match(name):
        raise ValueError(
            "ROI mask names may contain only letters, numbers, underscores, and hyphens."
        )
    return name


def _validate_mask(name: str, mask: np.ndarray, image_shape: Tuple[int, int]) -> np.ndarray:
    if not isinstance(mask, np.ndarray):
        raise ValueError(f"ROI mask {name!r} must be a NumPy array.")
    if mask.dtype != np.bool_:
        raise ValueError(f"ROI mask {name!r} must have dtype bool.")
    if tuple(mask.shape) != image_shape:
        raise ValueError(
            f"ROI mask {name!r} has shape {tuple(mask.shape)}, expected {image_shape}."
        )
    return mask.copy()


def build_roi_mask_metadata(
    masks: Mapping[str, np.ndarray],
    *,
    image_shape: Optional[Iterable[int]] = None,
    roi_kind: str = "manual",
    source_file: Optional[str] = None,
    notes: Optional[str] = None,
    extra_metadata: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Build JSON-compatible metadata for a saved ROI mask container."""

    if not masks:
        raise ValueError("At least one ROI mask must be provided.")

    if image_shape is None:
        first_mask = next(iter(masks.values()))
        image_shape_tuple = _normalize_image_shape(first_mask.shape)
    else:
        image_shape_tuple = _normalize_image_shape(image_shape)

    mask_entries = []
    for name, mask in masks.items():
        validated_name = _validate_mask_name(name)
        validated_mask = _validate_mask(validated_name, mask, image_shape_tuple)
        mask_entries.append(
            {
                "name": validated_name,
                "shape": list(validated_mask.shape),
                "pixel_count": int(validated_mask.sum()),
            }
        )

    metadata: dict[str, Any] = {
        "format": ROI_MASK_FORMAT,
        "version": ROI_MASK_VERSION,
        "image_shape": list(image_shape_tuple),
        "roi_kind": roi_kind,
        "source_file": source_file,
        "notes": notes,
        "masks": mask_entries,
    }
    if extra_metadata:
        metadata["extra"] = dict(extra_metadata)
    return metadata


def save_roi_masks(
    path: Path | str,
    masks: Mapping[str, np.ndarray],
    *,
    image_shape: Optional[Iterable[int]] = None,
    roi_kind: str = "manual",
    source_file: Optional[str] = None,
    notes: Optional[str] = None,
    extra_metadata: Optional[Mapping[str, Any]] = None,
) -> None:
    """Save named boolean ROI masks in the FRAP-Toolbox ``.npz`` container."""

    metadata = build_roi_mask_metadata(
        masks,
        image_shape=image_shape,
        roi_kind=roi_kind,
        source_file=source_file,
        notes=notes,
        extra_metadata=extra_metadata,
    )
    normalized_shape = _normalize_image_shape(metadata["image_shape"])

    arrays: dict[str, np.ndarray] = {
        "metadata_json": np.asarray(json.dumps(metadata, sort_keys=True)),
    }
    for name, mask in masks.items():
        validated_name = _validate_mask_name(name)
        arrays[f"{_MASK_KEY_PREFIX}{validated_name}"] = _validate_mask(
            validated_name,
            mask,
            normalized_shape,
        )

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        np.savez_compressed(handle, **arrays)


def save_roi_mask(
    path: Path | str,
    name: str,
    mask: np.ndarray,
    **kwargs: Any,
) -> None:
    """Save one named boolean ROI mask."""

    save_roi_masks(path, {_validate_mask_name(name): mask}, **kwargs)


def _mask_source_label(path: Any) -> str:
    if isinstance(path, (str, Path)):
        return str(path)
    return str(getattr(path, "name", "<stream>"))


def _reset_mask_source(path: Any) -> None:
    if hasattr(path, "seek"):
        path.seek(0)


def load_roi_masks(
    path: Path | str | Any,
    *,
    expected_shape: Optional[Tuple[int, int]] = None,
) -> ROIMaskSet:
    """Load and validate a saved FRAP-Toolbox ROI mask container."""

    input_label = _mask_source_label(path)
    input_source = Path(path) if isinstance(path, (str, Path)) else path
    _reset_mask_source(input_source)
    with np.load(input_source, allow_pickle=False) as data:
        if "metadata_json" not in data.files:
            raise ValueError(f"{input_label} is not a FRAP-Toolbox ROI mask file.")
        metadata = json.loads(str(data["metadata_json"].item()))

        if metadata.get("format") != ROI_MASK_FORMAT:
            raise ValueError(f"Unsupported ROI mask file format: {metadata.get('format')!r}.")
        if int(metadata.get("version", 0)) != ROI_MASK_VERSION:
            raise ValueError(f"Unsupported ROI mask file version: {metadata.get('version')!r}.")

        image_shape = _normalize_image_shape(metadata["image_shape"])
        if expected_shape is not None and image_shape != tuple(expected_shape):
            expected_shape_tuple = tuple(expected_shape)
            raise ValueError(
                f"ROI mask file declares image shape {image_shape}, "
                f"expected {expected_shape_tuple}."
            )

        masks: dict[str, np.ndarray] = {}
        for key in data.files:
            if not key.startswith(_MASK_KEY_PREFIX):
                continue
            name = _validate_mask_name(key[len(_MASK_KEY_PREFIX) :])
            mask = np.asarray(data[key])
            masks[name] = _validate_mask(name, mask, image_shape)

    metadata_names = {entry["name"] for entry in metadata.get("masks", [])}
    loaded_names = set(masks)
    if metadata_names and metadata_names != loaded_names:
        raise ValueError(
            "ROI mask metadata does not match stored arrays: "
            f"metadata={sorted(metadata_names)}, arrays={sorted(loaded_names)}."
        )
    if not masks:
        raise ValueError(f"{input_label} does not contain any ROI mask arrays.")

    _reset_mask_source(input_source)
    return ROIMaskSet(masks=masks, metadata=metadata)


def load_roi_mask(
    path: Path | str | Any,
    name: str,
    *,
    expected_shape: Optional[Tuple[int, int]] = None,
    allow_single: bool = False,
) -> np.ndarray:
    """Load one named boolean ROI mask from a saved container."""

    mask_set = load_roi_masks(path, expected_shape=expected_shape)
    return mask_set.require_mask(name, expected_shape=expected_shape, allow_single=allow_single)
