from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from typing import Any, Iterable, Optional


MODEL_KEYS = {"diffusion", "reaction1", "reaction2"}
ROI_KINDS = {"circle", "polygon", "freehand", "mask"}
ROI_ROLES = {"bleach", "cell", "adjacent"}


class RunStatus:
    """Known lifecycle states for a cloud FRAP analysis run."""

    CREATED = "created"
    UPLOADING = "uploading"
    UPLOADED = "uploaded"
    PREVIEW_QUEUED = "preview_queued"
    PREVIEW_READY = "preview_ready"
    SUBMITTED = "submitted"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"

    ALL = {
        CREATED,
        UPLOADING,
        UPLOADED,
        PREVIEW_QUEUED,
        PREVIEW_READY,
        SUBMITTED,
        RUNNING,
        SUCCEEDED,
        FAILED,
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require_finite_number(value: Any, name: str) -> float:
    numeric = float(value)
    if numeric != numeric or numeric in {float("inf"), float("-inf")}:
        raise ValueError(f"{name} must be finite.")
    return numeric


@dataclass
class UploadedFile:
    """A raw microscopy file selected by the user and staged in object storage."""

    name: str
    object_uri: str
    size_bytes: int
    content_type: str = "application/octet-stream"
    checksum: Optional[str] = None

    def validate(self) -> None:
        if not self.name:
            raise ValueError("Uploaded file name is required.")
        if not self.object_uri:
            raise ValueError(f"Uploaded file {self.name!r} is missing object_uri.")
        if self.size_bytes < 0:
            raise ValueError(f"Uploaded file {self.name!r} has a negative size.")


@dataclass
class RoiGeometry:
    """Browser-authored ROI geometry in image pixel coordinates."""

    role: str
    kind: str
    center_x: Optional[float] = None
    center_y: Optional[float] = None
    radius: Optional[float] = None
    points: list[tuple[float, float]] = field(default_factory=list)
    mask_uri: Optional[str] = None
    source_file: Optional[str] = None
    notes: Optional[str] = None

    def validate(self) -> None:
        if self.role not in ROI_ROLES:
            raise ValueError(f"Unsupported ROI role {self.role!r}.")
        if self.kind not in ROI_KINDS:
            raise ValueError(f"Unsupported ROI kind {self.kind!r}.")
        if self.kind == "circle":
            if self.center_x is None or self.center_y is None or self.radius is None:
                raise ValueError("Circle ROI requires center_x, center_y, and radius.")
            _require_finite_number(self.center_x, "center_x")
            _require_finite_number(self.center_y, "center_y")
            radius = _require_finite_number(self.radius, "radius")
            if radius <= 0:
                raise ValueError("Circle ROI radius must be positive.")
        elif self.kind in {"polygon", "freehand"}:
            if len(self.points) < 3:
                raise ValueError(f"{self.kind} ROI requires at least three points.")
            for index, point in enumerate(self.points):
                if len(point) != 2:
                    raise ValueError(f"ROI point {index} must contain x and y.")
                _require_finite_number(point[0], f"points[{index}].x")
                _require_finite_number(point[1], f"points[{index}].y")
        elif self.kind == "mask" and not self.mask_uri:
            raise ValueError("Mask ROI requires mask_uri.")

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RoiGeometry":
        points = data.get("points") or []
        return cls(
            role=data["role"],
            kind=data["kind"],
            center_x=data.get("center_x"),
            center_y=data.get("center_y"),
            radius=data.get("radius"),
            points=[(float(point[0]), float(point[1])) for point in points],
            mask_uri=data.get("mask_uri"),
            source_file=data.get("source_file"),
            notes=data.get("notes"),
        )


@dataclass
class AnalysisSettings:
    """User-controlled analysis settings shared by frontend, API, and worker."""

    model: str = "diffusion"
    post_bleach_frame: int = 21
    pre_bleach_count: int = 10
    background: float = 0.0
    normalize_by_cell: bool = False
    fit_mode: Optional[str] = None
    use_adjacent_roi: bool = False
    adjacent_offset: float = 2.5
    max_profile_radius: Optional[float] = None

    def validate(self) -> None:
        if self.model not in MODEL_KEYS:
            raise ValueError(f"Unsupported model {self.model!r}.")
        if self.post_bleach_frame < 1:
            raise ValueError("post_bleach_frame is 1-based and must be at least 1.")
        if self.pre_bleach_count < 1:
            raise ValueError("pre_bleach_count must be at least 1.")
        if self.background < 0:
            raise ValueError("background must be non-negative.")
        if self.adjacent_offset <= 0:
            raise ValueError("adjacent_offset must be positive.")
        if self.max_profile_radius is not None and self.max_profile_radius <= 0:
            raise ValueError("max_profile_radius must be positive when provided.")

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "AnalysisSettings":
        return cls(**data)


@dataclass
class RunManifest:
    """The durable input contract for one cloud analysis job."""

    run_id: str
    user_id: str
    settings: AnalysisSettings
    files: list[UploadedFile]
    rois: list[RoiGeometry]
    output_prefix: str
    analysis_engine_version: str
    created_at: str = field(default_factory=utc_now)
    storage_mode: str = "managed"

    def validate(self) -> None:
        if not self.run_id:
            raise ValueError("run_id is required.")
        if not self.user_id:
            raise ValueError("user_id is required.")
        self.settings.validate()
        if not self.files:
            raise ValueError("At least one input file is required.")
        for uploaded_file in self.files:
            uploaded_file.validate()
        if not any(roi.role == "bleach" for roi in self.rois):
            raise ValueError("A bleach ROI is required.")
        if self.settings.normalize_by_cell and not any(roi.role == "cell" for roi in self.rois):
            raise ValueError("A cell ROI is required when whole-cell normalization is enabled.")
        for roi in self.rois:
            roi.validate()
        if not self.output_prefix:
            raise ValueError("output_prefix is required.")
        if not self.analysis_engine_version:
            raise ValueError("analysis_engine_version is required.")

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["model"] = self.settings.model
        return data

    def to_json(self) -> str:
        self.validate()
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RunManifest":
        settings_data = data.get("settings")
        if settings_data is None:
            settings_data = {key: data[key] for key in AnalysisSettings.__dataclass_fields__ if key in data}
        manifest = cls(
            run_id=data["run_id"],
            user_id=data["user_id"],
            settings=AnalysisSettings.from_dict(settings_data),
            files=[UploadedFile(**entry) for entry in data.get("files", [])],
            rois=[RoiGeometry.from_dict(entry) for entry in data.get("rois", [])],
            output_prefix=data["output_prefix"],
            analysis_engine_version=data["analysis_engine_version"],
            created_at=data.get("created_at", utc_now()),
            storage_mode=data.get("storage_mode", "managed"),
        )
        manifest.validate()
        return manifest

    @classmethod
    def from_json(cls, raw: str) -> "RunManifest":
        return cls.from_dict(json.loads(raw))


@dataclass
class OutputManifest:
    """The durable output contract written by the worker."""

    run_id: str
    status: str
    parameters: dict[str, float | None] = field(default_factory=dict)
    residual_summaries: dict[str, float] = field(default_factory=dict)
    result_tables: dict[str, str] = field(default_factory=dict)
    plots: dict[str, str] = field(default_factory=dict)
    roi_masks: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    container_version: str = "local"
    input_hashes: dict[str, str] = field(default_factory=dict)
    started_at: str = field(default_factory=utc_now)
    finished_at: Optional[str] = None

    def validate(self) -> None:
        if not self.run_id:
            raise ValueError("run_id is required.")
        if self.status not in RunStatus.ALL:
            raise ValueError(f"Unsupported run status {self.status!r}.")
        if self.status == RunStatus.FAILED and not self.errors:
            raise ValueError("Failed output manifests must include at least one error.")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self) -> str:
        self.validate()
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "OutputManifest":
        manifest = cls(**data)
        manifest.validate()
        return manifest

    @classmethod
    def from_json(cls, raw: str) -> "OutputManifest":
        return cls.from_dict(json.loads(raw))


def coerce_roi_list(rois: Iterable[dict[str, Any] | RoiGeometry]) -> list[RoiGeometry]:
    coerced = [roi if isinstance(roi, RoiGeometry) else RoiGeometry.from_dict(roi) for roi in rois]
    for roi in coerced:
        roi.validate()
    return coerced
