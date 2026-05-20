from __future__ import annotations

from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from typing import Any, Optional

from .manifest import RunStatus, utc_now


@dataclass
class RunRecord:
    """Metadata stored outside raw microscopy files."""

    run_id: str
    user_id: str
    status: str = RunStatus.CREATED
    settings: dict[str, Any] = field(default_factory=dict)
    files: list[dict[str, Any]] = field(default_factory=list)
    rois: list[dict[str, Any]] = field(default_factory=list)
    manifest_uri: Optional[str] = None
    output_manifest_uri: Optional[str] = None
    preview_job_id: Optional[str] = None
    analysis_job_id: Optional[str] = None
    error: Optional[str] = None
    created_at: str = field(default_factory=utc_now)
    updated_at: str = field(default_factory=utc_now)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RunRecord":
        return cls(**data)


class RunRepository(ABC):
    """Metadata persistence boundary for API handlers."""

    @abstractmethod
    def create(self, record: RunRecord) -> RunRecord:
        raise NotImplementedError

    @abstractmethod
    def get(self, run_id: str) -> Optional[RunRecord]:
        raise NotImplementedError

    @abstractmethod
    def update(self, run_id: str, **changes: Any) -> RunRecord:
        raise NotImplementedError


class InMemoryRunRepository(RunRepository):
    """Tiny repository for tests and local API development."""

    def __init__(self):
        self._records: dict[str, RunRecord] = {}

    def create(self, record: RunRecord) -> RunRecord:
        if record.run_id in self._records:
            raise ValueError(f"Run {record.run_id!r} already exists.")
        self._records[record.run_id] = deepcopy(record)
        return deepcopy(record)

    def get(self, run_id: str) -> Optional[RunRecord]:
        record = self._records.get(run_id)
        return deepcopy(record) if record is not None else None

    def update(self, run_id: str, **changes: Any) -> RunRecord:
        record = self._records.get(run_id)
        if record is None:
            raise KeyError(run_id)
        for key, value in changes.items():
            if not hasattr(record, key):
                raise AttributeError(f"RunRecord has no field {key!r}.")
            setattr(record, key, value)
        record.updated_at = utc_now()
        self._records[run_id] = deepcopy(record)
        return deepcopy(record)


class FirestoreRunRepository(RunRepository):
    """Firestore-backed repository for Cloud Run deployment."""

    def __init__(self, collection_name: str = "frap_runs"):
        try:
            from google.cloud import firestore
        except ImportError as exc:  # pragma: no cover - optional deployment dependency
            raise ImportError("Install the cloud extra to use FirestoreRunRepository.") from exc

        self.client = firestore.Client()
        self.collection = self.client.collection(collection_name)

    def create(self, record: RunRecord) -> RunRecord:
        self.collection.document(record.run_id).set(record.to_dict())
        return record

    def get(self, run_id: str) -> Optional[RunRecord]:
        snapshot = self.collection.document(run_id).get()
        if not snapshot.exists:
            return None
        return RunRecord.from_dict(snapshot.to_dict())

    def update(self, run_id: str, **changes: Any) -> RunRecord:
        changes["updated_at"] = utc_now()
        document = self.collection.document(run_id)
        document.update(changes)
        snapshot = document.get()
        return RunRecord.from_dict(snapshot.to_dict())
