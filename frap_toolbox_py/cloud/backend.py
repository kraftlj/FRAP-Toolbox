from __future__ import annotations

from dataclasses import asdict
import os
from pathlib import Path
import re
from typing import Any, Optional
import uuid

try:
    from fastapi import Depends, FastAPI, Header, HTTPException, Request, Response, status
    from fastapi.middleware.cors import CORSMiddleware
    from pydantic import BaseModel, Field
except ImportError as exc:  # pragma: no cover - optional app dependency
    raise ImportError("Install the cloud-api extra to run the FastAPI backend.") from exc

from .jobs import GCPBatchJobBackend, JobBackend, LocalJobBackend
from .manifest import AnalysisSettings, RoiGeometry, RunManifest, RunStatus, UploadedFile
from .repository import FirestoreRunRepository, InMemoryRunRepository, RunRecord, RunRepository
from .storage import GCSStorageBackend, LocalStorageBackend, StorageBackend, UploadSession
from .worker import join_uri, package_version


def _model_dump(model: BaseModel) -> dict[str, Any]:
    if hasattr(model, "model_dump"):
        return model.model_dump(exclude_none=True)
    return model.dict(exclude_none=True)


class AnalysisSettingsIn(BaseModel):
    model: str = "diffusion"
    post_bleach_frame: int = 21
    pre_bleach_count: int = 10
    background: float = 0.0
    normalize_by_cell: bool = False
    fit_mode: Optional[str] = None
    use_adjacent_roi: bool = False
    adjacent_offset: float = 2.5
    max_profile_radius: Optional[float] = None


class CreateRunRequest(BaseModel):
    settings: AnalysisSettingsIn = Field(default_factory=AnalysisSettingsIn)


class UploadFileRequest(BaseModel):
    file_name: str
    size_bytes: int = 0
    content_type: str = "application/octet-stream"
    checksum: Optional[str] = None


class CreateUploadsRequest(BaseModel):
    files: list[UploadFileRequest]


class RoiGeometryIn(BaseModel):
    role: str
    kind: str
    center_x: Optional[float] = None
    center_y: Optional[float] = None
    radius: Optional[float] = None
    points: list[tuple[float, float]] = Field(default_factory=list)
    mask_uri: Optional[str] = None
    source_file: Optional[str] = None
    notes: Optional[str] = None


class SubmitRunRequest(BaseModel):
    settings: Optional[AnalysisSettingsIn] = None
    rois: list[RoiGeometryIn]


class PreviewRequest(BaseModel):
    input_uri: Optional[str] = None


class AppState:
    def __init__(
        self,
        repository: RunRepository,
        storage: StorageBackend,
        jobs: JobBackend,
        allowed_users: set[str],
    ):
        self.repository = repository
        self.storage = storage
        self.jobs = jobs
        self.allowed_users = allowed_users
        self.local_uploads: dict[str, dict[str, Any]] = {}


def _allowed_users_from_env() -> set[str]:
    raw = os.getenv("FRAP_ALLOWED_USERS", "")
    return {entry.strip().lower() for entry in raw.split(",") if entry.strip()}


def _normalize_iap_email(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    if ":" in value:
        return value.rsplit(":", 1)[-1].lower()
    return value.lower()


_CONTENT_RANGE_PATTERN = re.compile(r"^bytes (?P<start>\d+)-(?P<end>\d+)/(?P<total>\d+)$")


def _local_upload_session(state: AppState, session: UploadSession, size_bytes: int) -> UploadSession:
    token = uuid.uuid4().hex
    state.local_uploads[token] = {
        "object_uri": session.object_uri,
        "size_bytes": size_bytes,
        "received_bytes": 0,
    }
    return UploadSession(
        file_name=session.file_name,
        object_uri=session.object_uri,
        upload_url=f"/local-uploads/{token}",
        upload_method="PUT",
        headers=session.headers,
        expires_at=session.expires_at,
    )


def _default_repository() -> RunRepository:
    if os.getenv("FRAP_METADATA_BACKEND", "local").lower() == "firestore":
        return FirestoreRunRepository(os.getenv("FRAP_FIRESTORE_COLLECTION", "frap_runs"))
    return InMemoryRunRepository()


def _default_storage() -> StorageBackend:
    if os.getenv("FRAP_STORAGE_BACKEND", "local").lower() == "gcs":
        bucket = os.environ["FRAP_GCS_BUCKET"]
        return GCSStorageBackend(bucket, upload_origin=os.getenv("FRAP_UPLOAD_ORIGIN"))
    return LocalStorageBackend(Path(os.getenv("FRAP_LOCAL_STORAGE_ROOT", ".frap-cloud-storage")))


def _default_jobs() -> JobBackend:
    if os.getenv("FRAP_JOB_BACKEND", "local").lower() == "google-batch":
        return GCPBatchJobBackend(
            project_id=os.environ["FRAP_GCP_PROJECT"],
            region=os.getenv("FRAP_GCP_REGION", "us-central1"),
            container_image=os.environ["FRAP_WORKER_IMAGE"],
            service_account=os.getenv("FRAP_BATCH_SERVICE_ACCOUNT"),
        )
    return LocalJobBackend()


def create_app(
    *,
    repository: Optional[RunRepository] = None,
    storage: Optional[StorageBackend] = None,
    jobs: Optional[JobBackend] = None,
    allowed_users: Optional[set[str]] = None,
) -> FastAPI:
    repository = repository or _default_repository()
    storage = storage or _default_storage()
    jobs = jobs or _default_jobs()
    state = AppState(
        repository=repository,
        storage=storage,
        jobs=jobs,
        allowed_users=allowed_users if allowed_users is not None else _allowed_users_from_env(),
    )

    app = FastAPI(title="FRAP Toolbox Cloud API", version=package_version())
    app.state.frap = state
    app.add_middleware(
        CORSMiddleware,
        allow_origins=[origin.strip() for origin in os.getenv("FRAP_CORS_ORIGINS", "*").split(",")],
        allow_credentials=True,
        allow_methods=["GET", "POST", "OPTIONS"],
        allow_headers=["*"],
    )

    def current_user(
        x_frap_user: Optional[str] = Header(default=None),
        x_goog_authenticated_user_email: Optional[str] = Header(default=None),
    ) -> str:
        email = _normalize_iap_email(x_goog_authenticated_user_email) or _normalize_iap_email(x_frap_user)
        if not email:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Authenticated user header required.")
        if state.allowed_users and email not in state.allowed_users:
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="User is not invited to this pilot.")
        return email

    def require_run(run_id: str, user_id: str) -> RunRecord:
        record = state.repository.get(run_id)
        if record is None:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Run not found.")
        if record.user_id != user_id:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Run not found.")
        return record

    @app.get("/healthz")
    def healthz() -> dict[str, str]:
        return {"status": "ok"}

    @app.post("/runs", status_code=status.HTTP_201_CREATED)
    def create_run(request: CreateRunRequest, user_id: str = Depends(current_user)) -> dict[str, Any]:
        settings = AnalysisSettings.from_dict(_model_dump(request.settings))
        try:
            settings.validate()
        except ValueError as exc:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(exc)) from exc

        run_id = uuid.uuid4().hex
        record = RunRecord(
            run_id=run_id,
            user_id=user_id,
            settings=asdict(settings),
        )
        return state.repository.create(record).to_dict()

    @app.post("/runs/{run_id}/uploads")
    def create_uploads(
        run_id: str,
        request: CreateUploadsRequest,
        user_id: str = Depends(current_user),
    ) -> dict[str, Any]:
        record = require_run(run_id, user_id)
        if not request.files:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="At least one file is required.")

        sessions = []
        file_records = []
        for index, file_request in enumerate(request.files):
            safe_name = Path(file_request.file_name).name
            object_name = f"runs/{user_id}/{run_id}/inputs/{index:04d}-{safe_name}"
            session = state.storage.create_upload_session(
                object_name,
                safe_name,
                content_type=file_request.content_type,
                size_bytes=file_request.size_bytes,
            )
            if isinstance(state.storage, LocalStorageBackend):
                session = _local_upload_session(state, session, file_request.size_bytes)
            sessions.append(session.to_dict())
            file_records.append(
                {
                    "name": safe_name,
                    "object_uri": session.object_uri,
                    "size_bytes": file_request.size_bytes,
                    "content_type": file_request.content_type,
                    "checksum": file_request.checksum,
                }
            )

        updated = state.repository.update(
            record.run_id,
            files=file_records,
            status=RunStatus.UPLOADING,
        )
        return {"run": updated.to_dict(), "upload_sessions": sessions}

    @app.put("/local-uploads/{token}")
    async def local_upload(token: str, request: Request) -> Response:
        if not isinstance(state.storage, LocalStorageBackend):
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Local uploads are disabled.")

        target = state.local_uploads.get(token)
        if target is None:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Upload session not found.")

        body = await request.body()
        content_range = request.headers.get("content-range")
        if content_range is None:
            state.storage.write_bytes(target["object_uri"], body)
            target["received_bytes"] = len(body)
            return Response(status_code=status.HTTP_200_OK)

        match = _CONTENT_RANGE_PATTERN.match(content_range)
        if match is None:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid Content-Range header.")

        start = int(match.group("start"))
        end = int(match.group("end"))
        total = int(match.group("total"))
        if end < start or len(body) != end - start + 1:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Upload chunk size does not match Content-Range.")
        if target["size_bytes"] and total != target["size_bytes"]:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Upload total does not match prepared session.")

        next_offset = state.storage.write_range(target["object_uri"], start, body, total)
        target["received_bytes"] = max(int(target["received_bytes"]), next_offset)
        if next_offset < total:
            return Response(status_code=308, headers={"Range": f"bytes=0-{next_offset - 1}"})

        return Response(status_code=status.HTTP_200_OK)

    @app.post("/runs/{run_id}/preview")
    def queue_preview(
        run_id: str,
        request: PreviewRequest,
        user_id: str = Depends(current_user),
    ) -> dict[str, Any]:
        record = require_run(run_id, user_id)
        if not record.files:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Upload files before preview.")
        input_uri = request.input_uri or record.files[0]["object_uri"]
        output_prefix = state.storage.object_uri(f"runs/{user_id}/{run_id}/preview")
        submission = state.jobs.submit_preview_job(record.run_id, input_uri, output_prefix)
        updated = state.repository.update(
            record.run_id,
            preview_job_id=submission.job_id,
            status=RunStatus.PREVIEW_QUEUED,
        )
        return {"run": updated.to_dict(), "job": submission.__dict__}

    @app.post("/runs/{run_id}/submit")
    def submit_run(
        run_id: str,
        request: SubmitRunRequest,
        user_id: str = Depends(current_user),
    ) -> dict[str, Any]:
        record = require_run(run_id, user_id)
        if not record.files:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Upload files before analysis.")

        settings = AnalysisSettings.from_dict(_model_dump(request.settings) if request.settings else record.settings)
        rois = [RoiGeometry.from_dict(_model_dump(roi)) for roi in request.rois]
        output_prefix = state.storage.object_uri(f"runs/{user_id}/{run_id}/outputs")
        manifest = RunManifest(
            run_id=record.run_id,
            user_id=user_id,
            settings=settings,
            files=[UploadedFile(**file_record) for file_record in record.files],
            rois=rois,
            output_prefix=output_prefix,
            analysis_engine_version=package_version(),
        )
        try:
            manifest.validate()
        except ValueError as exc:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(exc)) from exc

        manifest_uri = state.storage.object_uri(f"runs/{user_id}/{run_id}/manifest.json")
        state.storage.write_json(manifest_uri, manifest.to_dict())
        output_manifest_uri = join_uri(output_prefix, "output-manifest.json")
        submission = state.jobs.submit_analysis_job(record.run_id, manifest_uri, output_prefix)
        updated = state.repository.update(
            record.run_id,
            settings=asdict(settings),
            rois=[asdict(roi) for roi in rois],
            manifest_uri=manifest_uri,
            output_manifest_uri=output_manifest_uri,
            analysis_job_id=submission.job_id,
            status=RunStatus.SUBMITTED,
        )
        return {"run": updated.to_dict(), "job": submission.__dict__, "manifest_uri": manifest_uri}

    @app.get("/runs/{run_id}")
    def get_run(run_id: str, user_id: str = Depends(current_user)) -> dict[str, Any]:
        return require_run(run_id, user_id).to_dict()

    @app.get("/runs/{run_id}/results")
    def get_results(run_id: str, user_id: str = Depends(current_user)) -> dict[str, Any]:
        record = require_run(run_id, user_id)
        if not record.output_manifest_uri:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Results are not available yet.")
        try:
            return state.storage.read_json(record.output_manifest_uri)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Results are not available yet.") from exc

    return app


app = create_app()
