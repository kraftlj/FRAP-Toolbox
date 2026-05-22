from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import json
from pathlib import Path
import shutil
from typing import Any, Optional
from urllib.parse import quote, unquote


@dataclass
class UploadSession:
    """A browser-consumable upload target for one raw input object."""

    file_name: str
    object_uri: str
    upload_url: str
    upload_method: str = "PUT"
    headers: dict[str, str] = field(default_factory=dict)
    expires_at: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "file_name": self.file_name,
            "object_uri": self.object_uri,
            "upload_url": self.upload_url,
            "upload_method": self.upload_method,
            "headers": self.headers,
            "expires_at": self.expires_at,
        }


class StorageBackend(ABC):
    """Object-storage boundary used by the API and worker."""

    @abstractmethod
    def object_uri(self, object_name: str) -> str:
        raise NotImplementedError

    @abstractmethod
    def create_upload_session(
        self,
        object_name: str,
        file_name: str,
        content_type: str = "application/octet-stream",
        size_bytes: Optional[int] = None,
    ) -> UploadSession:
        raise NotImplementedError

    @abstractmethod
    def read_bytes(self, uri: str) -> bytes:
        raise NotImplementedError

    @abstractmethod
    def write_bytes(self, uri: str, data: bytes, content_type: str = "application/octet-stream") -> None:
        raise NotImplementedError

    def read_json(self, uri: str) -> dict[str, Any]:
        return json.loads(self.read_bytes(uri).decode("utf-8"))

    def write_json(self, uri: str, data: dict[str, Any]) -> None:
        self.write_bytes(uri, json.dumps(data, indent=2, sort_keys=True).encode("utf-8"), "application/json")

    def materialize(self, uri: str, destination: Path) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(self.read_bytes(uri))
        return destination


class LocalStorageBackend(StorageBackend):
    """Filesystem-backed storage for tests and local development."""

    def __init__(self, root: Path | str, bucket_name: str = "frap-local"):
        self.root = Path(root).resolve()
        self.bucket_name = bucket_name
        self.root.mkdir(parents=True, exist_ok=True)

    def object_uri(self, object_name: str) -> str:
        normalized = object_name.strip("/")
        return f"local://{self.bucket_name}/{quote(normalized)}"

    def create_upload_session(
        self,
        object_name: str,
        file_name: str,
        content_type: str = "application/octet-stream",
        size_bytes: Optional[int] = None,
    ) -> UploadSession:
        uri = self.object_uri(object_name)
        target = self._path_for_uri(uri)
        target.parent.mkdir(parents=True, exist_ok=True)
        return UploadSession(
            file_name=file_name,
            object_uri=uri,
            upload_url=target.as_uri(),
            upload_method="PUT",
            headers={"content-type": content_type},
        )

    def read_bytes(self, uri: str) -> bytes:
        if uri.startswith("file://"):
            return Path(unquote(uri.removeprefix("file://"))).read_bytes()
        path = self._path_for_uri(uri)
        return path.read_bytes()

    def write_bytes(self, uri: str, data: bytes, content_type: str = "application/octet-stream") -> None:
        path = self._path_for_uri(uri)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)

    def write_range(self, uri: str, start: int, data: bytes, total_size: int) -> int:
        """Write one sequential upload chunk and return the next byte offset."""

        if start < 0:
            raise ValueError("Upload range start must be non-negative.")
        if total_size < 0:
            raise ValueError("Upload total size must be non-negative.")
        if start + len(data) > total_size:
            raise ValueError("Upload chunk exceeds declared total size.")

        path = self._path_for_uri(uri)
        path.parent.mkdir(parents=True, exist_ok=True)
        mode = "r+b" if path.exists() else "wb"
        if start > 0 and not path.exists():
            raise ValueError("Cannot resume upload before the first chunk is written.")

        with path.open(mode) as handle:
            if start == 0:
                handle.truncate(total_size)
            handle.seek(start)
            handle.write(data)
        return start + len(data)

    def materialize(self, uri: str, destination: Path) -> Path:
        if uri.startswith("file://"):
            source = Path(unquote(uri.removeprefix("file://")))
            if source.resolve() == destination.resolve():
                return destination
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            return destination
        if uri.startswith("/") or uri.startswith("."):
            source = Path(uri)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            return destination
        return super().materialize(uri, destination)

    def _path_for_uri(self, uri: str) -> Path:
        prefix = f"local://{self.bucket_name}/"
        if uri.startswith(prefix):
            object_name = unquote(uri.removeprefix(prefix))
            return self.root / object_name
        if uri.startswith("local://"):
            _, _, rest = uri.partition("local://")
            bucket, _, object_name = rest.partition("/")
            if bucket != self.bucket_name:
                raise ValueError(f"URI bucket {bucket!r} does not match local bucket {self.bucket_name!r}.")
            return self.root / unquote(object_name)
        if uri.startswith("/") or uri.startswith("."):
            return Path(uri)
        raise ValueError(f"Unsupported local storage URI: {uri}")


class GCSStorageBackend(StorageBackend):
    """Google Cloud Storage adapter used in deployed environments."""

    def __init__(self, bucket_name: str, upload_origin: Optional[str] = None):
        try:
            from google.cloud import storage
        except ImportError as exc:  # pragma: no cover - optional deployment dependency
            raise ImportError("Install the cloud extra to use GCSStorageBackend.") from exc

        self.bucket_name = bucket_name
        self.upload_origin = upload_origin
        self.client = storage.Client()
        self.bucket = self.client.bucket(bucket_name)

    def object_uri(self, object_name: str) -> str:
        return f"gs://{self.bucket_name}/{object_name.strip('/')}"

    def create_upload_session(
        self,
        object_name: str,
        file_name: str,
        content_type: str = "application/octet-stream",
        size_bytes: Optional[int] = None,
    ) -> UploadSession:
        blob = self.bucket.blob(object_name.strip("/"))
        upload_url = blob.create_resumable_upload_session(
            content_type=content_type,
            size=size_bytes,
            origin=self.upload_origin,
        )
        return UploadSession(
            file_name=file_name,
            object_uri=self.object_uri(object_name),
            upload_url=upload_url,
            upload_method="PUT",
            headers={"content-type": content_type},
        )

    def read_bytes(self, uri: str) -> bytes:
        blob = self._blob_for_uri(uri)
        return blob.download_as_bytes()

    def write_bytes(self, uri: str, data: bytes, content_type: str = "application/octet-stream") -> None:
        blob = self._blob_for_uri(uri)
        blob.upload_from_string(data, content_type=content_type)

    def materialize(self, uri: str, destination: Path) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        self._blob_for_uri(uri).download_to_filename(destination)
        return destination

    def _blob_for_uri(self, uri: str):
        prefix = f"gs://{self.bucket_name}/"
        if not uri.startswith(prefix):
            raise ValueError(f"URI must start with {prefix!r}.")
        return self.bucket.blob(uri.removeprefix(prefix))
