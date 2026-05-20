from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional
import uuid


@dataclass
class JobSubmission:
    job_id: str
    provider: str
    status: str = "queued"
    details: dict[str, Any] | None = None


class JobBackend(ABC):
    """Batch job boundary used by the API."""

    @abstractmethod
    def submit_preview_job(self, run_id: str, input_uri: str, output_prefix: str) -> JobSubmission:
        raise NotImplementedError

    @abstractmethod
    def submit_analysis_job(self, run_id: str, manifest_uri: str, output_prefix: str) -> JobSubmission:
        raise NotImplementedError


class LocalJobBackend(JobBackend):
    """Records job submissions without executing them."""

    def __init__(self):
        self.submissions: list[JobSubmission] = []

    def submit_preview_job(self, run_id: str, input_uri: str, output_prefix: str) -> JobSubmission:
        submission = JobSubmission(
            job_id=f"local-preview-{uuid.uuid4().hex[:12]}",
            provider="local",
            details={"run_id": run_id, "input_uri": input_uri, "output_prefix": output_prefix},
        )
        self.submissions.append(submission)
        return submission

    def submit_analysis_job(self, run_id: str, manifest_uri: str, output_prefix: str) -> JobSubmission:
        submission = JobSubmission(
            job_id=f"local-analysis-{uuid.uuid4().hex[:12]}",
            provider="local",
            details={"run_id": run_id, "manifest_uri": manifest_uri, "output_prefix": output_prefix},
        )
        self.submissions.append(submission)
        return submission


def build_batch_job_body(
    *,
    job_id: str,
    container_image: str,
    manifest_uri: str,
    output_prefix: str,
    machine_type: str = "e2-standard-4",
    boot_disk_gb: int = 100,
    service_account: Optional[str] = None,
) -> dict[str, Any]:
    """Build a Google Batch REST job body for one FRAP worker run."""

    body: dict[str, Any] = {
        "name": job_id,
        "labels": {"app": "frap-toolbox", "kind": "analysis"},
        "taskGroups": [
            {
                "taskSpec": {
                    "runnables": [
                        {
                            "container": {
                                "imageUri": container_image,
                                "commands": [
                                    "frap-toolbox-run",
                                    "--manifest",
                                    manifest_uri,
                                    "--output-prefix",
                                    output_prefix,
                                ],
                            }
                        }
                    ],
                    "computeResource": {
                        "cpuMilli": 4000,
                        "memoryMib": 16384,
                        "bootDiskMib": boot_disk_gb * 1024,
                    },
                    "maxRetryCount": 1,
                    "maxRunDuration": "7200s",
                },
                "taskCount": 1,
                "parallelism": 1,
            }
        ],
        "allocationPolicy": {
            "instances": [
                {
                    "policy": {
                        "machineType": machine_type,
                    }
                }
            ]
        },
        "logsPolicy": {"destination": "CLOUD_LOGGING"},
    }
    if service_account:
        body["allocationPolicy"]["serviceAccount"] = {"email": service_account}
    return body


class GCPBatchJobBackend(JobBackend):
    """Google Batch adapter for deployed analysis jobs."""

    def __init__(
        self,
        project_id: str,
        region: str,
        container_image: str,
        service_account: Optional[str] = None,
    ):
        try:
            from google.cloud import batch_v1
            from google.protobuf import json_format
        except ImportError as exc:  # pragma: no cover - optional deployment dependency
            raise ImportError("Install the cloud extra to use GCPBatchJobBackend.") from exc

        self.project_id = project_id
        self.region = region
        self.container_image = container_image
        self.service_account = service_account
        self.batch_v1 = batch_v1
        self.json_format = json_format
        self.client = batch_v1.BatchServiceClient()

    def submit_preview_job(self, run_id: str, input_uri: str, output_prefix: str) -> JobSubmission:
        # The v1 backend queues preview extraction through the same container.
        manifest_uri = f"{output_prefix.rstrip('/')}/preview-manifest.json"
        return self._submit(run_id, manifest_uri, output_prefix, kind="preview")

    def submit_analysis_job(self, run_id: str, manifest_uri: str, output_prefix: str) -> JobSubmission:
        return self._submit(run_id, manifest_uri, output_prefix, kind="analysis")

    def _submit(self, run_id: str, manifest_uri: str, output_prefix: str, kind: str) -> JobSubmission:
        job_id = f"frap-{kind}-{run_id}".lower().replace("_", "-")
        body = build_batch_job_body(
            job_id=job_id,
            container_image=self.container_image,
            manifest_uri=manifest_uri,
            output_prefix=output_prefix,
            service_account=self.service_account,
        )
        parent = f"projects/{self.project_id}/locations/{self.region}"
        job = self.json_format.ParseDict(body, self.batch_v1.Job())
        response = self.client.create_job(parent=parent, job_id=job_id, job=job)
        return JobSubmission(
            job_id=response.name,
            provider="google-batch",
            details={"run_id": run_id, "manifest_uri": manifest_uri, "output_prefix": output_prefix},
        )
