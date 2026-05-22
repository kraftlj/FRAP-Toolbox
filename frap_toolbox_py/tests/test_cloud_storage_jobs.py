from __future__ import annotations

from frap_toolbox_py.cloud.jobs import build_batch_job_body
from frap_toolbox_py.cloud.storage import LocalStorageBackend


def test_local_storage_generates_uri_and_upload_session(tmp_path):
    storage = LocalStorageBackend(tmp_path)

    session = storage.create_upload_session(
        "runs/user/run-1/inputs/example.nd2",
        "example.nd2",
        size_bytes=42,
    )

    assert session.object_uri == "local://frap-local/runs/user/run-1/inputs/example.nd2"
    assert session.upload_method == "PUT"
    assert session.upload_url.startswith("file://")


def test_batch_body_runs_worker_with_manifest_and_output_prefix():
    body = build_batch_job_body(
        job_id="frap-analysis-run-1",
        container_image="us-docker.pkg.dev/frap/worker:latest",
        manifest_uri="gs://frap-runs/runs/run-1/manifest.json",
        output_prefix="gs://frap-runs/runs/run-1/outputs",
        service_account="worker@example.iam.gserviceaccount.com",
    )

    runnable = body["taskGroups"][0]["taskSpec"]["runnables"][0]["container"]

    assert runnable["imageUri"] == "us-docker.pkg.dev/frap/worker:latest"
    assert runnable["commands"] == [
        "frap-toolbox-run",
        "--manifest",
        "gs://frap-runs/runs/run-1/manifest.json",
        "--output-prefix",
        "gs://frap-runs/runs/run-1/outputs",
    ]
    assert body["allocationPolicy"]["serviceAccount"]["email"] == "worker@example.iam.gserviceaccount.com"
