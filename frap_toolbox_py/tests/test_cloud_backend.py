from __future__ import annotations

from fastapi.testclient import TestClient

from frap_toolbox_py.cloud.backend import create_app
from frap_toolbox_py.cloud.jobs import LocalJobBackend
from frap_toolbox_py.cloud.repository import InMemoryRunRepository
from frap_toolbox_py.cloud.storage import LocalStorageBackend


def _client(tmp_path):
    jobs = LocalJobBackend()
    app = create_app(
        repository=InMemoryRunRepository(),
        storage=LocalStorageBackend(tmp_path),
        jobs=jobs,
        allowed_users={"pilot@example.edu"},
    )
    return TestClient(app), jobs


def test_backend_requires_invited_user(tmp_path):
    client, _ = _client(tmp_path)

    response = client.post("/runs", headers={"x-frap-user": "stranger@example.edu"}, json={})

    assert response.status_code == 403


def test_create_upload_submit_and_fetch_run(tmp_path):
    client, jobs = _client(tmp_path)
    headers = {"x-frap-user": "pilot@example.edu"}

    created = client.post("/runs", headers=headers, json={"settings": {"model": "diffusion"}})
    assert created.status_code == 201
    run_id = created.json()["run_id"]

    uploads = client.post(
        f"/runs/{run_id}/uploads",
        headers=headers,
        json={"files": [{"file_name": "Venus_Cytoplasm_1.lsm", "size_bytes": 30_000_000}]},
    )
    assert uploads.status_code == 200
    upload_session = uploads.json()["upload_sessions"][0]
    assert upload_session["object_uri"].endswith("/0000-Venus_Cytoplasm_1.lsm")
    assert upload_session["upload_url"].startswith("/local-uploads/")

    submitted = client.post(
        f"/runs/{run_id}/submit",
        headers=headers,
        json={
            "rois": [
                {
                    "role": "bleach",
                    "kind": "circle",
                    "center_x": 256,
                    "center_y": 23,
                    "radius": 9,
                }
            ]
        },
    )

    assert submitted.status_code == 200
    body = submitted.json()
    assert body["run"]["status"] == "submitted"
    assert body["manifest_uri"].endswith("/manifest.json")
    assert jobs.submissions[-1].details["manifest_uri"] == body["manifest_uri"]

    fetched = client.get(f"/runs/{run_id}", headers=headers)
    assert fetched.status_code == 200
    assert fetched.json()["analysis_job_id"].startswith("local-analysis-")


def test_local_upload_endpoint_accepts_resumable_chunks(tmp_path):
    client, _ = _client(tmp_path)
    headers = {"x-frap-user": "pilot@example.edu"}
    created = client.post("/runs", headers=headers, json={})
    run_id = created.json()["run_id"]
    uploads = client.post(
        f"/runs/{run_id}/uploads",
        headers=headers,
        json={"files": [{"file_name": "tiny.lsm", "size_bytes": 10}]},
    )
    session = uploads.json()["upload_sessions"][0]

    first = client.put(
        session["upload_url"],
        content=b"abcd",
        headers={"content-range": "bytes 0-3/10"},
    )
    second = client.put(
        session["upload_url"],
        content=b"efghij",
        headers={"content-range": "bytes 4-9/10"},
    )

    assert first.status_code == 308
    assert first.headers["range"] == "bytes=0-3"
    assert second.status_code == 200
    assert (tmp_path / "runs" / "pilot@example.edu" / run_id / "inputs" / "0000-tiny.lsm").read_bytes() == b"abcdefghij"
