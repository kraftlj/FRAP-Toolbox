from __future__ import annotations

from frap_toolbox_py.cloud.manifest import AnalysisSettings, RoiGeometry, RunManifest, UploadedFile
from frap_toolbox_py.cloud.storage import LocalStorageBackend
from frap_toolbox_py.cloud.worker import run_analysis


def test_worker_dry_run_writes_output_manifest(tmp_path):
    storage = LocalStorageBackend(tmp_path)
    manifest = RunManifest(
        run_id="run-1",
        user_id="pilot@example.edu",
        settings=AnalysisSettings(model="diffusion"),
        files=[UploadedFile(name="example.lsm", object_uri=storage.object_uri("inputs/example.lsm"), size_bytes=1)],
        rois=[RoiGeometry(role="bleach", kind="circle", center_x=5, center_y=5, radius=2)],
        output_prefix=storage.object_uri("outputs/run-1"),
        analysis_engine_version="0.1.0",
    )
    manifest_uri = storage.object_uri("manifests/run-1.json")
    storage.write_json(manifest_uri, manifest.to_dict())

    output = run_analysis(manifest_uri, storage=storage, dry_run=True)

    assert output.status == "succeeded"
    stored = storage.read_json(storage.object_uri("outputs/run-1/output-manifest.json"))
    assert stored["warnings"] == ["Dry run completed without reading microscopy inputs."]
