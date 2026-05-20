from __future__ import annotations

import pytest

from frap_toolbox_py.cloud.manifest import (
    AnalysisSettings,
    OutputManifest,
    RoiGeometry,
    RunManifest,
    RunStatus,
    UploadedFile,
)


def test_run_manifest_round_trips_with_required_cloud_fields():
    manifest = RunManifest(
        run_id="run-1",
        user_id="scientist@example.edu",
        settings=AnalysisSettings(model="diffusion"),
        files=[
            UploadedFile(
                name="Venus_Cytoplasm_1.lsm",
                object_uri="gs://frap-runs/runs/run-1/inputs/Venus_Cytoplasm_1.lsm",
                size_bytes=30_000_000,
                checksum="sha256:example",
            )
        ],
        rois=[RoiGeometry(role="bleach", kind="circle", center_x=256, center_y=23, radius=9)],
        output_prefix="gs://frap-runs/runs/run-1/outputs",
        analysis_engine_version="0.1.0",
    )

    loaded = RunManifest.from_json(manifest.to_json())

    assert loaded.run_id == "run-1"
    assert loaded.settings.model == "diffusion"
    assert loaded.files[0].checksum == "sha256:example"
    assert loaded.rois[0].radius == pytest.approx(9)


def test_manifest_requires_cell_roi_when_whole_cell_normalization_is_enabled():
    manifest = RunManifest(
        run_id="run-1",
        user_id="scientist@example.edu",
        settings=AnalysisSettings(model="reaction1", normalize_by_cell=True),
        files=[UploadedFile(name="Venus_1002.nd2", object_uri="local://frap/run/input.nd2", size_bytes=1)],
        rois=[RoiGeometry(role="bleach", kind="polygon", points=[(0, 0), (1, 0), (1, 1)])],
        output_prefix="local://frap/run/outputs",
        analysis_engine_version="0.1.0",
    )

    with pytest.raises(ValueError, match="cell ROI"):
        manifest.validate()


def test_output_manifest_requires_failed_runs_to_explain_error():
    output = OutputManifest(run_id="run-1", status=RunStatus.FAILED)

    with pytest.raises(ValueError, match="Failed output"):
        output.validate()
