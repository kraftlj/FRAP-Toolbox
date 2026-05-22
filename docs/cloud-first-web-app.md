# Cloud-First Web App

This document describes the first implementation of the FRAP-Toolbox cloud workflow.
The design follows sequencing-style analysis: raw microscopy files are uploaded to
object storage, analysis runs as a containerized batch job, and the web app tracks
metadata and results.

## Components

- **Web client**: `web/` is a Vite/React/TypeScript app. It creates runs, asks the
  API for upload sessions, uploads files directly to storage, collects ROI geometry,
  submits analysis, and shows result-manifest values.
- **API service**: `frap_toolbox_py.cloud.backend` is a FastAPI app intended for
  Cloud Run. It handles invite-gated identity, run metadata, upload-session
  creation, preview/job queueing, and result lookup.
- **Worker**: `frap_toolbox_py.cloud.worker` provides the `frap-toolbox-run`
  command. It reads a run manifest, materializes inputs, calls the existing Python
  FRAP analysis engine, and writes result artifacts plus `output-manifest.json`.
- **Storage/job boundaries**: local adapters support tests and development. GCS,
  Firestore, and Google Batch adapters are selected by environment variables in
  deployment.

## Local Development

Run the API in local-adapter mode:

```bash
python -m pip install -e ".[test]"
uvicorn frap_toolbox_py.cloud.backend:app --reload
```

Run the web client:

```bash
cd web
npm install
npm run dev
```

The frontend proxies `/api` to `http://127.0.0.1:8000`.

The local API returns tokenized `/local-uploads/...` URLs that accept the same
chunked `PUT` flow as the browser client. Production uses GCS resumable upload
session URLs so raw microscopy files still bypass the Cloud Run API.

## Web Deployment

Build the static client with:

```bash
cd web
npm ci
VITE_FRAP_API_BASE=https://<cloud-run-api-url> npm run build
```

Serve `web/dist/` from a static host such as Firebase Hosting, Cloud Storage
static hosting behind HTTPS, or the same load-balanced domain as the API.

## Cloud Run API

Build and deploy the API image from `deploy/cloud-run-api.Dockerfile`. Required
environment:

- `FRAP_ALLOWED_USERS`: comma-separated invited user emails.
- `FRAP_STORAGE_BACKEND=gcs`
- `FRAP_GCS_BUCKET`: managed pilot bucket.
- `FRAP_UPLOAD_ORIGIN`: browser origin allowed to use resumable upload sessions.
- `FRAP_METADATA_BACKEND=firestore`
- `FRAP_FIRESTORE_COLLECTION=frap_runs`
- `FRAP_JOB_BACKEND=google-batch`
- `FRAP_GCP_PROJECT`
- `FRAP_GCP_REGION`, default `us-central1`
- `FRAP_WORKER_IMAGE`: Artifact Registry image for the Batch worker.
- `FRAP_BATCH_SERVICE_ACCOUNT`: service account used by worker jobs.
- `FRAP_CORS_ORIGINS`: deployed web origin.

Recommended Cloud Run front door:

- Put the service behind Google-authenticated access, such as Identity-Aware Proxy
  or an equivalent load balancer/auth layer.
- Pass the authenticated email to the API through `X-Goog-Authenticated-User-Email`.
- Keep `X-FRAP-USER` for local development and tests only.

Configure bucket CORS for the deployed web origin so the browser can upload
directly to GCS resumable session URLs. The upload client sends `PUT` requests
with `Content-Type` and `Content-Range` headers.

## Batch Worker

Build the worker image from `deploy/batch-worker.Dockerfile` and push it to Artifact
Registry. The API submits Google Batch jobs that run:

```bash
frap-toolbox-run --manifest gs://.../manifest.json --output-prefix gs://.../outputs
```

The worker writes:

- `result-parameters.csv`
- `frap-fit.png`
- `post-bleach-profile.png` for diffusion runs
- `roi-masks.npz`
- `output-manifest.json`

## Pilot Storage Policy

For the managed pilot path, raw uploads are retained by default under:

```text
gs://<bucket>/runs/<user-email>/<run-id>/inputs/
```

Results are written under:

```text
gs://<bucket>/runs/<user-email>/<run-id>/outputs/
```

Do not enable automatic raw-data deletion for the pilot unless the retention policy
changes. Use bucket labels, object prefixes, and Firestore run records for quota and
manual cleanup reporting.

## BYO GCP Path

The API, storage, repository, and job layers are already split behind small adapter
interfaces. A future BYO GCP version should add a per-lab storage/job profile and
select adapters per run instead of from process-wide environment variables.
