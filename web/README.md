# FRAP Toolbox Web Client

This is the cloud-first responsive TypeScript client for the FRAP-Toolbox pilot workflow.

```bash
npm install
npm run dev
```

The Vite dev server proxies `/api` to `http://127.0.0.1:8000`, where the FastAPI app can run with:

```bash
uvicorn frap_toolbox_py.cloud.backend:app --reload
```

For deployed builds, set `VITE_FRAP_API_BASE` to the Cloud Run API URL.
