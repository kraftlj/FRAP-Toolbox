import {
  Activity,
  CheckCircle2,
  CloudUpload,
  FileUp,
  FlaskConical,
  Loader2,
  Play,
  RefreshCw,
  ShieldCheck
} from "lucide-react";
import type React from "react";
import { useMemo, useState } from "react";
import {
  AnalysisSettings,
  OutputManifest,
  RoiGeometry,
  RunRecord,
  UploadSession,
  createRun,
  createUploadSessions,
  getResults,
  getRun,
  queuePreview,
  submitRun,
  uploadFileResumable
} from "./api";

type Stage = "setup" | "upload" | "roi" | "submit" | "results";
type UploadProgress = Record<string, number>;

const defaultSettings: AnalysisSettings = {
  model: "diffusion",
  post_bleach_frame: 21,
  pre_bleach_count: 10,
  background: 0,
  normalize_by_cell: false,
  fit_mode: "global",
  use_adjacent_roi: false,
  adjacent_offset: 2.5
};

function fitModes(model: AnalysisSettings["model"]): string[] {
  if (model === "diffusion") {
    return ["global", "individual", "average_curve", "simplified_kang", "simplified_kang_global"];
  }
  return ["individual", "average_curve"];
}

function bytesLabel(bytes: number): string {
  if (bytes < 1024 * 1024) return `${Math.round(bytes / 1024)} KB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(2)} GB`;
}

function Stepper({ stage }: { stage: Stage }) {
  const steps: Array<[Stage, string]> = [
    ["setup", "Setup"],
    ["upload", "Upload"],
    ["roi", "ROI"],
    ["submit", "Run"],
    ["results", "Results"]
  ];
  const activeIndex = steps.findIndex(([key]) => key === stage);
  return (
    <div className="stepper" aria-label="Run workflow">
      {steps.map(([key, label], index) => (
        <div className={`step ${index <= activeIndex ? "active" : ""}`} key={key}>
          <span>{index + 1}</span>
          <strong>{label}</strong>
        </div>
      ))}
    </div>
  );
}

function RoiEditor({
  rois,
  setRois,
  normalizeByCell
}: {
  rois: RoiGeometry[];
  setRois: (rois: RoiGeometry[]) => void;
  normalizeByCell: boolean;
}) {
  const bleach = rois.find((roi) => roi.role === "bleach") ?? {
    role: "bleach",
    kind: "circle",
    center_x: 256,
    center_y: 256,
    radius: 24
  };
  const cell = rois.find((roi) => roi.role === "cell") ?? {
    role: "cell",
    kind: "circle",
    center_x: 256,
    center_y: 256,
    radius: 160
  };

  function update(next: RoiGeometry) {
    const withoutRole = rois.filter((roi) => roi.role !== next.role);
    const nextRois = [...withoutRole, next];
    setRois(normalizeByCell ? nextRois : nextRois.filter((roi) => roi.role !== "cell"));
  }

  function placeBleach(event: React.MouseEvent<SVGSVGElement>) {
    const rect = event.currentTarget.getBoundingClientRect();
    const x = ((event.clientX - rect.left) / rect.width) * 512;
    const y = ((event.clientY - rect.top) / rect.height) * 512;
    update({ ...bleach, center_x: Math.round(x), center_y: Math.round(y) });
  }

  return (
    <div className="roi-grid">
      <section className="panel roi-canvas-panel">
        <div className="panel-title">
          <Activity size={18} />
          <h2>ROI Preview</h2>
        </div>
        <svg viewBox="0 0 512 512" className="roi-canvas" onClick={placeBleach} role="img" aria-label="ROI canvas">
          <defs>
            <pattern id="grid" width="32" height="32" patternUnits="userSpaceOnUse">
              <path d="M 32 0 L 0 0 0 32" fill="none" stroke="rgba(43, 54, 67, 0.12)" strokeWidth="1" />
            </pattern>
          </defs>
          <rect width="512" height="512" fill="#f5f8fa" />
          <rect width="512" height="512" fill="url(#grid)" />
          {normalizeByCell && (
            <circle
              cx={cell.center_x}
              cy={cell.center_y}
              r={cell.radius}
              className="roi-cell"
            />
          )}
          <circle cx={bleach.center_x} cy={bleach.center_y} r={bleach.radius} className="roi-bleach" />
          <line x1={bleach.center_x} x2={bleach.center_x} y1="0" y2="512" className="roi-guide" />
          <line y1={bleach.center_y} y2={bleach.center_y} x1="0" x2="512" className="roi-guide" />
        </svg>
      </section>
      <section className="panel">
        <div className="panel-title">
          <FlaskConical size={18} />
          <h2>ROI Geometry</h2>
        </div>
        <div className="field-grid two">
          <label>
            Bleach X
            <input
              type="number"
              value={bleach.center_x}
              onChange={(event) => update({ ...bleach, center_x: Number(event.target.value) })}
            />
          </label>
          <label>
            Bleach Y
            <input
              type="number"
              value={bleach.center_y}
              onChange={(event) => update({ ...bleach, center_y: Number(event.target.value) })}
            />
          </label>
          <label>
            Bleach radius
            <input
              type="number"
              min="1"
              value={bleach.radius}
              onChange={(event) => update({ ...bleach, radius: Number(event.target.value) })}
            />
          </label>
          {normalizeByCell && (
            <label>
              Cell radius
              <input
                type="number"
                min="1"
                value={cell.radius}
                onChange={(event) => update({ ...cell, radius: Number(event.target.value) })}
              />
            </label>
          )}
        </div>
      </section>
    </div>
  );
}

export function App() {
  const [stage, setStage] = useState<Stage>("setup");
  const [userEmail, setUserEmail] = useState("pilot@example.edu");
  const [settings, setSettings] = useState<AnalysisSettings>(defaultSettings);
  const [files, setFiles] = useState<File[]>([]);
  const [run, setRun] = useState<RunRecord | null>(null);
  const [uploadSessions, setUploadSessions] = useState<UploadSession[]>([]);
  const [uploadProgress, setUploadProgress] = useState<UploadProgress>({});
  const [rois, setRois] = useState<RoiGeometry[]>([
    { role: "bleach", kind: "circle", center_x: 256, center_y: 23, radius: 9 }
  ]);
  const [results, setResults] = useState<OutputManifest | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const totalSize = useMemo(() => files.reduce((sum, file) => sum + file.size, 0), [files]);

  function updateSetting<K extends keyof AnalysisSettings>(key: K, value: AnalysisSettings[K]) {
    const next = { ...settings, [key]: value };
    if (key === "model") {
      next.fit_mode = fitModes(value as AnalysisSettings["model"])[0];
    }
    setSettings(next);
  }

  async function handleCreateRun() {
    setBusy(true);
    setError(null);
    try {
      const created = await createRun(userEmail, settings);
      setRun(created);
      setStage("upload");
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setBusy(false);
    }
  }

  async function handlePrepareUploads() {
    if (!run || files.length === 0) return;
    setBusy(true);
    setError(null);
    try {
      const response = await createUploadSessions(userEmail, run.run_id, files);
      setRun(response.run);
      setUploadSessions(response.upload_sessions);
      setUploadProgress(Object.fromEntries(files.map((file) => [file.name, 0])));
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setBusy(false);
    }
  }

  async function handleUploadFiles() {
    setBusy(true);
    setError(null);
    try {
      for (const file of files) {
        const session = uploadSessions.find((entry) => entry.file_name === file.name);
        if (!session) throw new Error(`Missing upload session for ${file.name}`);
        await uploadFileResumable(file, session, (uploaded) => {
          setUploadProgress((current) => ({ ...current, [file.name]: uploaded / file.size }));
        });
      }
      if (run) {
        const preview = await queuePreview(userEmail, run.run_id);
        setRun(preview.run);
      }
      setStage("roi");
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setBusy(false);
    }
  }

  async function handleSubmit() {
    if (!run) return;
    setBusy(true);
    setError(null);
    try {
      const submitted = await submitRun(userEmail, run.run_id, settings, rois);
      setRun(submitted.run);
      setStage("submit");
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setBusy(false);
    }
  }

  async function refreshRun() {
    if (!run) return;
    setBusy(true);
    setError(null);
    try {
      const nextRun = await getRun(userEmail, run.run_id);
      setRun(nextRun);
      if (nextRun.output_manifest_uri) {
        const output = await getResults(userEmail, run.run_id);
        setResults(output);
        setStage("results");
      }
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught));
    } finally {
      setBusy(false);
    }
  }

  return (
    <main className="app-shell">
      <header className="topbar">
        <div>
          <p className="eyebrow">FRAP Toolbox</p>
          <h1>Cloud Analysis Workspace</h1>
        </div>
        <div className="identity">
          <ShieldCheck size={18} />
          <input
            aria-label="Pilot email"
            value={userEmail}
            onChange={(event) => setUserEmail(event.target.value)}
          />
        </div>
      </header>

      <Stepper stage={stage} />

      {error && <div className="alert">{error}</div>}

      <div className="workspace">
        <section className="panel control-panel">
          <div className="panel-title">
            <FlaskConical size={18} />
            <h2>Analysis Setup</h2>
          </div>
          <div className="field-grid">
            <label>
              Model
              <select value={settings.model} onChange={(event) => updateSetting("model", event.target.value as AnalysisSettings["model"])}>
                <option value="diffusion">Diffusion</option>
                <option value="reaction1">Reaction 1</option>
                <option value="reaction2">Reaction 2</option>
              </select>
            </label>
            <label>
              Fit mode
              <select value={settings.fit_mode} onChange={(event) => updateSetting("fit_mode", event.target.value)}>
                {fitModes(settings.model).map((mode) => (
                  <option key={mode} value={mode}>
                    {mode}
                  </option>
                ))}
              </select>
            </label>
            <label>
              Post-bleach frame
              <input
                type="number"
                min="1"
                value={settings.post_bleach_frame}
                onChange={(event) => updateSetting("post_bleach_frame", Number(event.target.value))}
              />
            </label>
            <label>
              Pre-bleach frames
              <input
                type="number"
                min="1"
                value={settings.pre_bleach_count}
                onChange={(event) => updateSetting("pre_bleach_count", Number(event.target.value))}
              />
            </label>
            <label>
              Background
              <input
                type="number"
                min="0"
                value={settings.background}
                onChange={(event) => updateSetting("background", Number(event.target.value))}
              />
            </label>
          </div>
          <div className="toggles">
            <label className="check">
              <input
                type="checkbox"
                checked={settings.normalize_by_cell}
                onChange={(event) => updateSetting("normalize_by_cell", event.target.checked)}
              />
              Whole-cell normalization
            </label>
            <label className="check">
              <input
                type="checkbox"
                checked={settings.use_adjacent_roi}
                disabled={settings.model !== "diffusion"}
                onChange={(event) => updateSetting("use_adjacent_roi", event.target.checked)}
              />
              Adjacent ROI correction
            </label>
          </div>
          <button className="primary" onClick={handleCreateRun} disabled={busy || !userEmail}>
            {busy ? <Loader2 className="spin" size={18} /> : <Play size={18} />}
            Create run
          </button>
        </section>

        <section className="panel upload-panel">
          <div className="panel-title">
            <CloudUpload size={18} />
            <h2>Data Upload</h2>
          </div>
          <label className="dropzone">
            <FileUp size={22} />
            <span>{files.length ? `${files.length} files selected, ${bytesLabel(totalSize)}` : "Select microscopy files"}</span>
            <input
              type="file"
              multiple
              onChange={(event) => setFiles(Array.from(event.target.files ?? []))}
            />
          </label>
          <div className="file-list">
            {files.map((file) => (
              <div className="file-row" key={file.name}>
                <span>{file.name}</span>
                <strong>{bytesLabel(file.size)}</strong>
                {uploadProgress[file.name] !== undefined && (
                  <progress value={uploadProgress[file.name]} max="1" aria-label={`${file.name} upload progress`} />
                )}
              </div>
            ))}
          </div>
          <div className="button-row">
            <button onClick={handlePrepareUploads} disabled={busy || !run || !files.length}>
              {busy ? <Loader2 className="spin" size={18} /> : <CloudUpload size={18} />}
              Prepare
            </button>
            <button onClick={handleUploadFiles} disabled={busy || !uploadSessions.length}>
              <CloudUpload size={18} />
              Upload
            </button>
          </div>
        </section>
      </div>

      <RoiEditor rois={rois} setRois={setRois} normalizeByCell={settings.normalize_by_cell} />

      <section className="panel run-panel">
        <div className="panel-title">
          <Activity size={18} />
          <h2>Run Status</h2>
        </div>
        <div className="status-grid">
          <div>
            <span>Run</span>
            <strong>{run?.run_id ?? "Not created"}</strong>
          </div>
          <div>
            <span>Status</span>
            <strong>{run?.status ?? "Idle"}</strong>
          </div>
          <div>
            <span>Job</span>
            <strong>{run?.analysis_job_id ?? "Not submitted"}</strong>
          </div>
        </div>
        <div className="button-row">
          <button className="primary" onClick={handleSubmit} disabled={busy || !run}>
            {busy ? <Loader2 className="spin" size={18} /> : <Play size={18} />}
            Submit analysis
          </button>
          <button onClick={refreshRun} disabled={busy || !run}>
            <RefreshCw size={18} />
            Refresh
          </button>
        </div>
      </section>

      <section className="panel results-panel">
        <div className="panel-title">
          <CheckCircle2 size={18} />
          <h2>Results</h2>
        </div>
        {results ? (
          <div className="results-grid">
            {Object.entries(results.parameters).map(([name, value]) => (
              <div className="metric" key={name}>
                <span>{name}</span>
                <strong>{value === null ? "n/a" : value.toPrecision(4)}</strong>
              </div>
            ))}
          </div>
        ) : (
          <p className="empty">Results will appear when the output manifest is available.</p>
        )}
      </section>
    </main>
  );
}
