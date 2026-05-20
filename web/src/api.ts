export type ModelKey = "diffusion" | "reaction1" | "reaction2";

export interface AnalysisSettings {
  model: ModelKey;
  post_bleach_frame: number;
  pre_bleach_count: number;
  background: number;
  normalize_by_cell: boolean;
  fit_mode?: string;
  use_adjacent_roi: boolean;
  adjacent_offset: number;
  max_profile_radius?: number;
}

export interface RoiGeometry {
  role: "bleach" | "cell" | "adjacent";
  kind: "circle" | "polygon" | "freehand" | "mask";
  center_x?: number;
  center_y?: number;
  radius?: number;
  points?: Array<[number, number]>;
  mask_uri?: string;
  source_file?: string;
  notes?: string;
}

export interface RunRecord {
  run_id: string;
  user_id: string;
  status: string;
  settings: AnalysisSettings;
  files: Array<Record<string, unknown>>;
  rois: RoiGeometry[];
  manifest_uri?: string;
  output_manifest_uri?: string;
  preview_job_id?: string;
  analysis_job_id?: string;
  error?: string;
  created_at: string;
  updated_at: string;
}

export interface UploadSession {
  file_name: string;
  object_uri: string;
  upload_url: string;
  upload_method: "PUT";
  headers: Record<string, string>;
  expires_at?: string;
}

export interface OutputManifest {
  run_id: string;
  status: string;
  parameters: Record<string, number | null>;
  result_tables: Record<string, string>;
  plots: Record<string, string>;
  roi_masks: Record<string, string>;
  warnings: string[];
  errors: string[];
}

const API_BASE = import.meta.env.VITE_FRAP_API_BASE ?? "/api";

function apiUrl(pathOrUrl: string): string {
  if (/^https?:\/\//.test(pathOrUrl)) return pathOrUrl;
  if (pathOrUrl.startsWith("/")) return `${API_BASE.replace(/\/$/, "")}${pathOrUrl}`;
  return pathOrUrl;
}

async function request<T>(path: string, userEmail: string, init: RequestInit = {}): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`, {
    ...init,
    headers: {
      "content-type": "application/json",
      "x-frap-user": userEmail,
      ...(init.headers ?? {})
    }
  });
  if (!response.ok) {
    const detail = await response.text();
    throw new Error(detail || `${response.status} ${response.statusText}`);
  }
  return response.json() as Promise<T>;
}

export async function createRun(userEmail: string, settings: AnalysisSettings): Promise<RunRecord> {
  return request<RunRecord>("/runs", userEmail, {
    method: "POST",
    body: JSON.stringify({ settings })
  });
}

export async function createUploadSessions(
  userEmail: string,
  runId: string,
  files: File[]
): Promise<{ run: RunRecord; upload_sessions: UploadSession[] }> {
  return request<{ run: RunRecord; upload_sessions: UploadSession[] }>(`/runs/${runId}/uploads`, userEmail, {
    method: "POST",
    body: JSON.stringify({
      files: files.map((file) => ({
        file_name: file.name,
        size_bytes: file.size,
        content_type: file.type || "application/octet-stream"
      }))
    })
  });
}

export async function queuePreview(userEmail: string, runId: string): Promise<{ run: RunRecord }> {
  return request<{ run: RunRecord }>(`/runs/${runId}/preview`, userEmail, {
    method: "POST",
    body: JSON.stringify({})
  });
}

export async function submitRun(
  userEmail: string,
  runId: string,
  settings: AnalysisSettings,
  rois: RoiGeometry[]
): Promise<{ run: RunRecord }> {
  return request<{ run: RunRecord }>(`/runs/${runId}/submit`, userEmail, {
    method: "POST",
    body: JSON.stringify({ settings, rois })
  });
}

export async function getRun(userEmail: string, runId: string): Promise<RunRecord> {
  return request<RunRecord>(`/runs/${runId}`, userEmail);
}

export async function getResults(userEmail: string, runId: string): Promise<OutputManifest> {
  return request<OutputManifest>(`/runs/${runId}/results`, userEmail);
}

export async function uploadFileResumable(
  file: File,
  session: UploadSession,
  onProgress: (uploaded: number) => void,
  signal?: AbortSignal
): Promise<void> {
  const chunkSize = 8 * 1024 * 1024;
  let offset = 0;

  while (offset < file.size) {
    const endExclusive = Math.min(offset + chunkSize, file.size);
    const chunk = file.slice(offset, endExclusive);
    const headers: Record<string, string> = {
      ...session.headers,
      "content-range": `bytes ${offset}-${endExclusive - 1}/${file.size}`
    };
    const response = await fetch(apiUrl(session.upload_url), {
      method: session.upload_method,
      headers,
      body: chunk,
      signal
    });

    if (response.status === 308) {
      const range = response.headers.get("range");
      offset = range ? Number(range.split("-").pop()) + 1 : endExclusive;
      onProgress(offset);
      continue;
    }

    if (!response.ok) {
      throw new Error(`Upload failed for ${file.name}: ${response.status} ${response.statusText}`);
    }

    offset = endExclusive;
    onProgress(offset);
  }
}
