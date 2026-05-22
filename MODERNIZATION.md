# FRAP-Toolbox modernization plan

## Direction

The Python port should make the analysis engine the durable center of the project:

1. Core library: image loading, ROI masks, normalization, fitting, exports.
2. CLI: reproducible batch analysis and regression checks.
3. Local browser app: low-friction deployment for lab workstations.
4. MATLAB parity tests: keep historical behavior explainable while the Python implementation improves.

## Technology choices

- Image I/O: BioIO as the preferred reader layer, with `bioio-tifffile` installed by default.
- Format plugins: add `bioio-nd2` or `bioio-bioformats` only when needed.
- GUI: Streamlit for the first modern app surface because it avoids native desktop packaging and runs locally with one command.
- Legacy support: AICSImageIO remains as an optional fallback while the port transitions.

## Near-term milestones

1. Keep diffusion CLI and app numerically aligned with MATLAB reference outputs.
2. Port Reaction 1 and Reaction 2 into `frap_toolbox_py.models`.
3. Add export writers that reproduce the MATLAB tab-delimited output files.
4. Add interactive ROI drawing or image-preview workflows.
5. Package test data intentionally, likely with Git LFS or downloadable fixtures rather than normal git blobs.

## Local-only assets

- `test-data/` is intentionally ignored for now because the microscopy fixtures are large.
- `scratch/` is intentionally ignored and can hold one-off parity scripts, debug probes, and generated fit artifacts.
- The current MATLAB-vs-Python diffusion fitting investigation is preserved in `docs/matlab-to-python-port-testing.md`.
