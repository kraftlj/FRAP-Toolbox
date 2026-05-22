# Modernization Workstreams

This file records the modernization split after the parity checkpoint at commit
`92355f5` and its integration status.

## Integration Branch

- Branch: `codex/modernization-cleanup`
- Worktree: `/Users/lewiskraft/Coding/FRAP-Toolbox`
- Role: checkpoint branch and integration point for reviewed workstream branches.

## Workstream Branches

| Branch | Worktree | Primary Goal | Integration Status |
| --- | --- | --- | --- |
| `codex/roi-mask-io` | `/Users/lewiskraft/Coding/FRAP-Toolbox-roi-mask-io` | Saved ROI/mask format, validation, and CLI mask inputs | Merged in `5a88557` |
| `codex/packaging-ci` | `/Users/lewiskraft/Coding/FRAP-Toolbox-packaging-ci` | Install extras, CI, public-facing docs, and repo hygiene | Merged in `0b19021` |
| `codex/gui-deployment` | `/Users/lewiskraft/Coding/FRAP-Toolbox-gui-deployment` | Practical Streamlit workflow for diffusion/reaction analysis and ROI use | Merged in `3ccb185` |
| `codex/reaction-raw-parity` | `/Users/lewiskraft/Coding/FRAP-Toolbox-reaction-raw-parity` | Quantify and tighten raw Reaction 1/2 ND2 parity | Merged in `0099fe1` |

## Acceptance Gates

Each workstream should finish with:

- `git diff --check` clean.
- `.venv/bin/python -m pytest -q` passing, with only documented expected skips or
  xfails.
- Clear final notes listing changed paths, test results, and any integration
  caveats.

## Current Technical Facts

- Reaction 1 model parity is achieved from MATLAB-exported guide vectors.
- Reaction 2 model support exists, but the archive still lacks a Reaction 2 FRAP
  vector export.
- Raw ND2 timestamp loading works for `Venus_1002.nd2`,
  `Venus-Atg5_1002.nd2`, and `Venus-Atg5_1003.nd2`.
- `Venus_1001.nd2` currently fails in Python ND2 readers with a legacy metadata
  header parse error.
- Saved `bleach` and `cell` masks now have a canonical `.npz` container and can
  be loaded through the CLI and app for reaction workflows.
- The raw reaction workflows still need the actual hand-drawn analysis bleach
  and whole-cell ROI masks that MATLAB collected interactively via `roipoly`.
- Embedded ND2 stimulation ROI and threshold-derived candidate masks for
  `Venus_1002.nd2` are documented diagnostics, but neither is exact raw-curve
  parity.

## Completed Merge Order

1. `codex/roi-mask-io`
2. `codex/packaging-ci`
3. `codex/gui-deployment`
4. `codex/reaction-raw-parity`

## Integrated Verification

- Local parity environment with ignored `test-data/`: `44 passed, 1 xfailed`.
- Clean-checkout-style run from outside the repository: `22 passed, 23 skipped`.
