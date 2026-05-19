# Modernization Workstreams

This file records the active modernization split after the parity checkpoint at
commit `92355f5`.

## Integration Branch

- Branch: `codex/modernization-cleanup`
- Worktree: `/Users/lewiskraft/Coding/FRAP-Toolbox`
- Role: checkpoint branch and integration point for reviewed workstream branches.

## Workstream Branches

| Branch | Worktree | Primary Goal | Integration Dependency |
| --- | --- | --- | --- |
| `codex/roi-mask-io` | `/Users/lewiskraft/Coding/FRAP-Toolbox-roi-mask-io` | Saved ROI/mask format, validation, and CLI mask inputs | First |
| `codex/reaction-raw-parity` | `/Users/lewiskraft/Coding/FRAP-Toolbox-reaction-raw-parity` | Quantify and tighten raw Reaction 1/2 ND2 parity | After mask IO where possible |
| `codex/gui-deployment` | `/Users/lewiskraft/Coding/FRAP-Toolbox-gui-deployment` | Practical Streamlit workflow for diffusion/reaction analysis and ROI use | After or alongside mask IO |
| `codex/packaging-ci` | `/Users/lewiskraft/Coding/FRAP-Toolbox-packaging-ci` | Install extras, CI, public-facing docs, and repo hygiene | Can merge independently if clean |

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
- The raw reaction workflows still need the hand-drawn analysis bleach and
  whole-cell ROI masks that MATLAB collected interactively via `roipoly`.

## Merge Order

1. Merge `codex/roi-mask-io` first so all other workstreams can consume the same
   saved-mask contract.
2. Merge `codex/packaging-ci` once it is clean, unless it changes installation
   docs that depend on the mask format.
3. Merge `codex/reaction-raw-parity` after reconciling any mask helper names or
   test fixtures with mask IO.
4. Merge `codex/gui-deployment` after the app can load the finalized mask
   format or has a clearly documented interim path.
