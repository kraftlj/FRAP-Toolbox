from __future__ import annotations

MODEL_KEYS = ("diffusion", "reaction1", "reaction2")
MODEL_LABELS = {
    "diffusion": "Diffusion",
    "reaction1": "Reaction 1",
    "reaction2": "Reaction 2",
}
MODEL_INDEX = {
    "diffusion": 1,
    "reaction1": 2,
    "reaction2": 3,
}
DIFFUSION_FIT_MODES = (
    "global",
    "individual",
    "average_curve",
    "simplified_kang",
    "simplified_kang_global",
)
REACTION_FIT_MODES = ("individual", "average_curve")
OPTIMIZER_MODES = ("modern", "legacy_matlab")
DEFAULT_FIT_MODES = {
    "diffusion": "global",
    "reaction1": "individual",
    "reaction2": "average_curve",
}


def default_fit_mode(model_key: str) -> str:
    return DEFAULT_FIT_MODES.get(model_key, DEFAULT_FIT_MODES["diffusion"])
