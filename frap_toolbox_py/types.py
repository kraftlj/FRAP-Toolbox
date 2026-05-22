from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional
import numpy as np


@dataclass
class FRAPDataset:
    """Container for a single FRAP dataset and derived quantities."""

    name: str
    time: np.ndarray
    frap: np.ndarray
    norm_frap: np.ndarray
    corrected_frap: np.ndarray
    cell: Optional[np.ndarray] = None
    adjacent: Optional[np.ndarray] = None
    corrected_mobile_fraction: Optional[float] = None
    radius: Optional[np.ndarray] = None
    post_bleach_profile: Optional[np.ndarray] = None
    post_bleach_fit: Optional[np.ndarray] = None
    post_bleach_residuals: Optional[np.ndarray] = None
    fit_time: Optional[np.ndarray] = None
    frap_fit: Optional[np.ndarray] = None
    frap_residuals: Optional[np.ndarray] = None
    sum_squared_residuals: Optional[float] = None
    k: Optional[float] = None
    r_effective: Optional[float] = None
    half_time: Optional[float] = None
    diffusion_coefficient: Optional[float] = None
    mobile_fraction: Optional[float] = None
    reaction_parameters: dict = field(default_factory=dict)
    voxel_size_x: Optional[float] = None
    voxel_size_y: Optional[float] = None
    metadata: dict = field(default_factory=dict)


@dataclass
class BasicInputs:
    """Inputs gathered from the main GUI that control downstream analysis."""

    file_paths: List[Path]
    model_index: int
    roi_mode: int  # 1 = Circle, 2 = User defined
    normalize_by_cell: bool
    background_intensity: float
    post_bleach_frame: int  # Zero-based index of the first post-bleach frame
    roi_definition: Iterable[float]
    pre_bleach_frame_count: int
    use_adjacent_roi: bool


@dataclass
class FitBounds:
    """Bounds and initial guess for a single fit parameter."""

    initial: float
    lower: float
    upper: float
    mode: str  # "Fixed" or "Adjustable"


@dataclass
class DiffusionFitConfig:
    """Parameter configuration for the diffusion model fitting."""

    k: FitBounds
    r_effective: FitBounds
    diffusion_coefficient: FitBounds
    mobile_fraction: FitBounds
    decay_rate: FitBounds
    profile_range: tuple[int, int]
    frap_range: tuple[int, int]
    decay_fit_range: tuple[int, int]
    mobile_fraction_range: tuple[int, int]
    fit_averaged_data: bool
    optimizer_mode: str = "modern"
    # "auto", "individual", "average_curve", "global", "simplified_kang", or "simplified_kang_global"
    fit_mode: str = "auto"


@dataclass
class DiffusionFitResult:
    """Aggregated outputs from a diffusion model fit."""

    datasets: List[FRAPDataset]
    k: float
    r_effective: float
    half_time: Optional[float]
    diffusion_coefficient: float
    mobile_fraction: float
    corrected_mobile_fraction: Optional[float]
    sum_squared_residuals: float
    decay_rate: float
    averaged_profile_radius: Optional[np.ndarray] = None
    averaged_profile: Optional[np.ndarray] = None
    averaged_profile_fit: Optional[np.ndarray] = None
    averaged_profile_residuals: Optional[np.ndarray] = None
    averaged_time: Optional[np.ndarray] = None
    averaged_frap: Optional[np.ndarray] = None
    averaged_frap_fit: Optional[np.ndarray] = None
    averaged_frap_residuals: Optional[np.ndarray] = None


@dataclass
class ReactionFitConfig:
    """Parameter configuration for reaction model fitting."""

    a: FitBounds
    b: FitBounds
    c: FitBounds
    decay_rate: FitBounds
    frap_range: tuple[int, int]
    decay_fit_range: tuple[int, int]
    fit_averaged_data: bool
    d: Optional[FitBounds] = None
    f: Optional[FitBounds] = None
    optimizer_mode: str = "modern"


@dataclass
class ReactionFitResult:
    """Aggregated outputs from a reaction model fit."""

    model_order: int
    datasets: List[FRAPDataset]
    parameters: dict[str, float]
    decay_rate: float
    sum_squared_residuals: float
    averaged_time: Optional[np.ndarray] = None
    averaged_frap: Optional[np.ndarray] = None
    averaged_frap_fit: Optional[np.ndarray] = None
    averaged_frap_residuals: Optional[np.ndarray] = None
