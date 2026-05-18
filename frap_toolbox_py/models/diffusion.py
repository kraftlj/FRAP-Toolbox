from __future__ import annotations

import math
from typing import Iterable, List, Tuple

import numpy as np
from scipy.optimize import curve_fit, least_squares

from ..photodecay import estimate_decay_rate
from ..types import (
    BasicInputs,
    DiffusionFitConfig,
    DiffusionFitResult,
    FitBounds,
    FRAPDataset,
)


def _is_adjustable(bounds: FitBounds) -> bool:
    return bounds.mode.lower() != "fixed"


def _slice_range(range_tuple: Tuple[int, int]) -> slice:
    start, end = range_tuple
    return slice(start, end)


def kang_frap(
    t: np.ndarray,
    r_effective: float,
    r_nominal: float,
    diffusion: float,
    bleach_depth: float,
    terms: int = 10,
) -> np.ndarray:
    m = np.arange(0, terms + 1)[:, None]
    factorials = np.array([math.factorial(i) for i in range(terms + 1)], dtype=float)[:, None]
    numerator = ((-bleach_depth) ** m) * (r_effective**2)
    numerator /= factorials
    denominator = r_effective**2 + m * (8.0 * diffusion * t[None, :] + r_nominal**2)
    partial = numerator / denominator
    return partial.sum(axis=0)


def _pbp_model(r: np.ndarray, bleach_depth: float, r_effective: float) -> np.ndarray:
    return np.exp(-bleach_depth * np.exp(-2.0 * r**2 / r_effective**2))


def _fit_curve(
    model,
    xdata: np.ndarray,
    ydata: np.ndarray,
    initial: List[float],
    lower: List[float],
    upper: List[float],
    legacy_matlab: bool,
) -> np.ndarray:
    if legacy_matlab:
        result = least_squares(
            lambda params: model(xdata, *params) - ydata,
            x0=initial,
            bounds=(lower, upper),
            ftol=1e-6,
            xtol=1e-6,
            gtol=1e-6,
            max_nfev=max(200, 100 * len(initial)),
        )
        return result.x

    popt, _ = curve_fit(model, xdata, ydata, p0=initial, bounds=(lower, upper))
    return popt


def _prepare_average_segments(
    datasets: List[FRAPDataset],
    frap_slice: slice,
    profile_slice: slice,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ref_length = len(datasets[0].corrected_frap[frap_slice])
    ref_time = datasets[0].time[frap_slice]

    frap_matrix = []
    for dataset in datasets:
        segment = dataset.corrected_frap[frap_slice]
        if len(segment) != ref_length:
            raise ValueError("FRAP segments must have consistent length for averaging.")
        time_segment = dataset.time[frap_slice]
        if not np.allclose(time_segment, ref_time, rtol=1e-3, atol=1e-6):
            if not np.all(np.isfinite(time_segment)):
                raise ValueError("Encountered non-finite timestamps; cannot interpolate.")
            adjusted_time = np.maximum.accumulate(time_segment.astype(float, copy=True))
            diffs = np.diff(adjusted_time)
            if np.any(diffs <= 0):
                eps = np.finfo(adjusted_time.dtype).eps
                for idx in range(1, adjusted_time.size):
                    if adjusted_time[idx] <= adjusted_time[idx - 1]:
                        adjusted_time[idx] = adjusted_time[idx - 1] + eps
            segment = np.interp(ref_time, adjusted_time, segment)
        frap_matrix.append(segment)
    averaged_frap = np.asarray(frap_matrix).mean(axis=0)

    radius_reference = datasets[0].radius[profile_slice]
    profile_matrix = []
    for dataset in datasets:
        r = dataset.radius[profile_slice]
        pbp = dataset.post_bleach_profile[profile_slice]
        interp = np.interp(radius_reference, r, pbp)
        profile_matrix.append(interp)
    averaged_profile = np.asarray(profile_matrix).mean(axis=0)

    return ref_time, averaged_frap, radius_reference, averaged_profile


def _fit_diffusion_components(
    fit_time: np.ndarray,
    corrected_frap: np.ndarray,
    radius: np.ndarray,
    post_bleach_profile: np.ndarray,
    nominal_radius: float,
    config: DiffusionFitConfig,
) -> dict:
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    k_bounds = config.k
    r_bounds = config.r_effective

    if _is_adjustable(k_bounds) and _is_adjustable(r_bounds):
        bleach_depth, r_effective = _fit_curve(
            _pbp_model,
            radius,
            post_bleach_profile,
            [k_bounds.initial, r_bounds.initial],
            [k_bounds.lower, r_bounds.lower],
            [k_bounds.upper, r_bounds.upper],
            legacy_matlab,
        )
    elif _is_adjustable(k_bounds):
        r_effective = r_bounds.initial
        bleach_depth = _fit_curve(
            lambda r, k: _pbp_model(r, k, r_effective),
            radius,
            post_bleach_profile,
            [k_bounds.initial],
            [k_bounds.lower],
            [k_bounds.upper],
            legacy_matlab,
        )[0]
    elif _is_adjustable(r_bounds):
        bleach_depth = k_bounds.initial
        r_effective = _fit_curve(
            lambda r, re: _pbp_model(r, bleach_depth, re),
            radius,
            post_bleach_profile,
            [r_bounds.initial],
            [r_bounds.lower],
            [r_bounds.upper],
            legacy_matlab,
        )[0]
    else:
        bleach_depth = k_bounds.initial
        r_effective = r_bounds.initial

    frap_initial = corrected_frap[0]
    if _is_adjustable(k_bounds):
        bleach_depth = _fit_curve(
            lambda t, k: kang_frap(t, r_effective, nominal_radius, 1.0, k),
            np.array([0.0]),
            np.array([frap_initial]),
            [bleach_depth],
            [k_bounds.lower],
            [k_bounds.upper],
            legacy_matlab,
        )[0]
        if _is_adjustable(r_bounds):
            r_effective = _fit_curve(
                lambda r, re: _pbp_model(r, bleach_depth, re),
                radius,
                post_bleach_profile,
                [r_effective],
                [r_bounds.lower],
                [r_bounds.upper],
                legacy_matlab,
            )[0]

    D_bounds = config.diffusion_coefficient
    MF_bounds = config.mobile_fraction
    weights = corrected_frap.sum()
    weighted_target = corrected_frap / (fit_time + weights)

    def weighted_model(t, diffusion, mobile_fraction):
        numerator = kang_frap(t, r_effective, nominal_radius, diffusion, bleach_depth) * mobile_fraction
        numerator += (1.0 - mobile_fraction) * frap_initial
        return numerator / (t + weights)

    if _is_adjustable(D_bounds) and _is_adjustable(MF_bounds):
        diffusion, mobile_fraction = _fit_curve(
            weighted_model,
            fit_time,
            weighted_target,
            [D_bounds.initial, MF_bounds.initial],
            [D_bounds.lower, MF_bounds.lower],
            [D_bounds.upper, MF_bounds.upper],
            legacy_matlab,
        )
    elif _is_adjustable(D_bounds):
        mobile_fraction = MF_bounds.initial
        diffusion = _fit_curve(
            lambda t, diffusion: weighted_model(t, diffusion, mobile_fraction),
            fit_time,
            weighted_target,
            [D_bounds.initial],
            [D_bounds.lower],
            [D_bounds.upper],
            legacy_matlab,
        )[0]
    elif _is_adjustable(MF_bounds):
        diffusion = D_bounds.initial
        mobile_fraction = _fit_curve(
            lambda t, mobile_fraction: weighted_model(t, diffusion, mobile_fraction),
            fit_time,
            weighted_target,
            [MF_bounds.initial],
            [MF_bounds.lower],
            [MF_bounds.upper],
            legacy_matlab,
        )[0]
    else:
        diffusion = D_bounds.initial
        mobile_fraction = MF_bounds.initial

    weighted_prediction = weighted_model(fit_time, diffusion, mobile_fraction)
    profile_prediction = _pbp_model(radius, bleach_depth, r_effective)
    return {
        "k": float(bleach_depth),
        "r_effective": float(r_effective),
        "diffusion_coefficient": float(diffusion),
        "mobile_fraction": float(mobile_fraction),
        "sum_squared_residuals": float(np.sum((weighted_target - weighted_prediction) ** 2)),
        "post_bleach_fit": profile_prediction,
        "post_bleach_residuals": post_bleach_profile - profile_prediction,
        "frap_fit": weighted_prediction * (fit_time + weights),
        "frap_residuals": weighted_target - weighted_prediction,
    }


def _evaluate_diffusion_components(
    fit_time: np.ndarray,
    corrected_frap: np.ndarray,
    radius: np.ndarray,
    post_bleach_profile: np.ndarray,
    nominal_radius: float,
    k: float,
    r_effective: float,
    diffusion: float,
    mobile_fraction: float,
) -> dict:
    frap_initial = corrected_frap[0]
    weights = corrected_frap.sum()
    weighted_target = corrected_frap / (fit_time + weights)
    weighted_prediction = (
        kang_frap(fit_time, r_effective, nominal_radius, diffusion, k) * mobile_fraction
        + (1.0 - mobile_fraction) * frap_initial
    ) / (fit_time + weights)
    profile_prediction = _pbp_model(radius, k, r_effective)
    return {
        "sum_squared_residuals": float(np.sum((weighted_target - weighted_prediction) ** 2)),
        "post_bleach_fit": profile_prediction,
        "post_bleach_residuals": post_bleach_profile - profile_prediction,
        "frap_fit": weighted_prediction * (fit_time + weights),
        "frap_residuals": weighted_target - weighted_prediction,
    }


def fit_diffusion_model(
    datasets: List[FRAPDataset],
    inputs: BasicInputs,
    config: DiffusionFitConfig,
) -> DiffusionFitResult:
    if not datasets:
        raise ValueError("At least one dataset is required for fitting.")

    roi_values = list(inputs.roi_definition)
    if len(roi_values) < 3:
        raise ValueError("ROI definition must include at least three values (x, y, radius).")

    voxel_size = next((d.voxel_size_x for d in datasets if d.voxel_size_x is not None), 1.0)
    nominal_radius = float(roi_values[2]) * (voxel_size or 1.0)
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"

    decay_bounds = config.decay_rate
    if _is_adjustable(decay_bounds):
        decay_rate = estimate_decay_rate(
            datasets,
            inputs.post_bleach_frame,
            config.decay_fit_range,
            decay_bounds.initial,
            decay_bounds.lower,
            decay_bounds.upper,
            legacy_matlab=legacy_matlab,
        )
    else:
        decay_rate = decay_bounds.initial

    mf_slice = _slice_range(config.mobile_fraction_range)
    for dataset in datasets:
        original_time = dataset.time.astype(float, copy=True)
        time_from_first = original_time - original_time[0]
        time_from_post = original_time - original_time[inputs.post_bleach_frame]
        correction = np.exp(-decay_rate * time_from_first)
        dataset.corrected_frap = dataset.norm_frap / correction
        dataset.time = time_from_post
        if dataset.adjacent is not None and inputs.use_adjacent_roi:
            dataset.corrected_mobile_fraction = float(
                1.0 - (np.mean(dataset.adjacent[mf_slice]) - np.mean(dataset.norm_frap[mf_slice]))
            )

    frap_slice = _slice_range(config.frap_range)
    profile_slice = _slice_range(config.profile_range)
    ref_time, averaged_frap, radius_reference, averaged_profile = _prepare_average_segments(
        datasets,
        frap_slice,
        profile_slice,
    )

    if not config.fit_averaged_data:
        for dataset in datasets:
            fit = _fit_diffusion_components(
                dataset.time[frap_slice],
                dataset.corrected_frap[frap_slice],
                dataset.radius[profile_slice],
                dataset.post_bleach_profile[profile_slice],
                nominal_radius,
                config,
            )
            dataset.k = fit["k"]
            dataset.r_effective = fit["r_effective"]
            dataset.diffusion_coefficient = fit["diffusion_coefficient"]
            dataset.mobile_fraction = fit["mobile_fraction"]
            dataset.post_bleach_fit = fit["post_bleach_fit"]
            dataset.post_bleach_residuals = fit["post_bleach_residuals"]
            dataset.fit_time = dataset.time[frap_slice]
            dataset.frap_fit = fit["frap_fit"]
            dataset.frap_residuals = fit["frap_residuals"]
            dataset.sum_squared_residuals = fit["sum_squared_residuals"]

        k = float(np.nanmean([dataset.k for dataset in datasets]))
        r_effective = float(np.nanmean([dataset.r_effective for dataset in datasets]))
        diffusion = float(np.nanmean([dataset.diffusion_coefficient for dataset in datasets]))
        mobile_fraction = float(np.nanmean([dataset.mobile_fraction for dataset in datasets]))
        aggregate = _evaluate_diffusion_components(
            ref_time,
            averaged_frap,
            radius_reference,
            averaged_profile,
            nominal_radius,
            k,
            r_effective,
            diffusion,
            mobile_fraction,
        )
        sum_squared_residuals = aggregate["sum_squared_residuals"]
        profile_prediction = aggregate["post_bleach_fit"]
        profile_residuals = aggregate["post_bleach_residuals"]
        frap_fit = aggregate["frap_fit"]
        frap_residuals = aggregate["frap_residuals"]
    else:
        fit = _fit_diffusion_components(
            ref_time,
            averaged_frap,
            radius_reference,
            averaged_profile,
            nominal_radius,
            config,
        )
        k = fit["k"]
        r_effective = fit["r_effective"]
        diffusion = fit["diffusion_coefficient"]
        mobile_fraction = fit["mobile_fraction"]
        sum_squared_residuals = fit["sum_squared_residuals"]
        profile_prediction = fit["post_bleach_fit"]
        profile_residuals = fit["post_bleach_residuals"]
        frap_fit = fit["frap_fit"]
        frap_residuals = fit["frap_residuals"]

    valid_cmf = [d.corrected_mobile_fraction for d in datasets if d.corrected_mobile_fraction is not None]
    corrected_mobile_fraction = float(np.nanmean(valid_cmf)) if valid_cmf else None

    return DiffusionFitResult(
        datasets=datasets,
        k=k,
        r_effective=r_effective,
        diffusion_coefficient=diffusion,
        mobile_fraction=mobile_fraction,
        corrected_mobile_fraction=corrected_mobile_fraction,
        sum_squared_residuals=sum_squared_residuals,
        decay_rate=float(decay_rate),
        averaged_profile_radius=radius_reference,
        averaged_profile=averaged_profile,
        averaged_profile_fit=profile_prediction,
        averaged_profile_residuals=profile_residuals,
        averaged_time=ref_time,
        averaged_frap=averaged_frap,
        averaged_frap_fit=frap_fit,
        averaged_frap_residuals=frap_residuals,
    )


def default_diffusion_config(
    frap_length: int,
    profile_length: int,
    post_bleach_frame: int,
) -> DiffusionFitConfig:
    post_start = min(max(post_bleach_frame, 0), max(frap_length - 1, 0))
    decay_fit_start = max(frap_length - 50, post_start)
    corrected_mf_start = max(int(frap_length * 0.9), post_start)
    return DiffusionFitConfig(
        k=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        r_effective=FitBounds(3.0, 0.0, np.inf, "Adjustable"),
        diffusion_coefficient=FitBounds(10.0, 0.0, np.inf, "Adjustable"),
        mobile_fraction=FitBounds(1.0, 0.0, 2.0, "Adjustable"),
        decay_rate=FitBounds(0.001, 0.0, np.inf, "Adjustable"),
        profile_range=(0, profile_length),
        frap_range=(post_start, frap_length),
        decay_fit_range=(decay_fit_start, frap_length),
        mobile_fraction_range=(corrected_mf_start, frap_length),
        fit_averaged_data=True,
    )
