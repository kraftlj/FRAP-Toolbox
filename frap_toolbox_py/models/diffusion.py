from __future__ import annotations

import math
import warnings
from typing import List, Tuple

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit, least_squares

from ..optimization import legacy_matlab_trf
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


def _resolve_fit_mode(config: DiffusionFitConfig) -> str:
    mode = config.fit_mode.lower()
    if mode == "auto":
        return "average_curve" if config.fit_averaged_data else "individual"
    if mode in {"global", "true_global", "global_shared_dm"}:
        return "global"
    if mode in {"simplified", "simplified_kang", "kang_simplified", "half_time"}:
        return "simplified_kang"
    if mode in {"simplified_global", "simplified_kang_global", "global_simplified_kang"}:
        return "simplified_kang_global"
    if mode in {"average", "averaged", "average_curve", "averaged_curve"}:
        return "average_curve"
    if mode in {"individual", "per_dataset"}:
        return "individual"
    raise ValueError(f"Unsupported diffusion fit mode: {config.fit_mode!r}")


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


def simplified_kang_diffusion_coefficient(
    r_nominal: float,
    r_effective: float,
    half_time: float,
) -> float:
    """Return the Kang et al. confocal FRAP half-time diffusion estimate."""

    if not np.isfinite(half_time) or half_time <= 0.0:
        raise ValueError("Half time must be finite and greater than zero.")
    return float((r_nominal**2 + r_effective**2) / (8.0 * half_time))


def simplified_kang_recovery(
    t: np.ndarray,
    r_effective: float,
    r_nominal: float,
    diffusion: float,
    mobile_fraction: float,
    initial: float,
    bleach_depth: float | None = None,
    prebleach: float = 1.0,
) -> np.ndarray:
    """Evaluate the simplified Kang et al. confocal FRAP recovery curve."""

    t = np.asarray(t, dtype=float)
    if r_effective <= 0.0 or diffusion <= 0.0:
        return np.full_like(t, np.nan, dtype=float)

    gamma_squared = (r_nominal / r_effective) ** 2
    tau_d = r_effective**2 / (4.0 * diffusion)
    if bleach_depth is None:
        bleach_depth = (1.0 + gamma_squared) * (1.0 - initial / prebleach)
    mobile_pool = prebleach * (
        1.0 - bleach_depth / (1.0 + gamma_squared + 2.0 * t / tau_d)
    )
    return mobile_pool * mobile_fraction + (1.0 - mobile_fraction) * initial


def estimate_recovery_half_time(
    time: np.ndarray,
    recovery: np.ndarray,
    steady_state: float | None = None,
) -> float:
    """Estimate recovery half time by linear interpolation through the midpoint."""

    time = np.asarray(time, dtype=float)
    recovery = np.asarray(recovery, dtype=float)
    mask = np.isfinite(time) & np.isfinite(recovery)
    time = time[mask]
    recovery = recovery[mask]
    if time.size < 2:
        raise ValueError("At least two finite recovery samples are required.")

    order = np.argsort(time)
    time = time[order]
    recovery = recovery[order]
    start = recovery[0]
    end = recovery[-1] if steady_state is None else float(steady_state)
    if not np.isfinite(end):
        raise ValueError("Steady-state recovery estimate must be finite.")

    amplitude = end - start
    if np.isclose(amplitude, 0.0):
        raise ValueError("Recovery amplitude is too small to estimate a half time.")

    target = start + 0.5 * amplitude
    signed_distance = np.sign(amplitude) * (recovery - target)
    crossing_indices = np.where((signed_distance >= 0.0) & (time > time[0]))[0]
    if crossing_indices.size == 0:
        raise ValueError("Recovery curve does not cross its half-recovery level.")

    idx = int(crossing_indices[0])
    prev_idx = idx - 1
    t0, t1 = time[prev_idx], time[idx]
    y0, y1 = recovery[prev_idx], recovery[idx]
    if np.isclose(y1, y0):
        return float(t1)

    fraction = (target - y0) / (y1 - y0)
    return float(t0 + fraction * (t1 - t0))


def _pbp_model(r: np.ndarray, bleach_depth: float, r_effective: float) -> np.ndarray:
    return np.exp(-bleach_depth * np.exp(-2.0 * r**2 / r_effective**2))


def _simplified_pbp_model(r: np.ndarray, bleach_depth: float, r_effective: float) -> np.ndarray:
    return 1.0 - bleach_depth * np.exp(-2.0 * r**2 / r_effective**2)


def _fit_curve(
    model,
    xdata: np.ndarray,
    ydata: np.ndarray,
    initial: List[float],
    lower: List[float],
    upper: List[float],
    legacy_matlab: bool,
    legacy_x_scale: List[float] | None = None,
) -> np.ndarray:
    if legacy_matlab:
        def residual(params):
            return model(xdata, *params) - ydata

        if legacy_x_scale is not None:
            return legacy_matlab_trf(
                residual,
                np.asarray(initial, dtype=float),
                np.asarray(lower, dtype=float),
                np.asarray(upper, dtype=float),
                np.asarray(legacy_x_scale, dtype=float),
            )

        result = least_squares(
            residual,
            x0=initial,
            bounds=(lower, upper),
            ftol=1e-6,
            xtol=1e-6,
            gtol=1e-6,
            max_nfev=max(200, 100 * len(initial)),
        )
        return result.x

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=OptimizeWarning)
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
        legacy_x_scale = None
        if legacy_matlab:
            # Empirically matches the MATLAB R2013-era D/MF trust-region path
            # against the archived guide exports while keeping modern fits
            # independent of the legacy shoulder-point emulation.
            legacy_x_scale = [
                max(abs(D_bounds.initial), 1.0),
                80.0 * max(abs(MF_bounds.initial), 1.0),
            ]
        diffusion, mobile_fraction = _fit_curve(
            weighted_model,
            fit_time,
            weighted_target,
            [D_bounds.initial, MF_bounds.initial],
            [D_bounds.lower, MF_bounds.lower],
            [D_bounds.upper, MF_bounds.upper],
            legacy_matlab,
            legacy_x_scale,
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


def _fit_profile_model(
    profile_model,
    radius: np.ndarray,
    post_bleach_profile: np.ndarray,
    config: DiffusionFitConfig,
) -> dict:
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    k_bounds = config.k
    r_bounds = config.r_effective

    if _is_adjustable(k_bounds) and _is_adjustable(r_bounds):
        bleach_depth, r_effective = _fit_curve(
            profile_model,
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
            lambda r, k: profile_model(r, k, r_effective),
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
            lambda r, re: profile_model(r, bleach_depth, re),
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

    profile_prediction = profile_model(radius, bleach_depth, r_effective)
    return {
        "k": float(bleach_depth),
        "r_effective": float(r_effective),
        "post_bleach_fit": profile_prediction,
        "post_bleach_residuals": post_bleach_profile - profile_prediction,
    }


def _fit_initial_conditions(
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
    initial_conditions = _fit_profile_model(_pbp_model, radius, post_bleach_profile, config)
    bleach_depth = initial_conditions["k"]
    r_effective = initial_conditions["r_effective"]

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

    profile_prediction = _pbp_model(radius, bleach_depth, r_effective)
    return {
        "k": float(bleach_depth),
        "r_effective": float(r_effective),
        "post_bleach_fit": profile_prediction,
        "post_bleach_residuals": post_bleach_profile - profile_prediction,
    }


def _fit_simplified_initial_conditions(
    radius: np.ndarray,
    post_bleach_profile: np.ndarray,
    config: DiffusionFitConfig,
) -> dict:
    return _fit_profile_model(_simplified_pbp_model, radius, post_bleach_profile, config)


def _concatenate_profile_segments(
    datasets: List[FRAPDataset],
    profile_slice: slice,
) -> Tuple[np.ndarray, np.ndarray]:
    radius = np.concatenate([dataset.radius[profile_slice] for dataset in datasets])
    post_bleach_profile = np.concatenate(
        [dataset.post_bleach_profile[profile_slice] for dataset in datasets]
    )
    return radius, post_bleach_profile


def _fit_global_initial_conditions(
    datasets: List[FRAPDataset],
    frap_slice: slice,
    profile_slice: slice,
    nominal_radius: float,
    config: DiffusionFitConfig,
) -> dict:
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    k_bounds = config.k
    r_bounds = config.r_effective
    radius, post_bleach_profile = _concatenate_profile_segments(datasets, profile_slice)
    initial_conditions = _fit_profile_model(_pbp_model, radius, post_bleach_profile, config)
    bleach_depth = initial_conditions["k"]
    r_effective = initial_conditions["r_effective"]

    if _is_adjustable(k_bounds):
        frap_initial = np.asarray(
            [dataset.corrected_frap[frap_slice][0] for dataset in datasets],
            dtype=float,
        )
        bleach_depth = _fit_curve(
            lambda t, k: kang_frap(t, r_effective, nominal_radius, 1.0, k),
            np.zeros_like(frap_initial),
            frap_initial,
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

    profile_prediction = _pbp_model(radius, bleach_depth, r_effective)
    return {
        "k": float(bleach_depth),
        "r_effective": float(r_effective),
        "post_bleach_fit": profile_prediction,
        "post_bleach_residuals": post_bleach_profile - profile_prediction,
    }


def _fit_global_residual(
    residual,
    initial: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    legacy_matlab: bool,
    legacy_x_scale: np.ndarray | None = None,
) -> np.ndarray:
    if legacy_matlab:
        return legacy_matlab_trf(
            residual,
            initial,
            lower,
            upper,
            legacy_x_scale if legacy_x_scale is not None else np.maximum(np.abs(initial), 1.0),
        )

    result = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=max(200, 100 * initial.size),
    )
    return result.x


def _fit_global_diffusion_components(
    datasets: List[FRAPDataset],
    frap_slice: slice,
    profile_slice: slice,
    nominal_radius: float,
    config: DiffusionFitConfig,
) -> dict:
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    global_fit = _fit_global_initial_conditions(
        datasets,
        frap_slice,
        profile_slice,
        nominal_radius,
        config,
    )
    bleach_depth = global_fit["k"]
    r_effective = global_fit["r_effective"]
    pieces = []

    for dataset in datasets:
        fit_time = dataset.time[frap_slice]
        corrected_frap = dataset.corrected_frap[frap_slice]
        radius = dataset.radius[profile_slice]
        post_bleach_profile = dataset.post_bleach_profile[profile_slice]
        profile_prediction = _pbp_model(radius, bleach_depth, r_effective)
        dataset.k = bleach_depth
        dataset.r_effective = r_effective
        dataset.post_bleach_fit = profile_prediction
        dataset.post_bleach_residuals = post_bleach_profile - profile_prediction
        dataset.fit_time = fit_time

        weights = corrected_frap.sum()
        pieces.append(
            {
                "dataset": dataset,
                "time": fit_time,
                "frap": corrected_frap,
                "initial": corrected_frap[0],
                "weights": weights,
                "target": corrected_frap / (fit_time + weights),
            }
        )

    D_bounds = config.diffusion_coefficient
    MF_bounds = config.mobile_fraction

    def residual(params: np.ndarray) -> np.ndarray:
        diffusion, mobile_fraction = params
        residuals = []
        for piece in pieces:
            numerator = kang_frap(
                piece["time"],
                r_effective,
                nominal_radius,
                diffusion,
                bleach_depth,
            ) * mobile_fraction
            numerator += (1.0 - mobile_fraction) * piece["initial"]
            prediction = numerator / (piece["time"] + piece["weights"])
            residuals.append(prediction - piece["target"])
        return np.concatenate(residuals)

    if _is_adjustable(D_bounds) and _is_adjustable(MF_bounds):
        initial = np.asarray([D_bounds.initial, MF_bounds.initial], dtype=float)
        lower = np.asarray([D_bounds.lower, MF_bounds.lower], dtype=float)
        upper = np.asarray([D_bounds.upper, MF_bounds.upper], dtype=float)
        diffusion, mobile_fraction = _fit_global_residual(
            residual,
            initial,
            lower,
            upper,
            legacy_matlab,
            np.asarray([max(abs(D_bounds.initial), 1.0), 80.0 * max(abs(MF_bounds.initial), 1.0)]),
        )
    elif _is_adjustable(D_bounds):
        mobile_fraction = MF_bounds.initial

        def diffusion_residual(params: np.ndarray) -> np.ndarray:
            return residual(np.asarray([params[0], mobile_fraction], dtype=float))

        diffusion = _fit_global_residual(
            diffusion_residual,
            np.asarray([D_bounds.initial], dtype=float),
            np.asarray([D_bounds.lower], dtype=float),
            np.asarray([D_bounds.upper], dtype=float),
            legacy_matlab,
            np.asarray([max(abs(D_bounds.initial), 1.0)]),
        )[0]
    elif _is_adjustable(MF_bounds):
        diffusion = D_bounds.initial

        def mf_residual(params: np.ndarray) -> np.ndarray:
            return residual(np.asarray([diffusion, params[0]], dtype=float))

        mobile_fraction = _fit_global_residual(
            mf_residual,
            np.asarray([MF_bounds.initial], dtype=float),
            np.asarray([MF_bounds.lower], dtype=float),
            np.asarray([MF_bounds.upper], dtype=float),
            legacy_matlab,
            np.asarray([80.0 * max(abs(MF_bounds.initial), 1.0)]),
        )[0]
    else:
        diffusion = D_bounds.initial
        mobile_fraction = MF_bounds.initial

    residual_vectors = []
    for piece in pieces:
        dataset = piece["dataset"]
        numerator = kang_frap(
            piece["time"],
            r_effective,
            nominal_radius,
            diffusion,
            bleach_depth,
        ) * mobile_fraction
        numerator += (1.0 - mobile_fraction) * piece["initial"]
        weighted_prediction = numerator / (piece["time"] + piece["weights"])
        residual_vector = piece["target"] - weighted_prediction
        dataset.diffusion_coefficient = float(diffusion)
        dataset.mobile_fraction = float(mobile_fraction)
        dataset.frap_fit = weighted_prediction * (piece["time"] + piece["weights"])
        dataset.frap_residuals = residual_vector
        dataset.sum_squared_residuals = float(np.dot(residual_vector, residual_vector))
        residual_vectors.append(residual_vector)

    return {
        "k": float(bleach_depth),
        "r_effective": float(r_effective),
        "diffusion_coefficient": float(diffusion),
        "mobile_fraction": float(mobile_fraction),
        "sum_squared_residuals": float(sum(np.dot(res, res) for res in residual_vectors)),
    }


def _fit_simplified_kang_components(
    fit_time: np.ndarray,
    corrected_frap: np.ndarray,
    radius: np.ndarray,
    post_bleach_profile: np.ndarray,
    nominal_radius: float,
    config: DiffusionFitConfig,
) -> dict:
    initial_conditions = _fit_simplified_initial_conditions(
        radius,
        post_bleach_profile,
        config,
    )
    bleach_depth = initial_conditions["k"]
    r_effective = initial_conditions["r_effective"]

    steady_state = float(corrected_frap[-1])
    half_time = estimate_recovery_half_time(fit_time, corrected_frap, steady_state)
    D_bounds = config.diffusion_coefficient
    if _is_adjustable(D_bounds):
        diffusion = simplified_kang_diffusion_coefficient(nominal_radius, r_effective, half_time)
        diffusion = float(np.clip(diffusion, D_bounds.lower, D_bounds.upper))
    else:
        diffusion = D_bounds.initial

    MF_bounds = config.mobile_fraction
    if _is_adjustable(MF_bounds):
        denominator = 1.0 - corrected_frap[0]
        if np.isclose(denominator, 0.0):
            mobile_fraction = MF_bounds.initial
        else:
            mobile_fraction = (steady_state - corrected_frap[0]) / denominator
            if not np.isfinite(mobile_fraction):
                mobile_fraction = MF_bounds.initial
        mobile_fraction = float(np.clip(mobile_fraction, MF_bounds.lower, MF_bounds.upper))
    else:
        mobile_fraction = MF_bounds.initial

    evaluation = _evaluate_simplified_kang_components(
        fit_time,
        corrected_frap,
        radius,
        post_bleach_profile,
        nominal_radius,
        bleach_depth,
        r_effective,
        diffusion,
        mobile_fraction,
    )
    return {
        **initial_conditions,
        "half_time": float(half_time),
        "diffusion_coefficient": float(diffusion),
        "mobile_fraction": float(mobile_fraction),
        "sum_squared_residuals": evaluation["sum_squared_residuals"],
        "frap_fit": evaluation["frap_fit"],
        "frap_residuals": evaluation["frap_residuals"],
    }


def _fit_global_simplified_kang_components(
    datasets: List[FRAPDataset],
    averaged_time: np.ndarray,
    averaged_frap: np.ndarray,
    frap_slice: slice,
    profile_slice: slice,
    nominal_radius: float,
    config: DiffusionFitConfig,
) -> dict:
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    radius, post_bleach_profile = _concatenate_profile_segments(datasets, profile_slice)
    global_fit = _fit_simplified_initial_conditions(radius, post_bleach_profile, config)
    k = global_fit["k"]
    r_effective = global_fit["r_effective"]

    pieces = []
    for dataset in datasets:
        fit_time = dataset.time[frap_slice]
        corrected_frap = dataset.corrected_frap[frap_slice]
        weights = corrected_frap.sum()
        pieces.append(
            {
                "dataset": dataset,
                "time": fit_time,
                "frap": corrected_frap,
                "initial": corrected_frap[0],
                "weights": weights,
                "target": corrected_frap / (fit_time + weights),
            }
        )

    D_bounds = config.diffusion_coefficient
    MF_bounds = config.mobile_fraction

    try:
        half_time_initial = estimate_recovery_half_time(
            averaged_time,
            averaged_frap,
            float(averaged_frap[-1]),
        )
        diffusion_initial = simplified_kang_diffusion_coefficient(
            nominal_radius,
            r_effective,
            half_time_initial,
        )
    except ValueError:
        diffusion_initial = D_bounds.initial
    if not np.isfinite(diffusion_initial) or diffusion_initial <= D_bounds.lower:
        diffusion_initial = D_bounds.initial
    if np.isfinite(D_bounds.upper):
        diffusion_initial = min(diffusion_initial, D_bounds.upper)
    diffusion_initial = max(diffusion_initial, max(D_bounds.lower, np.finfo(float).eps))

    denominator = 1.0 - averaged_frap[0]
    if np.isclose(denominator, 0.0):
        mobile_fraction_initial = MF_bounds.initial
    else:
        mobile_fraction_initial = (float(averaged_frap[-1]) - averaged_frap[0]) / denominator
    if not np.isfinite(mobile_fraction_initial):
        mobile_fraction_initial = MF_bounds.initial
    mobile_fraction_initial = float(
        np.clip(mobile_fraction_initial, MF_bounds.lower, MF_bounds.upper)
    )

    def residual(params: np.ndarray) -> np.ndarray:
        diffusion, mobile_fraction = params
        residuals = []
        for piece in pieces:
            prediction = simplified_kang_recovery(
                piece["time"],
                r_effective,
                nominal_radius,
                diffusion,
                mobile_fraction,
                piece["initial"],
                bleach_depth=k,
            )
            weighted_prediction = prediction / (piece["time"] + piece["weights"])
            residuals.append(weighted_prediction - piece["target"])
        return np.concatenate(residuals)

    if _is_adjustable(D_bounds) and _is_adjustable(MF_bounds):
        initial = np.asarray([diffusion_initial, mobile_fraction_initial], dtype=float)
        lower = np.asarray([max(D_bounds.lower, np.finfo(float).eps), MF_bounds.lower], dtype=float)
        upper = np.asarray([D_bounds.upper, MF_bounds.upper], dtype=float)
        diffusion, mobile_fraction = _fit_global_residual(
            residual,
            initial,
            lower,
            upper,
            legacy_matlab,
            np.asarray(
                [
                    max(abs(diffusion_initial), 1.0),
                    80.0 * max(abs(MF_bounds.initial), 1.0),
                ]
            ),
        )
    elif _is_adjustable(D_bounds):
        mobile_fraction = MF_bounds.initial

        def diffusion_residual(params: np.ndarray) -> np.ndarray:
            return residual(np.asarray([params[0], mobile_fraction], dtype=float))

        diffusion = _fit_global_residual(
            diffusion_residual,
            np.asarray([diffusion_initial], dtype=float),
            np.asarray([max(D_bounds.lower, np.finfo(float).eps)], dtype=float),
            np.asarray([D_bounds.upper], dtype=float),
            legacy_matlab,
            np.asarray([max(abs(diffusion_initial), 1.0)]),
        )[0]
    elif _is_adjustable(MF_bounds):
        diffusion = D_bounds.initial

        def mf_residual(params: np.ndarray) -> np.ndarray:
            return residual(np.asarray([diffusion, params[0]], dtype=float))

        mobile_fraction = _fit_global_residual(
            mf_residual,
            np.asarray([mobile_fraction_initial], dtype=float),
            np.asarray([MF_bounds.lower], dtype=float),
            np.asarray([MF_bounds.upper], dtype=float),
            legacy_matlab,
            np.asarray([80.0 * max(abs(MF_bounds.initial), 1.0)]),
        )[0]
    else:
        diffusion = D_bounds.initial
        mobile_fraction = MF_bounds.initial

    half_time = (
        simplified_kang_diffusion_coefficient(nominal_radius, r_effective, 1.0) / diffusion
        if diffusion > 0.0
        else np.inf
    )

    residual_vectors = []
    for piece in pieces:
        dataset = piece["dataset"]
        fit_time = piece["time"]
        corrected_frap = piece["frap"]
        dataset_radius = dataset.radius[profile_slice]
        dataset_profile = dataset.post_bleach_profile[profile_slice]
        evaluation = _evaluate_simplified_kang_components(
            fit_time,
            corrected_frap,
            dataset_radius,
            dataset_profile,
            nominal_radius,
            k,
            r_effective,
            diffusion,
            mobile_fraction,
        )
        dataset.k = k
        dataset.r_effective = r_effective
        dataset.half_time = half_time
        dataset.diffusion_coefficient = diffusion
        dataset.mobile_fraction = mobile_fraction
        dataset.post_bleach_fit = evaluation["post_bleach_fit"]
        dataset.post_bleach_residuals = evaluation["post_bleach_residuals"]
        dataset.fit_time = fit_time
        dataset.frap_fit = evaluation["frap_fit"]
        dataset.frap_residuals = evaluation["frap_residuals"]
        dataset.sum_squared_residuals = evaluation["sum_squared_residuals"]
        residual_vectors.append(evaluation["frap_residuals"])

    return {
        **global_fit,
        "half_time": float(half_time),
        "diffusion_coefficient": float(diffusion),
        "mobile_fraction": float(mobile_fraction),
        "sum_squared_residuals": float(sum(np.dot(res, res) for res in residual_vectors)),
    }


def _evaluate_simplified_kang_components(
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
    prediction = simplified_kang_recovery(
        fit_time,
        r_effective,
        nominal_radius,
        diffusion,
        mobile_fraction,
        frap_initial,
        bleach_depth=k,
    )
    weighted_prediction = prediction / (fit_time + weights)
    profile_prediction = _simplified_pbp_model(radius, k, r_effective)
    return {
        "sum_squared_residuals": float(np.sum((weighted_target - weighted_prediction) ** 2)),
        "post_bleach_fit": profile_prediction,
        "post_bleach_residuals": post_bleach_profile - profile_prediction,
        "frap_fit": prediction,
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
    fit_mode = _resolve_fit_mode(config)

    if fit_mode == "individual":
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
            dataset.half_time = None
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
        half_time = None
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
    elif fit_mode == "average_curve":
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
        half_time = None
        diffusion = fit["diffusion_coefficient"]
        mobile_fraction = fit["mobile_fraction"]
        sum_squared_residuals = fit["sum_squared_residuals"]
        profile_prediction = fit["post_bleach_fit"]
        profile_residuals = fit["post_bleach_residuals"]
        frap_fit = fit["frap_fit"]
        frap_residuals = fit["frap_residuals"]
    elif fit_mode == "simplified_kang":
        for dataset in datasets:
            fit = _fit_simplified_kang_components(
                dataset.time[frap_slice],
                dataset.corrected_frap[frap_slice],
                dataset.radius[profile_slice],
                dataset.post_bleach_profile[profile_slice],
                nominal_radius,
                config,
            )
            dataset.k = fit["k"]
            dataset.r_effective = fit["r_effective"]
            dataset.half_time = fit["half_time"]
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
        half_time = float(np.nanmean([dataset.half_time for dataset in datasets]))
        diffusion = float(np.nanmean([dataset.diffusion_coefficient for dataset in datasets]))
        mobile_fraction = float(np.nanmean([dataset.mobile_fraction for dataset in datasets]))
        aggregate = _evaluate_simplified_kang_components(
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
    elif fit_mode == "simplified_kang_global":
        fit = _fit_global_simplified_kang_components(
            datasets,
            ref_time,
            averaged_frap,
            frap_slice,
            profile_slice,
            nominal_radius,
            config,
        )
        k = fit["k"]
        r_effective = fit["r_effective"]
        half_time = fit["half_time"]
        diffusion = fit["diffusion_coefficient"]
        mobile_fraction = fit["mobile_fraction"]
        aggregate = _evaluate_simplified_kang_components(
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
        sum_squared_residuals = fit["sum_squared_residuals"]
        profile_prediction = aggregate["post_bleach_fit"]
        profile_residuals = aggregate["post_bleach_residuals"]
        frap_fit = aggregate["frap_fit"]
        frap_residuals = aggregate["frap_residuals"]
    elif fit_mode == "global":
        fit = _fit_global_diffusion_components(
            datasets,
            frap_slice,
            profile_slice,
            nominal_radius,
            config,
        )
        for dataset in datasets:
            dataset.half_time = None
        k = fit["k"]
        r_effective = fit["r_effective"]
        half_time = None
        diffusion = fit["diffusion_coefficient"]
        mobile_fraction = fit["mobile_fraction"]
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
        sum_squared_residuals = fit["sum_squared_residuals"]
        profile_prediction = aggregate["post_bleach_fit"]
        profile_residuals = aggregate["post_bleach_residuals"]
        frap_fit = aggregate["frap_fit"]
        frap_residuals = aggregate["frap_residuals"]
    else:
        raise ValueError(f"Unsupported diffusion fit mode: {fit_mode!r}")

    valid_cmf = [d.corrected_mobile_fraction for d in datasets if d.corrected_mobile_fraction is not None]
    corrected_mobile_fraction = float(np.nanmean(valid_cmf)) if valid_cmf else None

    return DiffusionFitResult(
        datasets=datasets,
        k=k,
        r_effective=r_effective,
        half_time=half_time,
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
        fit_mode="global",
    )
