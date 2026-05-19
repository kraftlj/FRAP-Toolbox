from __future__ import annotations

import warnings
from typing import List, Sequence, Tuple

import numpy as np
from scipy.optimize import OptimizeWarning, least_squares

from ..optimization import legacy_matlab_trf
from ..photodecay import estimate_decay_rate
from ..types import (
    BasicInputs,
    FitBounds,
    FRAPDataset,
    ReactionFitConfig,
    ReactionFitResult,
)


def _is_adjustable(bounds: FitBounds) -> bool:
    return bounds.mode.lower() != "fixed"


def _slice_range(range_tuple: Tuple[int, int]) -> slice:
    start, end = range_tuple
    return slice(start, end)


def reaction1_curve(t: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Evaluate the FRAP-Toolbox one-component reaction model."""

    return a - b * np.exp(-c * np.asarray(t, dtype=float))


def reaction2_curve(
    t: np.ndarray,
    a: float,
    b: float,
    c: float,
    d: float,
    f: float,
) -> np.ndarray:
    """Evaluate the FRAP-Toolbox two-component reaction model."""

    t = np.asarray(t, dtype=float)
    return a - b * np.exp(-c * t) - d * np.exp(-f * t)


def _parameter_bounds(config: ReactionFitConfig, model_order: int) -> tuple[list[str], list[FitBounds]]:
    if model_order == 1:
        return ["a", "b", "c"], [config.a, config.b, config.c]
    if model_order == 2:
        if config.d is None or config.f is None:
            raise ValueError("Reaction 2 fitting requires d and f parameter bounds.")
        return ["a", "b", "c", "d", "f"], [config.a, config.b, config.c, config.d, config.f]
    raise ValueError("Reaction model order must be 1 or 2.")


def _reaction_curve(model_order: int, t: np.ndarray, params: Sequence[float]) -> np.ndarray:
    if model_order == 1:
        return reaction1_curve(t, params[0], params[1], params[2])
    return reaction2_curve(t, params[0], params[1], params[2], params[3], params[4])


def _evaluate_reaction_components(
    fit_time: np.ndarray,
    corrected_frap: np.ndarray,
    model_order: int,
    params: Sequence[float],
) -> dict:
    weights = fit_time + corrected_frap.sum()
    weighted_target = corrected_frap / weights
    prediction = _reaction_curve(model_order, fit_time, params)
    weighted_prediction = prediction / weights
    residuals = weighted_target - weighted_prediction
    return {
        "sum_squared_residuals": float(np.dot(residuals, residuals)),
        "frap_fit": prediction,
        "frap_residuals": residuals,
    }


def _fit_reaction_components(
    fit_time: np.ndarray,
    corrected_frap: np.ndarray,
    model_order: int,
    config: ReactionFitConfig,
) -> dict:
    names, bounds = _parameter_bounds(config, model_order)
    initial = np.asarray([bound.initial for bound in bounds], dtype=float)
    lower = np.asarray([bound.lower for bound in bounds], dtype=float)
    upper = np.asarray([bound.upper for bound in bounds], dtype=float)
    fixed = np.asarray([not _is_adjustable(bound) for bound in bounds], dtype=bool)
    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    weights = fit_time + corrected_frap.sum()

    def normalized_params(raw_params: np.ndarray) -> np.ndarray:
        params = np.asarray(raw_params, dtype=float).copy()
        params[fixed] = initial[fixed]
        return params

    def residual(raw_params: np.ndarray) -> np.ndarray:
        params = normalized_params(raw_params)
        prediction = _reaction_curve(model_order, fit_time, params)
        unweighted = prediction - corrected_frap
        if model_order == 1:
            return unweighted / weights
        return unweighted

    if legacy_matlab:
        fitted = legacy_matlab_trf(
            residual,
            initial,
            lower,
            upper,
            np.ones_like(initial),
        )
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=OptimizeWarning)
            result = least_squares(
                residual,
                x0=initial,
                bounds=(lower, upper),
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
                max_nfev=max(200, 100 * initial.size),
            )
        fitted = result.x

    params = normalized_params(fitted)
    evaluation = _evaluate_reaction_components(
        fit_time,
        corrected_frap,
        model_order,
        params,
    )
    return {
        "parameters": {name: float(value) for name, value in zip(names, params)},
        **evaluation,
    }


def _prepare_average_segments(
    datasets: List[FRAPDataset],
    frap_slice: slice,
) -> Tuple[np.ndarray, np.ndarray]:
    ref_length = len(datasets[0].corrected_frap[frap_slice])
    ref_time = datasets[0].time[frap_slice]

    frap_matrix = []
    for dataset in datasets:
        segment = dataset.corrected_frap[frap_slice]
        if len(segment) != ref_length:
            raise ValueError("FRAP segments must have consistent length for averaging.")
        frap_matrix.append(segment)

    return ref_time, np.asarray(frap_matrix).mean(axis=0)


def fit_reaction_model(
    datasets: List[FRAPDataset],
    inputs: BasicInputs,
    config: ReactionFitConfig,
    model_order: int,
) -> ReactionFitResult:
    """Fit Reaction 1 or Reaction 2 FRAP models using MATLAB-compatible rules."""

    if not datasets:
        raise ValueError("At least one dataset is required for fitting.")
    if model_order not in {1, 2}:
        raise ValueError("Reaction model order must be 1 or 2.")

    legacy_matlab = config.optimizer_mode.lower() == "legacy_matlab"
    if inputs.normalize_by_cell:
        decay_rate = 0.0
    elif _is_adjustable(config.decay_rate):
        decay_rate = estimate_decay_rate(
            datasets,
            inputs.post_bleach_frame,
            config.decay_fit_range,
            config.decay_rate.initial,
            config.decay_rate.lower,
            config.decay_rate.upper,
            legacy_matlab=legacy_matlab,
        )
    else:
        decay_rate = config.decay_rate.initial

    for dataset in datasets:
        original_time = dataset.time.astype(float, copy=True)
        time_from_first = original_time - original_time[0]
        time_from_post = original_time - original_time[inputs.post_bleach_frame]
        if inputs.normalize_by_cell:
            dataset.corrected_frap = dataset.norm_frap.copy()
        else:
            dataset.corrected_frap = dataset.norm_frap / np.exp(-decay_rate * time_from_first)
        dataset.time = time_from_post

    frap_slice = _slice_range(config.frap_range)
    averaged_time, averaged_frap = _prepare_average_segments(datasets, frap_slice)

    if config.fit_averaged_data:
        fit = _fit_reaction_components(averaged_time, averaged_frap, model_order, config)
        parameters = fit["parameters"]
        sum_squared_residuals = fit["sum_squared_residuals"]
        averaged_fit = fit["frap_fit"]
        averaged_residuals = fit["frap_residuals"]
    else:
        parameter_rows = []
        for dataset in datasets:
            fit_time = dataset.time[frap_slice]
            corrected_frap = dataset.corrected_frap[frap_slice]
            fit = _fit_reaction_components(fit_time, corrected_frap, model_order, config)
            dataset.reaction_parameters = fit["parameters"]
            dataset.fit_time = fit_time
            dataset.frap_fit = fit["frap_fit"]
            dataset.frap_residuals = fit["frap_residuals"]
            dataset.sum_squared_residuals = fit["sum_squared_residuals"]
            parameter_rows.append(fit["parameters"])

        names, _ = _parameter_bounds(config, model_order)
        parameters = {
            name: float(np.mean([row[name] for row in parameter_rows]))
            for name in names
        }
        aggregate = _evaluate_reaction_components(
            averaged_time,
            averaged_frap,
            model_order,
            [parameters[name] for name in names],
        )
        sum_squared_residuals = aggregate["sum_squared_residuals"]
        averaged_fit = aggregate["frap_fit"]
        averaged_residuals = aggregate["frap_residuals"]

    return ReactionFitResult(
        model_order=model_order,
        datasets=datasets,
        parameters=parameters,
        decay_rate=float(decay_rate),
        sum_squared_residuals=sum_squared_residuals,
        averaged_time=averaged_time,
        averaged_frap=averaged_frap,
        averaged_frap_fit=averaged_fit,
        averaged_frap_residuals=averaged_residuals,
    )


def default_reaction1_config(
    frap_length: int,
    post_bleach_frame: int,
) -> ReactionFitConfig:
    decay_start = max(frap_length - 51, post_bleach_frame)
    return ReactionFitConfig(
        a=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        b=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        c=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        decay_rate=FitBounds(0.001, 0.0, np.inf, "Adjustable"),
        frap_range=(post_bleach_frame, frap_length),
        decay_fit_range=(decay_start, frap_length),
        fit_averaged_data=False,
    )


def default_reaction2_config(
    frap_length: int,
    post_bleach_frame: int,
) -> ReactionFitConfig:
    decay_start = max(frap_length - 51, post_bleach_frame)
    return ReactionFitConfig(
        a=FitBounds(1.0, 0.0, np.inf, "Adjustable"),
        b=FitBounds(0.5, 0.0, np.inf, "Adjustable"),
        c=FitBounds(0.05, 0.0, np.inf, "Adjustable"),
        d=FitBounds(0.5, 0.0, np.inf, "Adjustable"),
        f=FitBounds(0.0005, 0.0, np.inf, "Adjustable"),
        decay_rate=FitBounds(0.001, 0.0, np.inf, "Adjustable"),
        frap_range=(post_bleach_frame, frap_length),
        decay_fit_range=(decay_start, frap_length),
        fit_averaged_data=True,
    )
