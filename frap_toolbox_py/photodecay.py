from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np
from scipy.optimize import curve_fit

from .optimization import legacy_matlab_trf
from .types import FRAPDataset


def estimate_decay_rate(
    datasets: Iterable[FRAPDataset],
    post_bleach_frame: int,
    decay_fit_range: Tuple[int, int],
    initial_rate: float,
    lower_bound: float,
    upper_bound: float,
    legacy_matlab: bool = False,
) -> float:
    """Estimate photodecay using an exponential fit to the averaged curve."""

    aligned_fraps = []
    reference_time = None
    for dataset in datasets:
        time_offset = dataset.time - dataset.time[0]
        if reference_time is None:
            reference_time = time_offset
            aligned_fraps.append(dataset.norm_frap)
            continue

        if legacy_matlab:
            is_close = np.abs(reference_time - time_offset) <= (
                0.1 * np.maximum(np.abs(reference_time), np.abs(time_offset))
                + np.finfo(float).eps
            )
            if not np.all(is_close):
                raise ValueError("Time vectors differ beyond MATLAB legacy tolerance.")
            aligned_fraps.append(dataset.norm_frap)
            continue

        current_time = np.maximum.accumulate(time_offset.astype(float, copy=True))
        diffs = np.diff(current_time)
        if np.any(diffs <= 0):
            eps = np.finfo(current_time.dtype).eps
            for idx in range(1, current_time.size):
                if current_time[idx] <= current_time[idx - 1]:
                    current_time[idx] = current_time[idx - 1] + eps
        if not np.allclose(current_time, reference_time, rtol=1e-3, atol=1e-6):
            interpolated = np.interp(reference_time, current_time, dataset.norm_frap)
            aligned_fraps.append(interpolated)
        else:
            aligned_fraps.append(dataset.norm_frap)

    if reference_time is None:
        raise ValueError("No datasets provided for photodecay estimation.")

    mean_time = reference_time
    mean_frap = np.mean(aligned_fraps, axis=0)

    start, end = decay_fit_range
    fit_t = mean_time[start:end]
    fit_f = mean_frap[start:end]

    def model(x, amplitude, rate):
        return amplitude * np.exp(-rate * x)

    if legacy_matlab:
        popt = legacy_matlab_trf(
            lambda params: model(fit_t, params[0], params[1]) - fit_f,
            initial=np.asarray([0.9, initial_rate], dtype=float),
            lower=np.asarray([0.0, lower_bound], dtype=float),
            upper=np.asarray([2.0, upper_bound], dtype=float),
            x_scale=np.ones(2),
        )
    else:
        popt, _ = curve_fit(
            model,
            fit_t,
            fit_f,
            p0=[0.9, initial_rate],
            bounds=([0.0, lower_bound], [2.0, upper_bound]),
        )

    decay_rate = popt[1]
    return decay_rate
