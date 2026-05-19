from __future__ import annotations

import numpy as np
from numpy.linalg import norm
from scipy.linalg import svd
from scipy.optimize._lsq.common import CL_scaling_vector, make_strictly_feasible
from scipy.optimize._lsq.common import solve_lsq_trust_region, update_tr_radius
from scipy.optimize._lsq.trf import select_step
from scipy.optimize._numdiff import approx_derivative


def legacy_matlab_trf(
    residual,
    initial: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    x_scale: np.ndarray,
    ftol: float = 1e-6,
    xtol: float = 1e-6,
    gtol: float = 1e-6,
    max_nfev: int | None = None,
) -> np.ndarray:
    """Emulate the MATLAB-era trust-region path for archived FRAP fits.

    MATLAB's historical `lsqcurvefit` default was a Coleman-Li
    trust-region-reflective method with absolute TolFun/TolX-style stopping.
    SciPy's public `least_squares` uses a related method, but its small dense
    path and relative function-tolerance test continue past several of the
    shoulder points reported in the FRAP-Toolbox user guide.
    """

    x = make_strictly_feasible(initial.astype(float, copy=True), lower, upper)
    x_scale = np.asarray(x_scale, dtype=float)
    scale_inv = 1.0 / x_scale
    if max_nfev is None:
        max_nfev = max(200, 100 * x.size)

    f = residual(x)
    nfev = 1
    jacobian = approx_derivative(residual, x, method="2-point", bounds=(lower, upper))
    cost = 0.5 * float(np.dot(f, f))
    gradient = jacobian.T.dot(f)

    v, dv = CL_scaling_vector(x, gradient, lower, upper)
    v[dv != 0] *= scale_inv[dv != 0]
    trust_radius = norm(x * scale_inv / np.sqrt(v))
    if trust_radius == 0:
        trust_radius = 1.0
    alpha = 0.0

    while nfev < max_nfev:
        v, dv = CL_scaling_vector(x, gradient, lower, upper)
        optimality = norm(gradient * v, ord=np.inf)
        if optimality < gtol:
            break

        scaled_v = v.copy()
        scaled_v[dv != 0] *= scale_inv[dv != 0]
        d = np.sqrt(scaled_v) * x_scale
        diag_h = np.maximum(gradient * dv * x_scale, 0.0)
        g_h = d * gradient

        m, n = jacobian.shape
        f_augmented = np.zeros(m + n)
        f_augmented[:m] = f
        jacobian_augmented = np.empty((m + n, n))
        jacobian_augmented[:m] = jacobian * d
        jacobian_h = jacobian_augmented[:m]
        jacobian_augmented[m:] = np.diag(np.sqrt(diag_h))

        u, s, vt = svd(jacobian_augmented, full_matrices=False)
        uf = u.T.dot(f_augmented)
        v_svd = vt.T
        theta = max(0.995, 1.0 - optimality)

        actual_reduction = -1.0
        step = np.zeros_like(x)
        step_h = np.zeros_like(x)
        ratio = 0.0
        new_trust_radius = trust_radius
        while actual_reduction <= 0 and nfev < max_nfev:
            p_h, alpha, _ = solve_lsq_trust_region(
                n,
                m,
                uf,
                s,
                v_svd,
                trust_radius,
                initial_alpha=alpha,
            )
            p = d * p_h
            step, step_h, predicted_reduction = select_step(
                x,
                jacobian_h,
                diag_h,
                g_h,
                p,
                p_h,
                d,
                trust_radius,
                lower,
                upper,
                theta,
            )

            x_new = make_strictly_feasible(x + step, lower, upper, rstep=0)
            f_new = residual(x_new)
            nfev += 1
            step_h_norm = norm(step_h)
            if not np.all(np.isfinite(f_new)):
                trust_radius = 0.25 * step_h_norm
                continue

            cost_new = 0.5 * float(np.dot(f_new, f_new))
            actual_reduction = cost - cost_new
            new_trust_radius, ratio = update_tr_radius(
                trust_radius,
                actual_reduction,
                predicted_reduction,
                step_h_norm,
                step_h_norm > 0.95 * trust_radius,
            )
            if actual_reduction <= 0:
                alpha *= trust_radius / new_trust_radius
                trust_radius = new_trust_radius

        if actual_reduction <= 0:
            break

        previous_x = x
        x = x_new
        f = f_new
        cost = cost_new
        jacobian = approx_derivative(residual, x, method="2-point", bounds=(lower, upper))
        gradient = jacobian.T.dot(f)
        alpha *= trust_radius / new_trust_radius
        trust_radius = new_trust_radius

        residual_sum_square_reduction = 2.0 * actual_reduction
        if residual_sum_square_reduction < ftol and ratio > 0.75:
            break
        if norm(step) < xtol * (xtol + norm(previous_x)):
            break

    return x
