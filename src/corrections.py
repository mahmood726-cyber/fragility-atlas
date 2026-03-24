"""Publication bias corrections: trim-and-fill (Duval-Tweedie L0) and PET-PEESE."""

import math
import numpy as np
from src.estimators import meta_analysis, MetaResult


def trim_and_fill(yi: np.ndarray, sei: np.ndarray, estimator: str = 'DL',
                  ci_method: str = 'Wald', conf_level: float = 0.95,
                  max_iter: int = 20) -> MetaResult:
    """Duval-Tweedie trim-and-fill (L0 estimator, right side).

    Estimates k0 missing studies, imputes them, and re-runs the meta-analysis.
    """
    k = len(yi)
    if k < 3:
        return meta_analysis(yi, sei, estimator, ci_method, conf_level)

    # Iterative: estimate center, count asymmetry, impute, repeat
    yi_current = yi.copy()
    sei_current = sei.copy()

    for _ in range(max_iter):
        # Estimate center using current data
        result = meta_analysis(yi_current, sei_current, estimator, ci_method, conf_level)
        theta0 = result.theta

        # Rank studies by distance from center
        di = yi[:k] - theta0  # use original studies only
        ranks = np.argsort(np.abs(di))
        signed_ranks = np.zeros(k)
        for i, r in enumerate(ranks):
            signed_ranks[r] = (i + 1) * np.sign(di[r])

        # L0 estimator for k0
        T_plus = np.sum(signed_ranks[signed_ranks > 0])
        # k0 = max(0, round(T_plus - k*(k+1)/4) ... simplified L0)
        # Duval 2000: k0 = max(0, round((4*S+ - k*(k+1)) / (2*k + 1)))
        S_plus = T_plus
        k0_est = max(0, round((4 * S_plus - k * (k + 1)) / (2 * k + 1)))

        if k0_est == 0:
            break

        # Impute missing studies by reflecting the k0 most extreme studies
        order = np.argsort(di)
        # Most extreme on the positive side (assuming right-side asymmetry)
        extreme_idx = order[-k0_est:]

        yi_filled = np.concatenate([yi, 2 * theta0 - yi[extreme_idx]])
        sei_filled = np.concatenate([sei, sei[extreme_idx]])

        # Check convergence
        old_yi = yi_current
        yi_current = yi_filled
        sei_current = sei_filled

        if len(yi_current) == len(old_yi) and np.allclose(yi_current, old_yi, atol=1e-10):
            break

    # Final meta-analysis on filled data
    return meta_analysis(yi_current, sei_current, estimator, ci_method, conf_level)


def pet_peese(yi: np.ndarray, sei: np.ndarray, estimator: str = 'DL',
              ci_method: str = 'Wald', conf_level: float = 0.95) -> MetaResult:
    """Conditional PET-PEESE publication bias correction.

    PET: regress yi on sei (precision-effect test)
    PEESE: regress yi on sei² (precision-effect estimate with SE)
    Use PET if PET intercept p >= 0.05, else use PEESE.
    """
    k = len(yi)
    if k < 3:
        return meta_analysis(yi, sei, estimator, ci_method, conf_level)

    wi = 1.0 / (sei ** 2)

    # PET: WLS regression of yi on sei, weighted by 1/vi
    pet_result = _weighted_regression(yi, sei, wi, use_se_squared=False)
    pet_intercept, pet_se, pet_p = pet_result

    if pet_p >= 0.05:
        # Use PET intercept
        intercept, se_int, p_val = pet_result
    else:
        # Use PEESE intercept
        intercept, se_int, p_val = _weighted_regression(yi, sei, wi, use_se_squared=True)

    # Construct result (tau2 from base estimator for I²/Q reporting)
    base = meta_analysis(yi, sei, estimator, ci_method, conf_level)

    from src.utils import normal_quantile, normal_cdf
    alpha = 1.0 - conf_level
    z_crit = normal_quantile(1.0 - alpha / 2)
    ci_lo = intercept - z_crit * se_int
    ci_hi = intercept + z_crit * se_int

    return MetaResult(
        theta=float(intercept),
        se_theta=float(se_int),
        ci_lo=float(ci_lo),
        ci_hi=float(ci_hi),
        p_value=float(p_val),
        tau2=base.tau2,
        i2=base.i2,
        q_stat=base.q_stat,
        estimator=f'{estimator}+PETPEESE',
        ci_method=ci_method,
        k=k,
    )


def _weighted_regression(yi, sei, wi, use_se_squared=False):
    """WLS regression of yi on sei (or sei²), return (intercept, SE, p-value)."""
    from src.utils import normal_cdf

    x = sei ** 2 if use_se_squared else sei
    n = len(yi)

    # WLS: minimize sum(wi * (yi - a - b*x)^2)
    sw = np.sum(wi)
    sx = np.sum(wi * x)
    sy = np.sum(wi * yi)
    sxx = np.sum(wi * x ** 2)
    sxy = np.sum(wi * x * yi)

    denom = sw * sxx - sx ** 2
    if abs(denom) < 1e-300:
        # Degenerate: return simple weighted mean
        intercept = sy / sw
        se = 1.0 / math.sqrt(sw)
        z = intercept / se if se > 0 else 0
        p = 2.0 * (1.0 - normal_cdf(abs(z)))
        return intercept, se, p

    intercept = (sxx * sy - sx * sxy) / denom
    slope = (sw * sxy - sx * sy) / denom

    # Residual variance (sigma²)
    residuals = yi - intercept - slope * x
    sigma2 = np.sum(wi * residuals ** 2) / (n - 2) if n > 2 else 1.0

    # SE of intercept
    var_intercept = sigma2 * sxx / denom
    se = math.sqrt(max(0, var_intercept))

    z = intercept / se if se > 0 else 0
    p = 2.0 * (1.0 - normal_cdf(abs(z)))

    return intercept, se, p
