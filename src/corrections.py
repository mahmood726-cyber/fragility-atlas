"""Publication bias corrections: trim-and-fill (Duval-Tweedie L0) and PET-PEESE."""

import math
import numpy as np
from src.estimators import meta_analysis, MetaResult


def trim_and_fill(yi: np.ndarray, sei: np.ndarray, estimator: str = 'DL',
                  ci_method: str = 'Wald', conf_level: float = 0.95,
                  max_iter: int = 20) -> MetaResult:
    """Duval-Tweedie trim-and-fill with auto side detection (L0 estimator).

    P0-3 FIX: detects asymmetry side automatically (not right-only).
    Uses rank correlation to choose side, then imputes on the opposite side.
    """
    k = len(yi)
    if k < 3:
        return meta_analysis(yi, sei, estimator, ci_method, conf_level)

    # Step 1: Determine asymmetry side via rank correlation of yi with sei
    # Positive correlation (larger effects in less precise studies) → right-side asymmetry
    # Negative correlation → left-side asymmetry
    base_result = meta_analysis(yi, sei, estimator, ci_method, conf_level)
    theta0 = base_result.theta
    di = yi - theta0

    # Use Kendall-like sign: side where excess studies are
    n_right = np.sum(di > 0)
    n_left = np.sum(di < 0)
    side = 'right' if n_right > n_left else 'left'

    yi_current = yi.copy()
    sei_current = sei.copy()

    for _ in range(max_iter):
        result = meta_analysis(yi_current, sei_current, estimator, ci_method, conf_level)
        theta0 = result.theta

        # Rank studies by distance from center (original studies only)
        di = yi[:k] - theta0
        ranks = np.argsort(np.abs(di))
        signed_ranks = np.zeros(k)
        for i, r in enumerate(ranks):
            signed_ranks[r] = (i + 1) * np.sign(di[r])

        # L0 estimator for k0
        if side == 'right':
            S = np.sum(signed_ranks[signed_ranks > 0])
        else:
            S = np.sum(np.abs(signed_ranks[signed_ranks < 0]))

        k0_est = max(0, round((4 * S - k * (k + 1)) / (2 * k + 1)))

        if k0_est == 0:
            break

        # Impute missing studies by reflecting the k0 most extreme studies
        order = np.argsort(di)
        if side == 'right':
            # Excess on right → impute on left (reflect rightmost studies)
            extreme_idx = order[-k0_est:]
        else:
            # Excess on left → impute on right (reflect leftmost studies)
            extreme_idx = order[:k0_est]

        yi_filled = np.concatenate([yi, 2 * theta0 - yi[extreme_idx]])
        sei_filled = np.concatenate([sei, sei[extreme_idx]])

        old_yi = yi_current
        yi_current = yi_filled
        sei_current = sei_filled

        if len(yi_current) == len(old_yi) and np.allclose(yi_current, old_yi, atol=1e-10):
            break

    return meta_analysis(yi_current, sei_current, estimator, ci_method, conf_level)


def pet_peese(yi: np.ndarray, sei: np.ndarray, estimator: str = 'DL',
              ci_method: str = 'Wald', conf_level: float = 0.95) -> MetaResult:
    """Conditional PET-PEESE publication bias correction.

    P0-4 FIX: Uses t-distribution (df=k-2) instead of z for inference.
    """
    k = len(yi)
    if k < 3:
        return meta_analysis(yi, sei, estimator, ci_method, conf_level)

    wi = 1.0 / (sei ** 2)
    df = max(1, k - 2)

    # PET: WLS regression of yi on sei, weighted by 1/vi
    pet_result = _weighted_regression(yi, sei, wi, df, use_se_squared=False)
    pet_intercept, pet_se, pet_p = pet_result

    if pet_p >= 0.05:
        intercept, se_int, p_val = pet_result
    else:
        intercept, se_int, p_val = _weighted_regression(yi, sei, wi, df, use_se_squared=True)

    # Construct result (tau2 from base estimator for I²/Q reporting)
    base = meta_analysis(yi, sei, estimator, ci_method, conf_level)

    from src.utils import t_quantile
    alpha = 1.0 - conf_level
    t_crit = t_quantile(1.0 - alpha / 2, df)
    ci_lo = intercept - t_crit * se_int
    ci_hi = intercept + t_crit * se_int

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


def _weighted_regression(yi, sei, wi, df, use_se_squared=False):
    """WLS regression of yi on sei (or sei²), return (intercept, SE, p-value).

    P0-4 FIX: Uses t-distribution with df degrees of freedom.
    P0-5 FIX: Uses relative threshold for near-singular detection.
    """
    from src.utils import t_cdf

    x = sei ** 2 if use_se_squared else sei
    n = len(yi)

    # WLS: minimize sum(wi * (yi - a - b*x)^2)
    sw = np.sum(wi)
    sx = np.sum(wi * x)
    sy = np.sum(wi * yi)
    sxx = np.sum(wi * x ** 2)
    sxy = np.sum(wi * x * yi)

    denom = sw * sxx - sx ** 2

    # P0-5 FIX: Relative threshold for near-singular detection
    scale_ref = max(abs(sw * sxx), abs(sx ** 2), 1e-30)
    if abs(denom) < 1e-10 * scale_ref:
        # Degenerate: return simple weighted mean
        intercept = float(sy / sw) if sw > 0 else 0.0
        se = float(1.0 / math.sqrt(sw)) if sw > 0 else 1.0
        t_stat = intercept / se if se > 0 else 0
        p = 2.0 * (1.0 - t_cdf(abs(t_stat), max(1, df)))
        return intercept, se, p

    intercept = float((sxx * sy - sx * sxy) / denom)
    slope = float((sw * sxy - sx * sy) / denom)

    # Residual variance (sigma²)
    residuals = yi - intercept - slope * x
    sigma2 = float(np.sum(wi * residuals ** 2) / max(n - 2, 1))

    # SE of intercept
    var_intercept = sigma2 * sxx / denom
    se = math.sqrt(max(0, float(var_intercept)))

    t_stat = intercept / se if se > 0 else 0
    p = 2.0 * (1.0 - t_cdf(abs(t_stat), max(1, df)))

    return intercept, se, p
