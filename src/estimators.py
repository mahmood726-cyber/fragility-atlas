"""Meta-analysis estimators: 7 tau² methods + 3 CI methods.

All estimators take yi (effect sizes) and sei (standard errors) as numpy arrays
and return a MetaResult dataclass.
"""

import math
import numpy as np
from dataclasses import dataclass
from src.utils import normal_quantile, t_quantile


@dataclass
class MetaResult:
    """Result of a single meta-analysis specification."""
    theta: float          # pooled effect estimate
    se_theta: float       # standard error of pooled estimate
    ci_lo: float          # lower confidence interval bound
    ci_hi: float          # upper confidence interval bound
    p_value: float        # two-sided p-value
    tau2: float           # between-study variance estimate
    i2: float             # I² heterogeneity (%)
    q_stat: float         # Cochran's Q statistic
    estimator: str        # name of tau² estimator
    ci_method: str        # name of CI method
    k: int                # number of studies


def meta_analysis(yi: np.ndarray, sei: np.ndarray, estimator: str = 'DL',
                  ci_method: str = 'Wald', conf_level: float = 0.95) -> MetaResult:
    """Run a meta-analysis with specified estimator and CI method."""
    k = len(yi)
    wi = 1.0 / (sei ** 2)

    # Step 1: Estimate tau²
    tau2 = _estimate_tau2(yi, sei, wi, estimator)
    tau2 = max(0.0, tau2)  # truncate at 0

    # Step 2: Compute pooled estimate with random-effects weights
    wi_star = 1.0 / (sei ** 2 + tau2)
    theta = np.sum(wi_star * yi) / np.sum(wi_star)
    se_theta = 1.0 / math.sqrt(np.sum(wi_star))

    # Step 3: Compute Q statistic (always using FE weights)
    theta_fe = np.sum(wi * yi) / np.sum(wi)
    q_stat = float(np.sum(wi * (yi - theta_fe) ** 2))

    # Step 4: Compute I²
    i2 = max(0.0, (q_stat - (k - 1)) / q_stat * 100) if q_stat > 0 else 0.0

    # Step 5: Compute CI using specified method
    alpha = 1.0 - conf_level
    ci_lo, ci_hi, p_value = _compute_ci(theta, se_theta, yi, sei, tau2, wi_star,
                                         k, alpha, ci_method)

    return MetaResult(
        theta=float(theta),
        se_theta=float(se_theta),
        ci_lo=float(ci_lo),
        ci_hi=float(ci_hi),
        p_value=float(p_value),
        tau2=float(tau2),
        i2=float(i2),
        q_stat=float(q_stat),
        estimator=estimator,
        ci_method=ci_method,
        k=k,
    )


def _estimate_tau2(yi, sei, wi, estimator):
    """Estimate tau² using the specified method."""
    k = len(yi)

    if estimator == 'FE':
        return 0.0

    elif estimator == 'DL':
        return _dl_tau2(yi, wi, k)

    elif estimator == 'REML':
        return _reml_tau2(yi, sei, k)

    elif estimator == 'PM':
        return _pm_tau2(yi, sei, k)

    elif estimator == 'SJ':
        return _sj_tau2(yi, sei, k)

    elif estimator == 'HS':
        return _hs_tau2(yi, sei, k)

    elif estimator == 'HE':
        return _he_tau2(yi, sei, k)

    else:
        raise ValueError(f"Unknown estimator: {estimator}")


def _dl_tau2(yi, wi, k):
    """DerSimonian-Laird moment estimator (1986)."""
    theta_fe = np.sum(wi * yi) / np.sum(wi)
    Q = float(np.sum(wi * (yi - theta_fe) ** 2))
    C = float(np.sum(wi) - np.sum(wi ** 2) / np.sum(wi))
    if C <= 0:
        return 0.0
    return max(0.0, (Q - (k - 1)) / C)


def _reml_tau2(yi, sei, k, max_iter=100, tol=1e-8):
    """Restricted Maximum Likelihood via Fisher scoring (Viechtbauer 2005)."""
    vi = sei ** 2
    # Initial estimate: DL
    wi = 1.0 / vi
    tau2 = _dl_tau2(yi, wi, k)
    tau2 = max(0.0, tau2)

    for _ in range(max_iter):
        wi_star = 1.0 / (vi + tau2)
        theta = np.sum(wi_star * yi) / np.sum(wi_star)
        resid = yi - theta

        # REML log-likelihood gradient
        # d(ll)/d(tau2) = -0.5 * sum(wi_star) + 0.5 * sum(wi_star^2 * resid^2)
        #                 + 0.5 * sum(wi_star^2) / sum(wi_star)  [REML correction]
        gradient = (-0.5 * np.sum(wi_star)
                    + 0.5 * np.sum(wi_star ** 2 * resid ** 2)
                    + 0.5 * np.sum(wi_star ** 2) / np.sum(wi_star))

        # Fisher information
        fisher_info = 0.5 * (np.sum(wi_star ** 2)
                             - 2.0 * np.sum(wi_star ** 3) / np.sum(wi_star)
                             + (np.sum(wi_star ** 2) / np.sum(wi_star)) ** 2)

        if fisher_info <= 0:
            break

        tau2_new = tau2 + gradient / fisher_info
        tau2_new = max(0.0, tau2_new)

        if abs(tau2_new - tau2) < tol:
            tau2 = tau2_new
            break
        tau2 = tau2_new

    return tau2


def _pm_tau2(yi, sei, k, max_iter=100, tol=1e-8):
    """Paule-Mandel iterative estimator (1982)."""
    vi = sei ** 2
    # Initial guess: DL
    wi = 1.0 / vi
    tau2 = _dl_tau2(yi, wi, k)
    tau2 = max(0.0, tau2)

    for _ in range(max_iter):
        wi_star = 1.0 / (vi + tau2)
        theta = np.sum(wi_star * yi) / np.sum(wi_star)
        Q_star = float(np.sum(wi_star * (yi - theta) ** 2))

        # PM criterion: Q*(tau2) = k - 1
        if abs(Q_star - (k - 1)) < tol:
            break

        # Update using bisection-like step
        C = float(np.sum(wi_star) - np.sum(wi_star ** 2) / np.sum(wi_star))
        if C <= 0:
            break
        tau2_new = tau2 + (Q_star - (k - 1)) / C
        tau2_new = max(0.0, tau2_new)

        if abs(tau2_new - tau2) < tol:
            tau2 = tau2_new
            break
        tau2 = tau2_new

    return tau2


def _sj_tau2(yi, sei, k):
    """Sidik-Jonkman moment estimator (2005)."""
    if k < 2:
        return 0.0
    vi = sei ** 2
    # Unweighted mean
    theta_uw = np.mean(yi)
    # Initial tau2 from residuals
    tau2_0 = np.sum((yi - theta_uw) ** 2) / (k - 1)

    # One-step update
    wi = 1.0 / (vi + tau2_0)
    theta = np.sum(wi * yi) / np.sum(wi)
    Q = float(np.sum(wi * (yi - theta) ** 2))
    C = float(np.sum(wi) - np.sum(wi ** 2) / np.sum(wi))

    if C <= 0:
        return max(0.0, tau2_0)

    return max(0.0, (Q - (k - 1)) / C)


def _hs_tau2(yi, sei, k):
    """Hunter-Schmidt variance components estimator."""
    if k < 2:
        return 0.0
    vi = sei ** 2
    # Unweighted (sample-size weighted in original HS, but we use equal weights)
    theta = np.mean(yi)
    S2 = np.sum((yi - theta) ** 2) / (k - 1)
    v_bar = np.mean(vi)
    return max(0.0, S2 - v_bar)


def _he_tau2(yi, sei, k):
    """Hedges unweighted estimator."""
    if k < 2:
        return 0.0
    vi = sei ** 2
    theta = np.mean(yi)
    Q_uw = np.sum((yi - theta) ** 2)
    return max(0.0, (Q_uw - np.sum(vi)) / (k - 1))


def _compute_ci(theta, se_theta, yi, sei, tau2, wi_star, k, alpha, ci_method):
    """Compute confidence interval and p-value using specified method."""
    from src.utils import normal_cdf

    if ci_method == 'Wald':
        z_crit = normal_quantile(1.0 - alpha / 2)
        ci_lo = theta - z_crit * se_theta
        ci_hi = theta + z_crit * se_theta
        z = theta / se_theta if se_theta > 0 else 0
        p_value = 2.0 * (1.0 - normal_cdf(abs(z)))

    elif ci_method == 'HKSJ':
        # Hartung-Knapp-Sidik-Jonkman adjustment
        if k < 2:
            return theta, theta, 1.0

        resid = yi - theta
        q_hat = float(np.sum(wi_star * resid ** 2)) / (k - 1)
        q_hat = max(q_hat, 1.0)  # HKSJ: bound q_hat >= 1 to avoid anti-conservative CI

        se_adj = se_theta * math.sqrt(q_hat)
        t_crit = t_quantile(1.0 - alpha / 2, k - 1)
        ci_lo = theta - t_crit * se_adj
        ci_hi = theta + t_crit * se_adj

        from src.utils import t_cdf as _t_cdf
        t_stat = theta / se_adj if se_adj > 0 else 0
        p_value = 2.0 * (1.0 - _t_cdf(abs(t_stat), k - 1))

    elif ci_method == 't-dist':
        # Knapp-Hartung: use t-distribution with k-1 df, no variance adjustment
        t_crit = t_quantile(1.0 - alpha / 2, k - 1)
        ci_lo = theta - t_crit * se_theta
        ci_hi = theta + t_crit * se_theta

        from src.utils import t_cdf as _t_cdf
        t_stat = theta / se_theta if se_theta > 0 else 0
        p_value = 2.0 * (1.0 - _t_cdf(abs(t_stat), k - 1))

    else:
        raise ValueError(f"Unknown CI method: {ci_method}")

    return ci_lo, ci_hi, p_value
