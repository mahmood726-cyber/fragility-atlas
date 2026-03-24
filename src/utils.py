"""Math primitives for meta-analysis computations.

Uses scipy for statistical distributions (reliable, well-tested).
"""

import math
from scipy import stats


def normal_cdf(x):
    """Standard normal CDF."""
    return float(stats.norm.cdf(x))


def normal_quantile(p):
    """Inverse normal CDF."""
    return float(stats.norm.ppf(p))


def t_cdf(t_val, df):
    """Student's t CDF."""
    return float(stats.t.cdf(t_val, df))


def t_quantile(p, df):
    """Inverse t-distribution CDF."""
    return float(stats.t.ppf(p, df))


def chi2_quantile(p, df):
    """Chi-squared quantile."""
    return float(stats.chi2.ppf(p, df))
