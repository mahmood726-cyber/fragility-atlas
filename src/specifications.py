"""Specification grid generator and executor for multiverse meta-analysis."""

import numpy as np
from itertools import product
from dataclasses import dataclass
from typing import List
from src.estimators import meta_analysis, MetaResult
from src.corrections import trim_and_fill, pet_peese
from src.loader import ReviewData


ESTIMATORS = ['FE', 'DL', 'REML', 'PM', 'SJ', 'HS', 'HE']
CI_METHODS = ['Wald', 'HKSJ', 't-dist']
BIAS_CORRECTIONS = ['none', 'trim-and-fill', 'PET-PEESE']


@dataclass
class SpecResult:
    """Result of a single specification within the multiverse."""
    review_id: str
    estimator: str
    ci_method: str
    bias_correction: str
    leave_out: str          # '' for full set, study label for LOO
    theta: float
    se_theta: float
    ci_lo: float
    ci_hi: float
    p_value: float
    tau2: float
    i2: float
    is_significant: bool    # p < 0.05
    direction: int          # +1 or -1


def generate_specifications(review: ReviewData, conf_level: float = 0.95) -> List[SpecResult]:
    """Generate and execute all specifications for a single review.

    Dimensions:
    - 7 estimators × 3 CI methods × 3 bias corrections × (k+1) LOO = 63·(k+1) specs
    """
    results = []

    # Build LOO datasets: full set + k leave-one-out subsets
    loo_sets = [('', review.yi, review.sei)]  # full set
    for i in range(review.k):
        mask = np.ones(review.k, dtype=bool)
        mask[i] = False
        label = review.study_labels[i]
        loo_sets.append((label, review.yi[mask], review.sei[mask]))

    for (loo_label, yi, sei) in loo_sets:
        if len(yi) < 2:
            continue

        for estimator, ci_method, bias_corr in product(ESTIMATORS, CI_METHODS, BIAS_CORRECTIONS):
            try:
                result = _run_specification(yi, sei, estimator, ci_method,
                                            bias_corr, conf_level)
                if result is None:
                    continue

                is_sig = result.p_value < 0.05
                direction = 1 if result.theta >= 0 else -1

                results.append(SpecResult(
                    review_id=review.review_id,
                    estimator=estimator,
                    ci_method=ci_method,
                    bias_correction=bias_corr,
                    leave_out=loo_label,
                    theta=result.theta,
                    se_theta=result.se_theta,
                    ci_lo=result.ci_lo,
                    ci_hi=result.ci_hi,
                    p_value=result.p_value,
                    tau2=result.tau2,
                    i2=result.i2,
                    is_significant=is_sig,
                    direction=direction,
                ))
            except (ValueError, FloatingPointError, ZeroDivisionError):
                continue
            except Exception as e:
                import sys
                print(f"Warning: spec {estimator}/{ci_method}/{bias_corr} "
                      f"failed for {review.review_id}: {e}", file=sys.stderr)
                continue

    return results


def _run_specification(yi, sei, estimator, ci_method, bias_corr, conf_level):
    """Execute a single specification."""
    if bias_corr == 'none':
        return meta_analysis(yi, sei, estimator, ci_method, conf_level)
    elif bias_corr == 'trim-and-fill':
        return trim_and_fill(yi, sei, estimator, ci_method, conf_level)
    elif bias_corr == 'PET-PEESE':
        return pet_peese(yi, sei, estimator, ci_method, conf_level)
    return None
