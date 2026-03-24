"""Robustness classification and dimension attribution for the Fragility Atlas."""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict
from src.specifications import SpecResult
from src.loader import ReviewData


@dataclass
class ReviewClassification:
    """Robustness classification for a single Cochrane review."""
    review_id: str
    review_doi: str
    analysis_name: str
    k: int
    scale: str
    total_specs: int
    agreeing_specs: int
    robustness_score: float        # 0-100%
    classification: str            # Robust/Moderate/Fragile/Unstable
    cochrane_significant: bool
    cochrane_direction: int
    # Dimension attribution (eta² values)
    eta2_estimator: float
    eta2_ci_method: float
    eta2_bias_correction: float
    eta2_leave_out: float
    # Summary stats across specs
    median_theta: float
    iqr_theta: float
    frac_significant: float
    frac_reversed: float


def classify_review(review: ReviewData, specs: List[SpecResult]) -> ReviewClassification:
    """Classify a review based on its specification results."""
    if not specs:
        return _empty_classification(review)

    total = len(specs)

    # Determine the "reference" conclusion from the full-set DL-Wald specification
    ref_spec = _find_reference_spec(specs)
    if ref_spec is None:
        ref_direction = 1 if np.mean([s.theta for s in specs]) >= 0 else -1
        ref_significant = review.is_significant
    else:
        ref_direction = ref_spec.direction
        ref_significant = ref_spec.is_significant

    # Count agreements: spec agrees if same significance AND same direction
    agreeing = 0
    reversed_count = 0
    for s in specs:
        same_direction = (s.direction == ref_direction)
        same_significance = (s.is_significant == ref_significant)
        if same_direction and same_significance:
            agreeing += 1
        if s.is_significant and (s.direction != ref_direction):
            reversed_count += 1

    robustness = (agreeing / total * 100) if total > 0 else 0

    # Classification
    if robustness >= 90:
        classification = 'Robust'
    elif robustness >= 70:
        classification = 'Moderate'
    elif robustness >= 50:
        classification = 'Fragile'
    else:
        classification = 'Unstable'

    # Dimension attribution (eta²)
    agreement_vector = np.array([
        1 if (s.direction == ref_direction and s.is_significant == ref_significant) else 0
        for s in specs
    ], dtype=float)

    eta2 = _compute_eta2(specs, agreement_vector)

    # Summary stats
    thetas = np.array([s.theta for s in specs])
    sig_frac = np.mean([s.is_significant for s in specs])

    return ReviewClassification(
        review_id=review.review_id,
        review_doi=review.review_doi,
        analysis_name=review.analysis_name,
        k=review.k,
        scale=review.scale,
        total_specs=total,
        agreeing_specs=agreeing,
        robustness_score=round(robustness, 2),
        classification=classification,
        cochrane_significant=ref_significant,
        cochrane_direction=ref_direction,
        eta2_estimator=eta2.get('estimator', 0),
        eta2_ci_method=eta2.get('ci_method', 0),
        eta2_bias_correction=eta2.get('bias_correction', 0),
        eta2_leave_out=eta2.get('leave_out', 0),
        median_theta=float(np.median(thetas)),
        iqr_theta=float(np.percentile(thetas, 75) - np.percentile(thetas, 25)),
        frac_significant=float(sig_frac),
        frac_reversed=float(reversed_count / total) if total > 0 else 0,
    )


def _find_reference_spec(specs: List[SpecResult]):
    """Find the DL + Wald + no-correction + full-set specification as reference."""
    for s in specs:
        if (s.estimator == 'DL' and s.ci_method == 'Wald' and
                s.bias_correction == 'none' and s.leave_out == ''):
            return s
    # Fallback: any full-set spec
    for s in specs:
        if s.leave_out == '':
            return s
    return None


def _compute_eta2(specs: List[SpecResult], agreement: np.ndarray) -> Dict[str, float]:
    """Compute eta² (variance explained) for each specification dimension."""
    total_var = np.var(agreement)
    if total_var < 1e-15:
        return {'estimator': 0, 'ci_method': 0, 'bias_correction': 0, 'leave_out': 0}

    result = {}
    for dim_name, dim_getter in [
        ('estimator', lambda s: s.estimator),
        ('ci_method', lambda s: s.ci_method),
        ('bias_correction', lambda s: s.bias_correction),
        ('leave_out', lambda s: 'full' if s.leave_out == '' else 'loo'),
    ]:
        groups = {}
        for i, s in enumerate(specs):
            key = dim_getter(s)
            if key not in groups:
                groups[key] = []
            groups[key].append(agreement[i])

        # Between-group sum of squares
        grand_mean = np.mean(agreement)
        ss_between = sum(
            len(vals) * (np.mean(vals) - grand_mean) ** 2
            for vals in groups.values()
        )
        ss_total = total_var * len(agreement)

        result[dim_name] = round(ss_between / ss_total, 4) if ss_total > 0 else 0

    return result


def _empty_classification(review):
    return ReviewClassification(
        review_id=review.review_id,
        review_doi=review.review_doi,
        analysis_name=review.analysis_name,
        k=review.k,
        scale=review.scale,
        total_specs=0,
        agreeing_specs=0,
        robustness_score=0,
        classification='Insufficient',
        cochrane_significant=review.is_significant,
        cochrane_direction=0,
        eta2_estimator=0,
        eta2_ci_method=0,
        eta2_bias_correction=0,
        eta2_leave_out=0,
        median_theta=0,
        iqr_theta=0,
        frac_significant=0,
        frac_reversed=0,
    )
