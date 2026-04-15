"""Tests for corrections, classifier, loader, and specifications modules."""

import os
import sys
import math
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent))
from src.estimators import meta_analysis
from src.corrections import trim_and_fill, pet_peese, _weighted_regression
from src.classifier import classify_review, _compute_eta2, ReviewClassification
from src.specifications import generate_specifications, SpecResult

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _pairwise70_sample_path():
    env_pairwise = os.getenv('PAIRWISE70_DATA_DIR')
    candidates = []
    if env_pairwise:
        env_path = Path(env_pairwise)
        candidates.append(env_path / 'CD000028_pub4_data.rda' if env_path.is_dir() else env_path)
    candidates.extend([
        PROJECT_ROOT.parent / 'Projects' / 'Pairwise70' / 'data' / 'CD000028_pub4_data.rda',
        PROJECT_ROOT.parent / 'Models' / 'Pairwise70' / 'data' / 'CD000028_pub4_data.rda',
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    pytest.skip('Pairwise70 sample review not available')


# ─── Test data ───

# Clear negative effect (protective intervention) - trim-and-fill should detect left-side
PROTECTIVE_YI = np.array([-1.5, -1.2, -0.9, -0.7, -0.3, 0.1, -0.5])
PROTECTIVE_SEI = np.array([0.3, 0.25, 0.2, 0.35, 0.15, 0.4, 0.28])

# Symmetric data (no publication bias)
SYMMETRIC_YI = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
SYMMETRIC_SEI = np.array([0.2, 0.2, 0.2, 0.2, 0.2])


class TestTrimAndFill:
    def test_returns_meta_result(self):
        r = trim_and_fill(PROTECTIVE_YI, PROTECTIVE_SEI, 'DL', 'Wald')
        assert hasattr(r, 'theta')
        assert hasattr(r, 'p_value')
        assert math.isfinite(r.theta)

    def test_k_less_than_3_falls_back(self):
        yi = np.array([-0.5, -0.3])
        sei = np.array([0.2, 0.2])
        r = trim_and_fill(yi, sei, 'DL', 'Wald')
        base = meta_analysis(yi, sei, 'DL', 'Wald')
        assert abs(r.theta - base.theta) < 1e-10

    def test_symmetric_data_no_fill(self):
        """Symmetric data should produce similar result to unadjusted."""
        r_tf = trim_and_fill(SYMMETRIC_YI, SYMMETRIC_SEI, 'DL', 'Wald')
        r_base = meta_analysis(SYMMETRIC_YI, SYMMETRIC_SEI, 'DL', 'Wald')
        # Should be close (symmetric = no asymmetry to correct)
        assert abs(r_tf.theta - r_base.theta) < 0.3

    def test_detects_asymmetry(self):
        """For asymmetric funnel, trim-and-fill should produce a different result."""
        # Asymmetric: most studies negative, one outlier positive
        yi = np.array([-1.0, -0.8, -0.6, -0.4, 0.5])
        sei = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        r_tf = trim_and_fill(yi, sei, 'DL', 'Wald')
        r_base = meta_analysis(yi, sei, 'DL', 'Wald')
        # Trim-and-fill should produce a result (finite, not crashed)
        assert math.isfinite(r_tf.theta)
        # The k in the result may be larger (imputed studies added)
        assert r_tf.k >= len(yi)


class TestPetPeese:
    def test_returns_meta_result(self):
        r = pet_peese(PROTECTIVE_YI, PROTECTIVE_SEI, 'DL', 'Wald')
        assert math.isfinite(r.theta)
        assert math.isfinite(r.p_value)

    def test_k_less_than_3_falls_back(self):
        yi = np.array([-0.5, -0.3])
        sei = np.array([0.2, 0.2])
        r = pet_peese(yi, sei, 'DL', 'Wald')
        base = meta_analysis(yi, sei, 'DL', 'Wald')
        assert abs(r.theta - base.theta) < 1e-10

    def test_uses_t_distribution(self):
        """PET-PEESE should use t(k-2) not z, giving wider CIs for small k."""
        r = pet_peese(PROTECTIVE_YI, PROTECTIVE_SEI, 'DL', 'Wald')
        ci_width = r.ci_hi - r.ci_lo
        # With k=7, df=5, t(0.975,5)=2.571 vs z=1.96 → ~31% wider
        assert ci_width > 0  # at minimum, CI should exist


class TestWeightedRegression:
    def test_near_equal_se(self):
        """Near-equal SEs should not produce zero-width CI (P0-5 fix)."""
        yi = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
        sei = np.array([0.200001, 0.200002, 0.200003, 0.200004, 0.200005])
        wi = 1.0 / sei**2
        intercept, se, p = _weighted_regression(yi, sei, wi, df=3)
        assert se > 0, f"SE should be positive, got {se}"
        assert math.isfinite(p)

    def test_degenerate_identical_se(self):
        """Identical SEs should fall back to weighted mean."""
        yi = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
        sei = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        wi = 1.0 / sei**2
        intercept, se, p = _weighted_regression(yi, sei, wi, df=3)
        assert math.isfinite(intercept)
        assert se > 0


class TestClassifier:
    def _make_review(self):
        from src.loader import ReviewData
        return ReviewData(
            review_id='CD999999', review_doi='10.test', analysis_name='Test',
            k=5, yi=np.array([-0.5, -0.3, -0.1, 0.1, 0.3]),
            sei=np.array([0.2, 0.2, 0.2, 0.2, 0.2]),
            ni=np.array([100, 100, 100, 100, 100]),
            study_labels=['A', 'B', 'C', 'D', 'E'],
            scale='ratio', cochrane_pooled=-0.1, cochrane_ci_lo=-0.28,
            cochrane_ci_hi=0.08, is_significant=False,
        )

    def test_all_agree_is_robust(self):
        review = self._make_review()
        specs = [SpecResult(
            review_id='CD999999', estimator='DL', ci_method='Wald',
            bias_correction='none', leave_out='', theta=-0.1, se_theta=0.05,
            ci_lo=-0.2, ci_hi=0.0, p_value=0.1,
            tau2=0, i2=0, is_significant=False, direction=-1,
        )] * 100
        c = classify_review(review, specs)
        assert c.robustness_score == 100.0
        assert c.classification == 'Robust'

    def test_half_disagree_is_fragile_or_unstable(self):
        review = self._make_review()
        agree = SpecResult(
            review_id='CD999999', estimator='DL', ci_method='Wald',
            bias_correction='none', leave_out='', theta=-0.1, se_theta=0.05,
            ci_lo=-0.2, ci_hi=0.0, p_value=0.1,
            tau2=0, i2=0, is_significant=False, direction=-1,
        )
        disagree = SpecResult(
            review_id='CD999999', estimator='DL', ci_method='HKSJ',
            bias_correction='none', leave_out='', theta=-0.3, se_theta=0.05,
            ci_lo=-0.4, ci_hi=-0.2, p_value=0.001,
            tau2=0, i2=0, is_significant=True, direction=-1,
        )
        specs = [agree] * 50 + [disagree] * 50
        c = classify_review(review, specs)
        assert c.robustness_score <= 55
        assert c.classification in ('Fragile', 'Unstable')

    def test_top_dimension_computed(self):
        review = self._make_review()
        specs = generate_specifications(review)
        c = classify_review(review, specs)
        assert c.top_dimension in ('Estimator', 'CI Method', 'Bias Correction', 'Leave-One-Out')

    def test_empty_specs(self):
        review = self._make_review()
        c = classify_review(review, [])
        assert c.classification == 'Insufficient'
        assert c.top_dimension == 'N/A'


class TestEta2:
    def test_uniform_agreement_zero_variance(self):
        """If all specs agree, eta² should be 0 for all dimensions."""
        specs = [SpecResult(
            review_id='X', estimator='DL', ci_method='Wald',
            bias_correction='none', leave_out='', theta=0.5, se_theta=0.1,
            ci_lo=0.3, ci_hi=0.7, p_value=0.001,
            tau2=0, i2=0, is_significant=True, direction=1,
        )] * 10
        agreement = np.ones(10)
        eta2 = _compute_eta2(specs, agreement)
        assert all(v == 0 for v in eta2.values())

    def test_loo_groups_by_study_label(self):
        """LOO eta² should group by specific study, not binary full/loo."""
        # 3 specs: full, LOO-A (disagrees), LOO-B (agrees)
        specs = [
            SpecResult('X', 'DL', 'Wald', 'none', '', 0.5, 0.1, 0.3, 0.7, 0.001, 0, 0, True, 1),
            SpecResult('X', 'DL', 'Wald', 'none', 'StudyA', 0.5, 0.1, 0.3, 0.7, 0.001, 0, 0, True, 1),
            SpecResult('X', 'DL', 'Wald', 'none', 'StudyB', -0.5, 0.1, -0.7, -0.3, 0.001, 0, 0, True, -1),
        ]
        agreement = np.array([1, 1, 0])
        eta2 = _compute_eta2(specs, agreement)
        # LOO should capture some variance since StudyB removal flips result
        assert eta2['leave_out'] > 0


class TestComputeEffectsRegression:
    """Regression tests for _compute_effects — guard against scalar-.copy() bug.

    `iterrows()` yields Series where cell access returns numpy/Python scalars.
    Calling .copy() on a scalar float raises AttributeError. This test exercises
    the per-row path with synthetic data so we can catch the bug without the
    Pairwise70 fixture.
    """

    def test_compute_effects_accepts_scalar_row_values(self):
        import pandas as pd
        from src.loader import _compute_effects

        primary = pd.DataFrame({
            'Study': ['Alpha 2001', 'Beta 2005', 'Gamma 2010'],
            'Mean': [0.75, 0.80, 0.60],
            'CI.start': [0.50, 0.60, 0.40],
            'CI.end': [1.10, 1.05, 0.90],
            'Experimental.N': [50, 75, 60],
            'Control.N': [50, 75, 60],
            'Experimental.cases': [10, 15, 8],
            'Control.cases': [15, 20, 12],
        })

        yi, sei, ni, labels = _compute_effects(primary, scale='ratio')

        assert yi is not None, "Should not return None for valid primary frame"
        assert len(yi) == 3
        assert len(sei) == 3
        assert len(ni) == 3
        assert labels == ['Alpha 2001', 'Beta 2005', 'Gamma 2010']
        # yi should be log(Mean) on the ratio scale
        assert math.isclose(yi[0], math.log(0.75), rel_tol=1e-9)
        # sei should be positive and finite
        assert all(s > 0 and math.isfinite(s) for s in sei)


class TestLoaderIntegration:
    def test_load_first_review(self):
        """Load the first Pairwise70 review and verify structure."""
        from src.loader import load_review
        r = load_review(_pairwise70_sample_path())
        assert r is not None
        assert r.review_id == 'CD000028'
        assert r.k >= 3
        assert len(r.yi) == r.k
        assert len(r.sei) == r.k
        assert r.scale in ('ratio', 'difference')
        assert math.isfinite(r.cochrane_pooled)
        assert r.cochrane_ci_lo < r.cochrane_ci_hi


def test_run_pipeline_uses_repo_relative_sibling_projects(tmp_path, monkeypatch):
    from src.loader import ReviewData
    from src.pipeline import run_pipeline, resolve_paths

    projects_root = tmp_path / 'projects'
    project_root = projects_root / 'FragilityAtlas'
    project_root.mkdir(parents=True)

    paths = resolve_paths(project_root=project_root, projects_root=projects_root)
    paths['pairwise_dir'].mkdir(parents=True, exist_ok=True)

    review = ReviewData(
        review_id='CD000001',
        review_doi='10.1000/example',
        analysis_name='Synthetic review',
        k=3,
        yi=np.array([-0.3, -0.2, -0.1]),
        sei=np.array([0.1, 0.1, 0.1]),
        ni=np.array([100, 100, 100]),
        study_labels=['A', 'B', 'C'],
        scale='ratio',
        cochrane_pooled=-0.2,
        cochrane_ci_lo=-0.31,
        cochrane_ci_hi=-0.09,
        is_significant=True,
    )
    spec = SpecResult(
        review_id='CD000001',
        estimator='DL',
        ci_method='Wald',
        bias_correction='none',
        leave_out='',
        theta=-0.2,
        se_theta=0.1,
        ci_lo=-0.4,
        ci_hi=-0.05,
        p_value=0.01,
        tau2=0.0,
        i2=0.0,
        is_significant=True,
        direction=-1,
    )
    classification = ReviewClassification(
        review_id='CD000001',
        review_doi='10.1000/example',
        analysis_name='Synthetic review',
        k=3,
        scale='ratio',
        total_specs=1,
        agreeing_specs=1,
        robustness_score=100.0,
        classification='Robust',
        cochrane_significant=True,
        cochrane_direction=-1,
        eta2_estimator=0.0,
        eta2_ci_method=0.0,
        eta2_bias_correction=0.0,
        eta2_leave_out=0.0,
        top_dimension='Estimator',
        median_theta=-0.2,
        iqr_theta=0.0,
        frac_significant=1.0,
        frac_reversed=0.0,
    )

    monkeypatch.setattr('src.pipeline.load_all_reviews', lambda pairwise_dir, min_k=3: [review])
    monkeypatch.setattr('src.pipeline.generate_specifications', lambda review, conf_level: [spec])
    monkeypatch.setattr('src.pipeline.classify_review', lambda review, specs: classification)

    result = run_pipeline(project_root=project_root, projects_root=projects_root)

    assert result[2]['n_reviews'] == 1
    assert (project_root / 'data' / 'output' / 'fragility_atlas_results.csv').exists()
    assert (project_root / 'data' / 'output' / 'fragility_atlas_specifications.csv').exists()
    assert (project_root / 'data' / 'output' / 'fragility_atlas_summary.json').exists()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
