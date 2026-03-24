"""Tests for meta-analysis estimators — internal consistency and mathematical properties.

Validates estimator formulas, CI methods, and edge cases using algebraically
verifiable properties rather than hardcoded reference values.
"""

import sys
import math
import numpy as np
import pytest

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent))
from src.estimators import meta_analysis


# Test dataset: 13 studies with known heterogeneity
TEST_YI = np.array([
    -0.8893, -1.5854, -1.3481, -1.4415, -0.2175, -0.7861,
    -1.6209, -0.5077, 0.0198, -0.4292, -0.3398, -1.4061, -0.3397
])
TEST_SEI = np.array([
    0.4659, 0.6840, 0.2827, 0.4225, 0.1355, 0.2067,
    0.2136, 0.2800, 0.2266, 0.2108, 0.2109, 0.3181, 0.3225
])


class TestDL:
    """DerSimonian-Laird estimator tests."""

    def test_theta_negative(self):
        r = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        assert r.theta < 0, f"DL theta={r.theta} should be negative"

    def test_tau2_positive(self):
        r = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        assert r.tau2 > 0, f"DL tau2={r.tau2} should be positive for heterogeneous data"

    def test_significant(self):
        r = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        assert r.p_value < 0.05

    def test_dl_formula_consistency(self):
        """Verify DL tau2 = max(0, (Q - (k-1)) / C)."""
        wi = 1.0 / TEST_SEI**2
        theta_fe = np.sum(wi * TEST_YI) / np.sum(wi)
        Q = float(np.sum(wi * (TEST_YI - theta_fe)**2))
        C = float(np.sum(wi) - np.sum(wi**2) / np.sum(wi))
        expected_tau2 = max(0, (Q - 12) / C)
        r = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        assert abs(r.tau2 - expected_tau2) < 1e-10


class TestREML:
    """REML estimator tests."""

    def test_theta_close_to_dl(self):
        """REML theta should be close to DL theta (within ~10%)."""
        r_dl = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        r_reml = meta_analysis(TEST_YI, TEST_SEI, 'REML', 'Wald')
        assert abs(r_reml.theta - r_dl.theta) < 0.1

    def test_tau2_nonnegative(self):
        r = meta_analysis(TEST_YI, TEST_SEI, 'REML', 'Wald')
        assert r.tau2 >= 0


class TestPM:
    """Paule-Mandel estimator tests."""

    def test_tau2_reasonable(self):
        """PM tau2 should be in the same ballpark as DL."""
        r_dl = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        r_pm = meta_analysis(TEST_YI, TEST_SEI, 'PM', 'Wald')
        assert abs(r_pm.tau2 - r_dl.tau2) < r_dl.tau2  # within 100% of DL


class TestFE:
    """Fixed-effect model tests."""

    def test_bcg_theta(self):
        r = meta_analysis(TEST_YI, TEST_SEI, 'FE', 'Wald')
        # FE theta should be inverse-variance weighted mean
        wi = 1.0 / TEST_SEI**2
        expected = np.sum(wi * TEST_YI) / np.sum(wi)
        assert abs(r.theta - expected) < 1e-10

    def test_tau2_zero(self):
        r = meta_analysis(TEST_YI, TEST_SEI, 'FE', 'Wald')
        assert r.tau2 == 0.0


class TestCIMethods:
    """CI method comparison tests."""

    def test_hksj_wider_than_wald(self):
        """HKSJ CI should generally be wider than Wald for heterogeneous data."""
        r_wald = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        r_hksj = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'HKSJ')
        wald_width = r_wald.ci_hi - r_wald.ci_lo
        hksj_width = r_hksj.ci_hi - r_hksj.ci_lo
        assert hksj_width >= wald_width * 0.9, \
            f"HKSJ width={hksj_width:.4f} vs Wald width={wald_width:.4f}"

    def test_tdist_wider_than_wald(self):
        """t-dist CI should be wider than Wald for small k."""
        r_wald = meta_analysis(TEST_YI, TEST_SEI, 'DL', 'Wald')
        r_tdist = meta_analysis(TEST_YI, TEST_SEI, 'DL', 't-dist')
        wald_width = r_wald.ci_hi - r_wald.ci_lo
        tdist_width = r_tdist.ci_hi - r_tdist.ci_lo
        assert tdist_width >= wald_width


class TestEdgeCases:
    """Edge case tests."""

    def test_k2_all_methods(self):
        """All methods should work with k=2 (minimum for tau2 estimation)."""
        yi = np.array([-0.5, 0.5])
        sei = np.array([0.3, 0.3])
        # k=2 is too small for our pipeline (min_k=3) but estimators should handle it
        for est in ['FE', 'DL', 'REML', 'PM', 'SJ', 'HS', 'HE']:
            r = meta_analysis(yi, sei, est, 'Wald')
            assert math.isfinite(r.theta), f"{est} produced non-finite theta"
            assert math.isfinite(r.p_value), f"{est} produced non-finite p-value"

    def test_homogeneous_data(self):
        """When all studies agree, tau2 should be ~0 and I2 ~0."""
        yi = np.array([-0.5, -0.5, -0.5, -0.5, -0.5])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        r = meta_analysis(yi, sei, 'DL', 'Wald')
        assert r.tau2 < 0.01, f"tau2={r.tau2} for homogeneous data"
        assert r.i2 < 5, f"I2={r.i2}% for homogeneous data"

    def test_highly_heterogeneous(self):
        """High heterogeneity should produce large tau2 and I2."""
        yi = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        r = meta_analysis(yi, sei, 'DL', 'Wald')
        assert r.tau2 > 1.0, f"tau2={r.tau2} for heterogeneous data"
        assert r.i2 > 90, f"I2={r.i2}% for heterogeneous data"

    def test_all_estimators_same_direction(self):
        """All estimators should produce same direction for clear effect."""
        for est in ['FE', 'DL', 'REML', 'PM', 'SJ', 'HS', 'HE']:
            r = meta_analysis(TEST_YI, TEST_SEI, est, 'Wald')
            assert r.theta < 0, f"{est} produced positive theta={r.theta} for BCG data"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
