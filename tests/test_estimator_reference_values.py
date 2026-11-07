"""Regression checks against external metafor reference results.

The heterogeneous eight-study reference is recorded in AllMeta's
tests/test_ma_core_tau2.py (metafor 4.6.0, 2026-06-02). The CLASP example is
published by the metafor author at:
https://www.metafor-project.org/doku.php/analyses:dersimonian2007

These tests check numerical agreement, rather than just positive variances
or the direction of the pooled result. Published CLASP values have four
decimal places, so their tolerance reflects that rounding.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from src.estimators import meta_analysis


REML_CASES = json.loads(
    (Path(__file__).parent / 'fixtures' / 'reml_nonconvergence.json').read_text()
)['reviews']


@pytest.mark.parametrize('case', REML_CASES, ids=lambda case: case['review_id'])
def test_reml_nonconvergence_uses_profile_likelihood_solution(case, capsys):
    result = meta_analysis(np.array(case['yi']), np.sqrt(case['vi']), 'REML')
    assert result.tau2 == pytest.approx(case['expected_tau2'], abs=2e-7)
    assert 'did not converge' not in capsys.readouterr().err


@pytest.mark.parametrize("method, expected", [
    ("DL", 0.05450642),
    ("HE", 0.05538393),
    ("HS", 0.04368853),
    ("SJ", 0.05989549),
    ("REML", 0.05432882),
    ("PM", 0.05474718),
])
def test_heterogeneous_tau2_matches_metafor(method, expected):
    yi = np.array([0.20, 0.50, -0.10, 0.80, 0.30, 0.60, 0.10, 0.45])
    vi = np.array([0.020, 0.030, 0.015, 0.050, 0.025, 0.040, 0.018, 0.030])
    result = meta_analysis(yi, np.sqrt(vi), method, "Wald")
    assert result.tau2 == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize("method, expected", [
    ("PM", (0.3681, -0.3811, 0.2060)),
    ("HE", (0.4410, -0.4035, 0.2327)),
    ("DL", (0.2323, -0.3240, 0.1540)),
    ("REML", (0.1843, -0.2980, 0.1343)),
    ("HS", (0.1330, -0.2666, 0.1125)),
    ("SJ", (0.4572, -0.4079, 0.2386)),
])
def test_clasp_example_matches_published_metafor_result(method, expected):
    a = np.array([5, 5, 12, 69, 9, 313], dtype=float)
    c = np.array([8, 17, 9, 94, 11, 352], dtype=float)
    n1 = np.array([156, 303, 565, 1570, 103, 4659], dtype=float)
    n2 = np.array([74, 303, 477, 1565, 105, 4650], dtype=float)
    yi = np.log(a / (n1 - a)) - np.log(c / (n2 - c))
    vi = 1 / a + 1 / (n1 - a) + 1 / c + 1 / (n2 - c)
    result = meta_analysis(yi, np.sqrt(vi), method, "Wald")
    actual = (np.sqrt(result.tau2), result.theta, result.se_theta)
    assert actual == pytest.approx(expected, abs=5.1e-5)


@pytest.mark.parametrize("method", ["HS", "HE", "SJ"])
def test_corrected_estimators_handle_identical_effects(method):
    result = meta_analysis(np.full(4, 0.25), np.array([0.1, 0.2, 0.3, 0.4]),
                           method, "Wald")
    assert result.tau2 == 0.0
    assert result.theta == pytest.approx(0.25)


@pytest.mark.parametrize("method", ["HS", "HE", "SJ"])
def test_corrected_estimators_are_invariant_to_effect_direction(method):
    yi = np.array([-0.5, 0.3, 0.8, 1.2])
    sei = np.array([0.1, 0.2, 0.35, 0.5])
    forward = meta_analysis(yi, sei, method, "Wald")
    reverse = meta_analysis(-yi, sei, method, "Wald")
    assert reverse.tau2 == pytest.approx(forward.tau2, abs=1e-12)
    assert reverse.theta == pytest.approx(-forward.theta, abs=1e-12)
