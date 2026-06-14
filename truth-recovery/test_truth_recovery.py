"""Truth-recovery invariants for the Fragility Atlas multiverse classifier.

Feeds the repo's OWN generate_specifications + classify_review a seeded
known-truth corpus and asserts the robustness verdict behaves honestly:

  1. It does NOT manufacture robustness on no-signal data: among KNOWN-NULL
     reviews that look significant by chance, the "Robust" rate is ~0.
  2. CONDITIONAL on significance, robustness_score discriminates true effects
     from false positives (AUROC well above chance).
  3. Documented honest caveat: on the FULL corpus the score does NOT separate
     truth (AUROC << among-significant AUROC), because "robustly null" scores
     high too. The metric is only interpretable given the reference claim.

Two scenarios are each evaluated ONCE (module scope) and shared across asserts
to keep the multiverse cost bounded.

Run from repo root:
  PYTHONPATH=.:truth-recovery python -m pytest truth-recovery/test_truth_recovery.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from harness import run_scenario  # noqa: E402

# Evaluate each regime once (multiverse is expensive). Small n_each keeps the
# full file under the verify-time budget while staying statistically meaningful.
CLEAN = run_scenario(n_each=25, seed=2026, effect_mu=0.35, n_scale=200)
STRESS = run_scenario(n_each=40, seed=2026, effect_mu=0.35,
                      tau_null=0.3, tau_effect=0.3, n_scale=40)


def test_no_false_robustness_clean_regime():
    assert CLEAN["false_robust_rate_all_null"] == 0.0, CLEAN
    if CLEAN["n_null_significant"] > 0:
        assert CLEAN["false_robust_among_sig_null"] == 0.0, CLEAN


def test_stress_regime_produces_significant_nulls():
    # Low precision + heterogeneity must let some true-null reviews cross p<.05,
    # otherwise the over-detection test below is vacuous.
    assert STRESS["n_null_significant"] >= 1, STRESS


def test_no_false_robustness_stress_regime():
    assert STRESS["false_robust_among_sig_null"] == 0.0, STRESS


def test_robustness_discriminates_among_significant():
    assert STRESS["auroc_robustness_among_sig"] >= 0.75, STRESS


def test_significant_nulls_score_lower_than_significant_effects():
    assert STRESS["mean_rob_null_sig"] < STRESS["mean_rob_eff_sig"], STRESS


def test_full_corpus_auroc_documents_conflation_caveat():
    # Honest interpretability limit, not a pass/fail of the estimator: full-corpus
    # AUROC is materially weaker than the among-significant AUROC.
    assert STRESS["auroc_robustness_full"] < STRESS["auroc_robustness_among_sig"], STRESS
