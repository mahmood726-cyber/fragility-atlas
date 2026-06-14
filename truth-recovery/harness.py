"""Truth-recovery harness for the Fragility Atlas multiverse classifier.

Wires the repo's OWN functions:
    src.specifications.generate_specifications  (multiverse over 7x3x3x(k+1) specs)
    src.classifier.classify_review              (robustness_score + Robust/Fragile)

against a seeded known-truth corpus (truth-recovery/dgp.py) and measures whether
the "robustness" verdict tracks reality.

Key truth-recovery questions for an *atlas/multiverse* robustness metric:

  Q1 (over-detection / false robustness):
      Among KNOWN-NULL reviews that happen to look significant, how often does
      the classifier stamp them "Robust"? A high rate means the relative
      agreement threshold manufactures robustness on no-signal data.

  Q2 (discrimination):
      Does robustness_score separate true-effect reviews from null reviews?
      Report AUROC. An uncalibrated relative score (% of specs agreeing with a
      reference spec) can be high for BOTH null and effect reviews, giving
      AUROC ~ 0.5 -- i.e. the headline metric carries little truth signal once
      you condition on significance.
"""

import numpy as np

from src.classifier import classify_review
from src.specifications import generate_specifications


def auroc(scores, labels):
    """Rank-based AUROC (Mann-Whitney). scores: array; labels: 0/1 array."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    _assign_tie_ranks(scores, ranks)
    sum_pos = ranks[labels == 1].sum()
    n_pos, n_neg = len(pos), len(neg)
    return float((sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def _assign_tie_ranks(scores, ranks):
    order = np.argsort(scores, kind="mergesort")
    s = scores[order]
    i = 0
    n = len(s)
    while i < n:
        j = i
        while j + 1 < n and s[j + 1] == s[i]:
            j += 1
        if j > i:
            avg = (ranks[order[i:j + 1]]).mean()
            ranks[order[i:j + 1]] = avg
        i = j + 1


def evaluate(corpus):
    """Run repo's multiverse + classifier on each review; collect outcomes."""
    rows = []
    for review, truth in corpus:
        specs = generate_specifications(review)
        cls = classify_review(review, specs)
        rows.append({
            "truth": truth,
            "ref_significant": bool(cls.cochrane_significant),
            "robustness": float(cls.robustness_score),
            "classification": cls.classification,
            "k": review.k,
        })
    return rows


def summarize(rows):
    """Compute truth-recovery metrics from harness rows."""
    truth = np.array([r["truth"] for r in rows])
    rob = np.array([r["robustness"] for r in rows])
    sig = np.array([r["ref_significant"] for r in rows])
    cls = np.array([r["classification"] for r in rows])

    is_robust = (cls == "Robust")

    # Q1: false robustness on KNOWN-NULL reviews.
    null_mask = truth == 0
    null_sig = null_mask & sig
    null_sig_robust = null_sig & is_robust
    n_null_sig = int(null_sig.sum())
    fpr_robust_among_sig_null = (
        float(null_sig_robust.sum()) / n_null_sig if n_null_sig else float("nan")
    )
    # robust-AND-significant on null, as a fraction of ALL null reviews
    false_robust_rate_all_null = float(null_sig_robust.sum()) / int(null_mask.sum())

    # Q2: discrimination of robustness score (true effect vs null), full corpus.
    auc_full = auroc(rob, truth)

    # Q2b: discrimination CONDITIONAL on being significant (the honest test --
    # once a review is significant, does robustness tell true from false +ve?).
    sig_mask = sig
    auc_among_sig = (
        auroc(rob[sig_mask], truth[sig_mask])
        if sig_mask.sum() > 0 and len(np.unique(truth[sig_mask])) == 2
        else float("nan")
    )

    # Mean robustness by truth class (are nulls just as "robust"?).
    mean_rob_null = float(rob[null_mask].mean())
    mean_rob_eff = float(rob[truth == 1].mean())

    # Mean robustness among SIGNIFICANT nulls vs SIGNIFICANT effects.
    eff_sig = (truth == 1) & sig
    mean_rob_null_sig = float(rob[null_sig].mean()) if null_sig.sum() else float("nan")
    mean_rob_eff_sig = float(rob[eff_sig].mean()) if eff_sig.sum() else float("nan")

    return {
        "n_reviews": len(rows),
        "n_null_significant": n_null_sig,
        "n_eff_significant": int(eff_sig.sum()),
        "false_robust_among_sig_null": fpr_robust_among_sig_null,
        "false_robust_rate_all_null": false_robust_rate_all_null,
        "auroc_robustness_full": auc_full,
        "auroc_robustness_among_sig": auc_among_sig,
        "mean_rob_null": mean_rob_null,
        "mean_rob_eff": mean_rob_eff,
        "mean_rob_null_sig": mean_rob_null_sig,
        "mean_rob_eff_sig": mean_rob_eff_sig,
    }


def run_scenario(n_each, seed=2026, effect_mu=0.35, tau_null=0.0,
                 tau_effect=0.0, n_scale=200):
    from dgp import make_corpus
    corpus = make_corpus(seed=seed, n_each=n_each, effect_mu=effect_mu,
                         tau_null=tau_null, tau_effect=tau_effect)
    # n_scale is applied inside make_review via make_corpus defaults; to vary
    # precision we regenerate with a custom scale.
    if n_scale != 200:
        import numpy as _np
        from dgp import make_review
        rng = _np.random.default_rng(seed)
        corpus = []
        for i in range(n_each):
            k = int(rng.integers(5, 26))
            corpus.append((make_review(rng, f"NULL{i:04d}", 0.0, tau_null, k,
                                       n_scale=n_scale), 0))
        for i in range(n_each):
            k = int(rng.integers(5, 26))
            corpus.append((make_review(rng, f"EFF{i:04d}", effect_mu, tau_effect,
                                       k, n_scale=n_scale), 1))
    return summarize(evaluate(corpus))


if __name__ == "__main__":
    import sys
    from dgp import make_corpus

    n_each = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    # Stress regime: low precision + heterogeneity so nulls CAN cross
    # significance -- the real test of false-robustness over-detection.
    n_scale = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    tau = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    corpus = make_corpus(seed=2026, n_each=n_each, tau_null=tau, tau_effect=tau)
    if n_scale != 200:
        import numpy as _np
        from dgp import make_review
        rng = _np.random.default_rng(2026)
        corpus = []
        for i in range(n_each):
            k = int(rng.integers(5, 26))
            corpus.append((make_review(rng, f"NULL{i:04d}", 0.0, tau, k,
                                       n_scale=n_scale), 0))
        for i in range(n_each):
            k = int(rng.integers(5, 26))
            corpus.append((make_review(rng, f"EFF{i:04d}", 0.35, tau, k,
                                       n_scale=n_scale), 1))
    rows = evaluate(corpus)
    m = summarize(rows)

    print("=" * 64)
    print("FRAGILITY ATLAS -- TRUTH-RECOVERY HARNESS")
    print("=" * 64)
    print(f"Reviews (null+effect)          : {m['n_reviews']}")
    print(f"Significant NULL reviews       : {m['n_null_significant']}")
    print(f"Significant EFFECT reviews     : {m['n_eff_significant']}")
    print("-" * 64)
    print("Q1  False robustness (over-detection on no-signal data)")
    print(f"  P(verdict=Robust | null & significant) : "
          f"{m['false_robust_among_sig_null']:.3f}")
    print(f"  Robust&Sig among ALL null reviews      : "
          f"{m['false_robust_rate_all_null']:.3f}")
    print("-" * 64)
    print("Q2  Does robustness_score track TRUTH?")
    print(f"  AUROC robustness (effect vs null), full corpus : "
          f"{m['auroc_robustness_full']:.3f}")
    print(f"  AUROC robustness, among SIGNIFICANT only       : "
          f"{m['auroc_robustness_among_sig']:.3f}")
    print("-" * 64)
    print("Mean robustness_score by class")
    print(f"  null   (all)        : {m['mean_rob_null']:.1f}")
    print(f"  effect (all)        : {m['mean_rob_eff']:.1f}")
    print(f"  null   (significant): {m['mean_rob_null_sig']:.1f}")
    print(f"  effect (significant): {m['mean_rob_eff_sig']:.1f}")
    print("=" * 64)
