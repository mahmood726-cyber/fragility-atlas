"""Seeded known-truth data-generating process for the Fragility Atlas.

Builds synthetic ReviewData objects with a KNOWN ground truth so the repo's own
multiverse robustness classifier can be measured against reality.

Two review classes:
  * NULL reviews  : true pooled effect = 0 (theta=0, optional heterogeneity).
  * EFFECT reviews: true pooled effect != 0 (a genuine signal).

Each review is a set of k studies on the log/difference scale (yi, sei). We
generate study-level effects from a random-effects model:

    theta_i ~ Normal(mu, tau^2)
    yhat_i  ~ Normal(theta_i, sei_i^2)

sei_i is drawn from realistic Cochrane-like precision (a few large, many small
trials). No RDA files, no pyreadr -- this is a standalone generator that emits
the exact `ReviewData` dataclass the repo's specification/classifier code
consumes.
"""

import numpy as np

from src.loader import ReviewData


def make_review(rng, review_id, mu, tau, k, n_scale=200, scale="difference"):
    """Create one synthetic ReviewData with a known true pooled mean `mu`.

    mu   : true pooled effect (on the analysis scale; 0 == null).
    tau  : between-study SD (sqrt of tau^2). 0 == homogeneous.
    k    : number of studies.
    n_scale: controls per-study precision (larger -> tighter sei).
    """
    # Per-study sample sizes: heavy-tailed (a few big trials, many small).
    ni = rng.lognormal(mean=np.log(n_scale), sigma=0.8, size=k)
    ni = np.clip(ni, 20, None)

    # Within-study SE shrinks ~ 1/sqrt(n). Constant chosen so a k~10 null review
    # has a realistic chance of being non-significant.
    sei = 1.6 / np.sqrt(ni)

    # Random-effects study means then observed estimates.
    theta_i = rng.normal(mu, tau, size=k)
    yi = rng.normal(theta_i, sei)

    labels = [f"Study{j+1}" for j in range(k)]

    wi = 1.0 / sei ** 2
    theta_fe = float(np.sum(wi * yi) / np.sum(wi))
    se_fe = float(1.0 / np.sqrt(np.sum(wi)))
    ci_lo = theta_fe - 1.96 * se_fe
    ci_hi = theta_fe + 1.96 * se_fe

    return ReviewData(
        review_id=review_id,
        review_doi="",
        analysis_name="synthetic",
        k=k,
        yi=yi.astype(float),
        sei=sei.astype(float),
        ni=ni.astype(float),
        study_labels=labels,
        scale=scale,
        cochrane_pooled=theta_fe,
        cochrane_ci_lo=ci_lo,
        cochrane_ci_hi=ci_hi,
        is_significant=(ci_lo > 0) or (ci_hi < 0),
    )


def make_corpus(seed=12345, n_each=200, k_range=(5, 25),
                effect_mu=0.35, tau_null=0.0, tau_effect=0.0):
    """Generate a corpus of NULL and EFFECT reviews with known labels.

    Returns list of (ReviewData, truth) where truth in {0,1}:
      truth==0 : true mu == 0 (null)
      truth==1 : true mu == effect_mu (real effect)

    The two classes share the same k distribution and precision so any
    separation the classifier achieves is attributable to the signal, not to
    confounded sample size.
    """
    rng = np.random.default_rng(seed)
    corpus = []
    for i in range(n_each):
        k = int(rng.integers(k_range[0], k_range[1] + 1))
        r = make_review(rng, f"NULL{i:04d}", mu=0.0, tau=tau_null, k=k)
        corpus.append((r, 0))
    for i in range(n_each):
        k = int(rng.integers(k_range[0], k_range[1] + 1))
        r = make_review(rng, f"EFF{i:04d}", mu=effect_mu, tau=tau_effect, k=k)
        corpus.append((r, 1))
    return corpus
