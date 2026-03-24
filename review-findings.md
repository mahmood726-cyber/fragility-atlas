## REVIEW CLEAN — 2 ROUNDS, ALL ISSUES FIXED
## Multi-Persona Review: Fragility Atlas (all src/, tests/, dashboard/)
### Date: 2026-03-24
### Round 1: 9 P0 + 8 P1 + 6 P2 = 23 issues → ALL FIXED
### Round 2: 1 P0 + 2 P1 + 2 P2 = 5 issues → ALL FIXED
### Statistical Methodologist Round 2: CLEAN (all formulas verified)
### Tests: 31/31 pass, pipeline verified: 403 reviews, 394,569 specs

---

## Persona 1: Statistical Methodologist

#### P0 -- Critical

- **SM-P0-1** Cochrane pooled uses median of study-level, not actual pooled (loader.py:72-74)
  - cochrane_pooled = np.nanmedian(means) is NOT the inverse-variance pooled estimate
  - is_significant determined from these medians → wrong reference conclusion
  - Fix: compute FE pooled from yi/sei

- **SM-P0-2** Scale detection infers ratio/difference from sign, not outcome type (loader.py:53)
  - All-positive MD misclassified as ratio → log-transform corrupts effects
  - Fix: use Pairwise70 outcome metadata or raw-data column presence

- **SM-P0-3** Trim-and-fill only detects RIGHT-side asymmetry (corrections.py:36-48)
  - Protective effects (log-RR < 0) have LEFT-side bias → k0=0, no correction
  - Fix: auto-detect side per metafor's approach

- **SM-P0-4** PET-PEESE uses z-distribution, not t (corrections.py:148-149, 94-97)
  - WLS regression has k-2 df; z is anti-conservative (esp k<10)
  - Fix: use t_cdf(z, k-2) and t_quantile(p, k-2)

- **SM-P0-5** _weighted_regression denom threshold 1e-300 too strict (corrections.py:129)
  - Near-equal SEs produce denom ~ 1e-13, passes check, causes division by near-zero
  - Fix: use relative threshold

#### P1 -- Important

- **SM-P1-1** LOO eta² uses binary grouping, dilutes influential-study signal (classifier.py:134)
- **SM-P1-2** Silent exception swallowing in spec loop (specifications.py:82)
- **SM-P1-3** Test coverage ~30% — no tests for corrections, classifier, loader, or specifications
- **SM-P1-4** PET-PEESE CI ignores ci_method parameter (corrections.py:94-97)
- **SM-P1-5** embed_dashboard_data.py uses fragile brace-counting (embed_dashboard_data.py:42-79)

#### P2 -- Minor

- **SM-P2-1** Median/quartile uses integer indexing, not numpy (pipeline.py:106)
- **SM-P2-2** SJ is one-step, not iterative (differs from metafor) (estimators.py:181-198)
- **SM-P2-3** REML non-convergence returns silently (estimators.py:146)
- **SM-P2-4** Hardcoded z=1.96 for SE back-calc (loader.py:151)
- **SM-P2-5** conf_level not propagated to loader (loader.py/pipeline.py)

---

## Persona 2: Security Auditor + Software Engineer

#### P0 -- Critical

- **SE-P0-1** CRASH: undefined rng() in dashboard scatter plot (dashboard:1017-1019)
  - Fisher-Yates shuffle calls rng() but PRNG was removed when real data was injected
  - Fix: add seeded PRNG back, or remove shuffle code

- **SE-P0-2** CRASH: populate_manuscript.py divides by n=0 (populate_manuscript.py:26-58)
  - Fix: add early return for n=0

- **SE-P0-3** CRASH: meta_analysis() with k=1 → division by zero in SJ/HS/HE (estimators.py:187,206,216)
  - Fix: add k<2 guard in those functions

- **SE-P0-4** = SM-P1-2: bare except swallows important errors (specifications.py:82)

#### P1 -- Important

- **SE-P1-1** = SM-P2-4: hardcoded z=1.96 (loader.py:151)
- **SE-P1-2** trim_and_fill loop performance (corrections.py:23-62)
- **SE-P1-3** populate_manuscript.py unguarded dict indexing (populate_manuscript.py:26-34)
- **SE-P1-4** XSS (low risk): showTooltip innerHTML (dashboard:594)
- **SE-P1-5** No validation of conf_level range (pipeline.py:211-212)
- **SE-P1-6** top_dimension field missing from Python pipeline (classifier.py, pipeline.py)
- **SE-P1-7** = SM-P1-5: brittle brace parser (embed_dashboard_data.py)
- **SE-P1-8** Median uses n//2 not true median (pipeline.py:106)

#### P2 -- Minor

- **SE-P2-1** CSV injection not guarded (pipeline.py:122-155)
- **SE-P2-2** iterrows() slow (loader.py:132)
- **SE-P2-3** Quantile method differs from numpy (pipeline.py:176-178)
- **SE-P2-4** Duplicate significance logic (loader.py:78-82)
- **SE-P2-5** Div balance OK: 19/19
- **SE-P2-6** Histogram bars missing keyboard tooltip (dashboard:739)

---

## Deduplicated P0 List (for fixing)

| ID | Description | File | Status |
|----|-------------|------|--------|
| P0-1 | Cochrane pooled uses median not actual FE pooled | loader.py:72-74 | [FIXED] |
| P0-2 | Scale detection by sign not outcome type | loader.py:53 | [FIXED] |
| P0-3 | Trim-and-fill right-side only | corrections.py:36-48 | [FIXED] |
| P0-4 | PET-PEESE uses z not t distribution | corrections.py | [FIXED] |
| P0-5 | WLS denom threshold too strict | corrections.py:129 | [FIXED] |
| P0-6 | Dashboard rng() undefined after data embed | dashboard:1017 | [FIXED] |
| P0-7 | populate_manuscript.py n=0 crash | populate_manuscript.py | [FIXED] |
| P0-8 | k=1 division by zero in SJ/HS/HE | estimators.py | [FIXED] |
| P0-9 | Bare except swallows errors | specifications.py:82 | [FIXED] |
