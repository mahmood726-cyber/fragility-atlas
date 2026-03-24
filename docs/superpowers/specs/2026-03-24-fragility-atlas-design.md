# Fragility Atlas of Evidence-Based Medicine — Design Spec

## 1. Problem Statement

Published meta-analyses form the backbone of clinical guidelines, yet no study has systematically quantified how sensitive their conclusions are to reasonable analyst choices. Individual fragility studies exist (Walsh 2014, Ioannidis 2007), multiverse analysis exists for primary studies (Steegen 2016), and reproducibility concerns are well-documented — but nobody has combined all dimensions at the scale of hundreds of real Cochrane reviews.

## 2. Research Question

> "What proportion of Cochrane meta-analysis conclusions are robust to reasonable variation in statistical methodology, and which analyst degrees of freedom most frequently alter conclusions?"

## 3. Dataset

**Pairwise70**: 501 Cochrane systematic reviews (RDA files), each containing study-level data:
- Effect sizes (Mean) with 95% CIs (CI.start, CI.end)
- Raw data: Experimental.cases/N, Control.cases/N (binary); Experimental.mean/SD, Control.mean/SD (continuous)
- O-E/Variance (survival)
- ~30,000+ study-level estimates total

**Inclusion criteria**: Reviews with k ≥ 3 studies in the primary analysis (largest k, binary preferred, alphabetical tiebreak). Estimated ~300-400 eligible reviews.

**Data location**: `C:\Models\Pairwise70\data\` (501 .rda files)

## 4. Specification Dimensions

Each eligible review is analyzed across these analyst degrees of freedom:

### 4.1 Estimator (7 levels)
- Fixed-Effect (FE: inverse-variance)
- DerSimonian-Laird (DL)
- Restricted Maximum Likelihood (REML)
- Paule-Mandel (PM)
- Sidik-Jonkman (SJ)
- Hunter-Schmidt (HS)
- Hedges-Olkin (HE)

### 4.2 CI Method (3 levels)
- Wald (normal approximation)
- Hartung-Knapp-Sidik-Jonkman (HKSJ)
- t-distribution (Knapp-Hartung variant)

### 4.3 Outlier Sensitivity (k+1 levels)
- Full set
- Leave-one-out (k studies → k specifications)

### 4.4 Publication Bias Correction (3 levels)
- Unadjusted
- Trim-and-fill (Duval-Tweedie L0)
- PET-PEESE (conditional: PET if p≥0.05, PEESE if p<0.05)

### 4.5 Small-Study Correction (2 levels)
- Standard weights (1/variance)
- HKSJ-adjusted variance

**Total specifications per review**: 7 × 3 × (k+1) × 3 × 2 = 126·(k+1)
For a typical k=10 review: 1,386 specifications.

## 5. Classification System

For each review, examine all specifications and classify:

### 5.1 Per-Specification Decision
- **Significant**: p < 0.05 AND same direction as Cochrane published result
- **Non-significant**: p ≥ 0.05
- **Reversed**: p < 0.05 BUT opposite direction

### 5.2 Review-Level Robustness Score (0-100%)
`robustness = (# specs agreeing with published conclusion) / (total specs) × 100`

### 5.3 Review-Level Classification
| Category | Robustness Score | Meaning |
|----------|-----------------|---------|
| **Robust** | ≥90% | Conclusion holds across nearly all specifications |
| **Moderately Robust** | 70-89% | Conclusion mostly holds |
| **Fragile** | 50-69% | Substantial sensitivity to methods |
| **Unstable** | <50% | Majority of specifications disagree |

### 5.4 Dimension-Level Attribution
For each dimension, compute: what fraction of specification disagreements are attributable to THIS dimension?
- Use ANOVA-like decomposition (eta²) across the factorial design
- Report which dimension (estimator vs CI method vs outlier vs pub-bias) drives fragility

## 6. Architecture

### 6.1 Python Pipeline (`src/`)
```
src/
  loader.py          — Load RDA files, parse, validate, select primary analysis
  estimators.py      — All 7 tau² estimators + FE, with CI methods
  specifications.py  — Generate and execute full specification grid
  corrections.py     — Trim-and-fill, PET-PEESE implementations
  classifier.py      — Robustness scoring and classification
  attribution.py     — Dimension-level variance decomposition (eta²)
  pipeline.py        — Orchestrate: load → specify → compute → classify → export
  utils.py           — Shared math (normalCDF, tCDF, chi2Quantile, etc.)
```

### 6.2 Interactive Dashboard (`dashboard/`)
Single-file HTML (signature style):
- Overview: bar chart of Robust/Moderate/Fragile/Unstable counts
- Specification curve: horizontal dots sorted by effect size, colored by significance
- Dimension attribution: stacked bar of which dimensions drive fragility
- Drill-down: click any review → full specification curve + forest plot
- Export: CSV of all results, print-ready report

### 6.3 Outputs (`data/output/`)
- `fragility_atlas_results.csv` — One row per review with robustness score + classification
- `fragility_atlas_specifications.csv` — One row per specification (large: ~500K rows)
- `fragility_atlas_attribution.csv` — Dimension-level eta² per review
- `fragility_atlas_summary.json` — Aggregate statistics for dashboard

## 7. Statistical Methods

### 7.1 Effect Size Computation
- **Binary (OR/RR)**: Compute from 2×2 table (ai, bi, ci, di) with 0.5 continuity correction
- **Continuous (MD/SMD)**: From means, SDs, Ns. Hedges' g correction for SMD.
- **Survival (HR)**: From O-E and Variance columns
- **SE back-calculation**: If raw data unavailable, SE = (CI.end - CI.start) / (2 × 1.96)

### 7.2 Tau² Estimators
All implemented in pure Python (no R dependency for core pipeline):
- **DL**: Q-based moment estimator (DerSimonian & Laird 1986)
- **REML**: Fisher scoring Newton-Raphson (Viechtbauer 2005)
- **PM**: Iterative (Paule & Mandel 1982)
- **SJ**: Sidik-Jonkman moment estimator
- **HS**: Hunter-Schmidt variance components
- **HE**: Hedges estimator (unweighted)
- **FE**: tau² = 0 (fixed-effect)

### 7.3 Confidence Intervals
- **Wald**: θ̂ ± z_{α/2} × SE(θ̂)
- **HKSJ**: θ̂ ± t_{k-1,α/2} × √(q̂ × SE²(θ̂)), where q̂ = Q/(k-1)
- **t-dist**: θ̂ ± t_{k-1,α/2} × SE(θ̂)

### 7.4 Publication Bias Corrections
- **Trim-and-fill**: Duval-Tweedie L0 estimator, iterative
- **PET-PEESE**: Weighted regression of yi on sei (PET) or sei² (PEESE)
  - Conditional: use PET intercept if PET p ≥ 0.05, else PEESE intercept

### 7.5 Dimension Attribution
For each review with N total specifications:
- Compute eta² = SS_dimension / SS_total for each dimension
- SS_total = variance of binary agreement indicator across all specs
- SS_dimension = between-group variance for that dimension's levels

## 8. Primary Outcome Measures (for the manuscript)

1. **Distribution of robustness scores** across all eligible reviews (histogram)
2. **Proportion in each classification** (Robust/Moderate/Fragile/Unstable)
3. **Most influential dimension** — which analyst choice changes conclusions most often
4. **Risk factors for fragility** — regression of robustness on: k, tau², I², effect magnitude
5. **Specification curve exemplars** — 3-5 representative reviews (1 robust, 1 fragile, 1 unstable)

## 9. Validation Strategy

1. **R cross-validation**: For 10 randomly selected reviews, compare Python results against metafor::rma() with all methods
2. **Known-robust check**: Reviews with very large effects (e.g., k>20, I²<20%) should classify as Robust
3. **Known-fragile check**: Reviews with borderline significance and high heterogeneity should classify as Fragile
4. **Sensitivity of classification**: Bootstrap 95% CI on robustness scores

## 10. Non-Goals (Explicit Exclusions)

- NOT re-extracting effects from PDFs (we use Cochrane data directly)
- NOT computing network meta-analysis (pairwise only)
- NOT doing Bayesian analysis (frequentist multiverse only)
- NOT assessing individual study quality/RoB
- NOT creating a living/updating system (static one-time analysis)

## 11. Success Criteria

1. Pipeline processes ≥300 reviews without errors
2. R cross-validation matches within ±0.001 tolerance on 10/10 reviews
3. Dashboard loads in <3 seconds, all interactions responsive
4. Manuscript-ready tables and figures generated automatically
5. Full analysis reproducible from a single `python pipeline.py` command

## 12. Target

- **Journal**: BMJ (Research Article) — primary target
- **Backup**: Systematic Reviews, Research Synthesis Methods
- **Title**: "The Fragility Atlas: How Robust Are Cochrane Meta-Analysis Conclusions to Analyst Choices? A Large-Scale Multiverse Analysis of 465 Systematic Reviews"
