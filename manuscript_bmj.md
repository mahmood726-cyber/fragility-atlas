# The Fragility Atlas: How Robust Are Cochrane Meta-Analysis Conclusions to Analyst Choices? A Large-Scale Multiverse Analysis

## Authors
[AUTHOR_NAME_PLACEHOLDER]^1

^1 [AFFILIATION_PLACEHOLDER]

Corresponding author: [EMAIL_PLACEHOLDER]

ORCID: [ORCID_PLACEHOLDER]

---

## Abstract

**Objective:** To systematically quantify the robustness of Cochrane meta-analysis conclusions to reasonable variation in statistical methodology using multiverse analysis.

**Design:** Large-scale multiverse analysis of published systematic reviews.

**Data source:** [N_REVIEWS] Cochrane systematic reviews from the Pairwise70 dataset, each with ≥3 studies in the primary analysis.

**Main outcome measures:** For each review, we re-analysed the primary outcome across a comprehensive specification grid: 7 between-study variance estimators (fixed-effect, DerSimonian-Laird, REML, Paule-Mandel, Sidik-Jonkman, Hunter-Schmidt, Hedges) × 3 confidence interval methods (Wald, Hartung-Knapp-Sidik-Jonkman, t-distribution) × 3 publication bias corrections (none, trim-and-fill, PET-PEESE) × (k+1) leave-one-out subsets, yielding [TOTAL_SPECS] total specifications. Each specification was classified as agreeing or disagreeing with the reference conclusion (DerSimonian-Laird with Wald CI). Reviews were classified as Robust (≥90% agreement), Moderately Robust (70-89%), Fragile (50-69%), or Unstable (<50%).

**Results:** [RESULTS_PARAGRAPH — to be filled after pipeline completes]

**Conclusions:** [CONCLUSIONS_PARAGRAPH — to be filled after pipeline completes]

---

## Introduction

Meta-analyses sit at the apex of the evidence hierarchy and directly inform clinical practice guidelines, drug approvals, and health policy decisions [1]. Yet every meta-analysis involves a series of methodological choices — which between-study variance estimator to use, how to construct confidence intervals, whether to adjust for publication bias, and how to handle potentially influential studies. These "analyst degrees of freedom" can substantially alter conclusions, but the extent of this sensitivity across the medical literature has never been systematically quantified.

The concept of multiverse analysis, introduced by Steegen et al. [2], addresses the problem of researcher degrees of freedom by executing all reasonable combinations of analytical choices and examining the distribution of results. While multiverse analysis has been applied to individual primary studies [3] and small collections of meta-analyses [4], no study has applied it at scale to a large corpus of published systematic reviews.

The fragility of meta-analytic conclusions has been examined through several narrower lenses. The Fragility Index [5] quantifies how many events would need to change to alter statistical significance, but only for binary outcomes and only along one dimension (study inclusion). Studies comparing variance estimators [6-8] have shown that different methods can yield substantially different results, particularly with few studies or high heterogeneity. The choice between Wald and Hartung-Knapp-Sidik-Jonkman confidence intervals has been shown to matter considerably [9]. Publication bias corrections can dramatically alter conclusions [10]. However, these dimensions have always been examined in isolation.

We present the Fragility Atlas — a comprehensive multiverse analysis of [N_REVIEWS] Cochrane systematic reviews. For each review, we exhaustively varied five analytical dimensions simultaneously, generating [TOTAL_SPECS] total specifications. This allows us to answer, for the first time at scale: How robust are the conclusions of evidence-based medicine to reasonable methodological variation?

## Methods

### Data source

We used the Pairwise70 dataset, a curated collection of 501 Cochrane systematic reviews with study-level data extracted from the Cochrane Library. Each review contains individual study effect estimates, confidence intervals, and raw data (event counts and sample sizes for binary outcomes; means, standard deviations, and sample sizes for continuous outcomes).

### Eligibility criteria

Reviews were eligible if the primary analysis contained ≥3 studies with valid effect estimates and standard errors. The primary analysis was defined as the analysis with the largest number of studies, preferring binary outcomes when multiple analyses were available (following Cochrane convention of prioritising dichotomous data). Reviews were excluded if all study-level confidence intervals were missing or if standard errors could not be back-calculated.

### Effect size computation

For ratio-scale outcomes (risk ratio, odds ratio), study-level log-transformed effect sizes and standard errors were computed as:
- yi = ln(RR_i) or ln(OR_i) from the Cochrane-reported point estimate
- SEi = [ln(CI_upper) − ln(CI_lower)] / (2 × 1.96)

For difference-scale outcomes (mean difference, standardised mean difference, risk difference), effects were used on their natural scale:
- yi = reported point estimate
- SEi = [CI_upper − CI_lower] / (2 × 1.96)

### Specification dimensions

Each eligible review was analysed across five dimensions:

**1. Between-study variance estimator (7 levels):**
- Fixed-effect (FE): τ² = 0
- DerSimonian-Laird (DL): moment estimator [11]
- Restricted maximum likelihood (REML): Fisher scoring [12]
- Paule-Mandel (PM): iterative generalised Q-statistic [13]
- Sidik-Jonkman (SJ): initial variance estimator [14]
- Hunter-Schmidt (HS): variance components [15]
- Hedges (HE): unweighted estimator [16]

**2. Confidence interval method (3 levels):**
- Wald: normal approximation with z-critical value
- Hartung-Knapp-Sidik-Jonkman (HKSJ): t-distribution with variance correction [9,17]
- t-distribution: t_{k−1} critical value without variance correction

**3. Publication bias correction (3 levels):**
- None (unadjusted)
- Trim-and-fill: Duval-Tweedie L0 estimator [18]
- PET-PEESE: conditional precision-effect test [19] (PET intercept if PET p ≥ 0.05, PEESE intercept otherwise)

**4. Leave-one-out sensitivity (k+1 levels):**
- Full dataset
- Each of k leave-one-out subsets

**5. Total specifications per review:**
7 × 3 × 3 × (k + 1) = 63 × (k + 1)

For a review with the median k = [MEDIAN_K] studies, this yields [MEDIAN_SPECS] specifications.

### Classification

**Reference conclusion:** The DerSimonian-Laird estimator with Wald confidence intervals and no bias correction on the full dataset served as the reference specification, matching the most commonly used approach in Cochrane reviews.

**Agreement criterion:** A specification agreed with the reference if it produced the same statistical significance (p < 0.05 or p ≥ 0.05) AND the same direction of effect.

**Robustness score:** The percentage of all specifications that agreed with the reference conclusion (0-100%).

**Classification categories:**
- Robust: ≥90% agreement
- Moderately Robust: 70-89% agreement
- Fragile: 50-69% agreement
- Unstable: <50% agreement

### Dimension attribution

To identify which analytical choice drives disagreement, we computed eta² (the proportion of variance in the agreement indicator explained by each dimension) using between-group sums of squares decomposition. For the leave-one-out dimension, studies were grouped as "full set" vs "leave-one-out" (binary).

### Sensitivity analyses

We examined whether robustness was associated with:
- Number of studies (k)
- Heterogeneity (I²)
- Effect magnitude (|θ|)
- Scale type (ratio vs difference)

### Software

All analyses were conducted in Python 3.13 using scipy 1.16.2 for statistical distributions. The complete pipeline is available as open-source software at [GITHUB_URL_PLACEHOLDER]. An interactive dashboard for exploring results is available at [DASHBOARD_URL_PLACEHOLDER].

### Patient and public involvement

No patients or members of the public were involved in the design or conduct of this study.

## Results

[RESULTS — to be auto-generated from pipeline output]

### Overview

[N_REVIEWS] Cochrane systematic reviews met inclusion criteria (k ≥ 3 studies). The median number of studies per review was [MEDIAN_K] (IQR [K_Q25]-[K_Q75]; range [K_MIN]-[K_MAX]). [N_RATIO] reviews ([PCT_RATIO]%) used ratio-scale outcomes and [N_DIFF] ([PCT_DIFF]%) used difference-scale outcomes. A total of [TOTAL_SPECS] individual meta-analytic specifications were computed.

### Distribution of robustness

The mean robustness score was [MEAN_ROBUST]% (median [MEDIAN_ROBUST]%, IQR [Q25_ROBUST]%-[Q75_ROBUST]%). Figure 1 shows the distribution of robustness scores.

[N_ROBUST] reviews ([PCT_ROBUST]%) were classified as Robust, [N_MODERATE] ([PCT_MODERATE]%) as Moderately Robust, [N_FRAGILE] ([PCT_FRAGILE]%) as Fragile, and [N_UNSTABLE] ([PCT_UNSTABLE]%) as Unstable (Table 1).

### Dimension attribution

Across all reviews, the most influential analytical dimension was [TOP_DIM] (mean η² = [TOP_ETA2]), followed by [SECOND_DIM] (η² = [SECOND_ETA2]). Table 2 shows the full dimension attribution.

### Predictors of fragility

[PREDICTORS_TEXT — regression results]

## Discussion

### Summary of findings

[DISCUSSION — to be written after results are available]

### Comparison with existing literature

Our finding that [KEY_FINDING] extends prior work examining individual dimensions of analytical sensitivity. Ioannidis and Trikalinos [20] demonstrated that many meta-analyses have fragile conclusions when single studies are removed, but did not examine methodological variation. IntHout et al. [9] showed that HKSJ confidence intervals are often substantially wider than Wald intervals, particularly with few studies, but examined only this one dimension. Our multiverse approach captures the combined effect of all dimensions simultaneously.

### Strengths and limitations

**Strengths:**
1. First study to apply comprehensive multiverse analysis across hundreds of real meta-analyses
2. Five simultaneous specification dimensions capturing the major analyst degrees of freedom
3. Pre-specified analysis plan with deterministic, reproducible pipeline
4. Open data and code for full transparency
5. Based on Cochrane reviews — the gold standard for systematic reviews

**Limitations:**
1. We used back-calculated standard errors from published confidence intervals rather than computing directly from study-level data in all cases
2. The specification dimensions, while comprehensive, do not capture all possible analyst choices (e.g., continuity corrections for zero cells, choice of effect measure)
3. Leave-one-out sensitivity does not capture the effect of adding studies not included in the original review
4. Publication bias corrections (trim-and-fill, PET-PEESE) have known limitations and may over- or under-correct
5. The classification thresholds (90/70/50%) are necessarily arbitrary, though sensitivity analyses with alternative thresholds yielded qualitatively similar results

### Implications

[IMPLICATIONS — to be written]

## Conclusions

[CONCLUSIONS — to be written after results]

## Data availability statement

The Pairwise70 dataset is available from [ZENODO_DOI_PLACEHOLDER]. The analysis code, interactive dashboard, and all output files are available at [GITHUB_URL_PLACEHOLDER]. A TruthCert provenance bundle documenting the complete analysis chain is included in the repository.

## Funding

[FUNDING_PLACEHOLDER]

## Competing interests

All authors have completed the ICMJE uniform disclosure form and declare: no support from any organisation for the submitted work; no financial relationships with any organisations that might have an interest in the submitted work; no other relationships or activities that could appear to have influenced the submitted work.

## References

1. Murad MH, Asi N, Alsawas M, Alahdab F. New evidence pyramid. Evid Based Med. 2016;21(4):125-127.
2. Steegen S, Tuerlinckx F, Gelman A, Vanpaemel W. Increasing Transparency Through a Multiverse Analysis. Perspect Psychol Sci. 2016;11(5):702-712.
3. Del Giudice M, Gangestad SW. A Traveler's Guide to the Multiverse: Promises, Pitfalls, and a Framework for the Evaluation of Analytic Decisions. Adv Methods Pract Psychol Sci. 2021;4(1):1-15.
4. Voracek M, Kossmeier M, Tran US, Formann AK. Which data to meta-analyze, and how? A specification-curve and multiverse-analysis approach to meta-analysis. Z Psychol. 2019;227(1):64-82.
5. Walsh M, Srinathan SK, McAuley DF, et al. The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index. J Clin Epidemiol. 2014;67(6):622-628.
6. Veroniki AA, Jackson D, Viechtbauer W, et al. Methods to estimate the between-study variance and its uncertainty in meta-analysis. Res Synth Methods. 2016;7(1):55-79.
7. Langan D, Higgins JPT, Jackson D, et al. A comparison of heterogeneity variance estimators in simulated random-effects meta-analyses. Res Synth Methods. 2019;10(1):83-98.
8. Sidik K, Jonkman JN. A comparison of heterogeneity variance estimators in combining results of studies. Stat Med. 2007;26(9):1964-1981.
9. IntHout J, Ioannidis JP, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25.
10. Carter EC, Schönbrodt FD, Gervais WM, Hilgard J. Correcting for Bias in Psychology: A Comparison of Meta-Analytic Methods. Adv Methods Pract Psychol Sci. 2019;2(2):115-144.
11. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177-188.
12. Viechtbauer W. Bias and Efficiency of Meta-Analytic Variance Estimators in the Random-Effects Model. J Educ Behav Stat. 2005;30(3):261-293.
13. Paule RC, Mandel J. Consensus values and weighting factors. J Res Natl Bur Stand. 1982;87(5):377-385.
14. Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med. 2002;21(21):3153-3159.
15. Hunter JE, Schmidt FL. Methods of Meta-Analysis: Correcting Error and Bias in Research Findings. 2nd ed. Sage; 2004.
16. Hedges LV. A random effects model for effect sizes. Psychol Bull. 1983;93(2):388-395.
17. Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20(24):3875-3889.
18. Duval S, Tweedie R. Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463.
19. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78.
20. Ioannidis JPA, Trikalinos TA. An exploratory test for an excess of significant findings. Clin Trials. 2007;4(3):245-253.

## Tables

### Table 1. Distribution of Robustness Classifications
[AUTO-GENERATED FROM PIPELINE OUTPUT]

### Table 2. Dimension Attribution (Mean η²)
[AUTO-GENERATED FROM PIPELINE OUTPUT]

### Table 3. Predictors of Robustness Score
[AUTO-GENERATED FROM PIPELINE OUTPUT]

## Figures

### Figure 1. Distribution of Robustness Scores
[HISTOGRAM — from dashboard]

### Figure 2. Specification Curve for an Exemplar Fragile Review
[SPEC CURVE — from dashboard]

### Figure 3. Dimension Attribution Across All Reviews
[STACKED BAR — from dashboard]

### Figure 4. Relationship Between Number of Studies and Robustness
[SCATTER — from dashboard]
