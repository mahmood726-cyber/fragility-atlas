# The Fragility Atlas: How Robust Are Cochrane Meta-Analysis Conclusions to Analyst Choices? A Large-Scale Multiverse Analysis

## Authors
Mahmood Ahmad^1

^1 Royal Free Hospital, London, UK

Corresponding author: mahmood.ahmad2@nhs.net

ORCID: 0009-0003-7781-4478

---

## Abstract

**Objective:** To systematically quantify the robustness of Cochrane meta-analysis conclusions to reasonable variation in statistical methodology using multiverse analysis.

**Design:** Large-scale multiverse analysis of published systematic reviews.

**Data source:** 403 Cochrane systematic reviews from the Pairwise70 dataset, each with ≥3 studies in the primary analysis.

**Main outcome measures:** For each review, we re-analysed the primary outcome across a comprehensive specification grid: 7 between-study variance estimators (fixed-effect, DerSimonian-Laird, REML, Paule-Mandel, Sidik-Jonkman, Hunter-Schmidt, Hedges) × 3 confidence interval methods (Wald, Hartung-Knapp-Sidik-Jonkman, t-distribution) × 3 publication bias corrections (none, trim-and-fill, PET-PEESE) × (k+1) leave-one-out subsets, yielding 394,569 total specifications. Each specification was classified as agreeing or disagreeing with the reference conclusion (DerSimonian-Laird with Wald CI). Reviews were classified as Robust (≥90% agreement), Moderately Robust (70-89%), Fragile (50-69%), or Unstable (<50%).

**Results:** Of 403 eligible Cochrane reviews, 80 (19.9%) were classified as Robust, 89 (22.1%) as Moderately Robust, 145 (36.0%) as Fragile, and 89 (22.1%) as Unstable. The mean robustness score was 66.4% (median 66.7%, IQR 56.5-87.2%). The most influential analytical dimension was publication bias correction (mean eta-squared = 0.374), followed by CI method (eta-squared = 0.036). A total of 394,569 individual meta-analytic specifications were computed across all reviews.

**Conclusions:** In this large-scale multiverse analysis of 403 Cochrane meta-analyses, we found that 58% of reviews had conclusions that were sensitive to reasonable analytical choices. The publication bias correction dimension was the single most influential driver of disagreement, explaining ten times more variance than the choice of variance estimator. These findings suggest that the conclusions of evidence-based medicine are more contingent on methodological choices than is commonly appreciated, and highlight the need for routine multiverse sensitivity analysis in meta-analytic practice.

---

## Introduction

Meta-analyses sit at the apex of the evidence hierarchy and directly inform clinical practice guidelines, drug approvals, and health policy decisions [1]. Yet every meta-analysis involves a series of methodological choices — which between-study variance estimator to use, how to construct confidence intervals, whether to adjust for publication bias, and how to handle potentially influential studies. These "analyst degrees of freedom" can substantially alter conclusions, but the extent of this sensitivity across the medical literature has never been systematically quantified.

The concept of multiverse analysis, introduced by Steegen et al. [2], addresses the problem of researcher degrees of freedom by executing all reasonable combinations of analytical choices and examining the distribution of results. While multiverse analysis has been applied to individual primary studies [3] and small collections of meta-analyses [4], no study has applied it at scale to a large corpus of published systematic reviews.

The fragility of meta-analytic conclusions has been examined through several narrower lenses. The Fragility Index [5] quantifies how many events would need to change to alter statistical significance, but only for binary outcomes and only along one dimension (study inclusion). Studies comparing variance estimators [6-8] have shown that different methods can yield substantially different results, particularly with few studies or high heterogeneity. The choice between Wald and Hartung-Knapp-Sidik-Jonkman confidence intervals has been shown to matter considerably [9]. Publication bias corrections can dramatically alter conclusions [10]. However, these dimensions have always been examined in isolation.

We present the Fragility Atlas — a comprehensive multiverse analysis of 403 Cochrane systematic reviews. For each review, we exhaustively varied five analytical dimensions simultaneously, generating 394,569 total specifications. This allows us to answer, for the first time at scale: How robust are the conclusions of evidence-based medicine to reasonable methodological variation?

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

For a review with the median k = 8 studies, this yields 567 specifications.

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

All analyses were conducted in Python 3.13 using scipy 1.16.2 for statistical distributions. The complete pipeline is available as open-source software at https://github.com/mahmood726-cyber/fragility-atlas. An interactive dashboard for exploring results is available at https://github.com/mahmood726-cyber/fragility-atlas#dashboard.

### Patient and public involvement

No patients or members of the public were involved in the design or conduct of this study.

## Results

### Overview

403 Cochrane systematic reviews met inclusion criteria (k ≥ 3 studies). The median number of studies per review was 8 (IQR 5-16; range 3-180). 363 reviews (90.1%) used ratio-scale outcomes and 40 (9.9%) used difference-scale outcomes. A total of 394,569 individual meta-analytic specifications were computed.

### Distribution of robustness

The mean robustness score was 66.4% (median 66.7%, IQR 56.5%-87.2%). Figure 1 shows the distribution of robustness scores.

80 reviews (19.9%) were classified as Robust, 89 (22.1%) as Moderately Robust, 145 (36.0%) as Fragile, and 89 (22.1%) as Unstable (Table 1). Notably, the distribution of robustness scores was strongly bimodal, with peaks near 66.7% (the most common score, reflecting reviews where exactly one of three bias correction approaches disagreed) and near 100% (reviews with strong, unambiguous effects).

### Dimension attribution

Across all reviews, the most influential analytical dimension was publication bias correction (mean η² = 0.374), followed by leave-one-out sensitivity (η² = 0.156), CI method (η² = 0.036), and variance estimator (η² = 0.035) (Table 2). The dominance of bias correction was consistent across review sizes: for reviews with k = 3-5, bias correction η² = 0.38; for k > 50, η² = 0.35.

### Predictors of fragility

We examined whether review-level characteristics predicted robustness. In a multiple linear regression including log(k), absolute effect magnitude, effect size IQR, and outcome scale (ratio vs difference), the overall model explained only 1.8% of variance in robustness scores (R² = 0.018, adjusted R² = 0.006). Only absolute effect magnitude was a marginally significant predictor (beta = 0.64, p = 0.04); larger effects were associated with slightly higher robustness. The number of studies (log k: beta = 2.07, p = 0.12), effect size heterogeneity (IQR: beta = -0.87, p = 0.26), and outcome scale (ratio vs difference: beta = 3.14, p = 0.45) were not significant predictors.

Robustness did not meaningfully vary by number of studies: reviews with k = 3-5 had mean robustness of 64.3% compared to 72.1% for reviews with k > 50 (Table 3). The proportion classified as Fragile or Unstable was 54-69% across all k-bands, indicating that methodological sensitivity is pervasive regardless of the amount of available evidence.

## Discussion

### Summary of findings

In this large-scale multiverse analysis of 403 Cochrane systematic reviews encompassing 394,569 individual meta-analytic specifications, we found that the majority of meta-analysis conclusions are sensitive to reasonable analytical choices. Only one in five reviews (19.9%) produced conclusions that were robust across more than 90% of specifications. More than half (58.1%) were classified as Fragile or Unstable, meaning their conclusions changed under a substantial proportion of plausible analytical approaches.

The most striking finding was the dominant role of publication bias correction. Whether an analyst chose to apply trim-and-fill or PET-PEESE correction — versus reporting the unadjusted result — explained more than twice the variance of any other single dimension (mean η² = 0.374). Leave-one-out sensitivity was the second-largest driver (η² = 0.156), while the choice of variance estimator (η² = 0.035) and confidence interval method (η² = 0.036) contributed comparatively little. This suggests that the single most consequential methodological decision in meta-analysis is not how to pool the data, but whether and how to account for the selective publication of studies.

Equally notable was the finding that fragility was essentially unpredictable from observable review characteristics. The number of included studies, the magnitude of the pooled effect, the degree of heterogeneity, and the type of outcome scale collectively explained less than 2% of variance in robustness scores (R² = 0.018). A review with 50 studies was nearly as likely to have a fragile conclusion as one with 5 studies. This means that neither reviewers nor guideline panels can rely on surface-level features of a meta-analysis to judge the stability of its conclusions.

### Comparison with existing literature

Our finding that 58% of Cochrane meta-analysis conclusions are sensitive to analyst choices extends prior work examining individual dimensions of analytical sensitivity. Walsh et al. [5] reported that 53% of meta-analyses had a Fragility Index of 3 or fewer, indicating vulnerability to small changes in event counts. However, the Fragility Index operates along a single dimension (study inclusion for binary outcomes) and cannot capture methodological sensitivity. Ioannidis and Trikalinos [20] demonstrated that many meta-analyses have fragile conclusions when single studies are removed, but did not examine the interaction between study removal and methodological variation.

Several simulation studies have compared variance estimators and shown that different methods can yield substantially different results, particularly with fewer than 10 studies or high heterogeneity [6-8]. IntHout et al. [9] showed that HKSJ confidence intervals are often considerably wider than Wald intervals and recommended HKSJ as the default. Our empirical findings confirm this at scale: the choice between Wald and HKSJ/t-distribution does alter some conclusions (η² = 0.036), but this effect is dwarfed by the publication bias correction decision.

The dominant influence of publication bias correction aligns with the theoretical work of Carter et al. [10], who showed in simulation that different bias correction methods can produce widely divergent results. Our study demonstrates that this is not merely a theoretical concern — in real Cochrane reviews, the decision to apply (or not apply) publication bias correction is the single most consequential analytical choice.

The multiverse analysis framework itself has been applied to meta-analysis by Voracek et al. [4], but only to a small number of reviews and without the leave-one-out and bias correction dimensions. Our study is, to our knowledge, the first to apply a comprehensive five-dimensional multiverse analysis across hundreds of published meta-analyses.

### Strengths and limitations

**Strengths:**
1. First study to apply comprehensive multiverse analysis across hundreds of real meta-analyses simultaneously varying five analytical dimensions
2. Based on Cochrane reviews — the gold standard for systematic review methodology, providing a conservative estimate of fragility (less rigorous reviews may be more sensitive)
3. Pre-specified analysis plan with a fully deterministic, reproducible pipeline validated against R's metafor package (9/9 matched reviews showed agreement to four decimal places)
4. Open data and code for full transparency and reproducibility
5. The finding that fragility is unpredictable (R² = 0.018) emerged from the data rather than being assumed, and is itself methodologically important

**Limitations:**
1. Standard errors were back-calculated from published confidence intervals using the normal approximation (dividing the CI width by 2 × 1.96), which assumes Cochrane's standard 95% Wald CIs. While this is nearly always correct for Cochrane data, it introduces a small approximation error for reviews that used non-standard intervals
2. The specification dimensions, while comprehensive, do not capture all possible analyst choices. Continuity corrections for zero cells, choice of effect measure (e.g., OR vs RR), Bayesian approaches, and robust variance estimation were not included
3. Leave-one-out sensitivity removes individual studies but does not capture the effect of adding unpublished studies or studies excluded by the original reviewers
4. Both trim-and-fill and PET-PEESE have known limitations: trim-and-fill can over-correct when heterogeneity is high [18], and PET-PEESE assumes a specific functional form for the relationship between effect size and standard error [19]. These limitations may contribute to the large η² for the bias correction dimension
5. The classification thresholds (90/70/50%) are necessarily arbitrary. However, the key finding — that a majority of reviews are sensitive to analytical choices — is robust to alternative thresholds: using 95/80/60% cutoffs, 67% of reviews would be classified as Fragile or Unstable; using 85/65/45% cutoffs, 48% would be so classified
6. Our analysis used the largest analysis within each review (preferring binary outcomes), which may not always correspond to the primary outcome as designated by the original authors

### Implications for practice and policy

These findings have several implications for the conduct, reporting, and interpretation of meta-analyses.

**For systematic reviewers:** Our results support recent calls to move beyond single-method meta-analysis [4]. At minimum, reviewers should report sensitivity analyses across variance estimators (DL vs REML at a minimum), confidence interval methods (Wald vs HKSJ), and with and without publication bias correction. A multiverse-style specification curve could be provided as supplementary material to transparently show how conclusions vary across analytical choices.

**For clinical guideline panels:** The finding that 58% of meta-analysis conclusions are method-dependent has direct relevance for GRADE assessments of certainty of evidence. Currently, GRADE downgrades for inconsistency (heterogeneity across studies) but does not systematically assess inconsistency across analytical methods. Our findings suggest that "analytical inconsistency" deserves consideration as a dimension of evidence quality.

**For journal editors and peer reviewers:** The dominant role of publication bias correction (η² = 0.374) highlights a paradox in current practice: most Cochrane reviews report unadjusted pooled estimates as their primary result, with bias-corrected estimates relegated to sensitivity analyses (if reported at all). Yet this single decision drives more disagreement than all other methodological choices combined. Journals should consider requiring that bias-corrected estimates be reported alongside unadjusted results.

**For methodologists:** The near-zero predictive power of review characteristics (R² = 0.018) means that there is no shortcut for assessing robustness. Neither large sample sizes nor large effect sizes confer methodological stability. This argues for routine multiverse analysis as a standard component of meta-analytic practice, rather than selective application to "uncertain" reviews.

## Conclusions

In this comprehensive multiverse analysis of 403 Cochrane systematic reviews, we found that 58% of meta-analysis conclusions were sensitive to reasonable variation in statistical methodology. Only 20% of reviews produced conclusions that were robust across more than 90% of analytical specifications. The decision of whether to correct for publication bias was the single most influential analytical choice, explaining an order of magnitude more variation than the choice of variance estimator or confidence interval method. Fragility was pervasive and unpredictable from review characteristics, affecting reviews across all sizes, effect magnitudes, and clinical domains. These findings underscore the need for routine multiverse sensitivity analysis in meta-analytic practice and transparent reporting of how conclusions depend on methodological choices.

## Data availability statement

The Pairwise70 dataset is available from [TO BE DEPOSITED]. The analysis code, interactive dashboard, and all output files are available at https://github.com/mahmood726-cyber/fragility-atlas. A TruthCert provenance bundle documenting the complete analysis chain is included in the repository.

## Funding

None.

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

| Classification | Criteria | n | % |
|---|---|---|---|
| Robust | ≥90% specifications agree | 80 | 19.9 |
| Moderately Robust | 70-89% agree | 89 | 22.1 |
| Fragile | 50-69% agree | 145 | 36.0 |
| Unstable | <50% agree | 89 | 22.1 |
| **Total** | | **403** | **100.0** |

Mean robustness score: 66.4% (median 66.7%, IQR 56.5-87.2%, range 5.8-100.0%).

### Table 2. Dimension Attribution (Mean η² Across All Reviews)

| Analytical Dimension | Levels | Mean η² | Interpretation |
|---|---|---|---|
| Publication bias correction | None, Trim-and-fill, PET-PEESE | 0.374 | Dominant driver of disagreement |
| CI method | Wald, HKSJ, t-distribution | 0.036 | Minor influence |
| Leave-one-out | Full set, each study removed | 0.156 | Second-largest driver |
| Variance estimator | FE, DL, REML, PM, SJ, HS, HE | 0.035 | Minor influence |

η² represents the proportion of variance in the agreement indicator (agree/disagree with reference) explained by each dimension, averaged across all reviews.

### Table 3. Robustness by Number of Studies (k-bands)

| k-band | n reviews | Mean robustness (%) | Median (%) | % Fragile or Unstable |
|---|---|---|---|---|
| 3-5 | 143 | 64.3 | 66.7 | 54 |
| 6-10 | 110 | 65.2 | 66.7 | 55 |
| 11-20 | 67 | 67.4 | 66.7 | 69 |
| 21-50 | 64 | 70.6 | 66.7 | 61 |
| 51-180 | 19 | 72.1 | 66.7 | 58 |

### Table 4. Multiple Regression: Predictors of Robustness Score

| Variable | Beta | SE | t | p |
|---|---|---|---|---|
| Intercept | 58.99 | 5.19 | 11.37 | <0.001 |
| log(k) | 2.07 | 1.32 | 1.57 | 0.118 |
| |Effect magnitude| | 0.64 | 0.31 | 2.06 | 0.040 |
| IQR(effect) | -0.87 | 0.76 | -1.14 | 0.255 |
| Ratio scale (vs difference) | 3.14 | 4.15 | 0.76 | 0.449 |

Overall model: R² = 0.018, adjusted R² = 0.006, n = 403.

## Figures

### Figure 1. Distribution of Robustness Scores Across 403 Cochrane Meta-Analyses
Histogram of robustness scores (0-100%) with 5-percentage-point bins, colour-coded by classification (green = Robust, yellow = Moderate, red = Fragile, purple = Unstable). Available in the interactive dashboard.

### Figure 2. Specification Curve for an Exemplar Fragile Review
All specifications sorted by pooled effect estimate (x-axis), with dots coloured by statistical significance. The reference specification (DL + Wald + unadjusted) is highlighted. Panels below show which dimension levels apply to each specification.

### Figure 3. Dimension Attribution Across All Reviews
Horizontal bar chart showing mean η² for each analytical dimension. Publication bias correction (η² = 0.374) dominates, with CI method, variance estimator, and leave-one-out each contributing η² < 0.04.

### Figure 4. Relationship Between Number of Studies and Robustness
Scatter plot of k (log scale, x-axis) versus robustness score (y-axis), with points coloured by classification. The near-flat relationship (r = 0.081, p = 0.107) demonstrates that fragility is not resolved by including more studies.
