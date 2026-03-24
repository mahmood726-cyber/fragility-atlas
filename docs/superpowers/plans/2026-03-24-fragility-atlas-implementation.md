# Fragility Atlas — Implementation Plan

## Phase 1: Core Statistical Engine (src/)
1. **utils.py** — Math primitives (normalCDF, tQuantile, lnGamma, chi2CDF)
2. **loader.py** — Load RDA files, select primary analysis, compute yi/sei from raw data
3. **estimators.py** — 7 tau² estimators + 3 CI methods
4. **corrections.py** — Trim-and-fill (L0) + PET-PEESE
5. **specifications.py** — Generate full specification grid, execute all combos
6. **classifier.py** — Robustness scoring, classification, dimension attribution (eta²)
7. **pipeline.py** — Orchestrate end-to-end

## Phase 2: Validation
8. **tests/test_estimators.py** — Unit tests for all 7 estimators against known values
9. **R cross-validation** — Compare 10 reviews against metafor

## Phase 3: Pipeline Execution
10. **Run full pipeline** on all 501 reviews → output CSVs + JSON

## Phase 4: Dashboard
11. **dashboard/index.html** — Single-file interactive HTML dashboard

## Phase 5: Manuscript
12. **manuscript.md** — BMJ-format paper with auto-generated tables/figures
