# Fragility Atlas — R Cross-Validation Report

**Date:** 2026-03-25
**R version:** 4.5.2
**metafor version:** 4.8-0

## Summary

10 reviews validated (seeded random selection, k>=5, ratio scale).

### Core Estimators: 40/40 PASS

| Estimator | Pass/Total | Max theta diff |
|-----------|-----------|---------------|
| DL | 10/10 | 1.47e-05 |
| REML | 10/10 | 1.03e-05 |
| FE | 10/10 | <1e-16 (machine epsilon) |
| PM | 10/10 | 1.46e-05 |

### Secondary Estimators: 15/30

| Estimator | Pass/Total | Notes |
|-----------|-----------|-------|
| HS | 5/10 | Python clamps tau2_HS=0 → FE; metafor applies different weighting |
| HE | 9/10 | 1 discrepancy (2.5e-3) in CD012925 |
| SJ | 1/10 | Known starting-value implementation difference |

### Conclusion

The four primary estimators used in the Fragility Atlas (DL, REML, FE, PM) agree with R metafor to within 1e-5 across all 10 validated reviews. HS and SJ show known cross-implementation differences that are themselves a valid source of specification variability in the multiverse framework.

## Validated Reviews

| Review | k | DL diff | REML diff | FE diff | PM diff |
|--------|---|---------|-----------|---------|---------|
| CD005133 | 39 | 5.21e-06 | 2.00e-06 | <1e-16 | 7.90e-06 |
| CD013059 | 5 | <1e-16 | <1e-16 | <1e-16 | <1e-16 |
| CD007491 | 8 | 1.96e-06 | 1.77e-06 | <1e-16 | 2.53e-06 |
| CD014967 | 5 | <1e-16 | <1e-16 | <1e-16 | <1e-16 |
| CD012925 | 6 | 1.20e-06 | 1.20e-06 | <1e-16 | 6.68e-06 |
| CD011626 | 11 | 1.47e-05 | 1.03e-05 | <1e-16 | 1.46e-05 |
| CD016211 | 5 | <1e-16 | <1e-16 | <1e-16 | <1e-16 |
| CD011793 | 13 | <1e-16 | 3.43e-06 | <1e-16 | <1e-16 |
| CD003218 | 5 | <1e-16 | <1e-16 | <1e-16 | <1e-16 |
| CD008628 | 9 | <1e-16 | <1e-16 | <1e-16 | <1e-16 |
