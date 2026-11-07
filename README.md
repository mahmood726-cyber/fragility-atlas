# The Fragility Atlas

[![ci](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/ci.yml) [![codeql](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/codeql.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)

Fragility Atlas explores how pooled estimates and significance decisions vary across seven variance estimators, three interval methods, three bias-sensitivity options and leave-one-out subsets. It compares specifications with an internally recomputed DerSimonian-Laird/Wald reference, not with a review's published conclusion.

**Numerical correction (8 October 2026):** HS, HE and SJ implementation errors are fixed, and REML now uses a profile-likelihood fallback when Fisher scoring does not converge. External reference and real-review regression checks are documented in [the validation record](docs/numerical_validation_2026-10-08.md). Test failures now fail CI.

The regenerated exploratory run contains 473 eligible reviews and 469,350 specifications. The March manuscript and micro-paper use a different 403-review snapshot and are superseded pending scientific revalidation. Automatic outcome selection, scale detection, the specification grid and agreement thresholds require expert review before research or clinical reuse.

**Live dashboard:** <https://mahmood726-cyber.github.io/fragility-atlas/>

## Run

Open `index.html` (or `dashboard/index.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q
```

The suite includes external estimator references, real-review REML regressions and dashboard refresh checks. Set `PAIRWISE70_DATA_DIR` to include the real-data loader test.

## Repo layout

| Path | Purpose |
|---|---|
| `dashboard/index.html` | the dashboard (main artifact) |
| `index.html` | landing page |
| `tests/` | pytest tests |
| `e156-submission/` | E156 micro-paper bundle |
| `E156-PROTOCOL.md` | project metadata (E156 entry #60) |

## License

See `LICENSE` (MIT).
