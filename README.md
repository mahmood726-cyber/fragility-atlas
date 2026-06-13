# The Fragility Atlas

[![ci](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/ci.yml) [![codeql](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/fragility-atlas/actions/workflows/codeql.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)

How robust are Cochrane meta-analysis conclusions when analysts make different but equally reasonable methodological choices? We applied multiverse analysis to 403 Cochrane systematic reviews from the Pairwise70 dataset, encompassing 394,569 individual meta-analytic specifications across five analytical dimensions. Each review was re-analysed across seven variance estimators, three confidence interval methods, three publication bias corrections, and leave-one-out subsets, then classified by agreement with the standard DerSimonian-Laird reference conclusion. The mean robustness score was 66.4 percent (median 66.7 percent, IQR 56.5 to 87.2 percent), with 58 percent of all reviews classified as Fragile or Unstable. Publication bias correction was the single most influential dimension, explaining ten times more variance (eta-squared 0.374) than variance estimator choice (eta-squared 0.036). These findings suggest that most Cochrane conclusions are substantially sensitive to defensible analytical alternatives, particularly the decision whether to adjust for publication bias. This analysis is limited to binary agreement classification and cannot capture the magnitude of effect size changes across specifications.

**Live dashboard:** <https://mahmood726-cyber.github.io/fragilityatlas/>

## Run

Open `index.html` (or `index.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q
```

The suite under `tests/` includes 3 test file(s).

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
