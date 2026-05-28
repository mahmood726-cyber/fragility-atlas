"""Export yi/sei vectors from 10 selected reviews for R cross-validation."""
import csv
import json
import os
import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
from src.loader import load_review

OUTPUT_DIR = REPO_ROOT / "data" / "output"
OUTPUT = OUTPUT_DIR / "r_validation_inputs.json"
RESULTS_CSV = OUTPUT_DIR / "fragility_atlas_results.csv"


def resolve_pairwise_dir() -> Path:
    candidates = []
    env_root = os.environ.get("PAIRWISE70_DATA_DIR", "")
    if env_root:
        candidates.append(Path(env_root).expanduser())
    candidates.extend([
        REPO_ROOT.parent / "Projects" / "mahmood789" / "Pairwise70" / "data",
        REPO_ROOT.parent / "Projects" / "Pairwise70" / "data",
        REPO_ROOT.parent / "Models" / "Pairwise70" / "data",
        REPO_ROOT.parent / "Pairwise70" / "data",
    ])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(
        "Pairwise70 data directory not found. Set PAIRWISE70_DATA_DIR or restore one of: "
        + ", ".join(str(candidate) for candidate in candidates)
    )


PAIRWISE_DIR = resolve_pairwise_dir()

# Same seed as R script
random.seed(42)

# Load results to get eligible reviews (k>=5, ratio scale)
eligible = []
with open(RESULTS_CSV, encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f):
        if int(row["k"]) >= 5 and row["scale"] == "ratio":
            eligible.append(row["review_id"])

sample_ids = random.sample(eligible, min(10, len(eligible)))
print(f"Selected {len(sample_ids)} reviews: {sample_ids}")

exports = []
for rid in sample_ids:
    rda_files = list(PAIRWISE_DIR.glob(f"{rid}_*"))
    if not rda_files:
        print(f"  SKIP {rid}: no RDA file")
        continue

    review = load_review(str(rda_files[0]))
    if review is None:
        print(f"  SKIP {rid}: load failed")
        continue

    exports.append({
        "review_id": rid,
        "analysis_name": review.analysis_name,
        "k": review.k,
        "scale": review.scale,
        "yi": review.yi.tolist(),
        "sei": review.sei.tolist(),
    })
    print(f"  {rid}: k={review.k}, analysis={review.analysis_name}")

with open(OUTPUT, "w") as f:
    json.dump(exports, f, indent=2)

print(f"\nExported {len(exports)} reviews to {OUTPUT}")
