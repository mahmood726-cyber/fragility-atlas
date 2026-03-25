"""Export yi/sei vectors from 10 selected reviews for R cross-validation."""
import csv
import json
import random
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.loader import load_review

PAIRWISE_DIR = Path("C:/Models/Pairwise70/data")
OUTPUT = Path("C:/FragilityAtlas/data/output/r_validation_inputs.json")

# Same seed as R script
random.seed(42)

# Load results to get eligible reviews (k>=5, ratio scale)
results_csv = Path("C:/FragilityAtlas/data/output/fragility_atlas_results.csv")
eligible = []
with open(results_csv, encoding="utf-8", errors="replace") as f:
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
