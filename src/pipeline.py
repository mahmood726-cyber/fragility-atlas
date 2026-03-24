"""End-to-end Fragility Atlas pipeline.

Usage: python -m src.pipeline [--pairwise-dir DIR] [--output-dir DIR] [--max-reviews N]
"""

import sys
import json
import time
import argparse
import csv
from pathlib import Path
from dataclasses import asdict

from src.loader import load_all_reviews
from src.specifications import generate_specifications, SpecResult
from src.classifier import classify_review, ReviewClassification


DEFAULT_PAIRWISE_DIR = r'C:\Models\Pairwise70\data'
DEFAULT_OUTPUT_DIR = r'C:\FragilityAtlas\data\output'


def run_pipeline(pairwise_dir: str, output_dir: str, max_reviews: int = 0,
                 conf_level: float = 0.95):
    """Run the full Fragility Atlas pipeline."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Fragility Atlas Pipeline")
    print(f"========================")
    print(f"Data: {pairwise_dir}")
    print(f"Output: {output_dir}")
    print()

    # Phase 1: Load reviews
    print("Phase 1: Loading reviews...")
    reviews = list(load_all_reviews(pairwise_dir, min_k=3))
    if max_reviews > 0:
        reviews = reviews[:max_reviews]
    print(f"  Loaded {len(reviews)} eligible reviews (k >= 3)")
    print()

    # Phase 2: Generate and execute specifications
    print("Phase 2: Running multiverse analysis...")
    all_classifications = []
    all_specs = []
    total_specs = 0

    t0 = time.time()
    for i, review in enumerate(reviews):
        t_review = time.time()
        specs = generate_specifications(review, conf_level)
        classification = classify_review(review, specs)

        all_classifications.append(classification)
        all_specs.extend(specs)
        total_specs += len(specs)

        elapsed = time.time() - t_review
        if (i + 1) % 10 == 0 or (i + 1) == len(reviews):
            rate = (i + 1) / (time.time() - t0)
            eta = (len(reviews) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(reviews)}] {review.review_id}: k={review.k}, "
                  f"specs={len(specs)}, robustness={classification.robustness_score:.1f}% "
                  f"({classification.classification}) [{elapsed:.1f}s] ETA: {eta:.0f}s")
            sys.stdout.flush()

    total_time = time.time() - t0
    print(f"\n  Total: {total_specs:,} specifications in {total_time:.1f}s")
    print()

    # Phase 3: Export results
    print("Phase 3: Exporting results...")

    # 3a: Review-level classifications
    results_path = output_path / 'fragility_atlas_results.csv'
    _export_classifications(all_classifications, results_path)
    print(f"  Review results: {results_path}")

    # 3b: Specification-level detail (can be large)
    specs_path = output_path / 'fragility_atlas_specifications.csv'
    _export_specifications(all_specs, specs_path)
    print(f"  Specifications: {specs_path} ({len(all_specs):,} rows)")

    # 3c: Summary JSON for dashboard
    summary = _compute_summary(all_classifications, total_time)
    summary_path = output_path / 'fragility_atlas_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Summary: {summary_path}")

    # Print headline results
    print()
    print("=" * 60)
    print("HEADLINE RESULTS")
    print("=" * 60)
    n = len(all_classifications)
    for cat in ['Robust', 'Moderate', 'Fragile', 'Unstable']:
        count = sum(1 for c in all_classifications if c.classification == cat)
        pct = count / n * 100 if n > 0 else 0
        print(f"  {cat:12s}: {count:4d} ({pct:5.1f}%)")
    print()

    mean_robustness = sum(c.robustness_score for c in all_classifications) / n if n > 0 else 0
    print(f"  Mean robustness: {mean_robustness:.1f}%")
    print(f"  Median robustness: {sorted(c.robustness_score for c in all_classifications)[n//2]:.1f}%")

    # Most influential dimension
    eta2_means = {
        'Estimator': sum(c.eta2_estimator for c in all_classifications) / n if n > 0 else 0,
        'CI method': sum(c.eta2_ci_method for c in all_classifications) / n if n > 0 else 0,
        'Bias correction': sum(c.eta2_bias_correction for c in all_classifications) / n if n > 0 else 0,
        'Leave-one-out': sum(c.eta2_leave_out for c in all_classifications) / n if n > 0 else 0,
    }
    top_dim = max(eta2_means, key=eta2_means.get)
    print(f"  Most influential dimension: {top_dim} (mean eta²={eta2_means[top_dim]:.3f})")
    print()

    return all_classifications, all_specs, summary


def _export_classifications(classifications, path):
    """Export review-level results to CSV."""
    if not classifications:
        return
    fields = [
        'review_id', 'review_doi', 'analysis_name', 'k', 'scale',
        'total_specs', 'agreeing_specs', 'robustness_score', 'classification',
        'cochrane_significant', 'cochrane_direction',
        'eta2_estimator', 'eta2_ci_method', 'eta2_bias_correction', 'eta2_leave_out',
        'median_theta', 'iqr_theta', 'frac_significant', 'frac_reversed',
    ]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for c in classifications:
            row = asdict(c)
            writer.writerow({k: row[k] for k in fields})


def _export_specifications(specs, path):
    """Export specification-level results to CSV."""
    if not specs:
        return
    fields = [
        'review_id', 'estimator', 'ci_method', 'bias_correction', 'leave_out',
        'theta', 'se_theta', 'ci_lo', 'ci_hi', 'p_value',
        'tau2', 'i2', 'is_significant', 'direction',
    ]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for s in specs:
            row = asdict(s)
            writer.writerow({k: row[k] for k in fields})


def _compute_summary(classifications, elapsed_time):
    """Compute summary statistics for the dashboard."""
    n = len(classifications)
    if n == 0:
        return {'n_reviews': 0}

    counts = {}
    for cat in ['Robust', 'Moderate', 'Fragile', 'Unstable', 'Insufficient']:
        counts[cat] = sum(1 for c in classifications if c.classification == cat)

    robustness_scores = [c.robustness_score for c in classifications]
    k_values = [c.k for c in classifications]

    return {
        'n_reviews': n,
        'classification_counts': counts,
        'robustness_distribution': {
            'mean': round(sum(robustness_scores) / n, 2),
            'median': round(sorted(robustness_scores)[n // 2], 2),
            'q25': round(sorted(robustness_scores)[n // 4], 2),
            'q75': round(sorted(robustness_scores)[3 * n // 4], 2),
            'min': round(min(robustness_scores), 2),
            'max': round(max(robustness_scores), 2),
        },
        'k_distribution': {
            'mean': round(sum(k_values) / n, 1),
            'median': sorted(k_values)[n // 2],
            'min': min(k_values),
            'max': max(k_values),
        },
        'eta2_means': {
            'estimator': round(sum(c.eta2_estimator for c in classifications) / n, 4),
            'ci_method': round(sum(c.eta2_ci_method for c in classifications) / n, 4),
            'bias_correction': round(sum(c.eta2_bias_correction for c in classifications) / n, 4),
            'leave_out': round(sum(c.eta2_leave_out for c in classifications) / n, 4),
        },
        'scale_counts': {
            'ratio': sum(1 for c in classifications if c.scale == 'ratio'),
            'difference': sum(1 for c in classifications if c.scale == 'difference'),
        },
        'total_specifications': sum(c.total_specs for c in classifications),
        'elapsed_seconds': round(elapsed_time, 1),
    }


def main():
    parser = argparse.ArgumentParser(description='Fragility Atlas Pipeline')
    parser.add_argument('--pairwise-dir', default=DEFAULT_PAIRWISE_DIR,
                        help='Path to Pairwise70 RDA files')
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR,
                        help='Output directory')
    parser.add_argument('--max-reviews', type=int, default=0,
                        help='Max reviews to process (0 = all)')
    parser.add_argument('--conf-level', type=float, default=0.95,
                        help='Confidence level (default: 0.95)')
    args = parser.parse_args()

    run_pipeline(args.pairwise_dir, args.output_dir, args.max_reviews, args.conf_level)


if __name__ == '__main__':
    main()
