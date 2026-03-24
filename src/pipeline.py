"""End-to-end Fragility Atlas pipeline.

Usage: python -m src.pipeline [--pairwise-dir DIR] [--output-dir DIR] [--max-reviews N]
"""

import sys
import json
import time
import argparse
import csv
import numpy as np
from pathlib import Path
from dataclasses import asdict

from src.loader import load_all_reviews
from src.specifications import generate_specifications, SpecResult
from src.classifier import classify_review, ReviewClassification


DEFAULT_PAIRWISE_DIR = r'C:\Models\Pairwise70\data'
DEFAULT_OUTPUT_DIR = r'C:\FragilityAtlas\data\output'


def _process_review(args):
    """Process a single review (for multiprocessing)."""
    review, conf_level = args
    specs = generate_specifications(review, conf_level)
    classification = classify_review(review, specs)
    return classification, specs


def run_pipeline(pairwise_dir: str, output_dir: str, max_reviews: int = 0,
                 conf_level: float = 0.95, workers: int = 1):
    """Run the full Fragility Atlas pipeline."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Fragility Atlas Pipeline")
    print(f"========================")
    print(f"Data: {pairwise_dir}")
    print(f"Output: {output_dir}")
    if workers > 1:
        print(f"Workers: {workers}")
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

    if workers > 1:
        # Multiprocessing mode
        from concurrent.futures import ProcessPoolExecutor, as_completed
        tasks = [(review, conf_level) for review in reviews]
        done = 0
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_process_review, task): i
                       for i, task in enumerate(tasks)}
            # Collect results in submission order
            results_map = {}
            for future in as_completed(futures):
                idx = futures[future]
                classification, specs = future.result()
                results_map[idx] = (classification, specs)
                done += 1
                if done % 20 == 0 or done == len(reviews):
                    rate = done / (time.time() - t0)
                    eta = (len(reviews) - done) / rate if rate > 0 else 0
                    print(f"  [{done}/{len(reviews)}] ETA: {eta:.0f}s")
                    sys.stdout.flush()

        # Reassemble in original order
        for i in range(len(reviews)):
            classification, specs = results_map[i]
            all_classifications.append(classification)
            all_specs.extend(specs)
            total_specs += len(specs)
    else:
        # Sequential mode (original)
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

    results_path = output_path / 'fragility_atlas_results.csv'
    _export_classifications(all_classifications, results_path)
    print(f"  Review results: {results_path}")

    specs_path = output_path / 'fragility_atlas_specifications.csv'
    _export_specifications(all_specs, specs_path)
    print(f"  Specifications: {specs_path} ({len(all_specs):,} rows)")

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

    # P1-8 FIX: use numpy for proper median
    scores = np.array([c.robustness_score for c in all_classifications])
    mean_robustness = float(np.mean(scores)) if n > 0 else 0
    print(f"  Mean robustness: {mean_robustness:.1f}%")
    print(f"  Median robustness: {float(np.median(scores)):.1f}%")

    eta2_means = {
        'Estimator': float(np.mean([c.eta2_estimator for c in all_classifications])) if n > 0 else 0,
        'CI method': float(np.mean([c.eta2_ci_method for c in all_classifications])) if n > 0 else 0,
        'Bias correction': float(np.mean([c.eta2_bias_correction for c in all_classifications])) if n > 0 else 0,
        'Leave-one-out': float(np.mean([c.eta2_leave_out for c in all_classifications])) if n > 0 else 0,
    }
    top_dim = max(eta2_means, key=eta2_means.get)
    print(f"  Most influential dimension: {top_dim} (mean eta²={eta2_means[top_dim]:.3f})")
    print()

    return all_classifications, all_specs, summary


def _export_classifications(classifications, path):
    """Export review-level results to CSV."""
    if not classifications:
        return
    # P1-7 FIX: include top_dimension in export
    fields = [
        'review_id', 'review_doi', 'analysis_name', 'k', 'scale',
        'total_specs', 'agreeing_specs', 'robustness_score', 'classification',
        'cochrane_significant', 'cochrane_direction',
        'eta2_estimator', 'eta2_ci_method', 'eta2_bias_correction', 'eta2_leave_out',
        'top_dimension',
        'median_theta', 'iqr_theta', 'frac_significant', 'frac_reversed',
    ]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for c in classifications:
            row = asdict(c)
            # SE-P2-1 FIX: guard CSV injection for string fields
            for str_field in ('analysis_name', 'review_id', 'review_doi'):
                val = str(row.get(str_field, ''))
                if val and val[0] in ('=', '+', '@', '\t', '\r'):
                    row[str_field] = "'" + val
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

    # P1-8 FIX: use numpy for proper statistics
    robustness_scores = np.array([c.robustness_score for c in classifications])
    k_values = np.array([c.k for c in classifications])

    return {
        'n_reviews': n,
        'classification_counts': counts,
        'robustness_distribution': {
            'mean': round(float(np.mean(robustness_scores)), 2),
            'median': round(float(np.median(robustness_scores)), 2),
            'q25': round(float(np.percentile(robustness_scores, 25)), 2),
            'q75': round(float(np.percentile(robustness_scores, 75)), 2),
            'min': round(float(np.min(robustness_scores)), 2),
            'max': round(float(np.max(robustness_scores)), 2),
        },
        'k_distribution': {
            'mean': round(float(np.mean(k_values)), 1),
            'median': int(np.median(k_values)),
            'min': int(np.min(k_values)),
            'max': int(np.max(k_values)),
        },
        'eta2_means': {
            'estimator': round(float(np.mean([c.eta2_estimator for c in classifications])), 4),
            'ci_method': round(float(np.mean([c.eta2_ci_method for c in classifications])), 4),
            'bias_correction': round(float(np.mean([c.eta2_bias_correction for c in classifications])), 4),
            'leave_out': round(float(np.mean([c.eta2_leave_out for c in classifications])), 4),
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
    parser.add_argument('--workers', type=int, default=1,
                        help='Number of parallel workers (default: 1)')
    args = parser.parse_args()

    if not (0 < args.conf_level < 1):
        parser.error("conf-level must be between 0 and 1 (exclusive)")

    run_pipeline(args.pairwise_dir, args.output_dir, args.max_reviews,
                 args.conf_level, args.workers)


if __name__ == '__main__':
    main()
