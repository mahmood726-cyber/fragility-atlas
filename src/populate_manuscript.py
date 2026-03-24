"""Auto-populate the BMJ manuscript with results from the pipeline output."""

import json
import csv
import numpy as np
from pathlib import Path


def populate(output_dir: str = r'C:\FragilityAtlas\data\output',
             manuscript_path: str = r'C:\FragilityAtlas\manuscript_bmj.md'):
    """Read pipeline outputs and populate manuscript placeholders."""

    output = Path(output_dir)

    # Load summary
    with open(output / 'fragility_atlas_summary.json') as f:
        summary = json.load(f)

    # Load classifications
    reviews = []
    with open(output / 'fragility_atlas_results.csv', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            reviews.append(row)

    n = summary.get('n_reviews', 0)
    if n == 0:
        print("No reviews in pipeline output — cannot populate manuscript.")
        return
    counts = summary.get('classification_counts', {})
    rob = summary.get('robustness_distribution', {})
    k_dist = summary.get('k_distribution', {})
    eta2 = summary.get('eta2_means', {})

    # Compute scale breakdown
    n_ratio = summary['scale_counts']['ratio']
    n_diff = summary['scale_counts']['difference']

    # Find most/least robust reviews
    reviews_sorted = sorted(reviews, key=lambda r: float(r['robustness_score']))
    most_robust = reviews_sorted[-1]
    least_robust = reviews_sorted[0]

    # Dimension ranking
    dim_sorted = sorted(eta2.items(), key=lambda x: x[1], reverse=True)
    top_dim_name = dim_sorted[0][0].replace('_', ' ').title()
    top_dim_eta2 = dim_sorted[0][1]
    second_dim_name = dim_sorted[1][0].replace('_', ' ').title()
    second_dim_eta2 = dim_sorted[1][1]

    # Replacements
    replacements = {
        '[N_REVIEWS]': str(n),
        '[TOTAL_SPECS]': f"{summary['total_specifications']:,}",
        '[MEDIAN_K]': str(k_dist['median']),
        '[K_Q25]': str(int(np.percentile([int(r['k']) for r in reviews], 25))),
        '[K_Q75]': str(int(np.percentile([int(r['k']) for r in reviews], 75))),
        '[K_MIN]': str(k_dist['min']),
        '[K_MAX]': str(k_dist['max']),
        '[N_RATIO]': str(n_ratio),
        '[PCT_RATIO]': f"{n_ratio/n*100:.1f}",
        '[N_DIFF]': str(n_diff),
        '[PCT_DIFF]': f"{n_diff/n*100:.1f}",
        '[MEAN_ROBUST]': f"{rob['mean']:.1f}",
        '[MEDIAN_ROBUST]': f"{rob['median']:.1f}",
        '[Q25_ROBUST]': f"{rob['q25']:.1f}",
        '[Q75_ROBUST]': f"{rob['q75']:.1f}",
        '[N_ROBUST]': str(counts.get('Robust', 0)),
        '[PCT_ROBUST]': f"{counts.get('Robust', 0)/n*100:.1f}",
        '[N_MODERATE]': str(counts.get('Moderate', 0)),
        '[PCT_MODERATE]': f"{counts.get('Moderate', 0)/n*100:.1f}",
        '[N_FRAGILE]': str(counts.get('Fragile', 0)),
        '[PCT_FRAGILE]': f"{counts.get('Fragile', 0)/n*100:.1f}",
        '[N_UNSTABLE]': str(counts.get('Unstable', 0)),
        '[PCT_UNSTABLE]': f"{counts.get('Unstable', 0)/n*100:.1f}",
        '[TOP_DIM]': top_dim_name,
        '[TOP_ETA2]': f"{top_dim_eta2:.3f}",
        '[SECOND_DIM]': second_dim_name,
        '[SECOND_ETA2]': f"{second_dim_eta2:.3f}",
        '[MEDIAN_SPECS]': str(63 * (k_dist['median'] + 1)),
    }

    # Results paragraph
    results_para = (
        f"Of {n} eligible Cochrane reviews, "
        f"{counts.get('Robust', 0)} ({counts.get('Robust', 0)/n*100:.1f}%) were classified as Robust, "
        f"{counts.get('Moderate', 0)} ({counts.get('Moderate', 0)/n*100:.1f}%) as Moderately Robust, "
        f"{counts.get('Fragile', 0)} ({counts.get('Fragile', 0)/n*100:.1f}%) as Fragile, and "
        f"{counts.get('Unstable', 0)} ({counts.get('Unstable', 0)/n*100:.1f}%) as Unstable. "
        f"The mean robustness score was {rob['mean']:.1f}% (median {rob['median']:.1f}%, "
        f"IQR {rob['q25']:.1f}-{rob['q75']:.1f}%). "
        f"The most influential analytical dimension was {top_dim_name.lower()} "
        f"(mean eta-squared = {top_dim_eta2:.3f}), followed by {second_dim_name.lower()} "
        f"(eta-squared = {second_dim_eta2:.3f}). "
        f"A total of {summary['total_specifications']:,} individual meta-analytic specifications "
        f"were computed across all reviews."
    )
    replacements['[RESULTS_PARAGRAPH \u2014 to be filled after pipeline completes]'] = results_para

    conclusions_para = (
        f"In this large-scale multiverse analysis of {n} Cochrane meta-analyses, "
        f"we found that {(counts.get('Fragile', 0) + counts.get('Unstable', 0))/n*100:.0f}% of reviews "
        f"had conclusions that were sensitive to reasonable analytical choices. "
        f"The {top_dim_name.lower()} dimension was the most influential driver of disagreement. "
        f"These findings suggest that the conclusions of evidence-based medicine are more "
        f"contingent on methodological choices than is commonly appreciated, and highlight "
        f"the importance of conducting sensitivity analyses across multiple analytical dimensions."
    )
    replacements['[CONCLUSIONS_PARAGRAPH \u2014 to be filled after pipeline completes]'] = conclusions_para

    # Read and replace
    text = Path(manuscript_path).read_text(encoding='utf-8')
    for placeholder, value in replacements.items():
        text = text.replace(placeholder, value)

    Path(manuscript_path).write_text(text, encoding='utf-8')
    print(f"Manuscript populated with {len(replacements)} values.")
    print(f"Key findings:")
    print(f"  Reviews: {n}")
    print(f"  Robust: {counts.get('Robust', 0)} ({counts.get('Robust', 0)/n*100:.1f}%)")
    print(f"  Moderate: {counts.get('Moderate', 0)} ({counts.get('Moderate', 0)/n*100:.1f}%)")
    print(f"  Fragile: {counts.get('Fragile', 0)} ({counts.get('Fragile', 0)/n*100:.1f}%)")
    print(f"  Unstable: {counts.get('Unstable', 0)} ({counts.get('Unstable', 0)/n*100:.1f}%)")
    print(f"  Mean robustness: {rob['mean']:.1f}%")
    print(f"  Top dimension: {top_dim_name} (eta2={top_dim_eta2:.3f})")


if __name__ == '__main__':
    populate()
