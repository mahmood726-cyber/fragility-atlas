"""Embed pipeline results into the dashboard HTML file."""

import json
import csv
from pathlib import Path


def embed(output_dir: str = r'C:\FragilityAtlas\data\output',
          dashboard_path: str = r'C:\FragilityAtlas\dashboard\index.html'):
    """Replace placeholder data in dashboard with real pipeline results."""

    output = Path(output_dir)

    # Load summary
    with open(output / 'fragility_atlas_summary.json') as f:
        summary = json.load(f)

    # Load reviews CSV
    reviews = []
    with open(output / 'fragility_atlas_results.csv', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for key in ['k', 'total_specs', 'agreeing_specs']:
                row[key] = int(row[key])
            for key in ['robustness_score', 'eta2_estimator', 'eta2_ci_method',
                        'eta2_bias_correction', 'eta2_leave_out',
                        'median_theta', 'iqr_theta', 'frac_significant', 'frac_reversed']:
                row[key] = float(row[key])
            row['cochrane_significant'] = row['cochrane_significant'] == 'True'
            row['cochrane_direction'] = int(row['cochrane_direction'])
            reviews.append(row)

    # Generate JS data
    summary_js = json.dumps(summary, indent=2)
    reviews_js = json.dumps(reviews, indent=None)  # compact for size

    # Read dashboard HTML
    html = Path(dashboard_path).read_text(encoding='utf-8')

    # Replace placeholder data blocks
    # Look for: const SUMMARY_DATA = {...}; or similar patterns
    import re

    # Replace SUMMARY_DATA
    html = re.sub(
        r'const SUMMARY_DATA\s*=\s*\{[^;]*\};',
        f'const SUMMARY_DATA = {summary_js};',
        html, flags=re.DOTALL
    )

    # Replace REVIEWS_DATA
    html = re.sub(
        r'const REVIEWS_DATA\s*=\s*\[[^;]*\];',
        f'const REVIEWS_DATA = {reviews_js};',
        html, flags=re.DOTALL
    )

    Path(dashboard_path).write_text(html, encoding='utf-8')
    print(f"Dashboard data embedded: {len(reviews)} reviews, summary with {summary['total_specifications']:,} specs")


if __name__ == '__main__':
    embed()
