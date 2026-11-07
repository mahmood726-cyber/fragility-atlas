"""Embed pipeline results into the dashboard HTML file."""

import csv
import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / 'data' / 'output'
DEFAULT_DASHBOARD_PATH = REPO_ROOT / 'dashboard' / 'index.html'


def _replace_variable(html, name, value):
    """Replace one JS assignment without treating quoted brackets as syntax."""
    match = re.search(r'(?:const|var)\s+' + re.escape(name) + r'\s*=\s*', html)
    if not match:
        raise ValueError(f'Dashboard variable missing: {name}')
    depth = 0
    quote = None
    escaped = False
    for pos in range(match.end(), len(html)):
        char = html[pos]
        if quote:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == quote:
                quote = None
        elif char in ('"', "'"):
            quote = char
        elif char in '{[(':
            depth += 1
        elif char in '}])':
            depth -= 1
        elif char == ';' and depth == 0:
            # Prevent data strings from terminating an inline script element.
            encoded = json.dumps(value, ensure_ascii=True).replace('<', '\\u003c')
            return html[:match.start()] + f'var {name} = {encoded};' + html[pos + 1:]
    raise ValueError(f'Dashboard assignment is unterminated: {name}')


def embed(output_dir: str | Path = DEFAULT_OUTPUT_DIR,
          dashboard_path: str | Path = DEFAULT_DASHBOARD_PATH):
    """Replace placeholder data in dashboard with real pipeline results."""

    output = Path(output_dir)

    # Load summary
    with open(output / 'fragility_atlas_summary.json', encoding='utf-8') as f:
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

    # Keep the dashboard's display schema in sync with the pipeline's export
    # schema. Dumping the raw summary used to leave the UI fields undefined.
    counts = summary['classification_counts']
    n_reviews = summary['n_reviews']
    display_summary = {
        **summary,
        'total_reviews': n_reviews,
        'mean_robustness': summary['robustness_distribution']['mean'],
        'pct_fragile_unstable': round(
            100 * (counts['Fragile'] + counts['Unstable']) / n_reviews, 1
        ) if n_reviews else 0.0,
        'counts': counts,
        'mean_eta': {
            'Estimator': summary['eta2_means']['estimator'],
            'CI Method': summary['eta2_means']['ci_method'],
            'Bias Correction': summary['eta2_means']['bias_correction'],
            'Leave-One-Out': summary['eta2_means']['leave_out'],
        },
    }
    if len(reviews) != n_reviews:
        raise ValueError('Review count does not match summary; refusing partial refresh')

    # Read dashboard HTML
    dashboard = Path(dashboard_path)
    html = dashboard.read_text(encoding='utf-8')

    assignments = {
        'SUMMARY_DATA': display_summary,
        'REVIEWS_DATA': reviews,
        'TOTAL_REVIEWS': n_reviews,
        'counts': counts,
        'meanRobustness': display_summary['mean_robustness'],
        'pctFragileUnstable': display_summary['pct_fragile_unstable'],
        'totalSpecs': summary['total_specifications'],
    }
    for name, value in assignments.items():
        html = _replace_variable(html, name, value)

    dashboard.write_text(html, encoding='utf-8')
    print(f"Dashboard data embedded: {len(reviews)} reviews, summary with {summary['total_specifications']:,} specs")


if __name__ == '__main__':
    embed()
