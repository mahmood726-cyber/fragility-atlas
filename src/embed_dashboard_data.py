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

    # Find and replace data blocks using string search (not regex — avoids backslash issues)
    import re

    # Replace SUMMARY_DATA
    match = re.search(r'const SUMMARY_DATA\s*=\s*', html)
    if match:
        start = match.start()
        # Find the matching closing brace + semicolon
        brace_start = html.index('{', match.end())
        depth = 0
        pos = brace_start
        while pos < len(html):
            if html[pos] == '{':
                depth += 1
            elif html[pos] == '}':
                depth -= 1
                if depth == 0:
                    end = html.index(';', pos) + 1
                    html = html[:start] + f'const SUMMARY_DATA = {summary_js};' + html[end:]
                    break
            pos += 1

    # Replace REVIEWS_DATA
    match = re.search(r'const REVIEWS_DATA\s*=\s*', html)
    if match:
        start = match.start()
        bracket_start = html.index('[', match.end())
        depth = 0
        pos = bracket_start
        while pos < len(html):
            if html[pos] == '[':
                depth += 1
            elif html[pos] == ']':
                depth -= 1
                if depth == 0:
                    end = html.index(';', pos) + 1
                    html = html[:start] + f'const REVIEWS_DATA = {reviews_js};' + html[end:]
                    break
            pos += 1

    Path(dashboard_path).write_text(html, encoding='utf-8')
    print(f"Dashboard data embedded: {len(reviews)} reviews, summary with {summary['total_specifications']:,} specs")


if __name__ == '__main__':
    embed()
