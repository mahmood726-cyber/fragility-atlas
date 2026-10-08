"""Prevent a successful data refresh from leaving a broken or stale UI."""

import csv
import json
import re

import pytest

from src.embed_dashboard_data import _replace_variable, embed


def test_replacement_preserves_quoted_delimiters_and_escapes_script_end():
    html = 'var REVIEWS_DATA = [{"name": "] ; [ escaped \\\" quote"}];\nvar after = 1;'
    result = _replace_variable(html, 'REVIEWS_DATA', [{'name': '</script>'}])
    assert 'var after = 1;' in result
    assert '</script>' not in result
    assert json.loads(result.split(' = ', 1)[1].split(';\n', 1)[0]) == [{'name': '</script>'}]


def test_missing_variable_is_an_error():
    with pytest.raises(ValueError, match='missing'):
        _replace_variable('var other = 1;', 'SUMMARY_DATA', {})


def test_pipeline_summary_is_mapped_to_dashboard_display_fields(tmp_path):
    summary = {
        'n_reviews': 1,
        'classification_counts': {'Robust': 0, 'Moderate': 0, 'Fragile': 1,
                                  'Unstable': 0, 'Insufficient': 0},
        'robustness_distribution': {'mean': 66.67},
        'eta2_means': {'estimator': 0.1, 'ci_method': 0.2,
                       'bias_correction': 0.3, 'leave_out': 0.4},
        'total_specifications': 252,
    }
    (tmp_path / 'fragility_atlas_summary.json').write_text(json.dumps(summary))
    row = {
        'review_id': 'CD000001', 'analysis_name': 'A [bracket] and ; semicolon',
        'k': 3, 'total_specs': 252, 'agreeing_specs': 168,
        'robustness_score': 66.67, 'cochrane_significant': True,
        'cochrane_direction': -1,
        **{key: 0.1 for key in ['eta2_estimator', 'eta2_ci_method',
                                'eta2_bias_correction', 'eta2_leave_out',
                                'median_theta', 'iqr_theta', 'frac_significant',
                                'frac_reversed']},
    }
    with (tmp_path / 'fragility_atlas_results.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    dashboard = tmp_path / 'index.html'
    dashboard.write_text('\n'.join(f'var {name} = {{}};' for name in [
        'SUMMARY_DATA', 'REVIEWS_DATA', 'TOTAL_REVIEWS', 'counts',
        'meanRobustness', 'pctFragileUnstable', 'totalSpecs']))
    embed(tmp_path, dashboard)
    html = dashboard.read_text()
    match = re.search(r'var SUMMARY_DATA = (.*);', html)
    displayed = json.loads(match.group(1))
    assert displayed['total_reviews'] == 1
    assert displayed['mean_robustness'] == 66.67
    assert displayed['pct_fragile_unstable'] == 100.0
    assert displayed['mean_eta']['Bias Correction'] == 0.3
    assert 'var totalSpecs = 252;' in html
    assert 'var TOTAL_REVIEWS = 1;' in html

    summary['n_reviews'] = 2
    (tmp_path / 'fragility_atlas_summary.json').write_text(json.dumps(summary))
    previous = dashboard.read_bytes()
    with pytest.raises(ValueError, match='count does not match'):
        embed(tmp_path, dashboard)
    assert dashboard.read_bytes() == previous
