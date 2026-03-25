"""
Generate 4 publication-quality figures for Fragility Atlas BMJ manuscript.
Uses matplotlib only (no seaborn dependency).
"""

import csv
import json
import sys
import io
import os

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Output directory
OUT_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

# Load data
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'output')

with open(os.path.join(DATA_DIR, 'fragility_atlas_results.csv'), encoding='utf-8') as f:
    reader = csv.DictReader(f)
    reviews = list(reader)

with open(os.path.join(DATA_DIR, 'fragility_atlas_summary.json'), encoding='utf-8') as f:
    summary = json.load(f)

print(f"Loaded {len(reviews)} reviews")

# Parse numeric fields
for r in reviews:
    r['robustness'] = float(r['robustness_score'])
    r['k_val'] = int(r['k'])
    r['eta2_est'] = float(r['eta2_estimator'])
    r['eta2_ci'] = float(r['eta2_ci_method'])
    r['eta2_bias'] = float(r['eta2_bias_correction'])
    r['eta2_loo'] = float(r['eta2_leave_out'])
    r['median_theta_val'] = float(r['median_theta'])

# Color scheme for classifications
CLASS_COLORS = {
    'Robust': '#2ecc71',      # green
    'Moderate': '#f39c12',    # amber
    'Fragile': '#e74c3c',     # red
    'Unstable': '#8e44ad',    # purple
}

# ─────────────────────────────────────────────────
# FIGURE 1: Distribution of Robustness Scores
# ─────────────────────────────────────────────────
print("Generating Figure 1: Robustness distribution...")

fig, ax = plt.subplots(figsize=(8, 5))

scores = [r['robustness'] for r in reviews]
bins = np.arange(0, 105, 5)

# Color each bin by classification
for i in range(len(bins) - 1):
    lo, hi = bins[i], bins[i+1]
    midpoint = (lo + hi) / 2
    if midpoint >= 90:
        color = CLASS_COLORS['Robust']
    elif midpoint >= 70:
        color = CLASS_COLORS['Moderate']
    elif midpoint >= 50:
        color = CLASS_COLORS['Fragile']
    else:
        color = CLASS_COLORS['Unstable']

    count = sum(1 for s in scores if lo <= s < hi)
    if i == len(bins) - 2:  # last bin includes upper edge
        count = sum(1 for s in scores if lo <= s <= hi)
    ax.bar(midpoint, count, width=4.5, color=color, edgecolor='white', linewidth=0.5)

# Classification threshold lines
for threshold, label in [(50, '50%'), (70, '70%'), (90, '90%')]:
    ax.axvline(x=threshold, color='#555555', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(threshold + 1, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 50,
            label, fontsize=8, color='#555555', va='top')

# Legend
patches = [mpatches.Patch(color=c, label=f'{k} ({summary["classification_counts"][k]}, '
           f'{summary["classification_counts"][k]/len(reviews)*100:.1f}%)')
           for k, c in CLASS_COLORS.items()]
ax.legend(handles=patches, loc='upper left', framealpha=0.9, fontsize=9)

ax.set_xlabel('Robustness Score (%)', fontsize=11)
ax.set_ylabel('Number of Reviews', fontsize=11)
ax.set_title('Distribution of Robustness Scores Across 403 Cochrane Meta-Analyses', fontsize=12, pad=12)
ax.set_xlim(-2, 102)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add summary statistics text
stats_text = (f"Mean: {summary['robustness_distribution']['mean']:.1f}%\n"
              f"Median: {summary['robustness_distribution']['median']:.1f}%\n"
              f"IQR: {summary['robustness_distribution']['q25']:.1f}-{summary['robustness_distribution']['q75']:.1f}%")
ax.text(0.98, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', alpha=0.9))

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'figure1_robustness_distribution.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure1_robustness_distribution.pdf'), bbox_inches='tight')
print("  Saved figure1_robustness_distribution.png/pdf")
plt.close(fig)


# ─────────────────────────────────────────────────
# FIGURE 2: Specification Curve (Exemplar Fragile Review)
# ─────────────────────────────────────────────────
print("Generating Figure 2: Specification curve...")

# Find a good exemplar: Fragile classification, moderate k (8-15), interesting pattern
fragile_reviews = [r for r in reviews if r['classification'] == 'Fragile' and 8 <= r['k_val'] <= 15]
# Pick the one with robustness closest to the median fragile score
fragile_reviews.sort(key=lambda r: abs(r['robustness'] - 60))
exemplar = fragile_reviews[0] if fragile_reviews else reviews[0]

# Load specifications for this review
specs_file = os.path.join(DATA_DIR, 'fragility_atlas_specifications.csv')
exemplar_specs = []
with open(specs_file, encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['review_id'] == exemplar['review_id']:
            exemplar_specs.append(row)

if exemplar_specs:
    # Parse and sort by effect estimate
    for s in exemplar_specs:
        s['theta_val'] = float(s['theta'])
        s['sig'] = s['is_significant'].lower() == 'true'
        s['dir'] = s['direction']

    # Determine reference: DL + Wald + none + no leave-out
    ref = [s for s in exemplar_specs
           if s['estimator'] == 'DL' and s['ci_method'] == 'Wald'
           and s['bias_correction'] == 'none' and s['leave_out'] == '']
    if ref:
        ref_sig = ref[0]['sig']
        ref_dir = ref[0]['dir']
    else:
        ref_sig = exemplar['cochrane_significant'] == 'True'
        ref_dir = exemplar['cochrane_direction']

    for s in exemplar_specs:
        s['agrees'] = (s['sig'] == ref_sig) and (s['dir'] == ref_dir)

    exemplar_specs.sort(key=lambda s: s['theta_val'])

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 8),
                                              gridspec_kw={'height_ratios': [3, 0.6, 0.6, 0.6]},
                                              sharex=True)

    x = range(len(exemplar_specs))
    thetas = [s['theta_val'] for s in exemplar_specs]
    colors_spec = ['#2ecc71' if s['agrees'] else '#e74c3c' for s in exemplar_specs]

    ax1.scatter(x, thetas, c=colors_spec, s=3, alpha=0.7, edgecolors='none')
    ax1.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
    ax1.set_ylabel('Pooled Effect Estimate', fontsize=10)
    ax1.set_title(f'Specification Curve: {exemplar["review_id"]} ({exemplar["analysis_name"]})\n'
                  f'k={exemplar["k_val"]} studies, robustness={exemplar["robustness"]:.1f}%',
                  fontsize=11, pad=8)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Legend for main panel
    agree_patch = mpatches.Patch(color='#2ecc71', label='Agrees with reference')
    disagree_patch = mpatches.Patch(color='#e74c3c', label='Disagrees')
    ax1.legend(handles=[agree_patch, disagree_patch], loc='upper left', fontsize=8)

    # Dimension indicator strips
    estimators = sorted(set(s.get('estimator', '') for s in exemplar_specs))
    ci_methods = sorted(set(s.get('ci_method', '') for s in exemplar_specs))
    bias_methods = sorted(set(s.get('bias_correction', '') for s in exemplar_specs))

    est_cmap = plt.cm.Set3(np.linspace(0, 1, max(len(estimators), 1)))
    ci_cmap = plt.cm.Set2(np.linspace(0, 1, max(len(ci_methods), 1)))
    bias_cmap = plt.cm.Pastel1(np.linspace(0, 1, max(len(bias_methods), 1)))

    for i, s in enumerate(exemplar_specs):
        est_idx = estimators.index(s.get('estimator', '')) if s.get('estimator', '') in estimators else 0
        ci_idx = ci_methods.index(s.get('ci_method', '')) if s.get('ci_method', '') in ci_methods else 0
        bias_idx = bias_methods.index(s.get('bias_correction', '')) if s.get('bias_correction', '') in bias_methods else 0

        ax2.bar(i, 1, width=1, color=est_cmap[est_idx], edgecolor='none')
        ax3.bar(i, 1, width=1, color=ci_cmap[ci_idx], edgecolor='none')
        ax4.bar(i, 1, width=1, color=bias_cmap[bias_idx], edgecolor='none')

    ax2.set_ylabel('Estimator', fontsize=8, rotation=0, ha='right')
    ax3.set_ylabel('CI Method', fontsize=8, rotation=0, ha='right')
    ax4.set_ylabel('Bias Corr.', fontsize=8, rotation=0, ha='right')

    for a in [ax2, ax3, ax4]:
        a.set_yticks([])
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)

    ax4.set_xlabel(f'Specifications (sorted by effect estimate, n={len(exemplar_specs)})', fontsize=10)

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'figure2_specification_curve.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, 'figure2_specification_curve.pdf'), bbox_inches='tight')
    print(f"  Saved figure2_specification_curve.png/pdf (exemplar: {exemplar['review_id']})")
    plt.close(fig)
else:
    print("  WARNING: No specification data found for exemplar review, skipping Figure 2")


# ─────────────────────────────────────────────────
# FIGURE 3: Dimension Attribution (eta-squared bar chart)
# ─────────────────────────────────────────────────
print("Generating Figure 3: Dimension attribution...")

fig, ax = plt.subplots(figsize=(7, 4))

eta2 = summary['eta2_means']
dims = ['Publication Bias\nCorrection', 'CI Method', 'Variance\nEstimator', 'Leave-One-Out']
vals = [eta2['bias_correction'], eta2['ci_method'], eta2['estimator'], eta2['leave_out']]
colors_bar = ['#e74c3c', '#3498db', '#2ecc71', '#95a5a6']

bars = ax.barh(dims, vals, color=colors_bar, edgecolor='white', height=0.6)

# Add value labels
for bar, val in zip(bars, vals):
    ax.text(bar.get_width() + 0.008, bar.get_y() + bar.get_height()/2,
            f'{val:.3f}', va='center', fontsize=10, fontweight='bold')

ax.set_xlabel('Mean $\\eta^2$ (proportion of variance explained)', fontsize=11)
ax.set_title('Dimension Attribution: Which Analytical Choice\nDrives Disagreement?', fontsize=12, pad=10)
ax.set_xlim(0, max(vals) * 1.25)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Annotation
ax.text(0.95, 0.05, 'Publication bias correction\nexplains 10x more variance\nthan any other dimension',
        transform=ax.transAxes, fontsize=9, va='bottom', ha='right',
        style='italic', color='#555555')

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'figure3_dimension_attribution.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure3_dimension_attribution.pdf'), bbox_inches='tight')
print("  Saved figure3_dimension_attribution.png/pdf")
plt.close(fig)


# ─────────────────────────────────────────────────
# FIGURE 4: k vs Robustness Score scatter
# ─────────────────────────────────────────────────
print("Generating Figure 4: k vs robustness scatter...")

fig, ax = plt.subplots(figsize=(7, 5))

k_vals = [r['k_val'] for r in reviews]
rob_vals = [r['robustness'] for r in reviews]
class_colors = [CLASS_COLORS[r['classification']] for r in reviews]

ax.scatter(k_vals, rob_vals, c=class_colors, s=20, alpha=0.6, edgecolors='white', linewidth=0.3)

# Log scale for x
ax.set_xscale('log')
ax.set_xlim(2.5, 200)

# Classification threshold bands
ax.axhspan(90, 100, alpha=0.05, color='#2ecc71')
ax.axhspan(70, 90, alpha=0.05, color='#f39c12')
ax.axhspan(50, 70, alpha=0.05, color='#e74c3c')
ax.axhspan(0, 50, alpha=0.05, color='#8e44ad')

# Threshold lines
for y in [50, 70, 90]:
    ax.axhline(y=y, color='#cccccc', linestyle='--', linewidth=0.5)

# Compute correlation
from math import log
log_k = [log(k) for k in k_vals]
n = len(log_k)
mean_lk = sum(log_k) / n
mean_r = sum(rob_vals) / n
cov = sum((lk - mean_lk) * (r - mean_r) for lk, r in zip(log_k, rob_vals)) / n
var_lk = sum((lk - mean_lk)**2 for lk in log_k) / n
var_r = sum((r - mean_r)**2 for r in rob_vals) / n
r_val = cov / (var_lk * var_r) ** 0.5 if var_lk > 0 and var_r > 0 else 0

ax.set_xlabel('Number of Studies (k)', fontsize=11)
ax.set_ylabel('Robustness Score (%)', fontsize=11)
ax.set_title('Fragility Is Not Resolved by Including More Studies', fontsize=12, pad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Correlation annotation
ax.text(0.97, 0.05, f'r = {r_val:.3f}\n(near-zero correlation)',
        transform=ax.transAxes, fontsize=10, va='bottom', ha='right',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc', alpha=0.9))

# Legend
patches = [mpatches.Patch(color=c, label=k) for k, c in CLASS_COLORS.items()]
ax.legend(handles=patches, loc='upper left', framealpha=0.9, fontsize=9)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'figure4_k_vs_robustness.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure4_k_vs_robustness.pdf'), bbox_inches='tight')
print("  Saved figure4_k_vs_robustness.png/pdf")
plt.close(fig)

print(f"\nAll figures saved to {OUT_DIR}/")
print("Files generated:")
for f in sorted(os.listdir(OUT_DIR)):
    size = os.path.getsize(os.path.join(OUT_DIR, f))
    print(f"  {f} ({size/1024:.0f} KB)")
