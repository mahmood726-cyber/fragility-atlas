"""Load Pairwise70 RDA files and select primary analyses for the Fragility Atlas."""

import math
import pyreadr
import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ReviewData:
    """Processed data for a single Cochrane review's primary analysis."""
    review_id: str
    review_doi: str
    analysis_name: str
    k: int                        # number of studies
    yi: np.ndarray                # effect sizes (log scale for ratio, raw for diff)
    sei: np.ndarray               # standard errors
    ni: np.ndarray                # total sample sizes per study
    study_labels: list            # study identifiers
    scale: str                    # 'ratio' or 'difference'
    cochrane_pooled: float        # Cochrane's published pooled estimate (original scale)
    cochrane_ci_lo: float
    cochrane_ci_hi: float
    is_significant: bool          # Cochrane conclusion: CI excludes null


def load_review(rda_path: str) -> Optional[ReviewData]:
    """Load a single RDA file, select primary analysis, return ReviewData or None."""
    path = Path(rda_path)
    review_id = path.stem.split('_')[0]  # e.g., 'CD000028'

    result = pyreadr.read_r(str(path))
    df = list(result.values())[0].copy()

    # Standardize column names (handle both dot and space separators)
    df.columns = df.columns.str.replace(' ', '.', regex=False)

    # Get review metadata
    review_doi = str(df['review_doi'].iloc[0]) if 'review_doi' in df.columns else ''

    # Group analyses and find the primary one
    primary = _select_primary_analysis(df)
    if primary is None or len(primary) < 3:
        return None

    analysis_name = str(primary['Analysis.name'].iloc[0])

    # P0-2 FIX: Determine scale from raw data columns, not sign of means.
    # Binary outcomes (cases/N present and nonzero) → ratio scale (RR/OR).
    # Continuous outcomes (mean/SD present) → difference scale (MD/SMD).
    # Fallback: if all Mean > 0 → ratio; else → difference.
    has_binary = (primary['Experimental.cases'].notna() & (primary['Experimental.cases'] > 0)).any()
    exp_mean = pd.to_numeric(primary['Experimental.mean'], errors='coerce')
    exp_sd = pd.to_numeric(primary['Experimental.SD'], errors='coerce')
    has_continuous = (exp_mean.notna() & (exp_mean != 0)).any() and (exp_sd.notna() & (exp_sd != 0)).any()

    if has_continuous and not has_binary:
        scale = 'difference'
    elif has_binary:
        scale = 'ratio'
    else:
        # Fallback: infer from whether any Mean values are negative
        means = primary['Mean'].copy().dropna()
        scale = 'ratio' if (means > 0).all() else 'difference'

    # Compute yi and sei
    yi, sei, ni, labels = _compute_effects(primary, scale)

    if yi is None or len(yi) < 3:
        return None

    # Filter out studies with invalid/zero SE
    valid = (sei > 0) & np.isfinite(yi) & np.isfinite(sei)
    yi = yi[valid].copy()
    sei = sei[valid].copy()
    ni = ni[valid].copy()
    labels = [l for l, v in zip(labels, valid) if v]

    if len(yi) < 3:
        return None

    # P0-1 FIX: Compute actual FE inverse-variance pooled estimate from yi/sei
    wi = 1.0 / (sei ** 2)
    theta_fe = float(np.sum(wi * yi) / np.sum(wi))
    se_fe = float(1.0 / math.sqrt(np.sum(wi)))
    cochrane_pooled = theta_fe
    cochrane_ci_lo = theta_fe - 1.96 * se_fe
    cochrane_ci_hi = theta_fe + 1.96 * se_fe

    # Determine significance: CI excludes zero (log-scale for ratio, raw for difference)
    is_sig = (cochrane_ci_lo > 0) or (cochrane_ci_hi < 0)

    return ReviewData(
        review_id=review_id,
        review_doi=review_doi,
        analysis_name=analysis_name,
        k=len(yi),
        yi=yi,
        sei=sei,
        ni=ni,
        study_labels=labels,
        scale=scale,
        cochrane_pooled=cochrane_pooled,
        cochrane_ci_lo=cochrane_ci_lo,
        cochrane_ci_hi=cochrane_ci_hi,
        is_significant=is_sig,
    )


def _select_primary_analysis(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Select the primary analysis: largest k among binary outcomes, then largest k overall."""
    groups = []
    for (grp, num), sub in df.groupby(['Analysis.group', 'Analysis.number']):
        k = len(sub)
        has_binary = (sub['Experimental.cases'].notna() & (sub['Experimental.cases'] > 0)).any()
        groups.append({
            'grp': grp, 'num': num, 'k': k, 'binary': has_binary,
            'name': sub['Analysis.name'].iloc[0]
        })

    if not groups:
        return None

    groups_df = pd.DataFrame(groups)

    # Prefer binary outcomes with largest k
    binary = groups_df[groups_df['binary']].copy()
    if len(binary) > 0:
        best = binary.loc[binary['k'].idxmax()]
    else:
        best = groups_df.loc[groups_df['k'].idxmax()]

    grp, num = int(best['grp']), int(best['num'])
    return df[(df['Analysis.group'] == grp) & (df['Analysis.number'] == num)]


def _compute_effects(primary: pd.DataFrame, scale: str):
    """Compute yi (effect size) and sei (standard error) from raw data or back-calculation."""
    yi_list, sei_list, ni_list, labels = [], [], [], []

    for _, row in primary.iterrows():
        study = str(row['Study'])
        # iterrows yields a Series; cell access returns scalar floats, not
        # arrays/Series — so .copy() would raise AttributeError. Scalars are
        # immutable, so direct assignment is safe here.
        mean_val = row['Mean']
        ci_lo = row['CI.start']
        ci_hi = row['CI.end']

        # Total sample size
        n_exp = int(row['Experimental.N']) if pd.notna(row['Experimental.N']) else 0
        n_ctrl = int(row['Control.N']) if pd.notna(row['Control.N']) else 0
        n_total = n_exp + n_ctrl

        if pd.isna(mean_val) or pd.isna(ci_lo) or pd.isna(ci_hi):
            continue

        if scale == 'ratio':
            # Cochrane stores natural-scale RR/OR — log-transform
            if mean_val <= 0 or ci_lo <= 0 or ci_hi <= 0:
                continue
            y = math.log(mean_val)
            se = (math.log(ci_hi) - math.log(ci_lo)) / (2 * 1.96)
        else:
            # Difference scale (MD, SMD, RD) — already on additive scale
            y = mean_val
            se = (ci_hi - ci_lo) / (2 * 1.96)

        if se <= 0 or not math.isfinite(se) or not math.isfinite(y):
            continue

        yi_list.append(y)
        sei_list.append(se)
        ni_list.append(max(n_total, 1))
        labels.append(study)

    if not yi_list:
        return None, None, None, None

    return np.array(yi_list), np.array(sei_list), np.array(ni_list), labels


def load_all_reviews(pairwise_dir: str, min_k: int = 3):
    """Load all Pairwise70 reviews, yielding ReviewData for eligible ones."""
    pairwise_path = Path(pairwise_dir)
    rda_files = sorted(pairwise_path.glob('*.rda'))

    for rda in rda_files:
        try:
            review = load_review(str(rda))
            if review is not None and review.k >= min_k:
                yield review
        except Exception as e:
            print(f"Warning: failed to load {rda.name}: {e}")
            continue
