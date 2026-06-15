#!/usr/bin/env python3
"""
Манхэттен‑подобный график с категориальной осью X.
Ось Y: –log₁₀(скорректированное p-значение, FDR).
Чёрная обводка: q < α И значение в крайнем 0.1% хвосте своего SNV.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

try:
    from .stats_utils import calculate_z_score_p_value, apply_fdr_fast
except ImportError:
    try:
        from stats_utils import calculate_z_score_p_value, apply_fdr_fast
    except ImportError:
        from scipy.stats import norm
        def calculate_z_score_p_value(diffs):
            median = np.median(diffs)
            mad = np.median(np.abs(diffs - median))
            robust_std = mad / 0.6745 if mad != 0 else 1e-10
            z = (diffs - median) / robust_std
            p = 2 * (1 - norm.cdf(np.abs(z)))
            return z, p
        def apply_fdr_fast(pvals, alpha=0.05):
            from statsmodels.stats.multitest import multipletests
            reject, qvals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
            return reject, qvals

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
def find_master_json_files(data_dir: Path, recursive=True) -> List[Path]:
    return list(data_dir.rglob("*_master_data.json")) if recursive else list(data_dir.glob("*_master_data.json"))

def collect_all_points(json_files: List[Path], energy_type: str) -> pd.DataFrame:
    records = []
    for jf in json_files:
        try:
            with open(jf) as f:
                data = json.load(f)
        except Exception:
            continue
        if not data:
            continue
        pos = data[0].get('snp') or data[0].get('position')
        if pos is None:
            continue
        pos = int(pos)
        for rec in data:
            if rec.get('energy_type') != energy_type:
                continue
            ref = rec.get('ref_energy')
            alt = rec.get('alt_energy')
            if ref is None or alt is None:
                continue
            diff = float(ref) - float(alt)
            records.append({'position': pos, 'diff': diff})
    return pd.DataFrame(records)

def compute_pvalues(df: pd.DataFrame) -> pd.DataFrame:
    df['p_value'] = np.nan
    for pos, group in df.groupby('position'):
        diffs = group['diff'].values
        if len(diffs) < 3:
            continue
        _, p_vals = calculate_z_score_p_value(diffs)
        df.loc[group.index, 'p_value'] = p_vals
    return df

def add_fdr_and_extremes(df: pd.DataFrame, fdr_alpha: float = 0.01, extreme_percentile: float = 0.1) -> pd.DataFrame:
    valid = df.dropna(subset=['p_value'])
    if valid.empty:
        df['q_value'] = np.nan
        df['significant_fdr'] = False
        df['extreme'] = False
        return df
    res = apply_fdr_fast(valid['p_value'].values, alpha=fdr_alpha)
    if len(res) == 2:
        reject, qvals = res
    else:
        reject, qvals, _, _ = res
    valid['q_value'] = qvals
    valid['significant_fdr'] = reject
    df = df.merge(valid[['q_value', 'significant_fdr']], left_index=True, right_index=True, how='left')
    df['q_value'] = df['q_value'].fillna(1.0)
    df['significant_fdr'] = df['significant_fdr'].fillna(False)

    df['extreme'] = False
    low = extreme_percentile
    high = 100 - extreme_percentile
    for pos, group in df.groupby('position'):
        diffs = group['diff'].values
        if len(diffs) < 3:
            continue
        p_low = np.percentile(diffs, low)
        p_high = np.percentile(diffs, high)
        idx = group.index
        mask = (group['diff'] <= p_low) | (group['diff'] >= p_high)
        df.loc[idx[mask], 'extreme'] = True

    df['significant'] = df['significant_fdr'] & df['extreme']
    return df

def create_manhattan_plot(df: pd.DataFrame, energy_type: str, output_file: Path,
                          fdr_alpha: float = 0.01, png=True):
    df = add_fdr_and_extremes(df, fdr_alpha=fdr_alpha)
    df = df.dropna(subset=['q_value'])
    if df.empty:
        return
    df['neg_log10_q'] = -np.log10(df['q_value'])
    positions = sorted(df['position'].unique())
    df['pos_str'] = df['position'].astype(str)
    fig = px.scatter(
        df, x='pos_str', y='neg_log10_q',
        color_discrete_sequence=['gray'],
        category_orders={'pos_str': [str(p) for p in positions]},
        title=f"Манхэттен-график – {energy_type}<br>"
              f"<sup>Ось Y: –log₁₀(скорректированное p-значение, FDR). Чёрная обводка: q &lt; {fdr_alpha} и значение в крайнем 0.1% хвосте SNV.</sup>",
        labels={'pos_str': 'Позиция SNV', 'neg_log10_q': '–log₁₀(q-значение)'},
        hover_data={'position': True, 'diff': ':.3f', 'p_value': ':.2e', 'q_value': ':.2e', 'extreme': True}
    )
    fig.update_traces(marker=dict(size=6, opacity=0.7, line=dict(width=0)))
    sig = df[df['significant']]
    if not sig.empty:
        sig['pos_str'] = sig['position'].astype(str)
        fig.add_trace(go.Scatter(
            x=sig['pos_str'], y=sig['neg_log10_q'],
            mode='markers',
            marker=dict(color='gray', size=6, line=dict(width=1.5, color='black'), opacity=0.9),
            hoverinfo='skip',
            name=f'Значимо (q &lt; {fdr_alpha})'
        ))
    fig.add_hline(y=-np.log10(fdr_alpha), line_dash="dash", line_color="black",
                  annotation_text=f"q = {fdr_alpha}", annotation_position="top right")
    fig.update_layout(template='plotly_white', height=800, width=1200,
                      xaxis=dict(tickangle=-45, tickmode='array', tickvals=[str(p) for p in positions],
                                 ticktext=[str(p) for p in positions]),
                      legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1))
    save_plot(fig, output_file, png)

def save_plot(fig, path, png=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(fig, str(path))
    if png:
        try:
            pio.write_image(fig, str(path.with_suffix('.png')), scale=2)
        except Exception:
            pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("-o", "--output_dir", type=Path, default=None)
    parser.add_argument("--energy-type", required=True, choices=["Energy","EnergyLeft","EnergyRight"])
    parser.add_argument("--fdr-alpha", type=float, default=0.01, help="Уровень значимости FDR (например, 0.01)")
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument("--no-recurse", action="store_true")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Нет директории {data_dir}")
        sys.exit(1)

    out_dir = (args.output_dir or data_dir / "manhattan_extreme").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    json_files = find_master_json_files(data_dir, recursive=not args.no_recurse)
    logger.info(f"Найдено {len(json_files)} файлов")
    if not json_files:
        sys.exit(1)

    logger.info(f"Обработка {args.energy_type}...")
    df = collect_all_points(json_files, args.energy_type)
    if df.empty:
        logger.error("Нет данных")
        sys.exit(1)
    logger.info(f"Всего точек: {len(df)}")
    df = compute_pvalues(df)
    create_manhattan_plot(df, args.energy_type, out_dir / f"{args.energy_type}_manhattan.html",
                          fdr_alpha=args.fdr_alpha, png=not args.no_png)
    df.to_csv(out_dir / f"{args.energy_type}_data.csv", index=False)
    logger.info(f"Готово. Результаты в {out_dir}")

if __name__ == "__main__":
    main()