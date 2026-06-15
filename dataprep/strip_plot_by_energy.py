#!/usr/bin/env python3
"""
Strip plot: для каждого SNV показать все разности (alt - ref) по конструктам.
Ось X – позиция SNV, ось Y – разность (alt - ref). Цвет – значимость после FDR.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import wilcoxon, ttest_rel

# ----------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
# FDR (Benjamini-Hochberg) - собственная реализация, возвращает 2 значения
def fdr_bh(p_values: np.ndarray, alpha: float = 0.05):
    pvals = np.asarray(p_values)
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1)
    qvals = pvals[order] * n / ranks
    qvals = np.minimum.accumulate(qvals[::-1])[::-1]
    qvals = qvals[np.argsort(order)]
    qvals = np.minimum(qvals, 1.0)
    reject = qvals < alpha
    return reject, qvals

try:
    from .stats_utils import apply_fdr_fast as _apply_fdr_fast
    test_input = np.array([0.01, 0.02, 0.03])
    test_out = _apply_fdr_fast(test_input, alpha=0.05)
    if isinstance(test_out, tuple) and len(test_out) == 4:
        def apply_fdr(pvals, alpha=0.05):
            reject, qvals, _, _ = _apply_fdr_fast(pvals, alpha=alpha)
            return reject, qvals
    else:
        apply_fdr = _apply_fdr_fast
except ImportError:
    apply_fdr = fdr_bh
    logger.info("Используется встроенная FDR")

# ----------------------------------------------------------------------
# Гены
MT_GENE_RANGES = {
    'ND1': (3307, 4262), 'ND2': (4470, 5511), 'CO1': (5904, 7445),
    'CO2': (7586, 8269), 'ATP8': (8366, 8572), 'ATP6': (8527, 9207),
    'CO3': (9207, 9990), 'ND3': (10059, 10404), 'ND4L': (10470, 10766),
    'ND4': (10760, 12141), 'ND5': (12337, 14148), 'ND6': (14149, 14673),
    'CYB': (14747, 15887),
}
def infer_mt_gene(pos):
    if pos is None:
        return 'unknown'
    for name, (s, e) in MT_GENE_RANGES.items():
        if s <= pos <= e:
            return name
    return 'non-coding'

# ----------------------------------------------------------------------
def find_master_json_files(data_dir: Path, recursive=True):
    pattern = "*_master_data.json"
    return list(data_dir.rglob(pattern)) if recursive else list(data_dir.glob(pattern))

def collect_diff_and_pairs(json_files: List[Path]) -> Dict[str, Dict]:
    """
    Возвращает для каждого типа энергии:
        - diffs: список всех разностей (alt - ref) для каждого конструкта (с позицией)
        - pairs: ref и alt списки по позициям для расчёта p-value
    """
    diffs_by_energy = defaultdict(lambda: {'positions': [], 'diffs': [], 'genes': []})
    pairs_by_energy = defaultdict(lambda: defaultdict(lambda: {'ref': [], 'alt': [], 'gene': None}))
    
    for jf in json_files:
        try:
            with open(jf) as f:
                recs = json.load(f)
        except Exception as e:
            logger.warning(f"Ошибка {jf}: {e}")
            continue
        if not recs:
            continue
        pos = recs[0].get('snp') or recs[0].get('position')
        if pos is None:
            continue
        pos = int(pos)
        gene = infer_mt_gene(pos)
        for rec in recs:
            etype = rec.get('energy_type')
            if not etype:
                continue
            ref = rec.get('ref_energy')
            alt = rec.get('alt_energy')
            if ref is None or alt is None:
                continue
            try:
                ref_f = float(ref)
                alt_f = float(alt)
                diff = alt_f - ref_f   # alt - ref, как в вашем исходном коде
                diffs_by_energy[etype]['positions'].append(pos)
                diffs_by_energy[etype]['diffs'].append(diff)
                diffs_by_energy[etype]['genes'].append(gene)
                
                pairs_by_energy[etype][pos]['ref'].append(ref_f)
                pairs_by_energy[etype][pos]['alt'].append(alt_f)
                if pairs_by_energy[etype][pos]['gene'] is None:
                    pairs_by_energy[etype][pos]['gene'] = gene
            except (ValueError, TypeError):
                pass
    # Преобразуем diffs в DataFrame
    diff_dfs = {}
    for etype, d in diffs_by_energy.items():
        diff_dfs[etype] = pd.DataFrame({
            'position': d['positions'],
            'diff': d['diffs'],
            'gene': d['genes']
        })
    return diff_dfs, pairs_by_energy

def compute_pvalues_from_pairs(pairs_by_energy: Dict, test='wilcoxon') -> Dict[str, pd.DataFrame]:
    """
    pairs_by_energy: energy_type -> { pos: {'ref': list, 'alt': list, 'gene': str} }
    Возвращает словарь: energy_type -> DataFrame с колонками position, p_value, fdr_qvalue, significant, mean_diff, median_diff, n_pairs, gene
    """
    results = {}
    for etype, pos_dict in pairs_by_energy.items():
        rows = []
        for pos, data in pos_dict.items():
            ref = np.array(data['ref'])
            alt = np.array(data['alt'])
            if len(ref) < 3:
                continue
            try:
                if test == 'wilcoxon':
                    _, pval = wilcoxon(ref, alt)
                elif test == 'ttest':
                    _, pval = ttest_rel(ref, alt)
                else:
                    raise ValueError("test must be 'wilcoxon' or 'ttest'")
            except Exception:
                continue
            mean_diff = np.mean(alt - ref)   # alt - ref
            median_diff = np.median(alt - ref)
            n_pairs = len(ref)
            rows.append({
                'position': pos,
                'p_value': pval,
                'mean_diff': mean_diff,
                'median_diff': median_diff,
                'n_pairs': n_pairs,
                'gene': data['gene']
            })
        if not rows:
            results[etype] = pd.DataFrame()
            continue
        df = pd.DataFrame(rows)
        reject, qvals = apply_fdr(df['p_value'].values, alpha=0.05)
        df['fdr_qvalue'] = qvals
        df['significant'] = reject
        df['neg_log10_p'] = -np.log10(df['p_value'])
        results[etype] = df
    return results

def create_strip_plot(
    diff_df: pd.DataFrame,
    pval_df: pd.DataFrame,
    energy_type: str,
    opacity: float = 0.3,
    point_size: int = 3
) -> go.Figure:
    """
    diff_df: колонки position, diff, gene
    pval_df: колонки position, significant, mean_diff, median_diff
    """
    merged = diff_df.merge(pval_df[['position', 'significant', 'mean_diff', 'median_diff']], on='position', how='left')
    color_map = {True: 'red', False: 'gray'}
    fig = px.scatter(
        merged,
        x='position',
        y='diff',
        color='significant',
        color_discrete_map=color_map,
        hover_data={'position': True, 'diff': ':.3f', 'gene': True},
        title=f"Strip plot: разности (alt - ref) для всех конструктов – {energy_type}",
        labels={'position': 'Позиция SNV', 'diff': 'Δ = alt - ref (ккал/моль)', 'significant': 'Значимость (FDR<0.05)'},
        opacity=opacity
    )
    # Добавим аннотации средних и медиан
    for _, row in pval_df.iterrows():
        pos = row['position']
        mean_val = row['mean_diff']
        median_val = row['median_diff']
        fig.add_annotation(
            x=pos, y=mean_val,
            text=f"μ={mean_val:.2f}",
            showarrow=False,
            font=dict(size=9, color='blue'),
            bgcolor='white',
            opacity=0.8
        )
        fig.add_annotation(
            x=pos, y=median_val,
            text=f"med={median_val:.2f}",
            showarrow=False,
            font=dict(size=8, color='green'),
            bgcolor='white',
            opacity=0.7,
            yshift=10
        )
    fig.update_traces(marker=dict(size=point_size, line=dict(width=0)))
    fig.update_layout(
        template='plotly_white',
        height=800,
        width=1400,
        xaxis=dict(tickangle=-45),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
    )
    return fig

def save_plot(fig, path, png=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(fig, str(path))
    if png:
        try:
            pio.write_image(fig, str(path.with_suffix('.png')), scale=2)
        except Exception:
            pass

def main():
    parser = argparse.ArgumentParser(description="Strip plot всех разностей (alt - ref) по SNV")
    parser.add_argument("data_dir", type=Path, help="Директория с master_data.json")
    parser.add_argument("-o", "--output_dir", type=Path, default=None)
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument("--no-recurse", action="store_true")
    parser.add_argument("--test", choices=['wilcoxon','ttest'], default='wilcoxon')
    parser.add_argument("--opacity", type=float, default=0.3, help="Прозрачность точек (0-1)")
    parser.add_argument("--point-size", type=int, default=3, help="Размер точек")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Нет директории {data_dir}")
        sys.exit(1)
    out_dir = (args.output_dir or data_dir / "strip_plots").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    json_files = find_master_json_files(data_dir, recursive=not args.no_recurse)
    logger.info(f"Найдено {len(json_files)} файлов")
    if not json_files:
        sys.exit(1)

    logger.info("Сбор разностей (alt - ref) и пар ref/alt...")
    diff_dfs, pairs_by_energy = collect_diff_and_pairs(json_files)
    logger.info(f"Типы энергии: {list(diff_dfs.keys())}")

    logger.info("Вычисление p-value и FDR...")
    pval_dfs = compute_pvalues_from_pairs(pairs_by_energy, test=args.test)

    for etype in diff_dfs.keys():
        if etype not in pval_dfs or pval_dfs[etype].empty:
            logger.warning(f"Для {etype} нет p-value, пропускаем")
            continue
        diff_df = diff_dfs[etype]
        pval_df = pval_dfs[etype]
        logger.info(f"Создание strip plot для {etype}: {len(diff_df)} точек, {len(pval_df)} SNV")
        fig = create_strip_plot(diff_df, pval_df, etype, opacity=args.opacity, point_size=args.point_size)
        out_sub = out_dir / etype
        out_sub.mkdir(exist_ok=True)
        save_plot(fig, out_sub / f"{etype}_strip_plot.html", png=not args.no_png)
        pval_df.to_csv(out_sub / f"{etype}_pvalues_fdr.csv", index=False)
        logger.info(f"Сохранено для {etype}, значимых: {pval_df['significant'].sum()} из {len(pval_df)}")

    logger.info(f"Готово. Результаты в {out_dir}")

if __name__ == "__main__":
    main()