#!/usr/bin/env python3
"""
Volcano plot + Strip plot для всех SNV по каждому типу энергии.
Статистика: парный тест (Wilcoxon) по всем конструктам, FDR коррекция с заданным alpha.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

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
def fdr_bh(p_values: np.ndarray, alpha: float = 0.05) -> Tuple[np.ndarray, np.ndarray]:
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
    test_in = np.array([0.01, 0.02, 0.03])
    test_out = _apply_fdr_fast(test_in, alpha=0.05)
    if isinstance(test_out, tuple) and len(test_out) == 4:
        def apply_fdr(pvals, alpha=0.05):
            reject, qvals, _, _ = _apply_fdr_fast(pvals, alpha=alpha)
            return reject, qvals
    else:
        apply_fdr = _apply_fdr_fast
    logger.info("Используется apply_fdr_fast из stats_utils")
except ImportError:
    apply_fdr = fdr_bh
    logger.info("Используется встроенная FDR")

# ----------------------------------------------------------------------
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
    return list(data_dir.rglob("*_master_data.json")) if recursive else list(data_dir.glob("*_master_data.json"))

def collect_data_and_compute_pvalues(json_files: List[Path], test='wilcoxon', fdr_alpha=0.001) -> Dict[str, pd.DataFrame]:
    pairs = defaultdict(lambda: defaultdict(lambda: {'ref': [], 'alt': [], 'gene': None}))
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
                pairs[etype][pos]['ref'].append(float(ref))
                pairs[etype][pos]['alt'].append(float(alt))
                if pairs[etype][pos]['gene'] is None:
                    pairs[etype][pos]['gene'] = gene
            except:
                pass

    results = {}
    for etype, pos_dict in pairs.items():
        rows = []
        for pos, data in pos_dict.items():
            ref = np.array(data['ref'])
            alt = np.array(data['alt'])
            n_pairs = len(ref)
            if n_pairs < 3:
                continue
            try:
                if test == 'wilcoxon':
                    _, pval = wilcoxon(ref, alt)
                elif test == 'ttest':
                    _, pval = ttest_rel(ref, alt)
                else:
                    raise ValueError
            except Exception as e:
                logger.debug(f"Ошибка теста для {pos}: {e}")
                continue
            mean_diff = np.mean(alt - ref)
            median_diff = np.median(alt - ref)
            rows.append({
                'position': pos,
                'p_value': pval,
                'mean_diff': mean_diff,
                'median_diff': median_diff,
                'n_pairs': n_pairs,
                'gene': data['gene']
            })
        if not rows:
            continue
        df = pd.DataFrame(rows)
        reject, qvals = apply_fdr(df['p_value'].values, alpha=fdr_alpha)
        df['fdr_qvalue'] = qvals
        df['significant'] = reject
        df['neg_log10_p'] = -np.log10(df['p_value'])
        df['neg_log10_q'] = -np.log10(df['fdr_qvalue']).replace([np.inf, -np.inf], 0)
        results[etype] = df
    return results

def collect_all_diffs(json_files: List[Path]) -> Dict[str, pd.DataFrame]:
    diffs = defaultdict(lambda: {'positions': [], 'diffs': [], 'genes': []})
    for jf in json_files:
        try:
            with open(jf) as f:
                recs = json.load(f)
        except Exception:
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
                diff = float(alt) - float(ref)
                diffs[etype]['positions'].append(pos)
                diffs[etype]['diffs'].append(diff)
                diffs[etype]['genes'].append(gene)
            except:
                pass
    result = {}
    for etype, d in diffs.items():
        df = pd.DataFrame({'position': d['positions'], 'diff': d['diffs'], 'gene': d['genes']})
        result[etype] = df
    return result

def create_strip_plot(diff_df: pd.DataFrame, pval_df: pd.DataFrame, energy_type: str,
                      marker_size=3, opacity=0.3, fdr_alpha=0.001):
    # Убедимся, что в diff_df есть колонка 'position'
    if 'position' not in diff_df.columns:
        logger.error(f"В diff_df нет колонки 'position'. Столбцы: {diff_df.columns}")
        return px.scatter(title="Ошибка: нет колонки position")
    if 'position' not in pval_df.columns:
        logger.error(f"В pval_df нет колонки 'position'. Столбцы: {pval_df.columns}")
        return px.scatter(title="Ошибка: нет колонки position")
    merged = diff_df.merge(pval_df[['position', 'significant', 'mean_diff', 'neg_log10_p']], on='position', how='left')
    color_map = {True: 'red', False: 'gray'}
    fig = px.scatter(
        merged, x='position', y='diff', color='significant',
        color_discrete_map=color_map,
        hover_data={'position': True, 'diff': ':.3f', 'gene': True, 'neg_log10_p': ':.2f', 'mean_diff': ':.3f'},
        title=f"Strip plot: разности (alt - ref) – {energy_type}<br><sup>Красные точки: q < {fdr_alpha} (FDR).</sup>",
        labels={'position': 'Позиция SNV', 'diff': 'Δ (alt-ref), ккал/моль'},
        opacity=opacity
    )
    fig.update_traces(marker=dict(size=marker_size, line=dict(width=0)))
    for _, row in pval_df.iterrows():
        fig.add_annotation(
            x=row['position'], y=row['mean_diff'],
            text=f"μ={row['mean_diff']:.2f}",
            showarrow=False, font=dict(size=8, color='blue'), bgcolor='rgba(255,255,255,0.6)'
        )
    fig.update_layout(template='plotly_white', height=800, width=1200, xaxis_tickangle=-45)
    return fig

def create_volcano_plot(pval_df: pd.DataFrame, energy_type: str, color_map: Dict[int, str],
                        output_file: Path, fdr_alpha=0.001, effect_threshold=0.0,
                        marker_size=5, opacity=0.7, png=True):
    df = pval_df.copy()
    df['neg_log10_p'] = -np.log10(df['p_value'])
    df['neg_log10_q'] = -np.log10(df['fdr_qvalue']).replace([np.inf, -np.inf], 0)
    df['color'] = df['position'].map(color_map).fillna('gray')

    sig = df[df['significant']]
    non_sig = df[~df['significant']]

    fig = go.Figure()
    if not non_sig.empty:
        fig.add_trace(go.Scatter(
            x=non_sig['mean_diff'], y=non_sig['neg_log10_p'],
            mode='markers',
            marker=dict(color=non_sig['color'], size=marker_size, opacity=opacity, line=dict(width=0)),
            hovertext=non_sig.apply(lambda r: f"SNV {r['position']}<br>средняя разность={r['mean_diff']:.3f}<br>p-значение={r['p_value']:.2e}<br>скорректированное q={r['fdr_qvalue']:.2e}", axis=1),
            hoverinfo='text',
            name='не значимо'
        ))
    if not sig.empty:
        fig.add_trace(go.Scatter(
            x=sig['mean_diff'], y=sig['neg_log10_p'],
            mode='markers',
            marker=dict(color=sig['color'], size=marker_size, opacity=0.9, line=dict(width=1.5, color='black')),
            hovertext=sig.apply(lambda r: f"SNV {r['position']}<br>средняя разность={r['mean_diff']:.3f}<br>p-значение={r['p_value']:.2e}<br>скорректированное q={r['fdr_qvalue']:.2e}", axis=1),
            hoverinfo='text',
            name=f'значимо (q<{fdr_alpha})'
        ))
    fig.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='blue', annotation_text='p = 0,05')
    if df['significant'].any():
        min_q_log = df[df['significant']]['neg_log10_q'].min()
        fig.add_hline(y=min_q_log, line_dash='dot', line_color='red', annotation_text=f'FDR q = {fdr_alpha}')
    fig.add_vline(x=0, line_dash='dash', line_color='black', annotation_text='нет эффекта')
    if effect_threshold > 0:
        fig.add_vline(x=effect_threshold, line_dash='dash', line_color='green', annotation_text=f'+{effect_threshold}')
        fig.add_vline(x=-effect_threshold, line_dash='dash', line_color='green', annotation_text=f'-{effect_threshold}')

    fig.update_layout(
        title=f"График «вулкан» – {energy_type}<br><sup>Цвет – SNV. Чёрная обводка: q < {fdr_alpha} (FDR).</sup>",
        xaxis_title="Средняя разность (alt – ref), ккал/моль",
        yaxis_title="-log₁₀(p-значение)",
        template='plotly_white', height=700, width=1000,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
    )
    save_plot(fig, output_file, png)

def create_boxplot_by_gene(pval_df: pd.DataFrame, energy_type: str, output_file: Path, png=True):
    from scipy.stats import kruskal
    df_box = pval_df[['gene', 'mean_diff']].dropna()
    if df_box.empty:
        return
    genes = df_box['gene'].unique()
    groups = [df_box[df_box['gene'] == g]['mean_diff'].values for g in genes if len(df_box[df_box['gene'] == g]) >= 3]
    if len(groups) >= 2:
        h, p_kw = kruskal(*groups)
        title = f"Распределение средних разностей по генам – {energy_type}<br><sup>Тест Краскела–Уоллиса: χ² = {h:.2f}, p = {p_kw:.3e}</sup>"
    else:
        title = f"Распределение средних разностей по генам – {energy_type}"
    fig = px.box(df_box, x='gene', y='mean_diff', color='gene',
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 title=title, labels={'gene': 'Ген', 'mean_diff': 'Средняя разность (alt – ref), ккал/моль'},
                 points='outliers')
    fig.update_traces(marker=dict(size=4, line=dict(width=0)))
    fig.add_hline(y=0, line_dash="dash", line_color="black", annotation_text="Нет изменений")
    fig.update_layout(template='plotly_white', height=600, width=900,
                      legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1))
    save_plot(fig, output_file, png)

def save_plot(fig, path, png=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(fig, str(path))
    if png:
        try:
            pio.write_image(fig, str(path.with_suffix('.png')), scale=2)
        except Exception as e:
            logger.warning(f"Не удалось сохранить PNG: {e}")

def main():
    parser = argparse.ArgumentParser(description="Графики «вулкан» и стрип-плоты для SNV по типам энергии")
    parser.add_argument("data_dir", type=Path, help="Директория с master_data.json")
    parser.add_argument("-o", "--output_dir", type=Path, default=None)
    parser.add_argument("--test", choices=['wilcoxon','ttest'], default='wilcoxon')
    parser.add_argument("--fdr-alpha", type=float, default=0.001, help="Уровень значимости FDR")
    parser.add_argument("--effect-threshold", type=float, default=0.0)
    parser.add_argument("--marker-size", type=int, default=5)
    parser.add_argument("--opacity", type=float, default=0.7)
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument("--no-recurse", action="store_true")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Нет директории {data_dir}")
        sys.exit(1)

    out_dir = (args.output_dir or data_dir / "volcano_strip_plots").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    json_files = find_master_json_files(data_dir, recursive=not args.no_recurse)
    logger.info(f"Найдено {len(json_files)} файлов")
    if not json_files:
        sys.exit(1)

    logger.info("Вычисление p-значений и FDR...")
    pval_dfs = collect_data_and_compute_pvalues(json_files, test=args.test, fdr_alpha=args.fdr_alpha)
    logger.info("Сбор всех разностей для strip plot...")
    diff_dfs = collect_all_diffs(json_files)

    energy_types = set(pval_dfs.keys()) & set(diff_dfs.keys())
    if not energy_types:
        logger.error("Нет общих типов энергии")
        sys.exit(1)

    for etype in energy_types:
        pval_df = pval_dfs[etype]
        diff_df = diff_dfs[etype]
        logger.info(f"Обработка {etype}: {len(pval_df)} SNV, {len(diff_df)} точек")
        out_sub = out_dir / etype
        out_sub.mkdir(parents=True, exist_ok=True)

        pval_df.to_csv(out_sub / f"{etype}_statistics.csv", index=False)

        # Генерация цветовой карты для SNV
        snvs = sorted(pval_df['position'].unique())
        colors = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel + px.colors.qualitative.Plotly
        color_map = {pos: colors[i % len(colors)] for i, pos in enumerate(snvs)}

        # Strip plot
        fig_strip = create_strip_plot(diff_df, pval_df, etype,
                                      marker_size=args.marker_size,
                                      opacity=args.opacity,
                                      fdr_alpha=args.fdr_alpha)
        save_plot(fig_strip, out_sub / f"{etype}_strip_plot.html", png=not args.no_png)

        # Volcano plot
        volcano_file = out_sub / f"{etype}_volcano_plot.html"
        create_volcano_plot(pval_df, etype, color_map,
                            output_file=volcano_file,
                            fdr_alpha=args.fdr_alpha,
                            effect_threshold=args.effect_threshold,
                            marker_size=args.marker_size,
                            opacity=args.opacity,
                            png=not args.no_png)

        # Boxplot по генам
        if len(pval_df) > 1:
            create_boxplot_by_gene(pval_df, etype, out_sub / f"{etype}_boxplot_by_gene.html", png=not args.no_png)

        n_sig = pval_df['significant'].sum()
        logger.info(f"{etype}: значимых SNV после FDR (α={args.fdr_alpha}) = {n_sig} из {len(pval_df)}")

    logger.info(f"Готово. Результаты в {out_dir}")

if __name__ == "__main__":
    main()