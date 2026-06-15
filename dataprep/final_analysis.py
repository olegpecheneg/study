#!/usr/bin/env python3
"""
FINAL WORKING SCRIPT – все графики, исправленная версия.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from collections import defaultdict
from typing import List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import mannwhitneyu, kruskal, ttest_1samp

# ----------------------------------------------------------------------
# Импорт утилит с fallback-реализациями
try:
    from stats_utils import calculate_z_score_p_value, apply_fdr_fast
except ImportError:
    from scipy.stats import norm
    # fallback реализация Z-score
    def calculate_z_score_p_value(diffs):
        median = np.median(diffs)
        mad = np.median(np.abs(diffs - median))
        if mad == 0:
            z = np.zeros_like(diffs)
            p = np.ones_like(diffs)
        else:
            robust_std = mad / 0.6745
            z = (diffs - median) / robust_std
            p = 2 * (1 - norm.cdf(np.abs(z)))
        return z, p

    # fallback FDR – обязательно импортируем multipletests
    def apply_fdr_fast(pvals, alpha=0.05):
        from statsmodels.stats.multitest import multipletests
        reject, qvals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
        return reject, qvals

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

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

def find_master_json_files(data_dir: Path, recursive=True) -> List[Path]:
    return list(data_dir.rglob("*_master_data.json")) if recursive else list(data_dir.glob("*_master_data.json"))

def collect_all_points(json_files: List[Path], energy_type: str) -> pd.DataFrame:
    """Собирает все точки (разности ref-alt) для заданного типа энергии."""
    records = []
    for jf in json_files:
        try:
            with open(jf, encoding='utf-8') as f:
                data = json.load(f)
        except Exception as e:
            logger.warning(f"Не удалось прочитать {jf}: {e}")
            continue
        if not data:
            logger.debug(f"Пустой файл: {jf}")
            continue

        # Определяем позицию SNV из первой записи
        first = data[0]
        pos = first.get('snp') or first.get('position')
        if pos is None:
            logger.warning(f"В файле {jf.name} нет поля 'snp' или 'position', пропускаем")
            continue
        try:
            pos = int(pos)
        except (ValueError, TypeError):
            logger.warning(f"Некорректная позиция {pos} в {jf.name}, пропускаем")
            continue

        gene = infer_mt_gene(pos)

        for idx, rec in enumerate(data):
            if rec.get('energy_type') != energy_type:
                continue
            ref = rec.get('ref_energy')
            alt = rec.get('alt_energy')
            if ref is None or alt is None:
                continue
            try:
                diff = float(ref) - float(alt)
            except (ValueError, TypeError):
                continue
            records.append({
                'position': pos,
                'construct_index': idx,
                'diff': diff,
                'gene': gene,
            })

    if not records:
        logger.warning(f"Не найдено данных для типа энергии {energy_type}")
        return pd.DataFrame()
    df = pd.DataFrame(records)
    logger.info(f"Собрано {len(df)} точек для {energy_type} из {len(json_files)} файлов")
    return df

def compute_pvalues(df: pd.DataFrame) -> pd.DataFrame:
    """Вычисляет p-value для каждого SNV на основе разностей."""
    df['p_value'] = np.nan
    for pos, group in df.groupby('position'):
        diffs = group['diff'].values
        if len(diffs) < 3:
            continue
        _, p_vals = calculate_z_score_p_value(diffs)
        df.loc[group.index, 'p_value'] = p_vals
    return df

def add_fdr_and_extremes(df: pd.DataFrame, extreme_percentile: float = 0.1) -> pd.DataFrame:
    """Добавляет FDR q-value и флаг экстремальности (0.1% хвост)."""
    valid = df.dropna(subset=['p_value'])
    if valid.empty:
        df['q_value'] = np.nan
        df['significant_fdr'] = False
        df['extreme'] = False
        df['significant'] = False
        return df

    res = apply_fdr_fast(valid['p_value'].values, alpha=0.05)
    if len(res) == 4:
        reject, qvals = res[0], res[1]
    else:
        reject, qvals = res
    valid = valid.copy()
    valid['q_value'] = qvals
    valid['significant_fdr'] = reject

    # Объединяем обратно
    df = df.merge(valid[['q_value', 'significant_fdr']], left_index=True, right_index=True, how='left')
    df['q_value'] = df['q_value'].fillna(1.0)
    df['significant_fdr'] = df['significant_fdr'].fillna(False)

    # Экстремальные хвосты (указанный процентиль снизу и сверху)
    df['extreme'] = False
    for pos, group in df.groupby('position'):
        diffs = group['diff'].values
        if len(diffs) < 3:
            continue
        low = np.percentile(diffs, extreme_percentile)
        high = np.percentile(diffs, 100 - extreme_percentile)
        idx = group.index
        mask = (group['diff'] <= low) | (group['diff'] >= high)
        df.loc[idx[mask], 'extreme'] = True

    df['significant'] = df['significant_fdr'] & df['extreme']
    return df


def save_plot(fig, path, png=True):
    """Сохраняет график в HTML и, опционально, в PNG."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(fig, str(path))
    if png:
        try:
            # требуется kaleido или orca
            pio.write_image(fig, str(path.with_suffix('.png')), scale=2)
        except Exception as e:
            logger.warning(f"Не удалось сохранить PNG для {path.name}: {e}")

def plot_manhattan(df: pd.DataFrame, energy_type: str, extreme_percentile: float, out_file: Path, png=True):
    """Манхэттен-подобный график."""
    df = add_fdr_and_extremes(df, extreme_percentile)
    df = df.dropna(subset=['q_value'])
    if df.empty:
        logger.warning(f"Нет данных для Manhattan plot ({energy_type})")
        return

    df['neg_log10_q'] = -np.log10(df['q_value'])
    positions = sorted(df['position'].unique())
    df['pos_str'] = df['position'].astype(str)

    fig = px.scatter(
        df, x='pos_str', y='neg_log10_q', color='gene',
        category_orders={'pos_str': [str(p) for p in positions]},
        color_discrete_sequence=px.colors.qualitative.Set2,
        title=f"Manhattan plot – {energy_type}<br>"
              f"<sup>Y: -log10(q-value) (FDR). Чёрная обводка: q<0.05 & крайний хвост ({extreme_percentile}%).</sup>",
        labels={'pos_str': 'SNV position', 'neg_log10_q': '-log10(q-value)', 'gene': 'Gene'},
        hover_data={'position': True, 'diff': ':.3f', 'p_value': ':.2e', 'q_value': ':.2e'}
    )
    fig.update_traces(marker=dict(size=5, opacity=0.7, line=dict(width=0)))

    sig = df[df['significant']]
    if not sig.empty:
        gene_colors = {g: px.colors.qualitative.Set2[i % len(px.colors.qualitative.Set2)]
                       for i, g in enumerate(sig['gene'].unique())}
        fig.add_trace(go.Scatter(
            x=sig['pos_str'], y=sig['neg_log10_q'],
            mode='markers',
            marker=dict(color=[gene_colors[g] for g in sig['gene']],
                        size=5, line=dict(width=1.5, color='black'), opacity=0.9),
            hoverinfo='skip',
            name='значимые & экстремальные'
        ))

    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", annotation_text="q = 0.05")
    fig.update_layout(template='plotly_white', height=800, width=1200,
                      xaxis=dict(tickangle=-45, tickmode='array',
                                 tickvals=[str(p) for p in positions],
                                 ticktext=[str(p) for p in positions]),
                      legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1))
    save_plot(fig, out_file, png)

def plot_transition_heatmap(json_files: List[Path], energy_type: str, out_file: Path, png=True):
    """Тепловая карта по нуклеотидным заменам."""
    trans_data = defaultdict(list)
    for jf in json_files:
        try:
            with open(jf, encoding='utf-8') as f:
                data = json.load(f)
        except Exception:
            continue
        if not data:
            continue
        ref = data[0].get('ref_allele')
        alt = data[0].get('alt_allele')
        if not ref or not alt:
            continue
        trans = f"{ref}→{alt}"
        for rec in data:
            if rec.get('energy_type') != energy_type:
                continue
            r = rec.get('ref_energy')
            a = rec.get('alt_energy')
            if r is None or a is None:
                continue
            try:
                diff = float(r) - float(a)
                trans_data[trans].append(diff)
            except (ValueError, TypeError):
                pass

    if not trans_data:
        logger.warning(f"Нет данных для transition heatmap ({energy_type})")
        return

    transitions = sorted(trans_data.keys())
    means, pvals = [], []
    for t in transitions:
        vals = trans_data[t]
        if len(vals) < 3:
            means.append(np.nan)
            pvals.append(1.0)
        else:
            _, p = ttest_1samp(vals, 0)
            means.append(np.mean(vals))
            pvals.append(p)

    stars = ['***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '' for p in pvals]
    cell_text = [f"{means[i]:.2f}{stars[i]}" for i in range(len(means))]

    fig = go.Figure(data=go.Heatmap(
        z=[means], x=[energy_type], y=transitions,
        text=[cell_text], texttemplate="%{text}", textfont={"size": 12},
        colorscale='RdBu', zmid=0,
        colorbar=dict(title='Mean Δ (ref-alt), kcal/mol')
    ))
    fig.update_layout(title=f"Transition heatmap – {energy_type}<br><sup>Stars: t-test (H0: mean=0)</sup>",
                      xaxis_title="Energy type", yaxis_title="Nucleotide change",
                      height=600, width=500)
    save_plot(fig, out_file, png)

def plot_boxplot_by_gene(df: pd.DataFrame, energy_type: str, out_file: Path, png=True):
    """Boxplot средних разностей по генам."""
    snv_stats = df.groupby(['position', 'gene'])['diff'].mean().reset_index()
    if snv_stats.empty:
        logger.warning(f"Нет данных для boxplot по генам ({energy_type})")
        return

    groups = [snv_stats[snv_stats['gene'] == g]['diff'].values for g in snv_stats['gene'].unique()
              if len(snv_stats[snv_stats['gene'] == g]) >= 3]
    if len(groups) >= 2:
        h, p_kw = kruskal(*groups)
        title = f"Mean diff per gene – {energy_type}<br><sup>Kruskal-Wallis p = {p_kw:.3e}</sup>"
    else:
        title = f"Mean diff per gene – {energy_type}"

    fig = px.box(snv_stats, x='gene', y='diff', color='gene',
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 title=title, labels={'gene': 'Gene', 'diff': 'Mean (ref-alt), kcal/mol'},
                 points='outliers')
    fig.update_traces(marker=dict(size=4, line=dict(width=0)))
    fig.add_hline(y=0, line_dash="dash", line_color="black", annotation_text="No change")
    fig.update_layout(template='plotly_white', height=600, width=900,
                      legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1))
    save_plot(fig, out_file, png)

def plot_bar_extreme(df: pd.DataFrame, energy_type: str, out_dir: Path, top_n=5, png=True):
    """Столбчатые диаграммы для SNV с наиболее отрицательным и положительным эффектом."""
    snv_stats = df.groupby('position')['diff'].agg(mean='mean', std='std', count='size').reset_index()
    if snv_stats.empty:
        return
    snv_stats = snv_stats.sort_values('mean')
    neg = snv_stats.head(top_n)['position'].tolist()
    pos = snv_stats.tail(top_n)['position'].tolist()

    def _bar(positions, sign, out_file):
        sub = snv_stats[snv_stats['position'].isin(positions)].copy()
        if sub.empty:
            return
        sub = sub.sort_values('mean')
        sub['label'] = sub['position'].astype(str)

        all_other = df[~df['position'].isin(positions)]['diff'].values
        p_vals = []
        for _, row in sub.iterrows():
            vals = df[df['position'] == row['position']]['diff'].values
            if len(vals) < 3 or len(all_other) < 3:
                p_vals.append(1.0)
            else:
                _, p = mannwhitneyu(vals, all_other, alternative='two-sided')
                p_vals.append(p)
        # Bonferroni correction
        p_corr = [min(p * len(positions), 1.0) for p in p_vals]
        sub['p_corr'] = p_corr
        sub['star'] = sub['p_corr'].apply(lambda p: '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '')

        fig = px.bar(sub, x='label', y='mean', error_y='std', color='mean',
                     color_continuous_scale='RdBu',
                     title=f"Top {top_n} {sign} changes – {energy_type}<br>"
                           f"<sup>Stars: Mann-Whitney U vs all other SNV, Bonferroni corrected.</sup>",
                     labels={'label': 'SNV', 'mean': 'Mean diff (ref-alt), kcal/mol'},
                     hover_data={'position': True, 'count': True, 'p_corr': ':.3e'})
        for _, row in sub.iterrows():
            y_anno = row['mean'] + row['std'] + 0.1
            fig.add_annotation(x=row['label'], y=y_anno, text=row['star'],
                               showarrow=False, font=dict(size=12, color='red'))
        fig.add_hline(y=0, line_dash="dash", line_color="black", annotation_text="No change")
        fig.update_traces(marker_line_width=0, opacity=0.8)
        fig.update_layout(template='plotly_white', height=600, width=800)
        save_plot(fig, out_file, png)

    _bar(neg, "negative", out_dir / f"{energy_type}_bar_negative.html")
    _bar(pos, "positive", out_dir / f"{energy_type}_bar_positive.html")

def plot_strip_by_construct(df: pd.DataFrame, energy_type: str, out_file: Path, png=True):
    """Strip plot: разности для каждого конструкта, сгруппированные по SNV."""
    if 'construct_index' not in df.columns:
        df['construct_index'] = df.groupby('position').cumcount()
    max_len = df.groupby('position')['construct_index'].max().max() + 1 if not df.empty else 0
    logger.info(f"Strip plot для {energy_type}: max_len = {max_len}")

    snvs = sorted(df['position'].unique())
    colors = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel
    color_map = {pos: colors[i % len(colors)] for i, pos in enumerate(snvs)}

    fig = go.Figure()
    for pos, group in df.groupby('position'):
        fig.add_trace(go.Scatter(
            x=group['construct_index'], y=group['diff'],
            mode='markers',
            name=f"SNV {pos}",
            marker=dict(color=color_map[pos], size=3, opacity=0.5, line=dict(width=0)),
            hoverinfo='skip'
        ))
    fig.update_layout(
        title=f"Strip plot by construct – {energy_type}<br>"
              f"<sup>X: construct index (0..{max_len-1}). Y: ref - alt (kcal/mol). Color: SNV.</sup>",
        xaxis_title="Construct index",
        yaxis_title="ref - alt (kcal/mol)",
        template='plotly_white',
        height=800, width=1200,
        legend=dict(orientation='v', yanchor='top', y=1, xanchor='left', x=1.02)
    )
    fig.add_hline(y=0, line_dash="dash", line_color="black", annotation_text="No change")
    save_plot(fig, out_file, png)

def process_energy_type(json_files: List[Path], energy_type: str, out_dir: Path, extreme_percentile: float, png: bool):
    """Запускает все графики для одного типа энергии."""
    logger.info(f"Обработка {energy_type}...")
    df = collect_all_points(json_files, energy_type)
    if df.empty:
        logger.warning(f"Нет данных для {energy_type}, пропускаем")
        return
    logger.info(f"Всего точек: {len(df)}")
    df = compute_pvalues(df)

    # Сохраняем все данные в CSV
    df.to_csv(out_dir / f"{energy_type}_all_data.csv", index=False)

    # Генерируем графики
    plot_manhattan(df, energy_type, extreme_percentile, out_dir / f"{energy_type}_manhattan.html", png)
    plot_transition_heatmap(json_files, energy_type, out_dir / f"{energy_type}_transition_heatmap.html", png)
    plot_boxplot_by_gene(df, energy_type, out_dir / f"{energy_type}_boxplot_by_gene.html", png)
    plot_bar_extreme(df, energy_type, out_dir, top_n=5, png=png)
    plot_strip_by_construct(df, energy_type, out_dir / f"{energy_type}_strip_by_construct.html", png)

def main():
    parser = argparse.ArgumentParser(description="Final working analysis script")
    parser.add_argument("data_dir", type=Path, help="Директория с master_data.json (ищет рекурсивно)")
    parser.add_argument("-o", "--output_dir", type=Path, default=None,
                        help="Папка для сохранения результатов (по умолчанию data_dir/final_analysis)")
    parser.add_argument("--energy-type", default=None, choices=["Energy","EnergyLeft","EnergyRight"],
                        help="Тип энергии для анализа")
    parser.add_argument("--all-energy-types", action="store_true",
                        help="Обработать все три типа энергии")
    parser.add_argument("--extreme-percentile", type=float, default=0.1,
                        help="Процент для хвоста (по умолчанию 0.1%%)")
    parser.add_argument("--no-png", action="store_true", help="Не сохранять PNG")
    parser.add_argument("--no-recurse", action="store_true", help="Не искать рекурсивно в подпапках")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Директория не найдена: {data_dir}")
        sys.exit(1)

    out_dir = (args.output_dir or data_dir / "final_analysis").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    json_files = find_master_json_files(data_dir, recursive=not args.no_recurse)
    logger.info(f"Найдено {len(json_files)} master_data.json файлов")
    if not json_files:
        logger.error("Нет файлов для анализа")
        sys.exit(1)

    png = not args.no_png

    if args.energy_type:
        e_out = out_dir / args.energy_type
        e_out.mkdir(exist_ok=True)
        process_energy_type(json_files, args.energy_type, e_out, args.extreme_percentile, png)
    elif args.all_energy_types:
        for etype in ['Energy', 'EnergyLeft', 'EnergyRight']:
            e_out = out_dir / etype
            e_out.mkdir(exist_ok=True)
            process_energy_type(json_files, etype, e_out, args.extreme_percentile, png)
    else:
        logger.error("Укажите --energy-type или --all-energy-types")
        sys.exit(1)

    logger.info(f"Готово. Результаты в {out_dir}")

if __name__ == "__main__":
    main()