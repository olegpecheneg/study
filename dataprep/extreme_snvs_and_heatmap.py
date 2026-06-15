#!/usr/bin/env python3
"""
Экстремальные SNV (топ-5 положительных и отрицательных) и тепловая карта переходов.
Барплоты: средние значения и размах (min-max).
Тепловая карта: средние значения (без звёздочек).
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import ttest_1samp

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

# ----------------------------------------------------------------------
def find_master_json_files(data_dir: Path, recursive=True):
    return list(data_dir.rglob("*_master_data.json")) if recursive else list(data_dir.glob("*_master_data.json"))

def collect_all_diffs(json_files: List[Path], energy_type: str) -> pd.DataFrame:
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
            diff = float(alt) - float(ref)
            records.append({'position': pos, 'diff': diff})
    return pd.DataFrame(records)

def collect_snv_stats(json_files: List[Path], energy_type: str = 'Energy') -> pd.DataFrame:
    stats = {}
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
        diffs = []
        for rec in data:
            if rec.get('energy_type') != energy_type:
                continue
            ref = rec.get('ref_energy')
            alt = rec.get('alt_energy')
            if ref is None or alt is None:
                continue
            diffs.append(float(alt) - float(ref))
        if len(diffs) < 3:
            continue
        stats[pos] = {
            'mean_diff': np.mean(diffs),
            'median_diff': np.median(diffs),
            'std_diff': np.std(diffs),
            'min_diff': np.min(diffs),
            'max_diff': np.max(diffs),
            'n_pairs': len(diffs),
            'gene': infer_mt_gene(pos)
        }
    df = pd.DataFrame.from_dict(stats, orient='index').reset_index().rename(columns={'index': 'position'})
    return df

def create_bar_plot(df_stats: pd.DataFrame, df_all: pd.DataFrame, positions: List[int],
                    title: str, output_file: Path, png=True):
    df_sub = df_stats[df_stats['position'].isin(positions)].copy()
    if df_sub.empty:
        return
    df_sub = df_sub.sort_values('mean_diff')
    df_sub['label'] = df_sub['position'].apply(lambda x: f"SNV {x}")

    fig = px.bar(
        df_sub,
        x='label',
        y='mean_diff',
        color='mean_diff',
        color_continuous_scale='RdBu',
        title=title,
        labels={'label': 'SNV', 'mean_diff': 'Средняя разность (alt – ref), ккал/моль'},
        hover_data={'position': True, 'gene': True, 'median_diff': ':.3f', 'min_diff': ':.3f', 'max_diff': ':.3f', 'n_pairs': True}
    )
    for _, row in df_sub.iterrows():
        # Линия min-max
        fig.add_shape(type='line', x0=row['label'], x1=row['label'], y0=row['min_diff'], y1=row['max_diff'],
                      line=dict(color='black', width=2), xref='x', yref='y')
        # Чёрточки
        fig.add_shape(type='line', x0=row['label'], x1=row['label'], y0=row['min_diff']-0.2, y1=row['min_diff']+0.2,
                      line=dict(color='black', width=2), xref='x', yref='y')
        fig.add_shape(type='line', x0=row['label'], x1=row['label'], y0=row['max_diff']-0.2, y1=row['max_diff']+0.2,
                      line=dict(color='black', width=2), xref='x', yref='y')
    fig.update_traces(marker_line_width=0, opacity=0.8)
    fig.update_layout(template='plotly_white', height=600, width=800)
    save_plot(fig, output_file, png)

def build_transition_map(json_files: List[Path]) -> Dict[int, Tuple[str, str]]:
    trans_map = {}
    for jf in json_files:
        try:
            with open(jf) as f:
                data = json.load(f)
            if data:
                pos = data[0].get('snp') or data[0].get('position')
                if pos is not None:
                    ref = data[0].get('ref_allele')
                    alt = data[0].get('alt_allele')
                    if ref and alt:
                        trans_map[int(pos)] = (ref, alt)
        except Exception:
            continue
    return trans_map

def create_transition_heatmap(json_files: List[Path], output_file: Path, png=True):
    energy_types = ['Energy', 'EnergyLeft', 'EnergyRight']
    trans_map = build_transition_map(json_files)
    heat_data = defaultdict(lambda: defaultdict(list))

    for etype in energy_types:
        df = collect_snv_stats(json_files, energy_type=etype)
        for _, row in df.iterrows():
            pos = int(row['position'])
            ref_alt = trans_map.get(pos)
            if ref_alt:
                transition = f"{ref_alt[0]}→{ref_alt[1]}"
                heat_data[transition][etype].append(row['mean_diff'])

    transitions = sorted(heat_data.keys())
    matrix = []
    text_matrix = []
    for tr in transitions:
        row_vals = []
        row_text = []
        for etype in energy_types:
            vals = heat_data[tr].get(etype, [])
            if len(vals) < 3:
                row_vals.append(np.nan)
                row_text.append('')
                continue
            mean_val = np.mean(vals)
            _, p = ttest_1samp(vals, 0)
            star = ''
            if p < 0.001:
                star = '***'
            elif p < 0.01:
                star = '**'
            elif p < 0.05:
                star = '*'
            row_vals.append(mean_val)
            row_text.append(f"{mean_val:.1f}{star}")
        matrix.append(row_vals)
        text_matrix.append(row_text)

    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=energy_types,
        y=transitions,
        text=text_matrix,
        texttemplate="%{text}",
        textfont={"size": 10},
        colorscale='RdBu',
        zmid=0,
        colorbar=dict(title='Средняя разность (alt - ref), ккал/моль')
    ))
    fig.update_layout(
        title='Тепловая карта средних разностей по нуклеотидным заменам<br>'
              '<sup>Обозначения: * p &lt; 0,05; ** p &lt; 0,01; *** p &lt; 0,001 (t-тест, H<sub>0</sub>: среднее = 0).</sup>',
        xaxis_title='Тип энергии',
        yaxis_title='Нуклеотидная замена',
        height=600,
        width=800
    )
    save_plot(fig, output_file, png)


def save_plot(fig, path, png=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(fig, str(path))
    if png:
        try:
            pio.write_image(fig, str(path.with_suffix('.png')), scale=2)
        except Exception as e:
            logger.warning(f"PNG не сохранён: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("-o", "--output_dir", type=Path, default=None)
    parser.add_argument("--energy-type", default="Energy", choices=["Energy","EnergyLeft","EnergyRight"])
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument("--no-recurse", action="store_true")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Нет директории {data_dir}")
        sys.exit(1)

    out_dir = (args.output_dir or data_dir / "extreme_analysis").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    json_files = find_master_json_files(data_dir, recursive=not args.no_recurse)
    logger.info(f"Найдено {len(json_files)} файлов")
    if not json_files:
        sys.exit(1)

    logger.info("Сбор статистики по SNV...")
    df_stats = collect_snv_stats(json_files, energy_type=args.energy_type)
    if df_stats.empty:
        logger.error("Нет данных")
        sys.exit(1)

    df_all = collect_all_diffs(json_files, energy_type=args.energy_type)
    if df_all.empty:
        logger.error("Нет исходных данных для барплотов")
        sys.exit(1)

    df_sorted = df_stats.sort_values('mean_diff')
    top_negative = df_sorted.head(args.top)['position'].tolist()
    top_positive = df_sorted.tail(args.top)['position'].tolist()
    logger.info(f"Топ-{args.top} отрицательных: {top_negative}")
    logger.info(f"Топ-{args.top} положительных: {top_positive}")

    extreme_df = pd.concat([df_sorted.head(args.top), df_sorted.tail(args.top)])
    extreme_df.to_csv(out_dir / "extreme_snvs.csv", index=False)

    create_bar_plot(df_stats, df_all, top_negative,
                    f"Топ-{args.top} SNV с наибольшим отрицательным изменением (alt – ref)",
                    out_dir / "bar_plot_negative.html", png=not args.no_png)
    create_bar_plot(df_stats, df_all, top_positive,
                    f"Топ-{args.top} SNV с наибольшим положительным изменением (alt – ref)",
                    out_dir / "bar_plot_positive.html", png=not args.no_png)

    logger.info("Построение тепловой карты замен...")
    create_transition_heatmap(json_files, output_file=out_dir / "transition_heatmap.html", png=not args.no_png)

    logger.info(f"Готово. Результаты в {out_dir}")

if __name__ == "__main__":
    main()