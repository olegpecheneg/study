import os
import re
import logging
from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import json
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import colorsys
import mplcursors


def setup_matplotlib() -> None:
    """Настройка бэкенда matplotlib и шрифтов по умолчанию."""
    matplotlib.use('Agg')
    plt.rcParams['font.family'] = 'DejaVu Sans'


def setup_logger(name: str = "scatter_plus_n_std") -> logging.Logger:
    """Настройка и возврат логгера для модуля."""
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        log_path = Path(__file__).resolve().parent.parent / 'visualization.log'
        file_handler = logging.FileHandler(log_path, encoding='utf-8')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(stream_handler)
        logger.addHandler(file_handler)
        # Общий лог на INFO
        logger.setLevel(logging.INFO)
        # Отдельный файл только для ошибок (короткие записи)
        error_log_path = Path(__file__).resolve().parent.parent / 'errors.log'
        error_handler = logging.FileHandler(error_log_path, encoding='utf-8')
        error_handler.setLevel(logging.ERROR)
        # Короткий формат для ошибок
        error_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(error_handler)
    return logger


setup_matplotlib()
logger = setup_logger()


def load_snp_data(snp_file_path: str) -> set:
    """Загрузка позиций SNP из CSV файла."""
    snp_positions = set()
    try:
        with open(snp_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                _try_add_snp_position(line, snp_positions)
        logger.info(f"Загружено {len(snp_positions)} SNP из файла {snp_file_path}")
    except Exception as e:
        logger.error(f"Ошибка загрузки файла SNP {snp_file_path}: {e}")
    return snp_positions


def _try_add_snp_position(line: str, snp_positions: set) -> None:
    parts = line.strip().split(',')
    if parts:
        try:
            position = int(parts[0])
            snp_positions.add(position)
        except (ValueError, IndexError):
            logger.warning(f"Ошибка обработки строки SNP: {line.strip()}")


def parse_construct_id(construct_id: str) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
    """Парсинг ID конструкта и возврат параметров плеч."""
    try:
        arm_size = _parse_cgs(construct_id)
        center = _parse_cen(construct_id)
        arm3_start, arm4_start = _parse_con(construct_id)
        if None in (arm_size, center, arm3_start, arm4_start):
            return None, None, None, None
        return arm_size, center, arm3_start, arm4_start
    except Exception as e:
        logger.error(f"Ошибка парсинга ID конструкта {construct_id}: {e}")
        return None, None, None, None


def _parse_cgs(construct_id: str) -> Optional[int]:
    cgs_match = re.search(r'CGS-(\d+)-(\d+)-(\d+)-(\d+)-(\d+)-(\d+)', construct_id)
    if cgs_match:
        return int(cgs_match.group(6))
    logger.warning(f"Не найден блок CGS в ID: {construct_id}")
    return None


def _parse_cen(construct_id: str) -> Optional[int]:
    cen_match = re.search(r'CEN-(\d+)', construct_id)
    if cen_match:
        return int(cen_match.group(1))
    logger.warning(f"Не найден блок CEN в ID: {construct_id}")
    return None


def _parse_con(construct_id: str) -> Tuple[Optional[int], Optional[int]]:
    con_match = re.search(r'CON-(\d+)-(\d+)', construct_id)
    if con_match:
        return int(con_match.group(1)), int(con_match.group(2))
    logger.warning(f"Не найден блок CON в ID: {construct_id}")
    return None, None


def calculate_arm_ranges(arm_size: int, center: int, arm3_start: int, arm4_start: int, subseq_start: int = 10) -> List[Tuple[int, int]]:
    """Вычисление абсолютных диапазонов для четырех плеч конструкта."""
    arm1_start = center - arm_size
    arm2_start = center
    ranges = [
        (arm1_start + subseq_start, arm1_start + arm_size + subseq_start),
        (arm2_start + subseq_start, arm2_start + arm_size + subseq_start),
        (arm3_start + subseq_start, arm3_start + arm_size + subseq_start),
        (arm4_start + subseq_start, arm4_start + arm_size + subseq_start)
    ]
    logger.debug(f"Диапазоны конструкта (абсолютные): {ranges}")
    return ranges


def get_snps_in_construct(snp_positions: set, arm_size: int, center: int, arm3_start: int, arm4_start: int, subseq_start: int = 10) -> List[int]:
    """Возвращает SNP, попадающие в любое из плеч конструкта."""
    ranges = calculate_arm_ranges(arm_size, center, arm3_start, arm4_start, subseq_start)
    snps_in_construct = []
    for snp in snp_positions:
        for start, end in ranges:
            if start <= snp <= end:
                snps_in_construct.append(snp)
                break
    logger.debug(f"SNP в конструкте: {snps_in_construct}")
    return snps_in_construct


def generate_distinct_colors(n: int) -> List[Tuple[float, float, float]]:
    """Генерация списка визуально различимых RGB цветов."""
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        saturation = 0.9
        value = 0.9
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append(rgb)
    return colors


def calculate_outlier_stats(ref_data: List[float], alt_data: List[float]) -> Tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    """Вычисление статистики выбросов для сравнения референсных и альтернативных данных."""
    diff = np.array(ref_data) - np.array(alt_data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    upper_outliers = diff > mean_diff + 2 * std_diff
    lower_outliers = diff < mean_diff - 2 * std_diff
    normal_points = ~(upper_outliers | lower_outliers)
    return mean_diff, std_diff, upper_outliers, lower_outliers, normal_points


def plot_scatter_points(ax: Any, ref_data: List[float], alt_data: List[float], snp_values: List[Optional[int]], snp_colors: Dict[int, Tuple[float, float, float]], upper_outliers: np.ndarray, lower_outliers: np.ndarray, normal_points: np.ndarray) -> None:
    """Отрисовка точек scatter plot с цветовой кодировкой по SNP и выделением выбросов."""
    gray_rgb = (0.5, 0.5, 0.5)
    colors = [snp_colors.get(snp, gray_rgb) for snp in snp_values]
    ax.scatter(
        np.array(ref_data)[normal_points], 
        np.array(alt_data)[normal_points], 
        c=np.array(colors)[normal_points], 
        alpha=0.7, s=80, edgecolor='black', linewidth=0.5
    )
    ax.scatter(
        np.array(ref_data)[upper_outliers], 
        np.array(alt_data)[upper_outliers], 
        c=np.array(colors)[upper_outliers], 
        alpha=0.9, s=150, edgecolor='green', linewidth=3
    )
    ax.scatter(
        np.array(ref_data)[lower_outliers], 
        np.array(alt_data)[lower_outliers], 
        c=np.array(colors)[lower_outliers], 
        alpha=0.9, s=150, edgecolor='red', linewidth=3
    )


def add_diagonal_line(ax: Any, min_e: float, max_e: float) -> None:
    """Добавление пунктирной диагональной линии (y == x)."""
    ax.plot([min_e, max_e], [min_e, max_e], 'k--', linewidth=2, alpha=0.7)


def add_outlier_zones(ax: Any, x: np.ndarray, mean_diff: float, std_diff: float, min_e: float, max_e: float) -> None:
    """Затенение областей, соответствующих порогам выбросов."""
    line_upper = x - (mean_diff + 2 * std_diff)
    line_lower = x - (mean_diff - 2 * std_diff)
    line_upper_1 = x - (mean_diff + 3 * std_diff)
    line_lower_1 = x - (mean_diff - 3 * std_diff)
    line_upper_2 = x - (mean_diff + 4 * std_diff)
    line_lower_2 = x - (mean_diff - 4 * std_diff)
    ax.fill_between(x, min_e, line_upper, color='red', alpha=0.1, label='Верхние выбросы (+2std)')
    ax.fill_between(x, line_lower, max_e, color='green', alpha=0.1, label='Нижние выбросы (-2std)')
    ax.fill_between(x, min_e, line_upper_1, color='red', alpha=0.1, label='Верхние выбросы (+3std)')
    ax.fill_between(x, line_lower_1, max_e, color='green', alpha=0.1, label='Нижние выбросы (-3std)')
    ax.fill_between(x, min_e, line_upper_2, color='red', alpha=0.1, label='Верхние выбросы (+4std)')
    ax.fill_between(x, line_lower_2, max_e, color='green', alpha=0.1, label='Нижние выбросы (-4std)')


def create_legend_elements(snp_colors: Dict[int, Tuple[float, float, float]]) -> List[Line2D]:
    """Создание элементов легенды для графика на основе цветового отображения SNP."""
    legend_elements = []
    gray_rgb = (0.5, 0.5, 0.5)
    legend_snps = list(snp_colors.keys())
    if len(legend_snps) > 20:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
                   markersize=12, label=f'Другие SNP ({len(snp_colors)-20})', markeredgecolor='black')
        )
        legend_snps = legend_snps[:20]
    for snp in legend_snps:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor=snp_colors[snp], 
                   markersize=12, label=f'SNP {snp}', markeredgecolor='black')
        )
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
               markersize=12, label='Без SNP/Ошибка', markeredgecolor='black')
    )
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
               markeredgecolor='green', markersize=12, label='Верхние выбросы (+2std)', linewidth=3)
    )
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gray_rgb, 
               markeredgecolor='red', markersize=12, label='Нижние выбросы (-2std)', linewidth=3)
    )
    legend_elements.append(
        Line2D([0], [0], color='k', linestyle='--', linewidth=2, label='Диагональ (x=y)')
    )
    return legend_elements


def find_points_by_sequence(search_file: str, search_sequence: str, case_sensitive: bool = False) -> List[Dict]:
    """Поиск точек по последовательности в файле поиска."""
    if not case_sensitive:
        search_sequence = search_sequence.upper()
    
    try:
        df = pd.read_csv(search_file)
        results = []
        
        for _, row in df.iterrows():
            for window_num in range(1, 7):
                window_seq = str(row[f'window{window_num}'])
                if not case_sensitive:
                    window_seq = window_seq.upper()
                
                if search_sequence in window_seq:
                    results.append({
                        'point_index': row['point_index'],
                        'construct_id': row['construct_id'],
                        'ref_energy': row['ref_energy'],
                        'alt_energy': row['alt_energy'],
                        'snp': row['snp'],
                        'window_found': f'window{window_num}',
                        'sequence_found': window_seq
                    })
                    break
        
        return results
        
    except Exception as e:
        logger.error(f"Ошибка поиска в файле {search_file}: {e}")
        return []


def plot_energy_comparison(ref_data: List[float], alt_data: List[float], 
                         snp_values: List[Optional[int]], 
                         snp_colors: Dict[int, Tuple[float, float, float]], 
                         energy_type: str, output_dir: str, individual_id: str,
                         construct_details: List[Dict]) -> Optional[Dict[str, Any]]:
    """Создание scatter plot с интерактивными подсказками и экспортом детальных данных."""
    
    if not ref_data or not alt_data:
        logger.warning(f"Нет данных для построения графика {energy_type}")
        return

    min_len = min(len(ref_data), len(alt_data), len(snp_values), len(construct_details))
    ref_data = ref_data[:min_len]
    alt_data = alt_data[:min_len]
    snp_values = snp_values[:min_len]
    construct_details = construct_details[:min_len]
    
    mean_diff, std_diff, upper_outliers, lower_outliers, normal_points = calculate_outlier_stats(ref_data, alt_data)
    
    # Сохранение детальных данных в текстовый файл
    detailed_data_path = Path(output_dir) / f"{energy_type}_detailed_data.txt"
    with open(detailed_data_path, 'w', encoding='utf-8') as f:
        f.write(f"Детальные данные для {energy_type} - Индивидуум {individual_id}\n")
        f.write("=" * 80 + "\n\n")
        
        for i in range(min_len):
            details = construct_details[i]
            f.write(f"Точка {i}:\n")
            f.write(f"  Координаты: Ref={ref_data[i]:.2f}, Alt={alt_data[i]:.2f}\n")
            f.write(f"  ID конструкта: {details['construct_id']}\n")
            f.write(f"  SNP: {snp_values[i] if snp_values[i] else 'Нет'}\n")
            f.write(f"  Выброс: {'Верхний' if i < len(upper_outliers) and upper_outliers[i] else 'Нижний' if i < len(lower_outliers) and lower_outliers[i] else 'Норма'}\n")
            
            # ПРОВЕРКА НАЛИЧИЯ ДАННЫХ О КОНСТРУКТЕ
            if details['center_position'] != 'N/A':
                f.write("  Координаты плеч:\n")
                f.write(f"    Центр: {details['center_position']}\n")
                f.write(f"    Плечо1: {details['arm1_start']}\n")
                f.write(f"    Плечо2: {details['arm2_start']}\n")
                f.write(f"    Плечо3: {details['arm3_start']}\n")
                f.write(f"    Плечо4: {details['arm4_start']}\n")
            else:
                f.write("  Координаты плеч: Данные отсутствуют\n")
            
            # ПРОВЕРКА НАЛИЧИЯ ПОСЛЕДОВАТЕЛЬНОСТЕЙ
            if details['sequence_window1'] != 'N/A':
                f.write("  Последовательности Window:\n")
                f.write(f"    Window1: {details['sequence_window1']}\n")
                f.write(f"    Window2: {details['sequence_window2']}\n")
                f.write(f"    Window3: {details['sequence_window3']}\n")
                f.write(f"    Window4: {details['sequence_window4']}\n")
                f.write(f"    Window5: {details['sequence_window5']}\n")
                f.write(f"    Window6: {details['sequence_window6']}\n")
            else:
                f.write("  Последовательности Window: Данные отсутствуют\n")
                
            f.write("-" * 80 + "\n")
    
    logger.info(f"Детальные данные сохранены: {detailed_data_path}")

    # Создание CSV файла для поиска по последовательностям (только если есть данные о последовательностях)
    search_data_path = Path(output_dir) / f"{energy_type}_sequence_search.csv"
    search_data = []
    for i in range(min_len):
        details = construct_details[i]
        # ТОЛЬКО ЕСЛИ ЕСТЬ ДАННЫЕ О ПОСЛЕДОВАТЕЛЬНОСТЯХ
        if details['sequence_window1'] != 'N/A':
            search_data.append({
                'point_index': i,
                'ref_energy': ref_data[i],
                'alt_energy': alt_data[i],
                'construct_id': details['construct_id'],
                'snp': snp_values[i],
                'window1': details['sequence_window1'],
                'window2': details['sequence_window2'], 
                'window3': details['sequence_window3'],
                'window4': details['sequence_window4'],
                'window5': details['sequence_window5'],
                'window6': details['sequence_window6']
            })
    
    if search_data:
        search_df = pd.DataFrame(search_data)
        search_df.to_csv(search_data_path, index=False, encoding='utf-8')
        logger.info(f"Данные для поиска сохранены: {search_data_path}")
    else:
        logger.info(f"Данные для поиска не сохранены (отсутствуют последовательности Window)")

    # Создание графика
    fig, ax = plt.subplots(figsize=(16, 12))
    
    plot_scatter_points(ax, ref_data, alt_data, snp_values, snp_colors, upper_outliers, lower_outliers, normal_points)
    
    min_e = min(min(ref_data), min(alt_data))
    max_e = max(max(ref_data), max(alt_data))
    add_diagonal_line(ax, min_e, max_e)
    
    x = np.linspace(min_e, max_e, 100)
    add_outlier_zones(ax, x, mean_diff, std_diff, min_e, max_e)

    plt.title(f'Сравнение {energy_type} для индивидуума {individual_id}\nФайл данных: {detailed_data_path.name}', fontsize=16)
    plt.xlabel('Референсная энергия (kcal/mol)', fontsize=16)
    plt.ylabel('Альтернативная энергия (kcal/mol)', fontsize=16)
    
    legend_elements = create_legend_elements(snp_colors)
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, title="Легенда", title_fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.2)
    ax.set_xlim(min(ref_data), max(ref_data))
    ax.set_ylim(min(alt_data), max(alt_data))
    plt.tight_layout()

    # Интерактивные подсказки (только если mplcursors доступен)
    if mplcursors is not None:
        try:
            scatter = ax.collections[0]
            
            annotations = []
            for i in range(min_len):
                details = construct_details[i]
                snp_info = f"SNP: {snp_values[i]}" if snp_values[i] is not None else "Без SNP"
                outlier_info = ""
                if i < len(upper_outliers) and upper_outliers[i]:
                    outlier_info = " (Верхний выброс)"
                elif i < len(lower_outliers) and lower_outliers[i]:
                    outlier_info = " (Нижний выброс)"
                
                # ФОРМИРОВАНИЕ АННОТАЦИИ С УЧЕТОМ ОТСУТСТВУЮЩИХ ДАННЫХ
                if details['center_position'] != 'N/A':
                    coords_info = (
                        f"Координаты: Центр={details['center_position']}\n"
                        f"Плечо1={details['arm1_start']}, Плечо2={details['arm2_start']}\n"
                        f"Плечo3={details['arm3_start']}, Плечо4={details['arm4_start']}\n"
                    )
                else:
                    coords_info = "Координаты: Данные отсутствуют\n"
                
                annotation = (
                    f"Точка {i}\n"
                    f"Конструкт: {details['construct_id']}\n"
                    f"Ref: {ref_data[i]:.2f}, Alt: {alt_data[i]:.2f}\n"
                    f"{snp_info}{outlier_info}\n"
                    f"{coords_info}"
                    f"Полные данные в файле:\n{detailed_data_path.name}"
                )
                annotations.append(annotation)
            
            cursor = mplcursors.cursor(scatter, hover=True)
            
            @cursor.connect("add")
            def on_add(sel):
                sel.annotation.set_text(annotations[sel.index])
                sel.annotation.set_alpha(0.95)
                sel.annotation.set_backgroundcolor("white")
                sel.annotation.set_bbox(dict(
                    boxstyle="round,pad=0.5", 
                    facecolor="white", 
                    edgecolor="black", 
                    alpha=0.95
                ))
                sel.annotation.get_bbox_patch().set_linewidth(2)
            
            logger.info("Интерактивные подсказки активированы")
        except Exception as e:
            logger.warning(f"Ошибка при настройке интерактивных подсказок: {e}")
    else:
        logger.info("Интерактивные подсказки отключены (mplcursors недоступен)")

    output_path = Path(output_dir) / f"{energy_type}_comparison.png"
    plt.savefig(output_path, dpi=250, bbox_inches='tight')
    # Попробуем сохранить интерактивный HTML (mpld3) — если пакет доступен
    try:
        import mpld3 as _mpld3  # local name to avoid shadowing
        html_path = Path(output_dir) / f"{energy_type}_comparison.html"
        try:
            _mpld3.save_html(fig, str(html_path))
            logger.info(f"Интерактивный HTML сохранён: {html_path}")
        except Exception as e:
            logger.warning(f"Не удалось сохранить интерактивный HTML для {energy_type}: {e}")
    except Exception:
        # mpld3 отсутствует — не критично, PNG остаётся
        logger.info("Интерактивный HTML не сохранён (mpld3 недоступен)")

    plt.close()
    logger.info(f"График сохранён: {output_path}")
    # Попробуем сохранить интерактивный HTML через Plotly (более надёжно, чем mpld3)
    try:
        import plotly.graph_objects as go
        import plotly.io as pio

        def rgb_to_hex(rgb_tuple: Tuple[float, float, float]) -> str:
            return '#%02x%02x%02x' % tuple(int(255 * max(0, min(1, v))) for v in rgb_tuple)

        # Соберём данные для Plotly
        plot_data = []
        colors_plot = []
        hover_texts = []
        for i in range(min_len):
            details = construct_details[i]
            snp_val = snp_values[i]
            color_rgb = snp_colors.get(snp_val, (0.5, 0.5, 0.5))
            colors_plot.append(rgb_to_hex(color_rgb))

            outlier_info = ''
            if i < len(upper_outliers) and upper_outliers[i]:
                outlier_info = 'Верхний выброс'
            elif i < len(lower_outliers) and lower_outliers[i]:
                outlier_info = 'Нижний выброс'
            else:
                outlier_info = 'Норма'

            seq_block = ''
            if details.get('sequence_window1', 'N/A') != 'N/A':
                seq_block = f"Window1: {details.get('sequence_window1','')}"
                for w in range(2, 7):
                    seq_block += f"<br>Window{w}: {details.get(f'sequence_window{w}','')}"

            hover = (
                f"Точка {i}<br>"
                f"Конструкт: {details.get('construct_id')}<br>"
                f"Ref: {ref_data[i]:.2f}, Alt: {alt_data[i]:.2f}<br>"
                f"SNP: {snp_val if snp_val is not None else 'Нет'}<br>"
                f"Статус: {outlier_info}<br>"
                f"Центр: {details.get('center_position','N/A')}<br>"
                f"Плечо1: {details.get('arm1_start','N/A')}, Плечо2: {details.get('arm2_start','N/A')}<br>"
                f"Плечо3: {details.get('arm3_start','N/A')}, Плечо4: {details.get('arm4_start','N/A')}<br>"
                f"{seq_block}<br>"
                f"Полные данные в файле: {Path(detailed_data_path).name}"
            )
            hover_texts.append(hover)

        # Соберём структурированные записи для каждой точки — пригодится для встраивания JSON в HTML
        point_records = []
        for i in range(min_len):
            det = construct_details[i]
            outlier_status = 'Норма'
            if i < len(upper_outliers) and upper_outliers[i]:
                outlier_status = 'Верхний'
            elif i < len(lower_outliers) and lower_outliers[i]:
                outlier_status = 'Нижний'

            point_records.append({
                'index': i,
                'construct_id': det.get('construct_id'),
                'ref_energy': float(ref_data[i]),
                'alt_energy': float(alt_data[i]),
                'snp': snp_values[i],
                'outlier': outlier_status,
                'center_position': det.get('center_position', 'N/A'),
                'arm1_start': det.get('arm1_start', 'N/A'),
                'arm2_start': det.get('arm2_start', 'N/A'),
                'arm3_start': det.get('arm3_start', 'N/A'),
                'arm4_start': det.get('arm4_start', 'N/A'),
                'sequence_window1': det.get('sequence_window1', 'N/A'),
                'sequence_window2': det.get('sequence_window2', 'N/A'),
                'sequence_window3': det.get('sequence_window3', 'N/A'),
                'sequence_window4': det.get('sequence_window4', 'N/A'),
                'sequence_window5': det.get('sequence_window5', 'N/A'),
                'sequence_window6': det.get('sequence_window6', 'N/A'),
                'energy_type': energy_type,
                'detailed_file': Path(detailed_data_path).name
            })

        if len(ref_data) > 0:
            # Группируем точки по SNP для легенды
            unique_snps = sorted(set(snp for snp in snp_values if snp is not None))
            
            fig_p = go.Figure()
            
            # Добавляем точки по группам SNP для легенды
            for snp in unique_snps:
                snp_indices = [i for i, sv in enumerate(snp_values) if sv == snp]
                snp_ref = [ref_data[i] for i in snp_indices]
                snp_alt = [alt_data[i] for i in snp_indices]
                snp_colors_local = [colors_plot[i] for i in snp_indices]
                snp_hovers = [hover_texts[i] for i in snp_indices]

                fig_p.add_trace(go.Scatter(
                    x=snp_ref, y=snp_alt,
                    mode='markers',
                    name=f'SNP {snp}',
                    marker=dict(
                        color=snp_colors_local if snp_colors_local else ['#808080'] * len(snp_ref),
                        size=10,
                        line=dict(width=1, color='Black')
                    ),
                    hoverinfo='text',
                    text=snp_hovers,
                    customdata=snp_indices
                ))

            # Точки без SNP (если есть)
            mask_no_snp = [sv is None for sv in snp_values]
            if any(mask_no_snp):
                no_snp_indices = [i for i, sv in enumerate(snp_values) if sv is None]
                no_snp_ref = [ref_data[i] for i in no_snp_indices]
                no_snp_alt = [alt_data[i] for i in no_snp_indices]
                no_snp_hovers = [hover_texts[i] for i in no_snp_indices]
                no_snp_colors = [colors_plot[i] for i in no_snp_indices]

                fig_p.add_trace(go.Scatter(
                    x=no_snp_ref, y=no_snp_alt,
                    mode='markers',
                    name='Без SNP',
                    marker=dict(
                        color=no_snp_colors if no_snp_colors else ['#808080'] * len(no_snp_ref),
                        size=10,
                        line=dict(width=1, color='Black')
                    ),
                    hoverinfo='text',
                    text=no_snp_hovers,
                    customdata=no_snp_indices
                ))

            # Диагональ
            fig_p.add_trace(go.Scatter(
                x=[min_e, max_e], y=[min_e, max_e],
                mode='lines',
                line=dict(dash='dash', color='black'),
                name='Диагональ',
                hoverinfo='skip'
            ))

            # Добавим категорию "Выбросы" как отдельный trace (яркие, крупные маркеры)
            outlier_mask = [u or l for u, l in zip(upper_outliers, lower_outliers)]
            if any(outlier_mask):
                outlier_indices = [i for i, m in enumerate(outlier_mask) if m]
                outlier_ref = [ref_data[i] for i in outlier_indices]
                outlier_alt = [alt_data[i] for i in outlier_indices]
                outlier_hovers = [hover_texts[i] for i in outlier_indices]
                # Цвета для каждого выброса — соответствуют цветам исходных точек
                outlier_colors = [colors_plot[i] for i in outlier_indices]
                fig_p.add_trace(go.Scatter(
                    x=outlier_ref, y=outlier_alt,
                    mode='markers',
                    name='Выбросы',
                    marker=dict(
                        color=outlier_colors,
                        size=16,
                        symbol='diamond',
                        line=dict(width=1, color='black')
                    ),
                    hoverinfo='text',
                    text=outlier_hovers,
                    customdata=outlier_indices
                ))

            # Диагональ
            fig_p.add_trace(go.Scatter(
                x=[min_e, max_e], y=[min_e, max_e],
                mode='lines',
                line=dict(dash='dash', color='black'),
                name='Диагональ',
                hoverinfo='skip'
            ))

            # Подготовим видимость для кнопок управления (более надежный подход)
            total_traces = len(fig_p.data)
            visible_all = [True] * total_traces
            visible_snp = [getattr(t, 'name', '').startswith('SNP ') for t in fig_p.data]
            visible_no_snp = [getattr(t, 'name', '') == 'Без SNP' for t in fig_p.data]
            visible_outliers = [getattr(t, 'name', '') == 'Выбросы' for t in fig_p.data]

            # Настройка легенды и интерфейса
            fig_p.update_layout(
                title=f"Сравнение {energy_type} для индивидуума {individual_id}",
                xaxis_title='Ref energy (kcal/mol)',
                yaxis_title='Alt energy (kcal/mol)',
                width=1200, height=800,
                showlegend=True,
                legend=dict(
                    yanchor="top",
                    y=0.99,
                    xanchor="left",
                    x=0.01,
                    bgcolor='rgba(255, 255, 255, 0.8)',
                    bordercolor='rgba(0, 0, 0, 0.3)',
                    borderwidth=1
                ),
                updatemenus=[
                    dict(
                        buttons=list([
                            dict(
                                args=[{"visible": visible_all}],
                                label="Показать все",
                                method="update"
                            ),
                            dict(
                                args=[{"visible": visible_snp}],
                                label="Только с SNP",
                                method="update"
                            ),
                            dict(
                                args=[{"visible": visible_no_snp}],
                                label="Только без SNP",
                                method="update"
                            ),
                            dict(
                                args=[{"visible": visible_outliers}],
                                label="Выбросы",
                                method="update"
                            )
                        ]),
                        direction="down",
                        showactive=True,
                        x=0.1,
                        y=1.1,
                        xanchor="left",
                        yanchor="top"
                    )
                ],
                annotations=[
                    dict(
                        text="Фильтр:",
                        x=0.02, y=1.1,
                        xref="paper", yref="paper",
                        showarrow=False
                    )
                ],
                clickmode='event+select'
            )

            # Добавляем JavaScript/CSS/HTML для модального окна, панели поиска и sidebar фильтров
            modal_js = """
            <style>
            /* Top search panel (top-right, slight offset from top) */
            #topSearchPanel {
                position: absolute;
                z-index: 1002;
                right: 12px;
                top: 30px;
                background: rgba(255,255,255,0.97);
                padding: 10px;
                border-radius: 6px;
                box-shadow: 0 2px 12px rgba(0,0,0,0.12);
                font-family: Arial, sans-serif;
                font-size: 13px;
                width: 360px;
            }
            /* Bottom filters panel (bottom-right, slight offset from bottom) */
            #filterPanel {
                position: absolute;
                z-index: 1002;
                right: 12px;
                bottom: 30px;
                background: rgba(255,255,255,0.97);
                padding: 12px;
                border-radius: 6px;
                box-shadow: 0 2px 12px rgba(0,0,0,0.12);
                font-family: Arial, sans-serif;
                font-size: 13px;
                width: 360px;
                max-height: 48%;
                overflow: auto;
            }
            .panel-row { display:flex; align-items:center; gap:8px; margin-top:6px }
            .panel-row label { width:120px; font-size:12px; color:#333 }
            .panel-row input, .panel-row select, .panel-row textarea { padding:6px; border:1px solid #ccc; border-radius:4px; flex:1 }
            .panel-actions { display:flex; gap:8px; margin-top:8px }
            .panel-actions button { padding:6px 8px; border-radius:4px; border:1px solid #999; background:#f7f7f7; cursor:pointer }
            .panel-actions button:hover { background:#eee }

            #resultsSummary { margin-top:8px; font-size:12px }
            #resultsTable { margin-top:6px; border-top:1px solid #eee; padding-top:6px }
            #resultsTable table { width:100%; border-collapse:collapse; font-size:12px }
            #resultsTable th, #resultsTable td { padding:6px; border:1px solid #eee; text-align:left }

            .modal { display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.4); }
            .modal-content { background-color: #fff; margin: 10% auto; padding: 18px; border-radius: 6px; width: 60%; max-width: 700px; max-height: 70vh; overflow-y: auto; box-shadow: 0 6px 18px rgba(0,0,0,0.16) }
            .close { color: #777; float: right; font-size: 22px; font-weight: bold; cursor: pointer; position: sticky; top: 0; background: white; padding: 0 4px; }
            .close:hover { color: #111 }
            a.result-link { color:#0366d6; cursor:pointer; text-decoration:underline }
            </style>

            <div id="topSearchPanel">
                <div style="font-weight:bold; margin-bottom:6px">Быстрый поиск (по всем полям)</div>
                <div class="panel-row">
                    <label>Поиск (всё):</label>
                    <input id="searchAny" placeholder="введите текст, ищет по всем полям" />
                </div>
                <div class="panel-actions">
                    <button id="btnTopFind">Найти</button>
                    <button id="btnTopClear">Сброс</button>
                </div>
                <div id="topInfo" style="margin-top:8px;font-size:12px;color:#444">Результатов: <span id="resultsCount">0</span></div>
            </div>

            <div id="filterPanel">
                <div style="font-weight:bold">Фильтрация (подробно)</div>

                <div class="panel-row">
                    <label>ConstructID</label>
                    <input id="filterConstruct" placeholder="точное или подстрока" />
                </div>
                <div style="display:flex;gap:8px;margin-top:6px;align-items:center">
                    <label style="width:120px">Режим для ConstructID</label>
                    <label><input type="radio" name="constructMode" value="exact" checked/> Точное</label>
                    <label><input type="radio" name="constructMode" value="substring"/> Подстрока</label>
                </div>

                <div class="panel-row">
                    <label>SNP (multi)</label>
                    <select id="filterSnpMulti" multiple size="4"></select>
                </div>
                <div style="display:flex;gap:8px;margin-top:6px;align-items:center">
                    <label style="width:120px">Режим SNP</label>
                    <label><input type="radio" name="snpMode" value="exact" checked/> Точное</label>
                    <label><input type="radio" name="snpMode" value="substring"/> Подстрока</label>
                    <button id="btnSelectAllSnps" style="margin-left:auto">Выбрать все</button>
                </div>

                <div class="panel-row">
                    <label>Energy min / max</label>
                    <input id="filterMinEnergy" placeholder="min" />
                    <input id="filterMaxEnergy" placeholder="max" />
                </div>

                <div style="margin-top:8px;font-weight:600">Поля плеч (arms)</div>
                <div class="panel-row">
                    <label>Arm1</label>
                    <input id="filterArm1" placeholder="значение или часть" />
                </div>
                <div class="panel-row">
                    <label>Arm2</label>
                    <input id="filterArm2" placeholder="значение или часть" />
                </div>
                <div class="panel-row">
                    <label>Arm3</label>
                    <input id="filterArm3" placeholder="значение или часть" />
                </div>
                <div class="panel-row">
                    <label>Arm4</label>
                    <input id="filterArm4" placeholder="значение или часть" />
                </div>

                <div style="margin-top:8px;font-weight:600">Поля Window (последовательности)</div>
                <div class="panel-row"><label>Window1</label><input id="filterW1" placeholder="часть последовательности" /></div>
                <div class="panel-row"><label>Window2</label><input id="filterW2" placeholder="часть последовательности" /></div>
                <div class="panel-row"><label>Window3</label><input id="filterW3" placeholder="часть последовательности" /></div>
                <div class="panel-row"><label>Window4</label><input id="filterW4" placeholder="часть последовательности" /></div>
                <div class="panel-row"><label>Window5</label><input id="filterW5" placeholder="часть последовательности" /></div>
                <div class="panel-row"><label>Window6</label><input id="filterW6" placeholder="часть последовательности" /></div>

                <div class="panel-actions">
                    <button id="btnApply">Применить</button>
                    <button id="btnClear">Очистить</button>
                    <button id="btnDownloadCSV">Сохранить CSV</button>
                    <button id="btnDownloadJSON">Сохранить JSON</button>
                </div>

                <div id="resultsTable"></div>
            </div>

            <div id="dataModal" class="modal"><div class="modal-content"><span class="close" id="modalClose">&times;</span><div id="dataTable"></div></div></div>

            <script>
            // Ensure POINT_DATA is available by parsing a safe JSON script block if needed
            if(typeof POINT_DATA === 'undefined' || !Array.isArray(POINT_DATA)){
                try{
                    var _el = document.getElementById('POINT_DATA_JSON');
                    if(_el && _el.textContent){
                        window.POINT_DATA = JSON.parse(_el.textContent);
                    } else {
                        window.POINT_DATA = [];
                    }
                }catch(e){ window.POINT_DATA = []; console && console.error && console.error('POINT_DATA parse error', e); }
            }

            function getPlotDiv(){ var p = document.getElementById('plotly'); if(!p) p = document.getElementsByClassName('plotly-graph-div')[0]; return p; }

            var ORIGINAL_MARKER_SIZES = [];
            var ORIGINAL_MARKER_COLORS = [];
            var LAST_FILTERED = [];

            function storeOriginalMarkerState(){
                var plot = getPlotDiv(); if(!plot) return; var data = plot.data || [];
                for(var t=0;t<data.length;t++){ var m = data[t].marker||{}; var size = m.size||10; var color = m.color||'#808080';
                    if(!Array.isArray(size)){ var arr = new Array((data[t].x && data[t].x.length)||0).fill(size); size = arr }
                    if(!Array.isArray(color)){ var carr = new Array((data[t].x && data[t].x.length)||0).fill(color); color = carr }
                    ORIGINAL_MARKER_SIZES[t]=size.slice(); ORIGINAL_MARKER_COLORS[t]=color.slice(); }
            }

            function resetAllMarkers(){ var plot=getPlotDiv(); if(!plot) return; for(var t=0;t<plot.data.length;t++){ var s = ORIGINAL_MARKER_SIZES[t]||[]; var c = ORIGINAL_MARKER_COLORS[t]||[]; try{ Plotly.restyle(plot, {'marker.size':[s], 'marker.color':[c]}, [t]); }catch(e){} } document.getElementById('resultsCount').innerText='0'; document.getElementById('resultsTable').innerHTML=''; LAST_FILTERED=[] }

            function selectAllSnps(){ var sel = document.getElementById('filterSnpMulti'); if(!sel) return; for(var i=0;i<sel.options.length;i++){ sel.options[i].selected=true } }

            function applyFilters(){
                var construct = document.getElementById('filterConstruct').value.trim();
                var constructMode = document.querySelector('input[name="constructMode"]:checked') ? document.querySelector('input[name="constructMode"]:checked').value : 'exact';
                var snpSelect = document.getElementById('filterSnpMulti');
                var selectedSnps = [];
                if(snpSelect){ for(var i=0;i<snpSelect.options.length;i++){ if(snpSelect.options[i].selected) selectedSnps.push(snpSelect.options[i].value) } }
                var snpMode = document.querySelector('input[name="snpMode"]:checked') ? document.querySelector('input[name="snpMode"]:checked').value : 'exact';
                var minE = parseFloat(document.getElementById('filterMinEnergy').value);
                var maxE = parseFloat(document.getElementById('filterMaxEnergy').value);
                var arm1 = document.getElementById('filterArm1').value.trim();
                var arm2 = document.getElementById('filterArm2').value.trim();
                var arm3 = document.getElementById('filterArm3').value.trim();
                var arm4 = document.getElementById('filterArm4').value.trim();
                var w1 = document.getElementById('filterW1').value.trim();
                var w2 = document.getElementById('filterW2').value.trim();
                var w3 = document.getElementById('filterW3').value.trim();
                var w4 = document.getElementById('filterW4').value.trim();
                var w5 = document.getElementById('filterW5').value.trim();
                var w6 = document.getElementById('filterW6').value.trim();
                var any = document.getElementById('searchAny').value.trim().toLowerCase();

                var matches = [];
                for(var i=0;i<POINT_DATA.length;i++){
                    var r = POINT_DATA[i];
                    var ok = true;

                    // ConstructID filter
                    if(construct){
                        if(constructMode==='exact'){
                            if(String(r.construct_id)!==construct) ok=false;
                        } else {
                            if(!(String(r.construct_id||'').toLowerCase().indexOf(construct.toLowerCase())!==-1)) ok=false;
                        }
                    }

                    // SNP multi filter
                    if(selectedSnps.length>0){
                        var rSnp = r.snp===null? '': String(r.snp);
                        var matchedSnp = false;
                        for(var si=0; si<selectedSnps.length; si++){
                            var sval = selectedSnps[si];
                            if(snpMode==='exact'){
                                if(rSnp===sval) matchedSnp=true;
                            } else {
                                if(rSnp.toLowerCase().indexOf(sval.toLowerCase())!==-1) matchedSnp=true;
                            }
                        }
                        if(!matchedSnp) ok=false;
                    }

                    // Energy range (either ref or alt in range)
                    if(!isNaN(minE) || !isNaN(maxE)){
                        var lo = isNaN(minE)? -1e99 : minE;
                        var hi = isNaN(maxE)? 1e99 : maxE;
                        if(!((r.ref_energy>=lo && r.ref_energy<=hi) || (r.alt_energy>=lo && r.alt_energy<=hi))) ok=false;
                    }

                    // Arms filters (substring match)
                    if(arm1){ if(String(r.arm1_start||'').toLowerCase().indexOf(arm1.toLowerCase())===-1) ok=false }
                    if(arm2){ if(String(r.arm2_start||'').toLowerCase().indexOf(arm2.toLowerCase())===-1) ok=false }
                    if(arm3){ if(String(r.arm3_start||'').toLowerCase().indexOf(arm3.toLowerCase())===-1) ok=false }
                    if(arm4){ if(String(r.arm4_start||'').toLowerCase().indexOf(arm4.toLowerCase())===-1) ok=false }

                    // Windows filters (substring match)
                    if(w1){ if(String(r.sequence_window1||'').toLowerCase().indexOf(w1.toLowerCase())===-1) ok=false }
                    if(w2){ if(String(r.sequence_window2||'').toLowerCase().indexOf(w2.toLowerCase())===-1) ok=false }
                    if(w3){ if(String(r.sequence_window3||'').toLowerCase().indexOf(w3.toLowerCase())===-1) ok=false }
                    if(w4){ if(String(r.sequence_window4||'').toLowerCase().indexOf(w4.toLowerCase())===-1) ok=false }
                    if(w5){ if(String(r.sequence_window5||'').toLowerCase().indexOf(w5.toLowerCase())===-1) ok=false }
                    if(w6){ if(String(r.sequence_window6||'').toLowerCase().indexOf(w6.toLowerCase())===-1) ok=false }

                    // Any-field free-text search
                    if(any){ var j = JSON.stringify(r).toLowerCase(); if(j.indexOf(any)===-1) ok=false }

                    if(ok) matches.push(r);
                }

                LAST_FILTERED = matches.slice();
                highlightFiltered(matches.map(function(x){return x.index}));
                populateResults(matches);
            }

            function highlightFiltered(indices){ var plot = getPlotDiv(); if(!plot) return; for(var t=0;t<plot.data.length;t++){ var data = plot.data[t]; var sizes = ORIGINAL_MARKER_SIZES[t] ? ORIGINAL_MARKER_SIZES[t].slice() : new Array((data.x && data.x.length)||0).fill(10); var colors = ORIGINAL_MARKER_COLORS[t] ? ORIGINAL_MARKER_COLORS[t].slice() : new Array((data.x && data.x.length)||0).fill('#808080'); var custom = data.customdata||[]; for(var j=0;j<custom.length;j++){ if(indices.indexOf(custom[j])!==-1){ sizes[j] = Math.max(sizes[j]||10, 18); colors[j] = '#FFD54F'; } } try{ Plotly.restyle(plot, {'marker.size':[sizes], 'marker.color':[colors]}, [t]); }catch(e){} } if(indices.length>0){ var first = indices[0]; try{ Plotly.relayout(plot, {'xaxis.range':[POINT_DATA[first].ref_energy-1, POINT_DATA[first].ref_energy+1], 'yaxis.range':[POINT_DATA[first].alt_energy-1, POINT_DATA[first].alt_energy+1]}); }catch(e){} } document.getElementById('resultsCount').innerText = indices.length }

            function populateResults(matches){ var html = '<table><thead><tr><th>#</th><th>Construct</th><th>SNP</th><th>Ref</th><th>Alt</th><th>Energy</th></tr></thead><tbody>'; for(var i=0;i<matches.length;i++){ var r=matches[i]; html += '<tr><td>'+r.index+'</td><td><a class="result-link" data-idx="'+r.index+'">'+r.construct_id+'</a></td><td>'+(r.snp===null?'':r.snp)+'</td><td>'+r.ref_energy.toFixed(3)+'</td><td>'+r.alt_energy.toFixed(3)+'</td><td>'+r.energy_type+'</td></tr>' } html += '</tbody></table>'; document.getElementById('resultsTable').innerHTML = html }

            function openModalByIndex(idx){ var rec = (window.POINT_DATA && window.POINT_DATA[idx]) ? window.POINT_DATA[idx] : null; if(rec) openModal(rec) }

            function downloadCSV(){ var rows = LAST_FILTERED; if(!rows || !rows.length){ alert('Нет данных для сохранения'); return } var keys = Object.keys(rows[0]); var csv = keys.join(',') + '\\n'; for(var i=0;i<rows.length;i++){ var vals = keys.map(function(k){ return '"'+String(rows[i][k]).replace(/"/g,'""')+'"' }); csv += vals.join(',') + '\\n' } var blob = new Blob([csv], {type:'text/csv;charset=utf-8;'}); var url = URL.createObjectURL(blob); var a = document.createElement('a'); a.href = url; a.download = 'filtered_points.csv'; document.body.appendChild(a); a.click(); document.body.removeChild(a); URL.revokeObjectURL(url) }

            function downloadJSON(){ var rows = LAST_FILTERED; if(!rows || !rows.length){ alert('Нет данных для сохранения'); return } var s = JSON.stringify(rows, null, 2); var blob = new Blob([s], {type:'application/json;charset=utf-8;'}); var url = URL.createObjectURL(blob); var a = document.createElement('a'); a.href = url; a.download = 'filtered_points.json'; document.body.appendChild(a); a.click(); document.body.removeChild(a); URL.revokeObjectURL(url) }

            function clearFilters(){ document.getElementById('filterConstruct').value=''; var sel=document.getElementById('filterSnpMulti'); if(sel){ for(var i=0;i<sel.options.length;i++){ sel.options[i].selected=false }} document.getElementById('filterMinEnergy').value=''; document.getElementById('filterMaxEnergy').value=''; document.getElementById('filterArm1').value=''; document.getElementById('filterArm2').value=''; document.getElementById('filterArm3').value=''; document.getElementById('filterArm4').value=''; document.getElementById('filterW1').value=''; document.getElementById('filterW2').value=''; document.getElementById('filterW3').value=''; document.getElementById('filterW4').value=''; document.getElementById('filterW5').value=''; document.getElementById('filterW6').value=''; document.getElementById('searchAny').value=''; resetAllMarkers(); }

            function openModal(data){ var modal = document.getElementById('dataModal'); var table = document.getElementById('dataTable'); var content = '<table style="width:100%; border-collapse:collapse">'; content += '<tr style="background:#f6f6f6"><th style="padding:8px;border:1px solid #eee">Поле</th><th style="padding:8px;border:1px solid #eee">Значение</th></tr>'; for(var k in data){ var v = data[k]; if(Array.isArray(v)){ v = v.join(', ') } content += '<tr><td style="padding:8px;border:1px solid #eee;font-weight:bold">'+k+'</td><td style="padding:8px;border:1px solid #eee">'+String(v)+'</td></tr>' } content += '</table>'; table.innerHTML = content; modal.style.display='block' }

            function closeModal(){ document.getElementById('dataModal').style.display='none' }
            window.onclick = function(event){ var modal = document.getElementById('dataModal'); if(event.target==modal) modal.style.display='none' }

            // Populate SNP multi-select with unique SNP values from POINT_DATA
            function populateSnpOptions(){ try{ var sel = document.getElementById('filterSnpMulti'); if(!sel) return; var unique = {}; for(var i=0;i<POINT_DATA.length;i++){ var s = POINT_DATA[i].snp; if(s!==null && s!==undefined && String(s)!=='') unique[String(s)] = true } var keys = Object.keys(unique).sort(function(a,b){return Number(a)-Number(b)}); for(var k=0;k<keys.length;k++){ var opt=document.createElement('option'); opt.value=keys[k]; opt.text=keys[k]; sel.appendChild(opt) } }catch(e){} }

            // Expose core functions to window to ensure they are callable from UI elements
            window.applyFilters = applyFilters;
            window.clearFilters = clearFilters;
            window.downloadCSV = downloadCSV;
            window.downloadJSON = downloadJSON;
            window.openModalByIndex = openModalByIndex;

            window.addEventListener('load', function(){ setTimeout(function(){ try{ storeOriginalMarkerState(); populateSnpOptions(); }catch(e){ console && console.error && console.error('init error', e); } 
                    // Attach event listeners for buttons (safer than inline onclick)
                    var btnTopFind = document.getElementById('btnTopFind'); if(btnTopFind) btnTopFind.addEventListener('click', applyFilters);
                    var btnTopClear = document.getElementById('btnTopClear'); if(btnTopClear) btnTopClear.addEventListener('click', clearFilters);
                    var btnApply = document.getElementById('btnApply'); if(btnApply) btnApply.addEventListener('click', applyFilters);
                    var btnClear = document.getElementById('btnClear'); if(btnClear) btnClear.addEventListener('click', clearFilters);
                    var btnDownloadCSV = document.getElementById('btnDownloadCSV'); if(btnDownloadCSV) btnDownloadCSV.addEventListener('click', downloadCSV);
                    var btnDownloadJSON = document.getElementById('btnDownloadJSON'); if(btnDownloadJSON) btnDownloadJSON.addEventListener('click', downloadJSON);
                    var btnSelectAllSnps = document.getElementById('btnSelectAllSnps'); if(btnSelectAllSnps) btnSelectAllSnps.addEventListener('click', selectAllSnps);
                    // Delegate clicks from results table
                    var resultsTable = document.getElementById('resultsTable'); if(resultsTable){ resultsTable.addEventListener('click', function(ev){ var target = ev.target || ev.srcElement; if(target && target.classList && target.classList.contains('result-link')){ var idx = target.getAttribute('data-idx'); if(idx) openModalByIndex(Number(idx)); } }); }
                    // Modal close
                    var modalClose = document.getElementById('modalClose'); if(modalClose) modalClose.addEventListener('click', closeModal);
                }, 300); });
            </script>
            """

            # Настроим hover шаблон (customdata мы уже задавали per-trace)
            # Use %{text} in hovertemplate and supply trace.text (Plotly expects 'text' variable)
            fig_p.update_traces(
                hovertemplate="%{text}<extra></extra>",
                hoverlabel=dict(bgcolor="white"),
            )

            fig_p.for_each_trace(
                lambda t: t.update(
                    hoverlabel=dict(bgcolor="white"),
                    hovertemplate="%{text}<extra></extra>"
                ) if getattr(t, 'mode', None) == "markers" else t
            )

            html_path_plotly = Path(output_dir) / f"{energy_type}_comparison_plotly.html"
            try:
                config = {
                    'displayModeBar': True,
                    'scrollZoom': True,
                    'modeBarButtonsToAdd': ['select2d', 'lasso2d'],
                }
                
                # Сохраняем с дополнительным JavaScript
                html_content = pio.to_html(
                    fig_p,
                    config=config,
                    include_plotlyjs='cdn',
                    full_html=True,
                    validate=True
                )
                
                # Вставляем наш JSON-объект с детальными данными и JavaScript/CSS для модального окна перед </body>
                try:
                    # embed as application/json to avoid JS parsing/syntax issues when sequences contain problematic substrings
                    json_script = f"<script type=\"application/json\" id=\"POINT_DATA_JSON\">{json.dumps(point_records, ensure_ascii=False)}</script>\n"
                except Exception:
                    # fallback: empty JSON block
                    json_script = "<script type=\"application/json\" id=\"POINT_DATA_JSON\">[]</script>\n"

                html_content = html_content.replace('</body>', f'{json_script}{modal_js}</body>')
                
                # Добавляем обработчик кликов
                click_handler = """
                <script>
                var plot = document.getElementsByClassName('plotly-graph-div')[0];
                if(plot) plot.on('plotly_click', function(data) {
                    var point = data.points[0];
                    if(point && point.customdata !== undefined){
                        var gidx = point.customdata; // global index
                        var rec = (window.POINT_DATA && window.POINT_DATA[gidx]) ? window.POINT_DATA[gidx] : null;
                        if(rec){ openModal(rec); }
                    }
                });
                </script>
                """
                html_content = html_content.replace('</body>', f'{click_handler}</body>')
                
                with open(html_path_plotly, 'w', encoding='utf-8') as f:
                    f.write(html_content)
                    
                logger.info(f"Plotly HTML с интерактивной легендой и модальным окном сохранён: {html_path_plotly}")
            except Exception as e:
                logger.warning(f"Не удалось сохранить Plotly HTML для {energy_type}: {e}")
    except Exception as e:
        logger.info(f"Plotly не доступен или произошла ошибка при создании HTML: {e}")
    
    return {
        'mean_diff': mean_diff,
        'std_diff': std_diff,
        'upper_outliers': np.sum(upper_outliers),
        'lower_outliers': np.sum(lower_outliers),
        'total_points': len(ref_data),
        'detailed_data_file': str(detailed_data_path),
        'search_data_file': str(search_data_path) if search_data else "Не создан"
        , 'point_records': point_records
    }

def find_individual_dirs(base_dir: str) -> List[Tuple[str, str]]:
    """Поиск директорий тестовых индивидуумов в базовой директории энергий."""
    individual_dirs: List[Tuple[str, str]] = []
    for d in os.listdir(base_dir):
        if d.startswith("SEQ-g38_Mt-Short_Test-test_individual_") and os.path.isdir(os.path.join(base_dir, d)):
            individual_id = d.split('-')[-1].split('_')[2]
            individual_dirs.append((os.path.join(base_dir, d), individual_id))
    return individual_dirs


def process_individual(ref_dir: str, alt_dir: str, snp_file_path: str, output_base_dir: str, individual_id: str) -> None:
    """Обработка данных и создание графиков для одного индивидуума."""
    logger.info(f"Обработка индивидуума: {individual_id}")
    snp_positions = load_snp_data(snp_file_path)
    output_dir = os.path.join(output_base_dir, f"test_individual_{individual_id}")
    os.makedirs(output_dir, exist_ok=True)

    # ИСПРАВЛЕНИЕ: Ищем CF файлы в директории конструктов для конкретного индивидуума
    constructs_base_dir = "D:/pythonProject/MitoFragility/MitoFragilityScore/Constructs"
    individual_constructs_dir = os.path.join(constructs_base_dir, f"SEQ-g38_Mt-Short_Test-test_individual_{individual_id}")
    
    logger.info(f"Поиск CF файлов в директории: {individual_constructs_dir}")
    cf_files = [f for f in os.listdir(individual_constructs_dir) if f.endswith('-CF.csv')]
    logger.info(f"Найдено CF файлов: {len(cf_files)}")
    if cf_files:
        logger.info(f"Примеры CF файлов: {cf_files[:3]}")  # Покажем первые 3 файла
    
    # Контейнеры для данных графиков с деталями конструктов
    energy_data = {
        'EnergyLeft': {'ref': [], 'alt': [], 'snp_value': [], 'construct_details': []},
        'EnergyRight': {'ref': [], 'alt': [], 'snp_value': [], 'construct_details': []},
        'Energy': {'ref': [], 'alt': [], 'snp_value': [], 'construct_details': []}
    }
    
    snp_counter = defaultdict(int)
    construct_categories = {'with_snp': 0, 'without_snp': 0, 'error': 0}
    
    all_snps = sorted(snp_positions)
    snp_colors = {snp: generate_distinct_colors(len(all_snps))[i] for i, snp in enumerate(all_snps)} if all_snps else {}
    
    # ЗАГРУЖАЕМ ДАННЫЕ КОНСТРУКТОВ ИЗ CF ФАЙЛОВ (устойчиво к названиям колонок)
    construct_data_map = {}
    cf_files_processed = 0

    desired_fields = {
        'construct_id': ['constructid', 'construct_id', 'construct id', 'id'],
        'center_position': ['centerposition', 'center_position', 'center position', 'center'],
        'arm1_start': ['arm1start', 'arm1_start', 'arm1'],
        'arm2_start': ['arm2start', 'arm2_start', 'arm2'],
        'arm3_start': ['arm3start', 'arm3_start', 'arm3'],
        'arm4_start': ['arm4start', 'arm4_start', 'arm4'],
        'sequence_window1': ['sequencewindow1', 'sequence_window1', 'sequence window1', 'sequence1', 'window1'],
        'sequence_window2': ['sequencewindow2', 'sequence_window2', 'sequence window2', 'sequence2', 'window2'],
        'sequence_window3': ['sequencewindow3', 'sequence_window3', 'sequence window3', 'sequence3', 'window3'],
        'sequence_window4': ['sequencewindow4', 'sequence_window4', 'sequence window4', 'sequence4', 'window4'],
        'sequence_window5': ['sequencewindow5', 'sequence_window5', 'sequence window5', 'sequence5', 'window5'],
        'sequence_window6': ['sequencewindow6', 'sequence_window6', 'sequence window6', 'sequence6', 'window6']
    }

    if os.path.isdir(individual_constructs_dir):
        for cf_file in os.listdir(individual_constructs_dir):
            if cf_file.endswith('-CF.csv'):
                try:
                    cf_path = os.path.join(individual_constructs_dir, cf_file)
                    cf_df = pd.read_csv(cf_path, dtype=str)

                    if cf_df.empty:
                        logger.debug(f"Пустой CF-файл пропущен: {cf_path}")
                        continue

                    # карта lower->original column name (строки обрезаем, т.к. иногда есть ведущие пробелы)
                    lower_map = {col.lower().strip(): col for col in cf_df.columns}

                    # найти колонку ConstructID
                    construct_id_col = None
                    for variant in desired_fields['construct_id']:
                        if variant in lower_map:
                            construct_id_col = lower_map[variant]
                            break
                    if construct_id_col is None:
                        logger.warning(f"CF-файл {cf_file} пропущен: не найден столбец ConstructID (проверялись варианты).")
                        continue

                    # для каждой строки собрать record с безопасными значениями
                    for _, row in cf_df.iterrows():
                        construct_id = row.get(construct_id_col)
                        if not construct_id or pd.isna(construct_id):
                            continue

                        record = {'construct_id': construct_id}
                        for key, variants in desired_fields.items():
                            if key == 'construct_id':
                                continue
                            col_name = None
                            for v in variants:
                                if v in lower_map:
                                    col_name = lower_map[v]
                                    break
                            record[key] = row.get(col_name, 'N/A') if col_name is not None else 'N/A'

                        construct_data_map[construct_id] = {
                            'center_position': record['center_position'],
                            'arm1_start': record['arm1_start'],
                            'arm2_start': record['arm2_start'],
                            'arm3_start': record['arm3_start'],
                            'arm4_start': record['arm4_start'],
                            'sequence_window1': record['sequence_window1'],
                            'sequence_window2': record['sequence_window2'],
                            'sequence_window3': record['sequence_window3'],
                            'sequence_window4': record['sequence_window4'],
                            'sequence_window5': record['sequence_window5'],
                            'sequence_window6': record['sequence_window6']
                        }
                    cf_files_processed += 1
                    logger.debug(f"Загружен конструкт файл: {cf_file} с {len(cf_df)} конструктами")
                except Exception as e:
                    logger.warning(f"Ошибка загрузки конструкт файла {cf_file}: {e}")
    else:
        logger.warning(f"Директория конструктов не найдена: {individual_constructs_dir}")
    
    logger.info(f"Загружено {cf_files_processed} CF файлов, всего конструктов в памяти: {len(construct_data_map)}")
    
    if construct_data_map:
        sample_constructs = list(construct_data_map.keys())[:5]
        logger.info(f"Примеры загруженных ConstructID: {sample_constructs}")

    # ОБРАБАТЫВАЕМ ENERGY ФАЙЛЫ (EF)
    for file in os.listdir(alt_dir):
        if not file.endswith('-EF.csv'):
            continue
            
        alt_path = os.path.join(alt_dir, file)
        ref_file = file.replace(f"-test_individual_{individual_id}", "")
        ref_path = os.path.join(ref_dir, ref_file)
        
        if not os.path.exists(ref_path):
            logger.warning(f"Референсный файл не найден: {ref_path}")
            continue
        
        try:
            ref_df = pd.read_csv(ref_path)
            alt_df = pd.read_csv(alt_path)
            
            # ПРОВЕРКА НАЛИЧИЯ НЕОБХОДИМЫХ СТОЛБЦОВ
            required_columns = ['ConstructID', 'EnergyLeft', 'EnergyRight', 'Energy']
            missing_in_ref = [col for col in required_columns if col not in ref_df.columns]
            missing_in_alt = [col for col in required_columns if col not in alt_df.columns]
            
            if missing_in_ref:
                logger.warning(f"В референсном файле {ref_path} отсутствуют столбцы: {missing_in_ref}. Файл пропущен.")
                continue
            if missing_in_alt:
                logger.warning(f"В альтернативном файле {alt_path} отсутствуют столбцы: {missing_in_alt}. Файл пропущен.")
                continue
            
            for i in range(min(len(ref_df), len(alt_df))):
                construct_id = alt_df.iloc[i]['ConstructID']
                
                # ПОЛУЧАЕМ ДАННЫЕ КОНСТРУКТА ИЗ ПРЕДВАРИТЕЛЬНО ЗАГРУЖЕННОЙ MAP
                construct_details_from_cf = construct_data_map.get(construct_id, {})
                
                if construct_details_from_cf:
                    # Есть данные из CF файла
                    construct_details = {
                        'construct_id': construct_id,
                        'center_position': construct_details_from_cf['center_position'],
                        'arm1_start': construct_details_from_cf['arm1_start'],
                        'arm2_start': construct_details_from_cf['arm2_start'],
                        'arm3_start': construct_details_from_cf['arm3_start'],
                        'arm4_start': construct_details_from_cf['arm4_start'],
                        'sequence_window1': construct_details_from_cf['sequence_window1'],
                        'sequence_window2': construct_details_from_cf['sequence_window2'],
                        'sequence_window3': construct_details_from_cf['sequence_window3'],
                        'sequence_window4': construct_details_from_cf['sequence_window4'],
                        'sequence_window5': construct_details_from_cf['sequence_window5'],
                        'sequence_window6': construct_details_from_cf['sequence_window6']
                    }
                else:
                    logger.debug(f"ConstructID {construct_id} не найден в CF данных. Доступные ConstructID: {len(construct_data_map)}")
                    # Данных из CF файла нет, используем только то, что есть в EF файле
                    logger.warning(f"Данные конструкта не найдены для: {construct_id}")
                    construct_details = {
                        'construct_id': construct_id,
                        'center_position': 'N/A',
                        'arm1_start': 'N/A', 
                        'arm2_start': 'N/A',
                        'arm3_start': 'N/A',
                        'arm4_start': 'N/A',
                        'sequence_window1': 'N/A',
                        'sequence_window2': 'N/A',
                        'sequence_window3': 'N/A',
                        'sequence_window4': 'N/A',
                        'sequence_window5': 'N/A',
                        'sequence_window6': 'N/A'
                    }
                
                # Парсинг ID конструкта и проверка SNP
                arm_size, center, arm3_start, arm4_start = parse_construct_id(construct_id)
                if arm_size is None:
                    construct_categories['error'] += 1
                    logger.warning(f"Ошибка парсинга для конструкта: {construct_id}")
                    selected_snp = None
                else:
                    # Проверка SNP в плечах конструкта
                    snps = get_snps_in_construct(snp_positions, arm_size, center, arm3_start, arm4_start)
                    if snps:
                        construct_categories['with_snp'] += 1
                        selected_snp = min(snps)
                        snp_counter[selected_snp] += 1
                    else:
                        construct_categories['without_snp'] += 1
                        selected_snp = None
                
                # Сбор данных для построения графиков
                for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
                    ref_energy = ref_df.iloc[i][energy_type]
                    alt_energy = alt_df.iloc[i][energy_type]
                    
                    if not np.isnan(ref_energy) and not np.isnan(alt_energy):
                        energy_data[energy_type]['ref'].append(ref_energy)
                        energy_data[energy_type]['alt'].append(alt_energy)
                        energy_data[energy_type]['snp_value'].append(selected_snp)
                        energy_data[energy_type]['construct_details'].append(construct_details)
        
        except Exception as e:
            logger.error(f"Ошибка обработки файла {file}: {e}")
            construct_categories['error'] += 1
    
    total_constructs = sum(construct_categories.values())
    logger.info(f"Категории конструктов для {individual_id}: всего {total_constructs}, с SNP {construct_categories['with_snp']}, без SNP {construct_categories['without_snp']}, с ошибкой {construct_categories['error']}")
    
    sorted_snps = sorted(snp_counter.items(), key=lambda x: x[1], reverse=True)[:10]
    logger.info("Топ-10 самых частых SNP:")
    for snp, count in sorted_snps:
        logger.info(f"  SNP {snp}: {count} конструктов")
    
    # Генерация и сохранение графиков
    outliers_stats = {}
    master_records = []
    for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
        ref_vals = energy_data[energy_type]['ref']
        alt_vals = energy_data[energy_type]['alt']
        snp_vals = energy_data[energy_type]['snp_value']
        construct_details = energy_data[energy_type]['construct_details']
        
        if ref_vals and alt_vals:
            stats = plot_energy_comparison(
                ref_vals, alt_vals, snp_vals, snp_colors, 
                energy_type, output_dir, individual_id, construct_details
            )
            outliers_stats[energy_type] = stats
            # Сбор master-списка детальных записей для одного индивидуума
            if stats and isinstance(stats, dict) and stats.get('point_records'):
                # Включаем поле energy_type в каждую запись уже done в plot_energy_comparison
                master_records.extend(stats.get('point_records'))
        else:
            logger.warning(f"Нет данных для {energy_type}")
    
    # Запишем master JSON для всего индивидуума (все energy-type вместе)
    try:
        master_json_path = Path(output_dir) / f"{individual_id}_master_data.json"
        with open(master_json_path, 'w', encoding='utf-8') as mf:
            json.dump(master_records, mf, ensure_ascii=False, indent=2)
        logger.info(f"Master JSON с детальными данными сохранён: {master_json_path}")
    except Exception as e:
        logger.warning(f"Не удалось сохранить master JSON для {individual_id}: {e}")

    write_outlier_stats(outliers_stats, output_dir, individual_id)

def write_outlier_stats(outliers_stats: Dict[str, Dict[str, Any]], output_dir: str, individual_id: str) -> None:
    """Запись статистики выбросов в текстовый файл."""
    stats_path = Path(output_dir) / f"{individual_id}_outliers_statistics.txt"
    with open(stats_path, 'w', encoding='utf-8') as f:
        f.write("Статистика выбросов по типам энергии:\n")
        f.write("==================================================\n")
        for energy_type, stats in outliers_stats.items():
            f.write(f"{energy_type}:\n")
            f.write(f"  Всего точек: {stats['total_points']}\n")
            f.write(f"  Средняя разница (ref - alt): {stats['mean_diff']:.4f}\n")
            f.write(f"  Стандартное отклонение: {stats['std_diff']:.4f}\n")
            f.write(f"  Верхние выбросы (> +2std): {stats['upper_outliers']} ({stats['upper_outliers']/stats['total_points']*100:.2f}%)\n")
            f.write(f"  Нижние выбросы (< -2std): {stats['lower_outliers']} ({stats['lower_outliers']/stats['total_points']*100:.2f}%)\n")
            f.write("\n")
    logger.info(f"Статистика сохранена: {stats_path}")


def search_sequences_wrapper():
    """Утилита поиска по последовательностям через командную строку."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Поиск точек по последовательности')
    parser.add_argument('search_file', help='CSV файл с данными для поиска')
    parser.add_argument('sequence', help='Последовательность для поиска')
    parser.add_argument('--case-sensitive', action='store_true', help='Поиск с учетом регистра')
    
    args = parser.parse_args()
    
    results = find_points_by_sequence(args.search_file, args.sequence, args.case_sensitive)
    
    if results:
        print(f"Найдено {len(results)} точек с последовательностью '{args.sequence}':")
        for result in results:
            print(f"  Точка {result['point_index']}: {result['construct_id']}")
            print(f"    Энергии: Ref={result['ref_energy']:.2f}, Alt={result['alt_energy']:.2f}")
            print(f"    Найдено в: {result['window_found']}")
            print(f"    SNP: {result['snp']}")
            print()
    else:
        print(f"Последовательность '{args.sequence}' не найдена")


def main(base_dir: Optional[str] = None, output_base_dir: Optional[str] = None, 
         snp_base_dir: Optional[str] = None, search_mode: bool = False) -> None:
    """Основная функция с поддержкой режима поиска."""
    
    if search_mode:
        search_sequences_wrapper()
        return

    if base_dir is None:
        base_dir = "D:/pythonProject/MitoFragility/MitoFragilityScore/Energies"
    if output_base_dir is None:
        output_base_dir = "D:/pythonProject/MitoFragility/DataPreparing/output"
    if snp_base_dir is None:
        snp_base_dir = "D:/pythonProject/MitoFragility/MitoFragilityScore/Sequences/Relative"

    ref_dir = os.path.join(base_dir, "SEQ-g38_Mt-Short_Test")
    os.makedirs(output_base_dir, exist_ok=True)
    individual_dirs = find_individual_dirs(base_dir)
    if not individual_dirs:
        logger.warning("Не найдено ни одной директории теста!")
        return
    logger.info(f"Найдено директорий теста: {len(individual_dirs)}")
    for alt_dir, individual_id in individual_dirs:
        snp_file_path = os.path.join(snp_base_dir, f"test_individual_{individual_id}.csv")
        if not os.path.exists(alt_dir):
            logger.warning(f"Директория теста не существует: {alt_dir}")
            continue
        if not os.path.exists(snp_file_path):
            logger.warning(f"Файл SNP не найден: {snp_file_path}")
            continue
        try:
            process_individual(
                ref_dir,
                alt_dir,
                snp_file_path,
                output_base_dir,
                individual_id
            )
        except Exception as e:
            logger.error(f"Ошибка при обработке теста {individual_id}: {str(e)}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "search":
        search_sequences_wrapper()
    else:
        main()