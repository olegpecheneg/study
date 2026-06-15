from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
from logging import getLogger
import json
import numpy as np
import matplotlib.pyplot as plt
from .plot_helpers import plot_scatter_points, add_diagonal_line, add_outlier_zones, create_legend_elements
from .io import save_detailed_data, save_search_csv

logger = getLogger(__name__)


def plot_energy_comparison(ref_data: List[float], alt_data: List[float],
                           snp_values: List[Optional[int]],
                           snp_colors: Dict[int, Tuple[float, float, float]],
                           energy_type: str, output_dir: str, individual_id: str,
                           construct_details: List[Dict]) -> Optional[Dict[str, Any]]:
    """Создание scatter plot и экспорт данных (PNG + optional HTML).

    Возвращает словарь со статистикой и point_records.
    """
    if not ref_data or not alt_data:
        logger.warning(f"Нет данных для построения графика {energy_type}")
        return None

    min_len = min(len(ref_data), len(alt_data), len(snp_values), len(construct_details))
    ref_data = ref_data[:min_len]
    alt_data = alt_data[:min_len]
    snp_values = snp_values[:min_len]
    construct_details = construct_details[:min_len]

    mean_diff, std_diff, upper_outliers, lower_outliers, normal_points = __calculate_outliers(ref_data, alt_data)

    # Сохранение детальных данных (упрощённо через io.save_detailed_data)
    try:
        detailed_file = save_detailed_data(output_dir, energy_type, individual_id, construct_details, ref_data, alt_data, snp_values, upper_outliers, lower_outliers)
    except Exception as e:
        logger.warning(f"Ошибка при сохранении детальных данных: {e}")
        detailed_file = None

    # Создание CSV файла для поиска
    search_data = []
    for i in range(min_len):
        details = construct_details[i]
        if details.get('sequence_window1', 'N/A') != 'N/A':
            search_data.append({
                'point_index': i,
                'ref_energy': ref_data[i],
                'alt_energy': alt_data[i],
                'construct_id': details.get('construct_id'),
                'snp': snp_values[i],
                'window1': details.get('sequence_window1'),
                'window2': details.get('sequence_window2'),
                'window3': details.get('sequence_window3'),
                'window4': details.get('sequence_window4'),
                'window5': details.get('sequence_window5'),
                'window6': details.get('sequence_window6')
            })
    try:
        search_file = save_search_csv(output_dir, energy_type, search_data)
    except Exception as e:
        logger.warning(f"Ошибка при сохранении search CSV: {e}")
        search_file = None

    # Построение matplotlib графика
    fig, ax = plt.subplots(figsize=(16, 12))
    plot_scatter_points(ax, ref_data, alt_data, snp_values, snp_colors, upper_outliers, lower_outliers, normal_points)

    min_e = min(min(ref_data), min(alt_data))
    max_e = max(max(ref_data), max(alt_data))
    add_diagonal_line(ax, min_e, max_e)
    x = np.linspace(min_e, max_e, 100)
    add_outlier_zones(ax, x, mean_diff, std_diff, min_e, max_e)

    plt.title(f'Сравнение {energy_type} для индивидуума {individual_id}\nФайл данных: {Path(detailed_file).name if detailed_file else "N/A"}', fontsize=16)
    plt.xlabel('Референсная энергия (kcal/mol)', fontsize=16)
    plt.ylabel('Альтернативная энергия (kcal/mol)', fontsize=16)
    legend_elements = create_legend_elements(snp_colors)
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, title="Легенда", title_fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.2)
    ax.set_xlim(min(ref_data), max(ref_data))
    ax.set_ylim(min(alt_data), max(alt_data))
    plt.tight_layout()

    # Сохранение PNG
    output_path = Path(output_dir) / f"{energy_type}_comparison.png"
    plt.savefig(output_path, dpi=250, bbox_inches='tight')
    plt.close()
    logger.info(f"График сохранён: {output_path}")

    # Подготовка point_records для внешнего использования
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
            'detailed_file': Path(detailed_file).name if detailed_file else ''
        })

    return {
        'mean_diff': mean_diff,
        'std_diff': std_diff,
        'upper_outliers': int(np.sum(upper_outliers)),
        'lower_outliers': int(np.sum(lower_outliers)),
        'total_points': len(ref_data),
        'detailed_data_file': str(detailed_file) if detailed_file else '',
        'search_data_file': str(search_file) if search_file else '',
        'point_records': point_records
    }


def __calculate_outliers(ref_data: List[float], alt_data: List[float]):
    diff = np.array(ref_data) - np.array(alt_data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    upper_outliers = diff > mean_diff + 2 * std_diff
    lower_outliers = diff < mean_diff - 2 * std_diff
    normal_points = ~(upper_outliers | lower_outliers)
    return mean_diff, std_diff, upper_outliers, lower_outliers, normal_points
