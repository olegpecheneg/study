import os
import re
import logging
from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import colorsys


def setup_matplotlib() -> None:
    """Configure matplotlib backend and default font."""
    matplotlib.use('Agg')
    plt.rcParams['font.family'] = 'DejaVu Sans'


def setup_logger(name: str = "scatter_plus_n_std") -> logging.Logger:
    """Configures and returns a module logger.

    This mirrors the style used in `mt_DNA_builder.py` so the module can be
    imported dynamically and provide a predictable logger.
    """
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        log_path = Path(__file__).resolve().parent.parent / 'visualization.log'
        file_handler = logging.FileHandler(log_path, encoding='utf-8')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(stream_handler)
        logger.addHandler(file_handler)
        logger.setLevel(logging.INFO)
    return logger


setup_matplotlib()
logger = setup_logger()


def load_snp_data(snp_file_path: str) -> set:
    """Load SNP positions from a CSV-like file (first column: position).

    Returns a set of integer positions.
    """
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
    """Parse construct id and return (arm_size, center, arm3_start, arm4_start).

    Returns None for any value that cannot be parsed.
    """
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
    """Compute absolute ranges for the four arms of a construct.

    The function converts relative arm coordinates into absolute positions by
    adding `subseq_start` which represents the offset used in the dataset.
    """
    arm1_start = center - arm_size
    arm2_start = center
    ranges = [
        (arm1_start + subseq_start, arm1_start + arm_size + subseq_start),  # arm1 absolute
        (arm2_start + subseq_start, arm2_start + arm_size + subseq_start),  # arm2 absolute
        (arm3_start + subseq_start, arm3_start + arm_size + subseq_start),  # arm3 absolute
        (arm4_start + subseq_start, arm4_start + arm_size + subseq_start)   # arm4 absolute
    ]
    logger.debug(f"Construct ranges (absolute): {ranges}")
    return ranges

def get_snps_in_construct(snp_positions: set, arm_size: int, center: int, arm3_start: int, arm4_start: int, subseq_start: int = 10) -> List[int]:
    """Return SNPs that fall within any of the construct's arms.

    Uses the absolute ranges returned by `calculate_arm_ranges`.
    """
    ranges = calculate_arm_ranges(arm_size, center, arm3_start, arm4_start, subseq_start)
    snps_in_construct = []
    for snp in snp_positions:
        for start, end in ranges:
            if start <= snp <= end:
                snps_in_construct.append(snp)
                break
    logger.debug(f"SNPs in construct: {snps_in_construct}")
    return snps_in_construct

def generate_distinct_colors(n: int) -> List[Tuple[float, float, float]]:
    """Generate a list of visually distinct RGB colors.

    Uses the golden-ratio conjugate to space hues evenly in HSV space.
    """
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        saturation = 0.9
        value = 0.9
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append(rgb)
    return colors

def calculate_outlier_stats(ref_data: List[float], alt_data: List[float]) -> Tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    """Compute outlier statistics comparing reference and alternative arrays.

    Returns mean difference, std deviation and boolean masks for upper, lower
    and normal points (upper/lower defined as > +/- 2 std).
    """
    diff = np.array(ref_data) - np.array(alt_data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    upper_outliers = diff > mean_diff + 2 * std_diff
    lower_outliers = diff < mean_diff - 2 * std_diff
    normal_points = ~(upper_outliers | lower_outliers)
    return mean_diff, std_diff, upper_outliers, lower_outliers, normal_points

def plot_scatter_points(ax: Any, ref_data: List[float], alt_data: List[float], snp_values: List[Optional[int]], snp_colors: Dict[int, Tuple[float, float, float]], upper_outliers: np.ndarray, lower_outliers: np.ndarray, normal_points: np.ndarray) -> None:
    """Plot scatter points colored by SNP and highlight outliers.

    Normal points are plotted smaller; outliers are emphasized with larger
    markers and colored edges.
    """
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
    """Add a dashed diagonal reference line (y == x)."""
    ax.plot([min_e, max_e], [min_e, max_e], 'k--', linewidth=2, alpha=0.7)

def add_outlier_zones(ax: Any, x: np.ndarray, mean_diff: float, std_diff: float, min_e: float, max_e: float) -> None:
    """Shade regions corresponding to outlier thresholds (+/- 2,3,4 std).

    Red zones mark upper outliers and green zones mark lower outliers.
    """
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
    """Create legend elements for the plot from SNP color mapping."""
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

def plot_energy_comparison(ref_data: List[float], alt_data: List[float], snp_values: List[Optional[int]], snp_colors: Dict[int, Tuple[float, float, float]], energy_type: str, output_dir: str, individual_id: str) -> Optional[Dict[str, Any]]:
    """Create and save a scatterplot comparing reference and alternative energies.

    Points are colored by SNP, and outlier zones are shaded. Returns a small
    dictionary with outlier statistics for logging and file output.
    """
    if not ref_data or not alt_data:
        logger.warning(f"Нет данных для построения графика {energy_type}")
        return

    min_len = min(len(ref_data), len(alt_data), len(snp_values))
    ref_data = ref_data[:min_len]
    alt_data = alt_data[:min_len]
    snp_values = snp_values[:min_len]
    
    mean_diff, std_diff, upper_outliers, lower_outliers, normal_points = calculate_outlier_stats(ref_data, alt_data)
    
    logger.info(f"Для {energy_type}:")
    logger.info(f"  Средняя разница: {mean_diff:.2f}, Стандартное отклонение: {std_diff:.2f}")
    logger.info(f"  Верхние выбросы (> +2std): {np.sum(upper_outliers)} точек")
    logger.info(f"  Нижние выбросы (< -2std): {np.sum(lower_outliers)} точек")

    fig, ax = plt.subplots(figsize=(16, 12))
    
    plot_scatter_points(ax, ref_data, alt_data, snp_values, snp_colors, upper_outliers, lower_outliers, normal_points)
    
    min_e = min(min(ref_data), min(alt_data))
    max_e = max(max(ref_data), max(alt_data))
    add_diagonal_line(ax, min_e, max_e)
    
    x = np.linspace(min_e, max_e, 100)
    add_outlier_zones(ax, x, mean_diff, std_diff, min_e, max_e)

    plt.title(f'Comparison of {energy_type} for individual {individual_id}', fontsize=18)
    plt.xlabel('Reference energy (kcal/mol)', fontsize=16)
    plt.ylabel('Alternative energy (kcal/mol)', fontsize=16)
    
    legend_elements = create_legend_elements(snp_colors)
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, title="Legend", title_fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.2)
    ax.set_xlim(min(ref_data), max(ref_data))
    ax.set_ylim(min(alt_data), max(alt_data))
    plt.tight_layout()

    output_path = Path(output_dir) / f"{energy_type}_comparison.png"
    plt.savefig(output_path, dpi=250, bbox_inches='tight')
    plt.close()
    logger.info(f"График сохранён: {output_path}")
    
    return {
        'mean_diff': mean_diff,
        'std_diff': std_diff,
        'upper_outliers': np.sum(upper_outliers),
        'lower_outliers': np.sum(lower_outliers),
        'total_points': len(ref_data)
    }

def find_individual_dirs(base_dir: str) -> List[Tuple[str, str]]:
    """Find test individual directories in energies base directory."""
    individual_dirs: List[Tuple[str, str]] = []
    for d in os.listdir(base_dir):
        if d.startswith("SEQ-g38_Mt-Short_Test-test_individual_") and os.path.isdir(os.path.join(base_dir, d)):
            individual_id = d.split('-')[-1].split('_')[2]  # Извлекаем номер после 'test_individual_'
            individual_dirs.append((os.path.join(base_dir, d), individual_id))
    return individual_dirs

def process_individual(ref_dir: str, alt_dir: str, snp_file_path: str, output_base_dir: str, individual_id: str) -> None:
    """Process data and produce plots for a single individual."""
    logger.info(f"Обработка индивидуума: {individual_id}")
    # Load SNP positions
    snp_positions = load_snp_data(snp_file_path)
    # Create output directory
    output_dir = os.path.join(output_base_dir, f"test_individual_{individual_id}")
    os.makedirs(output_dir, exist_ok=True)
    # Initialize data containers for plots
    energy_data = {
        'EnergyLeft': {'ref': [], 'alt': [], 'snp_value': []},
        'EnergyRight': {'ref': [], 'alt': [], 'snp_value': []},
        'Energy': {'ref': [], 'alt': [], 'snp_value': []}
    }
    snp_counter = defaultdict(int)
    construct_categories = {'with_snp': 0, 'without_snp': 0, 'error': 0}
    
    # Generate colors for SNPs
    all_snps = sorted(snp_positions)
    snp_colors = {snp: generate_distinct_colors(len(all_snps))[i] for i, snp in enumerate(all_snps)} if all_snps else {}
    # Process energy files
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
            
            for i in range(min(len(ref_df), len(alt_df))):
                construct_id = alt_df.iloc[i]['ConstructID']
                
                # Parse construct ID and check for SNPs
                arm_size, center, arm3_start, arm4_start = parse_construct_id(construct_id)
                if arm_size is None:
                    construct_categories['error'] += 1
                    logger.warning(f"Ошибка парсинга для конструкта: {construct_id}")
                    selected_snp = None
                else:
                    # Check for SNPs in construct arms
                    snps = get_snps_in_construct(snp_positions, arm_size, center, arm3_start, arm4_start)
                    if snps:
                        construct_categories['with_snp'] += 1
                        selected_snp = min(snps)  # Выбираем минимальный SNP (из второй версии)
                        snp_counter[selected_snp] += 1
                    else:
                        construct_categories['without_snp'] += 1
                        selected_snp = None
                
                # Collect data for plotting
                for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
                    ref_energy = ref_df.iloc[i][energy_type]
                    alt_energy = alt_df.iloc[i][energy_type]
                    
                    if not np.isnan(ref_energy) and not np.isnan(alt_energy):
                        energy_data[energy_type]['ref'].append(ref_energy)
                        energy_data[energy_type]['alt'].append(alt_energy)
                        energy_data[energy_type]['snp_value'].append(selected_snp)
        
        except Exception as e:
            logger.error(f"Ошибка обработки файла {file}: {e}")
            construct_categories['error'] += 1
    
    # Log construct categories
    total_constructs = sum(construct_categories.values())
    logger.info(f"Категории конструктов для {individual_id}: всего {total_constructs}, с SNP {construct_categories['with_snp']}, без SNP {construct_categories['without_snp']}, с ошибкой {construct_categories['error']}")
    # Log top SNPs
    sorted_snps = sorted(snp_counter.items(), key=lambda x: x[1], reverse=True)[:10]
    logger.info("Топ-10 самых частых SNP:")
    for snp, count in sorted_snps:
        logger.info(f"  SNP {snp}: {count} конструктов")
    
    # Generate and save plots
    outliers_stats = {}
    for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
        ref_vals = energy_data[energy_type]['ref']
        alt_vals = energy_data[energy_type]['alt']
        snp_vals = energy_data[energy_type]['snp_value']
        
        if ref_vals and alt_vals:
            stats = plot_energy_comparison(
                ref_vals, alt_vals, snp_vals, snp_colors, 
                energy_type, output_dir, individual_id
            )
            outliers_stats[energy_type] = stats
        else:
            logger.warning(f"Нет данных для {energy_type}")
    
    # Save outlier statistics to file
    write_outlier_stats(outliers_stats, output_dir, individual_id)

def write_outlier_stats(outliers_stats: Dict[str, Dict[str, Any]], output_dir: str, individual_id: str) -> None:
    """Write outlier statistics to a text file in the output directory."""
    stats_path = Path(output_dir) / f"{individual_id}_outliers_statistics.txt"
    with open(stats_path, 'w', encoding='utf-8') as f:
        f.write("Outlier statistics by energy type:\n")
        f.write("==================================================\n")
        for energy_type, stats in outliers_stats.items():
            f.write(f"{energy_type}:\n")
            f.write(f"  Total points: {stats['total_points']}\n")
            f.write(f"  Mean difference (ref - alt): {stats['mean_diff']:.4f}\n")
            f.write(f"  Std deviation: {stats['std_diff']:.4f}\n")
            f.write(f"  Upper outliers (> +2std): {stats['upper_outliers']} ({stats['upper_outliers']/stats['total_points']*100:.2f}%)\n")
            f.write(f"  Lower outliers (< -2std): {stats['lower_outliers']} ({stats['lower_outliers']/stats['total_points']*100:.2f}%)\n")
            f.write("\n")
    logger.info(f"Статистика сохранена: {stats_path}")


def main(base_dir: Optional[str] = None, output_base_dir: Optional[str] = None, snp_base_dir: Optional[str] = None) -> None:
    """Main entry similar to `mt_DNA_builder` style but configurable for imports/tests.

    Paths default to the previous hard-coded values but can be overridden for tests.
    """
    if base_dir is None:
        base_dir = "D:/pythonProject/MitoFragility/MitoFragilityScore/Energies"
    if output_base_dir is None:
        output_base_dir = "D:/pythonProject/MitoFragility/DataPreparing/plots/output"
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

    main()
