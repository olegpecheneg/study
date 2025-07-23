import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pathlib import Path
import matplotlib
import logging
import re
from collections import defaultdict
from matplotlib.lines import Line2D
import colorsys

# Настройка бэкенда и шрифта
matplotlib.use('Agg')
plt.rcParams['font.family'] = 'DejaVu Sans'

# Настройка логирования с кодировкой UTF-8
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
file_handler = logging.FileHandler('visualization.log', encoding='utf-8')
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logging.basicConfig(level=logging.INFO, handlers=[stream_handler, file_handler])
logger = logging.getLogger(__name__)

def load_snp_data(snp_file_path):
    """
    Загружает данные о SNP из файла.
    Возвращает множество позиций SNP.
    """
    snp_positions = set()
    try:
        with open(snp_file_path, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if parts:
                    try:
                        position = int(parts[0])
                        snp_positions.add(position)
                    except (ValueError, IndexError):
                        logger.warning(f"Ошибка обработки строки SNP: {line.strip()}")
        logger.info(f"Загружено {len(snp_positions)} SNP из файла {snp_file_path}")
    except Exception as e:
        logger.error(f"Ошибка загрузки файла SNP {snp_file_path}: {e}")
    return snp_positions

def parse_construct_id(construct_id):
    """
    Извлекает параметры конструкта из ID.
    Возвращает arm_size, center, arm3_start, arm4_start.
    """
    try:
        cgs_match = re.search(r'CGS-(\d+)-(\d+)-(\d+)-(\d+)-(\d+)-(\d+)', construct_id)
        if cgs_match:
            arm_size = int(cgs_match.group(6))
        else:
            logger.warning(f"Не найден блок CGS в ID: {construct_id}")
            return None, None, None, None
        cen_match = re.search(r'CEN-(\d+)', construct_id)
        if cen_match:
            center = int(cen_match.group(1))
        else:
            logger.warning(f"Не найден блок CEN в ID: {construct_id}")
            return None, None, None, None
        con_match = re.search(r'CON-(\d+)-(\d+)', construct_id)
        if con_match:
            arm3_start = int(con_match.group(1))
            arm4_start = int(con_match.group(2))
        else:
            logger.warning(f"Не найден блок CON в ID: {construct_id}")
            return None, None, None, None
        return arm_size, center, arm3_start, arm4_start
    except Exception as e:
        logger.error(f"Ошибка парсинга ID конструкта {construct_id}: {e}")
        return None, None, None, None

def calculate_arm_ranges(arm_size, center, arm3_start, arm4_start):
    """
    Вычисляет диапазоны для всех четырёх плеч.
    """
    arm1_start = center - arm_size
    arm1_end = center
    arm2_start = center
    arm2_end = center + arm_size
    arm3_end = arm3_start + arm_size
    arm4_end = arm4_start + arm_size
    return [
        (arm1_start, arm1_end),
        (arm2_start, arm2_end),
        (arm3_start, arm3_end),
        (arm4_start, arm4_end)
    ]

def get_snps_in_construct(construct_id, snp_positions):
    """
    Возвращает список SNP, присутствующих в конструкте.
    """
    arm_size, center, arm3_start, arm4_start = parse_construct_id(construct_id)
    if None in (arm_size, center, arm3_start, arm4_start):
        return []
    arm_ranges = calculate_arm_ranges(arm_size, center, arm3_start, arm4_start)
    snps_in_construct = []
    for snp in snp_positions:
        for start, end in arm_ranges:
            if start <= snp <= end:
                snps_in_construct.append(snp)
                break
    return snps_in_construct

def generate_distinct_colors(n):
    """
    Генерирует набор максимально различимых цветов.
    """
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        saturation = 0.9
        value = 0.9
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append(rgb)
    return colors

def calculate_outlier_stats(ref_data, alt_data):
    """
    Рассчитывает статистику выбросов.
    """
    diff = np.array(ref_data) - np.array(alt_data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    upper_outliers = diff > mean_diff + 2 * std_diff
    lower_outliers = diff < mean_diff - 2 * std_diff
    normal_points = ~(upper_outliers | lower_outliers)
    return mean_diff, std_diff, upper_outliers, lower_outliers, normal_points

def plot_scatter_points(ax, ref_data, alt_data, snp_values, snp_colors, upper_outliers, lower_outliers, normal_points):
    """
    Рисует точки на графике.
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

def add_diagonal_line(ax, min_e, max_e):
    """
    Добавляет диагональную линию.
    """
    ax.plot([min_e, max_e], [min_e, max_e], 'k--', linewidth=2, alpha=0.7)

def add_outlier_zones(ax, x, mean_diff, std_diff, min_e, max_e):
    """
    Добавляет заполненные зоны для выбросов.
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

def create_legend_elements(snp_colors):
    """
    Создаёт элементы легенды.
    """
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

def plot_energy_comparison(ref_data, alt_data, snp_values, snp_colors, energy_type, output_dir, individual_id):
    """
    Строит scatterplot с раскраской точек по конкретным SNP и выделением выбросов.
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

    plt.title(f'Сравнение {energy_type} с выделением выбросов', fontsize=18)
    plt.xlabel('Референсная энергия (ккал/моль)', fontsize=16)
    plt.ylabel('Альтернативная энергия (ккал/моль)', fontsize=16)
    
    legend_elements = create_legend_elements(snp_colors)
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, title="Легенда", title_fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.2)
    ax.set_xlim(min(ref_data), max(ref_data))
    ax.set_ylim(min(alt_data), max(alt_data))
    plt.tight_layout()

    # Изменяем название файла для включения ID теста
    output_path = Path(output_dir) / f"test_individual_{individual_id}_{energy_type}_snp_outliers.png"
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

def find_individual_dirs(base_dir):
    """Находит все поддиректории, названия которых оканчиваются на цифру"""
    individual_dirs = []
    for entry in os.listdir(base_dir):
        entry_path = os.path.join(base_dir, entry)
        if os.path.isdir(entry_path) and re.search(r'\d$', entry):
            # Извлекаем ID теста из названия директории
            match = re.search(r'test_individual_(\d+)$', entry)
            if match:
                individual_id = match.group(1)
                individual_dirs.append((entry_path, individual_id))
            else:
                logger.warning(f"Не удалось извлечь ID из названия директории: {entry}")
    return individual_dirs

def initialize_data():
    """
    Инициализирует структуры данных для хранения энергий и статистики.
    """
    energy_data = {
        'EnergyLeft': {'ref': [], 'alt': [], 'snp_value': []},
        'EnergyRight': {'ref': [], 'alt': [], 'snp_value': []},
        'Energy': {'ref': [], 'alt': [], 'snp_value': []}
    }
    snp_counter = defaultdict(int)
    return energy_data, snp_counter

def process_individual(ref_dir, alt_dir, snp_file_path, output_dir, individual_id):
    """
    Обрабатывает данные для одного теста
    """
    logger.info(f"Обработка теста с ID: {individual_id}")
    logger.info(f"Директория теста: {alt_dir}")
    logger.info(f"Файл SNP: {snp_file_path}")
    
    snp_positions = set()
    if snp_file_path and os.path.exists(snp_file_path):
        snp_positions = load_snp_data(snp_file_path)
    else:
        logger.warning(f"Файл SNP не найден: {snp_file_path}. Все точки будут серыми.")
    
    energy_data, snp_counter = initialize_data()
    all_snps = sorted(snp_positions)
    snp_colors = {snp: generate_distinct_colors(len(all_snps))[i] for i, snp in enumerate(all_snps)} if all_snps else {}
    
    total_constructs, snp_constructs, error_constructs = process_constructs(
        ref_dir, alt_dir, snp_positions, energy_data, snp_counter, individual_id
    )
    log_statistics(total_constructs, snp_constructs, error_constructs, snp_counter)
    
    outliers_stats = {}
    for energy_type in energy_data:
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
    
    write_outlier_stats(outliers_stats, output_dir, individual_id)

def process_constructs(ref_dir, alt_dir, snp_positions, energy_data, snp_counter, individual_id):
    """
    Обрабатывает файлы и конструкты, собирая данные об энергиях и SNP.
    """
    total_constructs = 0
    snp_constructs = 0
    error_constructs = 0
    for alt_file in os.listdir(alt_dir):
        if not alt_file.endswith("EF.csv"):
            continue
        ref_file = alt_file.replace(f"SEQ-g38_Mt-Short_Test-test_individual_{individual_id}", "SEQ-g38_Mt-Short_Test")
        ref_path = os.path.join(ref_dir, ref_file)
        alt_path = os.path.join(alt_dir, alt_file)
        if not os.path.exists(ref_path):
            logger.warning(f"Референсный файл не найден: {ref_path}")
            continue
        try:
            ref_df = pd.read_csv(ref_path)
            alt_df = pd.read_csv(alt_path)
            if ref_df.empty or alt_df.empty:
                logger.warning(f"Один из файлов пуст: {alt_file}")
                continue
            for idx, (_, alt_row) in enumerate(alt_df.iterrows()):
                total_constructs += 1
                construct_id = alt_row['ConstructID']
                try:
                    snps_in_construct = get_snps_in_construct(construct_id, snp_positions)
                    if snps_in_construct:
                        snp_constructs += 1
                        selected_snp = min(snps_in_construct)
                        snp_counter[selected_snp] += 1
                    else:
                        selected_snp = None
                except Exception as e:
                    logger.error(f"Ошибка обработки SNP для конструкта {construct_id}: {e}")
                    error_constructs += 1
                    selected_snp = None
                if idx < len(ref_df):
                    ref_row = ref_df.iloc[idx]
                    for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
                        if energy_type in alt_row and energy_type in ref_row:
                            alt_val = alt_row[energy_type]
                            ref_val = ref_row[energy_type]
                            if not np.isnan(alt_val) and not np.isnan(ref_val):
                                energy_data[energy_type]['ref'].append(ref_val)
                                energy_data[energy_type]['alt'].append(alt_val)
                                energy_data[energy_type]['snp_value'].append(selected_snp)
                else:
                    logger.warning(f"Нет соответствующей строки в референсном файле для конструкта {construct_id}")
        except Exception as e:
            logger.error(f"Ошибка при обработке файла {alt_file}: {e}")
    return total_constructs, snp_constructs, error_constructs

def log_statistics(total_constructs, snp_constructs, error_constructs, snp_counter):
    """
    Логирует статистику обработки конструктов.
    """
    logger.info(f"Всего обработано конструктов: {total_constructs}")
    logger.info(f"Конструктов с SNP: {snp_constructs}")
    logger.info(f"Конструктов без SNP: {total_constructs - snp_constructs}")
    logger.info(f"Конструктов с ошибками обработки: {error_constructs}")
    sorted_snps = sorted(snp_counter.items(), key=lambda x: x[1], reverse=True)[:10]
    logger.info("Топ-10 самых частых SNP:")
    for snp, count in sorted_snps:
        logger.info(f"  SNP {snp}: {count} конструктов")

def write_outlier_stats(outliers_stats, output_dir, individual_id):
    """
    Записывает статистику выбросов в файл.
    """
    stats_path = Path(output_dir) / f"test_individual_{individual_id}_outliers_statistics.txt"
    with open(stats_path, 'w', encoding='utf-8') as f:
        f.write("Статистика выбросов по типам энергии:\n")
        f.write("=" * 50 + "\n")
        for energy_type, stats in outliers_stats.items():
            f.write(f"{energy_type}:\n")
            f.write(f"  Всего точек: {stats['total_points']}\n")
            f.write(f"  Средняя разница (ref - alt): {stats['mean_diff']:.4f}\n")
            f.write(f"  Стандартное отклонение: {stats['std_diff']:.4f}\n")
            f.write(f"  Верхние выбросы (> +2std): {stats['upper_outliers']} ({stats['upper_outliers']/stats['total_points']*100:.2f}%)\n")
            f.write(f"  Нижние выбросы (< -2std): {stats['lower_outliers']} ({stats['lower_outliers']/stats['total_points']*100:.2f}%)\n")
            f.write("\n")
    logger.info(f"Статистика по выбросам сохранена: {stats_path}")

def main():
    # Базовые пути
    base_dir = "D:/pythonProject/MitoFragility/MitoFragilityScore/Energies"
    output_base_dir = "D:/pythonProject/MitoFragility/DataPreparing/plots/output"
    snp_base_dir = "D:/pythonProject/MitoFragility/MitoFragilityScore/Sequences/Relative"
    
    # Фиксированная референсная директория
    ref_dir = os.path.join(base_dir, "SEQ-g38_Mt-Short_Test")
    
    # Создаем базовую выходную директорию
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Находим все директории тестов
    individual_dirs = find_individual_dirs(base_dir)
    
    if not individual_dirs:
        logger.warning("Не найдено ни одной директории теста!")
        return
    
    logger.info(f"Найдено директорий теста: {len(individual_dirs)}")
    
    # Обрабатываем каждую директорию теста
    for alt_dir, individual_id in individual_dirs:
        # Формируем путь к файлу SNP
        snp_file_path = os.path.join(snp_base_dir, f"test_individual_{individual_id}.csv")
        
        # Проверяем существование необходимых путей
        if not os.path.exists(alt_dir):
            logger.warning(f"Директория теста не существует: {alt_dir}")
            continue
            
        if not os.path.exists(snp_file_path):
            logger.warning(f"Файл SNP не найден: {snp_file_path}")
            continue
            
        try:
            # Обрабатываем текущий тест
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