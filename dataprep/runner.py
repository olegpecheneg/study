import os
from typing import List, Tuple, Dict, Any, Optional
from pathlib import Path
from logging import getLogger
import pandas as pd

from .config import resolve_project_paths
from .snputils import load_snp_data, get_snps_in_construct
from .parsers import parse_construct_id
from .colors import generate_distinct_colors
from .plot_export import plot_energy_comparison
from .io import save_master_json, write_outlier_stats

logger = getLogger(__name__)


def find_points_by_sequence(search_file: str, search_sequence: str, case_sensitive: bool = False) -> List[Dict]:
    if not case_sensitive:
        search_sequence = search_sequence.upper()
    try:
        df = pd.read_csv(search_file)
        results = []
        for _, row in df.iterrows():
            for window_num in range(1, 7):
                window_seq = str(row.get(f'window{window_num}', ''))
                if not case_sensitive:
                    window_seq = window_seq.upper()
                if search_sequence in window_seq:
                    results.append({
                        'point_index': row.get('point_index'),
                        'construct_id': row.get('construct_id'),
                        'ref_energy': row.get('ref_energy'),
                        'alt_energy': row.get('alt_energy'),
                        'snp': row.get('snp'),
                        'window_found': f'window{window_num}',
                        'sequence_found': window_seq
                    })
                    break
        return results
    except Exception as e:
        logger.error(f"Ошибка поиска в файле {search_file}: {e}")
        return []


def find_individual_dirs(base_dir: str) -> List[Tuple[str, str]]:
    individual_dirs: List[Tuple[str, str]] = []
    for d in os.listdir(base_dir):
        if d.startswith("SEQ-g38_Mt-Short_Test-test_individual_") and os.path.isdir(os.path.join(base_dir, d)):
            parts = d.split('_')
            if parts:
                individual_id = parts[-1]
                individual_dirs.append((os.path.join(base_dir, d), individual_id))
    return individual_dirs


def process_individual(ref_dir: str, alt_dir: str, snp_file_path: str, output_base_dir: str, individual_id: str) -> None:
    logger.info(f"Обработка индивидуума: {individual_id}")
    snp_positions = load_snp_data(snp_file_path)
    output_dir = os.path.join(output_base_dir, f"test_individual_{individual_id}")
    os.makedirs(output_dir, exist_ok=True)

    project_paths = resolve_project_paths()
    constructs_base_dir = project_paths.get('CONSTRUCTS_DIR')
    individual_constructs_dir = os.path.join(constructs_base_dir, f"SEQ-g38_Mt-Short_Test-test_individual_{individual_id}")

    energy_data = {
        'EnergyLeft': {'ref': [], 'alt': [], 'snp_value': [], 'construct_details': []},
        'EnergyRight': {'ref': [], 'alt': [], 'snp_value': [], 'construct_details': []},
        'Energy': {'ref': [], 'alt': [], 'snp_value': [], 'construct_details': []}
    }

    snp_counter = {}
    construct_categories = {'with_snp': 0, 'without_snp': 0, 'error': 0}

    all_snps = sorted(snp_positions)
    snp_colors = {snp: generate_distinct_colors(len(all_snps))[i] for i, snp in enumerate(all_snps)} if all_snps else {}

    # Загрузка CF файлов
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
                        continue
                    lower_map = {col.lower().strip(): col for col in cf_df.columns}
                    construct_id_col = None
                    for variant in desired_fields['construct_id']:
                        if variant in lower_map:
                            construct_id_col = lower_map[variant]
                            break
                    if construct_id_col is None:
                        logger.warning(f"CF-файл {cf_file} пропущен: не найден столбец ConstructID")
                        continue
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
                except Exception as e:
                    logger.warning(f"Ошибка загрузки конструкт файла {cf_file}: {e}")
    else:
        logger.warning(f"Директория конструктов не найдена: {individual_constructs_dir}")

    logger.info(f"Загружено {cf_files_processed} CF файлов, всего конструктов в памяти: {len(construct_data_map)}")

    # Обработка EF файлов
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
            required_columns = ['ConstructID', 'EnergyLeft', 'EnergyRight', 'Energy']
            missing_in_ref = [col for col in required_columns if col not in ref_df.columns]
            missing_in_alt = [col for col in required_columns if col not in alt_df.columns]
            if missing_in_ref or missing_in_alt:
                logger.warning(f"Отсутствуют нужные столбцы в {file}: {missing_in_ref + missing_in_alt}")
                continue
            for i in range(min(len(ref_df), len(alt_df))):
                construct_id = alt_df.iloc[i]['ConstructID']
                construct_details_from_cf = construct_data_map.get(construct_id, {})
                if construct_details_from_cf:
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

                arm_size, center, arm3_start, arm4_start = parse_construct_id(construct_id)
                if arm_size is None:
                    construct_categories['error'] += 1
                    selected_snp = None
                else:
                    snps = get_snps_in_construct(snp_positions, arm_size, center, arm3_start, arm4_start)
                    if snps:
                        construct_categories['with_snp'] += 1
                        selected_snp = min(snps)
                        snp_counter[selected_snp] = snp_counter.get(selected_snp, 0) + 1
                    else:
                        construct_categories['without_snp'] += 1
                        selected_snp = None

                for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
                    ref_energy = ref_df.iloc[i][energy_type]
                    alt_energy = alt_df.iloc[i][energy_type]
                    try:
                        if not pd.isna(ref_energy) and not pd.isna(alt_energy):
                            energy_data[energy_type]['ref'].append(float(ref_energy))
                            energy_data[energy_type]['alt'].append(float(alt_energy))
                            energy_data[energy_type]['snp_value'].append(selected_snp)
                            energy_data[energy_type]['construct_details'].append(construct_details)
                    except Exception:
                        pass

        except Exception as e:
            logger.error(f"Ошибка обработки файла {file}: {e}")
            construct_categories['error'] += 1

    total_constructs = sum(construct_categories.values())
    logger.info(f"Категории конструктов для {individual_id}: всего {total_constructs}, с SNP {construct_categories['with_snp']}, без SNP {construct_categories['without_snp']}, с ошибкой {construct_categories['error']}")

    # Генерация графиков
    outliers_stats = {}
    master_records = []
    for energy_type in ['EnergyLeft', 'EnergyRight', 'Energy']:
        ref_vals = energy_data[energy_type]['ref']
        alt_vals = energy_data[energy_type]['alt']
        snp_vals = energy_data[energy_type]['snp_value']
        construct_details = energy_data[energy_type]['construct_details']
        if ref_vals and alt_vals:
            stats = plot_energy_comparison(ref_vals, alt_vals, snp_vals, snp_colors, energy_type, output_dir, individual_id, construct_details)
            outliers_stats[energy_type] = stats
            if stats and isinstance(stats, dict) and stats.get('point_records'):
                master_records.extend(stats.get('point_records'))
        else:
            logger.warning(f"Нет данных для {energy_type}")

    try:
        master_json_path = save_master_json(output_dir, individual_id, master_records)
    except Exception as e:
        logger.warning(f"Не удалось сохранить master JSON для {individual_id}: {e}")

    try:
        write_outlier_stats(outliers_stats, output_dir, individual_id, job_prefix='')
    except Exception as e:
        logger.warning(f"Не удалось записать статистику выбросов: {e}")


def search_sequences_wrapper():
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
    if search_mode:
        search_sequences_wrapper()
        return

    paths = resolve_project_paths()
    if base_dir is None:
        base_dir = paths.get('ENERGIES_DIR')
    if output_base_dir is None:
        output_base_dir = paths.get('OUTPUT_DIR')
    if snp_base_dir is None:
        snp_base_dir = paths.get('SEQUENCES_RELATIVE_DIR')

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
            process_individual(ref_dir, alt_dir, snp_file_path, output_base_dir, individual_id)
        except Exception as e:
            logger.error(f"Ошибка при обработке теста {individual_id}: {str(e)}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "search":
        search_sequences_wrapper()
    else:
        main()
