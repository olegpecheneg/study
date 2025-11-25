from pathlib import Path
from typing import List, Dict, Any
from logging import getLogger
import json
import pandas as pd

logger = getLogger(__name__)


def save_detailed_data(output_dir: str, energy_type: str, individual_id: str, construct_details: List[Dict[str, Any]], ref_data: List[float], alt_data: List[float], snp_values: List[Any], upper_outliers=None, lower_outliers=None):
    detailed_data_path = Path(output_dir) / f"{energy_type}_detailed_data.txt"
    try:
        with open(detailed_data_path, 'w', encoding='utf-8') as f:
            f.write(f"Детальные данные для {energy_type} - Индивидуум {individual_id}\n")
            f.write("=" * 80 + "\n\n")
            min_len = min(len(ref_data), len(alt_data), len(snp_values), len(construct_details))
            for i in range(min_len):
                details = construct_details[i]
                f.write(f"Точка {i}:\n")
                f.write(f"  Координаты: Ref={ref_data[i]:.2f}, Alt={alt_data[i]:.2f}\n")
                f.write(f"  ID конструкта: {details['construct_id']}\n")
                f.write(f"  SNP: {snp_values[i] if snp_values[i] else 'Нет'}\n")
                f.write("-" * 80 + "\n")
        logger.info(f"Детальные данные сохранены: {detailed_data_path}")
        return str(detailed_data_path)
    except Exception as e:
        logger.error(f"Ошибка при записи детальных данных {detailed_data_path}: {e}")
        raise


def save_search_csv(output_dir: str, energy_type: str, search_data: List[Dict[str, Any]]):
    if not search_data:
        logger.info("Данные для поиска отсутствуют — CSV не создан")
        return None
    search_data_path = Path(output_dir) / f"{energy_type}_sequence_search.csv"
    try:
        df = pd.DataFrame(search_data)
        df.to_csv(search_data_path, index=False, encoding='utf-8')
        logger.info(f"Данные для поиска сохранены: {search_data_path}")
        return str(search_data_path)
    except Exception as e:
        logger.error(f"Ошибка при записи search CSV {search_data_path}: {e}")
        raise


def write_outlier_stats(outliers_stats: Dict[str, Dict[str, Any]], output_dir: str, individual_id: str, job_prefix: str = "") -> str:
    stats_path = Path(output_dir) / f"{job_prefix}{individual_id}_outliers_statistics.txt"
    try:
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
        return str(stats_path)
    except Exception as e:
        logger.error(f"Ошибка при записи статистики выбросов {stats_path}: {e}")
        raise


def save_master_json(output_dir: str, individual_id: str, master_records: List[Dict[str, Any]], job_prefix: str = "") -> str:
    master_json_path = Path(output_dir) / f"{job_prefix}{individual_id}_master_data.json"
    try:
        with open(master_json_path, 'w', encoding='utf-8') as mf:
            json.dump(master_records, mf, ensure_ascii=False, indent=2)
        logger.info(f"Master JSON с детальными данными сохранён: {master_json_path}")
        return str(master_json_path)
    except Exception as e:
        logger.warning(f"Не удалось сохранить master JSON для {individual_id}: {e}")
        raise
