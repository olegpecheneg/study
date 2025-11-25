"""Вывод и сохранение результатов обработки mt_DNA."""

from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
from logging import getLogger
import json
import pandas as pd

logger = getLogger(__name__)


def save_individual_results(
    output_dir: Path,
    individual_id: str,
    mutations: List[Dict[str, Any]],
    energies: Dict[str, Any],
    sequences: Dict[str, str]
) -> bool:
    """Сохранить результаты обработки индивидуума.

    Args:
        output_dir: Директория для вывода.
        individual_id: ID индивидуума.
        mutations: Список мутаций с их свойствами.
        energies: Словарь с энергиями.
        sequences: Словарь с последовательностями (ref, mutated, etc).

    Returns:
        True если успешно, False иначе.
    """
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Сохраняем мутации как JSON
        mutations_file = output_dir / f"{individual_id}_mutations.json"
        with open(mutations_file, 'w', encoding='utf-8') as f:
            json.dump(mutations, f, indent=2, ensure_ascii=False)
        logger.info(f"Мутации сохранены в {mutations_file}")
        
        # Сохраняем энергии как JSON
        energies_file = output_dir / f"{individual_id}_energies.json"
        with open(energies_file, 'w', encoding='utf-8') as f:
            json.dump(energies, f, indent=2, ensure_ascii=False)
        logger.info(f"Энергии сохранены в {energies_file}")
        
        # Сохраняем последовательности в FASTA
        sequences_file = output_dir / f"{individual_id}_sequences.fasta"
        with open(sequences_file, 'w', encoding='utf-8') as f:
            for seq_name, seq_data in sequences.items():
                f.write(f">{individual_id}_{seq_name}\n")
                # Записываем с переносом на каждые 60 символов
                for i in range(0, len(seq_data), 60):
                    f.write(f"{seq_data[i:i+60]}\n")
        logger.info(f"Последовательности сохранены в {sequences_file}")
        
        return True
    except Exception as e:
        logger.error(f"Ошибка сохранения результатов для {individual_id}: {e}")
        return False


def save_mutations_to_csv(
    output_file: Path,
    mutations: List[Dict[str, Any]]
) -> bool:
    """Сохранить мутации в CSV.

    Args:
        output_file: Путь к выходному файлу.
        mutations: Список мутаций.

    Returns:
        True если успешно, False иначе.
    """
    try:
        df = pd.DataFrame(mutations)
        df.to_csv(output_file, index=False, encoding='utf-8')
        logger.info(f"Мутации сохранены в CSV: {output_file}")
        return True
    except Exception as e:
        logger.error(f"Ошибка сохранения мутаций в CSV: {e}")
        return False


def save_energy_comparison(
    output_file: Path,
    comparison_data: Dict[str, Any]
) -> bool:
    """Сохранить сравнение энергий в JSON.

    Args:
        output_file: Путь к выходному файлу.
        comparison_data: Словарь со сравнением энергий.

    Returns:
        True если успешно, False иначе.
    """
    try:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(comparison_data, f, indent=2, ensure_ascii=False)
        logger.info(f"Сравнение энергий сохранено: {output_file}")
        return True
    except Exception as e:
        logger.error(f"Ошибка сохранения сравнения энергий: {e}")
        return False


def create_summary_report(
    output_file: Path,
    summary_data: Dict[str, Any]
) -> bool:
    """Создать итоговый отчёт в JSON.

    Args:
        output_file: Путь к выходному файлу.
        summary_data: Данные для итогового отчёта.

    Returns:
        True если успешно, False иначе.
    """
    try:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Добавляем метаинформацию
        report = {
            'metadata': {
                'version': '1.0',
                'description': 'mt_DNA analysis summary report'
            },
            'summary': summary_data
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        logger.info(f"Итоговый отчёт сохранён: {output_file}")
        return True
    except Exception as e:
        logger.error(f"Ошибка создания итогового отчёта: {e}")
        return False


def merge_results_to_dataframe(
    results_list: List[Dict[str, Any]]
) -> Optional[pd.DataFrame]:
    """Объединить результаты нескольких индивидуумов в DataFrame.

    Args:
        results_list: Список словарей с результатами.

    Returns:
        DataFrame или None в случае ошибки.
    """
    try:
        df = pd.DataFrame(results_list)
        logger.info(f"Объединено {len(df)} результатов в DataFrame")
        return df
    except Exception as e:
        logger.error(f"Ошибка объединения результатов: {e}")
        return None


def append_to_master_result_file(
    master_file: Path,
    individual_results: Dict[str, Any]
) -> bool:
    """Добавить результаты индивидуума к мастер-файлу.

    Args:
        master_file: Путь к мастер-файлу.
        individual_results: Результаты индивидуума.

    Returns:
        True если успешно, False иначе.
    """
    try:
        master_file = Path(master_file)
        
        # Загружаем существующие данные или создаём новый файл
        if master_file.exists():
            with open(master_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                if not isinstance(data, list):
                    data = [data]
        else:
            data = []
        
        data.append(individual_results)
        
        # Сохраняем обновленные данные
        master_file.parent.mkdir(parents=True, exist_ok=True)
        with open(master_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        logger.info(f"Результаты добавлены в мастер-файл: {master_file}")
        return True
    except Exception as e:
        logger.error(f"Ошибка добавления к мастер-файлу: {e}")
        return False
