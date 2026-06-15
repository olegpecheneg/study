import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Optional
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

# Импорты модулей DataPreparing
try:
    from .data_loader import extract_energy_data_from_master_json
    from .stats_utils import (
        calculate_full_distribution_stats,
        fit_gpd_to_tail,
        fit_gpd_to_left_tail,
        calculate_survival_and_cdf,
        calculate_wilcoxon_test,
    )
    from .plots import create_interactive_analysis_plot, create_statistics_table
    from .html_generator import generate_html_report, save_json_results
    from .utils import setup_logger, get_mutation_id
    from .html_builder import generate_plotly_html_from_master_json
except ImportError:
    from data_loader import extract_energy_data_from_master_json
    from stats_utils import (
        calculate_full_distribution_stats,
        fit_gpd_to_tail,
        fit_gpd_to_left_tail,
        calculate_survival_and_cdf,
        calculate_wilcoxon_test,
    )
    from plots import create_interactive_analysis_plot, create_statistics_table
    from html_generator import generate_html_report, save_json_results
    from utils import setup_logger, get_mutation_id
    from html_builder import generate_plotly_html_from_master_json


def extract_energy_data_from_prepared_json(prepared_data: Dict) -> Dict:
    """Извлекает энергии из подготовленных данных (альтернатива для подготовленных JSON)."""
    if not isinstance(prepared_data, dict) or 'energy_data' not in prepared_data:
        return {}
    result = {}
    for energy_type, energies in prepared_data['energy_data'].items():
        result[energy_type] = {
            'ref': np.array(energies['ref']),
            'alt': np.array(energies['alt']),
            'diff': np.array(energies['diff'])
        }
    return result


def analyze_mutation_energy_distribution(master_json_path: str, output_dir: str, is_prepared_data: bool = False):
    """Анализ распределения энергий для мутации."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(log_file=output_path / "mutation_analysis.log")
    mutation_id = get_mutation_id(master_json_path)

    logger.info(f"Начало анализа для мутации: {mutation_id}")
    logger.info(f"JSON файл: {master_json_path}")
    logger.info(f"Выходная директория: {output_dir}")
    logger.info(f"Подготовленные данные: {is_prepared_data}")

    # Извлекаем данные
    if is_prepared_data:
        try:
            with open(master_json_path, 'r', encoding='utf-8') as f:
                prepared_data = json.load(f)
            energy_data = extract_energy_data_from_prepared_json(prepared_data)
        except Exception as e:
            logger.error(f"Ошибка загрузки подготовленных данных: {e}")
            return {}
    else:
        energy_data = extract_energy_data_from_master_json(master_json_path)

    if not energy_data:
        logger.error("Не удалось извлечь данные")
        return {}

    all_results = {}

    for energy_type, data_dict in energy_data.items():
        logger.info(f"\n{'='*60}")
        logger.info(f"Анализ {energy_type}:")
        logger.info(f"{'='*60}")

        energy_dir = output_path / mutation_id / energy_type
        energy_dir.mkdir(parents=True, exist_ok=True)

        # Wilcoxon test один раз на тип энергии
        wilcoxon_results = None
        if 'ref' in data_dict and 'alt' in data_dict:
            ref_energies = data_dict['ref']
            alt_energies = data_dict['alt']
            if len(ref_energies) == len(alt_energies):
                logger.info(f"  Расчет Wilcoxon signed-rank test (ref vs alt)...")
                wilcoxon_results = calculate_wilcoxon_test(ref_energies, alt_energies)
                if 'error' not in wilcoxon_results:
                    p_val = wilcoxon_results['wilcoxon_p_value']
                    mean_diff = wilcoxon_results['mean_difference']
                    direction = "ref > alt" if mean_diff > 0 else "ref < alt"
                    significance = "ЗНАЧИМО" if p_val < 0.05 else "не значимо"
                    logger.info(f"    p = {p_val:.4f} ({significance})")
                    logger.info(f"    Средняя разность = {mean_diff:.3f} ({direction})")
                else:
                    logger.warning(f"    Ошибка Wilcoxon: {wilcoxon_results['error']}")

        results_for_type = {}

        for data_type, energies in data_dict.items():
            if len(energies) < 10:
                logger.warning(f"  {data_type}: только {len(energies)} точек")
                continue

            logger.info(f"\n  {data_type}: {len(energies)} точек")

            distribution_stats = calculate_full_distribution_stats(energies)
            gpd_result_right = fit_gpd_to_tail(energies)
            gpd_result_left = fit_gpd_to_left_tail(energies)
            survival_data = calculate_survival_and_cdf(energies)

            results_for_type[data_type] = {
                'distribution_stats': distribution_stats,
                'gpd_result_right': gpd_result_right,
                'gpd_result_left': gpd_result_left,
                'survival_data': survival_data,
                'wilcoxon_results': wilcoxon_results
            }

            plot_fig = create_interactive_analysis_plot(
                energy_type, data_type, energies, mutation_id,
                wilcoxon_results=wilcoxon_results
            )
            table_fig = create_statistics_table(
                energy_type, data_type, distribution_stats,
                gpd_result_right, gpd_result_left, survival_data,
                wilcoxon_results
            )

            html_path = generate_html_report(
                mutation_id, energy_type, data_type,
                plot_fig, table_fig, energies,
                distribution_stats, gpd_result_right, gpd_result_left, survival_data,
                wilcoxon_results, energy_dir
            )
            logger.info(f"    Сохранен HTML: {html_path}")

            json_path = save_json_results(
                mutation_id, energy_type, data_type,
                distribution_stats, gpd_result_right, gpd_result_left, survival_data,
                wilcoxon_results, energy_dir
            )
            logger.info(f"    Сохранен JSON: {json_path}")

        all_results[energy_type] = results_for_type

    logger.info(f"\n{'='*60}")
    logger.info(f"Анализ завершен для мутации {mutation_id}")
    logger.info(f"Результаты сохранены в: {output_dir}")
    logger.info(f"{'='*60}")
    return all_results


def prepare_comparison_data(ref_json_path: str, alt_json_path: str, logger) -> Dict:
    """Подготавливает данные для сравнения двух JSON файлов."""
    try:
        with open(ref_json_path, 'r', encoding='utf-8') as f:
            ref_data = json.load(f)
        with open(alt_json_path, 'r', encoding='utf-8') as f:
            alt_data = json.load(f)

        ref_id = Path(ref_json_path).stem.split('_')[0]
        alt_id = Path(alt_json_path).stem.split('_')[0]
        combined_id = f"{ref_id}-{alt_id}"

        logger.info(f"Создание сравнения {combined_id}")
        logger.info(f"Из файла {ref_id} всегда берем alt_energy как 'ref'")
        logger.info(f"Из файла {alt_id} всегда берем alt_energy как 'alt'")

        ref_by_type = {}
        alt_by_type = {}

        if isinstance(ref_data, list):
            for point in ref_data:
                if isinstance(point, dict):
                    energy_type = point.get('energy_type', 'Energy')
                    if energy_type not in ref_by_type:
                        ref_by_type[energy_type] = []
                    try:
                        alt_value = float(point.get('alt_energy', 0))
                        ref_by_type[energy_type].append(alt_value)
                    except (ValueError, TypeError):
                        continue
        else:
            logger.error(f"Неизвестный формат ref данных: {type(ref_data)}")
            return None

        if isinstance(alt_data, list):
            for point in alt_data:
                if isinstance(point, dict):
                    energy_type = point.get('energy_type', 'Energy')
                    if energy_type not in alt_by_type:
                        alt_by_type[energy_type] = []
                    try:
                        alt_value = float(point.get('alt_energy', 0))
                        alt_by_type[energy_type].append(alt_value)
                    except (ValueError, TypeError):
                        continue
        else:
            logger.error(f"Неизвестный формат alt данных: {type(alt_data)}")
            return None

        common_types = set(ref_by_type.keys()) & set(alt_by_type.keys())
        if not common_types:
            logger.warning(f"Нет общих типов энергии между файлами")
            return None

        logger.info(f"Общие типы энергии: {common_types}")

        new_energy_data = {}
        for energy_type in common_types:
            ref_values = ref_by_type[energy_type]
            alt_values = alt_by_type[energy_type]
            min_len = min(len(ref_values), len(alt_values))
            if min_len == 0:
                continue
            ref_values = ref_values[:min_len]
            alt_values = alt_values[:min_len]
            diff_values = [float(r) - float(a) for r, a in zip(ref_values, alt_values)]
            new_energy_data[energy_type] = {
                'ref': [float(v) for v in ref_values],
                'alt': [float(v) for v in alt_values],
                'diff': diff_values
            }
            mean_ref = np.mean(ref_values)
            mean_alt = np.mean(alt_values)
            mean_diff = np.mean(diff_values)
            logger.info(f"Тип {energy_type}: {min_len} пар")
            logger.info(f"  Среднее ref: {mean_ref:.3f}, alt: {mean_alt:.3f}, diff: {mean_diff:.3f}")

        return {
            'mutation_id': combined_id,
            'energy_data': new_energy_data
        }
    except Exception as e:
        logger.error(f"Ошибка подготовки данных для сравнения: {e}", exc_info=True)
        return None


def process_all_json_files(data_dir: Path, ref_data_dir: Optional[Path], output_base: Path, logger: logging.Logger, recursive: bool = True):
    """
    Полный пайплайн обработки всех JSON файлов.
    Параметры:
        data_dir: папка с файлами *_master_data.json (может содержать подпапки)
        ref_data_dir: папка с такими же файлами для сравнения (может быть None)
        output_base: базовая папка для сохранения результатов анализа
        logger: объект логгера
        recursive: искать ли рекурсивно в подпапках
    """
    logger.info("=" * 80)
    logger.info("ЗАПУСК ПОЛНОГО ПАЙПЛАЙНА")
    logger.info("=" * 80)

    # Часть 1: Обработка всех JSON из data_dir (рекурсивно)
    logger.info("\n" + "=" * 80)
    logger.info("ЧАСТЬ 1: Анализ всех JSON файлов из data/ (рекурсивно)")
    logger.info("=" * 80)

    if recursive:
        data_json_files = list(data_dir.rglob("*_master_data.json"))
    else:
        data_json_files = list(data_dir.glob("*_master_data.json"))

    logger.info(f"Найдено {len(data_json_files)} JSON файлов")

    for json_file in data_json_files:
        try:
            # Определяем относительный путь от data_dir
            rel_path = json_file.relative_to(data_dir).parent
            # ID мутации – имя файла без _master_data
            file_stem = json_file.stem
            if file_stem.endswith("_master_data"):
                mutation_id = file_stem[:-12]  # удаляем "_master_data"
            else:
                mutation_id = file_stem
            # Выходная папка: output_base/analysis_results/rel_path/mutation_id
            analysis_base = output_base / "analysis_results"
            output_dir = analysis_base / rel_path / mutation_id
            logger.info(f"\nОбработка файла: {json_file}")
            logger.info(f"ID мутации: {mutation_id}")
            logger.info(f"Выходная директория: {output_dir}")

            # Стандартный анализ
            results = analyze_mutation_energy_distribution(
                master_json_path=str(json_file),
                output_dir=str(output_dir),
                is_prepared_data=False
            )

            # Scatter plot (интерактивные графики)
            logger.info("Создание scatter plot...")
            scat_dir = output_dir / "scatter_plots"
            scat_dir.mkdir(exist_ok=True)
            generated_files = generate_plotly_html_from_master_json(
                master_json_path=str(json_file),
                output_dir=str(scat_dir)
            )
            logger.info(f"Создано {len(generated_files)} scatter plots")

        except Exception as e:
            logger.error(f"Ошибка при обработке {json_file}: {e}", exc_info=True)

    # Часть 2: Сравнение (только если задан ref_data_dir и он существует)
    if ref_data_dir and ref_data_dir.exists():
        logger.info("\n" + "=" * 80)
        logger.info("ЧАСТЬ 2: Сравнение JSON из ref_data/ и data/")
        logger.info("=" * 80)

        if recursive:
            ref_json_files = list(ref_data_dir.rglob("*_master_data.json"))
        else:
            ref_json_files = list(ref_data_dir.glob("*_master_data.json"))

        logger.info(f"Найдено {len(ref_json_files)} JSON файлов в ref_data/")

        for ref_json in ref_json_files:
            try:
                # Находим соответствующий файл в data_dir относительно структуры
                rel_path = ref_json.relative_to(ref_data_dir)
                data_json = data_dir / rel_path
                if not data_json.exists():
                    # Пробуем найти по имени файла в любом месте (не строго по структуре)
                    candidates = list(data_dir.rglob(ref_json.name))
                    if candidates:
                        data_json = candidates[0]
                    else:
                        logger.warning(f"Не найден соответствующий файл в data/: {ref_json.name}")
                        continue

                file_id = ref_json.stem.replace("_master_data", "")
                logger.info(f"\nСравнение файлов:")
                logger.info(f"  REF (из ref_data/): {ref_json}")
                logger.info(f"  ALT (из data/): {data_json}")
                logger.info(f"  ID: {file_id}")

                prepared_data = prepare_comparison_data(
                    ref_json_path=str(ref_json),
                    alt_json_path=str(data_json),
                    logger=logger
                )
                if prepared_data is None:
                    logger.error(f"Не удалось подготовить данные для сравнения {file_id}")
                    continue

                combined_id = f"{file_id}-{file_id}"
                comparison_base = output_base / "analysis_results_comparison"
                output_dir = comparison_base / rel_path.parent / combined_id
                output_dir.mkdir(parents=True, exist_ok=True)

                temp_json = output_dir / f"{combined_id}_prepared.json"
                with open(temp_json, 'w', encoding='utf-8') as f:
                    json.dump(prepared_data, f, indent=2, ensure_ascii=False)

                logger.info(f"Подготовленные данные сохранены: {temp_json}")

                results = analyze_mutation_energy_distribution(
                    master_json_path=str(temp_json),
                    output_dir=str(output_dir),
                    is_prepared_data=True
                )

                # Scatter plots для сравнения
                scat_ref_dir = output_dir / "scatter_ref"
                scat_ref_dir.mkdir(exist_ok=True)
                try:
                    generated_ref = generate_plotly_html_from_master_json(
                        master_json_path=str(ref_json),
                        output_dir=str(scat_ref_dir)
                    )
                    logger.info(f"Создано {len(generated_ref)} scatter plots для ref")
                except Exception as e:
                    logger.error(f"Ошибка создания scatter plot для ref: {e}")

                scat_alt_dir = output_dir / "scatter_alt"
                scat_alt_dir.mkdir(exist_ok=True)
                try:
                    generated_alt = generate_plotly_html_from_master_json(
                        master_json_path=str(data_json),
                        output_dir=str(scat_alt_dir)
                    )
                    logger.info(f"Создано {len(generated_alt)} scatter plots для alt")
                except Exception as e:
                    logger.error(f"Ошибка создания scatter plot для alt: {e}")

            except Exception as e:
                logger.error(f"Ошибка при сравнении {ref_json}: {e}", exc_info=True)
    else:
        logger.info("ref_data_dir не задан или не существует, сравнение пропущено.")

    logger.info("\n" + "=" * 80)
    logger.info("ВЕСЬ ПАЙПЛАЙН ЗАВЕРШЕН")
    logger.info("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Запуск полного пайплайна: рекурсивный анализ всех master_data.json в подпапках."
    )
    parser.add_argument(
        "data_dir",
        help="Путь к папке, содержащей JSON-файлы (*_master_data.json) в корне или подпапках."
    )
    parser.add_argument(
        "-r", "--ref_data_dir",
        default=None,
        help="(Опционально) Папка с ref JSON-файлами для сравнения. Если не указано, сравнение не выполняется."
    )
    parser.add_argument(
        "-o", "--output_dir",
        default=None,
        help="(Опционально) Базовая папка для сохранения результатов. По умолчанию внутри data_dir/analysis_results"
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Отключить рекурсивный поиск (искать только в корне data_dir)"
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        print(f"ОШИБКА: Папка с данными не найдена: {data_dir}")
        sys.exit(1)

    # Для логов используем либо указанную output_dir, либо data_dir
    log_dir = Path(args.output_dir) if args.output_dir else data_dir
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(log_file=log_dir / "full_pipeline.log")

    logger.info("=" * 80)
    logger.info("ЗАПУСК ПОЛНОГО ПАЙПЛАЙНА АНАЛИЗА (РЕКУРСИВНЫЙ ПОИСК)")
    logger.info("=" * 80)
    logger.info(f"Директория с данными: {data_dir}")
    if args.ref_data_dir:
        ref_data_dir = Path(args.ref_data_dir)
        logger.info(f"Директория ref_data: {ref_data_dir}")
    else:
        ref_data_dir = None
        logger.info("Сравнение с ref_data отключено (не указан -r)")

    recursive = not args.no_recursive
    logger.info(f"Рекурсивный поиск: {'ВКЛЮЧЕН' if recursive else 'ВЫКЛЮЧЕН'}")

    output_base_dir = Path(args.output_dir) if args.output_dir else data_dir
    output_base_dir.mkdir(parents=True, exist_ok=True)

    # Проверяем наличие хотя бы одного файла (рекурсивно)
    if recursive:
        found = list(data_dir.rglob("*_master_data.json"))
    else:
        found = list(data_dir.glob("*_master_data.json"))
    if not found:
        logger.warning(f"Не найдено JSON-файлов вида *_master_data.json в {data_dir}")
        print(f"ВНИМАНИЕ: Не найдено JSON-файлов в {data_dir}")
        return

    process_all_json_files(data_dir, ref_data_dir, output_base_dir, logger, recursive=recursive)

    print("\n" + "=" * 80)
    print("ГОТОВО!")
    print(f"Результаты сохранены в {output_base_dir / 'analysis_results'}")
    if ref_data_dir:
        print(f"Сравнение сохранено в {output_base_dir / 'analysis_results_comparison'}")
    print(f"Лог работы: {log_dir / 'full_pipeline.log'}")
    print("=" * 80)


if __name__ == "__main__":
    main()