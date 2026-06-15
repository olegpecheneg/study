#!/usr/bin/env python3
"""
Главный пайплайн: анализ JSON, аннотация, агрегация, volcano/strip, кластеризация.
Подробное логирование каждого SNV.
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from collections import defaultdict, Counter

import pandas as pd

# ----------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("FullPipeline")

# ----------------------------------------------------------------------
def find_master_json_files(data_dir: Path, recursive=True):
    return list(data_dir.rglob("*_master_data.json")) if recursive else list(data_dir.glob("*_master_data.json"))

def analyze_json_files(json_files):
    """Анализирует каждый master_data.json и возвращает DataFrame со статистикой."""
    snv_stats = defaultdict(lambda: {'Energy': 0, 'EnergyLeft': 0, 'EnergyRight': 0, 'total_records': 0})
    for jf in json_files:
        try:
            with open(jf) as f:
                data = json.load(f)
        except Exception as e:
            logger.error(f"Не удалось прочитать {jf}: {e}")
            continue
        if not data:
            logger.warning(f"Пустой файл: {jf}")
            continue
        pos = data[0].get('snp') or data[0].get('position')
        if pos is None:
            logger.warning(f"Нет позиции в {jf}")
            continue
        pos = int(pos)
        snv_stats[pos]['total_records'] += len(data)
        for rec in data:
            etype = rec.get('energy_type')
            if etype in ('Energy', 'EnergyLeft', 'EnergyRight'):
                snv_stats[pos][etype] += 1
    # Преобразуем в DataFrame
    rows = []
    for pos, d in snv_stats.items():
        rows.append({
            'position': pos,
            'Energy': d['Energy'],
            'EnergyLeft': d['EnergyLeft'],
            'EnergyRight': d['EnergyRight'],
            'total_records': d['total_records']
        })
    df = pd.DataFrame(rows).sort_values('position')
    return df

def run_script(script_name, args_list, step_name):
    """Запускает внешний скрипт и логирует результат."""
    cmd = [sys.executable, str(script_name)] + args_list
    logger.info(f"Запуск {step_name}: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"{step_name} завершён успешно")
        logger.debug(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} упал с кодом {e.returncode}")
        logger.error(e.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="Полный пайплайн обработки mtDNA данных")
    parser.add_argument("data_dir", type=Path, help="Корневая директория с master_data.json")
    parser.add_argument("-o", "--output_dir", type=Path, default=None, help="Базовая папка для всех результатов")
    parser.add_argument("--skip-annotate", action="store_true", help="Пропустить аннотацию")
    parser.add_argument("--skip-aggregate", action="store_true", help="Пропустить агрегатные графики")
    parser.add_argument("--skip-volcano", action="store_true", help="Пропустить volcano/strip plots")
    parser.add_argument("--skip-cluster", action="store_true", help="Пропустить кластеризацию")
    parser.add_argument("--fastas-dir", type=Path, help="Папка с FASTA для аннотации")
    parser.add_argument("--ref-fasta", type=Path, help="Референсный FASTA chrM")
    parser.add_argument("--effect-threshold", type=float, default=0.0)
    parser.add_argument("--max-constructs", type=int, default=500)
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Директория не найдена: {data_dir}")
        sys.exit(1)

    out_base = (args.output_dir or data_dir.parent / "full_pipeline_output").resolve()
    out_base.mkdir(parents=True, exist_ok=True)

    # Шаг 0: Анализ исходных JSON
    logger.info("=== Шаг 0: Анализ master_data.json ===")
    json_files = find_master_json_files(data_dir)
    logger.info(f"Найдено {len(json_files)} файлов")
    if not json_files:
        logger.error("Нет файлов, выход")
        sys.exit(1)

    snv_df = analyze_json_files(json_files)
    snv_df.to_csv(out_base / "snv_counts_per_energy.csv", index=False)
    logger.info(f"Статистика по SNV сохранена: {out_base / 'snv_counts_per_energy.csv'}")
    logger.info("Сводка по типам энергии:")
    logger.info(f"  Energy: SNV с >0 записей = {len(snv_df[snv_df['Energy']>0])}, всего записей = {snv_df['Energy'].sum()}")
    logger.info(f"  EnergyLeft: SNV с >0 записей = {len(snv_df[snv_df['EnergyLeft']>0])}, всего записей = {snv_df['EnergyLeft'].sum()}")
    logger.info(f"  EnergyRight: SNV с >0 записей = {len(snv_df[snv_df['EnergyRight']>0])}, всего записей = {snv_df['EnergyRight'].sum()}")
    logger.info("Первые 10 SNV с количеством пар по типам:")
    logger.info("\n" + snv_df.head(10).to_string())

    # Шаг 1: Аннотация
    if not args.skip_annotate:
        logger.info("=== Шаг 1: Аннотация SNV ===")
        annotate_args = [str(data_dir)]
        if args.fastas_dir:
            annotate_args.extend(["--fastas-dir", str(args.fastas_dir)])
        if args.ref_fasta:
            annotate_args.extend(["--ref-fasta", str(args.ref_fasta)])
        annotate_args.append("--force")
        script_path = Path(__file__).parent / "annotate_master_json.py"
        if not script_path.exists():
            logger.error(f"Не найден {script_path}, пропускаем аннотацию")
        else:
            run_script(script_path, annotate_args, "Аннотация")
    else:
        logger.info("Аннотация пропущена")

    # Шаг 2: Агрегатные графики
    if not args.skip_aggregate:
        logger.info("=== Шаг 2: Агрегатные графики распределений ===")
        agg_args = [str(data_dir), "-o", str(out_base / "aggregate_plots"), "--no-png"]
        script_path = Path(__file__).parent / "aggregate_plots.py"
        if not script_path.exists():
            logger.error(f"Не найден {script_path}, пропускаем")
        else:
            run_script(script_path, agg_args, "Агрегатные графики")
    else:
        logger.info("Агрегатные графики пропущены")

    # Шаг 3: Volcano + Strip plots
    if not args.skip_volcano:
        logger.info("=== Шаг 3: Volcano и Strip plots ===")
        volc_args = [
            str(data_dir), "-o", str(out_base / "volcano_strip"),
            "--effect-threshold", str(args.effect_threshold),
            "--no-png"
        ]
        script_path = Path(__file__).parent / "volcano_and_strip.py"
        if not script_path.exists():
            logger.error(f"Не найден {script_path}, пропускаем")
        else:
            run_script(script_path, volc_args, "Volcano/Strip")
    else:
        logger.info("Volcano/Strip plots пропущены")

    if not args.skip_cluster:
        logger.info("=== Шаг 4: Кластеризация SNV ===")
        for etype in ['Energy', 'EnergyLeft', 'EnergyRight']:
            if snv_df[etype].sum() == 0:
                logger.info(f"Тип {etype}: нет данных, пропускаем")
                continue
            # Дополнительная проверка: есть ли хотя бы два SNV с ненулевым числом записей
            snv_with_data = snv_df[snv_df[etype] > 0]['position'].values
            if len(snv_with_data) < 2:
                logger.warning(f"Тип {etype}: найдено всего {len(snv_with_data)} SNV с данными – кластеризация невозможна")
                continue
            logger.info(f"Запуск кластеризации для {etype} (SNV с данными: {len(snv_with_data)})")
            cluster_args = [
                str(data_dir), "-o", str(out_base / "clustering"),
                "--energy-type", etype,
                "--metric", "euclidean",
                "--max-constructs", str(args.max_constructs),
                "--no-png"
            ]
            script_path = Path(__file__).parent / "cluster_snvs.py"
            if not script_path.exists():
                logger.error(f"Не найден {script_path}, пропускаем")
                continue   # важно: перейти к следующему типу энергии, а не прерывать цикл break
            run_script(script_path, cluster_args, f"Кластеризация {etype}")
        else:
            logger.info("Кластеризация пропущена")

    logger.info("=== Пайплайн завершён ===")
    logger.info(f"Все результаты сохранены в {out_base}")

if __name__ == "__main__":
    main()

