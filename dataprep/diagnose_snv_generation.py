#!/usr/bin/env python3
"""
Диагностика генерации SNV: сравнивает Excel с реальными JSON.
Выводит причины отсутствия SNV.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

MT_GENE_RANGES = {  # для справки
    'ND1': (3307, 4262), 'ND2': (4470, 5511), 'CO1': (5904, 7445),
    'CO2': (7586, 8269), 'ATP8': (8366, 8572), 'ATP6': (8527, 9207),
    'CO3': (9207, 9990), 'ND3': (10059, 10404), 'ND4L': (10470, 10766),
    'ND4': (10760, 12141), 'ND5': (12337, 14148), 'ND6': (14149, 14673),
    'CYB': (14747, 15887),
}

def load_region(settings_path, ref_info_path):
    """Загрузить покрываемый регион из файлов настроек."""
    try:
        with open(settings_path) as f:
            settings = json.load(f)
        ref_key = settings["ReferenceSequence"][0]
        sub_key = settings["ReferenceSequence"][1]
        with open(ref_info_path) as f:
            ref_info = json.load(f)
        return ref_info[ref_key]["SubSequences"][sub_key]
    except Exception:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path, help="Папка с подпапками test_individual_1_snv_* (например, ../output/all_syn)")
    parser.add_argument("--excel", type=Path, default="D:/pythonProject/MitoFragility/DataPreparing/raw_data/MitoPhewas_associations.xlsx")
    parser.add_argument("--settings", type=Path, default="D:/pythonProject/MitoFragility/MitoFragilityScore/Settings.json")
    parser.add_argument("--ref-info", type=Path, default="D:/pythonProject/MitoFragility/MitoFragilityScore/ReferenceSequenceInformation.json")
    args = parser.parse_args()

    result_dir = args.result_dir.resolve()
    if not result_dir.is_dir():
        print(f"Ошибка: папка {result_dir} не существует")
        return

    # Загружаем Excel
    if not args.excel.exists():
        print(f"Excel не найден: {args.excel}")
        return
    df = pd.read_excel(args.excel)
    df = df[df['FDR'] < 0.056]
    positions_excel = set(df['Position'].unique())
    print(f"Всего SNV после FDR: {len(positions_excel)}")

    # Собираем существующие непустые SNV из JSON
    existing = set()
    json_files = list(result_dir.rglob("*_master_data.json"))
    print(f"Найдено JSON-файлов: {len(json_files)}")
    for jf in json_files:
        if jf.stat().st_size < 100:
            continue
        try:
            with open(jf) as f:
                data = json.load(f)
            if data:
                pos = data[0].get('snp') or data[0].get('position')
                if pos:
                    existing.add(int(pos))
        except Exception:
            pass
    print(f"Непустых SNV: {len(existing)}")

    missing = positions_excel - existing
    print(f"Отсутствуют {len(missing)} позиций.")
    if not missing:
        print("Все SNV на месте.")
        return

    # Загружаем регион, если возможно
    region = load_region(args.settings, args.ref_info)
    if region:
        print(f"Покрываемый регион конструктов: {region[0]} - {region[1]}")
    else:
        print("Не удалось загрузить регион (файлы настроек не найдены или повреждены)")

    # Анализируем отсутствующие
    print("\nАнализ первых 10 отсутствующих SNV:")
    for pos in sorted(missing)[:100]:
        row = df[df['Position'] == pos].iloc[0]
        fdr = row['FDR']
        minor = row['Minor allele']
        a1 = row['Allele1']
        a2 = row['Allele2']
        in_region = region and (region[0] <= pos <= region[1]) if region else "N/A"
        print(f"  {pos}: FDR={fdr:.4f}, minor={minor}, alleles={a1}/{a2}, in_region={in_region}")

if __name__ == "__main__":
    main()