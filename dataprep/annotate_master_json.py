#!/usr/bin/env python3
"""
Скрипт для аннотации master_data.json нуклеотидными переходами.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple, List

try:
    from .snvs_helper import parse_snv_header, load_reference_sequence
except ImportError:
    from snvs_helper import parse_snv_header, load_reference_sequence

logger = logging.getLogger(__name__)

def setup_logger(level=logging.INFO):
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

def find_master_json_files(data_dir: Path, recursive: bool = True) -> List[Path]:
    pattern = "*_master_data.json"
    return list(data_dir.rglob(pattern)) if recursive else list(data_dir.glob(pattern))

def find_fasta_for_master_json(master_path: Path, fastas_dir: Path) -> Optional[Path]:
    # Приоритет: имя родительской папки (обычно test_individual_1_snv_XXX)
    parent_name = master_path.parent.name
    if "snv" in parent_name.lower():
        candidate = fastas_dir / f"{parent_name}.fasta"
        if candidate.exists():
            return candidate
    # Fallback: имя файла без суффикса _master_data
    stem = master_path.stem
    if stem.endswith("_master_data"):
        fasta_stem = stem[:-12]
    else:
        fasta_stem = stem
    candidate = fastas_dir / f"{fasta_stem}.fasta"
    if candidate.exists():
        return candidate
    # Поиск по маске (на случай разных расширений)
    for f in fastas_dir.glob(f"*{fasta_stem}*.*"):
        if f.suffix.lower() in {".fasta", ".fa", ".fna"}:
            return f
    logger.warning(f"FASTA не найден для {master_path.name} (искали {candidate})")
    return None

def parse_fasta_header(fasta_path: Path) -> Optional[Tuple[int, str]]:
    try:
        with open(fasta_path, "r", encoding="utf-8") as f:
            header = f.readline().strip()
            if not header.startswith(">"):
                return None
            return parse_snv_header(header)
    except Exception as e:
        logger.error(f"Ошибка чтения {fasta_path}: {e}")
        return None

def update_master_json(master_path: Path, ref_map: Dict[int, str], fastas_dir: Path, force: bool = False) -> bool:
    try:
        with open(master_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        logger.error(f"Не удалось прочитать {master_path}: {e}")
        return False
    if not isinstance(data, list):
        logger.warning(f"Неожиданный формат {master_path}")
        return False

    fasta_path = find_fasta_for_master_json(master_path, fastas_dir)
    if fasta_path is None:
        return False

    header_info = parse_fasta_header(fasta_path)
    if header_info is None:
        logger.warning(f"Не удалось распарсить заголовок {fasta_path}")
        return False
    position, alt_allele = header_info

    ref_allele = ref_map.get(position)
    if ref_allele is None:
        logger.error(f"Референсный аллель для позиции {position} не найден")
        return False

    modified = False
    for record in data:
        if force or record.get("snp") is None:
            record["snp"] = position
            record["ref_allele"] = ref_allele
            record["alt_allele"] = alt_allele
            record["mutation_code"] = f"{position}{ref_allele}>{alt_allele}"
            modified = True
    if modified:
        with open(master_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        logger.info(f"Обновлён {master_path}: {position}{ref_allele}>{alt_allele}")
    else:
        logger.debug(f"В {master_path} уже есть аннотации (используйте --force для перезаписи)")
    return modified

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("--fastas-dir", type=Path, default=None)
    parser.add_argument("--ref-fasta", type=Path, default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-recursive", action="store_true")
    args = parser.parse_args()
    setup_logger()
    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Директория не существует: {data_dir}")
        sys.exit(1)

    # Определение fastas_dir
    if args.fastas_dir:
        fastas_dir = args.fastas_dir.resolve()
    else:
        candidates = [data_dir.parent / "fastas", Path(__file__).resolve().parent / "fastas"]
        fastas_dir = None
        for cand in candidates:
            if cand.exists():
                fastas_dir = cand
                break
        if fastas_dir is None:
            logger.error("Не найдена папка fastas, укажите --fastas-dir")
            sys.exit(1)

    # Референсный FASTA
    if args.ref_fasta:
        ref_fasta_path = args.ref_fasta.resolve()
    else:
        candidates = [
            Path(__file__).resolve().parent / "sequences" / "ref_seq" / "Homo_sapiens_assembly38.chrM.fasta",
            data_dir.parent / "sequences" / "ref_seq" / "Homo_sapiens_assembly38.chrM.fasta",
        ]
        ref_fasta_path = None
        for cand in candidates:
            if cand.exists():
                ref_fasta_path = cand
                break
        if ref_fasta_path is None:
            logger.error("Референсный FASTA не найден, укажите --ref-fasta")
            sys.exit(1)

    logger.info(f"Загрузка референса из {ref_fasta_path}")
    ref_map = load_reference_sequence(str(ref_fasta_path))
    if not ref_map:
        logger.error("Не удалось загрузить референс")
        sys.exit(1)

    json_files = find_master_json_files(data_dir, recursive=not args.no_recursive)
    if not json_files:
        logger.error(f"Нет файлов *_master_data.json в {data_dir}")
        sys.exit(1)
    logger.info(f"Найдено {len(json_files)} файлов")

    updated = 0
    for master_path in json_files:
        if update_master_json(master_path, ref_map, fastas_dir, force=args.force):
            updated += 1
    logger.info(f"Обновлено: {updated} из {len(json_files)}")

if __name__ == "__main__":
    main()