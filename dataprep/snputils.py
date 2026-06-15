from typing import List, Tuple, Set, Optional
from pathlib import Path
from logging import getLogger

logger = getLogger(__name__)


def load_snp_data(snp_file_path: str) -> Set[int]:
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


def calculate_arm_ranges(arm_size: int, center: int, arm3_start: int, arm4_start: int, subseq_start: int = 10) -> List[Tuple[int, int]]:
    arm1_start = center - arm_size
    arm2_start = center
    ranges = [
        (arm1_start + subseq_start, arm1_start + arm_size + subseq_start),
        (arm2_start + subseq_start, arm2_start + arm_size + subseq_start),
        (arm3_start + subseq_start, arm3_start + arm_size + subseq_start),
        (arm4_start + subseq_start, arm4_start + arm_size + subseq_start)
    ]
    logger.debug(f"Диапазоны конструкта (абсолютные): {ranges}")
    return ranges


def get_snps_in_construct(snp_positions: Set[int], arm_size: int, center: int, arm3_start: int, arm4_start: int, subseq_start: int = 10) -> List[int]:
    ranges = calculate_arm_ranges(arm_size, center, arm3_start, arm4_start, subseq_start)
    snps_in_construct = []
    for snp in snp_positions:
        for start, end in ranges:
            if start <= snp <= end:
                snps_in_construct.append(snp)
                break
    logger.debug(f"SNP в конструкте: {snps_in_construct}")
    return snps_in_construct
