from typing import Optional, Tuple
import re
from logging import getLogger

logger = getLogger(__name__)


def parse_construct_id(construct_id: str) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
    """Parse construct id and return (arm_size, center, arm3_start, arm4_start) or Nones on failure."""
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
    m = re.search(r'CGS-(\d+)-\d+-\d+-\d+-\d+-(\d+)', construct_id)
    if m:
        return int(m.group(2))
    logger.warning(f"Не найден блок CGS в ID: {construct_id}")
    return None


def _parse_cen(construct_id: str) -> Optional[int]:
    m = re.search(r'CEN-(\d+)', construct_id)
    if m:
        return int(m.group(1))
    logger.warning(f"Не найден блок CEN в ID: {construct_id}")
    return None


def _parse_con(construct_id: str) -> Tuple[Optional[int], Optional[int]]:
    m = re.search(r'CON-(\d+)-(?P<arm4>\d+)', construct_id)
    if m:
        return int(m.group(1)), int(m.group('arm4'))
    logger.warning(f"Не найден блок CON в ID: {construct_id}")
    return None, None
