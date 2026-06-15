"""Фильтрация и парсинг SNV из различных источников."""

from typing import List, Tuple, Optional, Set, Dict, Any
from pathlib import Path
import pandas as pd
import json
from logging import getLogger

logger = getLogger(__name__)


def filter_snv_by_fdr(snv_df: pd.DataFrame, fdr_threshold: float = 0.056) -> pd.DataFrame:
    """Фильтрация SNV по FDR.

    Args:
        snv_df: DataFrame с SNV.
        fdr_threshold: Порог FDR.

    Returns:
        Отфильтрованный DataFrame.
    """
    return snv_df[snv_df['FDR'] < fdr_threshold].copy()


def filter_by_positions(snv_df: pd.DataFrame, position_filter: str) -> pd.DataFrame:
    """Фильтрация SNV по позициям.

    Args:
        snv_df: DataFrame с SNV.
        position_filter: Позиции через запятую (например, "1000,2000,3000").

    Returns:
        Отфильтрованный DataFrame.
    """
    try:
        positions = [int(pos.strip()) for pos in position_filter.split(',')]
        return snv_df[snv_df['Position'].isin(positions)].copy()
    except Exception as e:
        logger.error(f"Ошибка фильтрации по позициям: {e}")
        return snv_df


def filter_by_region(snv_df: pd.DataFrame, region_filter: str) -> pd.DataFrame:
    """Фильтрация SNV по геномной области.

    Args:
        snv_df: DataFrame с SNV.
        region_filter: Область в формате "chr:start-end".

    Returns:
        Отфильтрованный DataFrame.
    """
    try:
        region_parts = region_filter.split(':')
        if len(region_parts) != 2:
            logger.warning(f"Неверный формат региона: {region_filter}")
            return snv_df
        
        chrom = region_parts[0]
        start, end = map(int, region_parts[1].split('-'))
        
        # Предполагаем, что есть столбец 'Chromosome' и 'Position'
        filtered = snv_df[
            (snv_df.get('Chromosome') == chrom) &
            (snv_df['Position'] >= start) &
            (snv_df['Position'] <= end)
        ]
        return filtered.copy()
    except ValueError as e:
        logger.error(f"Ошибка парсинга региона '{region_filter}': {e}")
        return snv_df


def filter_by_trait(snv_df: pd.DataFrame, trait_filter: str) -> pd.DataFrame:
    """Фильтрация SNV по признаку/черте.

    Args:
        snv_df: DataFrame с SNV.
        trait_filter: Название признака для поиска.

    Returns:
        Отфильтрованный DataFrame.
    """
    if 'Trait' not in snv_df.columns:
        logger.warning("Столбец 'Trait' не найден в DataFrame")
        return snv_df
    return snv_df[snv_df['Trait'].str.contains(trait_filter, case=False, na=False)].copy()


def extract_alleles(snv_df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """Извлечение референсного и альтернативного аллелей.

    Args:
        snv_df: DataFrame с SNV (должен содержать столбцы с нуклеотидами).

    Returns:
        Кортеж (список референсных аллелей, список альтернативных аллелей).
    """
    ref_alleles: List[str] = []
    alt_alleles: List[str] = []
    
    for _, row in snv_df.iterrows():
        # Поиск столбцов с аллелями (могут быть разные названия)
        ref = None
        alt = None
        
        # Варианты названий для референсного аллеля
        for col in ['ref', 'ref_allele', 'Ref', 'Ref_allele', 'REF']:
            if col in row and pd.notna(row[col]):
                ref = str(row[col]).upper()
                break
        
        # Варианты названий для альтернативного аллеля
        for col in ['alt', 'alt_allele', 'Alt', 'Alt_allele', 'ALT']:
            if col in row and pd.notna(row[col]):
                alt = str(row[col]).upper()
                break
        
        # Если не нашли явно, пытаемся использовать первую и последнюю букву
        if ref is None or alt is None:
            # Fallback: используем столбец 'Alleles' или похожий
            for col in ['Alleles', 'alleles', 'variants', 'Variants']:
                if col in row and pd.notna(row[col]):
                    alleles_str = str(row[col]).upper()
                    parts = alleles_str.split('/')
                    if len(parts) >= 2:
                        ref = ref or parts[0]
                        alt = alt or parts[1]
                    break
        
        ref_alleles.append(ref or 'N')
        alt_alleles.append(alt or 'N')
    
    return ref_alleles, alt_alleles


def get_covered_positions(settings_path: Path, ref_seq_info_path: Path) -> Set[int]:
    """Получить позиции, покрытые конструктами из settings.json и reference_sequence_information.json.

    Args:
        settings_path: Путь к settings.json.
        ref_seq_info_path: Путь к reference_sequence_information.json.

    Returns:
        Множество покрытых позиций.
    """
    covered_positions: Set[int] = set()
    
    try:
        # Загрузим settings.json
        with open(settings_path, 'r', encoding='utf-8') as f:
            settings = json.load(f)
        
        # Загрузим reference_sequence_information.json
        with open(ref_seq_info_path, 'r', encoding='utf-8') as f:
            ref_info = json.load(f)
        
        # Из settings берём регионы, из ref_info берём координаты
        # Это зависит от структуры JSON - упрощённо берём все числовые диапазоны
        if isinstance(settings, dict):
            for key, value in settings.items():
                if isinstance(value, dict):
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, (int, float)):
                            covered_positions.add(int(sub_value))
        
        logger.info(f"Загружено {len(covered_positions)} покрытых позиций")
    except Exception as e:
        logger.error(f"Ошибка загрузки покрытых позиций: {e}")
    
    return covered_positions
