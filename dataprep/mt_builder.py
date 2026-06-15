"""Главная логика для построения mt_DNA последовательностей с SNVs."""

from typing import Dict, List, Optional, Set, Any, Tuple
from pathlib import Path
from logging import getLogger
import json

import pandas as pd
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

from .snv_filter import (
    filter_snv_by_fdr,
    filter_by_positions,
    filter_by_region,
    filter_by_trait,
    extract_alleles
)
from .excel_parser import load_excel_file
from .mt_output import save_individual_results

logger = getLogger(__name__)


def get_covered_positions_from_settings(
    settings_path: Path,
    ref_seq_info_path: Path
) -> Set[int]:
    """Получить позиции, покрытые конструктами из settings.json.

    Args:
        settings_path: Путь к settings.json.
        ref_seq_info_path: Путь к reference_sequence_information.json.

    Returns:
        Множество покрытых позиций.
    """
    try:
        # Загружаем settings.json
        with open(settings_path, 'r', encoding='utf-8') as f:
            settings = json.load(f)
        
        ref_seq_key = settings["ReferenceSequence"][0]
        subseq_key = settings["ReferenceSequence"][1]
        
        logger.info(f"Using reference sequence: {ref_seq_key}, sub-sequence: {subseq_key}")
        
        # Загружаем reference_sequence_information.json
        with open(ref_seq_info_path, 'r', encoding='utf-8') as f:
            ref_seq_info = json.load(f)
        
        # Получаем регион из SubSequences
        region = ref_seq_info[ref_seq_key]["SubSequences"][subseq_key]
        start_pos, end_pos = region
        
        logger.info(f"Region from settings: {start_pos} to {end_pos}")
        
        # Генерируем все позиции в регионе
        covered_positions = set(range(start_pos, end_pos + 1))
        logger.info(f"Generated {len(covered_positions)} covered positions")
        
        return covered_positions
    except Exception as e:
        logger.error(f"Error reading settings or reference sequence info: {e}")
        return set()


def csv_constructor(
    excel_path: Path,
    output_path: Path,
    position_filter: Optional[str] = None,
    region_filter: Optional[str] = None,
    fdr_threshold: float = 0.056,
    trait_filter: Optional[str] = None
) -> pd.DataFrame:
    """Создать CSV файл с SNVs из Excel с различными фильтрами.

    Args:
        excel_path: Путь к входному Excel файлу.
        output_path: Путь к выходному CSV файлу.
        position_filter: Позиции через запятую для фильтрации.
        region_filter: Регион для фильтрации (например, "chrM:1000-2000").
        fdr_threshold: Порог FDR для фильтрации.
        trait_filter: Название признака для фильтрации.

    Returns:
        DataFrame с отфильтрованными SNVs.
    """
    snv_df = load_excel_file(excel_path)
    if snv_df is None:
        logger.error(f"Failed to load Excel file: {excel_path}")
        return pd.DataFrame()
    
    filtered_snv_df = filter_snv_by_fdr(snv_df, fdr_threshold)
    
    # Применяем дополнительные фильтры
    if position_filter:
        filtered_snv_df = filter_by_positions(filtered_snv_df, position_filter)
    if region_filter:
        filtered_snv_df = filter_by_region(filtered_snv_df, region_filter)
    if trait_filter:
        filtered_snv_df = filter_by_trait(filtered_snv_df, trait_filter)
    
    # Извлекаем аллели
    try:
        ref_alleles, alt_alleles = _extract_alleles_from_dataframe(filtered_snv_df)
        filtered_snv_df['ref_allele'] = ref_alleles
        filtered_snv_df['alt_allele'] = alt_alleles
    except Exception as e:
        logger.error(f"Error extracting alleles: {e}")
        return pd.DataFrame()
    
    final_df = filtered_snv_df[['Position', 'ref_allele', 'alt_allele']].rename(
        columns={'Position': 'position'}
    )
    
    # Сохраняем CSV
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        final_df.to_csv(output_path, index=False)
        logger.info(f"Saved {len(final_df)} unique SNVs to {output_path}")
    except Exception as e:
        logger.error(f"Error saving CSV: {e}")
    
    return final_df


def _extract_alleles_from_dataframe(
    filtered_snv_df: pd.DataFrame
) -> Tuple[List[str], List[str]]:
    """Извлечь референсный и альтернативный аллели из DataFrame.

    Args:
        filtered_snv_df: Отфильтрованный DataFrame с SNVs.

    Returns:
        Кортеж (список референсных аллелей, список альтернативных аллелей).
    """
    ref_alleles: List[str] = []
    alt_alleles: List[str] = []
    
    for _, row in filtered_snv_df.iterrows():
        alt_allele = row['Minor allele']
        alt_alleles.append(alt_allele)
        
        if row['Allele1'] == alt_allele:
            ref_alleles.append(row['Allele2'])
        elif row['Allele2'] == alt_allele:
            ref_alleles.append(row['Allele1'])
        else:
            raise ValueError(
                f"Error at position {row['Position']}: "
                f"Minor allele '{alt_allele}' does not match "
                f"either Allele1 ({row['Allele1']}) or Allele2 ({row['Allele2']})"
            )
    
    return ref_alleles, alt_alleles


def apply_single_snv(
    ref_record: SeqRecord,
    snv: Dict[str, Any],
    log_path: Path,
    covered_positions: Set[int],
    snv_index: int
) -> Optional[SeqRecord]:
    """Применить одиночную SNV к референсной последовательности.

    Args:
        ref_record: SeqRecord с референсной последовательностью.
        snv: Словарь SNV с position, ref_allele, alt_allele.
        log_path: Путь к файлу логирования мутаций.
        covered_positions: Множество позиций, покрытых конструктами.
        snv_index: Индекс SNV для именования вывода.

    Returns:
        Изменённый SeqRecord или None если не удалось.
    """
    original_seq: str = str(ref_record.seq).upper()
    mutable_seq: MutableSeq = MutableSeq(original_seq)
    
    position = snv['position']
    ref_allele = snv['ref_allele'].upper()
    alt_allele = snv['alt_allele'].upper()
    
    # Проверяем, что позиция покрыта
    if position not in covered_positions:
        logger.warning(f"Position {position} not covered by analysis region. Skipping.")
        return None
    
    # Применяем SNV
    applied, mismatch_log = _apply_single_snv_internal(
        mutable_seq, original_seq, position, ref_allele, alt_allele
    )
    
    _save_mismatch_log(mismatch_log, log_path)
    
    if applied:
        return SeqRecord(
            seq=mutable_seq,
            id=f"custom_mtDNA_snv_{snv_index}",
            description=(
                f"Modified from {ref_record.id} | "
                f"SNV: {position} {ref_allele}->{alt_allele}"
            )
        )
    return None


def _apply_single_snv_internal(
    mutable_seq: MutableSeq,
    original_seq: str,
    position: int,
    ref_allele: str,
    alt_allele: str
) -> Tuple[bool, List[Dict[str, Any]]]:
    """Применить одиночную SNV к mutable последовательности.

    Args:
        mutable_seq: Mutable последовательность для изменения.
        original_seq: Строка оригинальной последовательности.
        position: Позиция для мутации.
        ref_allele: Референсный аллель.
        alt_allele: Альтернативный аллель.

    Returns:
        Кортеж (успех, лог несоответствий).
    """
    mismatch_log: List[Dict[str, Any]] = []
    idx = position - 1
    
    if idx >= len(original_seq):
        mismatch_log.append({
            'position': position,
            'original_base': None,
            'ref_allele': ref_allele,
            'alt_allele': alt_allele,
            'status': 'SKIPPED',
            'notes': 'Position out of sequence bounds'
        })
        return False, mismatch_log
    
    current_base = original_seq[idx]
    
    if ref_allele == alt_allele:
        mismatch_log.append({
            'position': position,
            'original_base': current_base,
            'ref_allele': ref_allele,
            'alt_allele': alt_allele,
            'status': 'SKIPPED',
            'notes': 'Ref and Alt alleles are identical'
        })
        return False, mismatch_log
    
    if current_base == ref_allele:
        mutable_seq[idx] = alt_allele
        mismatch_log.append({
            'position': position,
            'original_base': current_base,
            'ref_allele': ref_allele,
            'alt_allele': alt_allele,
            'status': 'APPLIED',
            'notes': ''
        })
        return True, mismatch_log
    else:
        notes = []
        if current_base == alt_allele:
            notes.append('ALT allele already present')
        else:
            notes.append(f"Expected ref: {ref_allele}, found: {current_base}")
        mismatch_log.append({
            'position': position,
            'original_base': current_base,
            'ref_allele': ref_allele,
            'alt_allele': alt_allele,
            'status': 'SKIPPED',
            'notes': '; '.join(notes)
        })
        return False, mismatch_log


def _save_mismatch_log(mismatch_log: List[Dict[str, Any]], log_path: Path) -> None:
    """Сохранить лог несоответствий в CSV файл.

    Args:
        mismatch_log: Список записей логирования.
        log_path: Путь к файлу логирования.
    """
    if mismatch_log:
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            mismatch_df = pd.DataFrame(mismatch_log)
            mismatch_df.to_csv(log_path, index=False)
            logger.info(f"Mutation log saved to {log_path}")
        except Exception as e:
            logger.error(f"Error saving mismatch log: {e}")


def process_snvs(
    snv_df: pd.DataFrame,
    input_fasta: Path,
    output_dir: Path,
    settings_path: Path,
    ref_seq_info_path: Path,
    run_id: int
) -> None:
    """Обработать SNVs и генерировать последовательности.

    Args:
        snv_df: DataFrame с SNVs.
        input_fasta: Путь к входному FASTA файлу.
        output_dir: Выходная директория для последовательностей.
        settings_path: Путь к settings.json.
        ref_seq_info_path: Путь к reference_sequence_information.json.
        run_id: Run ID для именования.
    """
    try:
        ref_record = SeqIO.read(input_fasta, "fasta")
        logger.info(f"Loaded reference sequence: {ref_record.id}")
        logger.info(f"Length: {len(ref_record.seq)} bp")
    except Exception as e:
        logger.error(f"Error loading reference sequence: {e}")
        return
    
    # Получаем покрытые позиции из settings
    covered_positions = get_covered_positions_from_settings(settings_path, ref_seq_info_path)
    
    if not covered_positions:
        logger.error("No covered positions found. Cannot generate sequences.")
        return
    
    # Создаём выходную директорию
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Error creating output directory: {e}")
        return
    
    # Обрабатываем каждую SNV индивидуально
    successful_count = 0
    for i, (_, snv_row) in enumerate(snv_df.iterrows()):
        snv = snv_row.to_dict()
        log_path = output_dir / f"snv_log_{run_id}_{i}.csv"
        output_fasta = output_dir / f"test_individual_{run_id}_snv_{i}.fasta"
        
        custom_record = apply_single_snv(
            ref_record, snv, log_path, covered_positions, i
        )
        
        if custom_record:
            try:
                SeqIO.write(custom_record, output_fasta, "fasta")
                successful_count += 1
                logger.info(f"Created sequence for SNV {i}: {output_fasta}")
            except Exception as e:
                logger.error(f"Error writing sequence for SNV {i}: {e}")
        else:
            logger.warning(f"Failed to apply SNV {i}: position {snv['position']}")
    
    logger.info(f"Successfully created {successful_count} out of {len(snv_df)} sequences")
