"""Построение и применение мутаций в последовательности."""

from typing import List, Tuple, Optional, Dict, Any
from pathlib import Path
from logging import getLogger

logger = getLogger(__name__)

# Генетический код (ДНК -> Белок)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def apply_mutation(seq: str, position: int, ref: str, alt: str) -> Optional[str]:
    """Применить мутацию к последовательности.

    Args:
        seq: Исходная последовательность.
        position: Позиция мутации (0-indexed).
        ref: Референсный аллель.
        alt: Альтернативный аллель.

    Returns:
        Мутированная последовательность или None в случае ошибки.
    """
    try:
        if position < 0 or position >= len(seq):
            logger.warning(f"Позиция {position} вне границ последовательности длины {len(seq)}")
            return None
        
        # Проверим, что на позиции находится ожидаемый аллель
        if seq[position].upper() != ref.upper():
            logger.warning(
                f"На позиции {position} найден '{seq[position].upper()}', "
                f"а не ожидаемый '{ref.upper()}'"
            )
            # Всё равно применим мутацию
        
        # Применим мутацию
        mutated = seq[:position] + alt.upper() + seq[position + 1:]
        return mutated
    except Exception as e:
        logger.error(f"Ошибка применения мутации на позиции {position}: {e}")
        return None


def apply_multiple_mutations(
    seq: str,
    mutations: List[Tuple[int, str, str]]
) -> Tuple[Optional[str], List[int]]:
    """Применить несколько мутаций к последовательности.

    Args:
        seq: Исходная последовательность.
        mutations: Список кортежей (позиция, ref, alt).

    Returns:
        Кортеж (мутированная последовательность, список индексов не применённых мутаций).
    """
    current_seq = seq
    failed_indices = []
    
    for i, (pos, ref, alt) in enumerate(mutations):
        result = apply_mutation(current_seq, pos, ref, alt)
        if result is None:
            failed_indices.append(i)
        else:
            current_seq = result
    
    if failed_indices:
        logger.warning(f"Не удалось применить {len(failed_indices)} мутаций")
    
    return current_seq, failed_indices


def translate_dna_to_protein(seq: str) -> str:
    """Перевести ДНК последовательность в белок.

    Args:
        seq: ДНК последовательность.

    Returns:
        Аминокислотная последовательность.
    """
    seq_upper = seq.upper()
    
    # Если длина не кратна 3, логируем предупреждение
    if len(seq_upper) % 3 != 0:
        logger.warning(f"Длина последовательности {len(seq_upper)} не кратна 3")
    
    protein = []
    for i in range(0, len(seq_upper) - 2, 3):
        codon = seq_upper[i:i + 3]
        amino_acid = GENETIC_CODE.get(codon, 'X')  # X для неизвестных кодонов
        protein.append(amino_acid)
    
    return ''.join(protein)


def calculate_mutation_impact(ref_aa: str, alt_aa: str) -> str:
    """Определить тип мутации на уровне белка.

    Args:
        ref_aa: Референсная аминокислота.
        alt_aa: Альтернативная аминокислота.

    Returns:
        Тип мутации: 'synonymous', 'missense', 'nonsense' или 'unknown'.
    """
    if ref_aa == alt_aa:
        return 'synonymous'
    elif alt_aa == '*':
        return 'nonsense'
    elif ref_aa == '*':
        return 'unknown'  # Редкий случай
    else:
        return 'missense'


def extract_region(seq: str, start: int, end: int) -> Optional[str]:
    """Извлечь регион из последовательности.

    Args:
        seq: Полная последовательность.
        start: Начало региона (0-indexed).
        end: Конец региона (exclusive).

    Returns:
        Регион последовательности или None.
    """
    try:
        if start < 0 or end > len(seq) or start >= end:
            logger.warning(f"Неверные границы: start={start}, end={end}, len={len(seq)}")
            return None
        return seq[start:end]
    except Exception as e:
        logger.error(f"Ошибка извлечения региона: {e}")
        return None


def calculate_gc_content(seq: str) -> float:
    """Рассчитать содержание G+C в последовательности.

    Args:
        seq: ДНК последовательность.

    Returns:
        Процент содержания G+C.
    """
    if not seq:
        return 0.0
    
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    return (gc_count / len(seq_upper)) * 100


def reverse_complement(seq: str) -> str:
    """Получить комплементарную обратную последовательность.

    Args:
        seq: ДНК последовательность.

    Returns:
        Комплементарная обратная последовательность.
    """
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                      'N': 'N', 'n': 'n'}
    
    try:
        complement = ''.join(complement_map.get(base, 'N') for base in seq)
        return complement[::-1]
    except Exception as e:
        logger.error(f"Ошибка вычисления комплемента: {e}")
        return seq
