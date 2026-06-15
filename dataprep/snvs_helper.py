import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def parse_snv_header(header: str) -> Optional[Tuple[int, str]]:
    if not header:
        return None
    match = re.search(r"SNV[_-](\d+)[_-]([ACGT])", header, re.IGNORECASE)
    if not match:
        return None
    try:
        position = int(match.group(1))
        if position <= 0:
            return None
        return position, match.group(2).upper()
    except ValueError:
        return None


def _locate_default_reference_fasta() -> Optional[Path]:
    candidates = [
        Path(__file__).resolve().parents[2] / 'sequences' / 'ref_seq' / 'Homo_sapiens_assembly38.chrM.fasta',
        Path(__file__).resolve().parents[3] / 'MitoFragilityScore' / 'Sequences' / 'Reference' / 'Homo_sapiens_assembly38.chrM.fasta',
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    search_root = Path(__file__).resolve().parents[3]
    for candidate in search_root.rglob('Homo_sapiens_assembly38.chrM.fasta'):
        return candidate

    return None


def load_reference_sequence(ref_fasta_path: Optional[str] = None) -> Dict[int, str]:
    if ref_fasta_path is None:
        ref_fasta_path = _locate_default_reference_fasta()
    if ref_fasta_path is None:
        logger.warning('Reference FASTA file not found for mtDNA lookup.')
        return {}

    ref_fasta_path = Path(ref_fasta_path)
    if not ref_fasta_path.exists():
        logger.warning(f'Reference FASTA path does not exist: {ref_fasta_path}')
        return {}

    sequence = []
    with ref_fasta_path.open('r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence.append(line.strip())
    sequence = ''.join(sequence).upper()
    return {i + 1: base for i, base in enumerate(sequence) if base in {'A', 'C', 'G', 'T'}}


def _normalize_snv_pair(position: int, ref: str, alt: str, ref_map: Dict[int, str]) -> Optional[Tuple[str, str]]:
    bases = {'A', 'C', 'G', 'T'}
    if len(ref) != 1 or len(alt) != 1:
        return None
    ref = ref.strip().upper()
    alt = alt.strip().upper()
    if ref not in bases or alt not in bases:
        return None
    genome_ref = ref_map.get(position)
    if genome_ref:
        if ref == genome_ref:
            pass
        elif alt == genome_ref:
            ref, alt = alt, ref
        else:
            logger.warning(
                f"SNV at position {position} is inconsistent with reference genome: csv ref={ref}, alt={alt}, genome={genome_ref}. Skipping."
            )
            return None
    if ref == alt:
        logger.warning(
            f"SNV at position {position} has identical ref/alt after normalization: {ref}->{alt}. Skipping."
        )
        return None
    return ref, alt


def load_fasta_snv_map(fastas_dir: Optional[str] = None, ref_fasta_path: Optional[str] = None) -> Dict[int, Tuple[str, str]]:
    fastas_dir = Path(fastas_dir or Path(__file__).parent.parent / 'fastas')
    ref_map = load_reference_sequence(ref_fasta_path)
    snv_map: Dict[int, Tuple[str, str]] = {}
    if not fastas_dir.exists():
        return snv_map

    for fasta_path in sorted(fastas_dir.glob('test_individual_*_snv_*.fasta')):
        try:
            with fasta_path.open('r', encoding='utf-8') as f:
                header = f.readline().strip()
                sequence = ''.join(line.strip() for line in f if line and not line.startswith('>')).upper()
        except Exception:
            continue
        if not header.startswith('>'):
            continue
        parsed = parse_snv_header(header)
        if parsed is None:
            continue
        pos, alt = parsed
        ref_allele = ref_map.get(pos)
        if not ref_allele or alt not in {'A', 'C', 'G', 'T'}:
            continue
        if ref_allele == alt:
            if len(sequence) >= pos:
                actual_alt = sequence[pos - 1]
                if actual_alt in {'A', 'C', 'G', 'T'} and actual_alt != ref_allele:
                    alt = actual_alt
                else:
                    logger.debug(
                        f"FASTA-derived SNV at position {pos} has ref==alt ({ref_allele}) and sequence does not show an alternative allele. Skipping."
                    )
                    continue
            else:
                logger.debug(
                    f"FASTA sequence too short to infer actual allele at position {pos} for {fasta_path}."
                )
                continue
        if ref_allele == alt:
            logger.debug(
                f"FASTA-derived SNV at position {pos} has ref==alt ({ref_allele}). Data input is invalid."
            )
            continue
        snv_map[pos] = (ref_allele, alt)
    return snv_map


def load_snv_nucleotide_data(snv_csv_path: Optional[str] = None, fastas_dir: Optional[str] = None, ref_fasta_path: Optional[str] = None) -> Dict[int, Tuple[str, str]]:
    snv_map: Dict[int, Tuple[str, str]] = {}
    if snv_csv_path is None:
        snv_csv_path = Path(__file__).parent.parent / 'snv_csv' / 'snvs.csv'

    ref_map = load_reference_sequence(ref_fasta_path)
    csv_counts: Dict[int, Counter[Tuple[str, str]]] = defaultdict(Counter)
    invalid_rows = 0

    try:
        df = pd.read_csv(snv_csv_path)
        for _, row in df.iterrows():
            try:
                position = int(row['position'])
                ref = str(row['ref_allele']).strip().upper()
                alt = str(row['alt_allele']).strip().upper()
                normalized = _normalize_snv_pair(position, ref, alt, ref_map)
                if normalized is not None:
                    csv_counts[position][normalized] += 1
                else:
                    invalid_rows += 1
            except Exception:
                invalid_rows += 1
                continue
    except Exception as e:
        logger.warning(f"Unable to read SNV CSV at {snv_csv_path}: {e}")

    for position, counter in csv_counts.items():
        if len(counter) > 1:
            logger.warning(
                f"Conflicting allele assignments for position {position}: {dict(counter)}. Using the most frequent pair."
            )
        canonical_pair, _ = counter.most_common(1)[0]
        snv_map[position] = canonical_pair

    if invalid_rows > 0:
        logger.info(f"Skipped {invalid_rows} invalid or inconsistent SNV rows from {snv_csv_path}")

    fasta_map = load_fasta_snv_map(fastas_dir, ref_fasta_path)
    # FASTA header is authoritative: override CSV-derived values when both exist
    for position, allele_pair in fasta_map.items():
        snv_map[position] = allele_pair

    return snv_map


def _find_fasta_for_master_json(master_json_path: str, fastas_dir: Optional[str] = None) -> Optional[Path]:
    fastas_dir = Path(fastas_dir or Path(__file__).parent.parent / 'fastas')
    if not fastas_dir.exists():
        return None

    dir_name = Path(master_json_path).parent.name
    # 1) Try exact mapping from metadata file
    info_path = fastas_dir / 'individuals_information.json'
    if info_path.exists():
        try:
            import json
            mapping = json.loads(info_path.read_text(encoding='utf-8'))
            if dir_name in mapping:
                fasta_name = mapping[dir_name][1]
                candidate = fastas_dir / fasta_name
                if candidate.exists():
                    return candidate
        except Exception:
            pass

    # 2) Try direct filename match by directory name
    candidate = fastas_dir / f"{dir_name}.fasta"
    if candidate.exists():
        return candidate

    # 3) Try to match by snv number in the directory name
    match = re.search(r'snv_(\d+)', dir_name, re.IGNORECASE)
    if match:
        target_pos = match.group(1)
        for fasta_path in fastas_dir.glob(f"*snv_{target_pos}*.fasta"):
            if fasta_path.exists():
                return fasta_path

    # 4) No direct candidate found
    return None


def infer_position_from_master_json_path(master_json_path: str) -> Optional[int]:
    stem = Path(master_json_path).stem
    match = re.search(r'snv[_-](\d+)', stem, re.IGNORECASE)
    if match:
        try:
            position = int(match.group(1))
            if position > 0:
                return position
            return None
        except ValueError:
            return None
    match = re.search(r'(\d+)_snv', stem, re.IGNORECASE)
    if match:
        try:
            position = int(match.group(1))
            if position > 0:
                return position
        except ValueError:
            pass
    return None


def infer_position_from_fasta_for_master_json(master_json_path: str, fastas_dir: Optional[str] = None) -> Optional[int]:
    fasta_path = _find_fasta_for_master_json(master_json_path, fastas_dir)
    if fasta_path is None:
        return infer_position_from_master_json_path(master_json_path)
    try:
        with fasta_path.open('r', encoding='utf-8') as f:
            header = f.readline().strip()
    except Exception:
        return infer_position_from_master_json_path(master_json_path)
    parsed = parse_snv_header(header)
    if parsed is not None:
        return parsed[0]
    return infer_position_from_master_json_path(master_json_path)


def infer_alt_allele_for_master_json(master_json_path: str, fastas_dir: Optional[str] = None) -> Optional[str]:
    fasta_path = _find_fasta_for_master_json(master_json_path, fastas_dir)
    if fasta_path is None:
        return None
    try:
        with fasta_path.open('r', encoding='utf-8') as f:
            header = f.readline().strip()
    except Exception:
        return None
    parsed = parse_snv_header(header)
    return parsed[1] if parsed is not None else None


def infer_snv_pair_for_master_json(master_json_path: str, fastas_dir: Optional[str] = None, ref_fasta_path: Optional[str] = None) -> Optional[Tuple[int, str, str]]:
    position = infer_position_from_fasta_for_master_json(master_json_path, fastas_dir)
    alt_allele = infer_alt_allele_for_master_json(master_json_path, fastas_dir)

    if position is not None and alt_allele is None:
        snv_map = load_snv_nucleotide_data(None, fastas_dir, ref_fasta_path)
        if position in snv_map:
            alt_allele = snv_map[position][1]

    if position is None or alt_allele is None:
        return None

    ref_map = load_reference_sequence(ref_fasta_path)
    ref_allele = ref_map.get(position)
    if ref_allele is None:
        logger.warning(f"Reference base not found for position {position} while parsing FASTA header for {master_json_path}")
        return None

    ref_allele = ref_allele.strip().upper()
    alt_allele = alt_allele.strip().upper()
    if ref_allele == alt_allele:
        logger.warning(
            f"FASTA-derived SNV at position {position} has identical ref/alt alleles ({ref_allele}); this pair will be ignored for {master_json_path}.")
        return None

    return position, ref_allele, alt_allele


def resolve_point_position(point: dict, master_json_path: str, fastas_dir: Optional[str] = None) -> Optional[int]:
    # Check point-level SNP annotation first, then fall back to FASTA header lookup.
    for key in ('snp', 'position', 'pos'):
        value = point.get(key)
        if value is None:
            continue
        if isinstance(value, str):
            norm = value.strip().lower()
            if norm in {'', 'none', 'n/a', 'na'}:
                continue
        try:
            position = int(value)
            if position > 0:
                return position
        except (TypeError, ValueError):
            continue

    return infer_position_from_fasta_for_master_json(master_json_path, fastas_dir)
