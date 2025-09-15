
import os
import random
import logging
from pathlib import Path
from typing import Tuple, List, Set, Dict, Any, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import importlib.util



def import_scatter_module() -> Any:
    """
    Dynamically imports the scatter_plus_n_std.py module from the plots directory.

    Returns:
        The imported module object.
    Raises:
        FileNotFoundError: If the module file does not exist.
        ImportError: If the module cannot be loaded.
    """
    plots_path = Path(__file__).resolve().parent.parent / "plots"
    scatter_path = plots_path / "scatter_plus_n_std.py"
    if not scatter_path.exists():
        raise FileNotFoundError(f"Module not found: {scatter_path}")
    spec = importlib.util.spec_from_file_location("scatter_plus_n_std", scatter_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load spec for {scatter_path}")
    scatter_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scatter_mod)
    return scatter_mod


scatter_mod = import_scatter_module()
parse_construct_id = scatter_mod.parse_construct_id
calculate_arm_ranges = scatter_mod.calculate_arm_ranges



def setup_logger() -> logging.Logger:
    """
    Configures and returns a logger for this script.

    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger("mt_DNA_builder")
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


logger = setup_logger()



def csv_constructor(excel_path: Path, output_path: Path) -> pd.DataFrame:
    """
    Creates a CSV file with SNVs from an Excel file.

    Args:
        excel_path (Path): Path to the input Excel file.
        output_path (Path): Path to the output CSV file.

    Returns:
        pd.DataFrame: DataFrame with filtered SNVs.
    """
    snv_df = pd.read_excel(excel_path)
    filtered_snv_df = _filter_snv_by_fdr(snv_df)
    ref_alleles, alt_alleles = _extract_alleles(filtered_snv_df)
    filtered_snv_df['ref_allele'] = ref_alleles
    filtered_snv_df['alt_allele'] = alt_alleles
    final_df = filtered_snv_df[['Position', 'ref_allele', 'alt_allele']].rename(columns={'Position': 'position'})
    final_df.to_csv(output_path, index=False)
    logger.info(f"Saved {len(final_df)} unique SNVs to {output_path}")
    return final_df



def _filter_snv_by_fdr(snv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters SNVs by FDR value.

    Args:
        snv_df (pd.DataFrame): DataFrame with SNVs.

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    return snv_df[snv_df['FDR'] < 0.056].copy()



def _extract_alleles(filtered_snv_df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """
    Extracts reference and alternative alleles from filtered SNV DataFrame.

    Args:
        filtered_snv_df (pd.DataFrame): Filtered SNV DataFrame.

    Returns:
        Tuple[List[str], List[str]]: Lists of reference and alternative alleles.
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



def get_covered_positions(ref_constructs_dir: Path) -> Set[int]:
    """
    Returns all positions covered by reference constructs.

    Args:
        ref_constructs_dir (Path): Path to the directory with reference constructs.

    Returns:
        Set[int]: Set of covered positions.
    """
    if not ref_constructs_dir.exists():
        logger.error(f"Constructs directory does not exist: {ref_constructs_dir}")
        return set()
    if not ref_constructs_dir.is_dir():
        logger.error(f"Path is not a directory: {ref_constructs_dir}")
        return set()
    logger.info(f"Scanning constructs directory: {ref_constructs_dir}")
    ref_constructs = _load_ref_constructs(ref_constructs_dir)
    covered_positions = _collect_covered_positions(ref_constructs)
    return covered_positions



def _load_ref_constructs(ref_constructs_dir: Path) -> List[str]:
    """
    Loads all ConstructIDs from reference construct files in the directory.

    Args:
        ref_constructs_dir (Path): Path to the directory with reference constructs.

    Returns:
        List[str]: List of ConstructIDs.
    """
    ref_constructs: List[str] = []
    file_count = 0
    files = list(ref_constructs_dir.iterdir())
    logger.info(f"Found {len(files)} files in directory")
    for file in files:
        if file.name.endswith("-EF.csv"):
            file_count += 1
            try:
                df = pd.read_csv(file)
                construct_count = len(df)
                ref_constructs.extend(df['ConstructID'].tolist())
                logger.debug(f"File {file.name}: loaded {construct_count} constructs")
            except Exception as e:
                logger.error(f"Error reading file {file.name}: {str(e)}")
    logger.info(f"Checked {file_count} construct files")
    logger.info(f"Total loaded reference constructs: {len(ref_constructs)}")
    return ref_constructs



def _collect_covered_positions(ref_constructs: List[str]) -> Set[int]:
    """
    Collects all positions covered by the given reference constructs.

    Args:
        ref_constructs (List[str]): List of ConstructIDs.

    Returns:
        Set[int]: Set of covered positions.
    """
    covered_positions: Set[int] = set()
    processed_constructs = 0
    if not ref_constructs:
        logger.warning("No reference constructs found!")
        return covered_positions
    for construct_id in ref_constructs:
        arm_size, center, arm3_start, arm4_start = parse_construct_id(construct_id)
        if None in (arm_size, center, arm3_start, arm4_start):
            continue
        try:
            arm_ranges = calculate_arm_ranges(arm_size, center, arm3_start, arm4_start)
            for start, end in arm_ranges:
                covered_positions.update(range(start, end + 1))
            processed_constructs += 1
        except Exception as e:
            logger.error(f"Error processing construct {construct_id}: {str(e)}")
    logger.info(f"Processed {processed_constructs} constructs")
    if covered_positions:
        logger.info(f"Found {len(covered_positions)} covered positions")
        logger.info(f"Coverage range: {min(covered_positions)} to {max(covered_positions)}")
    else:
        logger.warning("No covered positions found!")
    return covered_positions



def apply_snvs(
    ref_record: SeqRecord,
    snv_df: pd.DataFrame,
    log_path: Path,
    covered_positions: Set[int],
    num: int
) -> SeqRecord:
    """
    Applies 2 random SNVs to the reference sequence, ensuring their presence in the reference constructs.

    Args:
        ref_record (SeqRecord): Reference sequence record.
        snv_df (pd.DataFrame): DataFrame with SNVs.
        log_path (Path): Path to the mutation log file.
        covered_positions (Set[int]): Set of positions covered by constructs.
        num (int): Sequence number for output naming.

    Returns:
        SeqRecord: Modified sequence record.
    """
    original_seq: str = str(ref_record.seq).upper()
    mutable_seq: MutableSeq = MutableSeq(original_seq)
    filtered_snvs = _filter_snvs_by_covered_positions(snv_df, covered_positions)
    logger.info(f"Found {len(filtered_snvs)} SNVs covered by constructs")
    unique_positions, selected_positions = _select_snv_positions(filtered_snvs)
    logger.info(f"Selected positions for mutation: {selected_positions}")
    applied_count, mismatch_log, selected_snvs = _apply_selected_snvs(
        mutable_seq, original_seq, filtered_snvs, selected_positions
    )
    _save_mismatch_log(mismatch_log, log_path)
    _log_snv_stats(num, filtered_snvs, unique_positions, selected_positions, applied_count)
    return SeqRecord(
        seq=mutable_seq,
        id=f"custom_mtDNA_{num}",
        description=(
            f"Modified from {ref_record.id} | Applied {applied_count} "
            f"of {len(selected_positions)} selected SNVs"
        )
    )



def _filter_snvs_by_covered_positions(
    snv_df: pd.DataFrame, covered_positions: Set[int]
) -> List[Dict[str, str]]:
    """
    Filters SNVs to only those covered by constructs.

    Args:
        snv_df (pd.DataFrame): DataFrame with SNVs.
        covered_positions (Set[int]): Set of covered positions.

    Returns:
        List[Dict[str, str]]: List of SNV dicts.
    """
    filtered_snvs: List[Dict[str, str]] = []
    for _, row in snv_df.iterrows():
        position = int(row['position'])
        if position in covered_positions:
            filtered_snvs.append({
                'position': position,
                'ref_allele': row['ref_allele'].upper(),
                'alt_allele': row['alt_allele'].upper()
            })
    return filtered_snvs



def _select_snv_positions(
    filtered_snvs: List[Dict[str, str]]
) -> Tuple[List[int], List[int]]:
    """
    Selects up to 2 unique SNV positions for mutation.

    Args:
        filtered_snvs (List[Dict[str, str]]): List of SNV dicts.

    Returns:
        Tuple[List[int], List[int]]: (all unique positions, selected positions)
    """
    unique_positions: List[int] = list({snv['position'] for snv in filtered_snvs})
    if len(unique_positions) >= 2:
        selected_positions = random.sample(unique_positions, 2)
    else:
        selected_positions = unique_positions
    return unique_positions, selected_positions



def _apply_selected_snvs(
    mutable_seq: MutableSeq,
    original_seq: str,
    filtered_snvs: List[Dict[str, str]],
    selected_positions: List[int]
) -> Tuple[int, List[Dict[str, Any]], List[Dict[str, str]]]:
    """
    Applies selected SNVs to the mutable sequence.

    Args:
        mutable_seq (MutableSeq): Mutable sequence to modify.
        original_seq (str): Original sequence string.
        filtered_snvs (List[Dict[str, str]]): List of SNV dicts.
        selected_positions (List[int]): Positions to mutate.

    Returns:
        Tuple[int, List[Dict[str, Any]], List[Dict[str, str]]]:
            (number applied, mismatch log, selected SNVs)
    """
    applied_count = 0
    mismatch_log: List[Dict[str, Any]] = []
    selected_snvs: List[Dict[str, str]] = []
    for position in selected_positions:
        position_snvs = [snv for snv in filtered_snvs if snv['position'] == position]
        idx = position - 1
        if idx >= len(original_seq):
            mismatch_log.append({
                'position': position,
                'original_base': None,
                'ref_allele': '|'.join([s['ref_allele'] for s in position_snvs]),
                'alt_allele': '|'.join([s['alt_allele'] for s in position_snvs]),
                'status': 'SKIPPED',
                'notes': 'Position out of sequence bounds'
            })
            continue
        current_base = original_seq[idx]
        applied = False
        for snv in position_snvs:
            ref_allele = snv['ref_allele']
            alt_allele = snv['alt_allele']
            if ref_allele == alt_allele:
                continue
            if current_base == ref_allele:
                mutable_seq[idx] = alt_allele
                applied_count += 1
                applied = True
                selected_snvs.append(snv)
                mismatch_log.append({
                    'position': position,
                    'original_base': current_base,
                    'ref_allele': ref_allele,
                    'alt_allele': alt_allele,
                    'status': 'APPLIED',
                    'notes': ''
                })
                break
        if not applied:
            notes = []
            for snv in position_snvs:
                if current_base == snv['alt_allele']:
                    notes.append('ALT allele already present')
                else:
                    notes.append(f"Expected ref: {snv['ref_allele']}, found: {current_base}")
            mismatch_log.append({
                'position': position,
                'original_base': current_base,
                'ref_allele': '|'.join([s['ref_allele'] for s in position_snvs]),
                'alt_allele': '|'.join([s['alt_allele'] for s in position_snvs]),
                'status': 'SKIPPED',
                'notes': '; '.join(notes)
            })
    return applied_count, mismatch_log, selected_snvs



def _save_mismatch_log(mismatch_log: List[Dict[str, Any]], log_path: Path) -> None:
    """
    Saves the mismatch log to a CSV file if there are any entries.

    Args:
        mismatch_log (List[Dict[str, Any]]): List of mismatch log entries.
        log_path (Path): Path to the log file.
    """
    if mismatch_log:
        mismatch_df = pd.DataFrame(mismatch_log)
        mismatch_df.to_csv(log_path, index=False)
        logger.info(f"Mutation log saved to {log_path}")



def _log_snv_stats(
    num: int,
    filtered_snvs: List[Dict[str, str]],
    unique_positions: List[int],
    selected_positions: List[int],
    applied_count: int
) -> None:
    """
    Logs statistics about SNV application for a sequence.

    Args:
        num (int): Sequence number.
        filtered_snvs (List[Dict[str, str]]): Filtered SNVs.
        unique_positions (List[int]): All unique positions.
        selected_positions (List[int]): Selected positions for mutation.
        applied_count (int): Number of successfully applied SNVs.
    """
    logger.info(f"\nSNV application statistics for sequence #{num}:")
    logger.info(f"Total SNVs covered by constructs: {len(filtered_snvs)}")
    logger.info(f"Unique positions: {len(unique_positions)}")
    logger.info(f"Selected positions: {len(selected_positions)}")
    logger.info(f"Successfully applied: {applied_count}")



def main(num: int) -> None:
    """
    Main function for generating a sequence with SNVs.

    Args:
        num (int): Sequence number for output naming.
    """
    base_dir = Path(__file__).resolve().parent.parent
    snv_csv_path = base_dir / "snv_csv" / "snvs.csv"
    input_fasta = base_dir / "sequences" / "ref_seq" / "Homo_sapiens_assembly38.chrM.fasta"
    xlsx_path = base_dir / "raw_data" / "MitoPhewas_associations.xlsx"
    log_path = base_dir / "snv_log" / f"snv_log_{num}.csv"
    output_fasta = base_dir / "sequences" / "relative_seq" / f"test_individual_{num+4}.fasta"
    ref_constructs_dir = base_dir.parent / "MitoFragilityScore" / "Energies" / "SEQ-g38_Mt-Short_Test"

    if not input_fasta.exists():
        logger.error(f"Reference FASTA not found: {input_fasta}")
        return
    if not xlsx_path.exists():
        logger.error(f"Excel file not found: {xlsx_path}")
        return

    snv_df = _get_snv_df(num, snv_csv_path, xlsx_path)
    ref_record = SeqIO.read(input_fasta, "fasta")
    logger.info(f"Loaded reference sequence: {ref_record.id}")
    logger.info(f"Length: {len(ref_record.seq)} bp")
    covered_positions = get_covered_positions(ref_constructs_dir)
    custom_record = apply_snvs(ref_record, snv_df, log_path, covered_positions, num)
    SeqIO.write(custom_record, output_fasta, "fasta")
    logger.info(f"Result saved to {output_fasta}")
    logger.info(f"ID: {custom_record.id}")
    logger.info(f"Description: {custom_record.description}")



def _get_snv_df(num: int, snv_csv_path: Path, xlsx_path: Path) -> pd.DataFrame:
    """
    Gets the SNV DataFrame, creating it from Excel if needed.

    Args:
        num (int): Sequence number.
        snv_csv_path (Path): Path to the SNV CSV file.
        xlsx_path (Path): Path to the Excel file.

    Returns:
        pd.DataFrame: SNV DataFrame.
    """
    if num == 0 and not snv_csv_path.exists():
        return csv_constructor(xlsx_path, snv_csv_path)
    return pd.read_csv(snv_csv_path)


if __name__ == "__main__":
    for i in range(5):
        logger.info(f"\n{'=' * 50}")
        logger.info(f"Generating sequence #{i}")
        logger.info(f"{'=' * 50}")
        main(i)
