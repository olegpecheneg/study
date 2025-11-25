"""Entry-point script for building mt_DNA sequences with SNVs."""

import logging
import argparse
from pathlib import Path

from dataprep.config import read_pipeline_config
from dataprep.mt_builder import (
    csv_constructor,
    get_covered_positions_from_settings,
    process_snvs
)
from dataprep.excel_parser import load_excel_file
import pandas as pd

# Setup logger
logger = logging.getLogger("mt_DNA_builder")
if not logger.hasHandlers():
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def main_from_excel(
    num: int,
    position_filter: str = None,
    region_filter: str = None,
    fdr_threshold: float = 0.056,
    trait_filter: str = None
) -> None:
    """Generate sequences with SNVs from Excel.

    Args:
        num: Sequence number for output naming.
        position_filter: Comma-separated positions to filter.
        region_filter: Region to filter (e.g., "chrM:1000-2000").
        fdr_threshold: FDR threshold for filtering.
        trait_filter: Trait name to filter.
    """
    base_dir = Path(__file__).resolve().parent.parent
    snv_csv_path = base_dir / "snv_csv" / "snvs.csv"
    input_fasta = base_dir / "sequences" / "ref_seq" / "Homo_sapiens_assembly38.chrM.fasta"
    xlsx_path = base_dir / "raw_data" / "MitoPhewas_associations.xlsx"
    output_dir = base_dir / "sequences" / "relative_seq"
    
    settings_path = base_dir.parent / "MitoFragilityScore" / "src" / "settings.json"
    ref_seq_info_path = base_dir.parent / "MitoFragilityScore" / "Sequences" / "reference_sequence_information.json"

    if not input_fasta.exists():
        logger.error(f"Reference FASTA not found: {input_fasta}")
        return
    if not xlsx_path.exists():
        logger.error(f"Excel file not found: {xlsx_path}")
        return

    # Read DEFAULT_FDR from config if available
    cfg = read_pipeline_config()
    try:
        if 'DEFAULT_FDR' in cfg:
            fdr_threshold = float(cfg['DEFAULT_FDR'])
    except Exception:
        pass

    snv_df = csv_constructor(
        xlsx_path, snv_csv_path,
        position_filter=position_filter,
        region_filter=region_filter,
        fdr_threshold=fdr_threshold,
        trait_filter=trait_filter
    )
    
    process_snvs(snv_df, input_fasta, output_dir, settings_path, ref_seq_info_path, num)


def main_from_csv(num: int, snps_csv_path: Path) -> None:
    """Generate sequences with SNVs from existing CSV.

    Args:
        num: Sequence number for output naming.
        snps_csv_path: Path to the SNPs CSV file.
    """
    base_dir = Path(__file__).resolve().parent.parent
    input_fasta = base_dir / "sequences" / "ref_seq" / "Homo_sapiens_assembly38.chrM.fasta"
    output_dir = base_dir / "sequences" / "relative_seq"
    
    settings_path = base_dir.parent / "MitoFragilityScore" / "src" / "settings.json"
    ref_seq_info_path = base_dir.parent / "MitoFragilityScore" / "Sequences" / "reference_sequence_information.json"

    if not input_fasta.exists():
        logger.error(f"Reference FASTA not found: {input_fasta}")
        return
    if not snps_csv_path.exists():
        logger.error(f"SNPs CSV file not found: {snps_csv_path}")
        return

    snv_df = pd.read_csv(snps_csv_path)
    logger.info(f"Loaded {len(snv_df)} SNVs from {snps_csv_path}")
    
    required_columns = ['position', 'ref_allele', 'alt_allele']
    if not all(col in snv_df.columns for col in required_columns):
        logger.error(f"SNPs CSV must contain columns: {required_columns}")
        return

    process_snvs(snv_df, input_fasta, output_dir, settings_path, ref_seq_info_path, num)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mtDNA sequences with SNVs")
    parser.add_argument("--mode", choices=["excel", "csv"], required=True,
                       help="Input mode: excel (parse Excel file) or csv (use existing CSV)")
    parser.add_argument("--run-id", type=int, default=0,
                       help="Run ID for output naming")
    parser.add_argument("--csv-path", type=Path,
                       help="Path to SNPs CSV file (required for csv mode)")
    
    # Excel-specific filters
    parser.add_argument("--positions", type=str,
                       help="Comma-separated positions to filter (Excel mode only)")
    parser.add_argument("--region", type=str,
                       help="Genomic region to filter, e.g., 'chrM:1000-2000' (Excel mode only)")
    parser.add_argument("--fdr", type=float, default=0.056,
                       help="FDR threshold (Excel mode only)")
    parser.add_argument("--trait", type=str,
                       help="Trait name to filter (Excel mode only)")
    
    args = parser.parse_args()
    
    if args.mode == "excel":
        main_from_excel(
            num=args.run_id,
            position_filter=args.positions,
            region_filter=args.region,
            fdr_threshold=args.fdr,
            trait_filter=args.trait
        )
    elif args.mode == "csv":
        if not args.csv_path:
            parser.error("--csv-path is required for csv mode")
        main_from_csv(args.run_id, args.csv_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mtDNA sequences with SNVs")
    parser.add_argument("--mode", choices=["excel", "csv"], required=True,
                       help="Input mode: excel (parse Excel file) or csv (use existing CSV)")
    parser.add_argument("--run-id", type=int, default=0,
                       help="Run ID for output naming")
    parser.add_argument("--csv-path", type=Path,
                       help="Path to SNPs CSV file (required for csv mode)")
    
    # Excel-specific filters
    parser.add_argument("--positions", type=str,
                       help="Comma-separated positions to filter (Excel mode only)")
    parser.add_argument("--region", type=str,
                       help="Genomic region to filter, e.g., 'chrM:1000-2000' (Excel mode only)")
    parser.add_argument("--fdr", type=float, default=0.056,
                       help="FDR threshold (Excel mode only)")
    parser.add_argument("--trait", type=str,
                       help="Trait name to filter (Excel mode only)")
    
    args = parser.parse_args()
    
    if args.mode == "excel":
        main_from_excel(
            num=args.run_id,
            position_filter=args.positions,
            region_filter=args.region,
            fdr_threshold=args.fdr,
            trait_filter=args.trait
        )
    elif args.mode == "csv":
        if not args.csv_path:
            parser.error("--csv-path is required for csv mode")
        main_from_csv(args.run_id, args.csv_path)