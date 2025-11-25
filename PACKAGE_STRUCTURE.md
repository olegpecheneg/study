# DataPreparing Package Structure

## Overview
The `DataPreparing` project has been refactored into a modular package-based structure with clean separation of concerns. All business logic is moved into focused modules, leaving entry-point scripts clean and readable.

## Module Organization

### Core Package: `dataprep/`

#### Configuration & Paths
- **`config.py`**
  - `read_pipeline_config()` — Load config from `pipeline_config.conf`
  - `resolve_project_paths()` — Resolve absolute paths for project resources
  - Safe JOB_ID handling with validation

#### Data Input/Output
- **`excel_parser.py`**
  - `load_excel_file()` — Load data from Excel sheets
  - `get_sheet_names()` — List available sheets
  - `load_all_sheets()` — Load all sheets at once
  - `extract_column_values()` — Extract column data
  - `filter_dataframe_by_column_value()` — DataFrame filtering
  - `merge_dataframes()` — Combine multiple DataFrames
  - `validate_required_columns()` — Column validation
  - `save_dataframe_to_excel()` — Export to Excel

- **`io.py`**
  - `save_detailed_data()` — Save mutation data (CSV/JSON)
  - `save_search_csv()` — Save search results
  - `save_master_json()` — Consolidate results
  - `write_outlier_stats()` — Write statistical summaries

#### Data Processing & SNVs
- **`snputils.py`**
  - `load_snp_data()` — Load SNP data from CSV
  - `calculate_arm_ranges()` — Calculate genomic arm regions
  - `get_snps_in_construct()` — Filter SNPs by construct region
  - `_try_add_snp_position()` — Position lookup with fallback

- **`snv_filter.py`** (New: mt_DNA support)
  - `filter_snv_by_fdr()` — Filter by FDR threshold
  - `filter_by_positions()` — Filter by specific positions
  - `filter_by_region()` — Filter by genomic region
  - `filter_by_trait()` — Filter by phenotype trait
  - `extract_alleles()` — Extract ref/alt alleles
  - `get_covered_positions()` — Identify covered regions

- **`sequence_builder.py`** (New: mutation logic)
  - `apply_mutation()` — Apply single mutation to sequence
  - `apply_multiple_mutations()` — Batch mutation application
  - `translate_dna_to_protein()` — DNA to protein translation
  - `calculate_mutation_impact()` — Classify mutation type
  - `extract_region()` — Extract subsequence
  - `calculate_gc_content()` — Compute GC content
  - `reverse_complement()` — Get reverse complement

#### Analysis & Statistics
- **`stats.py`**
  - `calculate_outlier_stats()` — Compute mean, std, outliers
  - Handles multi-dimensional energy data

- **`parsers.py`**
  - `parse_construct_id()` — Parse construct naming (CGS/CEN/CON)
  - Position extraction from construct names
  - Fallback for non-standard names

#### Visualization
- **`colors.py`**
  - `generate_distinct_colors()` — HSV-based color generation
  - Ensures color distinctness for visualization

- **`plot_helpers.py`**
  - `plot_scatter_points()` — Scatter plot creation
  - `add_diagonal_line()` — Reference line addition
  - `add_outlier_zones()` — Mark statistical zones
  - `create_legend_elements()` — Legend generation

- **`plot_export.py`**
  - `plot_energy_comparison()` — Large multi-format export
  - Generates PNG, HTML (interactive), Plotly JSON
  - Handles energy scatter with filtering

- **`setup.py`** (New: initialization)
  - `setup_matplotlib()` — Configure matplotlib backend
  - `setup_logger()` — Configure logging (stream + file)
  - Separate error log for critical issues

#### Pipeline Orchestration
- **`runner.py`**
  - `process_individual()` — Main processing pipeline
  - `find_individual_dirs()` — Discover test individual directories
  - `find_points_by_sequence()` — Search functionality
  - `search_sequences_wrapper()` — Interactive search wrapper
  - `main()` — Top-level orchestration

#### mt_DNA Specific (New)
- **`mt_builder.py`** (New: mt_DNA sequence construction)
  - `get_covered_positions_from_settings()` — Load coverage from config
  - `csv_constructor()` — Build SNV CSV from Excel with filters
  - `apply_single_snv()` — Apply SNV to reference sequence
  - `process_snvs()` — Batch SNV processing
  - Integration with Bio.SeqIO for FASTA I/O

- **`mt_output.py`** (New: mt_DNA output)
  - `save_individual_results()` — Save mutations + energies + sequences
  - `save_mutations_to_csv()` — Export mutations as CSV
  - `save_energy_comparison()` — Save energy analysis JSON
  - `create_summary_report()` — Generate summary JSON
  - `merge_results_to_dataframe()` — Consolidate results
  - `append_to_master_result_file()` — Update master file

## Entry-Point Scripts

### `script/scatter_plus_n_std.py`
**Status**: Clean entry-point (20 lines)
```python
# Imports only:
from dataprep.setup import setup_logger, setup_matplotlib
from dataprep.runner import main, search_sequences_wrapper

# if __name__ == "__main__":
#   if "search" in argv → search_sequences_wrapper()
#   else → main()
```

**Functionality**:
- Main scatter + n_std pipeline (default mode)
- Interactive sequence search mode (`--search`)
- All logic delegated to `dataprep.runner.main()`

### `script/mt_DNA_builder.py`
**Status**: Clean entry-point (~150 lines)
```python
# Imports only:
from dataprep.mt_builder import csv_constructor, process_snvs
from dataprep.config import read_pipeline_config

# if __name__ == "__main__":
#   argparse for --mode (excel|csv) and filters
#   → main_from_excel() or main_from_csv()
```

**Functionality**:
- Build mt_DNA sequences with SNVs from Excel or CSV
- Support for FDR, region, position, trait filters
- Config-driven DEFAULT_FDR overrides
- Generates FASTA + mutation logs

## Data Flow

### Scatter + n_std Pipeline
```
find_individual_dirs()
  ↓
for each individual:
  process_individual()
    ├─ load SNP data (snputils.load_snp_data)
    ├─ parse construct IDs (parsers.parse_construct_id)
    ├─ compute energy stats (stats.calculate_outlier_stats)
    ├─ generate colors (colors.generate_distinct_colors)
    ├─ plot scatter (plot_export.plot_energy_comparison)
    └─ save results (io.save_*)
  ↓
search_sequences_wrapper() [optional]
  └─ find_points_by_sequence()
```

### mt_DNA Builder Pipeline
```
csv_constructor() [if Excel mode]
  ├─ load Excel (excel_parser.load_excel_file)
  ├─ filter SNVs (snv_filter.filter_*)
  ├─ extract alleles (snv_filter.extract_alleles)
  └─ save CSV
  ↓
process_snvs()
  ├─ load FASTA reference
  ├─ get covered positions (mt_builder.get_covered_positions_from_settings)
  └─ for each SNV:
      ├─ apply mutation (sequence_builder.apply_mutation)
      ├─ log result (mt_builder._save_mismatch_log)
      └─ write sequence (Bio.SeqIO)
```

## Configuration

### `pipeline_config.conf`
Located at: `d:\pythonProject\MitoFragility\config\pipeline_config.conf`

Example:
```ini
PROJECT_ROOT=d:\pythonProject\MitoFragility
CONSTRUCTS_DIR=%(PROJECT_ROOT)s\MitoFragilityScore\Constructs
ENERGIES_DIR=%(PROJECT_ROOT)s\MitoFragilityScore\Energies
OUTPUT_DIR=%(PROJECT_ROOT)s\DataPreparing\output
DEFAULT_FDR=0.056
JOB_ID_PREFIX=analysis
```

### Environment Setup
- **`setup_env.ps1`** — Windows PowerShell venv setup
- **`setup_env.sh`** — Unix/Linux bash venv setup
- **`requirements.txt`** — Dependencies (numpy, pandas, matplotlib, plotly, biopython, mplcursors, mpld3)

## Key Design Principles

1. **Single Responsibility**: Each module handles one logical domain
2. **Config-Driven**: Paths resolved from `pipeline_config.conf`, fallback to defaults
3. **Error Resilience**: Try/except blocks with detailed logging
4. **Clean Entry-Points**: Scripts are 90% imports, 10% orchestration
5. **Logging**: Dual-stream (console + file) with separate error log
6. **Testing Paths**: All paths use `Path` object for cross-platform compatibility

## Usage Examples

### Run Main Pipeline
```bash
cd d:\pythonProject\MitoFragility\DataPreparing
python script/scatter_plus_n_std.py
```

### Search Sequences
```bash
python script/scatter_plus_n_std.py search
```

### Build mt_DNA from Excel
```bash
python script/mt_DNA_builder.py --mode excel --run-id 1 --fdr 0.05
```

### Build mt_DNA from CSV
```bash
python script/mt_DNA_builder.py --mode csv --run-id 1 --csv-path path/to/snvs.csv
```

## Module Dependencies Graph

```
dataprep/
├─ config.py
│  └─ [no internal deps]
├─ excel_parser.py
│  └─ config.py (optional)
├─ io.py
│  └─ [no internal deps]
├─ snputils.py
│  └─ [no internal deps]
├─ snv_filter.py
│  └─ [no internal deps]
├─ sequence_builder.py
│  └─ [no internal deps]
├─ stats.py
│  └─ [no internal deps]
├─ parsers.py
│  └─ [no internal deps]
├─ colors.py
│  └─ [no internal deps]
├─ plot_helpers.py
│  └─ [no internal deps]
├─ setup.py
│  └─ [no internal deps]
├─ plot_export.py
│  └─ plot_helpers, colors
├─ runner.py
│  └─ config, snputils, parsers, stats, colors, plot_export, io
├─ mt_builder.py
│  └─ snv_filter, excel_parser, sequence_builder (indirect)
└─ mt_output.py
   └─ [no internal deps]

script/
├─ scatter_plus_n_std.py
│  └─ setup, runner
└─ mt_DNA_builder.py
   └─ config, mt_builder, excel_parser
```

## Testing & Validation

All modules have been syntax-validated with `py_compile`. To verify:

```bash
python -m py_compile dataprep/*.py
python -m py_compile script/*.py
```

## Future Enhancements

- [ ] Add type hints validation (mypy)
- [ ] Create unit tests for each module
- [ ] Add doctest examples to modules
- [ ] Implement CLI with Click instead of argparse
- [ ] Add progress bars for batch operations
- [ ] Implement caching for repeated calculations
