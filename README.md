# DataPreparing - SNV Generation and Analysis

Part of the MitoFragility pipeline. This module generates Single Nucleotide Variants (SNVs) from input FASTA files and performs statistical analysis.

## Quick Start

### 1. Setup Python Environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Prepare Input Data

Place FASTA files in `fastas/` directory:

```bash
cp your_sequences/*.fasta fastas/
```

### 3. Run Pipeline

```bash
python run_pipeline.py
```

Or via SLURM:

```bash
cd ../scripts
sbatch step1_dataprep_snv_generation.slurm
```

## Input Structure

```
fastas/
├── individual_1.fasta    # Input genome sequences
├── individual_2.fasta
└── individual_3.fasta

raw_data/                 # Raw SNP variant data
├── snps_list.csv
└── mutations.txt
```

## Output Structure

```
output/
├── *.csv              # Generated SNV files
├── *.json             # Analysis data
└── *.html             # Visualization reports

logs/
├── *.log              # Execution logs
└── *.txt              # Error logs
```

## Configuration

Edit `config/pipeline_config.conf`:

```ini
[PATHS]
OUTPUT_DIR = ./output
LOG_DIR = ./logs

[PARAMETERS]
MIN_SNV = 5
MAX_ENERGY = 10.0
```

## Module Descriptions

### dataprep/

- **data_loader.py** - Load and parse FASTA/JSON files
- **snv_filter.py** - Filter and validate SNVs
- **stats_utils.py** - Statistical analysis (Wilcoxon, etc.)
- **plots.py** - Visualization functions
- **html_builder.py** - Build interactive HTML reports
- **html_generator.py** - Generate summary HTML
- **runner.py** - Main execution orchestrator
- **utils.py** - Helper utilities

## Key Functions

### run_pipeline.py

Main entry point. Processes all individuals and generates outputs:

```python
python run_pipeline.py --help
```

### runner.py

Core processing logic:

```python
from dataprep.runner import main
main(base_dir="./fastas", output_base_dir="./output")
```

## Requirements

- Python 3.8+
- See `requirements.txt` for dependencies:
  - numpy
  - pandas
  - matplotlib
  - plotly
  - biopython
  - And others...

## Troubleshooting

### No FASTA files found

```bash
ls -la fastas/
# Ensure .fasta files are present
```

### Python module import errors

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

### Memory issues

Increase available memory or process smaller datasets.

## Integration

This module is the first step of the MitoFragility pipeline:

1. **Step 1** ← You are here: Generate SNVs
2. Step 2: Prepare input for Docker
3. Step 3: Run MitoFragilityScore Docker analysis
4. Step 4: Visualize and archive results

See `../README.md` for full pipeline documentation.

## Repository Info

- **GitHub**: https://github.com/olegpecheneg/study
- **Branch**: refactoring/dataprep-modularization
- **Remote**: origin

## Notes

- Generated output files are excluded from git (see .gitignore)
- Keep raw_data/ for reference but don't commit large data files
- Logs are automatically cleaned after 30 days
