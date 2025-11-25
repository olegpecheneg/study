"""Entry-point script for scatter + n_std visualization pipeline.

This script orchestrates the full pipeline:
1. Find individual test directories
2. Process each individual (mutations, energy, visualization)
3. Generate scatter plots with energy comparisons
4. Search for points by sequence (optional)
"""

import sys
import logging
from pathlib import Path

from dataprep.setup import setup_logger, setup_matplotlib
from dataprep.runner import main, search_sequences_wrapper

# Initialize matplotlib and logger
setup_matplotlib()
logger = setup_logger("scatter_plus_n_std")


def run_main_pipeline() -> None:
    """Run the main scatter + n_std pipeline."""
    try:
        main()
    except Exception as e:
        logger.error(f"Error in main pipeline: {e}", exc_info=True)
        sys.exit(1)


def run_search_pipeline() -> None:
    """Run the search sequences pipeline."""
    try:
        search_sequences_wrapper()
    except Exception as e:
        logger.error(f"Error in search pipeline: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "search":
        run_search_pipeline()
    else:
        run_main_pipeline()
