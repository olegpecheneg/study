#!/usr/bin/env python3
"""
Оркестратор для запуска всех финальных анализов.
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("Orchestrator")

def run_script(script_path: Path, args_list: list, step_name: str) -> bool:
    cmd = [sys.executable, str(script_path)] + args_list
    logger.info(f"Запуск {step_name}: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"{step_name} завершён успешно")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} упал с кодом {e.returncode}")
        logger.error(f"STDERR:\n{e.stderr}")
        return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("-o", "--output_dir", type=Path, default=None)
    parser.add_argument("--fdr-alpha", type=float, default=0.01, help="Уровень значимости FDR")
    parser.add_argument("--skip-manhattan", action="store_true")
    parser.add_argument("--skip-volcano", action="store_true")
    parser.add_argument("--skip-extreme", action="store_true")
    parser.add_argument("--test", choices=["wilcoxon","ttest"], default="wilcoxon")
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument("--no-recurse", action="store_true")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        logger.error(f"Директория не найдена: {data_dir}")
        sys.exit(1)

    out_base = (args.output_dir or data_dir / "final_analysis").resolve()
    out_base.mkdir(parents=True, exist_ok=True)

    script_dir = Path(__file__).parent
    manhattan_script = script_dir / "manhattan_extreme.py"
    volcano_script = script_dir / "volcano_and_strip.py"
    extreme_script = script_dir / "extreme_snvs_and_heatmap.py"

    common_args = []
    if args.no_png:
        common_args.append("--no-png")
    if args.no_recurse:
        common_args.append("--no-recurse")

    # 1. Manhattan
    if not args.skip_manhattan and manhattan_script.exists():
        for etype in ["Energy", "EnergyLeft", "EnergyRight"]:
            e_out = out_base / "manhattan_extreme" / etype
            e_out.mkdir(parents=True, exist_ok=True)
            cmd_args = [str(data_dir), "-o", str(e_out), "--energy-type", etype,
                        "--fdr-alpha", str(args.fdr_alpha), *common_args]
            run_script(manhattan_script, cmd_args, f"Manhattan {etype}")

    # 2. Volcano + Strip
    if not args.skip_volcano and volcano_script.exists():
        v_out = out_base / "volcano_strip"
        v_out.mkdir(parents=True, exist_ok=True)
        cmd_args = [str(data_dir), "-o", str(v_out), "--test", args.test,
                    "--fdr-alpha", str(args.fdr_alpha), *common_args]
        run_script(volcano_script, cmd_args, "Volcano+Strip (all types)")

    # 3. Extreme + Heatmap
    if not args.skip_extreme and extreme_script.exists():
        for etype in ["Energy", "EnergyLeft", "EnergyRight"]:
            e_out = out_base / "extreme_analysis" / etype
            e_out.mkdir(parents=True, exist_ok=True)
            cmd_args = [str(data_dir), "-o", str(e_out), "--energy-type", etype,
                        "--top", "5", *common_args]
            run_script(extreme_script, cmd_args, f"Extreme+Heatmap {etype}")

    logger.info("Оркестратор завершил работу")

if __name__ == "__main__":
    main()