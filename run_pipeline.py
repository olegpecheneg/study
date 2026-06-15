#!/usr/bin/env python3
"""
Единый пайплайн:
- Генерация FASTA (scripts/mt_DNA_builder.py)
- Анализ JSON (dataprep/main.py)

Поддерживает три режима:
  generate   – только генерация FASTA
  analyze    – только анализ JSON
  full       – генерация + анализ
"""

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
MT_BUILDER = PROJECT_ROOT / "script" / "mt_DNA_builder.py"
MAIN_SCRIPT = PROJECT_ROOT / "dataprep" / "main.py"

DEFAULT_EXCLUDE_SNPS = [
    "11719A>G", "10238T>C", "11778A>G", "1555A>G", "13617T>C",
    "11899T>C", "10081T>C", "7559A>G", "13708G>A", "8655C>T"
]

def run_generate(args):
    if not MT_BUILDER.exists():
        print(f"❌ Ошибка: не найден {MT_BUILDER}")
        return 1

    cmd = [
        sys.executable, str(MT_BUILDER),
        "--mode", "excel",
        "--run-id", str(args.run_id),
        "--fdr", str(args.fdr),
    ]
    if args.exclude_snps:
        cmd.extend(["--exclude-snps", args.exclude_snps])
    if args.positions:
        cmd.extend(["--positions", args.positions])
    if args.region:
        cmd.extend(["--region", args.region])
    if args.trait:
        cmd.extend(["--trait", args.trait])
    if getattr(args, 'synonymous_only', False):
        cmd.append("--synonymous-only")

    print("Запуск генерации FASTA...")
    print(">", " ".join(cmd))
    result = subprocess.run(cmd, cwd=PROJECT_ROOT)
    return result.returncode

def run_analyze(args):
    data_dir = Path(args.data_dir).resolve()
    ref_data_dir = Path(args.ref_data_dir).resolve() if args.ref_data_dir else None
    analysis_output_dir = Path(args.analysis_output_dir).resolve() if args.analysis_output_dir else None

    cmd = [
        sys.executable,
        "-m",
        "dataprep.main",
        str(data_dir)
    ]
    if ref_data_dir:
        cmd.extend(["-r", str(ref_data_dir)])
    if analysis_output_dir:
        cmd.extend(["-o", str(analysis_output_dir)])
    if getattr(args, 'summary_only', False):
        cmd.append("--summary-only")

    print("📊 Запуск анализа готовых мастер-данных...")
    print(">", " ".join(cmd))
    result = subprocess.run(cmd, cwd=PROJECT_ROOT)
    return result.returncode

def main():
    parser = argparse.ArgumentParser(description="Полный пайплайн: генерация FASTA и анализ JSON")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # generate
    gen_parser = subparsers.add_parser("generate", help="Только генерация FASTA")
    gen_parser.add_argument("--run-id", type=int, default=1)
    gen_parser.add_argument("--fdr", type=float, default=0.056)
    gen_parser.add_argument("--exclude-snps", type=str, default=",".join(DEFAULT_EXCLUDE_SNPS))
    gen_parser.add_argument("--positions", type=str)
    gen_parser.add_argument("--region", type=str)
    gen_parser.add_argument("--trait", type=str)
    gen_parser.add_argument("--synonymous-only", action="store_true")

    # analyze
    ana_parser = subparsers.add_parser("analyze", help="Только анализ JSON")
    ana_parser.add_argument("data_dir", type=str)
    ana_parser.add_argument("-r", "--ref_data_dir", type=str, default=None)
    ana_parser.add_argument("-o", "--analysis-output-dir", type=str, default=None)
    ana_parser.add_argument("--summary-only", action="store_true", help="Только создать агрегированный сводный отчет без анализа.")

    # full
    full_parser = subparsers.add_parser("full", help="Полный пайплайн: генерация + анализ")
    full_parser.add_argument("--run-id", type=int, default=1)
    full_parser.add_argument("--fdr", type=float, default=0.056)
    full_parser.add_argument("--exclude-snps", type=str, default=",".join(DEFAULT_EXCLUDE_SNPS))
    full_parser.add_argument("--positions", type=str)
    full_parser.add_argument("--region", type=str)
    full_parser.add_argument("--trait", type=str)
    full_parser.add_argument("--synonymous-only", action="store_true")
    full_parser.add_argument("data_dir", type=str)
    full_parser.add_argument("-r", "--ref_data_dir", type=str, default=None)
    full_parser.add_argument("-o", "--analysis-output-dir", type=str, default=None)

    args = parser.parse_args()

    if args.command == "generate":
        return run_generate(args)
    elif args.command == "analyze":
        return run_analyze(args)
    elif args.command == "full":
        ret = run_generate(args)
        if ret != 0:
            print("❌ Генерация FASTA завершилась с ошибкой, анализ не запущен.")
            return ret
        ret = run_analyze(args)
        return ret
    else:
        parser.print_help()
        return 1

if __name__ == "__main__":
    sys.exit(main())