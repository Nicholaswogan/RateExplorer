"""Regenerate repaired YAMLs and numerical validation artifacts from v0.3.2.

Run: python run_repair.py
Use --skip-independent --skip-plots to avoid fetching NASA/Burcat sources.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
from source_snapshot import REPAIRED_CONDENSATE, REPAIRED_GAS, RESULTS, ensure_source, original_text

HERE = Path(__file__).resolve().parent
GAS_SOURCE = "photochem_clima_data/data/reaction_mechanisms/zahnle_earth.yaml"
CONDENSATE_SOURCE = "photochem_clima_data/data/reaction_mechanisms/condensate_thermo.yaml"
GAS_OUTPUT = REPAIRED_GAS
CONDENSATE_OUTPUT = REPAIRED_CONDENSATE


def run(*arguments):
    subprocess.run([sys.executable, "-B", *map(str, arguments)], cwd=HERE, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-independent", action="store_true",
                        help="skip NASA9 metrics (plots still fetch NASA9 unless --skip-plots)")
    parser.add_argument("--skip-case-inputs", action="store_true",
                        help="skip Photochem case-input generation")
    parser.add_argument("--skip-plots", action="store_true",
                        help="skip all gas-fit figures")
    args = parser.parse_args()
    ensure_source()
    with tempfile.TemporaryDirectory(prefix="zahnle-repair-baseline-", dir=RESULTS) as directory:
        gas_baseline = Path(directory) / "zahnle_earth.yaml"
        condensate_baseline = Path(directory) / "condensate_thermo.yaml"
        gas_baseline.write_text(original_text(GAS_SOURCE), encoding="utf-8")
        condensate_baseline.write_text(original_text(CONDENSATE_SOURCE), encoding="utf-8")
        run("match_shomate_joins.py", gas_baseline, GAS_OUTPUT)
        run("derive_condensates.py", GAS_OUTPUT, condensate_baseline, CONDENSATE_OUTPUT)
    print(f"Repaired reaction mechanism: {GAS_OUTPUT}", flush=True)
    print(f"Repaired condensate thermodynamics: {CONDENSATE_OUTPUT}", flush=True)
    run("validate.py")
    if not args.skip_plots:
        run("plot_corrections.py")
        run("plot_join_repair_histogram.py")
    if not args.skip_independent:
        run("validate_independent.py")
    if not args.skip_case_inputs:
        run("photochem_case/generate_case_inputs.py")


if __name__ == "__main__":
    main()
