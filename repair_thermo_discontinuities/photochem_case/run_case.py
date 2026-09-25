"""Regenerate and run the bundled TOI-1231 b diagnostic."""

import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
LOG = HERE / "photochem_case_run.jsonl"


def main():
    subprocess.run(
        [sys.executable, "-B", str(HERE.parent / "run_repair.py"),
         "--skip-independent", "--skip-case-inputs", "--skip-plots"],
        check=True,
    )
    subprocess.run([sys.executable, "-B", str(HERE / "generate_case_inputs.py")], check=True)
    with LOG.open("w", encoding="utf-8") as stream:
        subprocess.run(
            [sys.executable, "-B", "-u", str(HERE / "diagnose_initial.py"),
             "--steps", "1600", "--interval", "200"],
            cwd=HERE, stdout=stream, check=True,
        )
    records = [json.loads(line) for line in LOG.read_text(encoding="utf-8").splitlines()]
    steps = [record for record in records if record.get("label") == "stepped"]
    if not steps or not steps[-1]["converged"]:
        raise RuntimeError(f"Case did not converge; see {LOG}")
    print(f"Converged at step {steps[-1]['step']}; log: {LOG}")


if __name__ == "__main__":
    main()
