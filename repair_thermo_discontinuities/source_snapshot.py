"""Retrieve and verify the immutable photochem_clima_data v0.3.2 source."""

import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PACKAGE = RESULTS / "photochem_clima_data_v0.3.2"
REPAIRED_GAS = HERE / "zahnle_earth.yaml"
REPAIRED_CONDENSATE = HERE / "condensate_thermo.yaml"
REMOTE = "https://github.com/Nicholaswogan/photochem_clima_data.git"
TAG = "v0.3.2"
COMMIT = "66789aba287a4eb8a846e3913ce346c0ceb5215f"


def ensure_source():
    RESULTS.mkdir(exist_ok=True)
    if not PACKAGE.exists():
        subprocess.run(
            ["git", "clone", "--depth", "1", "--branch", TAG, "--single-branch", REMOTE, str(PACKAGE)],
            check=True,
        )
    actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PACKAGE, text=True).strip()
    if actual != COMMIT:
        raise ValueError(f"Expected {TAG} at {COMMIT}, found {actual}")
    return PACKAGE


def original_text(relative_path):
    ensure_source()
    return subprocess.check_output(["git", "show", f"HEAD:{relative_path}"], cwd=PACKAGE, text=True)
