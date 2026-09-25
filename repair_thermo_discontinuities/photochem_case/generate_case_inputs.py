"""Generate the bundled TOI-1231 b Photochem case inputs from repaired YAMLs.

Run with the photochem environment from this directory. The Photochem
installation is not edited; its data root is redirected only in this process.
"""

import shutil
import sys
import tempfile
from pathlib import Path

from photochem.utils import _format

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from source_snapshot import REPAIRED_CONDENSATE, REPAIRED_GAS


HERE = Path(__file__).resolve().parent
CASE = HERE

def main():
    with tempfile.TemporaryDirectory(prefix="photochem-repaired-data-", dir=HERE) as directory:
        source = Path(directory) / "reaction_mechanisms"
        source.mkdir()
        shutil.copy2(REPAIRED_GAS, source / "zahnle_earth.yaml")
        shutil.copy2(REPAIRED_CONDENSATE, source / "condensate_thermo.yaml")
        previous_data_dir = _format.DATA_DIR
        try:
            _format.DATA_DIR = directory
            _format.zahnle_rx_and_thermo_files(
                atoms_names=["H", "He", "N", "O", "C", "S"],
                exclude_species=["S3", "S4", "S8", "S8aer"],
                rxns_filename=str(CASE / "photochem_rxns.yaml"),
                thermo_filename=str(CASE / "photochem_thermo.yaml"),
                remove_reaction_particles=True,
            )
        finally:
            _format.DATA_DIR = previous_data_dir
    print("Generated TOI-1231 b inputs in", CASE)


if __name__ == "__main__":
    main()
