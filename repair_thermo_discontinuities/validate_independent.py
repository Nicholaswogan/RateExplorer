"""Compare six unambiguously matched gas species with NASA9 source fits.

The reference files are fetched at validation time; their SHA256 hashes are
checked so changed upstream files cannot silently change this comparison.
Burcat's source license does not permit bundling database excerpts here.
"""

import csv
import hashlib
import math
import re
import urllib.request
from pathlib import Path

import numpy as np
import yaml

from match_shomate_joins import enthalpy_entropy, heat_capacity
from source_snapshot import REPAIRED_GAS, RESULTS, original_text


HERE = Path(__file__).resolve().parent
GAS_SOURCE = "photochem_clima_data/data/reaction_mechanisms/zahnle_earth.yaml"
R = 8.31446261815324
SOURCES = {
    "NASA CEA": (
        "https://raw.githubusercontent.com/nasa/cea/3f4441d28a02fccbb140e1a028d9902390981389/data/thermo.inp",
        "fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363",
    ),
    "Burcat": (
        "https://respecth.elte.hu/burcat/NEWNASA.TXT",
        "3e6717ed4e0e266eeb961c116a7589d8e8d71f63c9a935eebea7e263eb25c440",
    ),
}
SPECIES = {
    "C3H6": ("NASA CEA", "C3H6,propylene", 1400.0),
    "CH2CO": ("NASA CEA", "CH2CO,ketene", 1500.0),
    "CH3CN": ("NASA CEA", "CH3CN             Acetonitrile", 1500.0),
    "CH3CO": ("NASA CEA", "CH3CO,acetyl", 1500.0),
    "CH3O2": ("Burcat", "CH3O2 Methyl Peroxy Rad  ", 1500.0),
    "CH2CN": ("Burcat", "CH2CN Methyl-Cyanid Radical", 1500.0),
}
NUMBER = re.compile(r"[+-]?\d+\.\d+[DE][+-]\d+")


def fetch(label):
    url, expected_sha = SOURCES[label]
    with urllib.request.urlopen(url, timeout=40) as response:
        data = response.read()
    actual = hashlib.sha256(data).hexdigest()
    if actual != expected_sha:
        raise ValueError(f"{label} source changed: {actual} != {expected_sha}")
    return data.decode("utf-8").splitlines()


def parse_nasa9(lines, marker):
    matches = [i for i, line in enumerate(lines) if line.startswith(marker)]
    if len(matches) != 1:
        raise ValueError(f"Expected one source row for {marker!r}, found {len(matches)}")
    start = matches[0]
    count = int(lines[start + 1][:2])
    fits = []
    for i in range(count):
        row = start + 2 + 3 * i
        lower, upper = float(lines[row][:11]), float(lines[row][11:22])
        values = [float(x.replace("D", "E")) for x in NUMBER.findall(lines[row + 1] + lines[row + 2])]
        if len(values) not in (9, 10):
            raise ValueError(f"Malformed NASA9 row for {marker}")
        fits.append((lower, upper, values[:7] + values[-2:]))
    return fits


def nasa_properties(row, temperature):
    a, b, c, d, e, f, g, h, i = row
    t = temperature
    enthalpy = R * t * (-a / t**2 + b * math.log(t) / t + c + d * t / 2
                         + e * t**2 / 3 + f * t**3 / 4 + g * t**4 / 5 + h / t)
    entropy = R * (-a / (2 * t**2) - b / t + c * math.log(t) + d * t
                   + e * t**2 / 2 + f * t**3 / 3 + g * t**4 / 4 + i)
    cp = R * (a / t**2 + b / t + c + d * t + e * t**2 + f * t**3 + g * t**4)
    return enthalpy, entropy, enthalpy - t * entropy, cp


def shomate_properties(row, temperature):
    h, s = enthalpy_entropy(row, temperature)
    return h, s, h - temperature * s, heat_capacity(row, temperature)


def main():
    sources = {label: fetch(label) for label in SOURCES}
    old = {x["name"]: x for x in yaml.safe_load(original_text(GAS_SOURCE))["species"]}
    new = {x["name"]: x for x in yaml.safe_load(REPAIRED_GAS.read_text())["species"]}
    rows = []
    for name, (label, marker, join) in SPECIES.items():
        fits = parse_nasa9(sources[label], marker)
        temperatures = np.linspace(join + 1, 6000, 500)
        differences = {version: [] for version in ("old", "repaired")}
        for t in temperatures:
            reference = next((values for lo, hi, values in fits if lo <= t <= hi), None)
            if reference is None:
                raise ValueError(f"No NASA9 source range for {name} at {t}")
            reference_values = nasa_properties(reference, float(t))
            for version, entry in (("old", old[name]), ("repaired", new[name])):
                row = entry["thermo"]["data"][2]
                differences[version].append(np.subtract(shomate_properties(row, float(t)), reference_values))
        row = {"species": name, "source": label, "min_K": join + 1, "max_K": 6000}
        for version in differences:
            values = np.array(differences[version])
            for index, (property_name, unit) in enumerate((("H_kJ_mol", 1000), ("S_J_mol_K", 1), ("G_kJ_mol", 1000), ("Cp_J_mol_K", 1))):
                residual = values[:, index] / unit
                row[f"rms_{property_name}_{version}"] = float(np.sqrt(np.mean(residual**2)))
                row[f"max_abs_{property_name}_{version}"] = float(np.max(np.abs(residual)))
        rows.append(row)
    with (RESULTS / "independent_nasa9_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row["species"], "G RMS old/repaired kJ/mol:",
              f"{row['rms_G_kJ_mol_old']:.2f}/{row['rms_G_kJ_mol_repaired']:.2f}")


if __name__ == "__main__":
    main()
