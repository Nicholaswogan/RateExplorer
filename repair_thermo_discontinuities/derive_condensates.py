"""Derive condensate Shomate fits from repaired gas and particle saturation.

Usage: python derive_condensates.py GAS.yaml TEMPLATE.yaml OUTPUT.yaml

The template fixes species order, compositions, phase boundaries and formatting.
Saturation parameters come only from GAS.yaml, never from separate fit files.
The historical fitting windows below come from saturation_thermo/*.py.
"""

import math
import sys
from pathlib import Path

import numpy as np
import yaml

from match_shomate_joins import enthalpy_entropy, replace_coefficients


R_J = 8.31446261815324
R_CGS = R_J * 1e7

# Original windows used by saturation_thermo/fitting.py callers.
WINDOWS = {
    "H2O": (100.0, 1000.0),
    "CO2": (150.0, 1000.0),
    "S8": (None, 1000.0),
    "H2SO4": (None, 1500.0),
    "NH3": (160.0, 1000.0),
    "N2O": (None, 1000.0),
    "C2H2": (120.0, 1000.0),
    "C2H4": (90.0, 1000.0),
    "C2H6": (None, 1000.0),
    "CH3CN": (None, 1000.0),
    "HCCCN": (None, 1000.0),
    "HCN": (140.0, 1000.0),
    "CH4": (None, 1000.0),
}


def gas_row(thermo, temperature):
    for index, upper in enumerate(thermo["temperature-ranges"][1:]):
        if temperature <= upper:
            return thermo["data"][index]
    return thermo["data"][-1]


def gibbs(row, temperature):
    enthalpy, entropy = enthalpy_entropy(row, temperature)
    return enthalpy - temperature * entropy


def sat_log_pressure(saturation, temperature):
    """Natural log of saturation pressure in bar, using LinearLatentHeat."""
    parameters = saturation["parameters"]
    triple = parameters["T-triple"]
    critical = parameters["T-critical"]
    reference = parameters["T-ref"]

    def integral(branch, temp):
        a, b = saturation[branch]["a"], saturation[branch]["b"]
        return -a / temp + b * math.log(temp)

    base = integral("vaporization", reference)
    if temperature <= triple:
        delta = (integral("vaporization", triple) - base
                 + integral("sublimation", temperature)
                 - integral("sublimation", triple))
    elif temperature <= critical:
        delta = integral("vaporization", temperature) - base
    else:
        delta = (integral("vaporization", critical) - base
                 + integral("super-critical", temperature)
                 - integral("super-critical", critical))
    return math.log(parameters["P-ref"] / 1e6) + parameters["mu"] / R_CGS * delta


def fit_branch(gas_thermo, saturation, lower, upper):
    temperatures = np.linspace(lower, upper, 100)
    # For [A, 0, 0, 0, 0, F, G], the condensed Gibbs energy is
    # A*T*(1-ln(T/1000)) + 1000*F - T*G.
    design = np.column_stack((
        temperatures * (1 - np.log(temperatures / 1000)),
        np.full(temperatures.shape, 1000.0),
        -temperatures,
    ))
    target = np.array([
        gibbs(gas_row(gas_thermo, float(t)), float(t))
        + R_J * t * sat_log_pressure(saturation, float(t))
        for t in temperatures
    ])
    # Fitting log pressure gives each temperature the same weight, as in the
    # original saturation_thermo/fitting.py objective.
    weights = 1 / (R_J * temperatures)
    a, f, g = np.linalg.lstsq(design * weights[:, None], target * weights, rcond=None)[0]
    return [float(a), 0.0, 0.0, 0.0, 0.0, float(f), float(g)]


def derive(gas, template):
    gases = {entry["name"]: entry for entry in gas["species"]}
    particles = {
        entry["name"]: entry for entry in gas["particles"]
        if entry.get("formation") == "saturation"
    }
    entries = template["species"]
    if set(particles) != {entry["name"] for entry in entries}:
        raise ValueError("Condensate template does not match saturation particles")
    replacements = {}
    for entry in entries:
        particle_name = entry["name"]
        gas_name = particle_name.removesuffix("aer")
        if gas_name not in WINDOWS or gas_name not in gases:
            raise ValueError(f"Missing gas/window for {particle_name}")
        saturation = particles[particle_name]["saturation"]
        if saturation["model"] != "LinearLatentHeat":
            raise ValueError(f"Unexpected saturation model for {particle_name}")
        triple = saturation["parameters"]["T-triple"]
        critical = saturation["parameters"]["T-critical"]
        if entry["thermo"]["temperature-ranges"] != [0.0, triple, critical, 6000]:
            raise ValueError(f"Unexpected phase boundaries for {particle_name}")
        low, high = WINDOWS[gas_name]
        if low is None:
            low = triple - 30
        gas_thermo = gases[gas_name]["thermo"]
        rows = [
            fit_branch(gas_thermo, saturation, low, triple),
            fit_branch(gas_thermo, saturation, triple, critical),
            fit_branch(gas_thermo, saturation, critical, high),
        ]
        # Preserve Gibbs continuity at each phase boundary while leaving
        # phase-transition enthalpy, entropy and heat capacity free to jump.
        for index, boundary in enumerate((triple, critical), start=1):
            rows[index][5] += (gibbs(rows[index - 1], boundary)
                               - gibbs(rows[index], boundary)) / 1000
        for index, row in enumerate(rows):
            replacements[(particle_name, index)] = row
    return replacements


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: python derive_condensates.py GAS.yaml TEMPLATE.yaml OUTPUT.yaml")
    gas_path, template_path, output_path = map(Path, sys.argv[1:])
    gas = yaml.safe_load(gas_path.read_text(encoding="utf-8"))
    template_text = template_path.read_text(encoding="utf-8")
    replacements = derive(gas, yaml.safe_load(template_text))
    output = replace_coefficients(template_text, replacements)
    yaml.safe_load(output)
    output_path.write_text(output, encoding="utf-8")
    print(f"Derived {len(replacements) // 3} condensates in {output_path}")


if __name__ == "__main__":
    main()
