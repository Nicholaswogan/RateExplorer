"""Numerically validate both generated YAMLs and write reviewable metrics."""

import csv
import json
import math
from pathlib import Path

import numpy as np
import yaml

from derive_condensates import R_J, WINDOWS, gas_row, gibbs, sat_log_pressure
from match_shomate_joins import enthalpy_entropy, heat_capacity
from source_snapshot import REPAIRED_CONDENSATE, REPAIRED_GAS, RESULTS, original_text


HERE = Path(__file__).resolve().parent
GAS_SOURCE = "photochem_clima_data/data/reaction_mechanisms/zahnle_earth.yaml"
CONDENSATE_SOURCE = "photochem_clima_data/data/reaction_mechanisms/condensate_thermo.yaml"


def properties(row, temperature):
    h, s = enthalpy_entropy(row, temperature)
    return h, s, h - temperature * s, heat_capacity(row, temperature)


def write_csv(name, rows):
    with (RESULTS / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    original = yaml.safe_load(original_text(GAS_SOURCE))
    repaired = yaml.safe_load(REPAIRED_GAS.read_text())
    condensate_old = yaml.safe_load(original_text(CONDENSATE_SOURCE))
    condensate_new = yaml.safe_load(REPAIRED_CONDENSATE.read_text())
    old_gases = {entry["name"]: entry for entry in original["species"]}
    gases = {entry["name"]: entry for entry in repaired["species"]}
    if original["particles"] != repaired["particles"] or original["reactions"] != repaired["reactions"]:
        raise AssertionError("Gas generation altered particle or reaction records")
    if old_gases.keys() != gases.keys():
        raise AssertionError("Gas species inventory changed")

    joins = []
    curves = []
    for name, species in gases.items():
        old_thermo = old_gases[name].get("thermo", {})
        thermo = species.get("thermo", {})
        if thermo.get("model") != "Shomate":
            if old_thermo != thermo:
                raise AssertionError(f"Unexpected non-Shomate change in {name}")
            continue
        if old_thermo["temperature-ranges"] != thermo["temperature-ranges"]:
            raise AssertionError(f"Changed temperature ranges in {name}")
        anchor_index = 1 if len(thermo["data"]) > 1 else 0
        for index, (before, after) in enumerate(zip(old_thermo["data"], thermo["data"])):
            if before[1:5] != after[1:5] or (index == anchor_index and before != after):
                raise AssertionError(f"Changed anchored/shape coefficients in {name}")
            lower, upper = thermo["temperature-ranges"][index:index + 2]
            samples = np.linspace(max(lower, 10.0), upper, 101)
            effects = np.array([
                np.subtract(properties(after, float(t)), properties(before, float(t)))
                for t in samples
            ])
            curves.append({
                "species": name, "segment_index": index, "low_K": lower, "high_K": upper,
                "max_abs_delta_H_J_mol": float(np.max(np.abs(effects[:, 0]))),
                "max_abs_delta_S_J_mol_K": float(np.max(np.abs(effects[:, 1]))),
                "max_abs_delta_G_J_mol": float(np.max(np.abs(effects[:, 2]))),
                "max_abs_delta_Cp_J_mol_K": float(np.max(np.abs(effects[:, 3]))),
            })
        for index, boundary in enumerate(thermo["temperature-ranges"][1:-1]):
            left = properties(thermo["data"][index], boundary)
            right = properties(thermo["data"][index + 1], boundary)
            original_left = properties(old_thermo["data"][index], boundary)
            original_right = properties(old_thermo["data"][index + 1], boundary)
            changed = left if index == 0 else right
            original_changed = original_left if index == 0 else original_right
            joins.append({
                "species": name, "boundary_K": boundary,
                **{f"delta_{property}": right[i] - left[i]
                   for i, property in enumerate(("H_J_mol", "S_J_mol_K", "G_J_mol", "Cp_J_mol_K"))},
                "original_delta_G_J_mol": original_right[2] - original_left[2],
                "repair_delta_G_J_mol": changed[2] - original_changed[2],
                "original_delta_Cp_J_mol_K": original_right[3] - original_left[3],
                "repair_delta_Cp_J_mol_K": changed[3] - original_changed[3],
            })
    if len(joins) != 181:
        raise AssertionError(f"Expected 181 gas joins; found {len(joins)}")
    max_gas = {key: max(abs(row[key]) for row in joins) for key in joins[0] if key.startswith("delta_")}
    if max_gas["delta_G_J_mol"] > 1e-6 or max_gas["delta_H_J_mol"] > 1e-6 or max_gas["delta_S_J_mol_K"] > 1e-8 or max_gas["delta_Cp_J_mol_K"] > 1e-8:
        raise AssertionError(f"Gas joins not continuous: {max_gas}")
    write_csv("gas_joins.csv", joins)
    write_csv("gas_curve_changes.csv", curves)

    old_condensates = {entry["name"]: entry for entry in condensate_old["species"]}
    particles = {entry["name"]: entry for entry in repaired["particles"] if entry.get("formation") == "saturation"}
    results = []
    for entry in condensate_new["species"]:
        particle_name = entry["name"]
        gas_name = particle_name.removesuffix("aer")
        saturation = particles[particle_name]["saturation"]
        thermo = entry["thermo"]
        old_thermo = old_condensates[particle_name]["thermo"]
        if entry["composition"] != old_condensates[particle_name]["composition"] or thermo["temperature-ranges"] != old_thermo["temperature-ranges"]:
            raise AssertionError(f"Condensate metadata changed for {particle_name}")
        gas_thermo = gases[gas_name]["thermo"]
        triple, critical = thermo["temperature-ranges"][1:3]
        low, high = WINDOWS[gas_name]
        low = low if low is not None else triple - 30
        boundaries = []
        for index, temperature in enumerate((triple, critical)):
            left = properties(thermo["data"][index], temperature)
            right = properties(thermo["data"][index + 1], temperature)
            boundaries.append(right[2] - left[2])
        errors = []
        latent_errors = []
        old_errors = []
        for index, (start, stop) in enumerate(((low, triple), (triple, critical), (critical, high))):
            samples = np.linspace(start, stop, 301)
            log_errors = []
            old_log_errors = []
            latent = []
            for t in samples:
                t = float(t)
                gas_row_at_t = gas_row(gas_thermo, t)
                gas_g = gibbs(gas_row_at_t, t)
                goal = sat_log_pressure(saturation, t)
                log_errors.append((gibbs(thermo["data"][index], t) - gas_g) / (R_J * t) - goal)
                old_log_errors.append((gibbs(old_thermo["data"][index], t) - gas_g) / (R_J * t) - goal)
                if index < 2:
                    gas_h = properties(gas_row_at_t, t)[0]
                    condensed_h = properties(thermo["data"][index], t)[0]
                    branch = "sublimation" if index == 0 else "vaporization"
                    sat_heat = (saturation[branch]["a"] + saturation[branch]["b"] * t) * saturation["parameters"]["mu"] / 1e7
                    latent.append((gas_h - condensed_h) - sat_heat)
            errors.append(float(np.max(np.abs(log_errors))) / math.log(10))
            old_errors.append(float(np.max(np.abs(old_log_errors))) / math.log(10))
            if latent:
                latent_errors.append(float(np.max(np.abs(latent))))
        # The historical supercritical fit window ends at 1000 or 1500 K.
        # Audit its extrapolation separately rather than treating it as a
        # physically measured saturation curve at higher temperatures.
        extrapolated = []
        for t in np.linspace(high, 6000, 501):
            t = float(t)
            log_p = (gibbs(thermo["data"][2], t) - gibbs(gas_row(gas_thermo, t), t)) / (R_J * t)
            extrapolated.append((log_p - sat_log_pressure(saturation, t)) / math.log(10))
        results.append({
            "species": particle_name,
            "triple_K": triple, "critical_K": critical,
            "G_jump_triple_J_mol": boundaries[0],
            "G_jump_critical_J_mol": boundaries[1],
            "max_abs_log10_P_error_solid": errors[0],
            "max_abs_log10_P_error_liquid": errors[1],
            "max_abs_log10_P_error_supercritical": errors[2],
            "old_max_abs_log10_P_error_solid": old_errors[0],
            "old_max_abs_log10_P_error_liquid": old_errors[1],
            "old_max_abs_log10_P_error_supercritical": old_errors[2],
            "max_abs_log10_P_error_supercritical_extrapolated_to_6000K": float(np.max(np.abs(extrapolated))),
            "max_abs_latent_heat_error_below_critical_J_mol": max(latent_errors),
        })
    if len(results) != 13:
        raise AssertionError(f"Expected 13 condensates; found {len(results)}")
    max_condensate_jump = max(abs(row[key]) for row in results for key in ("G_jump_triple_J_mol", "G_jump_critical_J_mol"))
    max_subcritical_pressure_error = max(row[key] for row in results for key in ("max_abs_log10_P_error_solid", "max_abs_log10_P_error_liquid"))
    if max_condensate_jump > 1e-6 or max_subcritical_pressure_error > 0.03:
        raise AssertionError("Condensate Gibbs or physical saturation fit failed")
    write_csv("condensate_metrics.csv", results)
    summary = {
        "gas_join_count": len(joins),
        "gas_max_abs_join_deltas": max_gas,
        "gas_max_abs_curve_G_change_J_mol": max(row["max_abs_delta_G_J_mol"] for row in curves),
        "condensate_count": len(results),
        "condensate_max_abs_G_join_J_mol": max_condensate_jump,
        "condensate_max_abs_log10_P_error_below_critical": max_subcritical_pressure_error,
        "condensate_max_abs_log10_P_error_supercritical": max(row["max_abs_log10_P_error_supercritical"] for row in results),
        "condensate_max_abs_log10_P_error_supercritical_extrapolated_to_6000K": max(row["max_abs_log10_P_error_supercritical_extrapolated_to_6000K"] for row in results),
        "condensate_max_abs_latent_heat_error_below_critical_J_mol": max(row["max_abs_latent_heat_error_below_critical_J_mol"] for row in results),
    }
    (RESULTS / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
