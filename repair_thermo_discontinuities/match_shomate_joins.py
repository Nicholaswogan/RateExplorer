"""Make gas-phase Shomate H, S, G, and Cp continuous at every join.

The original 298 K-to-higher-temperature segment anchors each species. Match
its lower-temperature neighbor at 298 K, then match successive higher segments.
Only A, F, and G change in the adjusted segments; all other coefficients stay fixed.
Run on a gas-phase master YAML; phase-transition fits need separate treatment.

Usage: python match_shomate_joins.py INPUT.yaml OUTPUT.yaml
"""

import math
import re
import sys
from pathlib import Path

import yaml


def enthalpy_entropy(coefficients, temperature):
    a, b, c, d, e, f, g = coefficients
    t = temperature / 1000.0
    enthalpy = 1000.0 * (a * t + b * t**2 / 2 + c * t**3 / 3 + d * t**4 / 4 - e / t + f)
    entropy = a * math.log(t) + b * t + c * t**2 / 2 + d * t**3 / 3 - e / (2 * t**2) + g
    return enthalpy, entropy


def heat_capacity(coefficients, temperature):
    a, b, c, d, e = coefficients[:5]
    t = temperature / 1000.0
    return a + b * t + c * t**2 + d * t**3 + e / t**2


def corrected_rows(mechanism):
    replacements = {}
    for species in mechanism["species"]:
        thermo = species.get("thermo", {})
        if thermo.get("model") != "Shomate":
            continue
        rows = thermo["data"]
        ranges = thermo["temperature-ranges"]
        if len(rows) == 1:
            continue
        if ranges[1] != 298.0:
            raise ValueError(f"Expected a 298 K join for {species['name']}")

        # Preserve the polynomial that the 10-298 K approximation was built to meet.
        lower, anchor = rows[:2]
        temperature = ranges[1]
        lower[0] += heat_capacity(anchor, temperature) - heat_capacity(lower, temperature)
        anchor_h, anchor_s = enthalpy_entropy(anchor, temperature)
        lower_h, _ = enthalpy_entropy(lower, temperature)
        lower[5] += (anchor_h - lower_h) / 1000.0
        _, lower_s = enthalpy_entropy(lower, temperature)
        lower[6] += anchor_s - lower_s
        replacements[(species["name"], 0)] = lower

        for index in range(1, len(rows) - 1):
            temperature = ranges[index + 1]
            left, right = thermo["data"][index : index + 2]
            right[0] += heat_capacity(left, temperature) - heat_capacity(right, temperature)
            left_h, left_s = enthalpy_entropy(left, temperature)
            right_h, _ = enthalpy_entropy(right, temperature)
            right[5] += (left_h - right_h) / 1000.0
            _, right_s = enthalpy_entropy(right, temperature)
            right[6] += left_s - right_s
            replacements[(species["name"], index + 1)] = right
    return replacements


def replace_coefficients(text, replacements):
    output = []
    species = None
    data_index = 0
    in_data = False
    changed = set()
    lines = iter(text.splitlines(keepends=True))
    for line in lines:
        match = re.match(r"- name: (.+)", line)
        if match:
            species = match.group(1).strip().strip("'\"")
            in_data = False
        if line.startswith("    data:"):
            in_data = True
            data_index = 0
        elif in_data and line.startswith("    - ["):
            key = (species, data_index)
            if key in replacements:
                while "]" not in line:
                    try:
                        line += next(lines)
                    except StopIteration as error:
                        raise ValueError(f"Unclosed Shomate data for {key}") from error
                before, rest = line.split("[", 1)
                contents, after = rest.split("]", 1)
                old = contents.split(",")
                if len(old) != 7:
                    raise ValueError(f"Expected one-line Shomate data for {key}")
                old[0] = format(replacements[key][0], ".15g")
                old[5] = " " + format(replacements[key][5], ".15g")
                old[6] = " " + format(replacements[key][6], ".15g")
                line = before + "[" + ",".join(old) + "]" + after
                changed.add(key)
            data_index += 1
        output.append(line)
    if changed != replacements.keys():
        raise ValueError(f"Could not serialize {set(replacements) - changed}")
    return "".join(output)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python match_shomate_joins.py INPUT.yaml OUTPUT.yaml")
    source, destination = map(Path, sys.argv[1:])
    if source.resolve() == destination.resolve():
        raise SystemExit("Input and output must differ")
    text = source.read_text(encoding="utf-8")
    mechanism = yaml.safe_load(text)
    replacements = corrected_rows(mechanism)
    revised = replace_coefficients(text, replacements)
    yaml.safe_load(revised)  # Check serialization before writing.
    destination.write_text(revised, encoding="utf-8")
    print(f"Matched Cp, H, S, and G at {len(replacements)} joins in {destination}")


if __name__ == "__main__":
    main()
