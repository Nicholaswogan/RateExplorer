"""Plot the ten largest gas-fit corrections for each thermodynamic property.

Run after validate.py, or use run_repair.py to regenerate everything.
"""

import csv
import os
from pathlib import Path

_matplotlib_cache = Path(__file__).resolve().parent / "results" / ".matplotlib"
_matplotlib_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_matplotlib_cache))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from source_snapshot import REPAIRED_GAS, RESULTS, original_text
from validate_independent import SPECIES, fetch, nasa_properties, parse_nasa9


GAS_SOURCE = "photochem_clima_data/data/reaction_mechanisms/zahnle_earth.yaml"
PLOTS = {
    "enthalpy": ("Enthalpy", "H", "kJ/mol", "H_J_mol", 1000.0, 80.0),
    "entropy": ("Entropy", "S", "J/(mol K)", "S_J_mol_K", 1.0, 80.0),
    "heat_capacity": ("Heat capacity", "Cp", "J/(mol K)", "Cp_J_mol_K", 1.0, 35.0),
    "gibbs": ("Gibbs energy", "G", "kJ/mol", "G_J_mol", 1000.0, 80.0),
}
REFERENCE_SPECIES = {
    name: (label, marker) for name, (label, marker, _) in SPECIES.items()
}
REFERENCE_SPECIES.update({
    "C2H3": ("NASA CEA", "C2H3,vinyl"),
    "C2H4": ("NASA CEA", "C2H4              TRC(4/88)"),
    "C2H5": ("NASA CEA", "C2H5              Chen,1990."),
    "CN": ("NASA CEA", "CN                Hf: Huang,1992."),
    "CH": ("NASA CEA", "CH                Gurvich,1979"),
    "N2H4": ("NASA CEA", "N2H4              Gurvich,1989"),
    "HO2": ("NASA CEA", "HO2               Hf:Hills,1984"),
    "H2SO4": ("NASA CEA", "H2SO4             Gurvich,1989"),
    "S8": ("NASA CEA", "S8                Gurvich,1989"),
    "NH": ("NASA CEA", "NH                Hf:Anderson,1989"),
})
REFERENCE_COLORS = {"NASA CEA": "#743ca0", "Burcat": "#17865a"}
PROPERTY_INDEX = {"H": 0, "S": 1, "G": 2, "Cp": 3}


def property_values(row, temperatures, symbol):
    a, b, c, d, e, f, g = row
    t = temperatures / 1000.0
    h = 1000.0 * (a * t + b * t**2 / 2 + c * t**3 / 3 + d * t**4 / 4 - e / t + f)
    s = a * np.log(t) + b * t + c * t**2 / 2 + d * t**3 / 3 - e / (2 * t**2) + g
    if symbol == "H":
        return h
    if symbol == "S":
        return s
    if symbol == "G":
        return h - temperatures * s
    return a + b * t + c * t**2 + d * t**3 + e / t**2


def top_corrections(key):
    with (RESULTS / "gas_curve_changes.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    largest = {}
    column = f"max_abs_delta_{key}"
    for row in rows:
        name = row["species"]
        if name not in largest or float(row[column]) > float(largest[name][column]):
            largest[name] = row
    return sorted(largest.values(), key=lambda row: (-float(row[column]), row["species"]))[:10]


def draw_segments(ax, thermo, symbol, scale, color, linestyle, zorder):
    ranges = thermo["temperature-ranges"]
    for index, row in enumerate(thermo["data"]):
        low, high = ranges[index:index + 2]
        low = max(low, 200.0)
        if low >= high:
            continue
        temperatures = np.linspace(low, high, 220)
        ax.plot(temperatures, property_values(row, temperatures, symbol) / scale,
                color=color, linestyle=linestyle, linewidth=1.7, zorder=zorder)


def draw_nasa9(ax, reference, symbol, scale):
    if reference is None:
        return
    label, fits = reference
    for low, high, row in fits:
        low = max(low, 298.0)
        high = min(high, 6000.0)
        if low >= high:
            continue
        temperatures = np.linspace(low, high, 220)
        values = np.array([nasa_properties(row, float(t))[PROPERTY_INDEX[symbol]]
                           for t in temperatures]) / scale
        ax.plot(temperatures, values, color=REFERENCE_COLORS[label],
                linestyle=":", linewidth=1.7, zorder=4)


def nasa9_visible_in_zoom(reference, boundary, half_width, symbol, scale, limits):
    _, fits = reference
    for low, high, row in fits:
        low = max(low, boundary - half_width)
        high = min(high, boundary + half_width)
        if low >= high:
            continue
        temperatures = np.linspace(low, high, 80)
        values = np.array([nasa_properties(row, float(t))[PROPERTY_INDEX[symbol]]
                           for t in temperatures]) / scale
        if np.any((limits[0] <= values) & (values <= limits[1])):
            return True
    return False


def local_limits(old, new, boundary, segment_index, half_width, symbol, scale):
    values = []
    for thermo in (old, new):
        ranges = thermo["temperature-ranges"]
        for index in (segment_index - 1, segment_index):
            low = max(boundary - half_width, ranges[index], 200.0)
            high = min(boundary + half_width, ranges[index + 1])
            if low < high:
                temperatures = np.linspace(low, high, 100)
                values.extend(property_values(thermo["data"][index], temperatures, symbol) / scale)
    low, high = min(values), max(values)
    padding = max((high - low) * 0.12, abs(high) * 0.002, 0.02)
    return low - padding, high + padding


def render_plot(original, repaired, references, stem, config):
    title, symbol, units, key, scale, half_width = config
    largest = top_corrections(key)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(5, 4, figsize=(20, 18))
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.075, top=0.89,
                        wspace=0.26, hspace=0.48)
    for rank, item in enumerate(largest):
        name = item["species"]
        segment_index = int(item["segment_index"])
        old = original[name]["thermo"]
        new = repaired[name]["thermo"]
        reference = references.get(name)
        boundary = float(new["temperature-ranges"][segment_index])
        full = axes[rank // 2, 2 * (rank % 2)]
        zoom = axes[rank // 2, 2 * (rank % 2) + 1]
        for ax in (full, zoom):
            draw_segments(ax, new, symbol, scale, "#1764a0", "-", 2)
            draw_segments(ax, old, symbol, scale, "#d17a14", "--", 3)
            draw_nasa9(ax, reference, symbol, scale)
            ax.axvline(boundary, color="0.55", linewidth=0.8, linestyle=":", zorder=1)
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel(f"{title} ({units})")
            ax.grid(alpha=0.2, linewidth=0.5)
            ax.tick_params(labelsize=8)
        full.set_xlim(200, 6000)
        zoom.set_xlim(boundary - half_width, boundary + half_width)
        limits = local_limits(old, new, boundary, segment_index, half_width, symbol, scale)
        zoom.set_ylim(limits)
        full.set_title(f"{rank + 1}. {name} - full range", fontweight="bold")
        zoom.set_title(f"{name} - join at {boundary:g} K", fontweight="bold")
        before = old["data"]
        at_boundary = np.array([boundary])
        jump = (property_values(before[segment_index], at_boundary, symbol)[0]
                - property_values(before[segment_index - 1], at_boundary, symbol)[0]) / scale
        max_change = float(item[f"max_abs_delta_{key}"]) / scale
        zoom.text(0.03, 0.04,
                  f"Max correction: {max_change:.3g} {units}\nOriginal join: {jump:+.3g} {units}",
                  transform=zoom.transAxes, ha="left", va="bottom", fontsize=8,
                  zorder=10, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
        reference_note = None
        if reference is None:
            reference_note = "No matched NASA9 fit"
        elif not nasa9_visible_in_zoom(reference, boundary, half_width, symbol, scale, limits):
            reference_note = "NASA9 outside zoom"
        if reference_note:
            zoom.text(0.97, 0.96, reference_note, transform=zoom.transAxes,
                      ha="right", va="top", fontsize=8,
                      zorder=10, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
    handles = [
        plt.Line2D([0], [0], color="#d17a14", linestyle="--", linewidth=2, label="Original"),
        plt.Line2D([0], [0], color="#1764a0", linestyle="-", linewidth=2, label="Repaired"),
        plt.Line2D([0], [0], color=REFERENCE_COLORS["NASA CEA"], linestyle=":",
                   linewidth=2, label="NASA CEA NASA9"),
        plt.Line2D([0], [0], color=REFERENCE_COLORS["Burcat"], linestyle=":",
                   linewidth=2, label="Burcat NASA9"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=4,
               frameon=False)
    fig.suptitle(f"Ten species with the largest {title.lower()} corrections",
                 fontsize=16, fontweight="bold", y=0.957)
    fig.text(0.5, 0.025,
             "Ranked by maximum |repaired - original| over each species' fitted range. "
             "Dotted vertical lines mark the corrected join. NASA9 curves use SHA-256-verified "
             "NASA CEA and Burcat sources. Absolute properties; no offsets applied.",
             ha="center", va="top", fontsize=8)
    output = RESULTS / f"{stem}_corrections.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(f"Wrote {output}")


def render_nasa9_gibbs_residuals(original, repaired, references):
    fig, axes = plt.subplots(3, 2, figsize=(12, 11), sharex=True)
    for ax, (name, (_, _, join)) in zip(axes.flat, SPECIES.items()):
        label, fits = references[name]
        temperatures = np.linspace(join + 1, 6000, 500)
        reference = np.array([
            nasa_properties(next(row for low, high, row in fits if low <= t <= high), float(t))[2]
            for t in temperatures
        ]) / 1000.0
        for version, mechanism, color, style in (
            ("Original", original, "#d17a14", "--"),
            ("Repaired", repaired, "#1764a0", "-"),
        ):
            values = property_values(mechanism[name]["thermo"]["data"][2],
                                     temperatures, "G") / 1000.0
            residual = values - reference
            rms = np.sqrt(np.mean(residual**2))
            ax.plot(temperatures, residual, color=color, linestyle=style,
                    linewidth=1.8, label=f"{version}: RMS {rms:.2f} kJ/mol")
        ax.axhline(0, color="0.45", linewidth=0.8)
        ax.set_title(f"{name} ({label})")
        ax.set_ylabel("Gibbs residual (kJ/mol)")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("Temperature (K)")
    fig.suptitle("Gibbs energy minus independent NASA9 fit", fontsize=15, weight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    output = RESULTS / "nasa9_gibbs_residuals.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(f"Wrote {output}")


def main():
    original = {entry["name"]: entry for entry in yaml.safe_load(original_text(GAS_SOURCE))["species"]}
    repaired = {entry["name"]: entry for entry in yaml.safe_load(
        REPAIRED_GAS.read_text(encoding="utf-8"))["species"]}
    needed = ({row["species"] for config in PLOTS.values()
               for row in top_corrections(config[3])} | SPECIES.keys())
    matched = {name: REFERENCE_SPECIES[name] for name in needed if name in REFERENCE_SPECIES}
    sources = {label: fetch(label) for label in {label for label, _ in matched.values()}}
    references = {name: (label, parse_nasa9(sources[label], marker))
                  for name, (label, marker) in matched.items()}
    for stem, config in PLOTS.items():
        render_plot(original, repaired, references, stem, config)
    render_nasa9_gibbs_residuals(original, repaired, references)


if __name__ == "__main__":
    main()
