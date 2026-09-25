"""Plot Gibbs-energy and heat-capacity correction magnitudes at gas joins."""

import csv
import math
import os
from pathlib import Path

RESULTS = Path(__file__).resolve().parent / "results"
cache = RESULTS / ".matplotlib"
cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator


def plot_histogram(rows, property_name, column, unit, output_name):
    lower = []
    higher = []
    for row in rows:
        magnitude = abs(float(row[column]))
        if not math.isfinite(magnitude):
            raise ValueError(f"Nonfinite {property_name} repair for {row['species']}")
        (lower if float(row["boundary_K"]) == 298.0 else higher).append(magnitude)

    positive = np.array([value for value in lower + higher if value > 0])
    if not len(positive):
        raise ValueError(f"No positive {property_name} repairs to plot")
    first = math.floor(math.log10(positive.min()))
    last = math.ceil(math.log10(positive.max()))
    edges = 10.0 ** np.arange(first, last + 0.5, 0.5)
    zero_count = len(rows) - len(positive)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist([lower, higher], bins=edges, stacked=True,
            color=["#7ea6c5", "#d9863b"], edgecolor="white", linewidth=0.7,
            label=[f"298 K joins (n={len(lower)})", f"Higher-temperature joins (n={len(higher)})"])
    ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel(f"Absolute {property_name} change to the adjusted fit at the join ({unit})")
    ax.set_ylabel("Number of joins per half-decade bin")
    ax.set_title(f"Size of {property_name} repairs at all {len(rows)} gas joins", weight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")

    largest = max(rows, key=lambda row: abs(float(row[column])))
    largest_value = abs(float(largest[column]))
    largest_label = (f"{largest_value / 1000:.1f} kJ/mol" if property_name == "Gibbs-energy"
                     else f"{largest_value:.2f} {unit}")
    ax.text(0.98, 0.74,
            f"Largest: {largest['species']} at {float(largest['boundary_K']):g} K\n"
            f"{largest_label}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
    fig.text(0.5, 0.015,
             f"Each value is |{property_name}(repaired adjusted segment) - "
             f"{property_name}(original adjusted segment)| at that join. Source: pinned v0.3.2 gas fits.",
             ha="center", va="bottom", fontsize=9)
    if zero_count:
        fig.text(0.02, 0.015, f"Exact zeros: {zero_count} (outside log axis)",
                 ha="left", va="bottom", fontsize=9)
    fig.subplots_adjust(bottom=0.17, top=0.89, left=0.11, right=0.97)
    output = RESULTS / output_name
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Wrote {output} ({len(rows)} joins)")


def main():
    with (RESULTS / "gas_joins.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 181:
        raise ValueError(f"Expected 181 gas joins, found {len(rows)}")
    plot_histogram(rows, "Gibbs-energy", "repair_delta_G_J_mol", "J/mol",
                   "gibbs_join_repair_histogram.png")
    plot_histogram(rows, "heat-capacity", "repair_delta_Cp_J_mol_K", "J/(mol K)",
                   "heat_capacity_join_repair_histogram.png")


if __name__ == "__main__":
    main()
