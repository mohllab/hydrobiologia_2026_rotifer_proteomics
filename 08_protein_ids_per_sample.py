"""
08_protein_ids_per_sample.py
============================

Supplemental Figure S1 - bar chart of unique proteins identified per
sample, grouped by species with gaps between species blocks and dashed
separator lines. Species names are placed beneath the x-axis using
ax.get_xaxis_transform().

Recipe taken verbatim from the original ChatGPT analysis chat.

A protein is identified for a sample when its weighted spectral count
in that sample column is > 0.

INPUTS  : processed_protein_data.xlsx (Single Proteins), samples_metadata.csv
OUTPUTS : Supplemental_Figure_S1_rotifer_protein_ids_by_species_location.*
          rotifer_protein_ids_by_species.csv
"""

from __future__ import annotations

import sys
from pathlib import Path
from textwrap import fill

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (load_full_singles, load_metadata, sample_columns,
                    strip_group_prefix, PROJECT_ROOT, italic_species)


# Chat's species ordering for Sup S1
S1_SPECIES_ORDER = [
    "Notommata copeus",
    "Sinantherina socialis",
    "Euchlanis kingi",
    "Conochilus hippocrepis",
    "Lacinularia flosculosa",
    "Plationus patulus",
    "Cupelopagis vorax",
    "Hexarthra sp.",
]


def main():
    meta = load_metadata()
    df = load_full_singles()
    cols = sample_columns(df, meta)

    rows = []
    name_to_full = dict(zip(meta["sample_name"], meta["species_full"]))
    name_to_loc = dict(zip(meta["sample_name"], meta["location"]))
    name_to_num = dict(zip(meta["sample_name"], meta["sample_number"]))
    for c in cols:
        nm = strip_group_prefix(c)
        if nm not in name_to_full:
            continue
        col = pd.to_numeric(df[c], errors="coerce").fillna(0)
        rows.append({
            "Sample": nm,
            "Species": name_to_full[nm],
            "Location": name_to_loc[nm],
            "SampleNumber": int(name_to_num[nm]),
            "Proteins Identified": int((col > 0).sum()),
        })
    df_c = pd.DataFrame(rows)

    # Order: by S1 species order, then by sample number
    df_c["_species_order"] = df_c["Species"].map(
        {s: i for i, s in enumerate(S1_SPECIES_ORDER)})
    df_c = (df_c.sort_values(["_species_order", "Location", "Sample"])
            .reset_index(drop=True))

    # Build grouped x positions (gap=1 between species blocks)
    x_positions, x_labels = [], []
    species_centers = {}
    current_x = 0
    for sp in S1_SPECIES_ORDER:
        sub = df_c[df_c["Species"] == sp]
        if sub.empty:
            continue
        positions = list(range(current_x, current_x + len(sub)))
        x_positions.extend(positions)
        # Two-line tick label: collection site, then the sample number as
        # an explicit "sample N" sub-label so the number cannot be read as
        # a count. N is the biological-sample number used in the PRIDE
        # deposit and in the CSV written below.
        x_labels.extend([f"{r['Location']}\nsample {r['SampleNumber']}"
                         for _, r in sub.iterrows()])
        species_centers[sp] = sum(positions) / len(positions)
        current_x += len(sub) + 1

    plot_df = (pd.DataFrame({"x": x_positions, "_label": x_labels})
                .assign(_idx=lambda d: d.index)
                .reset_index(drop=True))
    plot_df = pd.concat([plot_df, df_c.reset_index(drop=True)], axis=1)

    csv_path = PROJECT_ROOT / "rotifer_protein_ids_by_species.csv"
    df_c[["Sample", "Species", "Location", "Proteins Identified"]] \
        .to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(17, 7.5), dpi=300)
    ax.bar(plot_df["x"], plot_df["Proteins Identified"],
           color="#1976D2", edgecolor="white", linewidth=0.6)
    ax.set_ylabel("Proteins identified", fontsize=14)
    ax.set_xlabel("Collection site (sample number)", fontsize=14)
    # No in-figure title: it collided with the species labels above the
    # axes, and the caption carries the description.
    ax.set_xticks(plot_df["x"])
    ax.set_xticklabels(plot_df["_label"], rotation=75, ha="right",
                        fontsize=12)
    ax.tick_params(axis="y", labelsize=12)

    # Species labels above the axes
    for sp, center in species_centers.items():
        # Every species label on two lines — genus above, epithet (or an
        # upright "sp.") below — so the labels are uniform across the panel.
        genus, epithet = sp.split(" ", 1)
        label = italic_species(genus) + "\n" + (
            epithet if epithet == "sp." else italic_species(epithet))
        ax.text(center, 1.01, label,
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=13)

    # Dashed separator lines between species
    last_species = [s for s in S1_SPECIES_ORDER if s in species_centers]
    for sp in last_species[:-1]:
        boundary = plot_df.loc[plot_df["Species"] == sp, "x"].max() + 0.5
        ax.axvline(boundary, linestyle="--", linewidth=0.8,
                   color="#1976D2", alpha=0.6)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.subplots_adjust(bottom=0.35)
    plt.tight_layout()

    stem = (PROJECT_ROOT /
            "Supplemental_Figure_S1_rotifer_protein_ids_by_species_location")
    for ext in ("png", "pdf", "svg", "tif"):
        fig.savefig(stem.with_suffix("." + ext), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem.name}.* (n={len(df_c)}, "
          f"range={df_c['Proteins Identified'].min()}-"
          f"{df_c['Proteins Identified'].max()})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
