"""
10_mucin_by_species.py
======================

Auxiliary chart (NOT a manuscript figure): per-species count of unique
mucin / glycan-associated / surface proteins, one bar per species. These
are the same eight values that the stacked Figure 4 prints above its
bars. An earlier revision carried this chart as Supplemental Figure S4;
it was dropped as redundant with Figure 4, and the contractile and
cytoskeleton panels became Supplemental Figures S4 and S5. The script is
kept because it asserts the eight counts independently of script 06.

Species run sessile -> motile, the same order used by Figure 4 and by
Supplemental Figures S4 and S5.

The previous version of this script called `species_subclass_table()`
and `GLYCAN_RULES` on 06_subclass_panels.py; neither name exists there,
so the script raised AttributeError and could not run. It now reads the
same computed subclass table that 06 uses.

Input:
  - rotifer_species_unique_protein_presence_with_mucin_contractile.csv
    (written by 09_top10_functional.py)

Outputs: Auxiliary_mucin_category_by_species.{png,pdf,svg,tif}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _subclass_rules import subclass_table  # noqa: E402
from common import italic_species  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

GLYCAN_CAT = "Mucin / glycan-associated / surface"

# Sessile -> motile, matching Figure 4 and Supplemental S4/S5
SPECIES_ORDER = [
    "Cupelopagis vorax",
    "Sinantherina socialis",
    "Lacinularia flosculosa",
    "Notommata copeus",
    "Conochilus hippocrepis",
    "Plationus patulus",
    "Euchlanis kingi",
    "Hexarthra sp.",
]

EXPECTED = [0, 8, 10, 4, 2, 3, 1, 1]


def main() -> int:
    table = subclass_table()[GLYCAN_CAT]
    vals = [sum(table.get(sp, {}).values()) for sp in SPECIES_ORDER]
    assert vals == EXPECTED, f"computed {vals}, published values are {EXPECTED}"

    fig, ax = plt.subplots(figsize=(12, 5.8), dpi=300)
    x = np.arange(len(SPECIES_ORDER))
    ax.bar(x, vals, color="#1976D2", width=0.7, edgecolor="white", linewidth=0.6)
    for xi, v in zip(x, vals):
        ax.text(xi, v + max(vals) * 0.02, str(int(v)),
                ha="center", va="bottom", fontsize=13)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [italic_species(sp) for sp in SPECIES_ORDER],
        rotation=20, ha="right", fontsize=13)
    ax.set_xlabel("Species", fontsize=14)
    ax.set_ylabel("Unique protein count", fontsize=14)
    ax.tick_params(axis="y", labelsize=12)
    # Title names the class as it is defined in the Methods, which is
    # broader than mucins alone — only one of the 22 unique proteins in
    # this class is a bona fide mucin.
    # No in-figure title: the caption carries it.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    stem = PROJECT_ROOT / "Auxiliary_mucin_category_by_species"
    for ext in ("png", "pdf", "svg", "tif"):
        fig.savefig(stem.with_suffix("." + ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem.name}.*  counts="
          f"{dict(zip(SPECIES_ORDER, vals))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
