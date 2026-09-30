"""
06_subclass_panels.py
=====================

Stacked bar charts of protein subclass composition within the three
focal classes:

  - Figure 4:        Mucin / glycan-associated / surface subclasses
  - Supplemental S4: Contractile / muscle-associated subclasses
  - Supplemental S5: Cytoskeleton / motility subclasses

The mucin/glycan panel is the main-text Figure 4: the per-species totals
are printed above its bars, so the separate per-species count chart
(script 10) is no longer a manuscript figure.

Figure numbers are NOT drawn into the artwork — they belong in the
typeset caption, and burning a supplemental-figure number into what is
now Figure 4 was how the previous version got mislabelled.

Counts are COMPUTED at run time from
`rotifer_species_unique_protein_presence_with_mucin_contractile.csv`
using the keyword rules in 06b_subclassify_proteins.py, then asserted
against the published per-species totals before anything is drawn.
(The previous version hard-coded the counts; computing them keeps the
figures, the supplementary CSV and the classifier rules in lockstep, as
promised in the response to Reviewer 2.)

Input:
  - rotifer_species_unique_protein_presence_with_mucin_contractile.csv
    (written by 09_top10_functional.py)

Outputs (per figure): .png (400 dpi), .pdf (vector), .tiff (300 dpi LZW)

Dependencies:
  pip install pandas matplotlib pillow
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _subclass_rules import SUBCLASS_ORDER, subclass_table  # noqa: E402
from common import italic_species  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'pdf.fonttype': 42,      # embed TrueType so text stays selectable
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
})

# Type sizes. The published version used 8 pt tick labels and an 11 pt
# title on a 10 in figure; Reviewer 2 asked for larger labels in the
# supplemental graphs, so these are scaled up on a wider canvas.
FS_TICK, FS_AXIS, FS_TITLE, FS_LEGEND, FS_VALUE = 14, 17, 19, 14, 14

SPECIES_ORDER = [
    'Cupelopagis vorax', 'Sinantherina socialis', 'Lacinularia flosculosa',
    'Notommata copeus', 'Conochilus hippocrepis', 'Plationus patulus',
    'Euchlanis kingi', 'Hexarthra sp.',
]
SPECIES_SHORT = [
    'C. vorax', 'S. socialis', 'L. flosculosa', 'N. copeus',
    'C. hippocrepis', 'P. patulus', 'E. kingi', 'Hexarthra sp.',
]
LIFESTYLE = [
    'Sessile\nSolitary', 'Sessile\nColonial', 'Sessile\nColonial',
    'Motile\nSolitary', 'Motile\nColonial', 'Motile\nSolitary',
    'Motile\nSolitary', 'Motile\nSolitary',
]

# Published per-species totals — asserted, never used to draw.
EXPECTED_TOTALS = {
    'Mucin / glycan-associated / surface': [0, 8, 10, 4, 2, 3, 1, 1],
    'Contractile / muscle-associated':     [11, 21, 12, 22, 15, 14, 9, 3],
    'Cytoskeleton / motility':             [105, 120, 96, 111, 109, 101, 93, 82],
}

COLORS = {
    'Mucin / glycan-associated / surface':
        ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628',
         '#f781bf', '#999999'],
    'Contractile / muscle-associated':
        ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f',
         '#8c564b', '#e377c2', '#bcbd22', '#999999'],
    'Cytoskeleton / motility':
        ['#009E73', '#D55E00', '#7570b3', '#E7298A', '#F0E442', '#66a61e',
         '#E6AB02', '#8B4513', '#666666', '#1f78b4', '#b3de69'],
}

TITLES = {
    'Mucin / glycan-associated / surface':
        'Mucin / glycan-associated / surface protein subclasses '
        'across rotifer species',
    'Contractile / muscle-associated':
        'Contractile / muscle-associated protein subclasses '
        'across rotifer species',
    'Cytoskeleton / motility':
        'Cytoskeleton / motility protein subclasses across rotifer species',
}

STEMS = {
    'Mucin / glycan-associated / surface': 'Figure_4_Glycan_Subclasses',
    'Contractile / muscle-associated':     'Supplemental_Figure_S4_Contractile_Subclasses',
    'Cytoskeleton / motility':             'Supplemental_Figure_S5_Cytoskeleton_Subclasses',
}

FIGSIZE = {
    'Mucin / glycan-associated / surface': (17.0, 8.5),
    'Contractile / muscle-associated':     (17.0, 8.5),
    'Cytoskeleton / motility':             (17.0, 9.5),
}


def _species_labels() -> list[str]:
    """Italic abbreviated species name over a roman lifestyle label."""
    out = []
    for short, life in zip(SPECIES_SHORT, LIFESTYLE):
        out.append(italic_species(short) + "\n" + life)
    return out


def make_stacked_bar(data, subclasses, colors, title, stem, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    bottoms = [0] * len(SPECIES_ORDER)
    for sub, col in zip(subclasses, colors):
        vals = [data.get(sp, {}).get(sub, 0) for sp in SPECIES_ORDER]
        if sum(vals) == 0:
            continue                     # keep empty subclasses out of the legend
        ax.bar(range(len(SPECIES_ORDER)), vals, bottom=bottoms, label=sub,
               color=col, width=0.62, edgecolor='white', linewidth=0.6)
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    for i, total in enumerate(bottoms):
        ax.text(i, total + max(bottoms) * 0.015, str(int(total)),
                ha='center', va='bottom', fontsize=FS_VALUE, fontweight='bold')

    ax.set_xticks(range(len(SPECIES_ORDER)))
    ax.set_xticklabels(_species_labels(), fontsize=FS_TICK)
    ax.set_ylabel('Number of unique proteins', fontsize=FS_AXIS)
    ax.tick_params(axis='y', labelsize=FS_TICK)
    # No in-figure title (journal style: the caption carries it); `title`
    # is kept in the signature so callers are unchanged.
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=FS_LEGEND,
              frameon=True, borderpad=0.7, labelspacing=0.6)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.set_ylim(0, max(bottoms) * 1.10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    png = PROJECT_ROOT / f'{stem}.png'
    fig.savefig(png, dpi=400, bbox_inches='tight')
    fig.savefig(PROJECT_ROOT / f'{stem}.pdf', bbox_inches='tight')
    plt.close(fig)
    Image.open(png).save(PROJECT_ROOT / f'{stem}.tiff', format='TIFF',
                         dpi=(300, 300), compression='tiff_lzw')
    print(f"  saved {stem}.png / .pdf / .tiff   totals={[int(b) for b in bottoms]}")


def main() -> int:
    table = subclass_table()

    for cat, expected in EXPECTED_TOTALS.items():
        got = [sum(table[cat].get(sp, {}).values()) for sp in SPECIES_ORDER]
        assert got == expected, (
            f"{cat}: computed {got} but published totals are {expected}")
    print("All subclass totals reproduce the published per-species values.\n")

    for cat in ('Mucin / glycan-associated / surface',
                'Contractile / muscle-associated',
                'Cytoskeleton / motility'):
        make_stacked_bar(table[cat], SUBCLASS_ORDER[cat], COLORS[cat],
                         TITLES[cat], STEMS[cat], FIGSIZE[cat])
    print("\nAll figures generated successfully.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
