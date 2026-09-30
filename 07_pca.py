"""
07_pca.py
=========

PCA for Figure 1 (species + lifestyle) and Supplemental Figure S2
(species + location). Recipe taken verbatim from the original ChatGPT
analysis chat that produced the published figures.

Method:
    X = samples x proteins matrix of spectral counts
    X_log = log2(X + 1)
    X_centered = X_log - X_log.mean(axis=0)   # protein-wise mean center
    sklearn PCA(n_components=2).fit_transform(X_centered)

Aesthetics:
    species color = matplotlib tab10 colormap, keyed to alphabetical
                    species order so the same species gets the same
                    color in both figures.
    lifestyle shape = o (Sessile colonial), s (Sessile non-colonial),
                      ^ (Motile colonial), D (Motile non-colonial).
                      Matches the published Figure 1 legend.
    location shape  = 'o', 's', '^', 'D', 'P', 'X', 'v', '<', '>',
                      'h', '8', 'p', 'd', '*', 'H' in the location
                      order from the chat (Album first).

Lifestyle labels for Figure 1 use "non-colonial" (matching the
published figure legend); the supplementals S4/S5/S6 use "Solitary"
in the published figures. samples_metadata.csv carries the "Solitary"
form in lifestyle_4 and we translate to "non-colonial" inline here so
both label conventions coexist without diverging from the manuscript.

INPUTS  : processed_protein_data.xlsx (Single Proteins), samples_metadata.csv
OUTPUTS : Figure_1_rotifer_pca_species_lifestyle.{png,pdf,svg,tif}
          Supplemental_Figure_S2_rotifer_pca_species_location.{png,pdf,svg,tif}
          rotifer_pca_species_lifestyle_coordinates.csv
          rotifer_pca_coordinates_corrected_colors.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (load_full_singles, load_metadata, sample_columns,
                    strip_group_prefix, PROJECT_ROOT, italic_species)


# Alphabetical species ordering used for color assignment (chat's recipe)
PCA_SPECIES_ORDER = [
    "Conochilus hippocrepis",
    "Cupelopagis vorax",
    "Euchlanis kingi",
    "Hexarthra sp.",
    "Lacinularia flosculosa",
    "Notommata copeus",
    "Plationus patulus",
    "Sinantherina socialis",
]

LIFESTYLE_ORDER = [
    "Sessile colonial",
    "Sessile non-colonial",
    "Motile colonial",
    "Motile non-colonial",
]
LIFESTYLE_MARKERS = {
    "Sessile colonial":     "o",
    "Sessile non-colonial": "s",
    "Motile colonial":      "^",
    "Motile non-colonial":  "D",
}

# Location order + markers (chat's recipe verbatim)
LOCATION_ORDER = [
    "Album", "Lake Nockamixon", "Mantua", "Moon Lake", "Pines Lake",
    "Mixed", "Luisa", "Tip Top Hueco", "Laguna Prieta", "Rattlesnake",
    "McKinney", "Ojo de la Casa", "Ojo de la Punta", "Lab", "Krystal Lake",
]
MARKER_CYCLE = ['o', 's', '^', 'D', 'P', 'X', 'v', '<', '>',
                'h', '8', 'p', 'd', '*', 'H']
LOCATION_MARKERS = dict(zip(LOCATION_ORDER, MARKER_CYCLE))

# tab10 colormap keyed to PCA_SPECIES_ORDER
TAB10 = plt.get_cmap("tab10").colors
SPECIES_COLOR_MAP = {sp: TAB10[i] for i, sp in enumerate(PCA_SPECIES_ORDER)}


def _lifestyle_pca_label(row: pd.Series) -> str:
    """Translate metadata's 'Solitary' terminology to the PCA-figure
    'non-colonial' terminology used in the published Figure 1 legend."""
    sessile = row["lifestyle_2"]       # "Sessile" or "Motile"
    colonial = row["colonial"]         # "Colonial" or "Solitary"
    cn = "colonial" if colonial == "Colonial" else "non-colonial"
    return f"{sessile} {cn}"


def compute_pca(df: pd.DataFrame, meta: pd.DataFrame):
    cols = sample_columns(df, meta)
    sample_names = [strip_group_prefix(c) for c in cols]
    X = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).T.values
    X_log = np.log2(X + 1.0)
    X_centered = X_log - X_log.mean(axis=0)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_centered)
    return sample_names, coords, pca.explained_variance_ratio_


def _species_full_from_meta(meta: pd.DataFrame) -> dict[str, str]:
    return dict(zip(meta["sample_name"], meta["species_full"]))


def _location_from_meta(meta: pd.DataFrame) -> dict[str, str]:
    return dict(zip(meta["sample_name"], meta["location"]))


def _lifestyle_pca_from_meta(meta: pd.DataFrame) -> dict[str, str]:
    return {row["sample_name"]: _lifestyle_pca_label(row)
            for _, row in meta.iterrows()}


def plot_pca(*, sample_names, coords, var_ratio, species_lookup,
             shape_lookup, shape_map, shape_order, shape_legend_title,
             title, out_stem, out_csv):
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    rows = []
    for nm, (x, y) in zip(sample_names, coords):
        sp = species_lookup.get(nm)
        sh = shape_lookup.get(nm)
        if sp is None or sh is None:
            continue
        rows.append({"Sample": nm, "Species": sp, "ShapeKey": sh,
                     "PCA1": x, "PCA2": y})
        ax.scatter(x, y,
                   color=SPECIES_COLOR_MAP.get(sp, "#444"),
                   marker=shape_map.get(sh, "o"),
                   s=100, alpha=0.9,
                   edgecolor="k", linewidth=0.5)
    pd.DataFrame(rows).to_csv(out_csv, index=False)

    # No in-figure title: the figure number and description live in the
    # typeset caption (journal style). `title` is kept in the signature so
    # callers are unchanged.
    ax.set_xlabel(f"PCA1 ({var_ratio[0]*100:.1f}% variance)", fontsize=14)
    ax.set_ylabel(f"PCA2 ({var_ratio[1]*100:.1f}% variance)", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.3)

    species_handles = [
        Line2D([0], [0], marker='o', linestyle='None',
               markerfacecolor=SPECIES_COLOR_MAP[s],
               markeredgecolor=SPECIES_COLOR_MAP[s],
               markersize=9, label=italic_species(s))
        for s in PCA_SPECIES_ORDER
    ]
    legend1 = ax.legend(handles=species_handles, title="Species (Color)",
                         bbox_to_anchor=(1.02, 1.0), loc="upper left",
                         borderaxespad=0.0, fontsize=12, title_fontsize=13)
    ax.add_artist(legend1)
    # legend1 is re-registered below via bbox_extra_artists so that
    # bbox_inches="tight" reserves room for it instead of clipping it.

    shape_handles = [
        Line2D([0], [0], marker=shape_map[v], linestyle='None',
               color='gray', markerfacecolor='gray',
               markeredgecolor='gray', markersize=9, label=v)
        for v in shape_order
        if any(s["ShapeKey"] == v for s in rows)
    ]
    legend2 = ax.legend(handles=shape_handles, title=shape_legend_title,
                        bbox_to_anchor=(1.02, 0.52), loc="upper left",
                        borderaxespad=0.0, fontsize=12, title_fontsize=13)

    # No tight_layout() here: it resizes the axes without accounting for
    # the two out-of-axes legends, which is what truncated the species
    # names in the previously published Figure 1 / Sup S2.
    for ext in ("png", "pdf", "svg", "tif"):
        fig.savefig(out_stem.with_suffix("." + ext), dpi=300,
                    bbox_inches="tight",
                    bbox_extra_artists=(legend1, legend2))
    plt.close(fig)
    print(f"  wrote {out_stem.name}.*")


def main():
    meta = load_metadata()
    df = load_full_singles()
    names, coords, vr = compute_pca(df, meta)

    species_lookup = _species_full_from_meta(meta)
    lifestyle_lookup = _lifestyle_pca_from_meta(meta)
    location_lookup = _location_from_meta(meta)

    plot_pca(
        sample_names=names, coords=coords, var_ratio=vr,
        species_lookup=species_lookup,
        shape_lookup=lifestyle_lookup, shape_map=LIFESTYLE_MARKERS,
        shape_order=LIFESTYLE_ORDER,
        shape_legend_title="Lifestyle (Shape)",
        title="PCA of Rotifer Species\nColor = Species, Shape = Lifestyle",
        out_stem=PROJECT_ROOT / "Figure_1_rotifer_pca_species_lifestyle",
        out_csv=PROJECT_ROOT
                / "rotifer_pca_species_lifestyle_coordinates.csv",
    )

    plot_pca(
        sample_names=names, coords=coords, var_ratio=vr,
        species_lookup=species_lookup,
        shape_lookup=location_lookup, shape_map=LOCATION_MARKERS,
        shape_order=LOCATION_ORDER,
        shape_legend_title="Location (Shape)",
        title="PCA of Rotifer Species\nColor = Species, Shape = Location",
        out_stem=(PROJECT_ROOT
                  / "Supplemental_Figure_S2_rotifer_pca_species_location"),
        out_csv=PROJECT_ROOT
                / "rotifer_pca_coordinates_corrected_colors.csv",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
