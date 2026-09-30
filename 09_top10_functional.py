"""
09_top10_functional.py
=======================

Figure 2 - Broad functional characterization of detected proteins by
rotifer species. 8-panel horizontal-bar layout, top-11 functional
categories per species + "Collapsed minor categories", total per
species labeled in the panel title.

Classifier matches the manuscript Methods text verbatim:
    "stress-response proteins are commonly discussed in rotifer biology,
    heat shock/chaperones and oxidative stress/antioxidant defense were
    treated as separate focal classes"
    "mucin/glycan-associated and contractile/muscle-associated were
    added to capture surface and mechanical functions"
    "top 11 most represented classes per species were displayed
    individually, whereas the remaining low-count classes were combined
    into Collapsed minor categories"

13 focal classes (first match wins) + "Annotated proteins outside focal
classes" (default) + "Unknown / unannotated" (no name / gene-model-style
name / "uncharacterized protein" / "predicted protein").

Recipe ported verbatim from the original ChatGPT analysis chat
(lines 1238-1430 of the captured chat scripts).

INPUTS  : processed_protein_data.xlsx (Single Proteins), samples_metadata.csv
OUTPUTS : Figure_2_rotifer_species_functional_annotation_panels_top10.*
          rotifer_species_functional_summary_full.csv
          rotifer_species_functional_summary_top_categories.csv
          rotifer_species_unique_protein_presence_with_mucin_contractile.csv
            (one row per species x detected protein, with the focal-class
             assignment; consumed by 06/06a/06b/10)
"""

from __future__ import annotations

import re
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
from common import (load_full_singles, load_metadata,
                    sample_columns, strip_group_prefix, PROJECT_ROOT, italic_species)


# Per-species bar colour. Keyed to the SAME alphabetical species order
# and tab10 colormap that 07_pca.py uses for Figure 1 and Sup S2, so a
# species keeps one colour across every figure that encodes species by
# colour. (Defined here rather than imported because 07_pca.py's module
# name starts with a digit.)
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

# Manuscript Fig 2 species ordering (verbatim from chat)
FIG2_SPECIES_ORDER = [
    "Cupelopagis vorax",
    "Lacinularia flosculosa",
    "Sinantherina socialis",
    "Hexarthra sp.",
    "Conochilus hippocrepis",
    "Euchlanis kingi",
    "Plationus patulus",
    "Notommata copeus",
]


def clean_name(name: str) -> str:
    """Strip Scaffold-export decorators that pollute the name."""
    s = str(name).strip()
    s = re.sub(r"^(PREDICTED:|LOW QUALITY PROTEIN:)\s*", "", s, flags=re.I)
    return s.strip()


def is_accession_like(name: str) -> bool:
    """Recognize gene-model-style names like 'asp13300.t1', 'notorw5207.t1'."""
    s = str(name).strip()
    patterns = [
        r"^hypothetical protein\b",
        r"^[A-Za-z]{1,6}\d+(?:\.\w+)?$",
        r"^[A-Za-z]{1,12}\d+\.t\d+$",
        r"^[A-Za-z]+DRAFT[_-]?\d+$",
        r"^[A-Za-z0-9]+_[0-9]{3,}$",
        r"^[A-Za-z]{2,}\d+_[0-9]{3,}$",
        r"^[a-z]{1,10}\d+\.t\d+$",
    ]
    return any(re.match(p, s) for p in patterns)


# 13-category classifier, verbatim from chat lines 1276-1359
def classify_protein(name: str) -> str:
    s = str(name).strip().lower()
    if s in {"", "-", "nan", "unknown"}:
        return "Unknown / unannotated"

    category_patterns = [
        ("Heat shock / chaperones", [
            r"heat shock", r"\bhsp\d*", r"stress-70", r"hsc70",
            r"chaperon", r"dnaj", r"glucose-regulated protein",
            r"\bgrp78\b", r"\bgrp94\b",
        ]),
        ("Oxidative stress / antioxidant defense", [
            r"glutathione", r"superoxide dismutase", r"\bsod\b",
            r"catalase", r"peroxiredoxin", r"thioredoxin",
            r"glutaredoxin", r"peroxidase", r"thioredoxin reductase",
            r"cytochrome p450", r"\bp450\b",
        ]),
        ("Mucin / glycan-associated / surface", [
            r"mucin", r"galectin", r"lectin", r"proteoglycan",
            r"heparan sulfate", r"chondroitin", r"ependymin",
            r"von willebrand", r"glycosyltransferase",
        ]),
        ("Contractile / muscle-associated", [
            r"titin", r"twitchin", r"tektin", r"rootletin", r"troponin",
            r"paramyosin", r"myomesin", r"obscurin", r"muscle lim",
            r"sarcoplasmic calcium-binding", r"striated muscle",
        ]),
        ("Translation / ribosome / RNA", [
            r"ribosomal", r"ribosome", r"\btrna\b", r"\brrna\b",
            r"translation", r"elongation factor", r"initiation factor",
            r"splic", r"rna binding", r"rna-binding", r"nucleolar",
            r"hnrnp", r"snrnp", r"dead-box", r"rna helicase",
        ]),
        ("Cytoskeleton / motility", [
            r"actin", r"tubulin", r"dynein", r"kinesin", r"myosin",
            r"microtubule", r"flagell", r"cilia", r"cili", r"centrin",
            r"filamin", r"profilin", r"cofilin", r"spectrin",
            r"intermediate filament", r"axonem", r"tropomyosin",
        ]),
        ("Metabolism / energy", [
            r"dehydrogenase", r"synthase", r"synthetase", r"isomerase",
            r"mutase", r"oxidase", r"reductase", r"lyase", r"aldolase",
            r"hydratase", r"fatty acid", r"beta-oxidation", r"glycol",
            r"mitochond", r"atp synthase", r"pyruvate", r"enolase",
            r"citrate synthase", r"fructose", r"glucose",
            r"phosphoglycerate", r"carboxylase", r"adenylate kinase",
            r"arginine kinase", r"transketolase",
            r"glycogen phosphorylase", r"adenosylhomocysteinase",
        ]),
        ("Proteostasis / degradation", [
            r"proteasome", r"ubiquitin", r"peptidase", r"protease",
            r"cathepsin", r"proteolysis", r"autophagy", r"calpain",
            r"lysosomal", r"foldase",
        ]),
        ("Signaling / regulation", [
            r"14-3-3", r"phosphatase", r"gtpase", r"\bras\b",
            r"\brho\b", r"\brab\b", r"\barf\b", r"calmodulin",
            r"wd repeat", r"ankyrin", r"leucine-rich",
            r"transcription factor", r"regulator", r"cyclin",
            r"checkpoint", r"adenylyl cyclase-associated",
            r"guanine nucleotide-binding",
        ]),
        ("Transport / membrane trafficking", [
            r"transporter", r"channel", r"\batpase\b", r"\babc\b",
            r"solute carrier", r"membrane", r"vesicle", r"golgi",
            r"endoplasmic reticulum", r"clathrin", r"coatomer",
            r"snare", r"vacuolar", r"endocyt", r"exocyt",
            r"proton pump", r"\bsec61\b",
        ]),
        ("DNA / chromatin / replication", [
            r"histone", r"dna ", r"dna-", r"chromatin", r"replication",
            r"repair", r"polymerase", r"nucleosome", r"topoisomerase",
            r"telomere",
        ]),
        ("Extracellular / structural", [
            r"extracellular", r"cuticle", r"chitin", r"collagen",
            r"laminin", r"fibronectin", r"cadherin", r"adhesion",
            r"matrix", r"keratin",
        ]),
    ]

    for category, patterns in category_patterns:
        if any(re.search(p, s) for p in patterns):
            return category

    if is_accession_like(s):
        return "Unknown / unannotated"
    return "Annotated proteins outside focal classes"


def main():
    meta = load_metadata()
    df = load_full_singles()

    full_to_cols: dict[str, list[str]] = {}
    name_to_full = dict(zip(meta["sample_name"], meta["species_full"]))
    for c in sample_columns(df, meta):
        sp_full = name_to_full.get(strip_group_prefix(c))
        if sp_full:
            full_to_cols.setdefault(sp_full, []).append(c)

    # Use Base Protein Name (split on '|') and clean prefixes
    df = df.copy()
    df["Base Protein Name"] = (
        df["Protein Name"].fillna("Unknown").astype(str)
          .str.split("|").str[0].str.strip()
    )
    df["Clean Protein Name"] = df["Base Protein Name"].apply(clean_name)
    df["Category"] = df["Clean Protein Name"].apply(classify_protein)

    full_rows = []
    plot_rows = []
    presence_rows = []
    for sp in FIG2_SPECIES_ORDER:
        cols = full_to_cols.get(sp, [])
        if not cols:
            print(f"[WARN] no samples for {sp}")
            continue
        spectral = (df[cols].apply(pd.to_numeric, errors="coerce")
                    .fillna(0).sum(axis=1))
        present_mask = spectral > 0
        spdf = df.loc[present_mask, ["Protein Name", "Accession Number",
                                      "Clean Protein Name", "Category"]]
        total = len(spdf)
        for _, r in spdf.iterrows():
            presence_rows.append({
                "Species": sp,
                "Accession Number": r["Accession Number"],
                "Clean Protein Name": r["Clean Protein Name"],
                "Category": r["Category"],
            })
        counts = spdf["Category"].value_counts()
        for cat, n in counts.items():
            full_rows.append({
                "Species": sp, "Category": cat,
                "Protein Count": int(n),
                "Percent of Species Proteome": round(n / total * 100, 2),
                "Total Proteins in Species": total,
            })
        top = counts.head(11)
        remainder = int(counts.iloc[11:].sum())
        for cat, n in top.items():
            plot_rows.append({
                "Species": sp, "Category": cat,
                "Protein Count": int(n),
                "Total Proteins in Species": total,
            })
        if remainder > 0:
            plot_rows.append({
                "Species": sp, "Category": "Collapsed minor categories",
                "Protein Count": remainder,
                "Total Proteins in Species": total,
            })

    summary_df = pd.DataFrame(full_rows)
    plot_df = pd.DataFrame(plot_rows)
    summary_df.to_csv(
        PROJECT_ROOT / "rotifer_species_functional_summary_full.csv",
        index=False)
    plot_df.to_csv(
        PROJECT_ROOT / "rotifer_species_functional_summary_top_categories.csv",
        index=False)
    pd.DataFrame(presence_rows).to_csv(
        PROJECT_ROOT /
        "rotifer_species_unique_protein_presence_with_mucin_contractile.csv",
        index=False)

    tab10 = plt.get_cmap("tab10").colors
    species_color = {sp: tab10[i] for i, sp in enumerate(PCA_SPECIES_ORDER)}

    fig, axes = plt.subplots(4, 2, figsize=(19, 26), sharex=False, dpi=300)
    axes = axes.flatten()
    for ax, sp in zip(axes, FIG2_SPECIES_ORDER):
        sub = (plot_df[plot_df["Species"] == sp]
               .sort_values("Protein Count", ascending=True).copy())
        y = np.arange(len(sub))
        ax.barh(y, sub["Protein Count"],
                color=species_color.get(sp, "#1976D2"),
                edgecolor="white", linewidth=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels([fill(c, 34) for c in sub["Category"]],
                            fontsize=15)
        ax.set_xlabel("Unique protein count", fontsize=16)
        ax.tick_params(axis="x", labelsize=14)
        total = (int(sub["Total Proteins in Species"].iloc[0])
                 if len(sub) else 0)
        ax.set_title(f"{italic_species(sp)}\nTotal identified: {total}",
                     fontsize=19, pad=12)
        if len(sub):
            mx = sub["Protein Count"].max()
            for i, v in enumerate(sub["Protein Count"]):
                ax.text(v + max(2, mx * 0.01), i, str(int(v)),
                        va="center", fontsize=14)
            ax.set_xlim(0, mx * 1.22)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    # No figure-level title: the description lives in the typeset caption.
    plt.tight_layout()
    stem = (PROJECT_ROOT /
            "Figure_2_rotifer_species_functional_annotation_panels_top10")
    for ext in ("png", "pdf", "svg", "tif"):
        fig.savefig(stem.with_suffix("." + ext), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem.name}.*")
    return 0


if __name__ == "__main__":
    sys.exit(main())
