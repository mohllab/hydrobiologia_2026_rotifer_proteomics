"""
11_database_robustness.py
=========================

Analyses supporting the response to Reviewer 1 (database- and
annotation-dependent bias). Nothing here changes a manuscript figure;
the script documents the numbers quoted in the response letter and
writes the species-by-source-genome table (Supplemental Table S2).

Four questions are answered:

  1. Which reference genome supplied each identification?  For every
     study species, the number and share of its unique proteins drawn
     from each of the thirteen genomes in the search database.
  2. Is the number of proteins identified for a species a function of
     the size of its predicted proteome?  (Pearson / Spearman, n = 7
     species with a conspecific genome.)
  3. Are the BH-FDR-significant proteins of Figure 3 skewed towards the
     species whose own gene models are in the database?  Direction and
     source genome of every red point in panels A-D.
  4. Does the class-level enrichment survive when quantitative values
     are summed across ALL homologous entries of a class (which removes
     any effect of a peptide being assigned to one homologue rather
     than another), and when only functionally annotated entries are
     used?

Inputs (project root):
  - Total_Protein_Output_052725.xlsx                        (For_Patrick)
  - Code_Supplement/samples_metadata.csv
  - rotifer_species_unique_protein_presence_with_mucin_contractile.csv (09)
  - Volcano (1)/Volcano/<comparison>/processed_protein_data (NN).xlsx

Outputs (project root):
  - Supplemental_Table_S2_species_by_source_genome.csv
  - database_robustness_summary.txt
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MATRIX_XLSX = PROJECT_ROOT / "Total_Protein_Output_052725.xlsx"
META_CSV = Path(__file__).resolve().parent / "samples_metadata.csv"
PRESENCE_CSV = (PROJECT_ROOT /
                "rotifer_species_unique_protein_presence_with_mucin_contractile.csv")
VOLCANO_DIR = PROJECT_ROOT / "Volcano (1)" / "Volcano"
OUT_TABLE = PROJECT_ROOT / "Supplemental_Table_S2_species_by_source_genome.csv"
OUT_TXT = PROJECT_ROOT / "database_robustness_summary.txt"

# Accession prefix -> genome. The four marked "unpublished" are the
# laboratory genomes first used here; their entries carry no functional
# descriptor in the search FASTA (protein name == accession).
GENOMES = {
    "notorw": "Notommata copeus",
    "cono":   "Conochilus hippocrepis",     # unpublished
    "ssocau": "Sinantherina socialis",
    "euki":   "Euchlanis kingi",
    "laflo":  "Lacinularia flosculosa",
    "ppat":   "Plationus patulus",
    "hex":    "Hexarthra sp.",
    "bv":     "Brachionus variabilis",      # unpublished
    "pquad":  "Platyias quadricornis",      # unpublished
    "asp":    "Asplanchna girodi",
    "echi":   "Epiphanes chihuahuaensis",
    "fl":     "Filinia longiseta",          # unpublished
    "ebra":   "Epiphanes brachionus",
}
UNPUBLISHED = {"Conochilus hippocrepis", "Brachionus variabilis",
               "Platyias quadricornis", "Filinia longiseta"}

# Predicted protein-coding gene models per genome, counted from the
# deposited search FASTA (unique gene models, secondary isoforms
# collapsed). These are the denominators used in Table S1.
PREDICTED = {
    "Notommata copeus": 33943, "Conochilus hippocrepis": 18071,
    "Sinantherina socialis": 22563, "Euchlanis kingi": 28646,
    "Lacinularia flosculosa": 31550, "Plationus patulus": 36604,
    "Hexarthra sp.": 11093, "Brachionus variabilis": 33306,
    "Platyias quadricornis": 47699, "Asplanchna girodi": 22211,
    "Epiphanes chihuahuaensis": 18205, "Filinia longiseta": 20843,
    "Epiphanes brachionus": 21474,
}

SPECIES_ORDER = [
    "Cupelopagis vorax", "Sinantherina socialis", "Lacinularia flosculosa",
    "Notommata copeus", "Conochilus hippocrepis", "Plationus patulus",
    "Euchlanis kingi", "Hexarthra sp.",
]

# Figure 3 panels: (folder, (0) group, (1) group). log2FC is
# log2(mean(1)+1) - log2(mean(0)+1), so negative = enriched in the (0)
# group (the first-named species).
PANELS = {
    "A": ("Cupelopagis vorax vs Plationus patulus", "Cupelopagis vorax", "Plationus patulus"),
    "B": ("Cupelopagis vorax vs Sinantherina socialis", "Cupelopagis vorax", "Sinantherina socialis"),
    "C": ("Sinantherina socialis vs Plationus patulus", "Plationus patulus", "Sinantherina socialis"),
    "D": ("Lacinularia flosculosa vs Sinantherina socialis", "Lacinularia flosculosa", "Sinantherina socialis"),
}

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def genome_of(acc: pd.Series) -> pd.Series:
    return acc.astype(str).str.extract(r"^([A-Za-z]+)")[0].map(GENOMES)


def find_legacy(folder: Path) -> Path:
    cands = sorted(folder.glob("processed_protein_data (*).xlsx"))
    if not cands:
        sys.exit(f"no processed_protein_data (NN).xlsx in {folder}")
    return cands[0]


# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------
m = pd.read_excel(MATRIX_XLSX, sheet_name="For_Patrick")
m = m[m["Accession Number"].notna()].copy()          # drops Scaffold footer
meta = pd.read_csv(META_CSV)
sp_of = dict(zip(meta.sample_name, meta.species_full))
samples = [c for c in m.columns if c in sp_of]
assert len(samples) == 20, samples
m["genome"] = genome_of(m["Accession Number"])
assert m["genome"].notna().all(), "unmapped accession prefix"
X = m[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)

pres = pd.read_csv(PRESENCE_CSV).drop_duplicates("Accession Number")
m["Category"] = m["Accession Number"].map(
    dict(zip(pres["Accession Number"], pres["Category"])))
assert m["Category"].notna().all(), "protein missing from presence CSV"
# An entry is "unannotated" when the keyword classifier (script 09) put it
# in the Unknown / unannotated class: every entry from the four unpublished
# sets (name == accession), plus published-genome entries whose BLASTP
# descriptor is blank ("-") or "hypothetical protein ...".
m["unannotated_entry"] = m["Category"] == "Unknown / unannotated"

present = pd.DataFrame({
    sp: (X[[s for s in samples if sp_of[s] == sp]] > 0).any(axis=1)
    for sp in SPECIES_ORDER})

say("=" * 72)
say("DATABASE / ANNOTATION ROBUSTNESS  (Reviewer 1, major concern 1)")
say("=" * 72)
say(f"{len(m):,} proteins; {len(samples)} samples; "
    f"{len(GENOMES)} genomes in the search database")

# ----------------------------------------------------------------------
# 1. species x source genome
# ----------------------------------------------------------------------
say("\n--- 1. Source genome of each species' unique proteins ---")
ct = pd.DataFrame({sp: m.loc[present[sp].values, "genome"].value_counts()
                   for sp in SPECIES_ORDER}).reindex(list(GENOMES.values())).fillna(0).astype(int)
pct = ct / ct.sum() * 100
say("counts:")
say(ct.to_string())
say("\npercent of the species' unique proteins:")
say(pct.round(1).to_string())

table = ct.copy()
table.index.name = "Source genome"
table["Predicted gene models"] = [PREDICTED[g] for g in table.index]
table["Status"] = ["unpublished, no functional descriptor" if g in UNPUBLISHED
                   else "Mohl et al. (2025)" for g in table.index]
table.to_csv(OUT_TABLE)
say(f"\nwrote {OUT_TABLE.name}")

say("\nself-genome share (proteins identified from the species' own genome):")
rows = []
for sp in SPECIES_ORDER:
    n = int(ct[sp].sum())
    if sp in ct.index:
        k = int(ct.loc[sp, sp])
        rows.append((sp, n, k, 100 * k / n, PREDICTED[sp], ct[sp].idxmax(), 100 * ct[sp].max() / n))
    else:
        rows.append((sp, n, np.nan, np.nan, np.nan, ct[sp].idxmax(), 100 * ct[sp].max() / n))
self_df = pd.DataFrame(rows, columns=["species", "unique_proteins", "conspecific_hits",
                                      "self_share_pct", "predicted_gene_models",
                                      "largest_source", "largest_source_pct"])
say(self_df.round(1).to_string(index=False))
own = self_df.dropna()
say(f"\n{int((own.largest_source == own.species).sum())} of {len(own)} species with a "
    "conspecific genome draw their largest share from it.")

unpub_share = {sp: 100 * int(ct.loc[list(UNPUBLISHED), sp].sum()) / int(ct[sp].sum())
               for sp in SPECIES_ORDER}
say("\nshare of unique proteins identified from the four unpublished (unannotated) sets:")
for sp in SPECIES_ORDER:
    say(f"  {sp:<24} {ct.loc[list(UNPUBLISHED), sp].sum():>4} / {ct[sp].sum():<4} = {unpub_share[sp]:.1f}%")
n_unpub = int(m["genome"].isin(UNPUBLISHED).sum())
say(f"  all proteins: {n_unpub} / {len(m)} = {100 * n_unpub / len(m):.1f}%")
n_unk = int((m["Category"] == "Unknown / unannotated").sum())
n_dash = int((~m["genome"].isin(UNPUBLISHED) & m["unannotated_entry"]
              & (m["Protein Name"].astype(str).str.strip() == "-")).sum())
say(f"  'Unknown / unannotated' class: {n_unk} / {len(m)} = {100 * n_unk / len(m):.1f}% "
    f"({n_unpub} from the unpublished sets + {n_unk - n_unpub} published-genome entries "
    f"whose descriptor is 'hypothetical protein' [{n_unk - n_unpub - n_dash}] or blank [{n_dash}])")

# ----------------------------------------------------------------------
# 2. detected vs predicted
# ----------------------------------------------------------------------
say("\n--- 2. Proteins detected vs predicted proteome size (n = 7) ---")
for col, label in [("conspecific_hits", "conspecific identifications"),
                   ("unique_proteins", "all unique proteins of the species")]:
    r, p = stats.pearsonr(own[col], own["predicted_gene_models"])
    rho, ps = stats.spearmanr(own[col], own["predicted_gene_models"])
    say(f"  {label:<36} Pearson r = {r:.2f} (p = {p:.2f}); Spearman rho = {rho:.2f} (p = {ps:.2f})")

# ----------------------------------------------------------------------
# 3. Figure 3 red points: direction and source genome
# ----------------------------------------------------------------------
say("\n--- 3. BH-FDR + |log2FC| >= 2 proteins of Figure 3: direction and source ---")
tot_third = tot_red = tot_self = 0
for panel, (folder, s0, s1) in PANELS.items():
    legacy = find_legacy(VOLCANO_DIR / folder)
    sig = pd.read_excel(legacy, sheet_name="Sig with correction")
    full = pd.read_excel(legacy, sheet_name="Single Proteins")
    g0 = [c for c in full.columns if isinstance(c, str) and re.match(r"^\(0\)", c)]
    g1 = [c for c in full.columns if isinstance(c, str) and re.match(r"^\(1\)", c)]
    a = full[g0].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    b = full[g1].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    full["log2FC"] = np.log2(b + 1) - np.log2(a + 1)
    full["genome"] = genome_of(full["Accession Number"])
    red = full[full["Accession Number"].isin(set(sig["Accession Number"]))
               & (full["log2FC"].abs() >= 2)].copy()
    say(f"\nPanel {panel}: (0) {s0} vs (1) {s1}: BH-significant {sig['Accession Number'].nunique()}, "
        f"with |log2FC| >= 2: {len(red)}")
    if red.empty:
        continue
    red["enriched_in"] = np.where(red["log2FC"] < 0, s0, s1)
    red["third_party"] = ~red["genome"].isin({s0, s1})
    red["own_genome"] = red["genome"] == red["enriched_in"]
    for sp in (s0, s1):
        say(f"  enriched in {sp:<24} {int((red.enriched_in == sp).sum()):>3}")
    say(f"  identified from a genome of neither compared species: {int(red.third_party.sum())} of {len(red)}")
    say(f"  identified from the genome of the species it is enriched in: {int(red.own_genome.sum())} of {len(red)}")
    say("  source genome x direction:")
    say("  " + pd.crosstab(red["genome"], red["enriched_in"]).to_string().replace("\n", "\n  "))
    tot_red += len(red); tot_third += int(red.third_party.sum()); tot_self += int(red.own_genome.sum())
say(f"\nAll panels: {tot_red} proteins; {tot_third} from third-party genomes; "
    f"{tot_self} from the genome of the species in which they are enriched.")

# ----------------------------------------------------------------------
# 4. class-level shares per sample
# ----------------------------------------------------------------------
say("\n--- 4. Class share of each sample's total quantitative value ---")
say("(summed across every database entry of the class, so independent of which")
say(" homologous entry a peptide was assigned to)")
tot_all = X.sum()
tot_annot = X[~m["unannotated_entry"]].sum()
unannot_share = (X[m["unannotated_entry"]].sum() / tot_all * 100)
say("\nshare of signal in Unknown / unannotated entries (four unpublished sets + hypothetical/blank descriptors):")
u = pd.DataFrame({"species": [sp_of[s] for s in samples], "pct": unannot_share.values})
say(u.groupby("species")["pct"].agg(["mean", "min", "max"]).reindex(SPECIES_ORDER).round(1).to_string())


def welch(d, sp_a, sp_b):
    a = d[d.species == sp_a]["pct"]; b = d[d.species == sp_b]["pct"]
    t = stats.ttest_ind(a, b, equal_var=False)
    return f"Welch t = {t.statistic:.2f}, p = {t.pvalue:.4f}"


for cls in ["Cytoskeleton / motility", "Contractile / muscle-associated"]:
    for denom_name, denom in [("all entries", tot_all), ("entries with an informative descriptor", tot_annot)]:
        share = X[m["Category"] == cls].sum() / denom * 100
        d = pd.DataFrame({"species": [sp_of[s] for s in samples], "pct": share.values})
        g = d.groupby("species")["pct"].agg(["mean", "std", "count"]).reindex(SPECIES_ORDER).round(1)
        say(f"\n{cls} — % of {denom_name}:")
        say(g.to_string())
        for sp in ["Cupelopagis vorax", "Plationus patulus", "Sinantherina socialis"]:
            say(f"  {sp:<24} replicates: {np.round(d[d.species == sp].pct.values, 1)}")
        say(f"  C. vorax vs P. patulus:  {welch(d, 'Cupelopagis vorax', 'Plationus patulus')}")
        say(f"  C. vorax vs S. socialis: {welch(d, 'Cupelopagis vorax', 'Sinantherina socialis')}")

OUT_TXT.write_text("\n".join(lines) + "\n")
say(f"\nwrote {OUT_TXT.name}")
