"""
06a_validate_subclass_counts.py
================================

Cross-checks the generated data files against the numbers stated in the
manuscript. Run this after 09_top10_functional.py (which writes both
CSVs) and before 06_subclass_panels.py.

The previous version of this script read three files that the pipeline
never produces —
    rotifer_mucin_category_by_species.csv
    rotifer_species_functional_summary_with_mucin_contractile.csv
    rotifer_species_unique_protein_presence_with_mucin_contractile.csv
      (with a column "Total Unique Proteins in Species")
— so it could not run. It now reads the files 09 actually writes and the
column names they actually carry.

Inputs:
  - rotifer_species_functional_summary_full.csv          (09)
  - rotifer_species_unique_protein_presence_with_mucin_contractile.csv (09)
  - Total_Protein_Output_052725.xlsx                      (source matrix)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _subclass_rules import subclass_table  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SUMMARY_CSV = PROJECT_ROOT / "rotifer_species_functional_summary_full.csv"
PRESENCE_CSV = (PROJECT_ROOT /
                "rotifer_species_unique_protein_presence_with_mucin_contractile.csv")
MATRIX_XLSX = PROJECT_ROOT / "Total_Protein_Output_052725.xlsx"

SPECIES_ORDER = [
    'Cupelopagis vorax', 'Sinantherina socialis', 'Lacinularia flosculosa',
    'Notommata copeus', 'Conochilus hippocrepis', 'Plationus patulus',
    'Euchlanis kingi', 'Hexarthra sp.',
]

# Values as printed in the manuscript
MS_SPECIES_TOTALS = [540, 971, 624, 959, 957, 692, 672, 432]
MS_SUBCLASS_TOTALS = {
    'Mucin / glycan-associated / surface': [0, 8, 10, 4, 2, 3, 1, 1],
    'Contractile / muscle-associated':     [11, 21, 12, 22, 15, 14, 9, 3],
    'Cytoskeleton / motility':             [105, 120, 96, 111, 109, 101, 93, 82],
}
MS_TOTAL_PROTEINS = 1675
MS_SAMPLE_RANGE = (178, 959)

fails = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    if not ok:
        fails.append(label)


print("=" * 70)
print("VALIDATION: do the generated data files match the manuscript?")
print("=" * 70)

# ---- CHECK 1: per-species totals -------------------------------------
print("\n--- CHECK 1: unique proteins per species ---")
summary = pd.read_csv(SUMMARY_CSV)
totals = summary.groupby("Species")["Total Proteins in Species"].first()
got = [int(totals[sp]) for sp in SPECIES_ORDER]
print(f"{'Species':<25}{'computed':>10}{'manuscript':>12}")
for sp, g, m in zip(SPECIES_ORDER, got, MS_SPECIES_TOTALS):
    print(f"{sp:<25}{g:>10}{m:>12}")
check("per-species totals", got == MS_SPECIES_TOTALS)

# ---- CHECK 2: categories sum to the total ----------------------------
print("\n--- CHECK 2: category counts sum to the species total ---")
sums = summary.groupby("Species")["Protein Count"].sum()
ok = all(int(sums[sp]) == int(totals[sp]) for sp in SPECIES_ORDER)
check("category sums == totals", ok)

# ---- CHECK 3: presence CSV row counts --------------------------------
print("\n--- CHECK 3: presence CSV agrees with the summary ---")
pres = pd.read_csv(PRESENCE_CSV)
ok = True
for cat, expected in MS_SUBCLASS_TOTALS.items():
    rows = (pres[pres.Category == cat].groupby("Species").size()
            .reindex(SPECIES_ORDER).fillna(0).astype(int).tolist())
    check(f"presence rows, {cat}", rows == expected, f"{rows}")
    ok &= rows == expected

# ---- CHECK 4: subclass rules reproduce the published totals ----------
print("\n--- CHECK 4: keyword subclass rules reproduce published totals ---")
table = subclass_table(PRESENCE_CSV)
for cat, expected in MS_SUBCLASS_TOTALS.items():
    got_sc = [sum(table[cat].get(sp, {}).values()) for sp in SPECIES_ORDER]
    check(f"subclass totals, {cat}", got_sc == expected, f"{got_sc}")

# ---- CHECK 5: total non-redundant proteins ---------------------------
print("\n--- CHECK 5: total non-redundant proteins ---")
mat = pd.read_excel(MATRIX_XLSX, sheet_name="For_Patrick")
n_acc = int(mat["Accession Number"].notna().sum())
check(f"1,675 non-redundant proteins", n_acc == MS_TOTAL_PROTEINS,
      f"computed {n_acc} (sheet has {len(mat)} rows; the last row is "
      f"Scaffold's per-sample summary footer and carries no accession)")

# ---- CHECK 6: per-sample identification range ------------------------
print("\n--- CHECK 6: proteins identified per sample ---")
ids = pd.read_excel(MATRIX_XLSX, sheet_name="Proteins Identified")
vals = pd.to_numeric(ids.iloc[:, 1], errors="coerce").dropna().astype(int)
check("range 178-959", (vals.min(), vals.max()) == MS_SAMPLE_RANGE,
      f"computed {vals.min()}-{vals.max()}")
check("20 samples", len(vals) == 20, f"computed {len(vals)}")

print("\n" + "=" * 70)
if fails:
    print(f"VALIDATION FAILED — {len(fails)} check(s): " + "; ".join(fails))
    sys.exit(1)
print("VALIDATION COMPLETE — all checks passed")
print("=" * 70)
