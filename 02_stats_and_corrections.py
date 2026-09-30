"""
02_stats_and_corrections.py
============================

Per-protein differential abundance testing between two rotifer species
using **Scaffold's pre-computed t-test column** as the source of truth
for p-values. Scaffold writes this column to the protein table under
the header "t-test, CL: Species, BRL: Scaffold Category" (or
"Permutation Test, ..." when n<3) and is the authoritative statistic
used in the manuscript's volcano plots.

Pipeline:
    "Single Proteins" sheet
      -> read Scaffold's t-test / permutation-test column
      -> parse "< 0.0001" as the literal 0.0001 boundary
         (matches Patrick's original Code.txt parser)
      -> write back as "Test Column" so 04_volcano_plot.py picks it up
      -> compute a reference Bonferroni-style threshold (alpha/n_tested)
         and write summary stats; this is NOT the threshold used for
         the published figures
      -> filter "Sig with/without correction" sheets

Note on multiple-test correction: the manuscript's published volcano
plots (Figure 3) use Scaffold's built-in Benjamini-Hochberg FDR at
q=0.05 (selected in Scaffold's "Configure Sample Organization and
Statistical Analysis" dialog). Patrick's per-pair Scaffold exports
already contain BH-FDR-filtered "Sig with correction" sheets and the
empirical FWER value in the Raw Extract — script 04 reads those
authoritative outputs directly, so this script's Bonferroni-style
summary is informational rather than a control on the figures.

Why we trust Scaffold's stat: the original analysis (and the
manuscript's volcano plots) was built on Scaffold's per-protein test,
which uses a weighted spectrum count test internally tuned to the
Scaffold normalization. A fresh Welch t-test on log-transformed counts
gives different (and slightly less precise) p-values; reproducing the
published figures requires using Scaffold's values.

INPUTS  : Volcano (1)/Volcano/<comparison>/processed_protein_data.xlsx
OUTPUTS : same workbook with updated "Single Proteins", added
          "Sig without correction", "Sig with correction"
          + statistics_summary.csv per comparison
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
VOLCANO_ROOT = PROJECT_ROOT / "Volcano (1)" / "Volcano"

CANONICAL_COMPARISONS = [
    "Cupelopagis vorax vs Plationus patulus",
    "Cupelopagis vorax vs Sinantherina socialis",
    "Sinantherina socialis vs Plationus patulus",
    "Lacinularia flosculosa vs Sinantherina socialis",
]


def parse_scaffold_p(x):
    """Scaffold writes p-values like 0.0034, '< 0.0001', or '--'.

    Convention used here matches Patrick's original `Code.txt` parser:
    '< 0.0001' is parsed as the literal 0.0001 boundary, '--' is NaN.
    """
    if pd.isna(x):
        return np.nan
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s in ("--", "", "nan", "NaN"):
        return np.nan
    if s.startswith("<"):
        try:
            return float(s.lstrip("<").strip())  # literal boundary, matches Patrick Code.txt parser
        except Exception:
            return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def _stat_column(singles: pd.DataFrame) -> str | None:
    """Find the Scaffold-emitted stat column (t-test or permutation)."""
    for c in singles.columns:
        if not isinstance(c, str):
            continue
        if "t-test" in c or "Permutation Test" in c:
            return c
    return None


def run_one(comparison_dir: Path, alpha: float = 0.05) -> None:
    wb = comparison_dir / "processed_protein_data.xlsx"
    singles = pd.read_excel(wb, sheet_name="Single Proteins")
    stat_col = _stat_column(singles)
    if stat_col is None:
        print(f"  [WARN] no Scaffold stat column in {wb}")
        return

    p_vals = singles[stat_col].apply(parse_scaffold_p)
    test_label = ("Scaffold permutation" if "Permutation" in stat_col
                  else "Scaffold t-test")
    singles["Test Column"] = p_vals

    p_arr = np.asarray(p_vals, dtype=float)
    n_tests = int(np.sum(~np.isnan(p_arr)))
    bonf_thresh = alpha / max(n_tests, 1)

    sig_uncorr = singles[p_arr < alpha].copy()
    sig_corr = singles[p_arr < bonf_thresh].copy()

    with pd.ExcelWriter(wb, engine="openpyxl", mode="a",
                        if_sheet_exists="replace") as xl:
        singles.to_excel(xl, sheet_name="Single Proteins", index=False)
        sig_uncorr.to_excel(xl, sheet_name="Sig without correction",
                            index=False)
        sig_corr.to_excel(xl, sheet_name="Sig with correction", index=False)

    summary = pd.DataFrame({
        "comparison": [comparison_dir.name],
        "test": [test_label],
        "n_singles": [len(singles)],
        "n_tested": [n_tests],
        "alpha": [alpha],
        "bonferroni_threshold": [bonf_thresh],
        "n_sig_uncorrected": [int(np.sum(p_arr < alpha))],
        "n_sig_bonferroni": [int(np.sum(p_arr < bonf_thresh))],
    })
    summary.to_csv(comparison_dir / "statistics_summary.csv", index=False)
    print(f"  {summary.to_string(index=False)}")


def main(argv):
    comparisons = argv[1:] or CANONICAL_COMPARISONS
    for comp in comparisons:
        d = VOLCANO_ROOT / comp
        wb = d / "processed_protein_data.xlsx"
        if not wb.exists():
            print(f"[SKIP] no processed workbook in {d}")
            continue
        print(f"\n=== {comp} ===")
        run_one(d)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
