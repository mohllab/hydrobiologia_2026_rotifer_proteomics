"""
01_data_prep.py
================

Prepare a Scaffold DDA spectral-count export for downstream comparative
analysis between two rotifer species.

Pipeline reproduced:
    Scaffold "Samples Report" (CSV/XLSX) --> processed_protein_data.xlsx
    with the same sheet layout used throughout the manuscript:

      - Raw Extract              : headerless dump of the protein table region
      - All Protein Data         : header-aware table, all rows
      - Single Proteins          : rows that carry a Molecular Weight value
                                   (single, non-clustered protein entries)
      - Protein Groups           : cluster header rows (no Molecular Weight)

Sample columns are interpreted by the "(0)" / "(1)" group-assignment prefix
that Scaffold applies when two species are flagged in the categories tab.

----------------------------------------------------------------------
INPUTS  (relative to the project root)
----------------------------------------------------------------------
    Volcano (1)/Volcano/<comparison>/Samples Report of My Experiment.csv

OUTPUTS
    Volcano (1)/Volcano/<comparison>/processed_protein_data.xlsx
----------------------------------------------------------------------

Manuscript reference
    "Following protein identification in Proteome Discoverer and Scaffold,
    we identified 1675 proteins..."
    Used for: Tables of differentially abundant proteins; input to
    02_stats_and_corrections.py and 04_volcano_plot.py.

Run
    python 01_data_prep.py "Cupelopagis vorax vs Plationus patulus"
    # or, no arg = run all four canonical comparisons
"""

from __future__ import annotations

import csv
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

HEADER_KEYS = ("Protein Name", "Accession Number", "Molecular Weight")
GROUP_PREFIX_RE = re.compile(r"^\(([01])\)")


def _read_raw(samples_report: Path) -> pd.DataFrame:
    """Read the Scaffold samples report as a uniform-width DataFrame.

    The Scaffold export has a long preamble with rows of varying column
    count, which makes pandas' CSV sniffer infer a 1-column table.
    We parse it with the standard csv module and pad rows manually.
    """
    if samples_report.suffix.lower() == ".csv":
        with open(samples_report, "r", encoding="utf-8",
                  errors="ignore", newline="") as fh:
            rows = list(csv.reader(fh))
        width = max(len(r) for r in rows) if rows else 0
        rows = [r + [""] * (width - len(r)) for r in rows]
        return pd.DataFrame(rows)
    return pd.read_excel(samples_report, header=None)


def _find_header_row(raw: pd.DataFrame) -> int:
    for i in range(min(len(raw), 200)):
        row = raw.iloc[i].astype(str).tolist()
        if all(any(k in cell for cell in row) for k in HEADER_KEYS):
            return i
    raise RuntimeError("Could not locate the Scaffold protein-table header row.")


def parse_samples_report(samples_report: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = _read_raw(samples_report)
    hdr = _find_header_row(raw)
    table = raw.iloc[hdr:].reset_index(drop=True).copy()
    table.columns = table.iloc[0].astype(str).str.strip().tolist()
    table = table.iloc[1:].reset_index(drop=True)
    # drop empty trailing rows
    keep = table.apply(lambda r: any(str(v).strip() for v in r), axis=1)
    table = table[keep].reset_index(drop=True)
    # drop the END OF FILE marker row if present
    table = table[~table.iloc[:, 0].astype(str).str.startswith("END OF FILE")]
    return raw, table.reset_index(drop=True)


def discover_group_columns(columns: list[str]) -> dict[int, list[str]]:
    groups: dict[int, list[str]] = {0: [], 1: []}
    for c in columns:
        if not isinstance(c, str):
            continue
        m = GROUP_PREFIX_RE.match(c)
        if m:
            groups[int(m.group(1))].append(c)
    return groups


def coerce_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        out[c] = pd.to_numeric(
            out[c].replace({"--": np.nan, "": np.nan}),
            errors="coerce",
        )
    return out


def split_singles_and_groups(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mw = pd.to_numeric(
        table.get("Molecular Weight").astype(str).str.replace(" kDa", "", regex=False),
        errors="coerce",
    )
    singles = table[mw.notna()].reset_index(drop=True)
    groups = table[mw.isna()].reset_index(drop=True)
    return singles, groups


def _pick_samples_report(comparison_dir: Path) -> Path:
    """Pick the samples-report file that actually carries the (0)/(1)
    group-assignment prefixes. Some exports lost those prefixes when
    saved to CSV, but kept them in the original XLSX."""
    csv = comparison_dir / "Samples Report of My Experiment.csv"
    xlsx = comparison_dir / "Samples Report of My Experiment.xlsx"

    def _has_group_prefix(path: Path) -> bool:
        try:
            _, table = parse_samples_report(path)
            return any(
                isinstance(c, str) and GROUP_PREFIX_RE.match(c)
                for c in table.columns
            )
        except Exception:
            return False

    if xlsx.exists() and _has_group_prefix(xlsx):
        return xlsx
    if csv.exists() and _has_group_prefix(csv):
        return csv
    # Fall back to whichever exists; group-discovery will warn.
    if csv.exists():
        return csv
    if xlsx.exists():
        return xlsx
    raise FileNotFoundError(f"No samples report in {comparison_dir}")


def process_one(comparison_dir: Path) -> Path:
    samples_csv = _pick_samples_report(comparison_dir)
    print(f"  source: {samples_csv.name}")

    raw, all_data = parse_samples_report(samples_csv)
    groups_map = discover_group_columns(list(all_data.columns))
    sample_cols = groups_map[0] + groups_map[1]
    all_data = coerce_numeric(all_data, sample_cols)
    singles, groups = split_singles_and_groups(all_data)

    out_path = comparison_dir / "processed_protein_data.xlsx"
    with pd.ExcelWriter(out_path, engine="openpyxl") as xl:
        raw.to_excel(xl, sheet_name="Raw Extract", header=False, index=False)
        all_data.to_excel(xl, sheet_name="All Protein Data", index=False)
        singles.to_excel(xl, sheet_name="Single Proteins", index=False)
        groups.to_excel(xl, sheet_name="Protein Groups", index=False)

    print(f"  group 0 ({len(groups_map[0])} samples)")
    print(f"  group 1 ({len(groups_map[1])} samples)")
    print(f"  single proteins: {len(singles)}, protein groups: {len(groups)}")
    print(f"  wrote {out_path}")
    return out_path


def main(argv: list[str]) -> int:
    comparisons = argv[1:] or CANONICAL_COMPARISONS
    for comp in comparisons:
        d = VOLCANO_ROOT / comp
        if not d.is_dir():
            print(f"[SKIP] not a directory: {d}")
            continue
        print(f"\n=== {comp} ===")
        process_one(d)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
