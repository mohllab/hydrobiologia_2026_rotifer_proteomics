"""
common.py
=========

Shared helpers used by scripts 06–10:
  - load the project's sample metadata (species, location, lifestyle)
  - load the full Single Proteins table for one comparison folder
    (any of the four folders carries the same 1675-protein matrix
    across all 20 samples, only the (0)/(1) grouping differs)
  - infer per-species protein detection (count of proteins with
    spectral count > 0 in at least one replicate of that species)
  - canonical species order for figures (sessile -> motile, matching
    the manuscript Sup S4-S6 ordering)
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CODE_ROOT = Path(__file__).resolve().parent
VOLCANO_ROOT = PROJECT_ROOT / "Volcano (1)" / "Volcano"
META_CSV = CODE_ROOT / "samples_metadata.csv"

# Pick any folder that produced the full 20-sample export
DEFAULT_FOLDER = VOLCANO_ROOT / "Cupelopagis vorax vs Plationus patulus"

# Manuscript figure ordering: sessile (left) then motile (right),
# matching Sup S4/S5/S6.
SPECIES_ORDER = [
    "C. vorax",        # sessile solitary
    "S. socialis",     # sessile colonial
    "L. flosculosa",   # sessile colonial
    "N. copeus",       # motile solitary
    "C. hippocrepis",  # motile colonial
    "P. patulus",      # motile solitary
    "E. kingi",        # motile solitary
    "Hexarthra sp.",   # motile solitary
]


def load_metadata() -> pd.DataFrame:
    return pd.read_csv(META_CSV)


def load_full_singles(folder: Path = DEFAULT_FOLDER) -> pd.DataFrame:
    """Read the Single Proteins sheet (all 1675 proteins x 20 samples)."""
    wb = folder / "processed_protein_data.xlsx"
    if not wb.exists():
        raise FileNotFoundError(
            f"{wb} missing - run 01_data_prep.py first"
        )
    return pd.read_excel(wb, sheet_name="Single Proteins")


_GROUP_RE = re.compile(r"^\(([01])\)\s*")


def strip_group_prefix(c: str) -> str:
    if not isinstance(c, str):
        return c
    return _GROUP_RE.sub("", c)


def sample_columns(df: pd.DataFrame, meta: pd.DataFrame) -> list[str]:
    """Return the columns in df that correspond to known sample names."""
    keep = []
    known = set(meta["sample_name"])
    for c in df.columns:
        clean = strip_group_prefix(c)
        if clean in known:
            keep.append(c)
    return keep


def species_detection_matrix(
    df: pd.DataFrame, meta: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build two species-level matrices from the protein x sample table:

      - detection_per_species : (n_proteins x n_species) boolean
        True where the protein has spectral count > 0 in at least one
        replicate of that species.
      - abundance_per_species : (n_proteins x n_species) float
        Mean spectral count across replicates of the species
        (NaN-safe). Useful for PCA at species-mean level.

    The original sample-level matrix is also returned via
    abundance_per_sample (n_proteins x n_samples).
    """
    cols = sample_columns(df, meta)
    abund = (
        df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).copy()
    )
    abund.columns = [strip_group_prefix(c) for c in cols]

    name_to_species = dict(zip(meta["sample_name"], meta["species_short"]))
    species_to_cols: dict[str, list[str]] = {}
    for c in abund.columns:
        sp = name_to_species.get(c)
        if sp:
            species_to_cols.setdefault(sp, []).append(c)

    det = pd.DataFrame(index=df.index, columns=SPECIES_ORDER, dtype=bool)
    mean = pd.DataFrame(index=df.index, columns=SPECIES_ORDER, dtype=float)
    for sp in SPECIES_ORDER:
        cs = species_to_cols.get(sp, [])
        if not cs:
            det[sp] = False
            mean[sp] = 0.0
            continue
        block = abund[cs]
        det[sp] = (block > 0).any(axis=1)
        mean[sp] = block.mean(axis=1)
    return det, mean, abund


def lifestyle_label(meta: pd.DataFrame, species_short: str,
                    style: str = "4") -> str:
    """Look up the lifestyle label string for a species."""
    sub = meta[meta["species_short"] == species_short]
    if sub.empty:
        return ""
    col = "lifestyle_4" if style == "4" else "lifestyle_2"
    return sub[col].iloc[0]


def italic_species(name: str) -> str:
    """Return a matplotlib mathtext label with the species name in italics.

    Taxonomic convention: the Latin binomial is italicised, but the
    abbreviation "sp." is not. "Hexarthra sp." therefore renders as an
    italic genus followed by an upright "sp.". Works for full names
    ("Cupelopagis vorax") and abbreviated ones ("C. vorax").
    """
    name = str(name).strip()
    if name.endswith(" sp."):
        genus = name[:-4]
        return "$\\it{" + genus.replace(" ", "\\ ") + "}$ sp."
    return "$\\it{" + name.replace(" ", "\\ ") + "}$"
