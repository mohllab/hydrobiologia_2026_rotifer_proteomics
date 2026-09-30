"""
06b_subclassify_proteins.py
============================

Keyword sub-classification of proteins WITHIN the three focal classes
assigned by 09_top10_functional.py:

  - Contractile / muscle-associated
  - Cytoskeleton / motility
  - Mucin / glycan-associated / surface

Each cleaned protein name is tested against an ordered list of keyword
rules and assigned to the FIRST rule it matches. Because assignment is
order-dependent, the more specific term must always be tested before any
term it contains as a substring. Two such pairs exist here:

    "paramyosin"  contains "myosin"   -> paramyosin is tested first
    "tropomyosin" contains "myosin"   -> tropomyosin is tested first

Both orderings were inverted in the original version of this script,
which made the Paramyosin and Tropomyosin branches unreachable. They are
correct here, and 06a_validate_subclass_counts.py asserts that the
result reproduces the published per-species totals.

Input:
  - rotifer_species_unique_protein_presence_with_mucin_contractile.csv
    (written by 09_top10_functional.py)

Output:
  - prints subclass counts per species for each focal class
  - importable: CLASSIFIERS, SUBCLASS_ORDER, subclass_table()
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_FILE = (PROJECT_ROOT /
              "rotifer_species_unique_protein_presence_with_mucin_contractile.csv")

SPECIES_ORDER = [
    'Cupelopagis vorax', 'Sinantherina socialis', 'Lacinularia flosculosa',
    'Notommata copeus', 'Conochilus hippocrepis', 'Plationus patulus',
    'Euchlanis kingi', 'Hexarthra sp.',
]


def classify_contractile(name: str) -> str:
    """Sub-classify contractile / muscle-associated proteins.

    'paramyosin' is tested before 'myosin' because it contains it.
    """
    n = str(name).lower()
    if 'troponin' in n:
        return 'Troponin'
    if 'paramyosin' in n:                 # must precede 'myosin'
        return 'Paramyosin'
    if 'myosin' in n:
        return 'Myosin'
    if 'twitchin' in n:
        return 'Twitchin'
    if 'titin' in n:
        return 'Titin'
    if 'muscle lim' in n:
        return 'Muscle LIM'
    if 'tektin' in n:
        return 'Tektin'
    if 'rootletin' in n:
        return 'Rootletin'
    if 'sarcoplasmic' in n:
        return 'Sarcoplasmic Ca-binding'
    return 'Other'


def classify_cytoskeleton(name: str) -> str:
    """Sub-classify cytoskeleton / motility proteins.

    'tropomyosin' is tested first because it contains 'myosin' and,
    in annotations that carry a GO tail, also 'actin'.

    Note on myosins: 09_top10_functional.py assigns myosins whose
    annotation specifies striated or sarcomeric muscle to the
    contractile class, so anything reaching this function is a
    non-muscle or unspecified-type myosin.
    """
    n = str(name).lower()
    if 'tropomyosin' in n:                # must precede 'actin'/'myosin'
        return 'Tropomyosin'
    if 'actin' in n and 'alpha-actin' not in n:
        return 'Actin'
    if 'tubulin' in n:
        return 'Tubulin'
    if 'intermediate filament' in n or 'ifa-' in n:
        return 'Intermediate filament'
    if 'myosin' in n:
        return 'Myosin (non-muscle/unspec.)'
    if 'dynein' in n:
        return 'Dynein'
    if 'kinesin' in n:
        return 'Kinesin'
    if 'spectrin' in n or 'alpha-actin' in n:
        return 'Spectrin/Alpha-actinin'
    if 'filamin' in n:
        return 'Filamin'
    if 'cofilin' in n or 'profilin' in n or 'gelsolin' in n:
        return 'Actin-binding/regulatory'
    return 'Other'


def classify_glycan(name: str) -> str:
    """Sub-classify mucin / glycan-associated / surface proteins."""
    n = str(name).lower()
    if 'mucin' in n:
        return 'Mucin'
    if 'galectin' in n:
        return 'Galectin'
    if 'proteoglycan' in n or 'heparan' in n or 'chondroitin' in n:
        return 'Proteoglycan'
    if 'von willebrand' in n or 'vwd' in n:
        return 'VWD domain protein'
    if 'dolichyl' in n or 'glycosyltransferase' in n:
        return 'Glycosyltransferase'
    if 'ependymin' in n:
        return 'Ependymin'
    if 'lectin' in n:
        return 'Lectin'
    return 'Other glycan'


CLASSIFIERS = {
    'Mucin / glycan-associated / surface': classify_glycan,
    'Contractile / muscle-associated':      classify_contractile,
    'Cytoskeleton / motility':              classify_cytoskeleton,
}

# Plot/stack order for each focal class (used by 06_subclass_panels.py)
SUBCLASS_ORDER = {
    'Mucin / glycan-associated / surface': [
        'Mucin', 'Galectin', 'Proteoglycan', 'VWD domain protein',
        'Glycosyltransferase', 'Ependymin', 'Lectin', 'Other glycan'],
    'Contractile / muscle-associated': [
        'Troponin', 'Myosin', 'Twitchin', 'Titin', 'Muscle LIM',
        'Paramyosin', 'Tektin', 'Rootletin', 'Sarcoplasmic Ca-binding',
        'Other'],
    'Cytoskeleton / motility': [
        'Actin', 'Tubulin', 'Intermediate filament',
        'Myosin (non-muscle/unspec.)', 'Tropomyosin',
        'Spectrin/Alpha-actinin', 'Filamin', 'Dynein', 'Kinesin',
        'Actin-binding/regulatory', 'Other'],
}


def subclass_table(csv_path: Path = INPUT_FILE) -> dict:
    """Return {focal class: {species: {subclass: count}}}."""
    df = pd.read_csv(csv_path)
    out = {c: defaultdict(lambda: defaultdict(int)) for c in CLASSIFIERS}
    for _, row in df.iterrows():
        cat = row['Category']
        fn = CLASSIFIERS.get(cat)
        if fn is None:
            continue
        out[cat][row['Species']][fn(row['Clean Protein Name'])] += 1
    return {c: {sp: dict(v) for sp, v in d.items()} for c, d in out.items()}


def main() -> int:
    if not INPUT_FILE.exists():
        print(f"[ERROR] {INPUT_FILE.name} not found — run "
              f"09_top10_functional.py first.")
        return 1
    table = subclass_table()
    for cat in CLASSIFIERS:
        print("\n" + "=" * 62)
        print(f"  {cat.upper()}")
        print("=" * 62)
        for sp in SPECIES_ORDER:
            d = table[cat].get(sp, {})
            print(f"\n  {sp}: {sum(d.values())} total")
            for sub in SUBCLASS_ORDER[cat]:
                if d.get(sub):
                    print(f"    {sub}: {d[sub]}")
    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
