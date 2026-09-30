"""
_subclass_rules.py
==================

Import shim. `06b_subclassify_proteins.py` starts with a digit, so it
cannot be imported with a normal `import` statement. This module loads
it by path and re-exports its public names, so 06 and 10 can simply do:

    from _subclass_rules import SUBCLASS_ORDER, subclass_table

This replaces the ad-hoc importlib block that used to sit at the top of
10_mucin_by_species.py (and which referenced names that script 06 never
defined).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_PATH = Path(__file__).resolve().parent / "06b_subclassify_proteins.py"
_spec = importlib.util.spec_from_file_location("_subclass_rules_impl", _PATH)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

classify_contractile = _mod.classify_contractile
classify_cytoskeleton = _mod.classify_cytoskeleton
classify_glycan = _mod.classify_glycan
CLASSIFIERS = _mod.CLASSIFIERS
SUBCLASS_ORDER = _mod.SUBCLASS_ORDER
SPECIES_ORDER = _mod.SPECIES_ORDER
subclass_table = _mod.subclass_table

__all__ = [
    "classify_contractile", "classify_cytoskeleton", "classify_glycan",
    "CLASSIFIERS", "SUBCLASS_ORDER", "SPECIES_ORDER", "subclass_table",
]
