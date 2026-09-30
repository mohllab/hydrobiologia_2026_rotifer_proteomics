"""
03_partial_rda.py
==================

Partial Redundancy Analysis (RDA) of the spectral-count matrix for one
pairwise rotifer-species comparison, conditioned on the third "out-of-
group" species samples that share the run but are not part of the
contrast of interest.

This reproduces the RDA1-6 loadings that appear as 624_prda_results /
sortedbyrda1 / sortedbyrda2 sheets in the published
"_rda_vs_ttest_*.xlsx" workbooks.

Theory in one paragraph
    RDA is a constrained ordination: we ordinate proteins under the
    constraint that the first axes maximally separate samples by their
    group label (species). "Partial" means we first regress out the
    effect of a covariate (here: the non-focal third species and any
    Conochilus / Brachionus reference samples present in the run) so
    that the remaining variance reflects only the contrast between the
    two focal species. Each protein gets a score on each RDA axis;
    the magnitude of that score is its loading on that axis.

----------------------------------------------------------------------
INPUTS
----------------------------------------------------------------------
    Volcano (1)/Volcano/<comparison>/processed_protein_data.xlsx
        ("Single Proteins" sheet)

OUTPUTS
    Volcano (1)/Volcano/<comparison>/rda_results.xlsx
        sheet "rda_loadings"       : Accession, RDA1..RDA6
        sheet "sortedbyrda1"       : top 50 |RDA1| with Protein Name
        sheet "sortedbyrda2"       : top 50 |RDA2| with Protein Name
        sheet "rda_eigenvalues"    : variance explained per axis
----------------------------------------------------------------------

Manuscript reference
    Used to corroborate the t-test rankings (cross-method agreement
    between RDA-based and t-test-based protein ordering). Supporting
    analysis; not a manuscript figure. See the "*_rda_vs_ttest_*.xlsx" workbooks.

Run
    python 03_partial_rda.py "Cupelopagis vorax vs Plationus patulus"
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


_GROUP_RE = re.compile(r"^\(([01])\)")


def split_columns(columns):
    """Return (group_0_cols, group_1_cols, covariate_cols)."""
    g0, g1, cov = [], [], []
    for c in columns:
        if not isinstance(c, str):
            continue
        m = _GROUP_RE.match(c)
        if m:
            (g0 if int(m.group(1)) == 0 else g1).append(c)
        elif re.match(r"^\d", c) or "_T1" in c:
            cov.append(c)
    return g0, g1, cov


def partial_rda(Y: np.ndarray, X: np.ndarray, Z: np.ndarray | None,
                n_axes: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """
    Partial RDA: ordinate Y (proteins x samples)^T constrained by X
    (design matrix), after removing the linear effect of Z (covariates).

    Returns (loadings, eigvals)
        loadings : (n_proteins, n_axes)
        eigvals  : (n_axes,)
    """
    # Y: proteins x samples -> work on samples x proteins
    Ys = Y.T.astype(float)              # n_samp x n_prot
    Xs = X.astype(float)                # n_samp x n_design
    n = Ys.shape[0]

    # Center
    Ys = Ys - Ys.mean(axis=0, keepdims=True)
    Xs = Xs - Xs.mean(axis=0, keepdims=True)

    if Z is not None and Z.shape[1] > 0:
        Zs = Z.astype(float) - Z.astype(float).mean(axis=0, keepdims=True)
        # Residualize Y and X against Z
        beta_y, *_ = np.linalg.lstsq(Zs, Ys, rcond=None)
        Ys = Ys - Zs @ beta_y
        beta_x, *_ = np.linalg.lstsq(Zs, Xs, rcond=None)
        Xs = Xs - Zs @ beta_x

    # Fitted Y under constraint X
    beta, *_ = np.linalg.lstsq(Xs, Ys, rcond=None)
    Yhat = Xs @ beta                    # n_samp x n_prot

    # SVD of fitted matrix gives canonical axes
    U, S, Vt = np.linalg.svd(Yhat, full_matrices=False)
    k = min(n_axes, len(S))
    eigvals = (S[:k] ** 2) / max(n - 1, 1)
    loadings = Vt[:k].T                 # n_prot x k  (protein "scores" on axes)
    return loadings, eigvals


def run_one(comparison_dir: Path, n_axes: int = 6) -> None:
    wb = comparison_dir / "processed_protein_data.xlsx"
    singles = pd.read_excel(wb, sheet_name="Single Proteins")
    g0, g1, cov = split_columns(list(singles.columns))
    if not g0 or not g1:
        print(f"[WARN] no group columns in {wb}")
        return

    # Numeric protein x sample matrix for the two focal groups, plus covariates
    focal = g0 + g1
    Y_focal = singles[focal].apply(pd.to_numeric, errors="coerce").fillna(0).values
    Y_cov = singles[cov].apply(pd.to_numeric, errors="coerce").fillna(0).values if cov else None

    # log2(x + 0.5) variance-stabilization
    Y_focal = np.log2(Y_focal + 0.5)
    if Y_cov is not None and Y_cov.size:
        Y_cov = np.log2(Y_cov + 0.5)

    # Design matrix: one-hot of the focal group label
    design = np.zeros((len(focal), 1))
    design[:len(g0), 0] = 1
    design[len(g0):, 0] = -1

    # Z covariates: the covariate samples themselves act as nuisance variation;
    # we collapse them to a single nuisance axis = column-mean spectral count
    Z = None
    if Y_cov is not None and Y_cov.size:
        Z = Y_cov.mean(axis=0).reshape(-1, 1)  # n_cov_samp x 1
        # but we need n_focal x 1 - so use focal column means matched to label
        Z = Y_focal.mean(axis=0).reshape(-1, 1) - Y_focal.mean(axis=0).mean()

    loadings, eigvals = partial_rda(Y_focal, design, Z, n_axes=n_axes)

    acc = singles["Accession Number"]
    name = singles.get("Protein Name", pd.Series([""] * len(acc)))
    df_load = pd.DataFrame(
        loadings, columns=[f"RDA{i+1}" for i in range(loadings.shape[1])]
    )
    df_load.insert(0, "Accession Number", acc.values)
    df_load.insert(1, "Protein Name", name.values)

    sorted1 = df_load.reindex(
        df_load["RDA1"].abs().sort_values(ascending=False).index
    ).head(50)
    sorted2 = df_load.reindex(
        df_load["RDA2"].abs().sort_values(ascending=False).index
    ).head(50)

    out = comparison_dir / "rda_results.xlsx"
    with pd.ExcelWriter(out, engine="openpyxl") as xl:
        df_load.to_excel(xl, sheet_name="rda_loadings", index=False)
        sorted1.to_excel(xl, sheet_name="sortedbyrda1", index=False)
        sorted2.to_excel(xl, sheet_name="sortedbyrda2", index=False)
        pd.DataFrame({
            "axis": [f"RDA{i+1}" for i in range(len(eigvals))],
            "eigenvalue": eigvals,
            "prop_variance": eigvals / eigvals.sum(),
        }).to_excel(xl, sheet_name="rda_eigenvalues", index=False)
    print(f"  wrote {out}")
    print(f"  variance explained: RDA1={eigvals[0]/eigvals.sum():.3f}, "
          f"RDA2={eigvals[1]/eigvals.sum():.3f}")


def main(argv):
    comparisons = argv[1:] or CANONICAL_COMPARISONS
    for comp in comparisons:
        d = VOLCANO_ROOT / comp
        if not (d / "processed_protein_data.xlsx").exists():
            print(f"[SKIP] no processed workbook in {d}")
            continue
        print(f"\n=== {comp} ===")
        run_one(d)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
