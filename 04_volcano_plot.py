"""
04_volcano_plot.py
===================

EnhancedVolcano-style volcano plots for the four pairwise rotifer
comparisons (Figure 3 in the manuscript, panels A-D) and the
uncorrected supplement (Sup S3).

Authoritative data source: Patrick's legacy
`processed_protein_data (NN).xlsx` files, which contain pre-computed
"Sig with correction" and "Sig without correction" sheets generated
from per-pair Scaffold exports with Benjamini-Hochberg FDR at
q=0.05. (The manuscript caption labels this "Bonferroni-corrected";
the published `Sig with correction` Test Column ranges of [0.00001,
0.0059] for Panel A and [0.00001, 0.006] for Panel B are
diagnostic of BH-FDR with q=0.05, not strict Bonferroni.)

When the legacy file is present, the volcano red dots are drawn from
its `Sig with correction` sheet directly, with `|log2FC| >= 2` as the
effect-size filter. When the legacy file is absent, we fall back to
processed_protein_data.xlsx + a p<0.001 threshold (script 02's output).

Color categories (EnhancedVolcano convention):
    grey  - "NS"
    green - "Log2 FC"  (|log2FC| above threshold, p above)
    blue  - "p-value"  (p below, |log2FC| below)
    red   - "p-value and log2 FC"

INPUTS  : Volcano (1)/Volcano/<comparison>/processed_protein_data (NN).xlsx
          (or processed_protein_data.xlsx as fallback)
OUTPUTS : <comparison>_volcano_corrected.{png,svg,pdf}
          <comparison>_volcano_uncorrected.{png,svg,pdf}
          + composite Figure_3_volcano_panels_corrected.{pdf,png,tif}
          + composite Supplemental_Figure_S3_volcano_panels_uncorrected.{pdf,png,tif}

The composite is assembled as a true vector PDF by merging the four
panel PDFs; the PNG/TIF versions are rasterised from the same panels for
journals that ask for them.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"axes.labelsize": 16, "xtick.labelsize": 14,
                     "ytick.labelsize": 14, "pdf.fonttype": 42,
                     "ps.fonttype": 42, "svg.fonttype": "none"})

try:
    from adjustText import adjust_text
    _HAVE_ADJUSTTEXT = True
except ImportError:
    _HAVE_ADJUSTTEXT = False


PROJECT_ROOT = Path(__file__).resolve().parent.parent
VOLCANO_ROOT = PROJECT_ROOT / "Volcano (1)" / "Volcano"

CANONICAL_COMPARISONS = [
    "Cupelopagis vorax vs Plationus patulus",
    "Cupelopagis vorax vs Sinantherina socialis",
    "Sinantherina socialis vs Plationus patulus",
    "Lacinularia flosculosa vs Sinantherina socialis",
]

# Comparisons presented as "Effect size only in corrected analysis"
# (manuscript Figure 3, Panel C). For these, the "p-value and log2 FC"
# label category is suppressed — proteins that would otherwise be
# labeled red are presented as green "Log2 FC" hits instead.
# Panel D is the within-lifestyle control — manuscript caption marks
# it as "log2FC" (effect size only). Panel C still uses BH-FDR
# significance but happens to have zero red hits.
EFFECT_SIZE_ONLY = {
    "Sinantherina socialis vs Plationus patulus",
    "Lacinularia flosculosa vs Sinantherina socialis",
}

# Panels that should display "| log2FC" instead of "| BH-FDR" in the
# title. Per user request, the corrected Panel D is now labeled
# "| BH-FDR" (the test we actually ran, even though no proteins pass).
# This set is empty intentionally — uncomment Panel D to restore the
# log2FC-suffix appearance.
LOG2FC_TITLE_PANELS = set()

# Published figure title format per panel — uses lifestyle-prefix
# naming with (0) group named first (matches the published Figure 3).
# Verified against the data: the (0) group lands on the left side of
# every panel, so the published "X vs Y" title puts X first because
# X is the (0) group. This is why Panel C reads "P. patulus vs
# S. socialis" even though the folder name has them reversed.
PANEL_TITLES = {
    "Cupelopagis vorax vs Plationus patulus":
        r"Sessile Solitary ($\it{C.\ vorax}$) vs Motile Solitary ($\it{P.\ patulus}$)",
    "Cupelopagis vorax vs Sinantherina socialis":
        r"Sessile Solitary ($\it{C.\ vorax}$) vs Sessile Colonial ($\it{S.\ socialis}$)",
    "Sinantherina socialis vs Plationus patulus":
        r"Motile Solitary ($\it{P.\ patulus}$) vs Sessile Colonial ($\it{S.\ socialis}$)",
    "Lacinularia flosculosa vs Sinantherina socialis":
        r"Sessile Colonial Control: $\it{L.\ flosculosa}$ vs $\it{S.\ socialis}$",
}

_GROUP_RE = re.compile(r"^\(([01])\)")

COLOR_MAP = {
    "NS": "#9E9E9E",
    "Log2 FC": "#2E7D32",
    "p-value": "#1565C0",
    "p-value and log2 FC": "#C62828",
}


def group_columns(columns):
    g = {0: [], 1: []}
    for c in columns:
        if isinstance(c, str):
            m = _GROUP_RE.match(c)
            if m:
                g[int(m.group(1))].append(c)
    return g[0], g[1]


def find_legacy_file(comparison_dir):
    """Return path to legacy processed_protein_data (NN).xlsx if present."""
    for f in sorted(os.listdir(comparison_dir)):
        if f.startswith("processed_protein_data (") and f.endswith(".xlsx"):
            return comparison_dir / f
    return None


def read_fwer(legacy_path):
    """Read Scaffold's Familywise Error Rate from the Raw Extract sheet.
    This is the multiple-testing-corrected p-cutoff Scaffold itself
    computed for the comparison; it's what Patrick used as the
    "Sig with correction" threshold (legacy sheet's Test Column max
    matches this value)."""
    try:
        raw = pd.read_excel(legacy_path, sheet_name="Raw Extract",
                             header=None)
    except Exception:
        return None
    for i in range(len(raw)):
        for j in range(raw.shape[1]):
            v = str(raw.iat[i, j])
            if "Familywise" in v or "familywise" in v:
                for k in range(j + 1, raw.shape[1]):
                    val = raw.iat[i, k]
                    if pd.notna(val) and str(val).strip():
                        try:
                            return float(val)
                        except Exception:
                            return None
    return None


def compute_volcano_from_legacy(comparison_dir, fc_thresh=2.0,
                                 use_correction=True):
    """Read the legacy Sig sheet and the full Single Proteins sheet,
    classify each protein for the volcano plot.
    
    `use_correction=True` -> uses 'Sig with correction' (BH-FDR sheet)
                            and labels red dots as "p-value and log2 FC".
    `use_correction=False` -> uses 'Sig without correction' (raw p<0.05).
    """
    legacy = find_legacy_file(comparison_dir)
    if legacy is None:
        return None, None
    sig_sheet = "Sig with correction" if use_correction else "Sig without correction"
    sig = pd.read_excel(legacy, sheet_name=sig_sheet)
    full = pd.read_excel(legacy, sheet_name="Single Proteins")
    return sig, full


def compute_volcano(singles, g0, g1, p_thresh, fc_thresh,
                     sig_accessions=None, effect_size_only=False):
    """Build the volcano DataFrame.
    
    If `sig_accessions` is provided (set of Accession Numbers from a
    pre-computed Sig sheet), red dots are those proteins ALSO meeting
    the |log2FC|>=fc_thresh cutoff. The legacy BH-FDR sheets are the
    source of truth for "significance"; we just add the effect-size
    filter on top to mark the red category.
    """
    a = singles[g0].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    b = singles[g1].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    log2fc = np.log2(b + 1) - np.log2(a + 1)
    # Prefer the explicit "Test Column" if present (written by script 02),
    # else parse the raw Scaffold "t-test" / "Permutation Test" column.
    if "Test Column" in singles.columns:
        p = pd.to_numeric(singles["Test Column"], errors="coerce")
    else:
        stat_col = next((c for c in singles.columns
                          if isinstance(c, str)
                          and ("t-test" in c or "Permutation Test" in c)), None)
        if stat_col is None:
            p = pd.Series([np.nan] * len(singles))
        else:
            def _parse(x):
                if pd.isna(x):
                    return np.nan
                if isinstance(x, (int, float)):
                    return float(x)
                s = str(x).strip()
                if s in ("--", "", "nan", "NaN"):
                    return np.nan
                if s.startswith("<"):
                    try:
                        return float(s.lstrip("<").strip())
                    except Exception:
                        return np.nan
                try:
                    return float(s)
                except Exception:
                    return np.nan
            p = singles[stat_col].apply(_parse)
    # Build short, gene-name-style labels from Scaffold protein names.
    # Strip:
    #   "PREDICTED:" / "LOW QUALITY PROTEIN:" prefixes,
    #   GO-term tail after first '|',
    #   isoform suffixes ("isoform X1", "isoform A", etc.),
    #   trailing descriptive clauses after first ','.
    def _gene_label(name):
        s = str(name).split("|")[0].strip()
        s = re.sub(r"^(PREDICTED:|LOW QUALITY PROTEIN:)\s*", "", s,
                    flags=re.IGNORECASE)
        s = re.sub(r"\s+isoform\s+\S+.*$", "", s, flags=re.IGNORECASE)
        s = s.split(",")[0].strip()
        return s
    # Labels are gene-model accession numbers (e.g. "laflo19114.t1",
    # "ppat5539.t1"), matching the published Figure 3.
    df = pd.DataFrame({
        "label": singles.get("Accession Number", pd.Series([""] * len(singles))).astype(str),
        "protein_name_short": (
            singles.get("Protein Name", pd.Series([""] * len(singles)))
            .astype(str).map(_gene_label)),
        "Accession Number": singles.get("Accession Number"),
        "log2FC": log2fc,
        "pval": p,
    }).dropna(subset=["pval"])
    df["-log10(pval)"] = -np.log10(df["pval"].clip(lower=1e-300))

    if sig_accessions is not None:
        is_sig = df["Accession Number"].isin(sig_accessions)
    else:
        is_sig = df["pval"] < p_thresh
    is_fc = df["log2FC"].abs() >= fc_thresh

    df["category"] = "NS"
    df.loc[is_fc & ~is_sig, "category"] = "Log2 FC"
    df.loc[~is_fc & is_sig, "category"] = "p-value"
    # "p-value and log2 FC" category restricted to the top-10-down + top-10-up
    # proteins per panel (the labeling convention the manuscript uses).
    # Any remaining sig+|log2FC| proteins beyond that cap fall back to
    # "Log2 FC" so they're still rendered as effect-size hits.
    if effect_size_only:
        # Manuscript caption: "Effect size only in corrected analysis"
        # for this panel — render all sig + |log2FC| hits as green
        # (Log2 FC), with no red "p-value and log2 FC" labels.
        df.loc[is_fc & is_sig, "category"] = "Log2 FC"
        # Still label the BH-FDR-significant + |log2FC| hits in green so
        # the most biologically interesting proteins get called out.
        df["_label_me"] = is_sig & is_fc
    else:
        # ALL sig + |log2FC| proteins are red; only the labeling is
        # capped at top-10-down + top-10-up (matches published Figure 3).
        df.loc[is_fc & is_sig, "category"] = "p-value and log2 FC"
        candidate = df[is_sig & is_fc].copy()
        top_down = (candidate[candidate["log2FC"] <= -fc_thresh]
                    .sort_values("log2FC").head(10))
        top_up = (candidate[candidate["log2FC"] >= fc_thresh]
                  .sort_values("log2FC", ascending=False).head(10))
        top_idx = set(top_down.index).union(top_up.index)
        df["_label_me"] = df.index.isin(top_idx)
    return df


def render(df, *, title, p_thresh, fc_thresh, out_png, top_n=5):
    # Figure is wider than tall so the lifestyle-prefix title
    # ("Sessile Solitary (C. vorax) vs Motile Solitary (P. patulus)
    #   | BH-FDR p<6.51e-03") fits comfortably without being cut off.
    fig, ax = plt.subplots(figsize=(13, 8), dpi=300)
    for cat, color in COLOR_MAP.items():
        sub = df[df["category"] == cat]
        ax.scatter(sub["log2FC"], sub["-log10(pval)"],
                   label=f"{cat} (n={len(sub)})",
                   alpha=0.7, c=color, s=18, edgecolors="none")
    ax.axvline(fc_thresh, ls="--", c="black", lw=0.6)
    ax.axvline(-fc_thresh, ls="--", c="black", lw=0.6)
    if p_thresh is not None:
        ax.axhline(-np.log10(p_thresh), ls="--", c="black", lw=0.6)
    ax.set_xlabel("Log2 fold change")
    ax.set_ylabel("-Log10 P")
    ax.set_title(title, fontsize=17, pad=12)

    # Label only proteins explicitly marked _label_me (top-10-per-side
    # of the sig + |log2FC| set). Fallback to top hits by p when none.
    if "_label_me" in df.columns and df["_label_me"].any():
        top = df[df["_label_me"]]
    else:
        # Effect-size-only panel or no sig hits — label top |log2FC|
        # hits as a courtesy (manuscript Panel D shows a few).
        fc = df[df["category"] == "Log2 FC"]
        if len(fc):
            top = pd.concat([
                fc.sort_values("log2FC", ascending=False).head(top_n),
                fc.sort_values("log2FC").head(top_n),
            ])
        else:
            top = df[df["category"] == "p-value"] \
                .sort_values("-log10(pval)", ascending=False).head(top_n)

    texts = []
    for _, r in top.iterrows():
        texts.append(ax.text(r["log2FC"], r["-log10(pval)"], str(r["label"]),
                              fontsize=12, color="black"))
    if _HAVE_ADJUSTTEXT and texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="0.4", lw=0.5))

    # Put legend outside the data area so labels don't overlap it.
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
              fontsize=14, frameon=True, borderaxespad=0.)
    # Ensure FWER cutoff line stays visible: extend y-axis above it
    # if the data alone wouldn't reach it (Panel C, D case).
    if p_thresh is not None and p_thresh > 0:
        ymax_data = ax.get_ylim()[1]
        cutoff_y = -np.log10(p_thresh)
        if cutoff_y > ymax_data * 0.9:
            ax.set_ylim(top=cutoff_y * 1.08)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_png.with_suffix(".svg"))
    fig.savefig(out_png.with_suffix(".pdf"))
    plt.close(fig)
    print(f"  wrote {out_png.name}")


def run_one(comparison_dir, fc_thresh=2.0):
    legacy = find_legacy_file(comparison_dir)
    if legacy is not None:
        fwer = read_fwer(legacy)
        for use_corr, label in [(True, "corrected"), (False, "uncorrected")]:
            sheet = "Sig with correction" if use_corr else "Sig without correction"
            sig = pd.read_excel(legacy, sheet_name=sheet)
            full = pd.read_excel(legacy, sheet_name="Single Proteins")
            sig_accs = set(sig["Accession Number"].dropna()) \
                if "Accession Number" in sig.columns else set()
            g0, g1 = group_columns(list(full.columns))
            if not g0 or not g1:
                print(f"[WARN] no group columns in legacy file")
                continue
            # No p-value override: keep Scaffold's literal "< 0.0001"
            # boundary parsing (so all "< 0.0001" proteins plot at
            # y=4.0, matching Patrick's original Figure 3).
            df_v = compute_volcano(
                full, g0, g1, p_thresh=None, fc_thresh=fc_thresh,
                sig_accessions=sig_accs,
                effect_size_only=(use_corr and
                                   comparison_dir.name in EFFECT_SIZE_ONLY),
            )
            out = comparison_dir / f"{comparison_dir.name}_volcano_{label}.png"
            # Horizontal dashed line: max(FWER, 0.0001) for the
            # corrected panel — Scaffold's "< 0.0001" notation puts a
            # floor at 0.0001, so when FWER drops below that we use
            # the visible floor instead (matches Patrick's original).
            # Uncorrected panels use raw p<0.05.
            if use_corr and fwer is not None:
                line_p = max(fwer, 0.0001)
            elif use_corr:
                line_p = 0.0001
            else:
                line_p = 0.05
            base = PANEL_TITLES.get(comparison_dir.name, comparison_dir.name)
            if use_corr and comparison_dir.name in LOG2FC_TITLE_PANELS:
                title = f"{base} | log2FC"
            elif use_corr and fwer is not None:
                title = f"{base} | BH-FDR p<{fwer:.2e}"
            elif use_corr:
                title = f"{base} | BH-FDR"
            else:
                title = f"{base} | uncorrected"
            render(df_v, title=title, p_thresh=line_p,
                   fc_thresh=fc_thresh, out_png=out)
    else:
        # Fallback: use processed_protein_data.xlsx + p<0.001 + |log2FC|>=2
        wb = comparison_dir / "processed_protein_data.xlsx"
        if not wb.exists():
            print(f"[SKIP] no data in {comparison_dir}")
            return
        singles = pd.read_excel(wb, sheet_name="Single Proteins")
        g0, g1 = group_columns(list(singles.columns))
        if not g0 or not g1:
            print(f"[WARN] no group columns in {wb}")
            return
        for p_thresh, label in [(0.05, "uncorrected"), (0.001, "corrected")]:
            df_v = compute_volcano(singles, g0, g1, p_thresh, fc_thresh)
            out = comparison_dir / f"{comparison_dir.name}_volcano_{label}.png"
            render(df_v, title=f"{comparison_dir.name} (p<{p_thresh})",
                   p_thresh=p_thresh, fc_thresh=fc_thresh, out_png=out)


def main(argv):
    comparisons = argv[1:] or CANONICAL_COMPARISONS
    for comp in comparisons:
        d = VOLCANO_ROOT / comp
        if not d.is_dir():
            print(f"[SKIP] not a directory: {d}")
            continue
        print(f"\n=== {comp} ===")
        run_one(d)
    print("\n=== composites ===")
    for label in ("corrected", "uncorrected"):
        assemble_composite(label)
    return 0


# ----------------------------------------------------------------------
# 4-panel composite assembly (Figure 3 + Sup S3)
# ----------------------------------------------------------------------
PANELS = [
    ("A", "Cupelopagis vorax vs Plationus patulus"),
    ("B", "Cupelopagis vorax vs Sinantherina socialis"),
    ("C", "Sinantherina socialis vs Plationus patulus"),
    ("D", "Lacinularia flosculosa vs Sinantherina socialis"),
]

STEM = {
    "corrected": "Figure_3_volcano_panels_corrected",
    "uncorrected": "Supplemental_Figure_S3_volcano_panels_uncorrected",
}


def assemble_composite(label):
    """Merge the four panel PDFs into one vector 2x2 page, then
    rasterise PNG/TIF from the same panels.

    Vector output matters here: Reviewer 2 asked for a figure in which
    the accession labels on the volcano points can be read by zooming.
    """
    import io
    from pypdf import PdfReader, PdfWriter, PageObject, Transformation
    from reportlab.pdfgen import canvas as rlcanvas

    pdfs = [VOLCANO_ROOT / c / f"{c}_volcano_{label}.pdf" for _, c in PANELS]
    if not all(p.exists() for p in pdfs):
        print(f"  [SKIP] missing panel PDFs for {label}")
        return
    pages = [PdfReader(str(p)).pages[0] for p in pdfs]
    w = max(float(p.mediabox.width) for p in pages)
    h = max(float(p.mediabox.height) for p in pages)
    gutter = 26                      # white space above each panel for its letter
    W, H = w * 2, (h + gutter) * 2

    out = PageObject.create_blank_page(width=W, height=H)
    grid = [(0, 1), (1, 1), (0, 0), (1, 0)]      # A B / C D
    for page, (col, row) in zip(pages, grid):
        out.merge_transformed_page(
            page, Transformation().translate(tx=col * w, ty=row * (h + gutter)))

    buf = io.BytesIO()
    cv = rlcanvas.Canvas(buf, pagesize=(W, H))
    cv.setFont("Helvetica-Bold", 26)
    for (letter, _), (col, row) in zip(PANELS, grid):
        cv.drawString(col * w + 14, row * (h + gutter) + h + 2, letter)
    cv.save()
    buf.seek(0)
    out.merge_page(PdfReader(buf).pages[0])

    stem = PROJECT_ROOT / STEM[label]
    writer = PdfWriter()
    writer.add_page(out)
    with open(stem.with_suffix(".pdf"), "wb") as fh:
        writer.write(fh)
    print(f"  wrote {stem.name}.pdf (vector, {W:.0f}x{H:.0f} pt)")

    # Raster companions, rendered from the vector composite when
    # pdf2image/poppler is available; otherwise skipped (the PDF is the
    # figure of record).
    try:
        from pdf2image import convert_from_path
        img = convert_from_path(str(stem.with_suffix(".pdf")), dpi=300)[0]
        img.save(str(stem) + ".png", dpi=(300, 300))
        img.save(str(stem) + ".tif", dpi=(300, 300), compression="tiff_lzw")
        print(f"  wrote {stem.name}.png/.tif")
    except Exception as exc:
        print(f"  [note] raster companions skipped ({exc.__class__.__name__});"
              f" the vector PDF is the figure of record")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
