"""
05_consensus_annotation.py
===========================

Build a bitscore-weighted consensus protein annotation for every query
accession that appears in one or more of the per-comparison
"corrected" / "uncorrected" Bonferroni tables.

Pipeline reproduced:
    BLAST/DIAMOND hits file  +  per-comparison hit-list workbook
        --> consensus_per_query  (one row per query)
        --> unique_accessions    (one row per unique query, joined with
                                  source-table metadata)
        --> occurrences          (one row per query x comparison occurrence)

The consensus rule: for each query, group hits by their normalized
protein-name token set (everything after "RecName:" trimmed of organism
suffixes, lower-cased, stop-words removed). Each group's weight is the
sum of bitscores of its member hits. The annotation of the heaviest
group is the consensus. We also export:
    confidence_score_0to1     : (top group weight) / (total weight)
    support_fraction_bitscore : top group weight / max possible
    agreement_fraction_hits   : (top group hit count) / (total hits)

----------------------------------------------------------------------
INPUTS
----------------------------------------------------------------------
    Corrected_Bonferroni_completed_table_with_manual_uniprot_accessions_031226 (1).xlsx
        sheets used: "Corrected", "Uncorrected"
    (Optional) external BLAST tabular file, default name:
    "blast_hits_per_query.tsv"
        columns: qseqid, sseqid, salltitles, evalue, bitscore, length,
                  pident, qstart, qend, sstart, send, qlen, slen,
                  mismatch, gapopen
    (Optional) "subject_title_cache.tsv" with columns:
        sseqid, title, organism, length_aa

OUTPUTS
    All_Accession_Numbers_compiled_with_consensus.xlsx
      sheets:
        unique_accessions, occurrences, consensus_per_query,
        subject_title_cache, hits_used_for_consensus, missing_subjects
----------------------------------------------------------------------

Manuscript reference
    Methods, "Protein Characterization": query sequences were aligned
    against the NCBI non-redundant (nr) protein database (Galaxy index
    dated 03 Sep 2023) using DIAMOND in BLASTp mode (Galaxy v2.1.22),
    and a bitscore-weighted consensus annotation was derived per query.
    Source of the consensus_annotation column in the supplemental
    tables. (Genome annotation itself was done separately in OmicsBox
    with BLASTP and InterProScan; that step is not reproduced here.)

Run
    python 05_consensus_annotation.py
"""

from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent

CORRECTED_TBL = (PROJECT_ROOT /
                 "Corrected_Bonferroni_completed_table_with_manual_uniprot_accessions_031226 (1).xlsx")
HITS_TSV = PROJECT_ROOT / "blast_hits_per_query.tsv"
SUBJ_TSV = PROJECT_ROOT / "subject_title_cache.tsv"
OUT_PATH = PROJECT_ROOT / "All_Accession_Numbers_compiled_with_consensus.xlsx"


# ----------------------------------------------------------------------
# Name normalization
# ----------------------------------------------------------------------
_STOP = set("the of and a an in by from with for to ".split())
_ORG_SPLIT_RE = re.compile(r"\s+OS=|\s+\[.*?\]$|\s+\(.*?\)$")
_NONWORD = re.compile(r"[^a-z0-9]+")


def normalize_title(title: str) -> str:
    if not isinstance(title, str):
        return ""
    t = title.split("RecName:", 1)[-1]
    t = _ORG_SPLIT_RE.split(t)[0]
    t = t.lower()
    t = _NONWORD.sub(" ", t)
    toks = [w for w in t.split() if w and w not in _STOP]
    return " ".join(sorted(toks))


# ----------------------------------------------------------------------
# Consensus per query
# ----------------------------------------------------------------------
def consensus_per_query(hits: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = hits.groupby("qseqid", sort=False)
    for q, g in grouped:
        g = g.copy()
        g["norm"] = g["salltitles"].map(normalize_title)
        weight_by_group = g.groupby("norm")["bitscore"].sum().sort_values(
            ascending=False)
        count_by_group = g.groupby("norm").size()
        top_norm = weight_by_group.idxmax()
        top_weight = weight_by_group.iloc[0]
        total_weight = weight_by_group.sum()
        top_hits = g[g["norm"] == top_norm].sort_values(
            "bitscore", ascending=False)
        rows.append({
            "qseqid": q,
            "consensus_annotation": top_hits["salltitles"].iloc[0],
            "confidence_score_0to1": top_weight / max(total_weight, 1e-12),
            "support_fraction_bitscore": top_weight / max(total_weight, 1e-12),
            "agreement_fraction_hits": count_by_group[top_norm] / len(g),
            "top_hit_sseqid": top_hits["sseqid"].iloc[0],
            "top_hit_title": top_hits["salltitles"].iloc[0],
            "top_hit_organism": _extract_organism(top_hits["salltitles"].iloc[0]),
        })
    return pd.DataFrame(rows)


def _extract_organism(title: str) -> str:
    if not isinstance(title, str):
        return ""
    m = re.search(r"\[(.*?)\]", title)
    if m:
        return m.group(1)
    m = re.search(r"OS=([^=]+?)\s+\w+=", title)
    return m.group(1).strip() if m else ""


# ----------------------------------------------------------------------
# Source-table assembly
# ----------------------------------------------------------------------
def load_source_tables(path: Path) -> pd.DataFrame:
    xl = pd.ExcelFile(path)
    frames = []
    for sheet in ("Corrected", "Uncorrected"):
        if sheet not in xl.sheet_names:
            continue
        df = pd.read_excel(path, sheet_name=sheet)
        df["P_value_status"] = sheet
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main() -> int:
    if not CORRECTED_TBL.exists():
        print(f"[ERROR] missing input: {CORRECTED_TBL}")
        return 1

    sources = load_source_tables(CORRECTED_TBL)
    print(f"loaded {len(sources)} occurrence rows from "
          f"Corrected_Bonferroni table")

    # BLAST hits side
    if HITS_TSV.exists():
        hits = pd.read_csv(HITS_TSV, sep="\t")
        hits["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
        consensus = consensus_per_query(hits.dropna(subset=["salltitles"]))
        missing = (set(sources["Accession Number"].dropna().unique())
                   - set(consensus["qseqid"].unique()))
    else:
        print(f"[INFO] {HITS_TSV.name} not found; consensus will be loaded "
              f"from the source table's existing 'consensus_annotation' column")
        hits = pd.DataFrame()
        consensus = (sources[["Accession Number", "consensus_annotation",
                              "confidence_score_0to1",
                              "support_fraction_bitscore"]]
                     .drop_duplicates(subset=["Accession Number"])
                     .rename(columns={"Accession Number": "qseqid"}))
        missing = set()

    # Unique accessions: one row per query
    unique = (sources.drop_duplicates(subset=["Accession Number"])
              .merge(consensus, left_on="Accession Number", right_on="qseqid",
                     how="left", suffixes=("", "_blast")))
    if "qseqid" in unique.columns:
        unique = unique.drop(columns=["qseqid"])

    # Subject title cache
    if SUBJ_TSV.exists():
        subj = pd.read_csv(SUBJ_TSV, sep="\t")
    elif not hits.empty:
        subj = (hits[["sseqid", "salltitles"]]
                .drop_duplicates(subset=["sseqid"])
                .rename(columns={"salltitles": "title"}))
        subj["organism"] = subj["title"].map(_extract_organism)
        subj["length_aa"] = np.nan
    else:
        subj = pd.DataFrame()

    # Hits actually used for consensus (top-group hits per query)
    if not hits.empty:
        hits["norm"] = hits["salltitles"].map(normalize_title)
        top_norm_by_q = (hits.groupby("qseqid").apply(
            lambda g: g.groupby("norm")["bitscore"].sum().idxmax())
                         .to_dict())
        hits_used = hits[hits.apply(
            lambda r: top_norm_by_q.get(r["qseqid"]) == r["norm"], axis=1
        )].drop(columns=["norm"])
    else:
        hits_used = pd.DataFrame()

    with pd.ExcelWriter(OUT_PATH, engine="openpyxl") as xl:
        unique.to_excel(xl, sheet_name="unique_accessions", index=False)
        sources.to_excel(xl, sheet_name="occurrences", index=False)
        consensus.to_excel(xl, sheet_name="consensus_per_query", index=False)
        subj.to_excel(xl, sheet_name="subject_title_cache", index=False)
        hits_used.to_excel(xl, sheet_name="hits_used_for_consensus",
                            index=False)
        pd.DataFrame({"Accession Number": sorted(missing)}).to_excel(
            xl, sheet_name="missing_subjects", index=False)

    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
