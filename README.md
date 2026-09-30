# Code Supplement
## Proteomic Insights into Sessility and Lifestyle Evolution in Rotifers (Grajeda et al.)

Python scripts that generate every quantitative figure in the manuscript.
Running them in numerical order reproduces:

- Figure 1              (PCA, colour = species, shape = lifestyle)
- Figure 2              (top-11 functional classes per species, 8 panels)
- Figure 3 (panels A–D) (volcano plots, BH-FDR q = 0.05 + |log2FC| ≥ 2)
- Figure 4              (mucin / glycan-associated / surface subclasses, totals above bars)
- Supplemental Fig S1   (protein-ID count per sample, grouped by species)
- Supplemental Fig S2   (PCA, colour = species, shape = location)
- Supplemental Fig S3   (uncorrected volcano plots, panels A–D)
- Supplemental Fig S4   (contractile / muscle-associated subclasses)
- Supplemental Fig S5   (cytoskeleton / motility subclasses)

Every figure is written as a vector PDF as well as PNG/TIFF. No figure
carries an in-figure title or figure number: the description lives in the
typeset caption. Species names are italicised with an upright "sp."
(`common.italic_species`).

---

## Script map

| Script | Produces | Figure |
|---|---|---|
| `01_data_prep.py` | `processed_protein_data.xlsx` per comparison | inputs for all stats |
| `02_stats_and_corrections.py` | `Sig with/without correction` sheets, `statistics_summary.csv` | stats |
| `03_partial_rda.py` | `rda_results.xlsx` (RDA1–6 loadings) per comparison | supporting analysis |
| `04_volcano_plot.py` | per-panel `_volcano_{corrected,uncorrected}.{png,svg,pdf}` + both composites | **Figure 3, Sup S3** |
| `05_consensus_annotation.py` | `All_Accession_Numbers_compiled_with_consensus.xlsx` | annotation |
| `06_subclass_panels.py` | `Figure_4_*`, `Supplemental_Figure_S4/S5_*.{png,pdf,tiff}` | **Figure 4, Sup S4 / S5** |
| `06a_validate_subclass_counts.py` | cross-checks every published number | validation |
| `06b_subclassify_proteins.py` | keyword sub-classification rules | rules |
| `07_pca.py` | `Figure_1_*`, `Supplemental_Figure_S2_*` | **Figure 1, Sup S2** |
| `08_protein_ids_per_sample.py` | `Supplemental_Figure_S1_*` | **Sup S1** |
| `09_top10_functional.py` | `Figure_2_*` + the three summary CSVs | **Figure 2** |
| `10_mucin_by_species.py` | `Auxiliary_mucin_category_by_species.*` | auxiliary check (not a manuscript figure) |
| `11_database_robustness.py` | `Supplemental_Table_S2_species_by_source_genome.csv`, `database_robustness_summary.txt` | **Sup Table S2**; response to Reviewer 1 |

`common.py` holds shared helpers. `_subclass_rules.py` is a small import
shim that lets `06` and `10` import from `06b` (whose filename starts
with a digit and so cannot be imported directly).
`samples_metadata.csv` is the sample → species / location / lifestyle
lookup.

## Inputs expected in the project root

- `Volcano (1)/Volcano/<comparison>/Samples Report of My Experiment.csv` (or `.xlsx`)
- `Volcano (1)/Volcano/<comparison>/processed_protein_data (NN).xlsx` — the
  per-pair Scaffold processing carrying pre-filtered `Sig with correction`
  and `Sig without correction` sheets and the Raw Extract FWER. Script 04
  reads these directly when present.
- `Total_Protein_Output_052725.xlsx` — master per-sample protein matrix
- `Corrected_Bonferroni_completed_table_with_manual_uniprot_accessions_031226 (1).xlsx`

## How to run

```
pip install pandas numpy scipy scikit-learn matplotlib openpyxl \
            adjustText pillow pypdf reportlab
```

```bash
cd Code_Supplement

python 01_data_prep.py
python 02_stats_and_corrections.py
python 03_partial_rda.py
python 09_top10_functional.py       # Figure 2 + the three summary CSVs
python 06a_validate_subclass_counts.py   # asserts every published number
python 06_subclass_panels.py        # Figure 4 + Sup S4 / S5
python 10_mucin_by_species.py       # auxiliary per-species count check
python 07_pca.py                    # Figure 1 + Sup S2
python 08_protein_ids_per_sample.py # Sup S1
python 04_volcano_plot.py           # Figure 3 + Sup S3 (panels + composites)
python 05_consensus_annotation.py
python 11_database_robustness.py    # Sup Table S2 + reviewer-response numbers (after 09)
```

**Order note:** `09_top10_functional.py` writes
`rotifer_species_unique_protein_presence_with_mucin_contractile.csv`,
which `06`, `06a`, `06b` and `10` all read, so it must run before them.

## Ordering matters in the keyword classifiers

Assignment is first-match-wins, so any term that contains another as a
substring has to be tested first. Two such pairs occur:

```
"paramyosin"  contains "myosin"   -> tested before it
"tropomyosin" contains "myosin"   -> tested before it
```

`06a_validate_subclass_counts.py` asserts that the rules reproduce the
published per-species totals, so an ordering regression fails loudly
rather than silently redistributing proteins between subclasses.

## How Figure 2 is reproduced

Script 09 implements the 13-category classifier described in the Methods
(first match wins), displays the top 11 classes per species and collapses
the remainder into "Collapsed minor categories". It classifies the
**cleaned base protein name** — the Scaffold annotation truncated at the
first `|`, with `PREDICTED:` / `LOW QUALITY PROTEIN:` prefixes removed —
not the full annotation string. Classifying the full string would let GO
terms such as "striated muscle thin filament" pull actins into the
contractile class.

Per-species totals count proteins carrying an accession. The
`For_Patrick` sheet has 1,676 rows; the last is Scaffold's per-sample
summary footer and has no accession, so it is excluded — which is why the
totals are 540 / 971 / … and not 541 / 972 / ….

## How Figure 4 and Figures S4 / S5 are reproduced

Script 06 computes the subclass composition at run time from the presence
CSV using the rules in 06b, asserts the result against the published
per-species totals, and only then draws. (An earlier version hard-coded
the counts, which allowed the figures, the CSV and the rules to drift
apart.)

Myosins whose annotation specifies striated or sarcomeric muscle are
assigned to the contractile class (S4); myosins without a muscle-type
qualifier fall in the cytoskeleton class (S5) under
"Myosin (non-muscle/unspec.)". This is why *Hexarthra* sp. contributes a
single myosin to S4 while carrying many myosin-annotated proteins overall.

## How Figure 3 is reproduced

Each comparison was processed in **Scaffold DDA 6.5.0**: t-test
(Scaffold's weighted-spectrum-count test), multiple-test correction by
Benjamini–Hochberg at **q = 0.05**. Script 04 then:

1. reads the per-pair `processed_protein_data (NN).xlsx`;
2. pulls the **FWER** from its Raw Extract sheet and draws the horizontal
   cutoff at `−log10(FWER)` — A 6.51e-03, B 6.55e-03, C 6.35e-05, D 5.12e-05;
3. pulls the BH-significant accessions from `Sig with correction`;
4. applies the `|log2FC| ≥ 2` effect-size filter;
5. colours and labels by the EnhancedVolcano convention, labelling the
   ten most up- and ten most down-regulated red points per panel.

**log2FC direction.** `log2FC = log2(mean(group 1) + 1) − log2(mean(group 0) + 1)`.
The (0) group is the one named **first** in each panel title, so a
**negative** log2FC means enrichment in the first-named species.

Panels C and D are listed in `EFFECT_SIZE_ONLY`, which suppresses the red
"p-value and log2 FC" category for those panels; proteins that would
otherwise be red are drawn green. Panel C's legend therefore reports
`n = 0` red even though three proteins pass BH-FDR there (two are drawn
blue, one green).

The composite is assembled by merging the four panel **PDFs**, so
Figure 3 and Sup S3 are vector and the accession labels stay legible at
any zoom. PNG/TIFF companions are rasterised from the same page.

## How Figure 1 / Sup S2 are reproduced

Script 07 runs `sklearn.decomposition.PCA` on `log2(x + 1)` with
protein-wise mean centring. Colour = species (`matplotlib.tab10`, keyed
to alphabetical species order); shape = lifestyle (Figure 1) or location
(Sup S2). Both legends sit outside the axes and are passed to
`savefig(bbox_extra_artists=…)`; without that, `bbox_inches="tight"`
clips them and the longer species names lose their final characters.

## How Sup S1 is reproduced

Script 08 draws one bar per sample, grouped into species blocks separated
by a gap and a dashed line, with species names above the axes (two
lines: genus / epithet) and two-line tick labels below: collection site,
then "sample N" (the biological-sample number used in the PRIDE deposit
and in `rotifer_protein_ids_by_species.csv`).

## Conventions used throughout

- **Group encoding.** Scaffold tags sample columns with a leading `(0)` or
  `(1)`; unprefixed columns are out-of-group reference samples, used as
  covariates in the RDA only.
- **Single Proteins vs Protein Groups.** Scaffold's parsimonious cluster
  header rows carry no Molecular Weight; per-protein statistics use the
  "Single Proteins" sheet only.
- **Log2 transformation.** Statistics and ordination use `log2(x + 0.5)`;
  PCA uses `log2(x + 1)`.
- **`"< 0.0001"` parsing.** Treated as the literal 0.0001 boundary, so
  those proteins plot at y = 4.0.

## Numerical fidelity

`06a_validate_subclass_counts.py` checks all of the following and exits
non-zero on any mismatch:

| Quantity | Value |
|---|---|
| Non-redundant proteins | 1,675 |
| Identifications per sample | 178 – 959, across 20 samples |
| Unique proteins per species | 540, 971, 624, 959, 957, 692, 672, 432 |
| Mucin / glycan-associated / surface (Figure 4) | 0, 8, 10, 4, 2, 3, 1, 1 |
| Contractile / muscle-associated (S4) | 11, 21, 12, 22, 15, 14, 9, 3 |
| Cytoskeleton / motility (S5) | 105, 120, 96, 111, 109, 101, 93, 82 |
| Figure 1 / Sup S2 variance | PC1 21.7 %, PC2 16.7 % |
| Figure 3 red points (BH-FDR + \|log2FC\| ≥ 2) | A 26, B 28, C 0\*, D 0 |

Species order for the subclass rows is C. vorax, S. socialis,
L. flosculosa, N. copeus, C. hippocrepis, P. patulus, E. kingi,
Hexarthra sp.

\* Panel C has three BH-significant proteins; none is drawn red because
the panel is in `EFFECT_SIZE_ONLY` (see above).

## Database composition and the unannotated class

The search database concatenates thirteen predicted proteomes: nine from
Mohl et al. (2025), whose FASTA headers carry OmicsBox (BLASTP + InterProScan)
descriptors, and four unpublished laboratory genomes (`cono` = *Conochilus
hippocrepis*, `fl` = *Filinia longiseta*, `bv` = *Brachionus variabilis*,
`pquad` = *Platyias quadricornis*) whose entries carry **no descriptor** —
in the Scaffold export their protein name is the accession itself. The
keyword classifier in script 09 therefore places all 413 proteins
identified through those four sets (24.7 % of 1,675) in
"Unknown / unannotated", together with 120 published-genome entries
annotated only as "hypothetical protein" or with a blank descriptor
(533 in total, 31.8 %). Script 11 quantifies the consequences: source
genome of every identification by species, detected-versus-predicted
correlation, direction and source of every Figure 3 red point, and
class-level abundance shares that are independent of which homologous
entry a peptide was assigned to.

## Note on the RDA implementation

`03_partial_rda.py` is a plain NumPy partial RDA (centre → regress out
covariates → SVD of the fitted matrix). With a single binary group label
only RDA1 carries non-degenerate loadings. The published
`*_rda_vs_ttest_*.xlsx` workbooks were produced against the full
four-species matrix in R/vegan, which yields six non-degenerate axes:

```r
vegan::rda(t(spectral_counts) ~ species, data = meta)
```

The RDA1 ranking from script 03 reproduces the manuscript's RDA1 ranking;
the multi-axis decomposition does not.

Contact: Brian Grajeda (bgrajeda16@gmail.com), Walsh lab, UTEP.
