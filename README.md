# tcr-immunodominance

**Does TCR repertoire availability track how broadly an epitope is recognised?**

An R pipeline that profiles the T-cell receptor repertoire behind individual
epitopes, applied to SARS-CoV-2 structural proteins, with HSV-2 glycoprotein D
as a case study in what happens when the data run out.

![Depth control](figures/fig2_depth_control.png)

> **New to this topic?** [**WORKFLOW.md**](WORKFLOW.md) explains the whole
> project in plain language with diagrams — no immunology background needed.

---

## Why

Epitope prediction for vaccine design is dominated by the MHC side — *will this
peptide be presented?* The other half is usually left out: *is there a T-cell
repertoire available to recognise it?* Immunodominance is set on both sides. A
well-presented peptide seen by few available clones will not dominate.

This pipeline measures the second half from public TCR data.

## Findings

**1. Repertoire sharing tracks response breadth, and survives depth control.**

| Measure | Control | rho vs donor breadth | n |
|---|---|---|---|
| Exact-match publicness | none | +0.603 | 17 |
| Cluster publicness | none | +0.658 | 17 |
| **Rarefied to 3 donors** | depth-matched | **+0.553** | 17 |
| **Rarefied to 5 donors** | depth-matched | **+0.545** | 16 |
| Permutation z-score | depth-free | +0.333 | 10 |

Publicness rises with donor count for purely combinatorial reasons — with 81
donors there are 3,240 pairs a sequence could be shared across; with 3 donors
there are 3. Two independent controls (donor rarefaction, and a permutation
null on deduplicated records) both leave the association intact.

*Caveats, stated plainly: n = 17, HLA is uncontrolled, and donor count is a
proxy for breadth rather than a measurement of it. This is an association
consistent with the hypothesis, not a test of it. No p-values are reported
because the sample cannot support them.*

**2. Exact CDR3 matching makes publicness unmeasurable.**

The mean number of donors per unique CDR3 is 1.0–1.17 across every analysable
epitope. That is not evidence of private repertoires: people converge on
**motifs**, not identical sequences. Clustering similar CDR3s raises the median
public fraction from 0.007 to 0.018 (edit distance ≤ 1) and 0.044 (≤ 2).

![Exact vs cluster](figures/fig1_exact_vs_cluster.png)

**3. No repertoire is more convergent than its own sampling structure
predicts.** No epitope reaches z > 2. Two are significantly *donor-structured*
(z < −2) — clonal families inside individuals, which a naive analysis would
misread as public sharing.

**4. The method cannot be applied to the target that needs it.**

HSV-2 glycoprotein D, a leading genital herpes vaccine antigen, has **one
epitope in VDJdb: `SLPITVYYA`, HLA-A\*02:01, 4 TCRs, 1 study, no donor IDs.**
Not one metric above is computable for it.

![Data scarcity](figures/fig3_data_scarcity.png)

That is the result, not a failure. The approach works on SARS-CoV-2 spike
because a pandemic funded the sequencing. **Repertoire data exists in
proportion to attention already paid, not to clinical importance** — so
database-driven epitope ranking is unavailable precisely where vaccine design
would most benefit. The way around it is a model that needs no observed TCRs:
generation probability from a recombination model (OLGA/SONIA), computable for
any CDR3.

---

## Running it

```bash
conda env create -f environment.yml && conda activate tcr-immuno
Rscript tests/testthat.R               # 147 tests, fully offline

Rscript scripts/01_download_data.R     # VDJdb chunks + UniProt references
Rscript scripts/02_parse_vdjdb.R       # parse, annotate, QC
Rscript scripts/03_epitope_features.R  # per-epitope features
Rscript scripts/04_cdr3_networks.R     # similarity networks
Rscript scripts/05_depth_controlled.R  # rarefaction + permutation null
Rscript scripts/06_figures.R           # figures
# full write-up - either of:
Rscript -e 'rmarkdown::render("reports/report.Rmd")'   # no extra installs
quarto render reports/report.qmd                       # needs quarto
```

Already have an R/bioinformatics conda environment? Only three packages are
typically missing, all leaves:

```bash
conda install -c conda-forge r-stringdist r-testthat r-ggseqlogo
```

---

## Method notes

**VDJdb chunk files are read directly** rather than the assembled table. The
per-study chunks carry donor IDs as real columns, pair α/β chains, and expose
verified non-binders. No Python assembler toolchain needed.

**Two exclusions that change the answer:**

- *In vitro display.* 36,795 records against one epitope came from phage
  display — engineered sequences, not T cells from any person. One such
  experiment outweighs every genuine repertoire in the dataset.
- *Record counts ≠ response magnitude.* `LLLDRLNQL` has 1,753 records from 3
  donors; `TFEYVSQPFLMDLE` has 2,220 from 70. **Donor count is the breadth
  proxy throughout.**

**Publicness is `NA` with one donor, never 0.** With a single donor nothing
*can* be shared, and a 0 would read as "this repertoire is private".

**Rarefaction never extrapolates.** Epitopes shallower than the target depth
get `NA` rather than a biased estimate.

**`immunarch` is not a dependency.** Shannon, Simpson, Chao1, Hill numbers,
clonality, rarefaction (Hurlbert's analytic expectation) and overlap are
implemented in `R/diversity.R`.

## Known limitations

- β chain only; amino-acid CDR3s with no nucleotide column, so convergent
  recombination is not computable
- Rarefied richness is at ceiling (15/17 epitopes at 29–30 of a maximum 30) and
  contributes no variance
- Single-linkage clustering can chain; `largest_cluster_fraction` is reported
  as a diagnostic
- `LTDEMIAQY` is one deeply resampled expansion (3,582 records, 276 unique
  sequences, clonality 0.86) and distorts any statistic it enters
- Database bias is not fully removable: epitopes have records because they were
  studied, and were studied partly because they were thought dominant

---

## Layout

```
R/                      reusable functions
  utils.R               paths, logging, IO
  load_vdjdb.R          assembled-table / web-export loader
  load_vdjdb_chunks.R   chunk-format loader (preferred)
  annotate_antigen.R    epitope -> protein by sequence, not by label
  diversity.R           diversity and richness estimators
  features.R            per-epitope feature table
  cdr3_network.R        similarity networks, cluster publicness
  publicness_null.R     donor rarefaction + permutation null
scripts/                numbered, run in order
reports/report.qmd      full write-up
WORKFLOW.md             plain-language walkthrough with diagrams
tests/testthat/         147 tests, no network required
```

Functions live in `R/`; scripts only orchestrate.

## Data sources

| Source | Use |
|---|---|
| [VDJdb](https://vdjdb.cdr3.net) | TCR–epitope pairs (CC-BY) |
| [UniProt](https://www.uniprot.org) | Reference proteins (CC-BY 4.0) |

Shugay M *et al.* VDJdb: a curated database of T-cell receptor sequences with
known antigen specificity. *Nucleic Acids Research* 2018;46:D419–D427.

## Licence

MIT (code). Data retain their original licences.
