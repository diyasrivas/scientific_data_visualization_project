# 4. Analysis workflow

The pipeline is a linear sequence of numbered stages. Each stage reads from `data/` and
writes to `data/interim/` or `data/processed/`, is idempotent, is driven by a YAML config
rather than arguments scattered through the code, and ends with an **acceptance check**
that must pass before the next stage is attempted.

Entry points live in `scripts/` as thin wrappers (`scripts/01_acquire_sgvp.py` and so on).
All real logic lives in the importable package `src/sgdiv/` so that it can be tested.

```
Stage 0  feasibility gate
Stage 1  acquire SGVP                     ─┐
Stage 2  acquire and reduce AIDA           ├─> Stage 4 harmonise
Stage 3  acquire clinical annotation      ─┘
Stage 4  harmonise coordinates and identifiers
Stage 5  genetic divergence statistics    ─┐
Stage 6  expression divergence statistics  ├─> Stage 7 joint structure
                                          ─┘
Stage 7  joint embedding and structure
Stage 8  hypothesis testing
Stage 9  figures
Stage 10 report and talk
```

---

## Stage 0: feasibility gate

**Purpose.** Establish, before committing effort, that the data actually exists in
downloadable form. This stage exists because the SGVP download links are unverified (see
[`03-data-sources.md`](03-data-sources.md#32-sgvp-singapore-genome-variation-project)).

**Steps.**

1. Open the SGVP landing page in a normal browser session and record the literal target
   URL behind each of the listed file groups.
2. For each URL, issue a `HEAD` request and record status code, content type and content
   length. Priority order: hg19 genotypes, pharmacogenomic variants, allele frequencies,
   iHS, varLD.
3. Download the single smallest available file end to end and confirm it parses.
4. Record every result in `docs/stage0-link-audit.md`, including failures and dates.

**Decision rule.**

- hg19 genotypes **and** pharmacogenomic variants resolve: proceed as designed.
- Only frequencies resolve: proceed, but drop the LD and selection layers and say so in
  the proposal.
- Nothing resolves: switch to the gnomAD fallback, accept the loss of the Malay group,
  and email the SGVP contact authors to ask for an archived copy. Do this in week 1, not
  week 6.

**Acceptance check.** `docs/stage0-link-audit.md` exists and states one of the three
decisions above explicitly.

---

## Stage 1: acquire SGVP

**Inputs.** URLs confirmed in Stage 0.
**Outputs.** `data/raw/sgvp/`, entries in `data/MANIFEST.tsv`.

**Steps.**

1. Download each required file; compute SHA-256; append to the manifest with access date,
   source URL, declared genome build and licence.
2. Load genotypes with `pandas-plink` or convert with PLINK 2 to a compact format.
3. Record per-population sample counts and assert they match the published post-QC counts
   (96 `CHS`, 89 `MAS`, 83 `INS`). A mismatch means a different release and must be
   investigated, not absorbed.
4. Apply a documented variant filter: autosomes only, biallelic, call rate above 0.95,
   minor allele frequency above 0.01 in at least one population, Hardy-Weinberg
   p above 1e-6 within each population.

**Acceptance check.** Variant and sample counts before and after filtering are written to
a stage report, and the three pairwise FST values computed on the filtered set are within
a reasonable margin of the published values (Indian versus Chinese about 3.9 per cent,
Indian versus Malay about 2.7 per cent, Chinese versus Malay about 0.6 per cent). This is
the single most valuable check in the pipeline: it validates loading, filtering and the
FST implementation in one step against a published number.

---

## Stage 2: acquire and reduce AIDA

**Inputs.** CELLxGENE `.h5ad` (about 13.3 GB).
**Outputs.** `data/processed/aida_pseudobulk.parquet`, `data/processed/aida_donors.parquet`.

**Steps.**

1. Download the dataset object once; checksum it; never modify it in place.
2. Stream it with `anndata` in backed mode. Do not load the full matrix into memory.
3. Subset to donors with country Singapore and self-reported ethnicity in {Chinese,
   Malay, Indian}. Record resulting donor counts per group.
4. Harmonise the cell type annotation to a single level of granularity. Use the authors'
   own annotation column; do not invent a new one. Record which column was used.
5. Build **pseudobulk**: for each (donor, cell type) pair, sum raw counts across cells.
   Record the number of cells contributing to each pseudobulk sample.
6. Filter pseudobulk samples with fewer than a documented minimum number of cells
   (start at 10, state the threshold, and test sensitivity to it at Stage 8).
7. Normalise: counts per million, then `log1p`. Retain genes expressed above a minimum
   threshold in a minimum fraction of samples within at least one cell type.

**Why pseudobulk.** Donor is the unit of replication for a population-level question.
Treating individual cells as replicates would inflate significance by orders of magnitude,
and is the single most common error in this kind of analysis.

**Acceptance check.** Donor counts per population group are reported; the pseudobulk
matrix has no all-zero rows or columns; a principal component plot of pseudobulk samples
separates by cell type far more strongly than by population group, which is the expected
biology and confirms nothing has been scrambled.

---

## Stage 3: acquire clinical annotation

**Outputs.** `data/processed/gene_annotation.parquet`.

**Steps.**

1. Download PharmGKB gene-drug and clinical annotation files; record the release date.
2. Download the CPIC gene-drug pair table with assigned levels.
3. Download GWAS Catalog associations.
4. Resolve every gene identifier to a stable **Ensembl gene ID**, not a symbol. Gene
   symbols are not stable across releases and are the usual source of silent join errors.
5. Produce one tidy table: one row per Ensembl gene ID, with boolean and categorical
   flags for CPIC level, PharmGKB presence, and GWAS association count.

**Acceptance check.** A spot list of ten well-known pharmacogenes (`CYP2D6`, `CYP2C19`,
`CYP2C9`, `SLCO1B1`, `VKORC1`, `TPMT`, `DPYD`, `NUDT15`, `UGT1A1`, `G6PD`) is present and
correctly flagged. `CYP2D6` is a known hard case due to structural variation and should be
flagged for special handling rather than quietly trusted.

---

## Stage 4: harmonise coordinates and identifiers

This stage is where the project will most likely lose data silently, so it is treated as a
first-class stage with its own report rather than as glue code.

**Steps.**

1. Lift SGVP variant coordinates to **GRCh38** using `pyliftover` or the UCSC chain files.
   Record how many variants lift, fail, or lift to a different chromosome.
2. Resolve strand and allele representation. Harmonise to the GRCh38 reference allele.
   Flag ambiguous A/T and C/G variants explicitly; decide whether to exclude them
   (recommended for the primary analysis) and report the count either way.
3. Cross-check a random sample of lifted variants against gnomAD by rsID and allele
   frequency to confirm the lift did not systematically corrupt positions.
4. Assign each variant to genes: gene body plus a flanking window (default 5 kb, exposed
   in config, sensitivity-tested).
5. Join everything on Ensembl gene ID.

**Acceptance check.** A written report of lift success rate, ambiguous variant count, and
the gnomAD concordance check. If lift success falls below 95 per cent, stop and
investigate rather than proceeding.

---

## Stage 5: genetic divergence statistics

**Output.** `data/processed/genetic_divergence.parquet`, one row per gene.

**Per-variant measures.**

- **FST**, computed pairwise and across all three populations. Use the Hudson estimator
  with a documented formula; Weir and Cockerham as a sensitivity check. Report which was
  used and why, because the estimators differ non-trivially for small, unequal samples.
- **Maximum pairwise allele frequency difference**, which is less elegant than FST but far
  more interpretable when talking to a clinician.
- **iHS** and **varLD** scores, taken from the SGVP release rather than recomputed.

**Gene-level summarisation.** Summarise variants to genes using both the maximum and the
95th percentile of per-variant FST within the gene window. Report both. The maximum is
sensitive to single outlier variants; the percentile is more stable. Disagreement between
them is itself informative and worth a supplementary panel.

**Acceptance check.** The genome-wide FST distribution is unimodal with a long right tail;
the top differentiated regions recover known signals reported in the SGVP paper, including
`EDAR` and `VKORC1`. Failing to recover these means something upstream is wrong.

---

## Stage 6: expression divergence statistics

**Output.** `data/processed/expression_divergence.parquet`, one row per (gene, cell type).

**The core statistic.** For each gene within each cell type, fit a linear mixed model to
donor-level pseudobulk expression with population group as the effect of interest and age,
sex and batch as covariates, then compute the **variance partition coefficient** for
population: the proportion of residual expression variance attributable to population
group.

This deliberately mirrors FST. FST is the proportion of genetic variance attributable to
population structure; the variance partition coefficient is the proportion of expression
variance attributable to population structure. The two axes of the central figure are then
*the same quantity measured on two different modalities*, which is what makes the scatter
interpretable rather than arbitrary.

**Steps.**

1. Fit per gene and cell type. `variancePartition` in R is the reference implementation;
   `statsmodels` mixed models in Python are acceptable if the specification is documented.
2. Calibrate by permutation: shuffle population labels across donors (preserving cell type
   structure), refit, and build an empirical null for the variance partition coefficient.
   Minimum 1,000 permutations, seed recorded.
3. Compute an empirical p-value per gene and cell type; control the false discovery rate
   across genes within each cell type (Benjamini-Hochberg).
4. Retain only cell types with a documented minimum number of donors per population group.

**Negative control.** Housekeeping genes must show low population variance partition
coefficients. If they do not, the model is mis-specified or the covariates are wrong. Run
this check before looking at any result of interest.

**Acceptance check.** The permutation null is approximately uniform in p-value under label
shuffling, and the housekeeping control passes.

---

## Stage 7: joint embedding and structure

**Output.** `data/processed/gene_embedding.parquet`, stability diagnostics.

**Steps.**

1. Build a gene-level feature matrix: genetic divergence measures, plus per-cell-type
   expression divergence measures. Standardise features; document the transform.
2. Run PCA first and always. Report the explained variance curve. PCA is linear,
   deterministic and interpretable, and it is the honest baseline against which any
   non-linear embedding must justify itself.
3. Run UMAP, t-SNE and PHATE across a grid of hyperparameters and at least 10 random
   seeds each.
4. Quantify what each embedding preserves: trustworthiness, continuity, k-nearest-
   neighbour overlap with the original space, and Spearman correlation between embedded
   and original pairwise distances.
5. Measure stability of any apparent cluster across seeds and hyperparameters.

**This stage has a visualization ethic attached.** An embedding places points near each
other, and by the Gestalt principle of proximity a viewer will read that nearness as
group membership whether or not it is real. Any cluster that does not survive the seed
and hyperparameter sweep must not be drawn as a cluster. See
[`05-visualisation-and-gestalt.md`](05-visualisation-and-gestalt.md#53-the-proximity-problem-in-embeddings).

**Acceptance check.** Stability diagnostics exist for every embedding that appears in any
figure. No exceptions.

---

## Stage 8: hypothesis testing

**Steps.**

1. **H1.** Spearman rank correlation between gene-level genetic divergence and the maximum
   expression divergence across cell types. Report the confidence interval, not only the
   p-value.
2. **H2.** Compare pharmacogenes against a **matched background** rather than against all
   genes. Match on expression level, gene length and variant density, since all three
   confound differentiation statistics. Use at least 1,000 matched null sets.
3. **H3.** Compare the distribution of variance partition coefficients across cell types;
   test for concentration against a uniform expectation.
4. **H4.** Quantify the non-redundancy of the two axes: variance explained by one from the
   other, and the count of genes in each off-diagonal quadrant.
5. Sensitivity analyses, every one pre-declared: flanking window size, minimum cells per
   pseudobulk sample, FST estimator, inclusion of ambiguous variants, gene-level
   summarisation rule.

**Multiple testing.** Benjamini-Hochberg within each family of tests. State the family
definitions before running anything.

**Acceptance check.** Every reported number in the report traces to a specific stage
output file and a config hash. No number is typed by hand into the text.

---

## Stage 9: figures

Specified in full in [`05-visualisation-and-gestalt.md`](05-visualisation-and-gestalt.md).
Each figure is produced by one script, from processed data only, with no manual editing
afterwards, and is written as vector output.

---

## Stage 10: report and talk

- Report: question, data, method, results, limitations (the ecological caveat stated
  prominently), and what the student would do next with individual-level data.
- Talk: built around the central divergence landscape figure. Open with the question,
  show the map, then walk through what each quadrant means.
- Prepare one slide on the embedding stability audit. It is the slide most likely to be
  asked about by a lecturer who teaches these algorithms, and having run the experiment
  turns a hard question into the best moment of the presentation.
