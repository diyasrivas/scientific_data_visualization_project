# 2. Research question, hypotheses and prior art

## 2.1 Primary question

> Among the loci that differ most in allele frequency between Singapore's Chinese, Malay
> and Indian populations, which also show population differences in gene expression at
> the level of specific immune cell types, and are those loci enriched or depleted for
> clinically actionable pharmacogenes?

## 2.2 Sub-questions

1. **Divergence structure.** How is genetic differentiation distributed across the genome
   for these three populations, and does it concentrate in particular functional classes
   of gene?
2. **Functional echo.** Do genes with high genetic differentiation also show high
   population-attributable variance in expression, and is that relationship uniform
   across immune cell types or concentrated in a few?
3. **Clinical relevance.** Are pharmacogenes (PharmGKB, CPIC level A and B) more
   differentiated than a matched genomic background, on either axis?
4. **Cell-type specificity.** Which immune cell types contribute most of the
   population-attributable expression variance, and does that correspond to anything
   biologically coherent?

## 2.3 Hypotheses

| # | Hypothesis | Direction | How it would be falsified |
| --- | --- | --- | --- |
| H1 | Genetic and expression divergence are positively but weakly associated at gene level | positive, small effect | A null or negative rank correlation after background matching |
| H2 | Pharmacogenes are *not* systematically more differentiated than matched background, i.e. clinical actionability is largely orthogonal to population structure | null | Significant enrichment after matching on expression level, gene length and SNP density |
| H3 | Population-attributable expression variance is concentrated in a minority of cell types rather than spread uniformly | concentrated | A flat distribution of variance partition coefficients across cell types |
| H4 | A subset of genes shows high genetic but negligible expression divergence (and vice versa), making the two axes genuinely non-redundant | non-redundant | Near-collinearity of the two divergence measures |

H2 is deliberately stated as a null hypothesis. A clean, well-visualised null is a
legitimate and more honest result than a fished-for positive, and the project is designed
so that either outcome produces the same figure.

## 2.4 The critical caveat: this is an ecological comparison

SGVP and AIDA are **different cohorts of different people**. Genotype and expression are
therefore not observed in the same individuals, and no individual-level genotype to
expression inference is possible from this design. Every comparison in this project is at
the level of the *population group*, not the person.

Consequences, all of which must be stated in the report and on the relevant figures:

- The analysis cannot establish that a variant causes an expression difference. It can
  only establish that the two kinds of divergence coincide, or fail to coincide, at gene
  level across three population groups.
- With three population groups, the effective sample size for any correlation across
  populations is three. Gene-level statistics are therefore pooled across genes, and
  population-level claims are made descriptively.
- Confounding by cohort (recruitment, age structure, sequencing platform, batch) cannot
  be fully removed. Age, sex and batch are adjusted for within AIDA; cohort effects
  between AIDA and SGVP are acknowledged, not corrected.

Stating this clearly is a strength of the proposal. It is the kind of limitation the
lecturer works on directly in her own missing-data and integration research, and
pretending otherwise would be the fastest way to lose credibility in the viva.

## 2.5 Prior art

Work that already exists and must be cited:

- **Teo et al., *Genome Research* 19:2154 (2009).** The SGVP resource paper. Establishes
  population structure, FST between the three groups (Indian versus Chinese 3.9 per cent,
  Indian versus Malay 2.7 per cent, Chinese versus Malay 0.6 per cent), LD decay,
  tagging efficiency, varLD regions, and iHS and XP-EHH selection scans. Notably already
  flags `VKORC1` and `CDKAL1` as differentiated.
- **Chan et al., *Pharmacogenomics* (2010), doi:10.2217/pgs.10.79.** Uses SGVP and HapMap
  allele frequencies to assess population divergence across 42 PharmGKB genes. This is
  the closest existing work to sub-question 3 and the one to position against explicitly.
- **Nature *Scientific Reports* 14 (2024), doi:10.1038/s41598-024-67691-6.** CYP2D6 star
  allele landscape across Chinese, Malay and Indian individuals in Singapore, from
  whole-genome sequencing.
- **Kock et al., *Cell* (2025), doi:10.1016/j.cell.2025.02.017.** The AIDA resource paper.
  Reports that sub-continental diversity, age and sex pervasively affect immune cell
  properties, and identifies cell-type-specific expression quantitative trait loci that
  are under-represented in non-Asian populations.

## 2.6 What is new here

None of the above connects **genome-wide population differentiation** to **cell-type
resolved expression divergence** for the same three Singaporean population groups, and
none renders that relationship as a single interpretable visual space. Specifically:

1. The 2010 pharmacogenomics work predates single-cell data entirely and is restricted to
   42 genes. This project is genome-wide and adds a functional axis.
2. The AIDA paper analyses expression and eQTLs but does not relate them to classical
   population-differentiation or selection statistics.
3. The "divergence landscape" as a visual object, with clinical annotation overlaid, does
   not exist in the literature for this population triad.

The contribution is therefore **integrative and visual**, not a new algorithm. That is
the appropriate scope for a visualization course project, and it should be described that
way rather than overclaimed as methodological novelty.

## 2.7 Optional module: intrinsic dimension

If time allows, one additional question is worth asking because it connects directly to
the lecturer's own research:

> What is the intrinsic dimension of the gene-level divergence feature space, and does it
> differ between pharmacogenes and background genes?

This should be framed to her as a question, not a claim: which estimators would she
recommend, and how should small-sample-size effects be handled. Use at least two
estimators (for example TwoNN and MLE, with DANCo via `scikit-dimension`) with bootstrap
intervals, and report disagreement between them rather than a single number. Treat this
as a supporting lens, never as the headline result.
