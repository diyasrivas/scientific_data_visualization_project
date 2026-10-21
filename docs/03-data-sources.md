# 3. Data sources

Four sources, three of them genuinely multimodal views of the same gene universe and one
providing clinical annotation. Nothing in this project requires controlled access, which
was a hard selection criterion: a semester project cannot wait on a data access
committee.

## 3.1 Summary table

| ID | Resource | Modality | Unit | Access | Role |
| --- | --- | --- | --- | --- | --- |
| `SGVP` | Singapore Genome Variation Project | SNP genotypes, LD, selection | variant | Open | Genetic divergence axis |
| `AIDA` | Asian Immune Diversity Atlas | single-cell RNA-seq | cell, donor | Open (CC BY 4.0) | Expression divergence axis |
| `PGX` | PharmGKB and CPIC | curated annotation | gene, gene-drug pair | Open (registration for some files) | Clinical actionability layer |
| `GWAS` | NHGRI-EBI GWAS Catalog | curated annotation | variant-trait | Open | Disease context layer |
| `GNOMAD` | gnomAD v4 | population allele frequencies | variant | Open | Cross-check and fallback |

## 3.2 SGVP (Singapore Genome Variation Project)

**Citation.** Teo YY, Sim X, Ong RTH, et al. Singapore Genome Variation Project: a
haplotype map of three Southeast Asian populations. *Genome Research* 19:2154-2162 (2009).
doi:10.1101/gr.095000.109

**Content.** 292 individuals genotyped (99 Chinese, 98 Malay, 95 Indian) on both the
Affymetrix SNP6.0 and Illumina 1M arrays. After quality control, 268 individuals remain
(96 Chinese `CHS`, 89 Malay `MAS`, 83 Indian `INS`), with roughly 1.58 million autosomal
SNPs per population and a mean inter-SNP distance of about 2 kb.

**Files listed on the project download page:**

| File group | Contents | Notes |
| --- | --- | --- |
| Genotypes (hg18) | individual genotype data | original release build |
| Genotypes (hg19) | individual genotypes in PLINK and VCF | **preferred starting point** |
| Frequencies | allele and genotype frequencies | per population |
| LD Data | D', LOD, r-squared | per population |
| Pharmacogenomic variants | 3,972 SNPs in drug disposition and response, 253 individuals | **directly relevant to this project** |
| Phasing Data | PHASE output | |
| Recombination Rates | rates and genetic distances | |
| Selection | iHS scores | |
| varLD | genome-wide varLD scores, 17 population pairs | |
| Copy Number Variation | CNV loci | |
| HLA | HLA genotypes and phased haplotypes, MHC imputation reference panel | optional extension |

Landing page: `https://blog.nus.edu.sg/sshsphphg/singapore-genome-variation-project`

> **Known risk, resolve first.** As of October 2026 the landing page content is readable,
> but an automated browser request to it was redirected to an institutional sign-in page,
> and the individual file URLs have not been confirmed to resolve. The genome browser
> referenced in the 2009 paper (`www.nus-cme.org.sg/SGVP/`) is from a host that may no
> longer be live. **Stage 0 of the workflow exists solely to settle this before any other
> work begins.** See [`04-analysis-workflow.md`](04-analysis-workflow.md#stage-0-feasibility-gate).

**Fallback if SGVP files are unobtainable.** Substitute gnomAD v4 population allele
frequencies for East Asian (`eas`) and South Asian (`sas`) groups. This is a real loss:
gnomAD has no Malay group, so the three-population design collapses to two, and the LD,
varLD and iHS layers disappear. The project remains viable but the Malay population, which
is precisely the group most under-represented in global reference data, would be lost. The
fallback is documented so the decision is explicit, not so it is taken lightly.

**Build.** Original release is hg18; an hg19 release exists. AIDA and all modern
annotation are on GRCh38. Coordinate lifting is therefore a required, and error-prone,
pipeline stage. See workflow Stage 4.

## 3.3 AIDA (Asian Immune Diversity Atlas)

**Citation.** Kock KH, Tan LM, Han KY, et al. Asian diversity in human immune cells.
*Cell* (2025). doi:10.1016/j.cell.2025.02.017

**Content.** 1,265,624 circulating immune cells from 619 donors across 7 population
groups in 5 Asian countries, plus 6 controls. Some press material quotes 625 donors; use
the figure from the paper (619 donors and 6 controls) and note the discrepancy rather than
silently picking one. Includes Singaporean donors self-identifying as Chinese, Malay and
Indian, which is exactly the triad needed to align with SGVP.

**Access.**

- Processed object on CZ CELLxGENE: collection
  `ced320a1-29f3-47c1-a735-513c7084d508`. The "AIDA Phase 1 Data Freeze v2" dataset is a
  single `.h5ad` of 1,265,624 cells by 35,477 genes, roughly 13.3 GB.
- Human Cell Atlas Data Explorer project
  `f0f89c14-7460-4bab-9d42-22228a91f185` (raw FASTQ, 369 donors listed there, not needed).
- Analysis code: `https://github.com/prabhakarlab/AIDA_Phase1`
- Licence: CC BY 4.0 under the HCA Data Release Policy.

**Practical note.** Do not work with the 13.3 GB object directly. Subset to Singaporean
donors and collapse to donor-by-cell-type pseudobulk as early as possible (workflow Stage
2); the resulting matrices are small enough to analyse on a laptop and to version
sensibly.

**Donor metadata needed.** Self-reported ethnicity, country, age, sex, and any batch or
library identifier. Age and sex are confirmed by the paper to affect immune cell
properties, so they are covariates, not nuisance to be ignored.

## 3.4 Clinical annotation

**PharmGKB** (`https://www.pharmgkb.org/downloads`). Gene-drug relationships, clinical
annotations with evidence levels, and variant annotations. Downloads are free; some
require accepting a licence for academic use. Cite the current data release date.

**CPIC** (`https://cpicpgx.org`). Gene-drug pairs with assigned levels (A, B, C, D). CPIC
level A is the operational definition of "clinically actionable" used throughout this
project, because it is a published, auditable criterion rather than a judgement call.

**NHGRI-EBI GWAS Catalog** (`https://www.ebi.ac.uk/gwas/docs/file-downloads`). Used only
for context on the case-study loci, not as a primary analysis layer. Be aware of its
well-documented European ancestry bias, which is itself worth one sentence in the report
given the subject matter.

## 3.5 gnomAD v4

Used for two purposes: sanity-checking SGVP allele frequencies against an independent
modern reference, and serving as the documented fallback described in 3.2. Frequencies
for East Asian and South Asian groups are the relevant comparators.

## 3.6 Data governance for this repository

- **Nothing in `data/` is ever committed.** The directory is gitignored.
- Every acquired file is recorded in `data/MANIFEST.tsv` with: logical ID, source URL,
  access date, SHA-256 checksum, file size, genome build, and licence.
- Acquisition is performed by a script in `scripts/`, never by hand, so that the manifest
  and the download cannot drift apart.
- Re-running an acquisition script against an unchanged manifest entry is a no-op; a
  checksum mismatch is a hard failure, not a warning.

This matters beyond tidiness. SGVP is a 2009 resource hosted on university
infrastructure, and the probability that a link rots during the semester is not small. A
checksummed local copy with a recorded access date is the only defence.
