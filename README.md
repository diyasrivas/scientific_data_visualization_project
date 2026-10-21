# Scientific Data Visualization Project

Course project for **Scientific Visualization** (6 ECTS, INFO-01/A), MSc Computer Science
(LM-18), Università degli Studi di Milano. Lecturer: Prof. Elena Casiraghi.

## The project in one paragraph

Singapore's resident population comprises three major groups (Chinese, Malay and Indian)
whose genetic ancestries diverged long before they shared a city. This project asks
whether the parts of the genome that differ most between those groups are also the parts
that matter functionally and clinically, or whether most of that divergence is
biologically inert. It joins three public resources: allele frequencies and haplotype
structure from the Singapore Genome Variation Project (SGVP), cell-type-resolved immune
gene expression from the Asian Immune Diversity Atlas (AIDA), and clinical annotation
from PharmGKB, CPIC and the GWAS Catalog. The deliverable is a set of
publication-standard figures that place each gene on a shared "divergence landscape",
together with a written report and a presentation.

## Why this shape of project

Three things had to hold at once:

1. **Close to the lecturer's interests without competing with her.** Multimodal data
   integration, dimensionality reduction and translational biomedical analysis are her
   territory. Population genetics statistics (FST, iHS, varLD) and pharmacogenomics are
   not, which leaves room for the student to hold genuine expertise.
2. **Genuinely multimodal.** Genotype, single-cell expression and curated clinical
   annotation are semantically different views of the same underlying entities (genes).
3. **Visualization is the deliverable, not decoration.** The central claim of the project
   is a picture: where genes sit in a two-dimensional space of genetic versus functional
   divergence.

## Design doctrine

Every figure in this repository is designed against the Gestalt principles of
**proximity**, **similarity** and **continuity**, and every figure is expected to survive
the question "what grouping is this image asking the viewer to believe?". See
[`docs/05-visualisation-and-gestalt.md`](docs/05-visualisation-and-gestalt.md).

Code follows research software engineering practice: configuration over hard-coded paths,
deterministic seeds, pinned environments, tested numerical kernels, and no analysis step
that cannot be re-run from a clean checkout. See
[`docs/06-software-engineering.md`](docs/06-software-engineering.md).

## Documentation

| Document | Purpose |
| --- | --- |
| [`docs/01-context-and-motivation.md`](docs/01-context-and-motivation.md) | Course context, supervisor fit, why this topic was chosen |
| [`docs/02-research-question.md`](docs/02-research-question.md) | Question, hypotheses, prior art, novelty, limitations |
| [`docs/03-data-sources.md`](docs/03-data-sources.md) | Every dataset, how to obtain it, licence, known risks |
| [`docs/04-analysis-workflow.md`](docs/04-analysis-workflow.md) | The pipeline, stage by stage, with acceptance checks |
| [`docs/05-visualisation-and-gestalt.md`](docs/05-visualisation-and-gestalt.md) | Figure specifications and the perceptual reasoning behind them |
| [`docs/06-software-engineering.md`](docs/06-software-engineering.md) | Repository layout, environment, testing, reproducibility |
| [`docs/07-project-plan.md`](docs/07-project-plan.md) | Milestones, risk register, deliverables |
| [`docs/08-proposal-to-supervisor.md`](docs/08-proposal-to-supervisor.md) | Draft proposal text to send to Prof. Casiraghi |

## Status

Planning. No data has been acquired yet. The first executable task is
[Stage 0 of the workflow](docs/04-analysis-workflow.md#stage-0-feasibility-gate), which
verifies that the SGVP download links still resolve before any further work is committed.

## Repository layout

```
.
├── README.md
├── docs/                 documentation (this planning material)
├── config/               YAML configuration, one file per pipeline run
├── src/sgdiv/            importable package: all reusable logic lives here
├── scripts/              thin numbered entry points that call into src/
├── notebooks/            exploration only, never a dependency of a result
├── tests/                pytest suite over the numerical kernels
├── data/                 gitignored: raw/, interim/, processed/
└── figures/              gitignored: generated vector output
```

Directories beyond `docs/` are created as the corresponding workflow stage begins.
