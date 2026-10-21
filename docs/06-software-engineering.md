# 6. Research software engineering practice

The standard applied here is simple: **a competent stranger should be able to clone this
repository, run one command, and obtain every number and every figure in the report.**
Everything below exists to serve that sentence.

## 6.1 Repository layout

```
.
├── README.md
├── pyproject.toml            project metadata, dependencies, tool config
├── environment.yml           conda environment with pinned versions
├── Makefile                  the pipeline DAG as targets
├── config/
│   ├── default.yaml          every tunable parameter, no exceptions
│   └── sensitivity/*.yaml    variants used for sensitivity analyses
├── src/sgdiv/
│   ├── io/                   loading, checksums, manifest handling
│   ├── genetics/             FST, allele frequency, liftover, gene assignment
│   ├── expression/           pseudobulk, variance partition, permutation
│   ├── integrate/            joins, feature matrix construction
│   ├── embed/                dimensionality reduction and stability metrics
│   └── plotting/
│       ├── style.py          the single source of visual truth
│       └── figures/          one module per figure
├── scripts/                  00_..., 01_..., thin CLI entry points
├── notebooks/                exploration only
├── tests/
├── data/                     gitignored
└── figures/                  gitignored
```

**The rule that matters:** `scripts/` contains argument parsing and orchestration only.
If a script contains logic, that logic cannot be tested, and it will be wrong.

## 6.2 Environment and dependencies

- Pin everything. `environment.yml` with explicit versions, and a lockfile committed.
- Record the environment used for any result: the pipeline writes `conda list` output
  alongside each stage report.
- Core stack: `numpy`, `pandas`, `scipy`, `scikit-learn`, `statsmodels`, `anndata`,
  `scanpy`, `pandas-plink`, `pyliftover`, `umap-learn`, `openTSNE`, `phate`,
  `scikit-dimension`, `matplotlib`, `pyarrow`.
- R is used only if `variancePartition` is preferred for Stage 6. If so, pin it with
  `renv` and document the handoff format (parquet in, parquet out). Do not create a
  half-R, half-Python pipeline with ad hoc CSV exchange.

## 6.3 Configuration

- **No hard-coded paths or parameters anywhere in `src/`.** All of it comes from the YAML
  config, loaded once and passed down.
- Every stage writes the **hash of the config it ran under** into its output directory.
  Any figure or number can therefore be traced to the exact parameters that produced it.
- Sensitivity analyses are config variants, not edited code. Changing the flanking window
  from 5 kb to 10 kb means a new YAML file, not a modified line.

## 6.4 Reproducibility

- **Seeds are set explicitly and recorded** for every stochastic operation: permutations,
  UMAP, t-SNE, PHATE, matched background sampling.
- Stages are **idempotent**: re-running with unchanged inputs and config produces
  byte-identical output, or skips with a logged reason.
- The `Makefile` encodes dependencies between stages, so `make figures` rebuilds only what
  is stale. If the pipeline grows beyond what a Makefile handles cleanly, move to
  Snakemake rather than to a hand-written driver script.
- Raw data is **read-only**. Nothing is ever written back into `data/raw/`.

## 6.5 Code style

House rules, several of which reflect explicit preferences for this project:

- **No stray `print` statements.** Use the `logging` module with a configured level. Stage
  progress belongs in a log, not in stdout noise.
- **No commentary comments.** Do not write `# loop over genes` above a loop over genes.
  Comments explain *why* a non-obvious choice was made, never *what* the next line does.
  Code that needs a what-comment needs a better name instead.
- **Docstrings on every public function**, numpydoc style, stating units and the shape of
  array arguments. For numerical functions, cite the source of the formula.
- **Type hints throughout.** Checked with `mypy` where practical.
- **Formatting and linting are automated**, not argued about: `ruff` for both, run in a
  pre-commit hook.
- **Small, pure functions for anything numerical.** A function that takes arrays and
  returns arrays can be tested; a function that reads a file, computes and plots cannot.
- **Prefer vectorised operations and explicit dtypes.** Avoid premature abstraction: no
  class hierarchy where a function will do, and no configuration framework beyond the
  single YAML file.

## 6.6 Testing

`pytest`, focused where correctness is actually at risk:

| Target | Test |
| --- | --- |
| FST estimator | Known analytic cases, plus agreement with published SGVP values |
| Variance partition | Synthetic data with a known injected population effect |
| Liftover and strand harmonisation | Hand-checked fixture of 20 variants including ambiguous A/T and C/G cases |
| Pseudobulk aggregation | Small synthetic AnnData where the correct sums are known by construction |
| Gene assignment | Variants placed at known window boundaries, including off-by-one cases |
| Permutation machinery | Under a true null, p-values are approximately uniform |

Fixtures are tiny and synthetic, committed to the repository, and run in seconds. The test
suite must never depend on downloaded data.

## 6.7 Version control

- Work on `main` for documentation; use short-lived branches for pipeline stages, merged
  when their acceptance check passes.
- **Small commits with specific messages.** Prefix with the stage where it helps, for
  example `stage5: use Hudson FST estimator and document the choice`.
- Commit messages describe *why*, not only *what*.
- **Commit messages must not include co-authorship or generated-by attribution lines.**
- Never commit data, figures, credentials or notebook output. The `.gitignore` enforces
  the first two; a pre-commit hook strips notebook output.
- Tag the commit used for the final submission.

## 6.8 Provenance and honesty

- `data/MANIFEST.tsv` records every input file with source URL, access date, SHA-256,
  genome build and licence.
- Every stage writes a short report to `data/interim/<stage>/report.md` containing counts
  in and out, parameters and the acceptance check result.
- Known failures and dropped records are reported as counts, never silently discarded. A
  pipeline that quietly loses 12 per cent of its variants is worse than one that crashes.
- Any number appearing in the report or the talk must be traceable to a stage output. None
  are typed by hand.
