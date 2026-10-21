# 5. Visualisation design and the Gestalt contract

## 5.1 The contract

Every figure in this project must answer one question before it is drawn:

> **What grouping is this image asking the viewer to believe, and is that grouping real?**

The Gestalt principles describe how the visual system imposes structure whether or not the
data warrants it. Proximity, similarity and continuity are therefore not decorative
guidelines here. They are the mechanism by which a figure makes a claim, and so they are
treated as part of the analysis, not part of the styling.

Three rules follow, and they apply to every figure without exception:

1. **Nothing is placed near anything else by accident.** Spatial adjacency is reserved for
   actual relatedness.
2. **Nothing shares a visual attribute by accident.** One colour means one thing across the
   entire project, in every figure, including supplementary ones.
3. **No continuous path is implied across discontinuous data.** Lines connect things that
   are genuinely ordered. Nothing else.

## 5.2 The three principles, operationalised

### Proximity

Elements placed near one another are perceived as belonging together.

| Applied to | Decision |
| --- | --- |
| Small multiples | Gaps within a group of panels are set to roughly half the gap between groups, so the grouping is read from spacing before any label is consulted |
| Legends | Direct labelling at the end of a series or beside a marked point, in preference to a separate key the eye must travel to and return from |
| Annotation | Case-study gene labels sit adjacent to their points with leader lines only where unavoidable |
| Axis titles | Close to their axis, never floating at a distance that invites association with the wrong panel |
| Composite figures | Panels about the same gene are adjacent; panels about the same modality are adjacent. These two groupings must not conflict, which constrains panel order and is decided deliberately |

### Similarity

Elements sharing colour, shape, size or texture are perceived as a class.

| Applied to | Decision |
| --- | --- |
| Population identity | One fixed colour per population group (Chinese, Malay, Indian), defined once in `src/sgdiv/plotting/style.py` and never redefined. If this colour varies between two figures, the project has lied to the viewer |
| Gene class | Marker **shape**, not colour, distinguishes CPIC level A genes from background. Shape and colour then encode independent variables and can be read simultaneously |
| Data role | Same line weight and dash pattern for the same role (observed, fitted, null) throughout |
| Deliberate dissimilarity | Null distributions and permutation envelopes are rendered in neutral grey with no fill saturation, so they are read as context rather than as a competing series |

### Continuity

The eye follows smooth, uninterrupted paths and completes interrupted ones.

| Applied to | Decision |
| --- | --- |
| Heatmap ordering | Rows and columns are **seriated** (optimal leaf ordering on the dendrogram), never left in input order. A heatmap in arbitrary order actively fights continuity and hides structure |
| Shared axes | Small multiples share scales and baselines so the eye can travel across panels without recalibrating |
| Genome tracks | Position on the genome is continuous and is drawn as such, with chromosome boundaries marked as interruptions rather than hidden |
| Ranked displays | Where an ordering exists, it is used, so the figure has a direction to read in |
| False continuity | Lines are never drawn between categorical levels. This is the most common way continuity is abused and it will not appear in this project |

## 5.3 The proximity problem in embeddings

This deserves its own section because it is where the course content and the design
doctrine meet, and because it is the most defensible intellectual point in the project.

UMAP and t-SNE produce layouts whose entire communicative force comes from proximity.
A viewer looking at a UMAP sees clusters because the Gestalt system groups nearby points,
and the viewer then attributes that grouping to the biology. But inter-cluster distances
in t-SNE are not meaningful, and both methods can manufacture apparent separation from
hyperparameter choices alone.

The project therefore adopts the following rules:

1. **No embedding appears in any figure without a stability diagnostic available.**
   Trustworthiness, continuity and k-nearest-neighbour overlap are computed and reported
   for every embedding shown (Stage 7).
2. **Apparent clusters that do not survive the seed sweep are not drawn as clusters.** No
   convex hulls, no shading, no cluster colouring applied to unstable structure.
3. **PCA appears alongside** any non-linear embedding in supplementary material, because
   PCA's axes mean something and provide a sober reference.
4. **A seed and hyperparameter grid is shown as a small multiple** in at least one figure,
   so the reader can see for themselves how much the picture moves.

Stated positively: the figure should be as honest about what is uncertain as it is
persuasive about what is solid. This is the argument to put to the lecturer, and it is the
one most likely to interest someone who teaches both perception and manifold learning.

## 5.4 Colour, type and format

Constraints, all enforced through a single matplotlib style module:

- **Restrained palette.** At most three categorical colours in any figure (the three
  population groups). Everything else is greyscale, with a single accent reserved for
  highlighting. Journal figures are not decorated; they are legible.
- **Perceptually uniform continuous maps only.** `viridis`, `cividis` or `magma` for
  sequential data; a balanced diverging map centred at a meaningful zero for signed data.
  No rainbow, no jet.
- **Colour-blind safe, and never load-bearing alone.** Any distinction carried by colour is
  redundantly encoded by shape, position or direct label. Verify with a deuteranopia
  simulation before submission.
- **Typography.** One sans-serif family. Axis labels and tick labels at fixed point sizes
  that remain legible at final print size. Every axis carries units.
- **Format.** Vector output (PDF for the report, SVG where a web version is needed) at
  explicit physical dimensions, typically single column 89 mm or double column 183 mm.
  Figures are never resized after export, because that silently rescales all text.
- **Data-ink.** No chartjunk, no 3D effects on 2D data, no gridlines heavier than the data,
  no boxes around legends.

## 5.5 Figure specifications

| ID | Figure | Carries |
| --- | --- | --- |
| F1 | Study design schematic | The three modalities, the shared population triad, and the explicit statement that cohorts differ |
| F2 | Allele frequency comparison panels | Pairwise population frequency density, pharmacogenes overlaid on background |
| F3 | Genome-wide divergence track | FST and iHS along the genome, with drug-gene positions annotated |
| F4 | **The divergence landscape** | Genetic divergence versus expression divergence, per gene, with CPIC level A marked by shape |
| F5 | Cell-type divergence heatmap | Gene by cell type variance partition coefficients, seriated |
| F6 | Locus case studies | For 3 or 4 loci: LD structure, population frequencies, and cell-type expression, vertically aligned |
| F7 | Embedding stability audit | Seed and hyperparameter grid, with preservation metrics |

### F4 in detail, since it is the centrepiece

- **Position** encodes the two divergence measures. Position is the most accurately
  perceived visual channel and is therefore spent on the quantities the argument rests on.
- **Shape** encodes CPIC level A membership. Shape, not colour, so that it reads as a
  categorical class without competing with any continuous encoding.
- **Grey density or contours** encode the background gene distribution, placed behind the
  highlighted points, so that similarity (uniform grey) marks the whole background as one
  class and the eye separates figure from ground immediately.
- **Quadrants** are defined by the permutation-derived thresholds on each axis, drawn as
  thin reference lines, and labelled in words at the plot edge rather than in a legend.
  The four quadrants are the finding, so they are named: high genetic and high expression
  divergence, high genetic only, high expression only, and neither.
- **Direct labels** on a small, pre-declared set of case-study genes. The set is chosen
  before the plot is made, so that labelling is not used to tell a story the statistics do
  not support.

## 5.6 Production rules

1. **One script per figure**, named for the figure ID, reading only from
   `data/processed/`.
2. **No manual editing after export.** If a figure needs a tweak, the tweak goes in the
   script. A figure that cannot be regenerated from a clean checkout does not go in the
   report.
3. **Every figure script records provenance**: git commit SHA and config hash written to a
   sidecar file next to the output.
4. **Caption written in the same commit as the figure.** Captions live in
   `figures/captions/` as markdown, so the figure and its claim stay together.
5. **Style is centralised.** No `rcParams` set inside a figure script. All of it lives in
   the style module, so that the similarity contract in 5.2 cannot be broken by accident.
