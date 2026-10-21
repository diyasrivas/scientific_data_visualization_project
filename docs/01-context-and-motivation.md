# 1. Context and motivation

## 1.1 The course

**Scientific Visualization**, 6 ECTS, sector INFO-01/A, taught by Prof. Elena Casiraghi
in the MSc in Computer Science (LM-18) at Università degli Studi di Milano. In the
programme manifesto the course sits in the **Data Analytics** area of the *Artificial
Intelligence, Data Analytics and Big Data* pillar, alongside Network Science, Algorithms
for Massive Datasets and Advanced Multivariate Statistics.

There is no published syllabus. From the lecture content so far the course spans:

- principles of data visualization and human visual perception,
- dimensionality reduction and manifold learning (t-SNE, UMAP and relatives),
- omics and biomedical data as application domain.

The lecturer has stressed the **Gestalt principles of proximity, similarity and
continuity**. This project treats those three principles as a design contract rather than
as a topic to be mentioned once, which is set out in
[`05-visualisation-and-gestalt.md`](05-visualisation-and-gestalt.md).

Assessment is a project: select data, run a full analysis pipeline, produce exceptional
visualizations, and present.

## 1.2 Reading the supervisor

Prof. Casiraghi is an associate professor in the Department of Computer Science "Giovanni
degli Antoni", scientific lead of AnacletoLab, and an affiliate researcher at Lawrence
Berkeley National Laboratory. The relevant parts of her record, taken from her CV and
publication list:

**Her core expertise. Use her methods, do not try to extend them.**

- Multimodal data integration through patient similarity networks, including the case of
  entirely missing data sources (`miss-SNF`, Bioinformatics 2025).
- Missing-data imputation and comparison of imputation techniques, as group leader of the
  "Feature Selection and Missing Data" working group in the NIH N3C COVID enclave.
- Biomedical knowledge graphs, with the Monarch Initiative and the RNA-KG project.
- Intrinsic dimension estimation. She co-proposed **DANCo** (Pattern Recognition, 2014),
  a state-of-the-art estimator, and her 2025 paper in *Artificial Intelligence in
  Medicine* uses intrinsic dimension to guide dimensionality reduction and data fusion in
  multi-omics data.

**Adjacent interests. This is the sweet spot.**

- Environmental and nutritional exposure linked to disease, via a knowledge graph and
  random forests (*International Journal of Medical Informatics*, 2024).
- Explainable and interpretable models intended to reach clinical practice.
- Teaching Scientific Visualization and a separate course on models of perception.

**Outside her published work. This is where the student can hold the expertise.**

- Population genetics statistics: FST, iHS, XP-EHH, varLD, linkage disequilibrium
  structure.
- Pharmacogenomics and clinical variant actionability (PharmGKB, CPIC).
- Singapore-specific population health and genomic resources.

The project is deliberately positioned so that the lecturer is the expert on one half
(single-cell data, dimensionality reduction, integration) and the student is the expert
on the other half (population genetics, pharmacogenomics, the Singaporean context). That
makes the viva a conversation rather than an examination, and it gives her something she
has a genuine reason to be curious about.

## 1.3 Why this topic and not another

Personal constraints that shaped the choice:

- Interest in medical and healthcare data, and in genetics specifically.
- Singaporean context, using openly available Singaporean data.
- A multimodal element.
- Preference for a national or gene-level question rather than one requiring familiarity
  with Singapore's planning areas and subzones, which rules out the geospatial variants
  of this project.
- Appetite for something novel enough to produce an unexpected result.

## 1.4 The intellectual hook

Population genetics tells us *that* Singapore's three groups differ. It is much less
clear *whether it matters*. A variant can be 40 percentage points apart in frequency
between two groups and do nothing at all; another can differ by five points and change
how a patient metabolises clopidogrel.

The project makes that distinction visible. It puts genetic divergence on one axis and
functional divergence (population differences in gene expression, resolved by immune cell
type) on the other, then marks which genes carry clinically actionable drug guidelines.
The resulting scatter is a map of where population genetic difference is and is not
biologically consequential.

That is a visualization-first framing: the argument of the project *is* the figure.
