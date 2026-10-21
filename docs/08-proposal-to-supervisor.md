# 8. Proposal to the supervisor

Draft text for the proposal to Prof. Elena Casiraghi. Send **after** Stage 0 has
confirmed data availability, so that nothing promised here turns out to be unobtainable.

## 8.1 Positioning notes

Three things this proposal is trying to do, which should survive any rewriting:

1. **Invite her in as the expert on one half.** The intrinsic dimension question and the
   embedding stability question are posed as questions to her, not as claims. That turns
   the strongest potential line of hard questioning into a collaboration.
2. **State the limitation before she finds it.** The ecological caveat (SGVP and AIDA are
   different cohorts) is raised explicitly. Given that her own research is on missing data
   and multimodal integration, she will see it within thirty seconds. Naming it first is
   the difference between looking careful and looking caught out.
3. **Make the visualization the deliverable.** This is a Scientific Visualization course.
   The proposal should read as a project whose result is a figure, not as a bioinformatics
   project with figures attached.

Keep it to roughly one page. Lecturers read proposals in the gaps between other work.

## 8.2 Draft email

> **Subject:** Scientific Visualization project proposal: population divergence and immune
> cell expression in Singapore
>
> Dear Professor Casiraghi,
>
> I would like to propose the following project for the Scientific Visualization course.
>
> **Question.** Singapore's resident population comprises three groups (Chinese, Malay and
> Indian) with distinct genetic ancestries. I would like to ask which of the loci that
> differ most between these groups also show differences in gene expression at the level
> of specific immune cell types, and whether those loci are enriched for genes with
> clinical pharmacogenomic guidelines. Put simply: where is population genetic difference
> actually consequential, and where is it not?
>
> **Data.** Three open, complementary modalities covering the same three population groups:
>
> - the Singapore Genome Variation Project (1.6 million SNPs in 268 individuals, with
>   allele frequencies, linkage disequilibrium, varLD and selection scores, plus a set of
>   3,972 variants in genes of drug disposition and response);
> - the Asian Immune Diversity Atlas (Cell, 2025: 1.26 million immune cells from 619
>   donors, including Singaporean Chinese, Malay and Indian donors);
> - clinical annotation from PharmGKB and CPIC.
>
> **Approach.** I would summarise each gene on two comparable axes: FST, the proportion of
> genetic variance attributable to population structure, and a variance partition
> coefficient computed from donor-level pseudobulk expression, which is the proportion of
> expression variance attributable to population group after adjusting for age, sex and
> batch. Because the two measures are constructed to mean the same thing on different
> modalities, they can be placed on a single scatter plot. The central deliverable of the
> project is that figure: a map of genes in which the four quadrants distinguish loci that
> diverge genetically and functionally, genetically only, functionally only, or neither,
> with clinically actionable genes marked.
>
> **Visualization focus.** I want to treat the Gestalt principles of proximity, similarity
> and continuity as a design contract rather than as decoration: one fixed colour per
> population group across every figure, panel spacing that encodes grouping before any
> label is read, and seriated rather than arbitrarily ordered heatmaps. I am also
> interested in a point where the course material and the honesty of a figure meet
> directly. UMAP and t-SNE communicate almost entirely through proximity, and a viewer
> reads nearness as group membership whether or not it is stable. I would therefore like
> to run a seed and hyperparameter sweep with preservation metrics (trustworthiness,
> continuity and nearest-neighbour overlap) and draw as a cluster only what survives it.
>
> **A limitation I want to be explicit about.** The two cohorts are different sets of
> people, so no individual-level genotype to expression link is possible. Every comparison
> is at the level of the population group, and I plan to state this on the design figure
> itself rather than in a footnote.
>
> **Two questions for you.** First, I would like to characterise the intrinsic dimension
> of the gene-level divergence feature space, and I know this is an area you have worked
> on directly. Are there estimators you would recommend, and how would you suggest
> handling small-sample-size effects here? Second, would you consider the ecological
> framing above acceptable for this course, or would you prefer I restrict the scope to a
> single modality?
>
> I have set up a repository with the full workflow, data provenance plan and figure
> specifications, and I am happy to send it or to discuss in office hours.
>
> Thank you for your time.
>
> Kind regards,
> Diya Srivastava

## 8.3 If she asks for something smaller

Two pre-planned reductions, in order of preference:

1. **Drop AIDA.** Build the pharmacogenomic divergence atlas from SGVP alone: allele
   frequency divergence, LD variation and selection signals across PharmGKB and CPIC
   genes, with the clinical layer retained. Still genetic, still clinical, still a strong
   visual deliverable, and roughly half the engineering.
2. **Drop the clinical layer.** Keep SGVP and AIDA and present the divergence landscape
   without actionability annotation, as a purely descriptive integration.

Do not reduce by dropping the Gestalt and embedding-honesty component. That is the part
that makes it a visualization project rather than a bioinformatics one.
