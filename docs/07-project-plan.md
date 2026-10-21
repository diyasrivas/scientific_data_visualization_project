# 7. Project plan

## 7.1 Milestones

| # | Milestone | Depends on | Definition of done |
| --- | --- | --- | --- |
| M0 | Feasibility settled | none | Stage 0 link audit complete; go or fallback decision recorded |
| M1 | Proposal sent and accepted | M0 | Prof. Casiraghi has replied; scope agreed |
| M2 | Data in hand | M0 | SGVP and AIDA downloaded, checksummed, manifested |
| M3 | Both divergence axes computed | M2 | Stages 5 and 6 pass acceptance checks, including the published-FST and housekeeping controls |
| M4 | Integration complete | M3 | Joint feature matrix built; embeddings with stability diagnostics |
| M5 | Results frozen | M4 | Stage 8 run; sensitivity analyses complete; no further analysis changes |
| M6 | Figures final | M5 | All seven figures regenerate from a clean checkout |
| M7 | Report and talk delivered | M6 | Submitted |

Order matters more than dates. M0 before M1 in particular: do not promise the professor a
dataset that has not been confirmed to download.

## 7.2 Suggested sequencing

Roughly a semester, assuming part-time effort:

- **Weeks 1 to 2.** Stage 0 and the proposal. Set up the repository skeleton, environment
  and CI while waiting for a reply.
- **Weeks 3 to 4.** Stages 1 to 3. Acquisition and annotation. Get the FST reproduction
  check passing, because it validates a large part of the pipeline at once.
- **Weeks 5 to 6.** Stage 4. Harmonisation. Budget more time than feels necessary; this is
  the stage that silently eats projects.
- **Weeks 7 to 9.** Stages 5 and 6. The two divergence axes, with their controls.
- **Weeks 10 to 11.** Stages 7 and 8. Integration, embeddings, testing.
- **Weeks 12 to 13.** Stage 9. Figures, iterated properly rather than rushed.
- **Week 14.** Report and talk.

## 7.3 Risk register

| Risk | Likelihood | Impact | Mitigation |
| --- | --- | --- | --- |
| SGVP files unreachable | Medium | High | Stage 0 gate in week 1; gnomAD fallback documented; contact authors early |
| Liftover loses many variants | Medium | Medium | Explicit Stage 4 report with a 95 per cent threshold; gnomAD concordance check |
| Too few Singaporean donors per group in AIDA for stable mixed models | Medium | High | Check donor counts immediately at Stage 2; if thin, collapse to coarser cell types and say so |
| AIDA object too large to handle | Low | Medium | Backed-mode streaming and early pseudobulk reduction |
| Gene identifier mismatches corrupt joins | Medium | High | Ensembl IDs only, never symbols; spot-check list of ten pharmacogenes |
| Result is a clean null | Medium | Low | H2 is pre-registered as a null; the figure and the talk work either way |
| Scope creep into the lecturer's own methods | Medium | Medium | Patient similarity networks and knowledge graphs are explicitly out of scope |
| Over-interpreting an ecological comparison | Medium | High | Caveat stated in the proposal, the report, the talk and on figure F1 |

## 7.4 Scope boundaries

**In scope.** Population genetics statistics on SGVP; pseudobulk expression divergence on
AIDA Singaporean donors; clinical annotation overlay; dimensionality reduction with
stability auditing; seven figures; report and talk.

**Out of scope.** Individual-level genotype to expression inference; new eQTL calling;
patient similarity network methods; knowledge graph construction; anything requiring
SG10K_Health or other controlled-access data; spatial or subzone-level analysis.

The out-of-scope list is as important as the in-scope list. It is what keeps the project
finishable, and what keeps it from wandering into territory where the lecturer has a
decade of advantage.

## 7.5 Deliverables

1. This repository, reproducible from a clean checkout.
2. Seven figures as vector files with captions.
3. A written report.
4. A presentation built around figure F4.
5. `docs/stage0-link-audit.md` and the per-stage reports, as the audit trail.
