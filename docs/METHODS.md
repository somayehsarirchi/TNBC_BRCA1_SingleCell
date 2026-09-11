# Methods

## Overview

The project implements a modular single-cell RNA-seq analysis
workflow for TNBC using GEO dataset GSE161529.

## Analytical stages

1. Data loading and quality control
2. Normalization and data integration
3. PCA, clustering, and UMAP
4. Cell-type annotation
5. Epithelial compartment refinement
6. CopyKAT preparation
7. CopyKAT batch inference
8. CopyKAT integration and CNV classification
9. Tumor-like epithelial characterization
10. Differential expression
11. Gene-set enrichment analysis
12. Integrated pathway prioritization
13. Biological interpretation
14. CellChat communication analysis
15. Monocle3 trajectory analysis
16. Pseudotime-associated gene analysis
17. Temporal gene dynamics

## Scientific interpretation safeguards

CopyKAT-defined Tumor_like cells represent an operational
computational classification and are not treated as definitive
experimental proof of malignancy.

Monocle3 pseudotime represents a relative transcriptional ordering
within the analyzed trajectory and is not interpreted as direct
chronological time.

Trajectory structure is not interpreted as experimentally validated
lineage, developmental direction, irreversible transition, or
causal mechanism.

Temporal gene clusters describe gene-expression dynamics rather
than cellular subtypes.

Pathway enrichment represents statistical overrepresentation and
is not by itself evidence of pathway activation or flux.
