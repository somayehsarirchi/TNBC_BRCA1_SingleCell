# TNBC BRCA1 Single-Cell Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22723081.svg)](https://doi.org/10.5281/zenodo.22723081)
[![R](https://img.shields.io/badge/R-%3E%3D4.5.2-276DC3?logo=r\&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5.5.0-success)](https://satijalab.org/seurat/)
[![Workflow](https://img.shields.io/badge/Pipeline-17%20Stages%20Reproducible-orange)](#analytical-workflow)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Executive Summary

This repository provides a stage-organized single-cell RNA-sequencing (scRNA-seq) workflow for investigating cellular heterogeneity, epithelial cell states, copy-number variation, pathway activity, cell–cell communication, and transcriptional dynamics in triple-negative breast cancer (TNBC), using the publicly available GSE161529 dataset.

The analysis compares the study-defined **BRCA1_tumour** and **TotalCell** groups while progressively resolving cellular composition, epithelial identity, and CNV-associated epithelial states before downstream molecular analyses.

The workflow integrates:

* Cell-type annotation and epithelial compartment refinement
* CopyKAT-based inference of epithelial CNV states
* Differential expression and multi-database pathway analysis
* CellChat-based analysis of intercellular communication
* Monocle3 trajectory inference and pseudotime analysis
* Identification and clustering of genes with distinct expression dynamics along pseudotime

> **Important:** `Tumor_like`, `Normal_like`, `Diploid`, and `Aneuploid` are analytical labels derived from the computational workflow. In particular, **aneuploidy should not be interpreted as experimental proof of malignancy**, and pseudotime should not be interpreted as direct chronological time or lineage tracing.

---

## Core Highlights

### CNV-Associated Epithelial State Characterization

CopyKAT was applied to the purified epithelial compartment together with annotated non-epithelial reference populations to infer large-scale copy-number states.

Epithelial cells were classified into:

* **Aneuploid**
* **Diploid**
* **Unknown**

The aneuploid epithelial population was subsequently represented by the analytical label `Tumor_like` for downstream tumor-associated molecular characterization.

This label reflects the computational classification used in the workflow and does not by itself establish malignant identity.

### Cross-Database Pathway Integration

Differentially expressed genes were evaluated across complementary gene-set resources, including:

* Gene Ontology Biological Process (GO BP)
* Reactome
* KEGG
* MSigDB Hallmark

Recurring biological themes were examined across these resources rather than relying on a single pathway database.

### Ecosystem-Level Cell–Cell Communication

CellChat was used to characterize predicted intercellular communication across epithelial and non-epithelial cellular populations, including tumor-associated and diploid epithelial states.

The analysis examines incoming and outgoing signaling patterns and provides an ecosystem-level view of predicted cellular communication within the TNBC microenvironment.

### Pseudotemporal Gene-Expression Programs

Monocle3 trajectory inference was used to construct an epithelial transcriptional manifold and assign relative pseudotime.

Trajectory-associated genes were identified using graph-based spatial autocorrelation (Moran's I), followed by temporal clustering to identify distinct patterns of gene-expression dynamics along pseudotime.

The resulting clusters represent **gene-expression programs**, not independent cell types or experimentally defined cellular states.

---

## Representative Results

| Figure       | Description                                               |
| ------------ | --------------------------------------------------------- |
| **Figure 1** | Cell-Type Annotation Landscape                            |
| **Figure 2** | CopyKAT-Inferred CNV States in the Epithelial Compartment |
| **Figure 3** | Intercellular Communication Weights                       |
| **Figure 4** | Epithelial Trajectory & Pseudotime                        |

### Figure 1 — Cell-Type Annotation Landscape

![Cell-Type Annotation](docs/figures/UMAP_Manual_CellTypes_Split.png)

### Figure 2 — CopyKAT-Inferred CNV States

![CopyKAT CNV State](docs/figures/CopyKAT_CNV_Status_UMAP.png)

### Figure 3 — Intercellular Communication

![CellChat Network](docs/figures/CellChat_Overall_Network_Weight.png)

### Figure 4 — Epithelial Trajectory and Pseudotime

![Monocle3 Trajectory](docs/figures/Trajectory_Pseudotime_Final.png)

---

# Study Design

The analytical workflow progressively resolves cellular composition, epithelial identity, CNV-associated states, molecular programs, intercellular communication, and transcriptional dynamics.

```text
All Cells
   │
   ▼
QC → Normalization → Integration → Clustering
   │
   ▼
Cell-Type Annotation
   │
   ▼
Epithelial Compartment Refinement
   │
   ▼
CopyKAT CNV Inference
   │
   ├───────────────┬──────────────────────┐
   ▼               ▼                      ▼
Diploid         Aneuploid              Unknown
                   │
                   ▼
             Tumor-associated
             epithelial subset
                   │
          ┌────────┼─────────────┐
          ▼        ▼             ▼
         DE       GSEA       Monocle3
                              │
                              ▼
                         Pseudotime
                              │
                              ▼
                     Gene Dynamics

CellChat is analyzed as a separate ecosystem-level branch using
the defined epithelial and non-epithelial cellular populations.
```

### Analytical Rationale

To reduce cellular-composition confounding, differential expression and pathway analyses were performed within the defined epithelial/CNV-associated analytical subset rather than treating all annotated cell populations as a single transcriptional compartment.

CellChat was retained as a separate branch because its purpose is to characterize predicted communication across the broader cellular ecosystem.

---

# Analytical Workflow

The repository is organized into **17 analytical stages**, with each stage corresponding to a defined component of the computational workflow.

```text
Stage 01 → Data Loading & Quality Control
Stage 02 → Normalization & Integration
Stage 03 → PCA, Clustering & UMAP
Stage 04 → Cell-Type Annotation
Stage 05 → Epithelial Compartment Analysis
Stage 06 → CopyKAT Preparation
Stage 07 → CopyKAT Batch Inference
Stage 08 → CopyKAT Integration & CNV State Assignment
Stage 09 → Tumor-Associated Epithelial Characterization
Stage 10 → Differential Expression
Stage 11 → GSEA
Stage 12 → GSEA Integration
Stage 13 → Representative Pathways & Biological Interpretation
Stage 14 → CellChat
Stage 15 → Monocle3 Trajectory
Stage 16 → Pseudotime Gene Dynamics & Gene Discovery
Stage 17 → Temporal Gene-Expression Dynamics
```

---

## Stage-by-Stage Overview

### Stage 01 — Data Loading and Quality Control

Initial loading and quality-control procedures were applied to the GSE161529 single-cell dataset.

The analysis begins with eight samples:

* GSM4909281
* GSM4909282
* GSM4909283
* GSM4909284
* GSM4909285
* GSM4909286
* GSM4909287
* GSM4909288

Cells passing the defined quality-control criteria were retained for downstream analysis.

---

### Stage 02 — Normalization and Integration

The retained single-cell datasets were normalized and integrated to reduce technical variation between samples while preserving biologically relevant transcriptional structure.

The workflow uses Seurat v5 and Harmony-based integration.

---

### Stage 03 — PCA, Clustering and UMAP

Principal component analysis, graph-based clustering, and UMAP dimensionality reduction were used to characterize the major transcriptional structure of the dataset.

---

### Stage 04 — Cell-Type Annotation

Cell populations were annotated using marker-based biological interpretation together with reference-based annotation.

The resulting cellular landscape provides the basis for subsequent compartment-specific analyses.

---

### Stage 05 — Epithelial Compartment Analysis

The epithelial compartment was refined from the annotated cellular landscape to support downstream CNV inference and epithelial-state characterization.

---

### Stage 06 — CopyKAT Preparation

The epithelial compartment was prepared for CopyKAT-based CNV inference.

Annotated non-epithelial populations were retained as reference cells to support computational inference of epithelial copy-number states.

---

### Stage 07 — CopyKAT Batch Inference

CopyKAT was run using the following principal parameters:

```text
batch size: 2500
min_genes_per_cell: 200
min_cells_per_gene: 3
min_normal_cells_per_batch: 20
genome: hg20
id.type: S
distance: pearson
n.cores: 4
seed: 1234
ngene.chr: 5
win.size: 25
KS.cut: 0.1
```

---

### Stage 08 — CopyKAT Integration and CNV State Assignment

CopyKAT outputs were integrated across batches and mapped back to the epithelial compartment.

The resulting analytical CNV categories were:

```text
Aneuploid
Diploid
Unknown
```

For downstream tumor-associated epithelial characterization:

```text
Aneuploid → Tumor_like
Diploid   → Normal_like
Unknown   → Unknown
```

Here, `Tumor_like` and `Normal_like` are workflow-specific analytical labels.

In particular:

> **Diploid ≠ experimentally proven normal**

and

> **Aneuploid ≠ experimentally proven malignant**

The CNV state represents a computationally inferred genomic characteristic.

---

### Stage 09 — Tumor-Associated Epithelial Characterization

The `Tumor_like` analytical subset was isolated for tumor-associated molecular characterization.

This subset corresponds to epithelial cells classified as aneuploid by the preceding CopyKAT workflow.

The analysis therefore uses `Tumor_like` as a computational subset label rather than as an independent experimental annotation of malignancy.

---

### Stage 10 — Differential Expression

Differential expression analysis was performed within the defined epithelial/CNV-associated analytical framework.

A total of **6,115 genes** were tested:

* **5,162** genes were significant
* **2,685** showed higher expression in the BRCA1_tumour group
* **2,477** showed higher expression in the TotalCell group

The results provide the basis for downstream pathway-level interpretation.

---

### Stage 11 — Gene Set Enrichment Analysis

Gene set enrichment analysis was performed using complementary pathway resources, including:

* GO Biological Process
* Reactome
* KEGG
* MSigDB Hallmark

A total of **5,846 genes** were mapped to the pathway-analysis framework, corresponding to approximately **95.6% mapping coverage**.

---

### Stage 12 — GSEA Integration

Enrichment results from the different gene-set resources were integrated to identify recurring biological themes across databases.

This cross-resource comparison was used to distinguish broadly supported biological patterns from findings that were specific to a single pathway resource.

---

### Stage 13 — Representative Pathways and Biological Interpretation

Representative pathways were selected for visualization and biological interpretation based on statistical significance, recurrence across pathway resources, and relevance to the epithelial/TNBC analytical context.

The selected pathways are intended as interpretable representatives of broader enrichment patterns rather than as an exhaustive list of all significant terms.

---

### Stage 14 — CellChat

CellChat was used as an independent ecosystem-level branch of the workflow.

The analysis included:

* Tumor-associated epithelial cells
* Diploid epithelial cells
* Annotated non-epithelial populations

Unknown epithelial cells were excluded from the defined CellChat analysis.

The analysis evaluates predicted ligand–receptor-mediated communication, including:

* Overall communication networks
* Incoming signaling
* Outgoing signaling
* Relative pathway contributions
* Communication patterns involving epithelial populations

Because CellChat infers communication from expression-based ligand–receptor relationships, these results represent **predicted signaling interactions** rather than experimentally validated cell–cell communication events.

---

### Stage 15 — Monocle3 Trajectory

Monocle3 was used to infer a transcriptional trajectory within the epithelial compartment.

Graph learning was performed on the epithelial transcriptional manifold, with programmatic root-node assignment anchored to a **diploid epithelial reference region**.

Pseudotime therefore represents a relative ordering of transcriptional states along the inferred manifold.

Where diploid and aneuploid cells occupy different regions of the trajectory, the result can be interpreted as an association between transcriptional trajectory position and inferred CNV state.

Importantly:

> Monocle3 does not demonstrate that diploid cells physically become aneuploid cells.

The trajectory is an inferred transcriptional structure rather than experimental lineage tracing or direct longitudinal observation.

---

### Stage 16 — Pseudotime Gene Dynamics and Gene Discovery

Genes associated with trajectory position were identified using graph-based spatial autocorrelation statistics, including Moran's I.

A total of:

* **24,865 genes** were tested
* **14,356 genes** were significant
* **2,908 genes** were retained as strong trajectory-associated candidates

These genes provide a focused set of candidates for investigating transcriptional programs associated with progression along the inferred epithelial trajectory.

---

### Stage 17 — Temporal Gene-Expression Dynamics

Trajectory-associated genes were evaluated across **20 pseudotime bins** to characterize patterns of gene-expression change along the inferred epithelial trajectory.

Temporal clustering identified two distinct gene-expression programs:

| Temporal Cluster | Gene Count | Pattern                    | Peak Bin | Trough Bin | Mean Centroid Correlation | Mean Silhouette | Biological Description                                                    |
| ---------------- | ---------: | -------------------------- | -------: | ---------: | ------------------------: | --------------: | ------------------------------------------------------------------------- |
| **Cluster 1**    |         53 | Early-peaking / decreasing |        3 |          8 |                     0.817 |           0.519 | Early-peaking immune/antigen-presentation-associated temporal program     |
| **Cluster 2**    |        516 | Late-increasing            |       19 |          6 |                     0.893 |           0.632 | Broad late-emerging epithelial/cellular-state-associated temporal program |

Cluster 1 showed higher mean standardized expression in early pseudotime (**Early Mean Z = 1.531**) and lower mean standardized expression in late pseudotime (**Late Mean Z = −0.397**), with its peak at pseudotime bin 3.

Cluster 2 showed lower mean standardized expression in early pseudotime (**Early Mean Z = −0.948**) and higher mean standardized expression in late pseudotime (**Late Mean Z = 1.010**), with its peak at pseudotime bin 19.

The relatively high centroid correlations and positive silhouette values indicate coherent within-cluster temporal expression patterns.

These clusters represent **temporal gene-expression programs along the inferred pseudotime trajectory**. They do not represent independent cell subtypes, experimentally defined cellular states, or evidence of physical cell-type conversion.

The biological descriptions summarize the dominant functional interpretation associated with each temporal program and should therefore be understood as **trajectory-associated transcriptional patterns**, rather than direct evidence of chronological cellular transformation.

---

# Integrated Biological Interpretation

The integrated analysis supports a multi-layered view of epithelial heterogeneity in the dataset:

```text
Cellular heterogeneity
        ↓
Epithelial compartment
        ↓
Inferred CNV-associated states
        ↓
Tumor-associated epithelial molecular programs
        ↓
Pathway-level biological interpretation
        ↓
Cell–cell communication within the cellular ecosystem
        ↓
Transcriptional trajectory
        ↓
Pseudotime-associated genes
        ↓
Distinct gene-expression dynamics
```

### CNV-Associated Epithelial States

CopyKAT provides a computational distinction between epithelial cells with inferred diploid and aneuploid CNV profiles.

The aneuploid epithelial population is subsequently used as a tumor-associated analytical subset.

### Transcriptional Programs

Differential expression and pathway analyses identify transcriptional programs associated with the defined epithelial comparison.

Cross-database integration provides complementary evidence for recurring biological themes.

### Cellular Communication

CellChat extends the analysis beyond intracellular transcriptional programs by examining predicted communication between epithelial and non-epithelial populations.

### Transcriptional Dynamics

Monocle3 provides a low-dimensional representation of epithelial transcriptional organization and a relative pseudotemporal ordering.

Trajectory-associated genes and their expression dynamics provide candidate molecular programs associated with different regions of the inferred transcriptional manifold.

---

# Important Interpretation Safeguards

1. **CopyKAT CNV calls are computational inferences.**
   They should not be treated as equivalent to experimentally measured genomic copy-number profiles without independent validation.

2. **Aneuploidy is not synonymous with malignancy.**
   The `Tumor_like` label is an analytical workflow label for the aneuploid epithelial subset and should not be interpreted as definitive experimental proof of tumor identity.

3. **Diploid is not synonymous with normal.**
   A diploid CNV profile does not establish that a cell is biologically normal.

4. **CellChat results are predictions.**
   Ligand–receptor communication inferred by CellChat represents computationally predicted interactions and requires experimental validation.

5. **Pseudotime is not chronological time.**
   Pseudotime represents a relative ordering of transcriptional states along an inferred trajectory and should not be interpreted as direct chronological time or experimental lineage tracing.

6. **Trajectory does not establish physical cell conversion.**
   Even when diploid and aneuploid cells occupy distinct regions of an inferred trajectory, the analysis alone cannot demonstrate that one population physically transforms into the other.

7. **Temporal clusters are gene-expression programs.**
   Clusters identified from pseudotemporal gene dynamics describe groups of genes with similar expression patterns along pseudotime; they are not automatically cell subtypes or biological states.

8. **Computational analyses require biological validation.**
   The results are intended to generate and prioritize biological hypotheses and candidate molecular programs for further investigation.

---

# Reproducibility

The repository is organized as a stage-based computational workflow with explicit intermediate outputs, figures, logs, documentation, and environment specifications.

The computational environment used for the analysis includes:

```text
R 4.5.2
Bioconductor 3.22
Seurat 5.5.0
SingleR 2.12.0
CopyKAT 1.1.0
CellChat 1.6.1
Monocle3 1.4.27
DESeq2 1.50.2
edgeR 4.8.2
clusterProfiler 4.18.4
survival 3.8.6
```

The complete environment specification is provided in:

```text
environment/R_environment.txt
```

Raw GEO data are not redistributed in this repository and should be obtained separately according to the availability and terms of the original dataset.

---

# Repository Structure

```text
TNBC_BRCA1_SingleCell/
│
├── scripts/
│   ├── Stage_01_...
│   ├── Stage_02_...
│   ├── ...
│   └── Stage_17_...
│
├── results/
│   ├── 01_Data_Loading_QC/
│   ├── 02_Normalization_Integration/
│   ├── 03_PCA_Clustering_UMAP/
│   ├── 04_CellType_Annotation/
│   ├── 05_Epithelial_Compartment/
│   ├── 06_CopyKAT_Preparation/
│   ├── 07_CopyKAT_Batch_Inference/
│   ├── 08_CopyKAT_Integration/
│   ├── 09_TumorLike_Epithelial/
│   ├── 10_Differential_Expression/
│   ├── 11_GSEA/
│   ├── 12_GSEA_Integration/
│   ├── 13_Representative_Pathways/
│   ├── 14_CellChat/
│   ├── 15_Monocle3_Trajectory/
│   ├── 16_Pseudotime_Gene_Dynamics/
│   └── 17_Temporal_Gene_Dynamics/
│
├── docs/
│   ├── figures/
│   └── ...
│
├── environment/
│   └── R_environment.txt
│
├── objects/
│   └── ...
│
├── LICENSE
└── README.md
```

---

# Quick Start

## Clone the Repository

```bash
git clone https://github.com/somayehsarirchi/TNBC_BRCA1_SingleCell.git
cd TNBC_BRCA1_SingleCell
```

## Environment

Install the required R and Bioconductor packages according to:

```text
environment/R_environment.txt
```

## Workflow

The analysis is designed to be executed sequentially from Stage 01 through Stage 17.

Each stage contains its own script, outputs, and documentation where applicable.

Because raw GEO data are not included in the repository, the required input data must be obtained separately before executing the workflow.

---

# Dataset

**GEO accession:** GSE161529

**Samples analyzed:**

```text
GSM4909281
GSM4909282
GSM4909283
GSM4909284
GSM4909285
GSM4909286
GSM4909287
GSM4909288
```

The study-defined analytical groups used in this repository are:

```text
BRCA1_tumour
TotalCell
```

These labels are retained as defined by the dataset/workflow and should not be replaced by stronger biological labels unless independently supported by the original study metadata.

---

# Citation

If you use this repository, please cite the associated Zenodo record:

**Sarirchi, Somayeh. TNBC BRCA1 Single-Cell Analysis. Zenodo.**

DOI:

```text
10.5281/zenodo.22723081
```

BibTeX:

```bibtex
@software{sarirchi_tnbc_brca1_singlecell,
  author       = {Sarirchi, Somayeh},
  title        = {TNBC BRCA1 Single-Cell Analysis},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.22723081}
}
```

---

# License

This project is released under the **MIT License**.

See the `LICENSE` file for the complete license text.

---

# Author

**Somayeh Sarirchi, PhD**

Computational Cancer Biology & Single-Cell Genomics

Research focus: TNBC cellular heterogeneity, tumor-associated epithelial states, and single-cell transcriptomics.
