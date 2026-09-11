# Results

Analysis outputs are organized by analytical stage.

## Stage structure

- `01_Data_Loading_QC/`
- `02_Normalization_Integration/`
- `03_PCA_Clustering_UMAP/`
- `04_CellType_Annotation/`
- `05_Epithelial_Compartment/`
- `06_CopyKAT_Preparation/`
- `07_CopyKAT_Batch_Inference/`
- `08_CopyKAT_Integration/`
- `09_TumorLike_Epithelial_Characterization/`
- `10_Differential_Expression/`
- `11_GSEA/`
- `12_GSEA_Integration/`
- `13_GSEA_Biological_Interpretation/`
- `14_CellChat/`
- `15_Monocle3_Trajectory/`
- `16_Pseudotime_Gene_Dynamics/`
- `17_Temporal_Gene_Dynamics/`

Each stage contains curated lightweight outputs such as
tables, figures, summaries, and validation reports.

Large RDS objects and computational intermediates are excluded.
