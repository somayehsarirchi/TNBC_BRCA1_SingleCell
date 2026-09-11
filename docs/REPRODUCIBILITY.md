# Reproducibility

The project is organized as a sequential, stage-based analysis
pipeline.

Each major stage is represented by an R script and a corresponding
results directory.

The workflow emphasizes:

- explicit analytical parameters
- stage-wise input/output relationships
- intermediate validation
- reproducible metadata
- complete ranked gene lists where appropriate
- explicit gene universes for enrichment
- separation of discovery and interpretation
- frozen trajectory validation
- final export and reload validation

Large biological objects are intentionally excluded from GitHub.
They remain in the local analysis project and are documented by
the corresponding object-stage directories.

The final Stage 17G output is an archival and validation stage.
It does not introduce a new biological analysis.
