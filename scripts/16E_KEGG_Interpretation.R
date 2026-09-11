# ============================================================
# TNBC Single-Cell Analysis Pipeline
# Stage 16E: KEGG Pathway Interpretation and Prioritization
#
# Purpose:
#   Interpretation and prioritization of the significant KEGG
#   pathways identified in Stage 16D.
#
# Frozen upstream inputs:
#   Stage 15  -> Frozen Monocle3 trajectory
#   Stage 16A -> Frozen analysis population
#   Stage 16B -> graph_test trajectory-associated genes
#   Stage 16C -> Gene characterization
#   Stage 16D -> KEGG enrichment
#
# Important:
#   - Stage 15 remains completely frozen.
#   - Stage 16A/16B/16C/16D remain unchanged.
#   - KEGG enrichment is NOT rerun.
#   - No trajectory reconstruction is performed.
#   - No pseudotime recalculation is performed.
#   - No graph_test is rerun.
#   - This stage performs only pathway cleaning, ranking,
#     prioritization and reporting.
#   - The shared Stage 16 directory architecture is preserved.
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace(
  "dplyr",
  quietly = TRUE
)) {
  
  stop(
    "Package 'dplyr' is required."
  )
  
}


# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------

project_dir <- "YOUR_PROJECT_DIRECTORY"


stage16_dir <- file.path(
  project_dir,
  "results",
  "16_Pseudotime_Gene_Dynamics"
)


stage16d_input_file <- file.path(
  stage16_dir,
  "tables",
  "main",
  "02_Significant_KEGG_Pathways.csv"
)


figures_main_dir <- file.path(
  stage16_dir,
  "figures",
  "main"
)


figures_supplementary_dir <- file.path(
  stage16_dir,
  "figures",
  "supplementary"
)


tables_main_dir <- file.path(
  stage16_dir,
  "tables",
  "main"
)


tables_supplementary_dir <- file.path(
  stage16_dir,
  "tables",
  "supplementary"
)


tables_validation_dir <- file.path(
  stage16_dir,
  "tables",
  "validation"
)


logs_dir <- file.path(
  stage16_dir,
  "logs"
)


objects_dir <- file.path(
  project_dir,
  "objects"
)


# ------------------------------------------------------------
# 3. Create shared Stage 16 directories if necessary
# ------------------------------------------------------------

dir.create(
  figures_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  figures_supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tables_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tables_supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tables_validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Header
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "TNBC Single-Cell Analysis Pipeline\n"
)

cat(
  "STAGE 16E: KEGG PATHWAY INTERPRETATION\n"
)

cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------
# 5. Validate Stage 16D input
# ------------------------------------------------------------

if (!dir.exists(stage16_dir)) {
  
  stop(
    paste(
      "Stage 16 directory does not exist:",
      stage16_dir
    )
  )
  
}


if (!file.exists(stage16d_input_file)) {
  
  stop(
    paste(
      "Canonical Stage 16D KEGG input was not found:",
      stage16d_input_file
    )
  )
  
}


cat(
  "Canonical Stage 16D input:\n",
  stage16d_input_file,
  "\n\n"
)


# ------------------------------------------------------------
# 6. Load significant KEGG pathways from Stage 16D
# ------------------------------------------------------------

cat(
  "Loading significant KEGG pathways from Stage 16D...\n\n"
)


kegg_interpreted <- tryCatch(
  
  read.csv(
    stage16d_input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  
  error = function(e) {
    
    stop(
      paste0(
        "Could not read Stage 16D KEGG result table:\n",
        conditionMessage(e)
      )
    )
    
  }
  
)


if (nrow(kegg_interpreted) == 0) {
  
  stop(
    "Stage 16D significant KEGG pathway table is empty."
  )
  
}


cat(
  "Stage 16D pathways loaded:",
  nrow(kegg_interpreted),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Validate required columns
# ------------------------------------------------------------

required_columns <- c(
  "ID",
  "Description",
  "pvalue",
  "p.adjust"
)


missing_columns <- setdiff(
  required_columns,
  colnames(kegg_interpreted)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Required KEGG columns are missing:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
  
}


# ------------------------------------------------------------
# 8. Convert statistical columns safely
# ------------------------------------------------------------

kegg_interpreted$pvalue <-
  suppressWarnings(
    as.numeric(
      as.character(
        kegg_interpreted$pvalue
      )
    )
  )


kegg_interpreted$p.adjust <-
  suppressWarnings(
    as.numeric(
      as.character(
        kegg_interpreted$p.adjust
      )
    )
  )


if ("qvalue" %in% colnames(kegg_interpreted)) {
  
  kegg_interpreted$qvalue <-
    suppressWarnings(
      as.numeric(
        as.character(
          kegg_interpreted$qvalue
        )
      )
    )
  
}


# ------------------------------------------------------------
# 9. Clean KEGG pathway table
# ------------------------------------------------------------

kegg_interpreted_clean <-
  kegg_interpreted %>%
  
  dplyr::filter(
    !is.na(ID),
    ID != "",
    !is.na(Description),
    Description != "",
    !is.na(p.adjust),
    is.finite(p.adjust)
  ) %>%
  
  dplyr::arrange(
    p.adjust,
    pvalue
  )


if (nrow(kegg_interpreted_clean) == 0) {
  
  stop(
    "No valid KEGG pathways remained after cleaning."
  )
  
}


# ------------------------------------------------------------
# 10. Remove duplicate pathway IDs if present
# ------------------------------------------------------------

duplicate_pathway_ids <-
  sum(
    duplicated(
      kegg_interpreted_clean$ID
    )
  )


if (duplicate_pathway_ids > 0) {
  
  kegg_interpreted_clean <-
    kegg_interpreted_clean %>%
    
    dplyr::group_by(
      ID
    ) %>%
    
    dplyr::slice_min(
      order_by = p.adjust,
      n = 1,
      with_ties = FALSE
    ) %>%
    
    dplyr::ungroup()
  
}


# Re-sort after duplicate removal

kegg_interpreted_clean <-
  kegg_interpreted_clean %>%
  
  dplyr::arrange(
    p.adjust,
    pvalue
  )


cat(
  "Valid KEGG pathways:",
  nrow(kegg_interpreted_clean),
  "\n"
)


cat(
  "Duplicate pathway IDs removed:",
  duplicate_pathway_ids,
  "\n\n"
)


# ------------------------------------------------------------
# 11. Export cleaned pathway table
# ------------------------------------------------------------

write.csv(
  
  kegg_interpreted_clean,
  
  file.path(
    tables_supplementary_dir,
    "16E_01_KEGG_Pathways_Cleaned.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 12. Create ranked pathway table
# ------------------------------------------------------------

pathway_ranking <-
  kegg_interpreted_clean


pathway_ranking$Rank <-
  seq_len(
    nrow(
      pathway_ranking
    )
  )


pathway_ranking <-
  pathway_ranking[
    ,
    c(
      "Rank",
      setdiff(
        colnames(pathway_ranking),
        "Rank"
      )
    ),
    drop = FALSE
  ]


write.csv(
  
  pathway_ranking,
  
  file.path(
    tables_main_dir,
    "16E_02_KEGG_Pathway_Ranking.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 13. Define top 20 pathways
# ------------------------------------------------------------

top20_pathways <-
  utils::head(
    pathway_ranking,
    20
  )


write.csv(
  
  top20_pathways,
  
  file.path(
    tables_main_dir,
    "16E_03_Top20_KEGG_Pathways.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 14. Define top 10 pathways
# ------------------------------------------------------------

top10_pathways <-
  utils::head(
    pathway_ranking,
    10
  )


write.csv(
  
  top10_pathways,
  
  file.path(
    tables_main_dir,
    "16E_04_Top10_KEGG_Pathways.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 15. Create compact top-20 table
# ------------------------------------------------------------

priority_columns <-
  c(
    "Rank",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "RichFactor",
    "FoldEnrichment",
    "zScore",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count"
  )


priority_columns <-
  priority_columns[
    priority_columns %in%
      colnames(
        top20_pathways
      )
  ]


if (length(priority_columns) == 0) {
  
  stop(
    "None of the expected KEGG priority columns were found."
  )
  
}


top20_compact <-
  top20_pathways[
    ,
    priority_columns,
    drop = FALSE
  ]


write.csv(
  
  top20_compact,
  
  file.path(
    tables_main_dir,
    "16E_05_Top20_KEGG_Pathways_Compact.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 16. Create pathway priority summary
# ------------------------------------------------------------

priority_summary <- data.frame(
  
  Rank_Category = c(
    "Top_1",
    "Top_5",
    "Top_10",
    "Top_20"
  ),
  
  Number_of_Pathways = c(
    min(1, nrow(pathway_ranking)),
    min(5, nrow(pathway_ranking)),
    min(10, nrow(pathway_ranking)),
    min(20, nrow(pathway_ranking))
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  
  priority_summary,
  
  file.path(
    tables_supplementary_dir,
    "16E_06_Pathway_Priority_Summary.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 17. Create overall Stage 16E summary
# ------------------------------------------------------------

pathway_summary <- data.frame(
  
  Parameter = c(
    
    "Stage16D_Significant_KEGG_Pathways",
    
    "Stage16E_Valid_KEGG_Pathways",
    
    "Duplicate_Pathway_IDs_Removed",
    
    "Top20_Pathways_Exported",
    
    "Top10_Pathways_Exported",
    
    "Minimum_Adjusted_P_Value",
    
    "Maximum_Adjusted_P_Value"
    
  ),
  
  Value = c(
    
    nrow(kegg_interpreted),
    
    nrow(kegg_interpreted_clean),
    
    duplicate_pathway_ids,
    
    nrow(top20_pathways),
    
    nrow(top10_pathways),
    
    min(
      kegg_interpreted_clean$p.adjust,
      na.rm = TRUE
    ),
    
    max(
      kegg_interpreted_clean$p.adjust,
      na.rm = TRUE
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  
  pathway_summary,
  
  file.path(
    tables_main_dir,
    "16E_07_KEGG_Pathway_Interpretation_Summary.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 18. Stage 16E validation
# ------------------------------------------------------------

validation_table <- data.frame(
  
  Check = c(
    
    "Stage16D_Input_File_Exists",
    
    "Stage16D_Input_Nonempty",
    
    "Required_Columns_Present",
    
    "Valid_Pathways_Available",
    
    "Pathway_IDs_Available",
    
    "Adjusted_P_Values_Finite",
    
    "Ranking_Completed",
    
    "Top20_Exported",
    
    "Top10_Exported"
    
  ),
  
  Result = c(
    
    file.exists(stage16d_input_file),
    
    nrow(kegg_interpreted) > 0,
    
    length(missing_columns) == 0,
    
    nrow(kegg_interpreted_clean) > 0,
    
    all(
      !is.na(
        kegg_interpreted_clean$ID
      )
    ),
    
    all(
      is.finite(
        kegg_interpreted_clean$p.adjust
      )
    ),
    
    identical(
      pathway_ranking$Rank,
      seq_len(
        nrow(pathway_ranking)
      )
    ),
    
    nrow(top20_pathways) ==
      min(
        20,
        nrow(pathway_ranking)
      ),
    
    nrow(top10_pathways) ==
      min(
        10,
        nrow(pathway_ranking)
      )
    
  ),
  
  stringsAsFactors = FALSE
  
)


validation_table$Status <-
  ifelse(
    validation_table$Result,
    "PASS",
    "FAIL"
  )


write.csv(
  
  validation_table,
  
  file.path(
    tables_validation_dir,
    "Stage16E_Input_Validation.csv"
  ),
  
  row.names = FALSE
  
)


if (
  any(
    validation_table$Status == "FAIL"
  )
) {
  
  stop(
    "Stage 16E validation failed."
  )
  
}


# ------------------------------------------------------------
# 19. Print Top 20 pathways
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "TOP 20 KEGG PATHWAYS\n"
)

cat(
  "============================================================\n\n"
)


print(
  top20_compact
)


# ------------------------------------------------------------
# 20. Save Stage 16E RDS
# ------------------------------------------------------------

stage16e_results <- list(
  
  source_file =
    stage16d_input_file,
  
  kegg_interpreted =
    kegg_interpreted,
  
  kegg_interpreted_clean =
    kegg_interpreted_clean,
  
  pathway_ranking =
    pathway_ranking,
  
  top20 =
    top20_pathways,
  
  top10 =
    top10_pathways,
  
  top20_compact =
    top20_compact,
  
  priority_summary =
    priority_summary,
  
  pathway_summary =
    pathway_summary,
  
  validation =
    validation_table,
  
  stage15_status =
    "FROZEN",
  
  stage16a_status =
    "FROZEN",
  
  stage16b_status =
    "FROZEN",
  
  stage16c_status =
    "FROZEN",
  
  stage16d_status =
    "FROZEN",
  
  stage16e_status =
    "COMPLETED"
  
)


stage16e_rds <-
  file.path(
    objects_dir,
    "Stage16E_KEGG_Interpretation_Results.rds"
  )


saveRDS(
  
  stage16e_results,
  
  stage16e_rds
  
)


# ------------------------------------------------------------
# 21. Create Top-20 text report
# ------------------------------------------------------------

report_file <-
  file.path(
    tables_supplementary_dir,
    "16E_08_Top20_KEGG_Pathways_Report.txt"
  )


report_lines <- c(
  
  "============================================================",
  
  "TNBC Single-Cell Analysis Pipeline",
  
  "Stage 16E: Top 20 KEGG Pathways",
  
  "============================================================",
  
  "",
  
  paste(
    "Generated:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste(
    "Stage 16D significant pathways:",
    nrow(kegg_interpreted)
  ),
  
  paste(
    "Valid pathways:",
    nrow(kegg_interpreted_clean)
  ),
  
  "",
  
  "Top pathways:",
  
  ""
  
)


if (
  nrow(top20_compact) > 0
) {
  
  for (
    i in seq_len(
      nrow(top20_compact)
    )
  ) {
    
    pathway_id <-
      if (
        "ID" %in%
        colnames(top20_compact)
      ) {
        
        as.character(
          top20_compact$ID[i]
        )
        
      } else {
        
        ""
        
      }
    
    
    pathway_description <-
      if (
        "Description" %in%
        colnames(top20_compact)
      ) {
        
        as.character(
          top20_compact$Description[i]
        )
        
      } else {
        
        ""
        
      }
    
    
    padj_value <-
      if (
        "p.adjust" %in%
        colnames(top20_compact)
      ) {
        
        format(
          top20_compact$`p.adjust`[i],
          scientific = TRUE,
          digits = 4
        )
        
      } else {
        
        ""
        
      }
    
    
    report_lines <-
      c(
        
        report_lines,
        
        paste0(
          i,
          ". ",
          pathway_id,
          " | ",
          pathway_description,
          " | adjusted p = ",
          padj_value
        )
        
      )
    
  }
  
}


report_lines <-
  c(
    
    report_lines,
    
    "",
    
    "============================================================",
    
    "Stage 16E performs interpretation and prioritization only.",
    
    "KEGG enrichment was not rerun.",
    
    "Stage 15 through Stage 16D remained frozen.",
    
    "============================================================"
    
  )


writeLines(
  
  report_lines,
  
  con = report_file
  
)


# ------------------------------------------------------------
# 22. Completion log
# ------------------------------------------------------------

stage16e_log <-
  file.path(
    logs_dir,
    "Stage16E_Completion.log"
  )


log_lines <- c(
  
  "============================================================",
  
  "TNBC Single-Cell Analysis Pipeline",
  
  "Stage 16E: KEGG Pathway Interpretation",
  
  "============================================================",
  
  "",
  
  paste(
    "Completion time:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste(
    "Stage 16D input:",
    stage16d_input_file
  ),
  
  paste(
    "Stage 16D significant pathways:",
    nrow(kegg_interpreted)
  ),
  
  paste(
    "Valid pathways:",
    nrow(kegg_interpreted_clean)
  ),
  
  paste(
    "Duplicate pathway IDs removed:",
    duplicate_pathway_ids
  ),
  
  paste(
    "Top 20 pathways:",
    nrow(top20_pathways)
  ),
  
  paste(
    "Top 10 pathways:",
    nrow(top10_pathways)
  ),
  
  "",
  
  "Stage 15: FROZEN",
  
  "Stage 16A: FROZEN",
  
  "Stage 16B: FROZEN",
  
  "Stage 16C: FROZEN",
  
  "Stage 16D: FROZEN",
  
  "Stage 16E: COMPLETED",
  
  "",
  
  "KEGG enrichment: NOT RERUN",
  
  "Trajectory reconstruction: NOT PERFORMED",
  
  "Pseudotime recalculation: NOT PERFORMED",
  
  "graph_test: NOT RERUN",
  
  "",
  
  paste(
    "Complete Stage 16E RDS:",
    stage16e_rds
  ),
  
  "",
  
  "STATUS: STAGE 16E COMPLETED SUCCESSFULLY",
  
  "============================================================"
  
)


writeLines(
  
  log_lines,
  
  con = stage16e_log
  
)


# ------------------------------------------------------------
# 23. Final output
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "STAGE 16E COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)

cat(
  "Stage 16D significant pathways:",
  nrow(kegg_interpreted),
  "\n"
)

cat(
  "Valid pathways:",
  nrow(kegg_interpreted_clean),
  "\n"
)

cat(
  "Duplicate pathway IDs removed:",
  duplicate_pathway_ids,
  "\n"
)

cat(
  "Top 20 pathways:",
  nrow(top20_pathways),
  "\n"
)

cat(
  "Top 10 pathways:",
  nrow(top10_pathways),
  "\n"
)

cat(
  "\nComplete Stage 16E RDS:\n",
  stage16e_rds,
  "\n"
)

cat(
  "\nCompletion log:\n",
  stage16e_log,
  "\n"
)

cat(
  "\n============================================================\n"
)

