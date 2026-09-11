############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
#
# STAGE 13: GSEA Biological Interpretation
#           Representative GSEA Pathways
#
# Comparison:
# BRCA1_tumour vs TotalCell
#
# PURPOSE:
# - Interpret significant integrated GSEA results
# - Reduce biological redundancy among related pathways
# - Select representative pathways for each direction
# - Limit final representative pathways to a maximum of
#   10 terms per enrichment direction
#
# INPUT:
# Stage 12 integrated significant GSEA results
#
# OUTPUT:
# - Representative pathway tables
# - Redundancy-group summary
# - Database distribution
# - QC summary
# - Representative pathway figure
# - Stage 13 log
#
# INTERPRETATION:
# NES > 0  = enrichment toward BRCA1_tumour
# NES < 0  = enrichment toward TotalCell
############################################################


############################################################
# 13.1 Load libraries
############################################################

library(dplyr)
library(ggplot2)
library(stringr)


############################################################
# 13.2 Project and Stage directories
############################################################

project_dir <- "YOUR_PROJECT_DIRECTORY"

stage13_dir <- file.path(
  project_dir,
  "results",
  "13_GSEA_Biological_Interpretation"
)

figures_main_dir <- file.path(
  stage13_dir,
  "figures",
  "main"
)

figures_supplementary_dir <- file.path(
  stage13_dir,
  "figures",
  "supplementary"
)

tables_main_dir <- file.path(
  stage13_dir,
  "tables",
  "main"
)

tables_supplementary_dir <- file.path(
  stage13_dir,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  stage13_dir,
  "tables",
  "validation"
)

logs_dir <- file.path(
  stage13_dir,
  "logs"
)


############################################################
# Create required directories
############################################################

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


############################################################
# 13.3 Load Stage 12 integrated GSEA results
############################################################

integrated_file <- file.path(
  project_dir,
  "results",
  "12_GSEA_Integration",
  "tables",
  "main",
  "GSEA_Integrated_significant.csv"
)

if (!file.exists(integrated_file)) {
  
  stop(
    "Stage 12 integrated GSEA file was not found:\n",
    integrated_file
  )
}


gsea_integrated <- read.csv(
  integrated_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 13.4 Validate required columns
############################################################

required_columns <- c(
  "Database",
  "Description",
  "NES",
  "p.adjust",
  "Enrichment_Group"
)

missing_columns <- setdiff(
  required_columns,
  colnames(gsea_integrated)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


############################################################
# 13.5 Keep statistically significant terms
############################################################

gsea_sig <- gsea_integrated %>%
  
  filter(
    !is.na(Description),
    Description != "",
    !is.na(NES),
    !is.na(p.adjust),
    p.adjust < 0.05
  )


if (nrow(gsea_sig) == 0) {
  
  stop(
    "No statistically significant GSEA terms remain after filtering."
  )
}


############################################################
# 13.6 Validate enrichment directions
############################################################

unexpected_direction <- gsea_sig %>%
  
  filter(
    
    (
      Enrichment_Group == "BRCA1_tumour" &
        NES <= 0
    ) |
      
      (
        Enrichment_Group == "TotalCell" &
          NES >= 0
      )
    
  )


if (nrow(unexpected_direction) > 0) {
  
  stop(
    "Unexpected NES direction detected in Stage 12 input. ",
    "Please inspect the integrated GSEA results."
  )
}


############################################################
# 13.7 Remove exact duplicate descriptions
############################################################

gsea_sig <- gsea_sig %>%
  
  arrange(
    Enrichment_Group,
    p.adjust,
    desc(abs(NES))
  ) %>%
  
  distinct(
    Description,
    Enrichment_Group,
    .keep_all = TRUE
  )


############################################################
# 13.8 Biological redundancy grouping
############################################################
#
# Conservative rule-based grouping.
#
# Terms that clearly describe the same broad biological
# process are grouped together.
#
# Terms without a sufficiently clear shared biological
# theme remain independent.
############################################################

gsea_sig <- gsea_sig %>%
  
  mutate(
    
    Redundancy_Group = case_when(
      
      ######################################################
      # Translation / ribosome / protein synthesis
      ######################################################
      
      grepl(
        "translation|ribosome|peptide chain elongation|cotranslational|protein targeting|translation termination|40S|60S|GTP hydrolysis.*ribosom",
        Description,
        ignore.case = TRUE
      ) ~ "Translation_Ribosome",
      
      
      ######################################################
      # Nonsense-mediated decay
      ######################################################
      
      grepl(
        "nonsense.*mediated.*decay|NMD",
        Description,
        ignore.case = TRUE
      ) ~ "Nonsense_Mediated_Decay",
      
      
      ######################################################
      # Selenocysteine / selenoamino acid metabolism
      ######################################################
      
      grepl(
        "selenocysteine|selenoamino",
        Description,
        ignore.case = TRUE
      ) ~ "Selenocysteine_Selenoamino_Acid_Metabolism",
      
      
      ######################################################
      # Viral infection / viral response
      ######################################################
      
      grepl(
        "influenza|viral|virus",
        Description,
        ignore.case = TRUE
      ) ~ "Viral_Response_Infection",
      
      
      ######################################################
      # TNF / NF-kB signaling
      ######################################################
      
      grepl(
        "TNFA|TNF|NFKB|NF.KB",
        Description,
        ignore.case = TRUE
      ) ~ "TNF_NFkB_Signaling",
      
      
      ######################################################
      # Inflammatory response
      ######################################################
      
      grepl(
        "inflammatory response|inflammation",
        Description,
        ignore.case = TRUE
      ) ~ "Inflammatory_Response",
      
      
      ######################################################
      # Wound healing / tissue repair
      ######################################################
      
      grepl(
        "wound healing|response to wounding|tissue repair",
        Description,
        ignore.case = TRUE
      ) ~ "Wound_Healing_Tissue_Repair",
      
      
      ######################################################
      # Microbial / bacterial response
      ######################################################
      
      grepl(
        "response to bacterium|response to bacteria|bacterial|microbial",
        Description,
        ignore.case = TRUE
      ) ~ "Microbial_Response",
      
      
      ######################################################
      # Migration / chemotaxis / locomotion
      ######################################################
      
      grepl(
        "locomotion|chemotaxis|taxis|cell migration|cell motility|leukocyte migration|regulation of cell motility",
        Description,
        ignore.case = TRUE
      ) ~ "Migration_Chemotaxis",
      
      
      ######################################################
      # Cell adhesion
      ######################################################
      
      grepl(
        "cell adhesion|cell-cell adhesion|cell substrate adhesion",
        Description,
        ignore.case = TRUE
      ) ~ "Cell_Adhesion",
      
      
      ######################################################
      # Cell activation
      ######################################################
      
      grepl(
        "cell activation",
        Description,
        ignore.case = TRUE
      ) ~ "Cell_Activation",
      
      
      ######################################################
      # Lipid response / metabolism
      ######################################################
      
      grepl(
        "response to lipid|lipid metabolism",
        Description,
        ignore.case = TRUE
      ) ~ "Lipid_Response_Metabolism",
      
      
      ######################################################
      # Tissue morphogenesis
      ######################################################
      
      grepl(
        "morphogenesis|morphogenetic|trabecula formation|trabecular formation",
        Description,
        ignore.case = TRUE
      ) ~ "Tissue_Morphogenesis",
      
      
      ######################################################
      # Melanogenesis
      ######################################################
      
      grepl(
        "melanogenesis",
        Description,
        ignore.case = TRUE
      ) ~ "Melanogenesis",
      
      
      ######################################################
      # Independent biological term
      ######################################################
      
      TRUE ~ paste0(
        "Independent_",
        Description
      )
    )
  )


############################################################
# 13.9 Redundancy-group summary
############################################################

group_summary <- gsea_sig %>%
  
  count(
    Enrichment_Group,
    Redundancy_Group,
    name = "n_terms"
  ) %>%
  
  arrange(
    Enrichment_Group,
    desc(n_terms),
    Redundancy_Group
  )


write.csv(
  group_summary,
  file.path(
    tables_validation_dir,
    "Step13_Redundancy_Group_Summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 13.10 Representative pathway selection function
############################################################

select_representative <- function(
    df,
    max_terms = 10
) {
  
  if (nrow(df) == 0) {
    
    return(df)
    
  }
  
  
  df %>%
    
    arrange(
      p.adjust,
      desc(abs(NES))
    ) %>%
    
    group_by(
      Redundancy_Group
    ) %>%
    
    slice_head(
      n = 1
    ) %>%
    
    ungroup() %>%
    
    arrange(
      p.adjust,
      desc(abs(NES))
    ) %>%
    
    slice_head(
      n = max_terms
    )
}


############################################################
# 13.11 BRCA1_tumour representative pathways
############################################################

gsea_brca1 <- gsea_sig %>%
  
  filter(
    Enrichment_Group == "BRCA1_tumour"
  )


representative_brca1 <-
  select_representative(
    gsea_brca1,
    max_terms = 10
  )


############################################################
# 13.12 TotalCell representative pathways
############################################################

gsea_totalcell <- gsea_sig %>%
  
  filter(
    Enrichment_Group == "TotalCell"
  )


representative_totalcell <-
  select_representative(
    gsea_totalcell,
    max_terms = 10
  )


############################################################
# 13.13 Combine representative pathways
############################################################

representative_gsea <- bind_rows(
  representative_brca1,
  representative_totalcell
)


############################################################
# 13.14 Add wrapped descriptions
############################################################

representative_gsea <- representative_gsea %>%
  
  mutate(
    
    Description_wrapped = str_wrap(
      Description,
      width = 55
    )
    
  )


############################################################
# 13.15 Save main representative pathway table
############################################################

write.csv(
  representative_gsea,
  file.path(
    tables_main_dir,
    "Step13_Representative_GSEA_Pathways.csv"
  ),
  row.names = FALSE
)


############################################################
# 13.16 Save direction-specific tables
############################################################

write.csv(
  representative_brca1,
  file.path(
    tables_supplementary_dir,
    "Step13_Representative_GSEA_BRCA1_tumour.csv"
  ),
  row.names = FALSE
)


write.csv(
  representative_totalcell,
  file.path(
    tables_supplementary_dir,
    "Step13_Representative_GSEA_TotalCell.csv"
  ),
  row.names = FALSE
)


############################################################
# 13.17 Database distribution
############################################################

database_distribution <-
  representative_gsea %>%
  
  count(
    Database,
    Enrichment_Group
  )


write.csv(
  database_distribution,
  file.path(
    tables_supplementary_dir,
    "Step13_Representative_GSEA_database_distribution.csv"
  ),
  row.names = FALSE
)


############################################################
# 13.18 Representative pathway plot
############################################################

plot_gsea <- representative_gsea %>%
  
  mutate(
    
    Direction = factor(
      Enrichment_Group,
      levels = c(
        "BRCA1_tumour",
        "TotalCell"
      )
    ),
    
    Description_wrapped = str_wrap(
      Description,
      width = 55
    )
    
  )


if (nrow(plot_gsea) > 0) {
  
  plot_gsea <- plot_gsea %>%
    
    mutate(
      
      Plot_Label = paste0(
        Database,
        ": ",
        Description_wrapped
      ),
      
      Plot_Label = reorder(
        Plot_Label,
        NES
      )
      
    )
  
  
  p_representative <- ggplot(
    
    plot_gsea,
    
    aes(
      x = NES,
      y = Plot_Label,
      size = -log10(p.adjust)
    )
    
  ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    
    geom_point(
      alpha = 0.85
    ) +
    
    facet_wrap(
      ~ Direction,
      scales = "free_y"
    ) +
    
    theme_classic() +
    
    labs(
      
      title =
        "Representative GSEA Pathways",
      
      subtitle =
        "Non-redundant significant pathways",
      
      x =
        "Normalized Enrichment Score (NES)",
      
      y =
        NULL,
      
      size =
        "-log10 adjusted p-value"
      
    ) +
    
    theme(
      
      axis.text.y =
        element_text(
          size = 9
        ),
      
      strip.text =
        element_text(
          face = "bold"
        ),
      
      plot.title =
        element_text(
          face = "bold"
        )
      
    )
  
  
  ggsave(
    
    file.path(
      figures_main_dir,
      "Step13_Representative_GSEA_Pathways.pdf"
    ),
    
    p_representative,
    
    width = 13,
    
    height =
      max(
        7,
        0.45 * nrow(plot_gsea)
      )
    
  )
}


############################################################
# 13.19 Summary statistics
############################################################

n_representative_total <-
  nrow(
    representative_gsea
  )

n_representative_brca1 <-
  nrow(
    representative_brca1
  )

n_representative_totalcell <-
  nrow(
    representative_totalcell
  )


############################################################
# 13.20 Quality control
############################################################

############################################################
# QC 1 — No duplicate redundancy groups within direction
############################################################

duplicate_group_check <-
  representative_gsea %>%
  
  count(
    Enrichment_Group,
    Redundancy_Group,
    name = "n"
  ) %>%
  
  filter(
    n > 1
  )


qc_no_duplicate_groups <-
  nrow(
    duplicate_group_check
  ) == 0


############################################################
# QC 2 — Maximum 10 representative terms per direction
############################################################

representative_counts <-
  representative_gsea %>%
  
  count(
    Enrichment_Group,
    name = "n_terms"
  )


qc_max_10 <-
  
  all(
    representative_counts$n_terms <= 10
  )


############################################################
# QC 3 — Correct NES direction
############################################################

qc_direction <-
  
  all(
    representative_brca1$NES > 0
  ) &&
  
  all(
    representative_totalcell$NES < 0
  )


############################################################
# QC 4 — No missing redundancy groups
############################################################

qc_no_missing_groups <-
  
  all(
    !is.na(
      representative_gsea$Redundancy_Group
    )
  )


############################################################
# QC 5 — Both expected directions represented
############################################################

qc_both_directions <-
  
  all(
    c(
      "BRCA1_tumour",
      "TotalCell"
    ) %in%
      representative_gsea$Enrichment_Group
  )


############################################################
# QC summary
############################################################

qc_summary <- data.frame(
  
  QC_Check = c(
    
    "No duplicate redundancy groups within direction",
    
    "Maximum 10 representative terms per direction",
    
    "Correct NES direction",
    
    "No missing redundancy groups",
    
    "Both enrichment directions represented"
    
  ),
  
  Result = c(
    
    qc_no_duplicate_groups,
    
    qc_max_10,
    
    qc_direction,
    
    qc_no_missing_groups,
    
    qc_both_directions
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  qc_summary,
  file.path(
    tables_validation_dir,
    "Step13_QC_Summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 13.21 Mandatory QC stop
############################################################

if (!all(qc_summary$Result)) {
  
  stop(
    "Stage 13 QC FAILED. ",
    "Please inspect Step13_QC_Summary.csv."
  )
}


############################################################
# 13.22 Output validation
############################################################

expected_outputs <- c(
  
  file.path(
    tables_main_dir,
    "Step13_Representative_GSEA_Pathways.csv"
  ),
  
  file.path(
    tables_supplementary_dir,
    "Step13_Representative_GSEA_BRCA1_tumour.csv"
  ),
  
  file.path(
    tables_supplementary_dir,
    "Step13_Representative_GSEA_TotalCell.csv"
  ),
  
  file.path(
    tables_supplementary_dir,
    "Step13_Representative_GSEA_database_distribution.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Step13_Redundancy_Group_Summary.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Step13_QC_Summary.csv"
  ),
  
  file.path(
    figures_main_dir,
    "Step13_Representative_GSEA_Pathways.pdf"
  )
)


missing_outputs <- expected_outputs[
  !file.exists(expected_outputs)
]


if (length(missing_outputs) > 0) {
  
  stop(
    "Stage 13 finished with missing expected outputs:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}


############################################################
# 13.23 Console summary
############################################################

message("")
message("====================================================")
message("STAGE 13 COMPLETED SUCCESSFULLY")
message("====================================================")

message(
  "Significant input terms: ",
  nrow(gsea_sig)
)

message(
  "Distinct redundancy groups: ",
  n_distinct(
    gsea_sig$Redundancy_Group
  )
)

message(
  "Representative terms overall: ",
  n_representative_total
)

message(
  "  BRCA1_tumour: ",
  n_representative_brca1
)

message(
  "  TotalCell: ",
  n_representative_totalcell
)

message(
  "Redundancy reduction: PASSED"
)

message(
  "Maximum 10 terms per direction: PASSED"
)

message(
  "NES direction validation: PASSED"
)

message(
  "Missing redundancy groups: PASSED"
)

message(
  "Both enrichment directions represented: PASSED"
)

message(
  "Output validation: PASSED"
)

message(
  "Stage 13 output directory: ",
  stage13_dir
)


############################################################
# 13.24 Write Stage 13 log
############################################################

log_file <- file.path(
  logs_dir,
  "Step13_Representative_GSEA.log"
)


log_lines <- c(
  
  "====================================================",
  
  "TNBC Single-cell RNA-seq Analysis Pipeline",
  
  "STAGE 13 — GSEA Biological Interpretation",
  
  "Representative GSEA Pathways",
  
  "====================================================",
  
  "",
  
  paste(
    "Date/time:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  "INPUT",
  
  "----------------------------------------------------",
  
  paste(
    "Stage 12 integrated GSEA:",
    integrated_file
  ),
  
  "Significant input criterion: p.adjust < 0.05",
  
  "",
  
  "SELECTION STRATEGY",
  
  "----------------------------------------------------",
  
  "Exact duplicate descriptions were removed.",
  
  "Rule-based biological redundancy groups were assigned.",
  
  "One representative term was retained per redundancy group.",
  
  "A maximum of 10 representative terms was retained per direction.",
  
  "BRCA1_tumour and TotalCell were evaluated separately.",
  
  "",
  
  "RESULTS",
  
  "----------------------------------------------------",
  
  paste(
    "Significant input terms:",
    nrow(gsea_sig)
  ),
  
  paste(
    "Distinct redundancy groups:",
    n_distinct(
      gsea_sig$Redundancy_Group
    )
  ),
  
  paste(
    "Representative terms overall:",
    n_representative_total
  ),
  
  paste(
    "Representative BRCA1_tumour terms:",
    n_representative_brca1
  ),
  
  paste(
    "Representative TotalCell terms:",
    n_representative_totalcell
  ),
  
  "",
  
  "QUALITY CONTROL",
  
  "----------------------------------------------------",
  
  paste(
    "No duplicate redundancy groups:",
    ifelse(
      qc_no_duplicate_groups,
      "PASSED",
      "FAILED"
    )
  ),
  
  paste(
    "Maximum 10 terms per direction:",
    ifelse(
      qc_max_10,
      "PASSED",
      "FAILED"
    )
  ),
  
  paste(
    "Correct NES direction:",
    ifelse(
      qc_direction,
      "PASSED",
      "FAILED"
    )
  ),
  
  paste(
    "No missing redundancy groups:",
    ifelse(
      qc_no_missing_groups,
      "PASSED",
      "FAILED"
    )
  ),
  
  paste(
    "Both enrichment directions represented:",
    ifelse(
      qc_both_directions,
      "PASSED",
      "FAILED"
    )
  ),
  
  "Output validation: PASSED",
  
  "",
  
  "INTERPRETATION",
  
  "----------------------------------------------------",
  
  "NES > 0 = enrichment toward BRCA1_tumour",
  
  "NES < 0 = enrichment toward TotalCell",
  
  "Representative pathways summarize major biological themes",
  
  "rather than listing every significant enrichment term.",
  
  "",
  
  "STAGE 13 COMPLETED SUCCESSFULLY.",
  
  "===================================================="
)


writeLines(
  log_lines,
  con = log_file
)


############################################################
# END OF STAGE 13
############################################################


############################################################
# OPTIONAL INTERACTIVE QC
############################################################

representative_brca1 %>%
  
  select(
    Database,
    Description,
    NES,
    p.adjust,
    Redundancy_Group
  )


representative_totalcell %>%
  
  select(
    Database,
    Description,
    NES,
    p.adjust,
    Redundancy_Group
  )


qc_summary


table(
  representative_gsea$Database,
  representative_gsea$Enrichment_Group
)

