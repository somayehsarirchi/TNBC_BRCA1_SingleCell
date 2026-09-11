#!/usr/bin/env Rscript

############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 7: CopyKAT copy-number inference on Linux
#
# This script:
# - Loads the CopyKAT input object generated in Step 6
# - Validates the raw count matrix
# - Filters low-quality genes and cells
# - Splits the dataset into manageable batches
# - Runs CopyKAT independently on each batch
# - Collects and standardizes prediction tables
# - Merges batch-level predictions into a final result
# - Saves all outputs and run parameters
#
# Input:
#   objects/copykat_inputs.rds
#
# Outputs:
#   results/copykat_batches/
#   results/copykat_final/copykat_combined_predictions.csv
#   results/copykat_final/copykat_combined_predictions.rds
#   results/copykat_final/copykat_run_parameters.rds
#
# Notes:
# - This script is intended to be executed on Linux.
# - The input matrix must contain raw UMI counts.
############################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(copykat)
})

set.seed(1234)

############################################################
# Project directories
############################################################

project_dir <- "YOUR_PROJECT_DIRECTORY"

objects_dir <- file.path(project_dir, "objects")
results_dir <- file.path(project_dir, "results")

input_rds <- file.path(objects_dir, "copykat_inputs.rds")
output_dir <- file.path(results_dir, "copykat_batches")
final_output_dir <- file.path(results_dir, "copykat_final")

############################################################
# Batch execution
#
# This script is designed to be executed from the Linux terminal:
#
#   nohup Rscript run_copykat_batches.R > run_copykat_batches.log 2>&1 &
#
# The pipeline automatically:
# - loads copykat_inputs.rds
# - splits cells into batches
# - runs CopyKAT on each batch
# - saves intermediate results
# - combines prediction tables across all batches
############################################################

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_output_dir, recursive = TRUE, showWarnings = FALSE)

target_batch_size <- 2500
min_genes_per_cell <- 200
min_cells_per_gene <- 3
min_normal_cells_per_batch <- 20

genome <- "hg20"
id_type <- "S"
distance <- "pearson"
n_cores <- 4

log_msg <- function(...) {
  msg <- paste0(...)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n", sep = "")
  flush.console()
}

save_rds_safe <- function(object, path) {
  saveRDS(object, path)
  log_msg("Saved: ", path)
}

get_matrix_dim <- function(x) {
  d <- dim(x)

  if (!is.null(d) && length(d) == 2) {
    return(d)
  }

  if (inherits(x, "Matrix") && length(x@Dim) == 2) {
    return(x@Dim)
  }

  NULL
}
############################################################
# Extract a raw count matrix from the input object
#
# The input RDS may contain:
# - a matrix,
# - a sparse Matrix object,
# - or a list containing a matrix under names such as
#   counts, raw_counts, matrix, or rawmat.
############################################################
extract_counts_object <- function(x) {
  if (
    inherits(x, "matrix") ||
    inherits(x, "Matrix") ||
    inherits(x, "sparseMatrix") ||
    inherits(x, "dgCMatrix") ||
    inherits(x, "dgTMatrix")
  ) {
    log_msg("Input is already matrix-like")
    return(x)
  }

  if (is.data.frame(x)) {
    log_msg("Input is a data.frame; converting to matrix")
    return(as.matrix(x))
  }

  if (is.list(x)) {
    log_msg("Input is a list. Names: ", paste(names(x), collapse = ", "))

    preferred_names <- c(
      "counts",
      "raw_counts",
      "copykat_counts",
      "rawmat",
      "rawmat_counts",
      "matrix"
    )

    forbidden_names <- c(
      "data",
      "normalized",
      "norm_data",
      "scale.data",
      "scaled",
      "expr",
      "expression"
    )

    for (nm in preferred_names) {
      if (nm %in% names(x)) {
        candidate <- x[[nm]]
        candidate_dim <- get_matrix_dim(candidate)

        if (!is.null(candidate_dim)) {
          log_msg("Using list element as raw counts matrix: ", nm)
          return(candidate)
        }
      }
    }

    for (nm in forbidden_names) {
      if (nm %in% names(x)) {
        candidate <- x[[nm]]
        candidate_dim <- get_matrix_dim(candidate)

        if (!is.null(candidate_dim)) {
          stop(
            "Found matrix-like element named '", nm, "', but this name suggests normalized data. ",
            "CopyKAT must be run on raw counts. Please store raw counts as 'counts' or 'raw_counts'."
          )
        }
      }
    }

    for (nm in names(x)) {
      candidate <- x[[nm]]
      candidate_dim <- get_matrix_dim(candidate)

      if (!is.null(candidate_dim)) {
        log_msg("Using first matrix-like list element: ", nm)
        return(candidate)
      }
    }
  }

  stop("Could not extract a counts matrix from the input RDS.")
}
############################################################
# Validate that the input matrix contains raw integer counts
#
# CopyKAT should be run on raw UMI/count matrices rather than
# normalized or log-transformed expression values.
############################################################
validate_raw_counts <- function(counts) {
  sample_values <- counts

  if (inherits(counts, "Matrix")) {
    nonzero_values <- counts@x
    if (length(nonzero_values) > 100000) {
      sample_values <- sample(nonzero_values, 100000)
    } else {
      sample_values <- nonzero_values
    }
  } else {
    values <- as.vector(counts)
    values <- values[values != 0]
    if (length(values) > 100000) {
      sample_values <- sample(values, 100000)
    } else {
      sample_values <- values
    }
  }

  if (length(sample_values) == 0) {
    stop("Counts matrix appears to contain only zeros.")
  }

  if (any(sample_values < 0, na.rm = TRUE)) {
    stop("Counts matrix contains negative values. This is not valid raw UMI/count data.")
  }

  non_integer_fraction <- mean(abs(sample_values - round(sample_values)) > 1e-6, na.rm = TRUE)

  if (non_integer_fraction > 0.01) {
    stop(
      "Counts matrix contains many non-integer values. ",
      "This suggests normalized data, not raw counts. ",
      "CopyKAT should be run on raw counts."
    )
  }

  log_msg("Raw-count validation passed")
}

standardize_prediction <- function(pred, batch_id, batch_dir) {
  pred <- as.data.frame(pred, stringsAsFactors = FALSE)

  possible_cell_columns <- c(
    "cell.names",
    "cell_name",
    "cell",
    "barcode",
    "barcodes",
    "Cell",
    "cells"
  )

  found_cell_column <- intersect(possible_cell_columns, colnames(pred))

  if (length(found_cell_column) == 0) {
    if (!is.null(rownames(pred)) && !all(rownames(pred) %in% c("", seq_len(nrow(pred))))) {
      pred$cell.names <- rownames(pred)
      log_msg("Added cell.names column from prediction rownames")
    } else {
      stop("Prediction table has no recognizable cell ID column and no useful rownames.")
    }
  } else if (!"cell.names" %in% colnames(pred)) {
    pred$cell.names <- pred[[found_cell_column[1]]]
    log_msg("Copied cell ID column to cell.names from: ", found_cell_column[1])
  }

  pred$copykat_batch <- batch_id
  pred$copykat_batch_dir <- basename(batch_dir)

  pred
}

run_copykat_one_batch <- function(batch_counts, batch_normal_cells, batch_id, batch_dir) {
  dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)

  old_wd <- getwd()
  setwd(batch_dir)
  on.exit(setwd(old_wd), add = TRUE)

  log_msg("Running CopyKAT for batch ", batch_id)
  log_msg("Batch cells: ", ncol(batch_counts))
  log_msg("Batch genes: ", nrow(batch_counts))
  log_msg("Batch normal cells: ", length(batch_normal_cells))
  log_msg("Batch matrix class before dense conversion: ", paste(class(batch_counts), collapse = ", "))

  if (length(batch_normal_cells) < min_normal_cells_per_batch) {
    log_msg(
      "WARNING: batch ", batch_id, " has only ", length(batch_normal_cells),
      " normal cells. Recommended minimum: ", min_normal_cells_per_batch
    )
  }

  batch_counts <- as.matrix(batch_counts)
  storage.mode(batch_counts) <- "numeric"

  log_msg("Batch matrix class after dense conversion: ", paste(class(batch_counts), collapse = ", "))
  log_msg("Dense matrix size: ", format(object.size(batch_counts), units = "GB"))

  result <- copykat(
    rawmat = batch_counts,
    id.type = id_type,
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = paste0("copykat_batch_", batch_id),
    distance = distance,
    norm.cell.names = if (length(batch_normal_cells) > 0) batch_normal_cells else NULL,
    genome = genome,
    n.cores = n_cores
  )

  save_rds_safe(result, file.path(batch_dir, paste0("copykat_result_batch_", batch_id, ".rds")))

  pred <- NULL

  if (is.list(result) && "prediction" %in% names(result)) {
    pred <- result[["prediction"]]
  }

  if (!is.null(pred)) {
    pred <- standardize_prediction(pred, batch_id, batch_dir)

    write.csv(
      pred,
      file = file.path(batch_dir, paste0("copykat_prediction_batch_", batch_id, ".csv")),
      row.names = FALSE
    )

    log_msg("Saved prediction CSV for batch ", batch_id)
  } else {
    log_msg("WARNING: no prediction table found for batch ", batch_id)
  }

  rm(batch_counts)
  gc()

  result
}

log_msg("Starting CopyKAT batch script")
log_msg("Base directory: ", project_dir)
log_msg("Input RDS: ", input_rds)
log_msg("Output directory: ", output_dir)
log_msg("Final output directory: ", final_output_dir)
log_msg("Target batch size: ", target_batch_size)
log_msg("Random seed: 1234")

if (!file.exists(input_rds)) {
  stop("Input RDS does not exist: ", input_rds)
}

log_msg("Saving sessionInfo")
capture.output(
  sessionInfo(),
  file = file.path(final_output_dir, "sessionInfo_start.txt")
)

log_msg("Loading input RDS")
input_obj <- readRDS(input_rds)

normal_cells <- character(0)

if (is.list(input_obj) && "normal_cells" %in% names(input_obj)) {
  normal_cells <- as.character(input_obj[["normal_cells"]])
  normal_cells <- normal_cells[!is.na(normal_cells)]
  normal_cells <- unique(normal_cells)
  log_msg("Found normal_cells in input RDS: ", length(normal_cells))
}

counts <- extract_counts_object(input_obj)

rm(input_obj)
gc()

counts_dim <- get_matrix_dim(counts)

if (is.null(counts_dim) || length(counts_dim) != 2) {
  stop("Extracted counts object is not a valid 2-dimensional matrix.")
}

log_msg("Counts class: ", paste(class(counts), collapse = ", "))
log_msg("Counts dimensions: ", counts_dim[1], " genes x ", counts_dim[2], " cells")

if (is.null(rownames(counts))) {
  stop("Counts matrix has no gene names in rownames(counts).")
}

if (is.null(colnames(counts))) {
  stop("Counts matrix has no cell names in colnames(counts).")
}

if (anyDuplicated(rownames(counts)) > 0) {
  stop("Counts matrix has duplicated gene names.")
}

if (anyDuplicated(colnames(counts)) > 0) {
  stop("Counts matrix has duplicated cell names.")
}

if (counts_dim[1] == 0) {
  stop("Counts matrix has zero genes.")
}

if (counts_dim[2] == 0) {
  stop("Counts matrix has zero cells.")
}

validate_raw_counts(counts)

log_msg("Filtering empty/low-quality genes and cells")

cell_gene_counts <- Matrix::colSums(counts > 0)
gene_cell_counts <- Matrix::rowSums(counts > 0)

keep_cells <- cell_gene_counts >= min_genes_per_cell
keep_genes <- gene_cell_counts >= min_cells_per_gene

log_msg("Cells before filtering: ", ncol(counts))
log_msg("Cells after filtering: ", sum(keep_cells))
log_msg("Genes before filtering: ", nrow(counts))
log_msg("Genes after filtering: ", sum(keep_genes))

if (sum(keep_cells) == 0) {
  stop("No cells remain after filtering.")
}

if (sum(keep_genes) == 0) {
  stop("No genes remain after filtering.")
}

counts <- counts[keep_genes, keep_cells, drop = FALSE]
gc()

normal_cells <- intersect(normal_cells, colnames(counts))
log_msg("Normal cells matching filtered counts matrix: ", length(normal_cells))

if (length(normal_cells) == 0) {
  log_msg("WARNING: no normal_cells remain after filtering. CopyKAT will infer references automatically.")
}

total_cells <- ncol(counts)
log_msg("Total cells for CopyKAT: ", total_cells)

cell_names <- sort(colnames(counts))

if (total_cells <= target_batch_size) {
  batches <- list(cell_names)
} else {
  batch_ids <- ceiling(seq_along(cell_names) / target_batch_size)
  batches <- split(cell_names, batch_ids)
}

log_msg("Number of batches: ", length(batches))

batch_info <- data.frame(
  batch = seq_along(batches),
  cells = vapply(batches, length, integer(1)),
  normal_cells = vapply(batches, function(x) length(intersect(normal_cells, x)), integer(1)),
  stringsAsFactors = FALSE
)

write.csv(
  batch_info,
  file = file.path(output_dir, "batch_info.csv"),
  row.names = FALSE
)

save_rds_safe(
  list(
    target_batch_size = target_batch_size,
    min_genes_per_cell = min_genes_per_cell,
    min_cells_per_gene = min_cells_per_gene,
    min_normal_cells_per_batch = min_normal_cells_per_batch,
    genome = genome,
    id_type = id_type,
    distance = distance,
    n_cores = n_cores,
    total_cells = total_cells,
    total_genes = nrow(counts),
    batch_info = batch_info
  ),
  file.path(final_output_dir, "copykat_run_parameters.rds")
)

for (i in seq_along(batches)) {
  batch_cells <- batches[[i]]
  batch_dir <- file.path(output_dir, paste0("batch_", i))

  done_file <- file.path(batch_dir, paste0("batch_", i, ".done"))
  failed_file <- file.path(batch_dir, paste0("batch_", i, ".failed"))
  result_file <- file.path(batch_dir, paste0("copykat_result_batch_", i, ".rds"))

  if (file.exists(done_file) && file.exists(result_file)) {
    log_msg("Skipping batch ", i, " because it is already complete")
    next
  }

  if (file.exists(failed_file)) {
    file.remove(failed_file)
  }

  log_msg("Preparing batch ", i, " of ", length(batches))

  batch_counts <- counts[, batch_cells, drop = FALSE]
  batch_normal_cells <- intersect(normal_cells, batch_cells)

  tryCatch(
    {
      result <- run_copykat_one_batch(
        batch_counts = batch_counts,
        batch_normal_cells = batch_normal_cells,
        batch_id = i,
        batch_dir = batch_dir
      )

      rm(result)
      gc()

      writeLines("done", done_file)
      log_msg("Finished batch ", i)
    },
    error = function(e) {
      log_msg("ERROR in batch ", i, ": ", conditionMessage(e))
      dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)
      writeLines(conditionMessage(e), failed_file)
    }
  )

  rm(batch_counts)
  gc()
}

log_msg("Finished batch loop")

############################################################
# Combine batch-level predictions
#
# After all batches finish successfully, prediction tables are
# merged into a single annotation file containing one row per cell.
#
# Final outputs are written to:
#   copykat_final/copykat_combined_predictions.csv
#   copykat_final/copykat_combined_predictions.rds
############################################################

prediction_files <- list.files(
  output_dir,
  pattern = "^copykat_prediction_batch_.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(prediction_files) > 0) {
  log_msg("Combining prediction CSV files: ", length(prediction_files))

  prediction_list <- lapply(prediction_files, function(f) {
    x <- read.csv(f, stringsAsFactors = FALSE)
    x$source_file <- basename(f)
    x
  })

  combined_predictions <- do.call(rbind, prediction_list)

  if (!"cell.names" %in% colnames(combined_predictions)) {
    stop("Combined predictions do not contain cell.names column.")
  }

  duplicate_predictions <- combined_predictions$cell.names[duplicated(combined_predictions$cell.names)]

  if (length(duplicate_predictions) > 0) {
    log_msg("WARNING: duplicated cell.names in combined predictions: ", length(unique(duplicate_predictions)))
  }

  predicted_cells <- unique(combined_predictions$cell.names)
  missing_predictions <- setdiff(cell_names, predicted_cells)
  extra_predictions <- setdiff(predicted_cells, cell_names)

  log_msg("Input cells after filtering: ", length(cell_names))
  log_msg("Unique predicted cells: ", length(predicted_cells))
  log_msg("Missing predictions: ", length(missing_predictions))
  log_msg("Extra predictions: ", length(extra_predictions))

  combined_prediction_csv <- file.path(final_output_dir, "copykat_combined_predictions.csv")
  combined_prediction_rds <- file.path(final_output_dir, "copykat_combined_predictions.rds")
  missing_prediction_csv <- file.path(final_output_dir, "copykat_missing_predictions.csv")
  extra_prediction_csv <- file.path(final_output_dir, "copykat_extra_predictions.csv")

  write.csv(combined_predictions, combined_prediction_csv, row.names = FALSE)
  save_rds_safe(combined_predictions, combined_prediction_rds)

  write.csv(
    data.frame(cell.names = missing_predictions),
    missing_prediction_csv,
    row.names = FALSE
  )

  write.csv(
    data.frame(cell.names = extra_predictions),
    extra_prediction_csv,
    row.names = FALSE
  )

  log_msg("Combined predictions saved: ", combined_prediction_csv)
} else {
  log_msg("WARNING: no prediction CSV files found to combine")
}

failed_files <- list.files(
  output_dir,
  pattern = "\\.failed$",
  recursive = TRUE,
  full.names = TRUE
)

done_files <- list.files(
  output_dir,
  pattern = "\\.done$",
  recursive = TRUE,
  full.names = TRUE
)

log_msg("Completed batches: ", length(done_files))
log_msg("Failed batches: ", length(failed_files))

if (length(failed_files) > 0) {
  log_msg("Failed batch files:")
  for (f in failed_files) {
    log_msg("  ", f)
  }
}

capture.output(
  sessionInfo(),
  file = file.path(final_output_dir, "sessionInfo_end.txt")
)

log_msg("CopyKAT batch script finished")
