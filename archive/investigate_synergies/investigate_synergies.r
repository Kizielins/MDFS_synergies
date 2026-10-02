#!/usr/bin/env Rscript
#
# investigate_synergies.r
#
#Rscript investigate_synergies.r \
#  --prepared_data_dir prepared_data \
#  --data_type combined \
#  --pairs pairs.csv \
# --out ig_synergy_combined.csv
#
# Goal:
#   For each feature pair, compute the 2D information gain (IG) and compare it to the
#   1D IG of the base feature. The 2D IG from ComputeMaxInfoGains is ADDITIONAL IG —
#   it represents how much the contributing variable improves the base feature's IG.
#
#   For a pair (f1, f2), we report BOTH directions:
#     - base=f1, contributing=f2: total_IG = f1_IG_1D + IG_2D_added(f2→f1)
#     - base=f2, contributing=f1: total_IG = f2_IG_1D + IG_2D_added(f1→f2)
#
#   Optional: 1D IG for derived synthetic transforms (LR / MULT / GM)
#
# Based on: run_mdfs.r (same data reading conventions).
#
# Usage (single cohort X/y):
#   Rscript investigate_synergies.r <X_tsv> <y_tsv> <pairs_tsv_or_csv> <out_csv>
#
# Usage (MERGE ALL cohorts from prepared_data):
#   Rscript investigate_synergies.r --prepared_data_dir prepared_data --data_type <taxa|functions|combined> --pairs <pairs.csv> --out <out.csv>
#
# Expected inputs:
#   - X_tsv: tab-delimited table with header; first column may be sample IDs (will be dropped)
#   - y_tsv: tab-delimited table; must contain labels in column 2 by default (index + labels)
#   - pairs file: CSV/TSV with columns: (Feature1, Feature2) or (base, contributing)
#       Optional columns preserved in output: Data_Type, Feature_Type, Feature (or any extra cols)
#
# Output:
#   - out_csv: TWO rows per input pair (one per direction) with columns:
#       base, contributing, base_IG_1D, contributing_IG_1D, IG_2D_added, total_IG
#       + optional 1D IG for LR/MULT/GM derived features (if compute_synthetic_1d = TRUE)
#

suppressPackageStartupMessages({
  library(data.table)
  library(MDFS)
})

args <- commandArgs(trailingOnly = TRUE)

parse_flag_args <- function(argv) {
  # Very small flag parser: --key value
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (startsWith(a, "--")) {
      key <- substring(a, 3)
      if (i == length(argv)) stop(paste0("Missing value for flag: ", a))
      out[[key]] <- argv[[i + 1]]
      i <- i + 2
    } else {
      out[["__positional__"]] <- c(out[["__positional__"]], a)
      i <- i + 1
    }
  }
  out
}

opts <- parse_flag_args(args)

# Backwards compatible positional mode:
#   <X> <y> <pairs> <out>
pos <- opts[["__positional__"]]

use_merged <- !is.null(opts[["prepared_data_dir"]]) || !is.null(opts[["data_type"]])

if (!use_merged) {
  if (is.null(pos) || length(pos) < 4) {
    stop("Usage: Rscript investigate_synergies.r <X_tsv> <y_tsv> <pairs_csv_or_tsv> <out_csv>\n",
         "   or: Rscript investigate_synergies.r --prepared_data_dir prepared_data --data_type <taxa|functions|combined> --pairs <pairs.csv> --out <out.csv>")
  }
  input_X_file <- pos[[1]]
  input_y_file <- pos[[2]]
  pairs_file   <- pos[[3]]
  output_csv   <- pos[[4]]
} else {
  prepared_data_dir <- opts[["prepared_data_dir"]]
  data_type <- opts[["data_type"]]
  pairs_file <- opts[["pairs"]]
  output_csv <- opts[["out"]]
  if (is.null(prepared_data_dir) || is.null(data_type) || is.null(pairs_file) || is.null(output_csv)) {
    stop("Merged mode requires: --prepared_data_dir, --data_type, --pairs, --out")
  }
}

# ---- knobs (edit as needed) ----
drop_first_col_if_id <- TRUE
y_col_index <- 2              # default: y is in 2nd column (index + labels)
pairs_sep_guess <- TRUE       # auto-detect CSV vs TSV based on extension
compute_synthetic_1d <- TRUE  # also compute 1D IG of LR/MULT/GM derived features
pseudocount <- 1e-9           # for LR/log operations and products

load_single_X_y <- function(x_path, y_path, cohort_prefix = NULL) {
  X_raw <- fread(x_path, sep = "\t", data.table = FALSE)
  if (ncol(X_raw) < 2) stop(paste0("X has too few columns: ", x_path))

  id_col_x <- colnames(X_raw)[1]
  ids_x <- as.character(X_raw[[id_col_x]])

  X_feat <- X_raw[, -1, drop = FALSE]

  y_raw <- fread(y_path, sep = "\t", data.table = FALSE)
  if (ncol(y_raw) < y_col_index) {
    stop(paste0("y file must have at least ", y_col_index, " columns (index + labels): ", y_path))
  }
  id_col_y <- colnames(y_raw)[1]
  ids_y <- as.character(y_raw[[id_col_y]])
  y_vec <- y_raw[[y_col_index]]

  # Ensure unique IDs across cohorts (avoids accidental collisions when merging)
  if (!is.null(cohort_prefix)) {
    ids_x <- paste0(cohort_prefix, "__", ids_x)
    ids_y <- paste0(cohort_prefix, "__", ids_y)
  }

  # Align y to X row order
  m <- match(ids_x, ids_y)
  keep <- !is.na(m)
  if (sum(!keep) > 0) {
    warning(paste0("Dropping ", sum(!keep), " rows from X because no matching y label found (", x_path, ")"))
  }
  ids_x <- ids_x[keep]
  X_feat <- X_feat[keep, , drop = FALSE]
  y_aligned <- y_vec[m[keep]]

  # Return as data.table-like with ID column retained (for merging)
  X_feat[[id_col_x]] <- ids_x
  # Put ID column first
  X_feat <- X_feat[, c(id_col_x, setdiff(colnames(X_feat), id_col_x)), drop = FALSE]

  list(X = X_feat, y = y_aligned, id_col = id_col_x)
}

if (!use_merged) {
  message("Loading single-cohort X: ", input_X_file)
  message("Loading single-cohort y: ", input_y_file)

  loaded <- load_single_X_y(input_X_file, input_y_file, cohort_prefix = NULL)
  X_dt <- loaded$X
  y <- loaded$y
  id_col <- loaded$id_col

  # Drop ID column for MDFS matrix
  if (drop_first_col_if_id && ncol(X_dt) > 1) {
    X_dt <- X_dt[, -1, drop = FALSE]
  }

  X_mat <- as.matrix(X_dt)
  feature_names <- colnames(X_dt)
} else {
  message("Loading MERGED cohorts from: ", prepared_data_dir, " | data_type=", data_type)

  x_paths <- list.files(prepared_data_dir, pattern = paste0("^X_.*_", data_type, "\\.tsv$"), full.names = TRUE)
  if (length(x_paths) == 0) stop(paste0("No X files found in ", prepared_data_dir, " for data_type=", data_type))

  cohorts <- sub(paste0("^X_(.*)_", data_type, "\\.tsv$"), "\\1", basename(x_paths))
  y_paths <- file.path(prepared_data_dir, paste0("y_", cohorts, "_", data_type, ".tsv"))

  missing_y <- !file.exists(y_paths)
  if (any(missing_y)) {
    stop(paste0("Missing y files for cohorts: ", paste(cohorts[missing_y], collapse = ", ")))
  }

  loaded_list <- vector("list", length(x_paths))
  for (k in seq_along(x_paths)) {
    loaded_list[[k]] <- load_single_X_y(x_paths[[k]], y_paths[[k]], cohort_prefix = cohorts[[k]])
  }

  # Union columns across cohorts, fill missing with NA then set to 0
  id_col <- loaded_list[[1]]$id_col
  X_list <- lapply(loaded_list, function(z) z$X)
  y_list <- lapply(loaded_list, function(z) z$y)

  X_merged <- data.table::rbindlist(X_list, fill = TRUE)
  # Replace NAs in feature columns with 0
  feat_cols <- setdiff(colnames(X_merged), id_col)
  for (c in feat_cols) {
    if (anyNA(X_merged[[c]])) X_merged[[c]][is.na(X_merged[[c]])] <- 0
  }

  y <- unlist(y_list, use.names = FALSE)

  # Drop ID col for matrix
  X_dt <- as.data.frame(X_merged)
  X_dt <- X_dt[, feat_cols, drop = FALSE]
  X_mat <- as.matrix(X_dt)
  feature_names <- colnames(X_dt)

  message("Merged samples: ", nrow(X_dt), " | features (union): ", ncol(X_dt))
}

# ---- load pairs ----
guess_sep <- function(path) {
  if (!pairs_sep_guess) return("\t")
  if (grepl("\\.csv$", path, ignore.case = TRUE)) return(",")
  return("\t")
}

pairs_sep <- guess_sep(pairs_file)
message("Loading pairs: ", pairs_file, " (sep='", pairs_sep, "')")
pairs_dt <- fread(pairs_file, sep = pairs_sep, data.table = FALSE)

# Accept both old (Feature1, Feature2) and new (base, contributing) column names
if ("base" %in% colnames(pairs_dt) && "contributing" %in% colnames(pairs_dt)) {
  pair_col1 <- "base"
  pair_col2 <- "contributing"
} else if ("Feature1" %in% colnames(pairs_dt) && "Feature2" %in% colnames(pairs_dt)) {
  pair_col1 <- "Feature1"
  pair_col2 <- "Feature2"
} else {
  stop("Pairs file must have columns (base, contributing) or (Feature1, Feature2)")
}

# If the pairs file includes a data-type column, filter it automatically in merged mode.
# Accept common variants: Data_Type, Data Type, data_type
if (use_merged) {
  dt_col <- NULL
  for (cand in c("Data_Type", "Data Type", "data_type")) {
    if (cand %in% colnames(pairs_dt)) {
      dt_col <- cand
      break
    }
  }
  if (!is.null(dt_col)) {
    before_n <- nrow(pairs_dt)
    pairs_dt <- pairs_dt[as.character(pairs_dt[[dt_col]]) == data_type, , drop = FALSE]
    message("Filtered pairs by ", dt_col, " == ", data_type, ": ", before_n, " -> ", nrow(pairs_dt))
  }
}

# ---- helper: compute 1D IG for one variable ----
compute_ig_1d <- function(vec, y_vec) {
  m <- as.matrix(data.frame(v = as.numeric(vec)))
  res <- tryCatch(
    ComputeMaxInfoGains(m, y_vec, dimensions = 1, return.tuples = FALSE),
    error = function(e) ComputeMaxInfoGains(m, y_vec, dimensions = 1)
  )
  res_df <- as.data.frame(res)
  if (nrow(res_df) == 0) return(NA_real_)
  as.numeric(res_df$IG[1])
}

# ---- helper: compute 2D IG for a pair of variables ----
# Returns BOTH directions: IG_2D_added for each variable (contributed by its partner).
# ComputeMaxInfoGains(dimensions=2) on 2 variables returns 2 rows:
#   Row 1 = IG_2D_added for variable 1 (v2 contributes to v1)
#   Row 2 = IG_2D_added for variable 2 (v1 contributes to v2)
compute_ig_2d_pair <- function(vec1, vec2, y_vec) {
  m <- as.matrix(data.frame(v1 = as.numeric(vec1), v2 = as.numeric(vec2)))
  res <- ComputeMaxInfoGains(m, y_vec, dimensions = 2, return.tuples = TRUE)
  res_df <- as.data.frame(res)
  if (nrow(res_df) == 0) return(list(ig_added_v1 = NA_real_, ig_added_v2 = NA_real_))
  # Row 1: IG_2D_added for v1 (v2 contributes to v1)
  # Row 2: IG_2D_added for v2 (v1 contributes to v2)
  list(
    ig_added_v1 = as.numeric(res_df$IG[1]),
    ig_added_v2 = as.numeric(res_df$IG[2])
  )
}

# ---- helper: synthetic transforms ----
safe_lr <- function(x1, x2, pc) {
  log((x1 + pc) / (x2 + pc))
}
safe_mult <- function(x1, x2, pc) {
  (x1 + pc) * (x2 + pc)
}
safe_gm <- function(x1, x2, pc) {
  sqrt((x1 + pc) * (x2 + pc))
}

# ---- main loop ----
# For each input pair, produce TWO rows (one per direction):
#   Direction 1: base=f1, contributing=f2 → total_IG = f1_IG_1D + IG_2D_added(f2→f1)
#   Direction 2: base=f2, contributing=f1 → total_IG = f2_IG_1D + IG_2D_added(f1→f2)
idx_map <- setNames(seq_along(feature_names), feature_names)

out_rows <- vector("list", nrow(pairs_dt) * 2)

# Preserve extra columns from pairs file (Data_Type, Feature_Type, etc.)
extra_cols <- setdiff(colnames(pairs_dt), c(pair_col1, pair_col2))

for (i in seq_len(nrow(pairs_dt))) {
  f1 <- as.character(pairs_dt[[pair_col1]][i])
  f2 <- as.character(pairs_dt[[pair_col2]][i])

  # Build base row with extra columns
  extra_vals <- as.list(pairs_dt[i, extra_cols, drop = FALSE])

  if (!(f1 %in% names(idx_map)) || !(f2 %in% names(idx_map))) {
    # Missing features: produce NA rows for both directions
    for (d in 1:2) {
      row <- data.frame(
        base = if (d == 1) f1 else f2,
        contributing = if (d == 1) f2 else f1,
        base_IG_1D = NA_real_,
        contributing_IG_1D = NA_real_,
        IG_2D_added = NA_real_,
        total_IG = NA_real_,
        stringsAsFactors = FALSE
      )
      if (compute_synthetic_1d) {
        row$IG_1D_LR <- NA_real_
        row$IG_1D_MULT <- NA_real_
        row$IG_1D_GM <- NA_real_
      }
      if (length(extra_cols) > 0) row <- cbind(row, extra_vals)
      out_rows[[(i - 1) * 2 + d]] <- row
    }
    next
  }

  v1 <- X_mat[, idx_map[[f1]]]
  v2 <- X_mat[, idx_map[[f2]]]

  ig1 <- compute_ig_1d(v1, y)
  ig2 <- compute_ig_1d(v2, y)
  igp <- compute_ig_2d_pair(v1, v2, y)

  # Synthetic 1D IGs (same for both directions)
  ig_lr <- NA_real_; ig_mult <- NA_real_; ig_gm <- NA_real_
  if (compute_synthetic_1d) {
    ig_lr <- compute_ig_1d(safe_lr(v1, v2, pseudocount), y)
    ig_mult <- compute_ig_1d(safe_mult(v1, v2, pseudocount), y)
    ig_gm <- compute_ig_1d(safe_gm(v1, v2, pseudocount), y)
  }

  # Direction 1: base=f1, contributing=f2
  row1 <- data.frame(
    base = f1,
    contributing = f2,
    base_IG_1D = ig1,
    contributing_IG_1D = ig2,
    IG_2D_added = igp$ig_added_v1,
    total_IG = ig1 + igp$ig_added_v1,
    stringsAsFactors = FALSE
  )
  if (compute_synthetic_1d) {
    row1$IG_1D_LR <- ig_lr
    row1$IG_1D_MULT <- ig_mult
    row1$IG_1D_GM <- ig_gm
  }
  if (length(extra_cols) > 0) row1 <- cbind(row1, extra_vals)

  # Direction 2: base=f2, contributing=f1
  row2 <- data.frame(
    base = f2,
    contributing = f1,
    base_IG_1D = ig2,
    contributing_IG_1D = ig1,
    IG_2D_added = igp$ig_added_v2,
    total_IG = ig2 + igp$ig_added_v2,
    stringsAsFactors = FALSE
  )
  if (compute_synthetic_1d) {
    row2$IG_1D_LR <- ig_lr
    row2$IG_1D_MULT <- ig_mult
    row2$IG_1D_GM <- ig_gm
  }
  if (length(extra_cols) > 0) row2 <- cbind(row2, extra_vals)

  out_rows[[(i - 1) * 2 + 1]] <- row1
  out_rows[[(i - 1) * 2 + 2]] <- row2
}

out_df <- do.call(rbind, out_rows)

# Add derived columns
out_df$ig_gain <- out_df$IG_2D_added
denom <- ifelse(out_df$base_IG_1D == 0, NA_real_, out_df$base_IG_1D)
out_df$ig_gain_pct <- (out_df$ig_gain / denom) * 100

# Sort by IG gain descending
out_df <- out_df[order(-out_df$IG_2D_added, -out_df$total_IG), , drop = FALSE]

message("Writing output: ", output_csv)
write.table(out_df, file = output_csv, sep = ",", row.names = FALSE, quote = TRUE)

message("Done. Input pairs: ", nrow(pairs_dt), " | Output rows (both directions): ", nrow(out_df))

