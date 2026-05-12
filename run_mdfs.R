# Load necessary libraries
library(data.table)
library(MDFS)

# === STEP 0: Get command line arguments ===
# Usage: Rscript script.R Train_set.txt Train_y.txt Train_mdfs_2d.csv Train_mdfs_1d_2d_unique.txt Train_mdfs_1d_ig.csv [seed] [no_bonferroni]
# Optional 6th argument: random seed (integer). If omitted, uses 12 (original behaviour).
# Optional 7th argument: "no_bonferroni" to skip the strict pair-level IG threshold
#   (used by apply_framework for exploratory synergy enumeration).
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript script.R <Train_set.txt> <Train_y.txt> <2D_output.csv> <1D_2D_unique.txt> <1D_IG_output.csv> [seed]")
}

seed <- if (length(args) >= 6) as.numeric(args[6]) else 12
set.seed(seed)

input_X_file  <- args[1]  # Train_set.txt
input_y_file  <- args[2]  # Train_y.txt
output_csv    <- args[3]  # Train_mdfs_2d.csv
output_unique <- args[4]  # Train_mdfs_1d_2d_unique.txt
output_1d_ig  <- args[5]  # Train_mdfs_1d_ig.csv

# === STEP 1: Read the features (X) and labels (y) files ===
X <- fread(input_X_file, sep = "\t")
X <- X[, -1, with = FALSE]  # drop first column if it's ID

y_dt <- fread(input_y_file, sep = "\t")
if (ncol(y_dt) >= 2) {
  y <- y_dt[[2]]
} else {
  stop("Error: y file must have at least two columns (index + labels).")
}

# === STEP 2: Convert to matrix ===
X_matrix <- as.matrix(X)
feature_names <- colnames(X)

# === STEP 3: MDFS 1D — significance + IGs ===
result_1d <- MDFS(X_matrix, y, dimensions = 1)

ig_1d <- result_1d$statistic  # one IG per variable (baseline)
names(ig_1d) <- feature_names

ig_1d_df <- data.frame(
  Feature = feature_names,
  IG = ig_1d,
  p_value = result_1d$p.value,
  adjusted_p_value = result_1d$adjusted.p.value,
  stringsAsFactors = FALSE
)
write.table(ig_1d_df, output_1d_ig, row.names = FALSE, sep = "\t")

relevant_1d_features <- feature_names[result_1d$relevant.variables]
cat("1D significant features:", length(relevant_1d_features), "\n")

# === STEP 4: MDFS 2D — significance ===
result_2d <- MDFS(X_matrix, y, dimensions = 2)

sig_2d_vars <- result_2d$relevant.variables
relevant_2d_features <- feature_names[sig_2d_vars]
cat("2D significant features:", length(relevant_2d_features), "\n")

# === STEP 5: ComputeInterestingTuples for significant 2D variables ===
# Pair-level chi-squared cutoff (Bonferroni-corrected over the number of pairs
# tested). MDFS reports IG in bits, undoubled; ComputePValue internally uses
# 2*log(2)*IG ~ chi^2(df). Default MDFS params (divisions=1, response.divisions=1,
# dimensions=2) give df = response.divisions * divisions * (divisions+1)^(dim-1) = 2.
#
# Optional 7th argument: "no_bonferroni" — skip the Bonferroni ig.thr and use
# default (0) instead. Used by apply_framework for exploratory synergy enumeration,
# while the CRC pipeline (framework.ipynb) keeps the strict threshold.
use_bonferroni <- TRUE
if (length(args) >= 7 && tolower(args[7]) == "no_bonferroni") {
  use_bonferroni <- FALSE
}

if (length(sig_2d_vars) > 0) {
  ig_thr_pair <- 0
  if (use_bonferroni) {
    alpha <- 0.05
    df_2d <- 1 * 1 * (1 + 1)^(2 - 1)
    n_total <- ncol(X_matrix)
    n_all_pairs <- max(1, n_total * (n_total - 1) / 2)
    alpha_corr <- alpha / n_all_pairs
    ig_thr_pair <- qchisq(1 - alpha_corr, df = df_2d) / (2 * log(2))
    cat(sprintf(
      "Pair-level IG cutoff: %.4f bits (chi^2 df=%d, alpha=%.3g Bonferroni / %d total pairs)\n",
      ig_thr_pair, df_2d, alpha, n_all_pairs))
  } else {
    cat("Bonferroni pair-level threshold SKIPPED (exploratory mode)\n")
  }

  tuples <- ComputeInterestingTuples(
    X_matrix, y,
    dimensions = 2,
    interesting.vars = sig_2d_vars,
    I.lower = ig_1d,
    ig.thr = ig_thr_pair
  )
  cat(sprintf("ComputeInterestingTuples returned %d rows for %d significant 2D variables\n",
              if (!is.null(tuples)) nrow(tuples) else 0L, length(sig_2d_vars)))

  if (!is.null(tuples) && nrow(tuples) > 0) {
    # Var = contributing variable, partner = base variable
    base_idx <- ifelse(tuples$Var == tuples$Tuple.1, tuples$Tuple.2, tuples$Tuple.1)

    out_df <- data.frame(
      base = feature_names[base_idx],
      contributing = feature_names[tuples$Var],
      base_IG_1D = ig_1d[base_idx],
      contributing_IG_1D = ig_1d[tuples$Var],
      IG_2D_added = tuples$IG,
      total_IG = ig_1d[base_idx] + tuples$IG,
      stringsAsFactors = FALSE
    )

    # True-synergy filter (symmetric): keep pairs where the joint 2D IG exceeds
    # the larger of the two 1D IGs — i.e. the pair carries more info than either
    # variable alone. Removes 'one strong + noise partner' artifacts.
    n_before <- nrow(out_df)
    max_1d <- pmax(out_df$base_IG_1D, out_df$contributing_IG_1D)
    out_df <- out_df[out_df$total_IG > max_1d, ]
    cat(sprintf("True-synergy filter (total_IG > max(base_1D, contrib_1D)): %d -> %d rows\n",
                n_before, nrow(out_df)))

    write.table(out_df, output_csv, row.names = FALSE, sep = "\t")
    cat("2D tuples written:", nrow(out_df), "rows\n")
  } else {
    cat("No interesting 2D tuples found.\n")
    write.table(data.frame(), output_csv, row.names = FALSE, sep = "\t")
  }
} else {
  cat("No significant 2D features — skipping ComputeInterestingTuples.\n")
  write.table(data.frame(), output_csv, row.names = FALSE, sep = "\t")
}

# === STEP 6: Unique features list (union of 1D + 2D significant) ===
unique_features <- unique(c(relevant_1d_features, relevant_2d_features))

if (length(unique_features) > 0) {
  writeLines(unique_features, output_unique)
  cat("Unique features saved:", length(unique_features), "to", output_unique, "\n")
} else {
  cat("No significant features found in 1D or 2D.\n")
}
