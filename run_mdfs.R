# One MDFS run: synergistic feature pairs of a feature matrix X and binary labels y.
#
# Usage: Rscript run_mdfs.R <X.tsv> <y.tsv> <pairs_output.tsv> [seed]
#   X.tsv   tab-separated, first column = sample IDs, other columns = features
#   y.tsv   tab-separated with a header, columns sample_id and label
#   seed    integer seed of the MDFS discretization (default: 42)
#
# Called by mdfs_pairs.py; see README.md for the method.

suppressPackageStartupMessages({
  library(data.table)
  library(MDFS)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_mdfs.R <X.tsv> <y.tsv> <pairs_output.tsv> [seed]")
}
input_X_file <- args[1]
input_y_file <- args[2]
output_pairs <- args[3]
seed <- if (length(args) >= 4) as.integer(args[4]) else 42L
set.seed(seed)

# === Read X (drop the sample ID column) and y ===
X <- fread(input_X_file, sep = "\t", header = TRUE)
X <- X[, -1, with = FALSE]
y_dt <- fread(input_y_file, sep = "\t", header = TRUE)
if (ncol(y_dt) < 2) {
  stop("y file must have two columns (sample_id, label).")
}
y <- y_dt[[2]]

X_matrix <- as.matrix(X)
feature_names <- colnames(X)

# === Shared discretization ===
# MDFS discretizes each feature at random. The 1D IGs, the 2D significance test
# and ComputeInterestingTuples all use the same seed and the 2D range, so each
# feature is discretized identically in every step and
#   total_IG = base_IG_1D + IG_2D_added
# is the joint IG of the pair (identical in both directions).
range_2d <- GetRange(n = nrow(X_matrix), dimensions = 2)

# === 1D IG of each feature ===
ig_1d <- ComputeMaxInfoGains(X_matrix, y, dimensions = 1, seed = seed, range = range_2d)$IG
names(ig_1d) <- feature_names

# === MDFS 2D: features significant jointly with at least one partner ===
result_2d <- MDFS(X_matrix, y, dimensions = 2, seed = seed, range = range_2d)
sig_2d_vars <- result_2d$relevant.variables
cat("2D significant features:", length(sig_2d_vars), "\n")

out_df <- data.frame()
if (length(sig_2d_vars) > 0) {
  # Pair-level chi-squared cutoff, Bonferroni-corrected over all pairs.
  # MDFS reports IG in bits, undoubled; 2*log(2)*IG ~ chi^2(df). Default MDFS
  # parameters (divisions = 1, response.divisions = 1, dimensions = 2) give
  # df = response.divisions * divisions * (divisions + 1)^(dimensions - 1) = 2.
  alpha <- 0.05
  df_2d <- 1 * 1 * (1 + 1)^(2 - 1)
  n_total <- ncol(X_matrix)
  n_all_pairs <- max(1, n_total * (n_total - 1) / 2)
  ig_thr_pair <- qchisq(1 - alpha / n_all_pairs, df = df_2d) / (2 * log(2))
  cat(sprintf("Pair-level IG cutoff: %.4f bits (chi^2 df=%d, alpha=%.3g Bonferroni / %d pairs)\n",
              ig_thr_pair, df_2d, alpha, n_all_pairs))

  # IG of each tuple is the joint IG minus I.lower of the partner (base) variable
  tuples <- ComputeInterestingTuples(
    X_matrix, y,
    dimensions = 2,
    seed = seed,
    range = range_2d,
    interesting.vars = sig_2d_vars,
    I.lower = ig_1d,
    ig.thr = ig_thr_pair
  )

  if (!is.null(tuples) && nrow(tuples) > 0) {
    # Var = contributing variable, the other tuple member = base variable
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

    # True-synergy filter: the pair carries more information than either feature
    # alone. Removes 'one strong feature + noise partner' artifacts.
    n_before <- nrow(out_df)
    out_df <- out_df[out_df$total_IG > pmax(out_df$base_IG_1D, out_df$contributing_IG_1D), ]
    cat(sprintf("True-synergy filter: %d -> %d pair directions\n", n_before, nrow(out_df)))
  }
}

if (nrow(out_df) == 0) cat("No synergistic pairs found.\n")
write.table(out_df, output_pairs, row.names = FALSE, sep = "\t")
