# Load necessary libraries
library(data.table)
library(MDFS)

# === STEP 0: Get command line arguments ===
# Usage: Rscript run_mdfs.R <X.tsv> <y.tsv> <2D_output.tsv> <1D_2D_unique.txt> <1D_IG_output.tsv> [seed]
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript run_mdfs.R <X.tsv> <y.tsv> <2D_output.tsv> <1D_2D_unique.txt> <1D_IG_output.tsv> [seed]")
}

seed <- if (length(args) >= 6) as.numeric(args[6]) else 12
set.seed(seed)

input_X_file  <- args[1]
input_y_file  <- args[2]
output_csv    <- args[3]
output_unique <- args[4]
output_1d_ig  <- args[5]

# === STEP 1: Read features (X) and labels (y) ===
X <- fread(input_X_file, sep = "\t")
X <- X[, -1, with = FALSE]  # drop sample ID column

y_dt <- fread(input_y_file, sep = "\t")
if (ncol(y_dt) >= 2) {
  y <- y_dt[[2]]
} else {
  stop("Error: y file must have at least two columns (index + labels).")
}

# === STEP 2: Convert to matrix ===
X_matrix <- as.matrix(X)

# === STEP 3: Compute 2D information gains ===
result <- ComputeMaxInfoGains(X_matrix, y, dimensions = 2, return.tuples = TRUE)

feature_names <- colnames(X)
result_df <- as.data.frame(result)
result_df$Feature1 <- feature_names[result_df$Tuple.1]
result_df$Feature2 <- feature_names[result_df$Tuple.2]
result_df <- result_df[, c("Feature1", "Feature2", "IG", "Discretization.nr")]

write.table(result_df, output_csv, row.names = FALSE, sep = "\t")

# === STEP 4: Compute 1D information gains ===
result_1d <- ComputeMaxInfoGains(X_matrix, y, dimensions = 1, return.tuples = FALSE)

if (is.list(result_1d)) {
  ig_values <- sapply(result_1d, function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    else if (is.list(x) && "IG" %in% names(x)) return(as.numeric(x$IG))
    else return(NA)
  })
  if (length(ig_values) == 0 || all(is.na(ig_values))) {
    ig_values <- as.numeric(unlist(result_1d))
  }
} else if (is.data.frame(result_1d)) {
  if ("IG" %in% colnames(result_1d)) {
    ig_values <- as.numeric(result_1d$IG)
  } else {
    numeric_cols <- sapply(result_1d, is.numeric)
    if (any(numeric_cols)) {
      ig_values <- as.numeric(result_1d[, which(numeric_cols)[1]])
    } else {
      ig_values <- as.numeric(unlist(result_1d))
    }
  }
} else {
  ig_values <- as.numeric(as.vector(result_1d))
}

if (length(ig_values) != length(feature_names)) {
  if (length(ig_values) > length(feature_names)) {
    ig_values <- ig_values[1:length(feature_names)]
  } else {
    ig_values <- c(ig_values, rep(NA, length(feature_names) - length(ig_values)))
  }
}

ig_1d_df <- data.frame(Feature = feature_names, IG = ig_values, stringsAsFactors = FALSE)
write.table(ig_1d_df, output_1d_ig, row.names = FALSE, sep = "\t")

# === STEP 5: Find relevant features (1D + 2D MDFS) and write unique list ===
result1 <- MDFS(X_matrix, y, dimensions = 1)
result2 <- MDFS(X_matrix, y, dimensions = 2)

relevant_1d <- feature_names[result1$relevant.variables]
relevant_2d <- feature_names[result2$relevant.variables]
unique_features <- unique(c(relevant_1d, relevant_2d))

if (length(unique_features) > 0) {
  writeLines(unique_features, output_unique)
} else {
  writeLines(character(0), output_unique)
}
