# Figure 1b input (your methodology, all points): for EVERY directed pair (base feature, contributing
# feature): base 1D IG and base 1D IG + 2D IG added by the contributing feature (current method,
# I.lower = 1D IG), with no threshold. MDFS calls are made in the same order as run_mdfs.R, so with the
# same seed the IGs are identical to the revision rerun (verified: identical rows, IG difference 0).
#
# status: 0 = not significant (added IG <= Bonferroni pair threshold, or neither feature 2D-relevant)
#         1 = significant but removed by the true-synergy filter (total IG <= max 1D IG of the two)
#         2 = synergistic (kept by run_mdfs.R)
# Usage: Rscript all_directed.R X.tsv y.tsv seed out.tsv
suppressMessages({library(data.table); library(MDFS)})
args <- commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[3]))
X <- fread(args[1], sep = "\t"); X <- as.matrix(X[, -1, with = FALSE])
y <- fread(args[2], sep = "\t")[[2]]
fn <- colnames(X); n <- ncol(X)

result_1d <- MDFS(X, y, dimensions = 1)
ig_1d <- result_1d$statistic
result_2d <- MDFS(X, y, dimensions = 2)
sig <- result_2d$relevant.variables
thr <- qchisq(1 - 0.05 / max(1, n * (n - 1) / 2), df = 2) / (2 * log(2))

tp <- ComputeInterestingTuples(X, y, dimensions = 2, interesting.vars = seq_len(n), I.lower = ig_1d, ig.thr = 0)
base <- ifelse(tp$Var == tp$Tuple.1, tp$Tuple.2, tp$Tuple.1)
total <- ig_1d[base] + tp$IG
significant <- (tp$Var %in% sig | base %in% sig) & (tp$IG > thr)   # rows run_mdfs.R receives
status <- ifelse(!significant, 0L, ifelse(total > pmax(ig_1d[base], ig_1d[tp$Var]), 2L, 1L))
fwrite(data.table(base = fn[base], contributing = fn[tp$Var], base_1D = ig_1d[base], total = total, status = status),
       args[4], sep = "\t")
