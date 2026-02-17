## Load required packages
if (!requireNamespace("expm", quietly = TRUE)) {
  stop("Package 'expm' is required for this example. Install with install.packages('expm')")
}
if (!requireNamespace("philentropy", quietly = TRUE)) {
  stop("Package 'philentropy' is required for this example. Install with install.packages('philentropy')")
}
if (!requireNamespace("speedyBBT", quietly = TRUE)) {
  stop("Package 'speedyBBT' is required for this example. Install with install.packages('speedyBBT')")
}


## Source project functions
source("R/rbd.R")
source("R/brute_force_functions.R")

library(speedyBBT)
library(expm)
library(philentropy)

## Read in Comparisons
tick <- Sys.time()
comparisons <- read.csv("data/Jones2017.csv", header = TRUE, check.names = FALSE)
pre_comparisons <- comparisons[1:700, ]  # appox first 20%
post_comparisons <- comparisons[-c(1:700), ]  # and approx last 80%

## Convert player names to numeric IDs (speedyBBTm requires numeric player IDs)
all_players <- unique(c(pre_comparisons$candidate_chosen, pre_comparisons$candidate_not_chosen))
player_ids <- setNames(seq_along(all_players), all_players)

## Fit BT Model
pre_model <- speedyBBTm(outcome = rep(0, length(pre_comparisons$candidate_chosen)),
                                  player1 = player_ids[pre_comparisons$candidate_chosen], 
                                  player2 = player_ids[pre_comparisons$candidate_not_chosen])

## Extract posterior covariance matrix for lambdas
C <- cov(pre_model$lambda)
colnames(C) <- rownames(C) <- all_players
tock <- Sys.time()
cat("Time taken to fit BT model and extract covariance matrix:", tock - tick, "\n")

## Compute design probabilities via RBD
tick <- Sys.time()
design_probs_rbd <- compute_design_probs_rbd(C)
tock <- Sys.time()
cat("Time taken for RBD computation:", tock - tick, "\n")


tick <- Sys.time()
B_brute <- bruteforceB(C)
design_probs_brute <- compute_design_probs(B_brute, item_labels = all_players)
tock <- Sys.time()
cat("Time taken for brute-force computation:", tock - tick, "\n")