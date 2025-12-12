## Load required packages
if (!requireNamespace("expm", quietly = TRUE)) {
  stop("Package 'expm' is required for this example. Install with install.packages('expm')")
}
if (!requireNamespace("philentropy", quietly = TRUE)) {
  stop("Package 'philentropy' is required for this example. Install with install.packages('philentropy')")
}

## Source project functions
source("R/rbd.R")
source("R/brute_force_functions.R")

## Create covariance matrix
# Read adjacency matrix and coerce to a 0/1 numeric matrix (not TRUE/FALSE)
# Keep original column names if present; allow logical, numeric or textual values
A_df <- read.csv("data/Seymour2022AdjacencyMatrix.csv", header = TRUE, check.names = FALSE)
A_mat <- as.matrix(A_df)
if (is.logical(A_mat) || is.numeric(A_mat)) {
  A <- matrix(as.integer(A_mat), nrow = nrow(A_mat), ncol = ncol(A_mat),
              dimnames = dimnames(A_mat))
} else {
  A <- matrix(0L, nrow = nrow(A_mat), ncol = ncol(A_mat), dimnames = dimnames(A_mat))
  true_idx <- which(toupper(as.character(A_mat)) %in% c("TRUE", "T", "1"))
  if (length(true_idx)) A[true_idx] <- 1L
}
A <- A[, -1] #Remove id column

C <- expm::expm(A)  # Exponentiate adjacency matrix to get covariance matrix
C <- diag(diag(C)^{-0.5}) %*% C %*% diag(diag(C)^{-0.5})  # Normalize to correlation matrix



## Compute design probabilities via RBD
tick <- Sys.time()
design_probs_rbd <- compute_design_probs_rbd(C)
tock <- Sys.time()
cat("Time taken for RBD computation:", tock - tick, "\n")
