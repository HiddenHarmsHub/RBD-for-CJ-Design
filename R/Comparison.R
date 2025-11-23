## svd_on_ECET_example.R
## Example: compute RBD and brute-force design probabilities, compare with KL
## Run with: Rscript R/Comparison.R

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

## Construct a small SPD test covariance matrix C (5x5)
C <- matrix(c(
  2, 0.3, 0.1, 0.0, 0.2,
  0.3, 2,   0.4, 0.1, 0.0,
  0.1, 0.4, 2,   0.2, 0.1,
  0.0, 0.1, 0.2, 2,   0.3,
  0.2, 0.0, 0.1, 0.3, 2
), nrow = 5, byrow = TRUE)
N <- nrow(C)

## Compute design probabilities via RBD
design_probs_rbd <- compute_design_probs_rbd(C)

## Compute brute-force design probabilities
B_brute <- bruteforceB(C)
design_probs_brute <- compute_design_probs(N, B_brute)

## Compute KL divergence (brute || rbd)
kl_divergence <- as.numeric(philentropy::KL(rbind(t(design_probs_brute), t(design_probs_rbd))))
cat("KL divergence (brute || rbd):", format(kl_divergence, digits = 6), "\n")

