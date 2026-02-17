## Construct brute-force B matrix from covariance C
##
## Given an N x N covariance matrix C, construct the matrix B indexed by
## unordered node pairs (i,j) with i < j. Entry B[p, q] corresponds to
## C[i_p, i_q] - C[i_p, j_q] - C[j_p, i_q] + C[j_p, j_q] where (i_p, j_p)
## and (i_q, j_q) are the node pairs for rows p and q respectively.
##
## @param C Numeric matrix, symmetric covariance (N x N).
## @return A matrix of size choose(N, 2) x choose(N, 2).
## @examples
## C <- diag(3)
## bruteforceB(C)
bruteforceB <- function(C) {
  if (missing(C)) stop("C is required")
  if (!is.matrix(C)) stop("C must be a matrix")
  N <- nrow(C)
  if (ncol(C) != N) stop("C must be square")
  if (N < 2) return(matrix(numeric(0), 0, 0))

  pairs <- utils::combn(N, 2)
  m <- ncol(pairs)
  i1 <- pairs[1, ]
  i2 <- pairs[2, ]

  B <- matrix(NA_real_, nrow = m, ncol = m)

  for (p in seq_len(m)) {
    for (q in seq_len(m)) {
      B[p, q] <- C[i1[p], i1[q]] - C[i1[p], i2[q]] - C[i2[p], i1[q]] + C[i2[p], i2[q]]
    }
  }

  B
}

## Compute design probabilities from B (brute-force method)
##
## Compute the design probabilities by spectral decomposition of B.
## The returned vector has length choose(N,2) and sums to 1 (up to numerical error).
##
## @param B Numeric matrix produced by `bruteforceB`.
## @param N Integer number of nodes; inferred from `nrow(B)` if not supplied.
## @param item_labels Optional character vector of length N with item names for labelling output.
## @return Numeric vector of length choose(N,2) with design probabilities.
## @examples
## C <- diag(4)
## B <- bruteforceB(C)
## compute_design_probs(B)
compute_design_probs <- function(B, N = NULL, item_labels = NULL) {

  if (is.null(N)) N <- (1 + sqrt(1 + 8 * nrow(B))) / 2
  N <- as.integer(round(N))

  eig <- eigen(B)
  values <- eig$values
  vectors <- eig$vectors

  weights <- sum(values)

  probs <- (vectors^2) %*% values
  probs <- as.numeric(probs / weights)

  if (is.null(item_labels)) item_labels <- seq_len(N)
  pairs <- combn(N, 2)
  names(probs) <- paste0(item_labels[pairs[1, ]], "-", item_labels[pairs[2, ]])

  probs
}
