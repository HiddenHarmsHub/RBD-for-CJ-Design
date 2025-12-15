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
## @param N Integer number of nodes.
## @param B Numeric matrix produced by `bruteforceB`.
## @return Numeric vector of length choose(N,2) with design probabilities.
## @examples
## C <- diag(4)
## B <- bruteforceB(C)
## compute_design_probs(4, B)
compute_design_probs <- function(N, B) {

  eig <- eigen(B)
  values <- eig$values
  vectors <- eig$vectors

  weights <- sum(values)

  probs <- (vectors^2) %*% values
  probs <- as.numeric(probs / weights)

  probs
}
