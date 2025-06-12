verbose <- FALSE
seed <- FALSE
start <- 10

rbd <- function(data, tol = 1e-6, max_cols = min(10, ncol(data))) {

  if (missing(data)) {
    stop("RBD needs at least one input, a matrix")
  }
  # preparation of the algorithm
  col_count <- ncol(data)
  bases <- matrix(0, nrow(data), max_cols)
  trans_mat <- matrix(0, max_cols, col_count)
  # select our first column (either using the seed or at random)
  if (seed) {
    xi_flag <- start
  } else {
    xi_flag <- sample(col_count, 1)
  }
  cur_err <- tol + 1  # just to get the loop going

  # Preparation for efficient error evaluation
  ftf <- rowSums(t(data) * t(data))
  ataatxi <- matrix(0, col_count, max_cols)

  # The RBD greedy algorithm
  i <- 1
  while (i <= max_cols && cur_err > tol) {
    if (verbose) {
      cat("we are on loop interation ", i, "\n")
      cat("using column in position ", xi_flag, "\n")
    }
    bi_cand <- data[, xi_flag]

    # Inside: Gram-Schmidt orthonormalization of the current candidate with all previously chosen basis vectors
    for (j in if (i > 1) 1:(i-1) else integer(0)) {
      bi_cand <- bi_cand - as.vector(t(bi_cand) %*% bases[, j, drop = FALSE]) * bases[, j]
    }
    normi <- sqrt(t(bi_cand) %*% bi_cand)
    if (verbose) {
      cat('normi is ', normi, '\n')
    }
    if (normi < 1e-7) { # should this be tol?? it seems to be in paper write up
      cat("Reduced system getting singular - to stop with ", i - 1, "basis functions\n")
      bases <- bases[, 1 : i - 1]
      trans_mat <- trans_mat[1 : i - 1, ]
      break
    } else {
      bases[, i] <- bi_cand / as.vector(normi)
    }
    trans_mat[i, ] <- t(bases[, i]) %*% data
    ataatxi[, i] <- t(data) %*% bases[, i]

    # Inside: Efficiently go through all the columns to identify where the error would be the largest if we were to use the current space for compression.
    tmm <- trans_mat[1 : i, , drop=FALSE]
    te1 <- rowSums(ataatxi[, 1 : i] * t(tmm))
    tm <- t(tmm)
    te2 <- rowSums(tm * t(tmm))
    errord <- ftf - 2 * te1 + te2
    cur_err <- max(errord)
    temp_pos <- which.max(errord)

    # Mark this location for the next round
    if (i < max_cols) {
      xi_flag <- temp_pos
    }
    if (verbose) {
      cat("cur_err is ", cur_err, "\n")
    }
    if (cur_err <= tol) {
      cat(sprintf("Reduced system getting accurate enough - to stop with %d basis functions\n", i))
      bases <- bases[, 1 : i]
      trans_mat <- trans_mat[1 : i, ]
    } else {
      i <- i + 1
    }

  }
  return(list(bases = bases, trans_mat = trans_mat))
}

test_df <- read.csv("data/C_Barnet.csv", header = TRUE)
test_matrix <- as.matrix(test_df)
result <- rbd(test_matrix, 1e-13, 22)
