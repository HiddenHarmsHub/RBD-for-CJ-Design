verbose <- FALSE
seed <- TRUE
start <- 10

rbd <- function(data, tol = 1e-6, col = min(10, ncol(data))) {

  if (missing(data)) {
    stop("RBD needs at least one input, a matrix")
  }

  # preparation of the algorithm
  nr <- nrow(data)
  nc <- ncol(data)
  bases <- matrix(0, nr, col)
  trans_mat <- matrix(0, col, nc)
  xi_flag <- integer(col)

  if (seed) {
    xi_flag[1] <- start
  } else {
    xi_flag[1] <- sample(nc, 1)
  }
  i <- 1
  cur_err <- tol + 1

  # Preparation for efficient error evaluation
  ftf <- colSums(t(data)^2)  # could we just use rowSums without the t()?
  ataatxi <- matrix(0, nc, col)

  # The RBD greedy algorithm
  while (i <= col && cur_err > tol) {
    if (verbose) {
      cat("we are on loop interation ", i, "\n")
      cat("starting column position is ", xi_flag[i], "\n")
    }
    bi_cand <- data[, xi_flag[i]]

    # Inside: Gram-Schmidt orthonormalization of the current candidate with all previously chosen basis vectors
    for (j in if (i > 1) 1:(i-1) else integer(0)) {
      bi_cand <- bi_cand - (t(bi_cand) %*% bases[, j]) * bases[, j]
    }
    normi <- sqrt(t(bi_cand) %*% bi_cand)
    if (verbose) {
      cat('normi is ', normi, '\n')
    }
    if (normi < 1e-7) {
      cat("Reduced system getting singular - to stop with ", i - 1, "basis functions\n")
      bases <- bases[, 1 : i - 1]
      trans_mat <- trans_mat[1 : i - 1, ]
      break
    } else {
      bases[, i] <- bi_cand / normi
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
    if (i < col) {
      xi_flag[i + 1] <- temp_pos
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
