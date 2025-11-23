verbose <- FALSE
seed <- FALSE
start <- 10

block_function <- function(i, N) {
  I <- diag(N - i)
  Z <- matrix(0, N-i, i-1)
  ones <- matrix(1, N-i, 1)
  cbind(Z, ones, -1*I)
}

transformation_matrix <- function(N) {
  E <- block_function(1, N)
  for(i in 2:N){
    ith_block <- block_function(i, N)
    E <- rbind(E, ith_block)
  }
  E
}

rbd <- function(data, tol = NULL, max_cols = NULL) {
  if (missing(data)) {
    stop("RBD needs at least one input, a matrix")
  }
  if (is.null(tol)) {
    tol = 1e-6
    max_cols = min(10, ncol(data))
  } else if (is.null(max_cols)) {
    max_cols = ncol(data)
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
    if (normi < tol) {
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
      #cat(sprintf("Reduced system getting accurate enough - to stop with %d basis functions\n", i))
      bases <- bases[, 1 : i]
      trans_mat <- trans_mat[1 : i, ]
    } else {
      i <- i + 1
    }

  }
  list(bases = bases, trans_mat = trans_mat)
}

compute_B <- function(C, tol = NULL, max_cols = NULL){
  N <- dim(C)[1]
  
  E <- transformation_matrix(N)
  
  # Set defaults for exact decomposition
  # E has rank N-1 or N, so we need max_cols >= rank(E)
  if (is.null(tol)) {
    tol <- 1e-13
  }
  if (is.null(max_cols)) {
    max_cols <- N  # E has at most rank N
  }
  
  E_tilde <- rbd(E, tol, max_cols)
  E_rbd <- E_tilde$trans_mat
  E_bases <- E_tilde$bases
  B <- E_bases %*% (E_rbd %*% C %*% t(E_rbd)) %*% t(E_bases)
  
  B
}

##' Compute q vector directly from C using the explicit pipeline
##'
##' This function implements the sequence:
##' E <- transformation_matrix(N)
##' E_rbd <- rbd(E, tol = tol, max_cols = N)
##' Y <- E_rbd$bases; E_tilde <- E_rbd$trans_mat
##' C_tilde <- E_tilde %*% C %*% t(E_tilde)
##' eig <- eigen(C_tilde)
##' V <- Y %*% eig$vectors; Lambda <- eig$values
##' q <- rowSums((V^2) * Lambda) / sum(Lambda)

compute_design_probs_rbd <- function(C, tol = 1e-13, max_cols = NULL){

  N <- dim(C)[1]

  if (is.null(max_cols)) max_cols <- N

  E <- transformation_matrix(N)
  E_rbd_res <- rbd(E, tol = tol, max_cols = max_cols)
  Y <- E_rbd_res$bases
  E_tilde <- E_rbd_res$trans_mat

  C_tilde <- E_tilde %*% C %*% t(E_tilde)
  eig <- eigen(C_tilde, symmetric = TRUE)

  V <- Y %*% eig$vectors
  Lambda <- eig$values

  denom <- sum(Lambda)
  if (denom == 0) stop("sum of eigenvalues is zero; cannot normalise")

  contributions <- sweep(V^2, 2, Lambda, FUN = "*")
  q <- rowSums(contributions) / denom

  q
}










