library(pacman)
pacman::p_load(tidyverse, igraph, parallel, parallelly, MCMCpack, philentropy)
source("R/rbd.R")
source("R/brute_force_functions.R")

# ============================================================================
# Matrix Generation Functions
# ============================================================================

## Generate network-based covariance matrix using sample_gnp
generate_matrix_gnp <- function(n, prob, signal.variance = 1) {
  g_np <- sample_gnp(n, prob, directed = FALSE, loops = FALSE)
  A <- as_adjacency_matrix(g_np, sparse = FALSE)
  C <- expm::expm(A)
  C <- signal.variance * diag(diag(C)^-0.5) %*% C %*% diag(diag(C)^-0.5)
  C
}

## Generate regularized graph Laplacian
generate_matrix_laplacian <- function(n, prob, lambda = 0.1, signal.variance = 1) {
  g_np <- sample_gnp(n, prob, directed = FALSE, loops = FALSE)
  A <- as_adjacency_matrix(g_np, sparse = FALSE)
  D <- diag(colSums(A))
  L <- D - A  # Graph Laplacian
  # Regularize: L_reg = (1-lambda)*I + lambda*L_normalized
  L_norm <- L / (max(abs(L)) + 1e-10)
  C <- (1 - lambda) * diag(n) + lambda * L_norm
  C <- signal.variance * diag(diag(C)^-0.5) %*% C %*% diag(diag(C)^-0.5)
  C
}

## Generate Toeplitz matrix
generate_matrix_toeplitz <- function(n, signal.variance = 1, rho = NULL) {
  # Create Toeplitz matrix with correlation structure
  if(is.null(rho)) rho <- runif(1, min = 0.3, max = 0.9)
  acf_values <- rho^(0:(n-1))
  C <- toeplitz(acf_values)
  C <- signal.variance * diag(diag(C)^-0.5) %*% C %*% diag(diag(C)^-0.5)
  C
}

## Generate inverse Wishart matrix
generate_matrix_inverse_wishart <- function(n, df = NULL, signal.variance = 1) {
  if (is.null(df)) df <- n + 1
  # Generate from inverse Wishart via Wishart
  
  S <- diag(n)  # Scale matrix
  W <- riwish(df, S)
  C <- signal.variance * diag(diag(W)^-0.5) %*% W %*% diag(diag(W)^-0.5)
  C
}

# ================================================
# Worker Function for Parallel Execution
# ================================================
run_single_iteration <- function(
  iteration,
  N,
  matrix_type = "gnp",
  rbd_tol = 1e-13,
  aux = NULL,
  output_dir = file.path("results", "simulation_study")
) {
  # Generate output file path
  file_name <- file.path(
    output_dir, 
    paste0(
      "N", N, 
      "_", matrix_type, 
      "_tol", rbd_tol, 
      "_aux", aux*10, 
      "_iteration", iteration,
      ".csv"
    )
  )

  # Check if file exists
  if (file.exists(file_name)) {
    message("File exists: ", file_name, ". Skipping computation.")
    return(NULL)
  }

  # Generate Covariance Matrix -----------------------------------------------
  C <- switch(as.character(matrix_type),
    "gnp" = generate_matrix_gnp(N, prob = aux),
    "laplacian" = generate_matrix_laplacian(N, prob = aux),
    "toeplitz" = generate_matrix_toeplitz(N),
    "inverse_wishart" = generate_matrix_inverse_wishart(N),
    stop("Unknown matrix_type: ", matrix_type)
  )

  # RBD Method ---------------------------------------------------------------
  tic <- Sys.time()
  design_probs_rbd <- compute_design_probs_rbd(C, tol = rbd_tol)
  toc <- Sys.time()
  time.rbd <- as.numeric(difftime(toc, tic, units = "secs"))

  # Brute Force Method -------------------------------------------------------
  tic <- Sys.time()
  B_brute <- bruteforceB(C)
  design_probs_brute <- compute_design_probs(N, B_brute)
  toc <- Sys.time()
  time.brute <- as.numeric(difftime(toc, tic, units = "secs"))

  # Compute KL Divergence
  kl_divergence <- as.numeric(philentropy::KL(rbind(t(design_probs_brute), t(design_probs_rbd))))

  # Save results to file
  results <- data.frame(time.rbd = time.rbd, time.brute = time.brute, kl_divergence = kl_divergence)
  write.csv(results, file = file_name, row.names = FALSE)

  return(results)
}

# ================================================
# Parallel Execution of Iterations
# ================================================
output_dir <- file.path("results", "simulation_study")
dir.create(output_dir, showWarnings = FALSE)

# Create list of all iterations to run
ER_matrices <- expand.grid(
  matrix_type = "gnp",
  aux = seq(from = 0.1, to = 0.9, by = 0.1)
)

RL_matrices <- expand.grid(
  matrix_type = "laplacian",
  aux = seq(from = 0.1, to = 0.9, by = 0.1)
)

IW_matrices <- expand.grid(
  matrix_type = "inverse_wishart"
)

T_matrices <- expand.grid(
  matrix_type = "toeplitz"
)

matrix_types <- bind_rows(
  ER_matrices,
  RL_matrices,
  IW_matrices %>% mutate(aux = NA),
  T_matrices %>% mutate(aux = NA)
)

n.iter <- 2
N <- c(4, 8)

# Create iteration grid by expanding iterations and sizes, then merging with matrix type settings
iteration_params <- merge(
  expand.grid(
    iteration = 1:n.iter,
    size_idx = 1:length(N),
    tols = c(1e-10, 1e-12, 1e-14, 1e-16)
  ),
  matrix_types,
  all = TRUE
)

# Filter iteration_params to exclude existing files
iteration_params <- iteration_params %>% filter(
  !file.exists(paste0(
    output_dir, "N", N[size_idx], "_", matrix_type, "_tol", tols, "_aux", aux, ".csv"
  ))
)

# Run iterations in parallel
parallel::mclapply(1:nrow(iteration_params), function(i) {
  params <- iteration_params[i, ]
  run_single_iteration(
    iteration = params$iteration,
    N = N[params$size_idx],
    matrix_type = params$matrix_type,
    rbd_tol = params$tols,
    aux = params$aux,
    output_dir = output_dir
  )
}, mc.cores = availableCores() - 1)