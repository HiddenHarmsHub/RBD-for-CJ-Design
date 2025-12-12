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
run_scenario <- function(
  index, # Use index instead of iteration-specific parameters
  N,
  matrix_type,
  tol,
  aux = NA,
  n_iter, # Number of iterations to run
  output_dir = file.path("results", "simulation_study")
) {
  # Generate output file path
  file_name <- file.path(output_dir, paste0("results_", index, ".csv"))

  # Check if file exists
  if (file.exists(file_name)) {
    message("File exists: ", file_name, ". Skipping computation.")
    return(NULL)
  }

  output_list <- list()

  for(iter in 1:n_iter) {
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
    design_probs_rbd <- compute_design_probs_rbd(C, tol = tol)
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
    output_list[[iter]] <- data.frame(
      time.rbd = time.rbd,
      time.brute = time.brute,
      kl_divergence = kl_divergence
    )
  }

  write.csv(bind_rows(output_list), file = file_name, row.names = FALSE)
}

output_dir <- file.path("results", "simulation_study")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create list of all iterations to run
ER_matrices <- expand.grid(
  matrix_type = "gnp",
  aux = seq(from = 0.1, to = 0.8, by = 0.1)
)

RL_matrices <- expand.grid(
  matrix_type = "laplacian",
  aux = seq(from = 0.1, to = 0.8, by = 0.1)
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

n.iter <- 100
Ns <- 2^(3:7)
tols <- c(1e-10, 1e-12, 1e-14, 1e-16)

iteration_params <- merge(
  expand.grid(
    N = Ns,
    tol = tols
  ),
  matrix_types,
  all = TRUE
)

# Create parameter index file
param_index_file <- file.path("results", "parameter_index.csv")
if (!file.exists(param_index_file)) {
  param_index <- iteration_params %>%
    mutate(index = row_number()) %>%
    dplyr::select(index, N, matrix_type, tol, aux)
  write.csv(param_index, file = param_index_file, row.names = FALSE)
} else {
  param_index <- read.csv(param_index_file)
}

# Run iterations in parallel
parallel::mclapply(
  1:nrow(param_index),
  function(i) {
    run_scenario(
      index = param_index$index[i],
      N = param_index$N[i],
      matrix_type = param_index$matrix_type[i],
      tol = param_index$tol[i],
      aux = param_index$aux[i],
      n_iter = n.iter,
      output_dir = output_dir
    )
  },
  mc.cores = availableCores()
)

# Combine output files
data_files <- list.files(
  path = "results/simulation_study",
  pattern = "*.csv",
  full.names = TRUE
)

combined_data <- lapply(data_files, function(file) {
    # Extract index from the file name
    index <- as.numeric(str_match(basename(file), "results_(\\d+)\\.csv")[2])

    # Retrieve parameters from param_index
    params <- param_index %>% filter(index == !!index)

    # Read the data from the file
    data <- read.csv(file)

    # Add parameter details to the data
    data$N <- params$N
    data$matrix_type <- params$matrix_type
    data$rbd_tol <- params$tol
    data$aux <- params$aux
    data$index <- index

    data
  }) %>%
  bind_rows()

# Save combined data to a new file
write.csv(
  combined_data,
  file = "results/combined_simulation_study_results.csv",
  row.names = FALSE
)