if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
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
    
    if(N <= 128) {
      tic <- Sys.time()
      B_brute <- bruteforceB(C)
      design_probs_brute <- compute_design_probs(N, B_brute)
      toc <- Sys.time()
      time.brute <- as.numeric(difftime(toc, tic, units = "secs"))
  
      # Compute KL Divergence
      kl_divergence <- as.numeric(philentropy::KL(rbind(t(design_probs_brute), t(design_probs_rbd))))
    } else {
      time.brute <- NA
      kl_divergence <- NA
      
    }

    # Brute Force Method -------------------------------------------------------

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
tols <- c(1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16)

full_iteration_params <- merge(
  expand.grid(
    N = Ns,
    tol = tols
  ),
  matrix_types,
  all = TRUE
) %>% 
  mutate(
    matrix_type = as.character(matrix_type),
    aux = round(aux, 5)
  ) %>% 
  bind_rows(data.frame( ## Also include an this scenario for the laplcian
    "N" = 512,
    "tol" = 1e-6,
    "matrix_type" = "laplacian",
    "aux" = 0.5
  ))


param_index_file <- file.path("results", "parameter_index.csv")

if (!file.exists(param_index_file)) {
  # Case 1: First run - create fresh index
  message("Creating new parameter index file...")
  param_index <- full_iteration_params %>%
    mutate(index = row_number()) %>%
    dplyr::select(index, N, tol, matrix_type, aux)
  
  write.csv(param_index, file = param_index_file, row.names = FALSE)
  
} else {
  # Case 2: Update existing run - append new params only
  message("Reading existing parameter index...")
  existing_index <- read.csv(param_index_file) %>%
    dplyr::select(index, N, tol, matrix_type, aux) %>% 
    mutate(aux = round(aux, 5))
  
  # Identify new combinations by anti-joining the full set with the existing set
  # Note: we exclude the 'index' column for comparison
  new_params <- anti_join(
    full_iteration_params,
    existing_index,
    by = c("N", "tol", "matrix_type", "aux")
  )
  
  if (nrow(new_params) > 0) {
    message("Found ", nrow(new_params), " new parameter combinations. Appending to index.")
    
    # Assign new indices starting from the previous maximum
    start_index <- max(existing_index$index) + 1
    new_params$index <- seq(start_index, start_index + nrow(new_params) - 1)
    
    # Ensure column order matches
    new_params <- new_params %>% dplyr::select(index, N, matrix_type, tol, aux)
    
    # Combine and save
    param_index <- bind_rows(existing_index, new_params)
    write.csv(param_index, file = param_index_file, row.names = FALSE)
  } else {
    message("No new parameter combinations found.")
    param_index <- existing_index
  }
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



y <- existing_index %>% filter(matrix_type == "gnp", aux == 0.3, tol == 1e-6)
x <- full_iteration_params %>% filter(matrix_type == "gnp", aux == 0.3, tol == 1e-6)
full_iteration_params %>% filter(matrix_type == "gnp", aux == 0.3)
full_iteration_params %>% filter(matrix_type == "gnp", aux == full_iteration_params$aux[103])

str(full_iteration_params$aux[103])


full_iteration_params %>% filter(matrix_type == "gnp", aux == 0.3, tol == 1e-6)


y
anti_join(
  full_iteration_params,
  y,
  by = c("N", "tol", "matrix_type", "aux")
)
