library(igraph)
source("R/rbd.R")
source("R/brute_force_functions.R")

N <- 2^c(3:7)
n.iter <- 100
outputs <- vector("list", length(N))
for (k in seq_along(outputs)) {
  outputs[[k]] <- matrix(NA, nrow = n.iter, ncol = 3)
}

for(i in 1:length(N)){
  for(j in 1:n.iter){
    # Generate Network and Matrices -------------------------------------------

    g_np <- sample_gnp(N[i], 0.3, directed = FALSE, loops = FALSE)
    #plot(g_np, vertex.size=5, vertex.label=NA)
    
    signal.variance <- 1
    A <- as_adjacency_matrix(g_np, sparse = FALSE)
    C <- expm::expm(A)
    C <- signal.variance*diag(diag(C)^-0.5)%*%C%*%diag(diag(C)^-0.5) #C = D(C)^-0.5 * C * D(C)^-0.5, ie. normalise C so it has one on the diagonal
    
    
    # RBD Method --------------------------------------------------------------
    tic <- Sys.time()
    design_probs_rbd <- compute_design_probs_rbd(C)
    toc <- Sys.time()

    time.rbd <- as.numeric(difftime(toc, tic, units = "secs"))
    
    
    # Brute method --------------------------------------------------------------
    tic <- Sys.time()
    B_brute <- bruteforceB(C)
    design_probs_brute <- compute_design_probs(N[i], B_brute)
    toc <- Sys.time()

    time.brute <- as.numeric(difftime(toc, tic, units = "secs"))
    kl_divergence <- as.numeric(philentropy::KL(rbind(t(design_probs_brute), t(design_probs_rbd))))
  
    outputs[[i]][j, ] <- c(time.rbd, time.brute, kl_divergence)
    print(paste0("Completed iteration ", j, " of size ", N[i]))
  }
}
