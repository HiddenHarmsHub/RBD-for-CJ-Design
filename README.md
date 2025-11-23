# RBD in R



## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.bham.ac.uk/seymourg-bsbt/rbd-in-r.git
git branch -M main
# RBD in R

This repository contains R code for computing reduced-basis design probabilities for networks (RBD) and the corresponding brute-force computations for comparison. The implementation provides:

- A greedy reduced-basis decomposition (`rbd`) and helpers in `R/rbd.R`.
- Brute-force construction of the pairwise matrix and spectral design probabilities in `R/brute_force_functions.R`.
- Example scripts under `examples/` demonstrating usage and numeric checks.

**Goal:** provide an efficient pipeline to compute edge-level design probabilities from a node covariance matrix `C` using a reduced-basis spectral decomposition.

**Edge ordering:** edge rows follow `combn(N, 2)` ordering (pairs with `i < j`).

---

## Requirements

- R (>= 4.0 recommended)
- CRAN packages: `expm`, `igraph`, `philentropy` (for KL divergence in examples)

Install the R dependencies in an R session:

```r
install.packages(c("expm", "igraph", "philentropy"))
```

---

## Key functions

- `rbd(data, tol = NULL, max_cols = NULL)` — greedy reduced-basis decomposition used internally.
- `transformation_matrix(N)` — builds the edge incidence-like transformation matrix `E` (rows correspond to edges ordered by `combn(N,2)`).
- `compute_B(C, tol = NULL, max_cols = NULL)` — compute reduced B via RBD pipeline.
- `compute_design_probs_rbd(C, tol = 1e-13, max_cols = NULL)` — compute design probabilities `q` directly from covariance `C` using the RBD pipeline.
- `bruteforceB(C)` — build the full brute-force `B` matrix (size choose(N,2) x choose(N,2)).
- `compute_design_probs(N, B)` — compute brute-force design probabilities from `B` (keeps legacy name used in examples).

Files:

- `R/rbd.R` — reduced-basis routines and `compute_design_probs_rbd`.
- `R/brute_force_functions.R` — brute-force `B` construction and spectral probabilities.


---

## Quick examples
Within an R session you can run the core steps manually:

```r
source("R/rbd.R")
source("R/brute_force_functions.R")

# build a test covariance
C <- diag(5)  # replace with your SPD covariance (N x N)

# RBD design probabilities
q_rbd <- compute_design_probs_rbd(C)

# Brute-force design probabilities
B <- bruteforceB(C)
q_brute <- compute_design_probs(nrow(C), B)

# Compare (KL requires philentropy)
philentropy::KL(rbind(t(q_brute), t(q_rbd)))
```

The returned `q` vectors have length `choose(N,2)` and correspond to edges ordered by `combn(N,2)`.

---

## Development notes

- Performance: `bruteforceB` is O(m^2) where m = choose(N,2). For large N use the RBD pipeline (`compute_design_probs_rbd`) which avoids forming the full `B`.

---
Contact: create an issue in this repo or edit the files directly — I'm happy to help further.
