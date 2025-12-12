# A Reduced Basis Decomposition Approach to Efficient Data Collection in Pairwise Comparisons Studies

This repository contains R code for computing reduced-basis decomposition (RBD) design probabilities and the corresponding brute-force computations for comparative judgement studies. The implementation provides:

- A greedy reduced-basis decomposition (`rbd`) and helpers in `R/rbd.R`. The main RBD algorithm is based on the [MATLAB code written by Yanlai Chen (2015)](https://www.mathworks.com/matlabcentral/fileexchange/50125-reduced-basis-decomposition). 
- Brute-force construction of the pairwise matrix and spectral design probabilities in `R/brute_force_functions.R`.
- Example scripts under `examples/` demonstrating usage and numeric checks.

**Goal:** provide an efficient pipeline to compute design probabilities for a comparative judgement study from a covariance matrix `C` using a reduced-basis spectral decomposition.

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
- `bruteforceB(C)` — build the full brute-force `B` matrix (size choose(N,2) x choose(N,2)).
- `compute_design_probs_rbd(C, tol = 1e-13, max_cols = NULL)` — compute design probabilities `q` directly from covariance `C` using the RBD pipeline.
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


---
## Acknowlegdements
The development of this software was supported by a UKRI Future Leaders Fellowship [MR/X034992/1].

## Licence for the Original Matlab Code
The original MATLAB code for the RBD algorithm is avaliable at [https://www.mathworks.com/matlabcentral/fileexchange/50125-reduced-basis-decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50125-reduced-basis-decomposition). It is subject to follow conditions:

Copyright (c) 2015, Yanlai Chen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

