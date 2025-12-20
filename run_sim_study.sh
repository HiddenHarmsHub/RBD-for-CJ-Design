#!/bin/bash

#SBATCH --job-name=R_sim_study   # Job name
#SBATCH --output=slurm_output/R_sim_study_%j.out       # Standard output and error log
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of MPI tasks (and Julia workers)
#SBATCH --cpus-per-task=64              # Number of CPU cores per task
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r.g.seymour@bham.ac.uk



set -e

module purge; module load bluebear

module load bear-apps/2024a
module load R/4.5.0-gfbf-2024a

Rscript examples/size_sim_study.R
