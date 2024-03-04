#!/bin/bash
#SBATCH --job-name=known_vs_tested_simulation
#SBATCH --output=/home/tlasisi/PODFRIDGE_upload/logfiles/known_vs_tested_simulation_%j.out
#SBATCH --error=/home/tlasisi/PODFRIDGE_upload/logfiles/known_vs_tested_simulation_%j.err
#SBATCH --time=16:00:00
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4000MB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tlasisi@umich.edu

# Load necessary modules
module load gcc/9.3.0 openblas/0.3.13 r/4.0.5

# Ensure the directory structure exists
mkdir -p /home/tlasisi/PODFRIDGE_upload/logfiles

# Run the R script
Rscript /home/tlasisi/PODFRIDGE_upload/code/known-vs-tested_simulation_script.R 100000 5000
