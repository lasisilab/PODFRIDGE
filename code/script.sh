#!/bin/bash
#SBATCH --job-name=known_vs_tested_simulation
#SBATCH --output=/home/%u/%u/slurm/%x-%j.log
#SBATCH --time=14-00:00:000
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=4000MB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jakevanc@umich.edu

# stimulate conda
CONDA_HOME="/home/jakevanc/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh" || exit 1

# Load conda env created with conda env create -f env.yml
conda activate simulations

# Ensure the directory structure exists
mkdir -p /home/jakevanc/jakevanc/PODFRIDGE/logfiles

# Run the R script
Rscript /home/jakevanc/jakevanc/PODFRIDGE/code/known-vs-tested_simulation_script.R 100000 50000
