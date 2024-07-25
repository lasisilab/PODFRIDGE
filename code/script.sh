#!/bin/bash
#SBATCH --job-name=known_vs_tested_simulation
#SBATCH --output=/home/%u/%u/slurm/%x-%j.log
#SBATCH --time=14-00:00:00
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=4000MB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER}@umich.edu

# Get the unique username
UNIQNAME=$USER

# Activate conda
CONDA_HOME="/home/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh" || exit 1

# Create the environment using the environment.yml file
conda env create -f /path/to/environment.yml

# Activate the environment
conda activate rstats

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/logfiles

# Run the R script
Rscript /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/code/known-vs-tested_simulation_script.R 100000 50000

