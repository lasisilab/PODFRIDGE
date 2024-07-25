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

# Instructions for exporting the SSH passphrase
# ----------------------------------------------
# Ensure you export the SSH passphrase before running this script.
# Example:
# export SSH_PASSPHRASE="your-ssh-key-passphrase"
# ----------------------------------------------

# Get the unique username
UNIQNAME=$(whoami)

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/$UNIQNAME/slurm
mkdir -p /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/logfiles

# Activate conda
CONDA_HOME="/home/$UNIQNAME/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"

# Activate the environment
conda activate rstats

# Run the R script with command line arguments
Rscript /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/code/STR_sims.R 10 50 /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/data/sim_processed_genotypes.csv /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/data/sim-summary_genotypes.csv

# Start the SSH agent and add the SSH key
eval "$(ssh-agent -s)"
ssh-add /home/$UNIQNAME/.ssh/id_rsa <<< "$SSH_PASSPHRASE"

# Navigate to the directory
cd /home/$UNIQNAME/$UNIQNAME/PODFRIDGE

# Add all changes to git
git add .

# Commit the changes with a timestamp message
git commit -m "Auto commit: Updated results on $(date)"

# Push the changes to the repository
git push

# Kill the SSH agent
ssh-agent -k
