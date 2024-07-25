#!/bin/bash
#SBATCH --job-name=STR-sims
#SBATCH --output=/home/%u/%u/slurm/%x-%j.log
#SBATCH --time=14-00:00:00
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=4000MB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$(whoami)@umich.edu

########################################################################
# Reminder: Please set your SSH passphrase using the following command #
# before running this script:                                          #
# export SSH_PASSPHRASE='your_passphrase'                              #
########################################################################

# Get the unique username
UNIQNAME=$(whoami)

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/$UNIQNAME/slurm
mkdir -p /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/logfiles

# Change to the project directory
cd /home/$UNIQNAME/$UNIQNAME/PODFRIDGE

# Activate conda
CONDA_HOME="/home/$UNIQNAME/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"

# Activate the environment
conda activate rstats

# Run the R script with command line arguments
Rscript code/STR_sims.R 2 5 data/sim_processed_genotypes.csv data/sim-summary_genotypes.csv

# Set up SSH agent and add SSH key
eval "$(ssh-agent -s)"
echo "$SSH_PASSPHRASE" | ssh-add -

# Commit and push changes to GitHub
git add .
git commit -m "Automated commit of new results $(date)"
git push origin main
