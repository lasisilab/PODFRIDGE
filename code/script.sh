#!/bin/bash
#SBATCH --job-name=STR_sims
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
UNIQNAME=$(whoami)

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/$UNIQNAME/slurm
mkdir -p /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/logfiles

# Activate conda
CONDA_HOME="/home/$UNIQNAME/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"

# Activate the environment
conda activate rstats

# Start the SSH agent and add the SSH key
eval "$(ssh-agent -s)"
echo "$SSH_PASSPHRASE" | ssh-add /home/$UNIQNAME/.ssh/id_rsa

# Run the R script with command line arguments
Rscript /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/code/STR_sims.R 10 50 /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/data/sim_processed_genotypes.csv /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/data/sim-summary_genotypes.csv

# Commit and push the changes to GitHub
cd /home/$UNIQNAME/$UNIQNAME/PODFRIDGE
git add .
git commit -m "Auto commit: Updated results on $(date)"
GIT_SSH_COMMAND="ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" git push origin main
