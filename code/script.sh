#!/bin/bash
#SBATCH --job-name=STR_sims
#SBATCH --output=/home/%u/%u/slurm/%x-%j.log
#SBATCH --time=1:00:00
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$(whoami)@umich.edu

########################################################################
# Reminder: Please set your GitHub PAT using the following command     #
# before running this script:                                          #
# export GITHUB_PAT='your_personal_access_token'                       #
# export GITHUB_USER='your_github_username'                            #
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

# Function to log resource usage
log_resource_usage() {
    squeue --job=$SLURM_JOB_ID --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %C %m %N" | tee -a logfiles/resource_usage.log
}

# Log resource usage periodically
while true; do
    log_resource_usage
    sleep 3600  # Pause for 1 hour
done &

# Run the R script with command line arguments
Rscript code/STR_sims.R 50 100

# Kill the logging background process
kill %1

# Final log of resource usage
log_resource_usage

# Configure Git to use HTTPS and PAT
git remote set-url origin https://github.com/lasisilab/PODFRIDGE.git

# Commit and push changes to GitHub
git add .
git commit -m "Automated commit of new results $(date)"
echo $GITHUB_PAT | git push https://$GITHUB_USER:$GITHUB_PAT@github.com/lasisilab/PODFRIDGE.git
