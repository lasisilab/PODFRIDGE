#!/bin/bash
#SBATCH --job-name=STR_sims
#SBATCH --time=01:00:00
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=4000MB
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

# Define the timestamp and create output directory for this job
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
OUTPUT_DIR="/home/$UNIQNAME/$UNIQNAME/PODFRIDGE/output/simulation_${TIMESTAMP}"

# Ensure the directory structure exists
mkdir -p $OUTPUT_DIR/slurm
mkdir -p $OUTPUT_DIR/logfiles

# Update SBATCH output to write to the newly created directory
#SBATCH --output=$OUTPUT_DIR/slurm/%x-%j.log

# Change to the project directory
cd /home/$UNIQNAME/$UNIQNAME/PODFRIDGE

# Activate conda environment
CONDA_HOME="/home/$UNIQNAME/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"
conda activate rstats

# Log memory and CPU usage every 15 minutes
log_resource_usage() {
    squeue --job=$SLURM_JOB_ID --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %C %m %N" | tee -a $OUTPUT_DIR/logfiles/resource_usage_$SLURM_JOB_ID.log
}

# Run the logging function every 15 minutes in the background
while true; do
    log_resource_usage
    sleep 900  # Pause for 15 minutes
done &

# Run the R script with command line arguments for simulation counts and output directory
Rscript code/STR_sims.R 5 10 $OUTPUT_DIR

# Configure Git to use HTTPS and PAT
git remote set-url origin https://github.com/lasisilab/PODFRIDGE.git
git add .
git commit -m "Automated commit of new results $(date)"
echo $GITHUB_PAT | git push https://$GITHUB_USER:$GITHUB_PAT@github.com/lasisilab/PODFRIDGE.git

# Log final memory and CPU usage after job completion
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocCPUs,Elapsed,MaxRSS,MaxVMSize,State > $OUTPUT_DIR/logfiles/usage_$SLURM_JOB_ID.log
