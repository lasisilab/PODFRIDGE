#!/bin/bash
#SBATCH --output=/home/%u/PODFRIDGE/slurm/%x-%j.log
#SBATCH --time=07:00:00
#SBATCH --partition=standard
#SBATCH --account=tlasisi0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=3g
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=robheath@umich.edu
#SBATCH --verbose

current_dir=${PWD##*/}
echo "Current dir: $current_dir"

# Get the unique username
UNIQNAME=$(whoami)

# Navigate to the project directory
cd/home/$UNIQNAME/PODFRIDGE


# Activate conda
CONDA_HOME="/home/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"

# Activate the environment
conda init
conda activate rstats

# Log memory and CPU usage every 15 minutes
log_resource_usage() {
    squeue --job=$SLURM_JOB_ID --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %C %m %N" | tee -a /home/$UNIQNAME/PODFRIDGE/slurm/resource_usage_$SLURM_JOB_ID.log
}

# Run the logging function every 15 minutes in the background (Lets try without this temporarily to check the sleep command isn't interfering)
while true; do
    log_resource_usage
    sleep 900  # Pause for 15 minutes
done &


current_dir=${PWD##*/}
echo "Current dir: $current_dir"

# Run the R script(s) with command line arguments
Rscript --vanilla PODFRIDGE/code/lr_cutoff.R $SLURM_JOB_NAME
Rscript --vanilla PODFRIDGE/code/lr_plots.R $SLURM_JOB_NAME

# Log memory and CPU usage after job completion
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocCPUs,Elapsed,MaxRSS,MaxVMSize,State > /home/$UNIQNAME/PODFRIDGE/slurm/usage_$SLURM_JOB_ID.log
