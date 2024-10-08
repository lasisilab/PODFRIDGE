#!/bin/bash
#SBATCH --output=/home/%u/PODFRIDGE/slurm/%x-%j.log
#SBATCH --time=07:00:00
#SBATCH --partition=standard
#SBATCH --account=tlasisi0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=1g
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=robheath@umich.edu
#SBATCH --verbose
#SBATCH --array=0-50

current_dir=${PWD##*/}
echo "Current dir: $current_dir"

# Get the unique username
UNIQNAME=$(whoami)

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/PODFRIDGE
mkdir -p /home/$UNIQNAME/PODFRIDGE/slurm
mkdir -p /home/$UNIQNAME/PODFRIDGE/logfiles
mkdir -p /home/$UNIQNAME/PODFRIDGE/data
mkdir -p /home/$UNIQNAME/PODFRIDGE/data/sims
mkdir -p /home/$UNIQNAME/PODFRIDGE/output

# Navigate to the project directory
cd/home/$UNIQNAME/PODFRIDGE

# Activate conda environment
CONDA_HOME="/home/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"
conda init
conda activate rstats #where rstats is the environment name

# Log memory and CPU usage every 15 minutes
log_resource_usage() {
    squeue --job=$SLURM_JOB_ID --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %C %m %N" | tee -a /home/$UNIQNAME/PODFRIDGE/slurm/resource_usage_$SLURM_JOB_ID.log
}

# Run the logging function every 15 minutes in the background
while true; do
    log_resource_usage
    sleep 900  # Pause for 15 minutes
done &

# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

# Run the R script with command line arguments

current_dir=${PWD##*/}
echo "Current dir: $current_dir"

Rscript --vanilla PODFRIDGE/code/simulate_genotypes.R 400 800 script_${i}.out $SLURM_JOB_NAME 0
Rscript --vanilla PODFRIDGE/code/lr.R ${SLURM_ARRAY_TASK_ID} $SLURM_JOB_NAME

# Log memory and CPU usage after job completion
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocCPUs,Elapsed,MaxRSS,MaxVMSize,State > /home/$UNIQNAME/PODFRIDGE/slurm/usage_$SLURM_JOB_ID.log
