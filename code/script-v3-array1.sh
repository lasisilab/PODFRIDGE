#This is a script for an array submission in which multiple identical jobs are run on the cluster at the same time
#The outputs will write to folders with suffixes named after each of the array tasks
#They can then be joined into one big table before moving on to subsequent tasks

#!/bin/bash
#SBATCH --job-name=simulate_genotypes_0923
#SBATCH --output=/home/%u/%u/PODFRIDGE/slurm/%x-%j.log
#SBATCH --time=07:00:00
#SBATCH --partition=standard
#SBATCH --account=tlasisi0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=1000m
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=robheath@umich.edu
#SBATCH --verbose
#SBATCH --array=0-5

########################################################################
# Reminder: If pushing to github automatically, please set your GitHub PAT using the following command before running this script:                                          #
#export GITHUB_PAT='pat'                       #
#export GITHUB_USER='username'                            #
########################################################################

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
    squeue --job=$SLURM_JOB_ID --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R %C %m %N" | tee -a /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/slurm/resource_usage_$SLURM_JOB_ID.log
}

# Run the logging function every 15 minutes in the background
while true; do
    log_resource_usage
    sleep 900  # Pause for 15 minutes
done &

# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

# Run the R script with command line arguments
#R CMD BATCH --no-save --no-restore simulate_genotypes.R 200 400 script_$SLURM_ARRAY_TASK_ID.out 1
Rscript --vanilla simulate_genotypes.R 4000 8000 script_${i}.out 0

#If automatically pushing to github:
# Configure Git to use HTTPS and PAT
#git remote set-url origin https://github.com/lasisilab/PODFRIDGE.git

# Commit and push changes to GitHub
#git add .
#git commit -m "Automated commit of new results $(date)"
#echo $GITHUB_PAT | git push https://$GITHUB_USER:$GITHUB_PAT@github.com/lasisilab/PODFRIDGE.git

# Log memory and CPU usage after job completion
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocCPUs,Elapsed,MaxRSS,MaxVMSize,State > /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/slurm/usage_$SLURM_JOB_ID.log
