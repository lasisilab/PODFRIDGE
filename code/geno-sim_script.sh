#!/bin/bash
#SBATCH --job-name=STR_sims
#SBATCH --output=/home/%u/%u/PODFRIDGE/slurm/%x-%j.log
#SBATCH --time=07:00:00
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=4000MB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$(whoami)@umich.edu

# Get the unique username
UNIQNAME=$(whoami)

# Ensure the directory structure exists
mkdir -p /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/slurm
mkdir -p /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/logfiles

# Change to the project directory
cd /home/$UNIQNAME/$UNIQNAME/PODFRIDGE

# Activate conda
CONDA_HOME="/home/$UNIQNAME/$UNIQNAME/miniconda3"
source "$CONDA_HOME/etc/profile.d/conda.sh"

# Activate the environment
conda activate rstats

# Run the genotype simulation
Rscript code/simulate_genotypes.R 100 500 $SLURM_JOB_ID

# Configure Git to use HTTPS and PAT
git remote set-url origin https://github.com/lasisilab/PODFRIDGE.git

# Log memory and CPU usage after job completion
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocCPUs,Elapsed,MaxRSS,MaxVMSize,State > /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/slurm/usage_$SLURM_JOB_ID.log

# Commit and push changes to GitHub
git add .
git commit -m "Automated commit of new results $(date)"
echo $GITHUB_PAT | git push https://$GITHUB_USER:$GITHUB_PAT@github.com/lasisilab/PODFRIDGE.git
