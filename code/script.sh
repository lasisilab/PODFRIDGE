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

# Check and install packages from STR-sims-env.txt with verbose output
while read -r package; do
    echo "Processing package: r-$package"
    if conda list -n rstats | grep -q "^r-$package"; then
        echo "Package r-$package is already installed."
    else
        echo "Installing package r-$package..."
        mamba install -n rstats -c r -c conda-forge "r-$package" --verbose || {
            echo "Failed to install package r-$package"
            exit 1
        }
    fi
done < /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/code/STR-sims-env.txt

# Run the R script with command line arguments
Rscript /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/code/STR_sims.R 10 50 /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/data/sim_processed_genotypes.csv /home/$UNIQNAME/$UNIQNAME/PODFRIDGE/data/sim-summary_genotypes.csv
