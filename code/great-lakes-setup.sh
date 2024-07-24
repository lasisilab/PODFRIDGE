# Setting up simulation on Great Lakes
#
#
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# run the installer
#
bash Miniconda3-latest-Linux-x86_64.sh


# Ensure the R and Conda-Forge channels are added
conda config --add channels r
conda config --add channels conda-forge

# Create a new environment named 'rstats'
conda create -n rstats -c r -c conda-forge r -y

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rstats

# Read the STR-sims-env.txt file and install each package
while read -r package; do
    conda install -n rstats -c r -c conda-forge "r-$package" -y
done < data/STR-sims-env.txt
