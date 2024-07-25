# Setting up account first time
# In the future start with exporting and go to ssh
export UNIQNAME=tlasisi

ssh $UNIQNAME@greatlakes.arc-ts.umich.edu

cd /nfs/turbo/lsa-tlasisi1
mkdir $UNIQNAME
ln -s /nfs/turbo/lsa-tlasisi1/$UNIQNAME /home/$UNIQNAME
cd $UNIQNAME

# Setting up simulation on Great Lakes
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

#check install isn't corrputed
shasum -a 256 Miniconda3-latest-Linux-x86_64.sh


# run the installer
bash Miniconda3-latest-Linux-x86_64.sh


#install miniconda here
/home/$UNIQNAME/$UNIQNAME/miniconda3


# Install mamba in the base conda environment with verbose output
conda install -n base -c conda-forge mamba -y -v

# Create a new environment named 'rstats' with verbose output
mamba create -n rstats -c conda-forge r -y -v


# Ensure the R and Conda-Forge channels are added
conda config --add channels r
conda config --add channels conda-forge

# Create a new environment named 'rstats'
mamba create -n rstats -c r -c conda-forge r --verbose

# Activate the environment
conda activate rstats

# clone PODFRIDGE
cd /home/$UNIQNAME/$UNIQNAME && git clone https://github.com/lasisilab/PODFRIDGE.git

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
done < code/STR-sims-env.txt

echo "All packages processed."
