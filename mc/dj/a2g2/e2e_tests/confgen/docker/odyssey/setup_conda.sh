#!/bin/bash

echo "Installing miniconda"
CONDA_INSTALLER_PATH="$HOME/conda_installer.sh"
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > $CONDA_INSTALLER_PATH
chmod +x $CONDA_INSTALLER_PATH 
bash $CONDA_INSTALLER_PATH -b -p $HOME/miniconda3
rm $CONDA_INSTALLER_PATH
echo "PATH=$HOME/miniconda3/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc
conda install conda-build
