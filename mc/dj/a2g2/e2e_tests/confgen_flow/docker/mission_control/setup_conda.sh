CONDA_INSTALLER_PATH="$HOME/conda_installer.sh"
CONDA_ROOT="${CONDA_INSTALL_TARGET:=$HOME/conda}"
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > $CONDA_INSTALLER_PATH
chmod +x $CONDA_INSTALLER_PATH 
bash $CONDA_INSTALLER_PATH -b -p $CONDA_ROOT
rm $CONDA_INSTALLER_PATH
echo "PATH=$CONDA_ROOT/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc
conda install conda-build
