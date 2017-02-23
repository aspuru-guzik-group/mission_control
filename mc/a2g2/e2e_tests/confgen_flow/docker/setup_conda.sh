CONDA_INSTALLER="/tmp/conda_installer.sh"
CONDA_ROOT="${CONDA_ROOT:=$HOME/conda}"
curl $MINICONDA_URL > $CONDA_INSTALLER && chmod +x $CONDA_INSTALLER
bash $CONDA_INSTALLER -b -p $CONDA_ROOT
rm $CONDA_INSTALLER
source $CONDA_ROOT/bin/activate && conda install conda-build
chmod -R a+rwx $CONDA_ROOT
echo "PATH=$CONDA_ROOT/bin:\$PATH" >> ~/.bashrc
