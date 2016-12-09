#!/bin/bash
#SBATCH --partition=queue_1
#SBATCH -n 1
#SBATCH -N 1

# Load modules.
module load openbabel
module load hpc/python-2.7.3
module load centos6/rdkit-2013.09.1_with_inchi

# Setup environment variables.
AG_SOFTWARE_DIR=/n/aagfs1/software
MOLGEN_EXE="python $AG_SOFTWARE_DIR/confgen/run_generator.py"

# Body.
# Setup scratch dir.
CWD=$(pwd)
SCRATCH_DIR="/scratch/tmp/confgen.$$/"
mkdir -p $SCRATCH_DIR

cp smiles.smi $SCRATCH_DIR
cd $SCRATCH_DIR
$MOLGEN_EXE \
  --output_limit="output_limit_value" \
  --candidate_limit="candidate_limit_value" \
  --forcefield_id="forcefield_id_value" \
  --output_dir="output_dir_value" \
  --smiles="smiles_value" \
  --inchikey="inchikey_value" \
  --prune_rms_thresh="prune_rms_thresh_value" \
  --energy_range="energy_range_value" \
  --energy_cluster_candidate_radius="energy_cluster_candidate_radius_value" \
  --energy_cluster_auto_accept_radius="energy_cluster_auto_accept_radius_value" \
  --rms_cluster_radius="rms_cluster_radius_value" \
  --fallback_to_align="fallback_to_align_value"

cp -r $SCRATCH_DIR/* $CWD/
rm -r $SCRATCH_DIR
