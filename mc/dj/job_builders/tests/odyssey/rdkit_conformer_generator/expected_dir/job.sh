#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=my_partition

MOLGEN_EXE="python /n/aagfs01/software/a2g2-utils.repo/a2g2_utils/rdkit_conformer_generator/run_generator.py"

module load openbabel
module load hpc/python-2.7.3
module load centos6/rdkit-2013.09.1_with_inchi

CWD=$(pwd)
SCRATCH_DIR="/scratch/tmp/rdkit_conformer_generator.$$/"
mkdir -p $SCRATCH_DIR

cp smiles.smi $SCRATCH_DIR
cd $SCRATCH_DIR
$MOLGEN_EXE \
  --output_limit="output_limit_value" \
  --candidate_limit="candidate_limit_value" \
  --forcefield_id="forcefield_id_value" \
  --output_dir="output_dir_value" \
  --smiles="smiles_value" \
  --xyz_filename_tpl="inchikey_value_conf_{index}.xyz" \
  --prune_rms_thresh="prune_rms_thresh_value" \
  --energy_range="energy_range_value" \
  --energy_cluster_candidate_radius="energy_cluster_candidate_radius_value" \
  --energy_cluster_auto_accept_radius="energy_cluster_auto_accept_radius_value" \
  --rms_cluster_radius="rms_cluster_radius_value" \
  --fallback_to_align="fallback_to_align_value"

cp -r $SCRATCH_DIR/* $CWD/
rm -r $SCRATCH_DIR

