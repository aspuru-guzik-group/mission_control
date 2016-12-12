import copy
import textwrap


def alter_dir_spec(dir_spec=None, job_spec=None):
    altered_dir_spec = copy.deepcopy(dir_spec)
    altered_dir_spec['sbatch_params'].update({'nodes': 1, 'ntasks': 1})
    altered_dir_spec['modules'].update([
        ('openbabel', True),
        ('hpc/python-2.7.3', True),
        ('centos6/rdkit-2013.09.1_with_inchi', True)
    ])
    altered_dir_spec['env_vars'].update({
        'MOLGEN_EXE': ('"python /n/aagfs01/software/a2g2-utils.repo/a2g2_utils/'
                       'rdkit_conformer_generator/run_generator.py"')
    })
    altered_dir_spec['script_body'] = generate_script_body(job_spec=job_spec)
    return altered_dir_spec

def generate_script_body(job_spec=None):
    script_body = textwrap.dedent(
        """
        CWD=$(pwd)
        SCRATCH_DIR="/scratch/tmp/rdkit_conformer_generator.$$/"
        mkdir -p $SCRATCH_DIR

        cp smiles.smi $SCRATCH_DIR
        cd $SCRATCH_DIR
        $MOLGEN_EXE \\
          --output_limit="output_limit_value" \\
          --candidate_limit="candidate_limit_value" \\
          --forcefield_id="forcefield_id_value" \\
          --output_dir="output_dir_value" \\
          --smiles="smiles_value" \\
          --xyz_filename_tpl="{inchikey}_conf_{{index}}.xyz" \\
          --prune_rms_thresh="prune_rms_thresh_value" \\
          --energy_range="energy_range_value" \\
          --energy_cluster_candidate_radius="energy_cluster_candidate_radius_value" \\
          --energy_cluster_auto_accept_radius="energy_cluster_auto_accept_radius_value" \\
          --rms_cluster_radius="rms_cluster_radius_value" \\
          --fallback_to_align="fallback_to_align_value"

        cp -r $SCRATCH_DIR/* $CWD/
        rm -r $SCRATCH_DIR
        """).format(inchikey=job_spec['params']['inchikey']).lstrip()
    return script_body


