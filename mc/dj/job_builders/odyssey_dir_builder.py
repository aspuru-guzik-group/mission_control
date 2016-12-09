import collections
import os
import tempfile
import textwrap


class OdysseyDirBuilder(object):
    def build_job(self, job_spec=None, output_dir=None):
        if not output_dir:
            output_dir = tempfile.mkdtemp(prefix='odyssey_dir.')
        dir_spec = self.generate_dir_spec(job_spec=None)
        self.build_dir(dir_spec=dir_spec, output_dir=output_dir)
        return output_dir

    def generate_dir_spec(self, job_spec=None):
        dir_spec = self.generate_initial_dir_spec(job_spec=job_spec)
        if job_spec.get('type') == 'molgen':
            self.update_dir_spec_for_molgen_job(job_spec=job_spec)

    def generate_initial_dir_spec(self, job_spec=None):
        dir_spec = collections.defaultdict(dict)
        dir_spec['sbatch_params'] = self.generate_initial_sbatch_params(
            job_spec=job_spec)
        dir_spec['modules'] = set()

    def generate_initial_sbatch_params(self, job_spec=None):
        initial_sbatch_params = {
            'nodes': 1,
            'ntasks': 1,
            'partition': 'my_partion'
        }
        return initial_sbatch_params

    def update_dir_spec_for_molgen_job(self, job_spec=None, dir_spec=None):
        dir_spec['sbatch_params'].update({'nodes': 1, 'ntasks': 1})
        dir_spec['modules'].update(set([
            'openbabel',
            'hpc/python-2.7.3',
            'centos6/rdkit-2013.09.1_with_inchi',
        ]))
        dir_spec['env_vars'].update({
            'MOLGEN_EXE': 'python $AG_SOFTWARE_DIR/confgen/run_generator.py'
        })
        dir_spec['body'] = self.generate_body_for_molgen_job(job_spec=job_spec)

    ## @TODO: HERE!!!! CLEAN THIS UP!
    def generate_job_body_for_molgen_job(self, job_spec=None):
        job_body = textwrap.dedent(
            """
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
            """).lstrip()
        return job_body


