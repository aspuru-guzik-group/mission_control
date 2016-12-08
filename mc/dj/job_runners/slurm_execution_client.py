import os
import re


SLURM_JOB_STATES_TO_RUN_STATES = {
    'running': ['CONFIGURING', 'COMPLETING', 'PENDING', 'RUNNING']
}

class SlurmExecutionClient(object):
    def __init__(self, process_runner=None):
        self.process_runner = process_runner

    def start_execution(self, job=None):
        cmd = ['sbatch',
               '--workdir=%s' % job['dir']['dir'],
               os.path.join(job['dir']['dir'], job['dir']['entrypoint'])]
        completed_proc = self.process_runner.run_process(cmd=cmd, check=True)
        slurm_job_id = self.parse_sbatch_output(completed_proc.stdout)
        return {'job_id': slurm_job_id}

    def parse_sbatch_output(self, sbatch_output=None):
        match = re.match(r'Submitted batch job (\d+)', sbatch_output)
        if not match:
            raise Exception("Could not parse slurm job id."
                            " sbatch output was: '%s'" % sbatch_output)
        slurm_job_id = match.group(1)
        return slurm_job_id

    def get_execution_state(self, job=None):
        job_id = job['execution']['job_id']
        cmd = ['scontrol', 'show', '--details', '--oneliner', 'job', job_id]
        completed_proc = self.process_runner.run_process(cmd=cmd, check=True)
        slurm_job_meta = self.parse_scontrol_output(completed_proc.stdout)
        run_state = self.get_execution_state_for_slurm_job_meta(slurm_job_meta)
        execution_state = {
            'executing': (run_state == 'running'),
            'slurm_job_meta': slurm_job_meta,
        }
        return execution_state

    def get_execution_state_for_slurm_job_meta(self, slurm_job_meta=None):
        slurm_job_state = slurm_job_meta['JobState']
        if slurm_job_state in SLURM_JOB_STATES_TO_RUN_STATES['running']:
            run_state = 'running'
        else:
            run_state = 'completed'
        return run_state

    def parse_scontrol_output(self, scontrol_output=None):
        parsed = {}
        for key_value_str in scontrol_output.split():
            key, value = key_value_str.split('=', 1)
            parsed[key] = value
        return parsed
