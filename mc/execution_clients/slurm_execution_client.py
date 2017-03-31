import logging
import re


SLURM_JOB_STATES_TO_RUN_STATES = {
    'RUNNING': set(['CONFIGURING', 'COMPLETING', 'PENDING', 'RUNNING']),
    'COMPLETED': set(['COMPLETED']),
    'FAILED': set(['BOOT_FAIL', 'CANCELLED', 'FAILED', 'NODE_FAIL', 'PREEMPTED',
                   'TIMEOUT'])
}

class SlurmExecutionClient(object):
    def __init__(self, process_runner=None, logger=None):
        self.process_runner = process_runner
        self.logger = logger or logging

    def start_execution(self, submission=None):
        workdir = submission['dir']
        entrypoint = workdir + '/' + submission['entrypoint']
        cmd = ['sbatch', '--workdir=%s' % workdir, entrypoint]
        try:
            completed_proc = self.process_runner.run_process(cmd=cmd, check=True)
        except self.process_runner.CalledProcessError as called_process_error:
            error_msg = ("Submission error: {error}, stdout: {stdout},"
                         " stderr: {stderr}").format(
                             error=called_process_error,
                             stdout=called_process_error.stdout,
                             stderr=called_process_error.stderr,
                         )
            raise Exception(error_msg)
        slurm_job_id = self.parse_sbatch_output(completed_proc.stdout)
        execution_meta = {
            'job_id': slurm_job_id,
            'submission': submission,
        }
        return execution_meta

    def parse_sbatch_output(self, sbatch_output=None):
        match = re.match(r'Submitted batch job (\d+)', sbatch_output)
        if not match:
            raise Exception("Could not parse slurm job id."
                            " sbatch output was: '%s'" % sbatch_output)
        slurm_job_id = match.group(1)
        return slurm_job_id

    def get_execution_state(self, execution_meta=None):
        job_id = execution_meta['job_id']
        cmd = ['scontrol', 'show', '--details', '--oneliner', 'job', job_id]
        completed_proc = self.process_runner.run_process(cmd=cmd, check=True)
        slurm_job_meta = self.parse_scontrol_output(completed_proc.stdout)
        run_status = self.get_run_status_from_slurm_job_meta(slurm_job_meta)
        execution_state = {
            'run_status': run_status,
            'slurm_job_meta': slurm_job_meta,
        }
        return execution_state

    def parse_scontrol_output(self, scontrol_output=None):
        parsed = {}
        for key_value_str in scontrol_output.split():
            key, value = key_value_str.split('=', 1)
            parsed[key] = value
        return parsed

    def get_run_status_from_slurm_job_meta(self, slurm_job_meta=None):
        job_state = slurm_job_meta['JobState']
        for run_status, job_states in SLURM_JOB_STATES_TO_RUN_STATES.items():
            if job_state in job_states: return run_status
        return 'UNKNOWN'
