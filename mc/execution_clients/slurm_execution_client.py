import logging
import os
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
        if run_status != 'RUNNING':
            execution_state = self.get_final_execution_state(
                execution_meta=execution_meta,
                completed_execution_state=execution_state
            )
        return execution_state

    def get_final_execution_state(self, execution_meta=None,
                                  completed_execution_state=None):
        final_execution_state = {**completed_execution_state}
        submission = execution_meta['submission']
        completed_dir = submission['dir']
        final_execution_state['completed_dir'] = completed_dir
        if 'checkpoint_files' in submission:
            try:
                self.validate_checkpoint_files(
                    checkpoint_files=submission['checkpoint_files'],
                    completed_dir=completed_dir
                )
            except Exception as error:
                final_execution_state['run_status'] = 'FAILED'
                final_execution_state['error'] = repr(error)
        return final_execution_state

    def validate_checkpoint_files(self, checkpoint_files=None,
                                  completed_dir=None):
        if 'completed' in checkpoint_files:
            completed_checkpoint_path = os.path.join(
                completed_dir, checkpoint_files['completed'])
            if not self.file_exists(path=completed_checkpoint_path):
                if 'failed' in checkpoint_files:
                    failed_checkpoint_path = os.path.join(
                        completed_dir, checkpoint_files['failed'])
                    try:
                        error = ("Contents of failure checkpoint file:\n"
                                 + self.read_file(path=failed_checkpoint_path))
                    except Exception as read_error:
                        error = (
                            "Error unknown, unable to read checkpoint file"
                            " '{path}': {read_error}."
                        ).format(
                            path=failed_checkpoint_path,
                            read_error=repr(read_error)
                        )
                    raise Exception(error)

    def file_exists(self, path=None):
        try:
            self.process_runner.run_process(cmd=['ls', path], check=True)
            return True
        except:
            return False

    def read_file(self, path=None):
        completed_proc = self.process_runner.run_process(
            cmd=['cat', path], check=True)
        return completed_proc.stdout

    def get_run_status_from_slurm_job_meta(self, slurm_job_meta=None):
        job_state = slurm_job_meta['JobState']
        for run_status, job_states in SLURM_JOB_STATES_TO_RUN_STATES.items():
            if job_state in job_states: return run_status
        return 'UNKNOWN'

    def parse_scontrol_output(self, scontrol_output=None):
        parsed = {}
        for key_value_str in scontrol_output.split():
            key, value = key_value_str.split('=', 1)
            parsed[key] = value
        return parsed
