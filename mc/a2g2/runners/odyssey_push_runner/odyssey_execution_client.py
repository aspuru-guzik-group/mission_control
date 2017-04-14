import os

from mc.execution_clients.remote_slurm_execution_client import (
    RemoteSlurmExecutionClient)


class OdysseyExecutionClient(object):
    def __init__(self, ssh_client=None):
        self.ssh_client = ssh_client
        self.remote_slurm_client = RemoteSlurmExecutionClient(
            ssh_client=self.ssh_client)

    def start_execution(self, *args, **kwargs):
        return self.remote_slurm_client.start_execution(*args, **kwargs)

    def get_execution_state(self, execution_meta=None):
        slurm_execution_state = self.remote_slurm_client.get_execution_state(
            execution_meta=execution_meta)
        if slurm_execution_state['run_status'] == 'RUNNING':
            execution_state = slurm_execution_state
        else:
            execution_state = self.decorate_completed_slurm_execution_state(
                execution_meta=execution_meta,
                slurm_execution_state=slurm_execution_state
            )
        return execution_state

    def decorate_completed_slurm_execution_state(self, execution_meta=None,
                                                 slurm_execution_state=None):
        decorated_execution_state = {**slurm_execution_state}
        submission = execution_meta['submission']
        completed_dir = submission['dir']
        decorated_execution_state['artifact'] = self.generate_artifact_spec(
            completed_dir=completed_dir)
        if 'checkpoint_files' in submission:
            try:
                self.validate_checkpoint_files(
                    checkpoint_files=submission['checkpoint_files'],
                    completed_dir=completed_dir)
            except Exception as exception:
                decorated_execution_state['run_status'] = 'FAILED'
                decorated_execution_state['error'] = \
                        self.stringify_exception(exception)
        decorated_execution_state['stdout'] = self.get_stdout(
            submission=submission, completed_dir=completed_dir)
        return decorated_execution_state

    def generate_artifact_spec(self, completed_dir=None):
        artifact_spec = {
            'artifact_type': 'a2g2.artifacts.odyssey',
            'artifact_params': {'path': completed_dir}
        }
        return artifact_spec

    def stringify_exception(self, exception=None):
        return '[%s] %s)' % (type(exception), exception)

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
                    except Exception as read_exception:
                        error = (
                            "Error unknown, unable to read checkpoint file"
                            " '{path}': {read_error}."
                        ).format(
                            path=failed_checkpoint_path,
                            read_error=self.stringify_exception(read_exception)
                        )
                    raise Exception(error)

    def file_exists(self, path=None):
        try:
            self.ssh_client.run_process(cmd=['ls', path], check=True)
            return True
        except:
            return False

    def read_file(self, path=None):
        completed_proc = self.ssh_client.run_process(cmd=['cat', path],
                                                     check=True)
        return completed_proc.stdout

    def get_stdout(self, submission=None, completed_dir=None):
        stdout = ''
        stdout_log_name = submission.get('std_log_files', {}).get('stdout')
        if stdout_log_name:
            stdout_log_path = os.path.join(completed_dir, stdout_log_name)
            stdout = self.read_file(stdout_log_path)
        return stdout
