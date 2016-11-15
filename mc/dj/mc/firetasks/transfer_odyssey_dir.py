import tempfile
import os
from fireworks.core.firework import Firework, FWAction, FireTaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from ..utils.ssh_client import SSHControlSocketClient

from .process_completed_job import ProcessCompletedJobFireTask

@explicit_serialize
class TransferOdysseyDirFireTask(FireTaskBase):
    def run_task(self, fw_spec):
        target_dir = tempfile.mkdtemp(prefix='ody.output.')
        self.transfer_odyssey_dir(fw_spec['odyssey_dir'], target_dir)
        ftask = ProcessCompletedJobFireTask()
        fwork = Firework([ftask])
        fwork.spec['job_dir'] = target_dir
        fwork.spec['_category'] = 'completed_job_processing'
        return FWAction(additions=[fwork])

    def transfer_odyssey_dir(self, odyssey_dir, target_dir):
        ssh_client = self.get_ssh_client()
        ssh_client.scp_from(
            src=odyssey_dir,
            dst=target_dir,
            flags='-r')

    def get_ssh_client(self):
        ssh_client = SSHControlSocketClient(
            user=os.environ.get('ODYSSEY_USER', os.environ.get('USER')),
            host=os.environ.get('ODYSSEY_HOST'),
            control_socket=os.environ.get('ODYSSEY_SSH_CONTROL_SOCKET'),
        )
        ssh_client.connect()
        return ssh_client
