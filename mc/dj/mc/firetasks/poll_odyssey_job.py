import os
from fireworks.core.firework import Firework, FWAction, FireTaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from ..utils.ssh_client import SSHControlSocketClient
from ..utils.odyssey_client import OdysseyClient

from .transfer_odyssey_dir import TransferOdysseyDirFireTask

@explicit_serialize
class PollOdysseyJobFireTask(FireTaskBase):
    def run_task(self, fw_spec):
        job_state = self.get_job_state(fw_spec['odyssey_job_metadata'])
        print("job_state: ", job_state)
        if job_state in ['PENDING', 'RUNNING']:
            ftask = PollOdysseyJobFireTask()
            fwork = Firework([ftask])
            fwork.spec['odyssey_job_metadata'] = fw_spec['odyssey_job_metadata']
            fwork.spec['_category'] = 'odyssey_job_polls'
        else:
            ftask = TransferOdysseyDirFireTask()
            fwork = Firework([ftask])
            fwork.spec['odyssey_dir'] = \
                    fw_spec['odyssey_job_metadata']['remote_dir']
            fwork.spec['_category'] = 'odyssey_transfers'
        return FWAction(additions=[fwork])

    def get_job_state(self, job_metadata):
        odyssey_client = self.get_odyssey_client()
        return odyssey_client.get_job_state(job_metadata)['JobState']

    def get_odyssey_client(self):
        ssh_client = SSHControlSocketClient(
            user=os.environ.get('ODYSSEY_USER', os.environ.get('USER')),
            host=os.environ.get('ODYSSEY_HOST'),
            control_socket=os.environ.get('ODYSSEY_SSH_CONTROL_SOCKET'),
        )
        ssh_client.connect()
        odyssey_client = OdysseyClient(ssh_client=ssh_client)
        return odyssey_client
