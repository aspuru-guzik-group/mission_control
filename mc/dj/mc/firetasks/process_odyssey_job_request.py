import os
import tempfile
import textwrap

from fireworks.core.firework import Firework, FWAction, FireTaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from .poll_odyssey_job import PollOdysseyJobFireTask
from ..utils.ssh_client import SSHControlSocketClient
from ..utils.odyssey_client import OdysseyClient

@explicit_serialize
class ProcessOdysseyJobRequestFireTask(FireTaskBase):
   def run_task(self, fw_spec):
       odyssey_job_dir = self.build_odyssey_job_dir()
       odyssey_job_metadata = self.submit_odyssey_job_dir(odyssey_job_dir)
       ftask = PollOdysseyJobFireTask()
       fwork = Firework([ftask])
       fwork.spec['odyssey_job_metadata'] = odyssey_job_metadata
       fwork.spec['_category'] = 'odyssey_job_polls'
       return FWAction(additions=[fwork])

   def build_odyssey_job_dir(self):
       odyssey_job_dir = tempfile.mkdtemp(prefix='ody.')
       job_file_path = os.path.join(odyssey_job_dir, 'job.sh')
       with open(job_file_path, 'w') as f:
           f.write(textwrap.dedent(
               """
               #!/bin/bash
               date >> ~/ody_test.txt
               echo start >> ~/ody_test.txt
               sleep 20
               date >> ~/ody_test.txt
               echo end >> ~/ody_test.txt
               """).strip())
       return odyssey_job_dir

   def submit_odyssey_job_dir(self, odyssey_job_dir):
       odyssey_client = self.get_odyssey_client()
       odyssey_job_metadata = odyssey_client.start_job(odyssey_job_dir)
       return odyssey_job_metadata

   def get_odyssey_client(self):
       ssh_client = SSHControlSocketClient(
           user=os.environ.get('ODYSSEY_USER', os.environ.get('USER')),
           host=os.environ.get('ODYSSEY_HOST'),
           control_socket=os.environ.get('ODYSSEY_SSH_CONTROL_SOCKET'),
       )
       ssh_client.connect()
       odyssey_client = OdysseyClient(ssh_client=ssh_client)
       return odyssey_client
