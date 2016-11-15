import os
import unittest

from ..odyssey_client import OdysseyClient
from ..ssh import SSHControlSocketClient


@unittest.skipIf(os.environ.get('ODYSSEY_SSH_CONTROL_SOCKET') is None,
                 "'ODYSSEY_SSH_CONTROL_SOCKET' was not set")
@unittest.skipIf(os.environ.get('ODYSSEY_HOST') is None,
                 "'ODYSSEY_HOST' was not set")
class JobRunnerTestCase(unittest.TestCase):
    def setUp(self):
        self.ssh_client = SSHControlSocketClient(
            user=os.environ.get('ODYSSEY_USER', os.environ.get('USER')),
            host=os.environ.get('ODYSSEY_HOST'),
            control_socket=os.environ.get('ODYSSEY_SSH_CONTROL_SOCKET'),
        )
        self.ssh_client.connect()
        self.odyssey_client = OdysseyClient(ssh_client=self.ssh_client)
        _dir = os.path.dirname(os.path.abspath(__file__))
        self.snapshots_dir = os.path.join(_dir, 'snapshots')
        self.fixtures_dir = os.path.join(_dir, 'fixtures')
        self.job_dir = os.path.join(self.snapshots_dir, 'stub_job')

    def test_upload_job(self):
        remote_job_dir = self.odyssey_client._upload_job_dir(self.job_dir)
        remote_dir_contents = self.ssh_client.run(
            ['ls -1',  remote_job_dir]).stdout.split()
        expected_contents = os.listdir(self.job_dir)
        self.assertEqual(remote_dir_contents, expected_contents)

    def test_start_job(self):
        job_metadata = self.odyssey_client.start_job(self.job_dir)
        self.assertTrue(job_metadata['job_id'] is not None)

    def test_get_job_state(self):
        job_metadata = self.odyssey_client.start_job(self.job_dir)
        job_state = self.odyssey_client.get_job_state(job_metadata)
        self.assertTrue(job_state is not None)

    def test_parse_raw_job_state(self):
        with open(os.path.join(self.fixtures_dir, 'raw_job_state.txt')) as f:
            raw_job_state = f.read()
        job_state = self.odyssey_client._parse_raw_job_state(raw_job_state)
        self.assertEqual(job_state['JobId'], '71927594')

if __name__ == '__main__':
    unittest.main()
