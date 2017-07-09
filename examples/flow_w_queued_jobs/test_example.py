import os
import subprocess
import textwrap
import unittest


class TestExample(unittest.TestCase):
    def test_example(self):
        cmd = 'cd {this_dir} && python entrypoint.py'.format(
            this_dir=os.path.dirname(os.path.abspath(__file__))
        )
        stdout = subprocess.check_output(cmd, shell=True).decode()
        expected_stdout = textwrap.dedent(
            '''
            running job 'job.0'
            running job 'job.2'
            running job 'job.1'
            No more flows to claim.
            '''
        ).lstrip()
        self.assertEqual(set(stdout.split("\n")),
                         set(expected_stdout.split("\n")))
