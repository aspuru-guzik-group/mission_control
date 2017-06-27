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
        expected_last_line = 'No more flows to claimed.'
        expected_stdout_lines = textwrap.dedent(
            '''
            I am task 'task_0' in flow 'flow_0'
            I am task 'task_0' in flow 'flow_2'
            I am task 'task_0' in flow 'flow_1'
            I am task 'task_1' in flow 'flow_0'
            I am task 'task_1' in flow 'flow_2'
            I am task 'task_1' in flow 'flow_1'
            I am task 'task_2' in flow 'flow_0'
            I am task 'task_2' in flow 'flow_2'
            I am task 'task_2' in flow 'flow_1'
            '''
        ).lstrip().split("\n") + [expected_last_line]
        stdout_lines = stdout.split("\n")
        self.assertEqual(sorted(stdout_lines), sorted(expected_stdout_lines))
