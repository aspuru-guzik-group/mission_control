from collections import defaultdict
import textwrap
import unittest
from unittest.mock import call, MagicMock

from ..slurm_execution_client import SlurmExecutionClient


class SlurmExecutionClientBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.process_runner = MagicMock()
        self.slurm_client = SlurmExecutionClient(
            process_runner=self.process_runner
        )
        self.submission = defaultdict(MagicMock)
        self.execution_meta = defaultdict(MagicMock)

    def generate_failed_proc(self):
        proc = MagicMock()
        proc.returncode = 1
        proc.stderr = 'some error'
        return proc

class StartExecutionTestCase(SlurmExecutionClientBaseTestCase):
    def test_calls_sbatch(self):
        self.process_runner.run_process.return_value = \
                self.generate_successful_sbatch_proc()
        self.slurm_client.start_execution(submission=self.submission)
        workdir = self.submission['dir']
        entrypoint = workdir + '/' + self.submission['entrypoint']
        expected_cmd = ['sbatch', '--workdir=%s' % workdir, entrypoint]
        self.assertEqual(self.process_runner.run_process.call_args,
                         call(cmd=expected_cmd, check=True))

    def generate_successful_sbatch_proc(self, job_id='12345'):
        proc = MagicMock()
        proc.returncode = 0
        proc.stdout = 'Submitted batch job %s' % job_id
        return proc

    def test_returns_expected_execution_meta_for_successful_submission(self):
        job_id = '12345'
        self.process_runner.run_process.return_value = \
                self.generate_successful_sbatch_proc(job_id=job_id)
        execution_meta = self.slurm_client.start_execution(
            submission=self.submission)
        expected_execution_meta = {
            'job_id': job_id,
            'submission': self.submission,
        }
        self.assertEqual(execution_meta, expected_execution_meta)

    def test_handles_failed_submission(self):
        failed_proc = self.generate_failed_proc()
        self.process_runner.run_cmd.return_value = failed_proc
        with self.assertRaises(Exception) as context:
            self.slurm_client.start_execution(job=self.job)
            self.assertEqual(str(context.exception), failed_proc.stderr)

class GetExecutionStateTestCase(SlurmExecutionClientBaseTestCase):
    def test_calls_scontrol(self):
        self.process_runner.run_process.return_value = \
                self.generate_completed_scontrol_proc()
        self.slurm_client.get_execution_state(
            execution_meta=self.execution_meta)
        expected_cmd = ['scontrol', 'show', '--details', '--oneliner',
                       'job', self.execution_meta['job_id']]
        self.assertEqual(self.process_runner.run_process.call_args,
                         call(cmd=expected_cmd, check=True))

    def generate_completed_scontrol_proc(self, job_state='RUNNING'):
        proc = MagicMock()
        proc.returncode = 0
        proc.stdout = textwrap.dedent(
            """
            JobId=75945797 JobName=foo.sh UserId=user(5900165)
            GroupId=rc_admin(40273) MCS_label=N/A Priority=19999669
            Nice=0 Account=rc_admin QOS=normal JobState=COMPLETED
            Reason=None Dependency=(null) Requeue=1 Restarts=0
            BatchFlag=1 Reboot=0 ExitCode=0:0 RunTime=00:00:01
            TimeLimit=00:10:00 TimeMin=N/A
            SubmitTime=2016-11-23T13:44:31
            EligibleTime=2016-11-23T13:44:31
            StartTime=2016-11-23T13:44:32 EndTime=2016-11-23T13:44:33
            Deadline=N/A PreemptTime=None SuspendTime=None
            SecsPreSuspend=0 Partition=serial_requeue
            AllocNode:Sid=rclogin03:9811 ReqNodeList=(null)
            ExcNodeList=(null) NodeList=holy2a16206
            BatchHost=holy2a16206 NumNodes=1 NumCPUs=1 NumTasks=0
            CPUs/Task=1 ReqB:S:C:T=0:0:*:* TRES=cpu=1,mem=100M,node=1
            Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
            MinCPUsNode=1 MinMemoryCPU=100M MinTmpDiskNode=0
            Features=(null) Gres=(null) Reservation=(null)
            OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
            Command=/home/user/foo.sh WorkDir=/home/user/foo
            StdErr=/home/user/slurm-75945797.out StdIn=/dev/null
            StdOut=/home/user/slurm-75945797.out Power=
            """)
        return proc

    def test_completed_proc(self):
        self.process_runner.run_process.return_value = \
                self.generate_completed_scontrol_proc()
        execution_state = self.slurm_client.get_execution_state(
            execution_meta=self.execution_meta)
        self.assertEqual(execution_state['run_status'], 'COMPLETED')
        self.assertEqual(execution_state['completed_dir'],
                        self.execution_meta['submission']['dir'])
        self.assertTrue('slurm_job_meta' in execution_state)

    def test_failed_proc(self):
        failed_proc = self.generate_failed_proc()
        self.process_runner.run_cmd.return_value = failed_proc
        with self.assertRaises(Exception) as context:
            self.slurm_client.get_execution_sate(job=self.job)
            self.assertEqual(str(context.exception), failed_proc.stderr)

if __name__ == '__main__':
    unittest.main()
