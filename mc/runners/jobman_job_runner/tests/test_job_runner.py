from collections import defaultdict, OrderedDict
import unittest
from unittest.mock import call, MagicMock, patch

from .. import job_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.submission_factory = MagicMock()
        self.jobman = MagicMock()
        self.job_runner = job_runner.JobRunner(
            job_record_client=MagicMock(),
            submission_factory=MagicMock(),
            cfg_factory=MagicMock(),
            jobman=MagicMock(),
            jobman_source_name=MagicMock(),
            artifact_processor=MagicMock(),
        )

def TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        for attr in ['process_executed_jobs, fill_jobman_queue']:
            setattr(self.job_runner, attr, MagicMock())
        self.job_runner.tick()

    def test_processes_executed_jobs(self):
        self.assertEqual(self.job_runner.process_executed_jobs.call_args,
                         call())

    def test_fills_jobman_queue(self):
        self.assertEqual(self.job_runner.fill_jobman_queue.call_args, call())

    def test_returns_tick_stats(self):
        self.fail()

class ProcessExecutedJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        for attr in ['get_executed_jobman_jobs', 'parse_jobman_jobs',
                     'parsed_jobman_jobs_to_keyed_patches', 'patch_job_records',
                     'finalize_jobman_jobs']:
            setattr(self.job_runner, attr, MagicMock())
        self.result = self.job_runner.process_executed_jobs()

    def test_parses_executed_jobman_jobs(self):
        self.assertEqual(self.job_runner.get_executed_jobman_jobs.call_args,
                         call())
        self.assertEqual(
            self.job_runner.parse_jobman_jobs.call_args,
            call(jobman_jobs=\
                 self.job_runner.get_executed_jobman_jobs.return_value)
        )

    def test_posts_keyed_patches(self):
        expected_parsed_jobman_jobs = \
                self.job_runner.parse_jobman_jobs.return_value
        self.assertEqual(
            self.job_runner.parsed_jobman_jobs_to_keyed_patches.call_args,
            call(parsed_jobman_jobs=expected_parsed_jobman_jobs)
        )
        expected_keyed_patches = \
                self.job_runner.parsed_jobman_jobs_to_keyed_patches.return_value
        self.assertEqual(self.job_runner.patch_job_records.call_args,
                         call(keyed_patches=expected_keyed_patches))

    def test_finalizes_jobman_jobs(self):
        self.assertEqual(
            self.job_runner.finalize_jobman_jobs.call_args,
            call(jobman_jobs=\
                 self.job_runner.get_executed_jobman_jobs.return_value))

    def test_returns_stats(self):
        expected_stats = {
            'executed': len(
                self.job_runner.get_executed_jobman_jobs.return_value)
        }
        self.assertEqual(self.result, expected_stats)

class GetExecutedJobmanJobsTestCase(BaseTestCase):
    def test_dispatches_to_jobman(self):
        result = self.job_runner.get_executed_jobman_jobs()
        self.assertEqual(
            self.job_runner.jobman.get_jobs.call_args,
            call(query={
                'filters': [
                    {'field': 'status', 'operator': 'IN',
                     'value': ['EXECUTED', 'FAILED']},
                    {'field': 'source', 'operator': '=',
                     'value': self.job_runner.jobman_source_name},
                ]
            })
        )
        self.assertEqual(result, self.job_runner.jobman.get_jobs.return_value)

class ParseJobmanJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.jobman_jobs = [MagicMock() for i in range(3)]
        for attr in ['parse_jobman_job']:
            setattr(self.job_runner, attr, MagicMock())
        self.result = self.job_runner.parse_jobman_jobs(
            jobman_jobs=self.jobman_jobs)

    def test_dispatches_to_parse_jobman_job(self):
        self.assertEqual(self.job_runner.parse_jobman_job.call_args_list,
                         [call(jobman_job=jobman_job)
                          for jobman_job in self.jobman_jobs])
        self.assertEqual(self.result,
                         [self.job_runner.parse_jobman_job.return_value
                          for job in self.jobman_jobs])

class ParseJobmanJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.jobman_job = MagicMock()
        for attr in ['get_std_log_contents_for_jobman_job']:
            setattr(self.job_runner, attr, MagicMock())

    def _parse(self):
        return self.job_runner.parse_jobman_job(jobman_job=self.jobman_job)

    def test_has_jobman_job(self):
        result = self._parse()
        self.assertEqual(result['jobman_job'], self.jobman_job)

    def test_has_artifact_spec(self):
        result = self._parse()
        self.assertEqual(
            self.job_runner.artifact_processor.dir_to_artifact.call_args,
            call(dir_=self.jobman_job['submission']['dir']))
        self.assertEqual(
            result['artifact'],
            self.job_runner.artifact_processor.dir_to_artifact.return_value)

    def test_has_std_log_contents(self):
        result = self._parse()
        self.assertEqual(
            self.job_runner.get_std_log_contents_for_jobman_job.call_args,
            call(jobman_job=self.jobman_job)
        )
        self.assertEqual(
            result['std_logs'],
            self.job_runner.get_std_log_contents_for_jobman_job.return_value
        )

    def test_copies_failure_log_content_to_error_key(self):
        failure_msg = 'some_failure'
        self.job_runner.get_std_log_contents_for_jobman_job.return_value = \
                {'failure': failure_msg}
        result = self._parse()
        self.assertEqual(result['error'], failure_msg)

    def test_sets_status_to_failed_if_error(self):
        self.job_runner.get_std_log_contents_for_jobman_job.return_value = \
                {'failure': MagicMock()}
        result = self._parse()
        self.assertEqual(result['status'], 'FAILED')

    def test_sets_status_to_completed_if_no_error(self):
        self.job_runner.get_std_log_contents_for_jobman_job.return_value = {}
        result = self._parse()
        self.assertEqual(result['status'], 'COMPLETED')

class GetJobStdLogFileContentsForJobmanJob(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_runner.read_submission_logs = MagicMock()
        self.submission = defaultdict(MagicMock, **{
            'std_log_file_names': {'log_%s' % i: MagicMock() for i in range(3)}
        })
        self.std_logs_to_expose = None

    def _get_std_log_file_contents_for_jobman_job(self):
        jobman_job = defaultdict(MagicMock, **{
            'submission': self.submission,
            'source_meta': {
                'mc_job': {
                    'job_spec': {'std_logs_to_expose': self.std_logs_to_expose}
                }
            }
        })
        return self.job_runner.get_std_log_contents_for_jobman_job(
            jobman_job=jobman_job)

    def test_handles_all_keyword(self):
        self.std_logs_to_expose = 'all'
        result = self._get_std_log_file_contents_for_jobman_job()
        self.assertEqual(
            self.job_runner.read_submission_logs.call_args,
            call(submission=self.submission,
                 logs=self.submission['std_log_file_names'].keys())
        )
        self.assertEqual(result,
                         self.job_runner.read_submission_logs.return_value)

    def test_reads_specified_logs(self):
        self.std_logs_to_expose = \
                list(self.submission['std_log_file_names'].keys())[:-1]
        result = self._get_std_log_file_contents_for_jobman_job()
        self.assertEqual(
            self.job_runner.read_submission_logs.call_args,
            call(submission=self.submission, logs=self.std_logs_to_expose)
        )
        self.assertEqual(result,
                         self.job_runner.read_submission_logs.return_value)

class ReadSubmissionLogsTestCase(BaseTestCase):
    def test_dispatches_to_read_submission_log_file(self):
        self.job_runner.read_submission_log = MagicMock()
        submission = MagicMock()
        logs = [MagicMock() for i in range(3)]
        result = self.job_runner.read_submission_logs(submission=submission,
                                                      logs=logs)
        expected_call_args_list = [call(submission=submission, log=log)
                                   for log in logs]
        self.assertEqual(self.job_runner.read_submission_log.call_args_list,
                         expected_call_args_list)
        expected_result = {
            log: self.job_runner.read_submission_log.return_value
            for log in logs
        }
        self.assertEqual(result, expected_result)

class ReadSubmissionLogTestCase(BaseTestCase):
    @patch.object(job_runner, 'os')
    @patch.object(job_runner, 'open')
    def test_reads_log_file(self, mock_open, mock_os):
        submission = MagicMock()
        log = MagicMock()
        self.job_runner.read_submission_log(submission=submission, log=log)
        expected_job_file_path = mock_os.path.join(
            submission['dir'], submission['std_log_files'][log])
        self.assertEqual(mock_open.call_args, call(expected_job_file_path))

class ParsedJobmanJobsToKeyedPatchesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_runner.parsed_jobman_job_to_patch = MagicMock()
        self.parsed_jobman_jobs = [MagicMock() for i in range(3)]

    def test_keys_patches_by_mc_job_key(self):
        result = self.job_runner.parsed_jobman_jobs_to_keyed_patches(
            parsed_jobman_jobs=self.parsed_jobman_jobs)
        expected_result = {}
        for parsed_jobman_job in self.parsed_jobman_jobs:
            expected_result[parsed_jobman_job['mc_job']['key']] = \
                    self.job_runner.parsed_jobman_job_to_patch.return_value
        self.assertEqual(result, expected_result)

class ParsedJobmanJobToPatchTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.parsed_jobman_job = MagicMock()
        self.job_runner.get_std_log_contents = MagicMock()
        self.patch = self.job_runner.parsed_jobman_job_to_patch(
            parsed_jobman_job=self.parsed_jobman_job)

    def test_has_status(self):
        self.assertEqual(self.patch['status'], self.parsed_jobman_job['status'])

    def test_has_artifact(self):
        self.assertEqual(self.patch['data']['artifact'],
                         self.parsed_jobman_job.get('artifact_spec'))

    def test_has_std_log_contents(self):
        self.assertEqual(self.patch['data']['std_logs'],
                         self.parsed_jobman_job.get('std_logs'))


class PatchJobRecordsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.keyed_patches = OrderedDict()
        for i in range(3): self.keyed_patches['key_%s' % i] = MagicMock()
        self.job_runner.patch_job_records(keyed_patches=self.keyed_patches)

    def test_dispatches_to_job_record_client(self):
        self.assertEqual(
            self.job_runner.job_record_client.patch_job_records.call_args,
            call(keyed_patches=self.keyed_patches)
        )

class FinalizeJobmanJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.jobman_jobs = [defaultdict(MagicMock) for i in range(3)]
        self.result = self.job_runner.finalize_jobman_jobs(
            jobman_jobs=self.jobman_jobs)

    def test_saves_jobs_with_completed_status(self):
        expected_finalized_jobs = [{**jobman_job, 'status': 'COMPLETED'}
                                   for jobman_job in self.jobman_jobs]
        self.assertEqual(self.job_runner.jobman.save_jobs.call_args,
                         call(jobs=expected_finalized_jobs))

class FillJobmanQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        for attr in ['build_job_submission', 'get_claim_limit']:
            setattr(self.job_runner, attr, MagicMock())
        self.mc_jobs = [MagicMock() for i in range(3)]
        self.job_runner.job_record_client.claim_job_records.return_value = \
                self.mc_jobs
        self.result =self.job_runner.fill_jobman_queue()

    def test_claims_jobs(self):
        self.assertEqual(
            self.job_runner.job_record_client.claim_job_records.call_args,
            call(params={
                'limit': self.job_runner.get_claim_limit.return_value
            })
        )

    def test_builds_job_submissions(self):
        self.assertEqual(self.job_runner.build_job_submission.call_args_list,
                         [call(mc_job=mc_job) for mc_job in self.mc_jobs])

    def test_submits_submissions(self):
        self.assertEqual(
            self.job_runner.jobman.submit_job.call_args_list,
            [
                call(
                    submission=self.job_runner.build_job_submission.return_value,
                    source_meta={'mc_job': mc_job},
                    source=self.job_runner.jobman_source_name
                )
                for mc_job in self.mc_jobs
            ]
        )

    def test_returns_stats(self):
        expected_num_claimed = len(
            self.job_runner.job_record_client.claim_job_records.return_value)
        expected_claimed_stats = {'claimed': expected_num_claimed,
                                  'submitted': expected_num_claimed,
                                  'failed': 0}
        self.assertEqual(self.result, expected_claimed_stats)

class BuildMcSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_job = MagicMock()
        for attr in ['prepare_job_inputs', 'get_submission_dir_prefix']:
            setattr(self.job_runner, attr, MagicMock())
        patcher = patch.object(job_runner, 'tempfile')
        self.addCleanup(patcher.stop)
        self.mocks = {'tempfile': patcher.start()}
        self.result = self.job_runner.build_job_submission(mc_job=self.mc_job)

    def test_creates_submission_dir(self):
        self.assertEqual(
            self.mocks['tempfile'].mkdtemp.call_args,
            call(dir=self.job_runner.submissions_dir,
                 prefix=self.job_runner.get_submission_dir_prefix.return_value)
        )

    def test_prepares_cfg(self):
        self.assertEqual(self.job_runner.cfg_factory.get_cfg.call_args,
                         call(job=self.mc_job))

    def test_prepares_inputs(self):
        self.assertEqual(
            self.job_runner.prepare_job_inputs.call_args,
            call(mc_job=self.mc_job,
                 cfg=self.job_runner.cfg_factory.get_cfg.return_value,
                 submission_dir=job_runner.tempfile.mkdtemp.return_value)
        )

    def test_dispatches_to_submission_factory(self):
        submission_factory = self.job_runner.submission_factory
        self.assertEqual(
            submission_factory.build_job_submission.call_args,
            call(job=self.mc_job,
                 cfg=self.job_runner.cfg_factory.get_cfg.return_value,
                 output_dir=job_runner.tempfile.mkdtemp.return_value)
        )
        self.assertEqual(self.result,
                         submission_factory.build_job_submission.return_value)

class PrepareJobInputsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.artifacts = {'artifact_%s' % i: MagicMock() for i in range(3)}
        self.mc_job = {
            'job_spec': {
                'inputs': {
                    'artifacts': self.artifacts
                }
            }
        }
        self.submission_dir = 'some_submission_dir'
        patcher = patch.object(job_runner.os, 'makedirs')
        patcher.start()
        self.addCleanup(patcher.stop)
        self.expected_inputs_dir = job_runner.os.path.join(
            self.submission_dir, 'inputs')
        self.job_runner.prepare_job_inputs(
            mc_job=self.mc_job, submission_dir=self.submission_dir)

    def test_makes_inputs_dir(self):
        self.assertEqual(job_runner.os.makedirs.call_args,
                         call(self.expected_inputs_dir, exist_ok=True))

    def test_calls_artifact_to_dir_for_each_artifact(self):
        expected_call_args_list =[]
        for artifact_key, artifact in self.artifacts.items():
            dest = job_runner.os.path.join(
                self.expected_inputs_dir, artifact_key)
            expected_call_args_list.append(call(artifact=artifact, dest=dest))
        self.assertEqual(
            self.job_runner.artifact_processor.artifact_to_dir.call_args_list,
            expected_call_args_list
        )
