from collections import defaultdict, OrderedDict
import json
import unittest
from unittest.mock import call, MagicMock, patch

from .. import job_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.submission_factory = MagicMock()
        self.jobman = MagicMock()
        self.job_runner = job_runner.JobRunner(
            job_client=MagicMock(),
            submission_factory=MagicMock(),
            jobman=MagicMock(),
            jobman_source_name=MagicMock()
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

class ProcessExecutedJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        for attr in ['get_executed_jobman_jobs', 'parse_jobman_jobs',
                     'parsed_jobman_jobs_to_keyed_patches', 'patch_jobs',
                     'finalize_jobman_jobs']:
            setattr(self.job_runner, attr, MagicMock())
        self.job_runner.process_executed_jobs()

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
        self.assertEqual(self.job_runner.patch_jobs.call_args,
                         call(keyed_patches=expected_keyed_patches))

    def test_finalizes_jobman_jobs(self):
        self.assertEqual(
            self.job_runner.finalize_jobman_jobs.call_args,
            call(jobman_jobs=\
                 self.job_runner.get_executed_jobman_jobs.return_value))

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
        for attr in ['generate_artifact_spec_for_dir',
                     'get_std_log_contents_for_jobman_job']:
            setattr(self.job_runner, attr, MagicMock())

    def _parse(self):
        return self.job_runner.parse_jobman_job(jobman_job=self.jobman_job)

    def test_has_jobman_job(self):
        result = self._parse()
        self.assertEqual(result['jobman_job'], self.jobman_job)

    def test_has_artifact_spec(self):
        result = self._parse()
        self.assertEqual(
            self.job_runner.generate_artifact_spec_for_dir.call_args,
            call(_dir=self.jobman_job['submission']['dir']))
        self.assertEqual(
            result['artifact_spec'],
            self.job_runner.generate_artifact_spec_for_dir.return_value)

    def test_has_std_log_contents(self):
        result = self._parse()
        self.assertEqual(
            self.job_runner.get_std_log_contents_for_jobman_job.call_args,
            call(jobman_job=self.jobman_job)
        )
        self.assertEqual(
            result['std_log_contents'],
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

class GenerateArtifactSpecForDirTestCase(BaseTestCase):
    def test_generates_artifact_spec(self):
        _dir = MagicMock()
        result = self.job_runner.generate_artifact_spec_for_dir(_dir=_dir)
        expected_result = {
            'artifact_type': 'a2g2.artifacts.odyssey',
            'artifact_params': {'path': _dir}
        }
        self.assertEqual(result, expected_result)


class GetJobStdLogFileContentsForJobmanJob(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_runner.read_submission_logs = MagicMock()
        self.submission = defaultdict(MagicMock, **{
            'std_log_files': {'log_%s' % i: MagicMock() for i in range(3)}
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
                 logs=self.submission['std_log_files'].keys())
        )
        self.assertEqual(result,
                         self.job_runner.read_submission_logs.return_value)

    def test_reads_specified_logs(self):
        self.std_logs_to_expose = \
                list(self.submission['std_log_files'].keys())[:-1]
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

    def test_keys_patches_by_mc_job_uuid(self):
        result = self.job_runner.parsed_jobman_jobs_to_keyed_patches(
            parsed_jobman_jobs=self.parsed_jobman_jobs)
        expected_result = {}
        for parsed_jobman_job in self.parsed_jobman_jobs:
            mc_job = parsed_jobman_job['jobman_job']['source_meta']['mc_job']
            expected_result[mc_job['uuid']] = \
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
        self.assertEqual(self.patch['data']['std_log_contents'],
                         self.parsed_jobman_job.get('std_log_contents'))


class PatchJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_runner.serialize_job_patch = MagicMock()
        self.keyed_patches = OrderedDict()
        for i in range(3): self.keyed_patches['key_%s' % i] = MagicMock()
        self.job_runner.patch_jobs(keyed_patches=self.keyed_patches)

    def test_dispatches_serialized_patches_to_job_client(self):
        self.assertEqual(
            self.job_runner.serialize_job_patch.call_args_list,
            [call(patch=patch) for patch in self.keyed_patches.values()]
        )
        self.assertEqual(
            self.job_runner.job_client.patch_jobs.call_args,
            call(keyed_patches={
                key: self.job_runner.serialize_job_patch.return_value
                for  key in self.keyed_patches
            })
        )

class SerializeJobPatchTestCase(BaseTestCase):
    def test_serializes_data_field(self):
        patch = {'data': {'some': 'data'}}
        result = self.job_runner.serialize_job_patch(patch=patch)
        self.assertEqual(result, {'data': json.dumps(patch['data'])})

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
        for attr in ['build_submission', 'get_claim_limit']:
            setattr(self.job_runner, attr, MagicMock())
        self.mc_jobs = [MagicMock() for i in range(3)]
        self.job_runner.job_client.claim_jobs.return_value = self.mc_jobs
        self.job_runner.fill_jobman_queue()

    def test_claims_jobs(self):
        self.assertEqual(
            self.job_runner.job_client.claim_jobs.call_args,
            call(params={
                'limit': self.job_runner.get_claim_limit.return_value
            })
        )

    def test_builds_submissions(self):
        self.assertEqual(self.job_runner.build_submission.call_args_list,
                         [call(mc_job=mc_job) for mc_job in self.mc_jobs])

    def test_submits_submissions(self):
        self.assertEqual(
            self.job_runner.jobman.submit_job.call_args_list,
            [
                call(submission=self.job_runner.build_submission.return_value,
                     source_meta={'mc_job': mc_job},
                     source=self.job_runner.jobman_source_name)
                for mc_job in self.mc_jobs
            ]
        )

class BuildSubmissionTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_job = MagicMock()
        for attr in ['prepare_job_inputs']:
            setattr(self.job_runner, attr, MagicMock())
        patcher = patch.object(job_runner, 'tempfile')
        patcher.start()
        self.result = self.job_runner.build_submission(mc_job=self.mc_job)
        self.addCleanup(patcher.stop)

    def test_creates_submission_dir(self):
        self.assertEqual(job_runner.tempfile.mkdtemp.call_args,
                         call(prefix='sf.'))

    def test_prepares_inputs(self):
        self.assertEqual(
            self.job_runner.prepare_job_inputs.call_args,
            call(mc_job=self.mc_job,
                 submission_dir=job_runner.tempfile.mkdtemp.return_value)
        )

    def test_dispatches_to_submission_factory(self):
        submission_factory = self.job_runner.submission_factory
        self.assertEqual(
            submission_factory.build_submission.call_args,
            call(job=self.mc_job,
                 output_dir=job_runner.tempfile.mkdtemp.return_value)
        )

class PrepareJobInputsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.artifacts = {i: MagicMock() for i in range(3)}
        self.mc_job = {
            'job_spec': {
                'inputs': {
                    'artifacts': self.artifacts
                }
            }
        }
        self.submission_dir = MagicMock()
        patcher = patch.object(job_runner, 'os')
        patcher.start()
        self.addCleanup(patcher.stop)
        for attr in ['prepare_input_artifact']:
            setattr(self.job_runner, attr, MagicMock())
        self.job_runner.prepare_job_inputs(
            mc_job=self.mc_job, submission_dir=self.submission_dir)

    def test_makes_inputs_dir(self):
        self.assertEqual(job_runner.os.path.join.call_args,
                         call(self.submission_dir, 'inputs'))
        self.assertEqual(
            job_runner.os.makedirs.call_args,
            call(job_runner.os.path.join.return_value, exist_ok=True))

    def test_calls_prepare_input_artifact_for_each_artifact(self):
        self.assertEqual(
            self.job_runner.prepare_input_artifact.call_args_list,
            [call(artifact_key=artifact_key,
                  artifact=artifact,
                  inputs_dir=job_runner.os.path.join.return_value)
             for artifact_key, artifact in self.artifacts.items()
            ]
        )
