from collections import defaultdict, OrderedDict
import unittest
from unittest.mock import call, MagicMock, patch

from .. import job_runner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.jobdir_factory = MagicMock()
        self.jobman = MagicMock()
        self.job_runner = job_runner.JobRunner(
            job_record_client=MagicMock(),
            build_jobdir_fn=MagicMock(),
            jobman=MagicMock(),
            jobman_source_name=MagicMock(),
            artifact_handler=MagicMock(),
        )
        self.mod_mocks = {}

    def mockify_job_runner_attrs(self, attrs=None):
        for attr in attrs: setattr(self.job_runner, attr, MagicMock())

    def mockify_module_attrs(self, attrs=None, module=job_runner):
        for attr in attrs:
            patcher = patch.object(module, attr)
            self.addCleanup(patcher.stop)
            self.mod_mocks[attr] = patcher.start()

def TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(attrs=[
            'process_unprocessed_finished_jobs', 'fill_jobman_queue'])
        self.job_runner.tick()

    def test_processes_unprocessed_finished_jobs(self):
        self.assertEqual(
            self.job_runner.process_unprocessed_finished_jobs.call_args,
            call()
        )

    def test_fills_jobman_queue(self):
        self.assertEqual(self.job_runner.fill_jobman_queue.call_args, call())

    def test_returns_tick_stats(self):
        self.fail()

class ProcessFinishedJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(attrs=[
            'get_unprocessed_finished_jobman_jobs', 'parse_jobman_jobs',
            'parsed_jobman_jobs_to_keyed_patches', 'patch_job_records',
            'finalize_jobman_jobs'])
        self.result = self.job_runner.process_finished_jobs()

    def test_parses_jobman_jobs(self):
        self.assertEqual(
            self.job_runner.get_unprocessed_finished_jobman_jobs.call_args,
            call()
        )
        self.assertEqual(
            self.job_runner.parse_jobman_jobs.call_args,
            call(jobman_jobs=(self.job_runner
                              .get_unprocessed_finished_jobman_jobs
                              .return_value))
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
            call(jobman_jobs=(self.job_runner
                              .get_unprocessed_finished_jobman_jobs
                              .return_value))
        )

    def test_returns_stats(self):
        expected_stats = {
            'finished': len(self.job_runner
                            .get_unprocessed_finished_jobman_jobs.return_value)
        }
        self.assertEqual(self.result, expected_stats)

class GetUnprocessedFinishedJobmanJobsTestCase(BaseTestCase):
    def test_dispatches_to_jobman(self):
        result = self.job_runner.get_unprocessed_finished_jobman_jobs()
        self.assertEqual(
            self.job_runner.jobman.get_jobs.call_args,
            call(query={
                'filters': [
                    {'field': 'status', 'op': 'IN',
                     'arg': ['COMPLETED', 'FAILED']},
                    {'field': 'source', 'op': '=',
                     'arg': self.job_runner.jobman_source_name},
                    {'field': 'source_tag', 'op': '! =',
                     'arg': self.job_runner.PROCESSED_TAG},
                ]
            })
        )
        self.assertEqual(result, self.job_runner.jobman.get_jobs.return_value)

class ParseJobmanJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.jobman_jobs = [MagicMock() for i in range(3)]
        self.mockify_job_runner_attrs(attrs=['parse_jobman_job'])
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
        self.mockify_job_runner_attrs(attrs=['_parse_std_log_contents'])
        self.result = self.job_runner.parse_jobman_job(
            jobman_job=self.jobman_job)

    def test_handles_artifact(self):
        self.assertEqual(
            self.job_runner.artifact_handler.dir_to_artifact.call_args,
            call(dir_=self.jobman_job['job_spec']['dir'])
        )

    def test_parses_std_log_contents(self):
        self.assertEqual(self.job_runner._parse_std_log_contents.call_args,
                         call(jobman_job=self.jobman_job))

    def test_returns_expected_parse_results(self):
        expected_result = {
            'jobman_job': self.jobman_job,
            'mc_job': self.jobman_job['source_meta']['mc_job'],
            'artifact': (self.job_runner.artifact_handler.dir_to_artifact
                         .return_value),
            **(self.job_runner._parse_std_log_contents.return_value)
        }
        self.assertEqual(self.result, expected_result)

class ParseStdLogContentsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.jobman_job = MagicMock()
        self.mockify_job_runner_attrs(
            attrs=['get_std_log_contents_for_jobman_job'])

    def _parse_std_log_contents(self, jobman_job=None):
        return self.job_runner._parse_std_log_contents(
            jobman_job=self.jobman_job)

    def test_calls_get_std_log_contents(self):
        self._parse_std_log_contents()
        self.assertEqual(
            self.job_runner.get_std_log_contents_for_jobman_job.call_args,
            call(jobman_job=self.jobman_job)
        )

    def test_result_has_std_logs(self):
        result = self._parse_std_log_contents()
        self.assertEqual(
            result['std_logs'], 
            self.job_runner.get_std_log_contents_for_jobman_job.return_value
        )

    def test_copies_failure_log_content_to_error_key(self):
        failure_msg = 'some_failure'
        self.job_runner.get_std_log_contents_for_jobman_job.return_value = \
                {'failure': failure_msg}
        result = self._parse_std_log_contents()
        self.assertEqual(result['error'], failure_msg)

    def test_sets_status_to_failed_if_error(self):
        self.job_runner.get_std_log_contents_for_jobman_job.return_value = \
                {'failure': MagicMock()}
        result = self._parse_std_log_contents()
        self.assertEqual(result['status'], 'FAILED')

    def test_sets_status_to_completed_if_no_error(self):
        self.job_runner.get_std_log_contents_for_jobman_job.return_value = {}
        result = self._parse_std_log_contents()
        self.assertEqual(result['status'], 'COMPLETED')

class GetJobStdLogFileContentsForJobmanJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(attrs=['read_jobdir_logs'])
        self.job_spec = defaultdict(MagicMock, **{
            'std_log_file_names': {'log_%s' % i: MagicMock() for i in range(3)}
        })
        self.std_logs_to_expose = None

    def _get_std_log_file_contents_for_jobman_job(self):
        jobman_job = defaultdict(MagicMock, **{
            'job_spec': self.job_spec,
            'source_meta': {
                'mc_job': {
                    'cfg': {'std_logs_to_expose': self.std_logs_to_expose}
                }
            }
        })
        return self.job_runner.get_std_log_contents_for_jobman_job(
            jobman_job=jobman_job)

    def test_handles_all_keyword(self):
        self.std_logs_to_expose = 'all'
        result = self._get_std_log_file_contents_for_jobman_job()
        self.assertEqual(
            self.job_runner.read_jobdir_logs.call_args,
            call(job_spec=self.job_spec,
                 logs=[*self.job_spec['std_log_file_names'].keys(),
                       'failure'])
        )
        self.assertEqual(result,
                         self.job_runner.read_jobdir_logs.return_value)

    def test_reads_specified_logs(self):
        self.std_logs_to_expose = \
                list(self.job_spec['std_log_file_names'].keys())[:-1]
        result = self._get_std_log_file_contents_for_jobman_job()
        self.assertEqual(
            self.job_runner.read_jobdir_logs.call_args,
            call(job_spec=self.job_spec, logs=self.std_logs_to_expose)
        )
        self.assertEqual(result,
                         self.job_runner.read_jobdir_logs.return_value)

class ReadJobdirLogsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(attrs=['read_jobdir_log'])

    def test_dispatches_to_read_jobdir_log_file(self):
        job_spec = MagicMock()
        logs = [MagicMock() for i in range(3)]
        result = self.job_runner.read_jobdir_logs(job_spec=job_spec,
                                                      logs=logs)
        expected_call_args_list = [call(job_spec=job_spec, log=log)
                                   for log in logs]
        self.assertEqual(self.job_runner.read_jobdir_log.call_args_list,
                         expected_call_args_list)
        expected_result = {
            log: self.job_runner.read_jobdir_log.return_value
            for log in logs
        }
        self.assertEqual(result, expected_result)

class ReadJobdirLogTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_module_attrs(attrs=['os', 'open'])

    def test_reads_log_file(self):
        job_spec = MagicMock()
        log = MagicMock()
        self.job_runner.read_jobdir_log(job_spec=job_spec, log=log)
        expected_job_file_path = self.mod_mocks['os'].path.join(
            job_spec['dir'], job_spec['std_log_files'][log])
        self.assertEqual(self.mod_mocks['open'].call_args,
                         call(expected_job_file_path))

class ParsedJobmanJobsToKeyedPatchesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(attrs=['parsed_jobman_job_to_patch'])
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
        self.mockify_job_runner_attrs(attrs=['get_std_log_contents'])
        self.parsed_jobman_job = MagicMock()
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

    def test_saves_jobs_with_processed_tag_and_purgeable(self):
        expected_finalized_jobs = [
            {**jobman_job, 'source_tag': self.job_runner.PROCESSED_TAG,
             'purgeable': 1}
            for jobman_job in self.jobman_jobs
        ]
        self.assertEqual(self.job_runner.jobman.save_jobs.call_args,
                         call(jobs=expected_finalized_jobs))

class FillJobmanQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(
            attrs=['submit_mc_job', 'get_claim_limit'])
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

    def test_submits_mc_jobs(self):
        self.assertEqual(self.job_runner.submit_mc_job.call_args_list,
                         [call(mc_job=mc_job) for mc_job in self.mc_jobs])

    def test_returns_stats(self):
        expected_num_claimed = len(
            self.job_runner.job_record_client.claim_job_records.return_value)
        expected_claimed_stats = {'claimed': expected_num_claimed,
                                  'submitted': expected_num_claimed,
                                  'failed': 0}
        self.assertEqual(self.result, expected_claimed_stats)

class SubmitMcJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_job_runner_attrs(attrs=['build_jobdir'])
        self.mc_job = MagicMock()
        self.result = self.job_runner.submit_mc_job(mc_job=self.mc_job)

    def test_builds_jobdir(self):
        self.assertEqual(self.job_runner.build_jobdir.call_args,
                         call(mc_job=self.mc_job))

    def test_submits_job_spec(self):
        self.assertEqual(
            self.job_runner.jobman.submit_job_spec.call_args,
            call(job_spec=self.job_runner.build_jobdir.return_value,
                 source=self.job_runner.jobman_source_name,
                 source_meta={'mc_job': self.mc_job})
        )

class BuildMcJobdirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mc_job = MagicMock()
        self.mockify_job_runner_attrs(
            attrs=['prepare_job_inputs', 'get_jobdir_prefix'])
        self.mockify_module_attrs(attrs=['tempfile'])
        self.result = self.job_runner.build_jobdir(mc_job=self.mc_job)

    def test_creates_jobdir(self):
        self.assertEqual(
            self.mod_mocks['tempfile'].mkdtemp.call_args,
            call(dir=self.job_runner.jobdirs_dir,
                 prefix=self.job_runner.get_jobdir_prefix.return_value)
        )

    def test_prepares_inputs(self):
        self.assertEqual(
            self.job_runner.prepare_job_inputs.call_args,
            call(mc_job=self.mc_job,
                 jobdir=job_runner.tempfile.mkdtemp.return_value)
        )

    def test_dispatches_to_build_jobdir_fn(self):
        self.assertEqual(
            self.job_runner.build_jobdir_fn.call_args,
            call(job=self.mc_job,
                 output_dir=job_runner.tempfile.mkdtemp.return_value)
        )
        self.assertEqual(self.result,
                         self.job_runner.build_jobdir_fn.return_value)

class PrepareJobInputsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.artifacts = {'artifact_%s' % i: MagicMock() for i in range(3)}
        self.mc_job = {'job_inputs': {'artifacts': self.artifacts}}
        self.jobdir = 'some_jobdir'
        self.mockify_module_attrs(attrs=['os'])
        self.expected_inputs_dir = self.mod_mocks['os'].path.join(
            self.jobdir, 'inputs')
        self.job_runner.prepare_job_inputs(mc_job=self.mc_job,
                                           jobdir=self.jobdir)

    def test_makes_inputs_dir(self):
        self.assertEqual(self.mod_mocks['os'].makedirs.call_args,
                         call(self.expected_inputs_dir, exist_ok=True))

    def test_calls_artifact_to_dir_for_each_artifact(self):
        expected_call_args_list =[]
        for artifact_key, artifact in self.artifacts.items():
            dest = job_runner.os.path.join(
                self.expected_inputs_dir, artifact_key)
            expected_call_args_list.append(call(artifact=artifact, dest=dest))
        self.assertEqual(
            self.job_runner.artifact_handler.artifact_to_dir.call_args_list,
            expected_call_args_list
        )
