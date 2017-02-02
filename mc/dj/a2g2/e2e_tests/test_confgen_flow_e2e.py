import io
import json
import os
import tarfile
import tempfile
import unittest
from unittest.mock import call, MagicMock, Mock

from django.conf.urls import url, include
from django.conf import settings
from django.test import TestCase, override_settings

from mc_utils import test_utils
from flow_engines.flow import Flow
from a2g2.a2g2_client.a2g2_client import A2G2_Client
from a2g2.runners.odyssey_push_runner.odyssey_push_runner import (
    OdysseyPushRunner)
from a2g2.a2g2_dj import urls as _a2g2_dj_urls
from jobs import urls as _jobs_urls
from missions import urls as _missions_urls
from storage import urls as _storage_urls
from storage_client.storage_client import MissionControlStorageClient
from a2g2.task_engines.job_task_engine import JobTaskEngine
from job_runners.action_processor import ActionProcessor


class ConfgenFlowGenerator(object):
    flow_type = 'confgen'

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.data['flow_spec'] = flow_spec
        flow.add_task(
            key='confgen_run', 
            as_root=True,
            task={
                'task_engine': JobTaskEngine.__name__,
                'input': {
                    'job_spec': {
                        'job_type': 'confgen',
                        'confgen': {
                            'smiles': flow_spec['input']['smiles'],
                            'params': flow_spec['input']['confgen_params'],
                        },
                        'post_exec_actions': [
                            {
                                'action': 'storage:upload',
                                'params': {
                                    'src': {
                                        'template': '{{ctx.completed_dir}}'
                                    }
                                },
                                'output_to_ctx_target': (
                                    'data.output.raw_dir_storage_params')
                            }
                        ]
                    }
                },
                'status': 'PENDING',
            }
        )
        flow.add_task(
            key='confgen_load',
            precursor_keys=['confgen_run'],
            task={
                'pre_start_actions': [
                    {
                        'action': 'set_ctx_value',
                        'description': 'wire output from job task to job input',
                        'params': {
                            'value': {
                                'template': (
                                    '{{ctx.flow.tasks.confgen_run.output'
                                    '.raw_dir_storage_params}}'
                                ),
                            },
                            'target': ('task.input.job_spec.input'
                                       '.json_storage_params'),
                        }
                    }
                ],
                'task_engine': JobTaskEngine.__name__,
                'input': {
                    'job_spec': {
                        'job_type': 'confgen:load',
                        'pre_build_actions': [
                            {
                                'action': 'storage:download',
                                'params': {
                                    'json_src_params': {
                                        'template': (
                                            '{{ctx.job_spec.input.'
                                            'json_storage_params}}'
                                        ),
                                    },
                                    'dest': {
                                        'template': '{{ctx.job_dir}}/raw_dir'
                                    }
                                },
                                'output_to_ctx_target': 'data.input.raw_dir'
                            }
                        ]
                    }
                },
                'status': 'PENDING'
            }
        )
        return flow

    @classmethod
    def get_dependencies(cls):
        return {
            'task_engines': set([JobTaskEngine()]),
        }

BASE_PATH = 'test_api/'
BASE_URL = '/' + BASE_PATH
STORAGE_PATH = BASE_PATH + 'storage/'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include(_a2g2_dj_urls.__name__)),
    url(r'^%s' % BASE_PATH, include(_jobs_urls.__name__)),
    url(r'^%s' % BASE_PATH, include(_missions_urls.__name__)),
    url(r'^%s' % STORAGE_PATH, include(_storage_urls.__name__)),
]

@override_settings(ROOT_URLCONF=__name__)
class ConfgenFlow_E2E_TestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.setup_storage_basedir()
        test_utils.patch_request_client(
            request_client=self.client,
            methods_to_patch=['get', 'post', 'patch']
        )
        self.a2g2_client = self.generate_a2g2_client()
        self.storage_client = self.generate_storage_client()
        self.action_processor = self.generate_action_processor()
        self.execution_client = self.generate_execution_client()
        self.job_dir_factory = MagicMock()
        self.flow_and_job_runner = self.generate_flow_and_job_runner()
        self.flow_client = self.flow_and_job_runner.flow_client
        self.flow_and_job_runner.tick = Mock(
            side_effect=self.flow_and_job_runner.tick)

    def setup_storage_basedir(self):
        self.storage_base_dir = tempfile.mkdtemp()
        settings_attr = 'STORAGE_FILESYSTEM_BACKEND_BASEDIR'
        setattr(settings, settings_attr, self.storage_base_dir)

    def generate_a2g2_client(self):
        return A2G2_Client(base_url=BASE_URL, request_client=self.client)

    def generate_storage_client(self):
        storage_client = MissionControlStorageClient(
            base_url='/' + STORAGE_PATH,
            request_client=self.client)
        return storage_client

    def generate_action_processor(self):
        action_processor = ActionProcessor()
        action_processor.register_handler(key='storage:upload',
                                          handler=self._upload_wrapper)
        action_processor.register_handler(key='storage:download',
                                          handler=self._download_wrapper)
        return action_processor

    def _upload_wrapper(self, *args, params=None, ctx=None, **kwargs):
        return _upload_action_handler(self.storage_client, *args,
                                      params=params, ctx=ctx)

    def _download_wrapper(self, *args, params=None, ctx=None, **kwargs):
        return _download_action_handler(self.storage_client, *args,
                                        params=params, ctx=ctx)

    def generate_execution_client(self):
        execution_client = MagicMock()
        execution_client.get_execution_state = Mock(return_value={})
        def mock_start_execution(job=None):
            return {'dir': tempfile.mkdtemp()}
        execution_client.start_execution.side_effect = mock_start_execution
        return execution_client

    def generate_flow_and_job_runner(self):
        return OdysseyPushRunner(
            action_processor=self.action_processor,
            flow_generator_classes=[ConfgenFlowGenerator],
            request_client=self.client,
            job_server_url=BASE_URL,
            flow_server_url=BASE_URL,
            job_dir_factory=self.job_dir_factory,
            job_runner_kwargs={'execution_client': self.execution_client,
                               'action_processor': self.action_processor,}
        )

    def test_flow(self):
        self.generate_molecule_library()
        self.create_flows()
        self.assertTrue(len(self.flow_client.fetch_tickable_flows()) > 0)
        self.run_flows_to_completion()
        self.assertTrue(self.flow_and_job_runner.tick_counter > 0)
        self.assert_domain_db_has_expected_state()

    def generate_molecule_library(self):
        initial_smiles = ['smiles_%s' % i for i in range(3)]
        for smiles in initial_smiles:
            self.a2g2_client.create_mol({'props': {'smiles': smiles}})

    def create_flows(self):
        for mol in self.a2g2_client.query(q={'collection': 'mols'})[:1]:
            self.create_confgen_flow(mol=mol)

    def create_confgen_flow(self, mol=None):
        flow_spec = {
            'flow_type': ConfgenFlowGenerator.flow_type,
            'input': {
                'smiles': mol['props']['smiles'],
                'confgen_params': {},
            },
        }
        flow = {'spec': json.dumps(flow_spec)}
        self.flow_client.create_flow(flow=flow)

    def run_flows_to_completion(self):
        max_ticks = 20
        incomplete_flows = self.get_incomplete_flows()
        while len(incomplete_flows) > 0:
            self.flow_and_job_runner.tick()
            if self.flow_and_job_runner.tick_counter > max_ticks:
                raise Exception("Exceeded max_ticks of '%s'" % max_ticks)
            incomplete_flows = self.get_incomplete_flows()

    def get_incomplete_flows(self):
        complete_flows = self.fetch_and_key_flows(
            query_params={'status': 'COMPLETED'})
        all_flows = self.fetch_and_key_flows()
        incomplete_flows = [flow for flow_uuid, flow in all_flows.items()
                            if flow_uuid not in complete_flows]
        return incomplete_flows

    def fetch_and_key_flows(self, query_params=None):
        keyed_flows = {
            flow['uuid']: flow
            for flow in self.flow_client.fetch_flows(query_params=query_params)
        }
        return keyed_flows

    def assert_domain_db_has_expected_state(self):
        expected_counts = {'Mol': 10}
        actual_counts = self.a2g2_client.get_counts()
        self.assertEqual(actual_counts, expected_counts)

def _upload_action_handler(storage_client, *args, params=None, ctx=None,
                           **kwargs):
    dir_path = ctx.transform_value(params['src'])
    tgz_bytes = _dir_to_tgz_bytes(dir_path=dir_path)
    updated_params = storage_client.post_data(data=tgz_bytes)
    serialized_params = json.dumps(updated_params)
    return serialized_params

def _dir_to_tgz_bytes(dir_path=None):
    mem_file = io.BytesIO()
    tgz = tarfile.open(mode='w:gz', fileobj=mem_file)
    tgz.add(dir_path, arcname='.')
    tgz.close()
    return mem_file.getvalue()

class UploadActionHandlerTestCase(unittest.TestCase):
    def setUp(self):
        self.storage_client = Mock()
        self.storage_client.post_data.return_value = {'some': 'object'}
        self.src_dir = self.setup_src_dir()
        self.ctx = self.generate_ctx()
        self.params = {'src': self.src_dir}

    def generate_ctx(self):
        ctx = Mock()
        def transform_value(value): return value
        ctx.transform_value.side_effect = transform_value
        return ctx

    def setup_src_dir(self):
        src_dir = tempfile.mkdtemp()
        for i in range(3):
            path = os.path.join(src_dir, 'file_%s' % i)
            with open(path, 'w') as f: f.write(str(i))
        return src_dir

    def do_upload(self):
        return _upload_action_handler(storage_client=self.storage_client,
                                      params=self.params, ctx=self.ctx)

    def test_posts_tgz_of_transformed_src(self):
        self.do_upload()
        self.assertEqual(self.ctx.transform_value.call_args,
                         call(self.params['src']))
        expected_data = _dir_to_tgz_bytes(
            dir_path=self.ctx.transform_value(self.params['src']))
        self.assertEqual(self.storage_client.post_data.call_args,
                         call(data=expected_data))

    def test_returns_serialized_result_of_post_data(self):
        result = self.do_upload()
        self.assertEqual(result,
                         json.dumps(self.storage_client.post_data.return_value))

def _download_action_handler(storage_client, *args, params=None, ctx=None,
                           **kwargs):
    json_src_params = ctx.transform_value(params['json_src_params'])
    src_params = json.loads(json_src_params)
    transformed_dest = ctx.transform_value(params['dest'])
    tgz_bytes = storage_client.get_data(storage_params=src_params)
    _tgz_bytes_to_dir(tgz_bytes=tgz_bytes, dir_path=transformed_dest)
    return transformed_dest

def _tgz_bytes_to_dir(tgz_bytes=None, dir_path=None):
    mem_file = io.BytesIO(tgz_bytes)
    tgz = tarfile.open(mode='r:gz', fileobj=mem_file)
    tgz.extractall(path=dir_path)

class DownloadActionHandlerTestCase(unittest.TestCase):
    def setUp(self):
        self.storage_client = Mock()
        self.src_params = {'some': 'src params'}
        self.json_src_params = json.dumps(self.src_params)
        self.src_dir = self.setup_src_dir()
        self.src_tgz_bytes = _dir_to_tgz_bytes(dir_path=self.src_dir)
        self.storage_client.get_data.return_value = self.src_tgz_bytes
        self.dest_dir = tempfile.mkdtemp()
        self.dest = os.path.join(self.dest_dir, 'dest')
        self.params = {'json_src_params': self.json_src_params,
                       'dest': self.dest}
        self.ctx = self.generate_ctx()

    def setup_src_dir(self):
        src_dir = tempfile.mkdtemp()
        for i in range(3):
            path = os.path.join(src_dir, 'file_%s' % i)
            with open(path, 'w') as f: f.write(str(i))
        return src_dir

    def generate_ctx(self):
        ctx = Mock()
        def transform_value(value): return value
        ctx.transform_value.side_effect = transform_value
        return ctx

    def do_download(self):
        return _download_action_handler(storage_client=self.storage_client,
                                        params=self.params, ctx=self.ctx)

    def test_calls_get_data_with_deserialized_transformed_src_params(self):
        self.do_download()
        transformed_params = self.ctx.transform_value(
            self.params['json_src_params'])
        deserialized_params = json.loads(transformed_params)
        self.assertEqual(self.storage_client.get_data.call_args,
                         call(storage_params=deserialized_params))

    def test_writes_tgz_to_transformed_dest(self):
        self.do_download()
        self.assert_dirs_equal(self.src_dir,
                               self.ctx.transform_value(self.params['dest']))

    def assert_dirs_equal(self, left, right):
        from filecmp import dircmp
        dcmp = dircmp(left, right)
        self.assertEqual(dcmp.left_only, [])
        self.assertEqual(dcmp.right_only, [])

    def test_returns_transformed_dest(self):
        result = self.do_download()
        self.assertEqual(result, self.ctx.transform_value(self.params['dest']))
