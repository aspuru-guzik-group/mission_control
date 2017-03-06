from collections import defaultdict
import unittest
from unittest.mock import patch, MagicMock

from mc.a2g2.task_engines.job_task_engine import JobTaskEngine

from .. import compute_parse_load_flow_generator

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_generator = \
                compute_parse_load_flow_generator.ComputeParseLoadFlowGenerator

class GetDependenciesTestCase(BaseTestCase):
    def test_returns_expected_dependencies(self):
        expected_dependencies = {
            'task_engines': set([self.flow_generator.job_task_engine])
        }
        self.assertEqual(self.flow_generator.get_dependencies(),
                         expected_dependencies)

class GenerateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = self.setup_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def setup_patchers(self):
        flow_generator_patches = {}
        for job_name in ['compute', 'parse', 'load']:
            attr_name = 'generate_{job_name}_job_task'.format(job_name=job_name)
            flow_generator_patches[attr_name] = MagicMock(
                return_value=defaultdict(MagicMock))
        patchers = {'flow_generator': patch.multiple(self.flow_generator,
                                                     **flow_generator_patches)}
        return patchers

    def test_flow_has_expected_tasks(self):
        flow = self.flow_generator.generate_flow()
        for job_name in ['compute', 'parse', 'load']:
            task_name = '{job_name}_job_task'.format(job_name=job_name)
            self.assertEqual(
                flow.tasks[task_name],
                getattr(self.flow_generator,
                        'generate_%s' % task_name).return_value
            )
        self.assertEqual(flow.root_task_key, 'compute_job_task')

class GenerateComputeJobTaskTestCase(BaseTestCase):
    def test_generates_expected_compute_task(self):
        flow = MagicMock()
        job_spec = flow.data['flow_spec']['compute_job_spec']
        expected_task = {
            'task_engine': JobTaskEngine.__name__,
            'input': {
                'job_spec': {
                    **job_spec,
                    'post_exec_actions': job_spec.get(
                        'post_exec_actions', []).append({
                            'description': ('upload completed computation'
                                            ' dir to storage'),
                            'action': 'storage:upload',
                            'params': {
                                'src': {'template': '{{ctx.completed_dir}}'}
                            },
                            'output_to_ctx_target': 'data.output.storage_meta'
                        })
                }
            }
        }
        self.assertEqual(
            self.flow_generator.generate_compute_job_task(flow=flow),
            expected_task
        )

class GenerateParseJobTask(BaseTestCase):
    def test_generates_expected_parse_job_task(self):
        flow = MagicMock()
        job_spec = flow.data['flow_spec']['parse_job_spec']
        expected_task = {
            'task_engine': JobTaskEngine.__name__,
            'pre_start_actions': [
                {
                    'action': 'set_ctx_value',
                    'description': 'wire output from job task to job input',
                    'params': {
                        'value': {
                            'template': (
                                '{{ctx.flow.tasks.confgen_run.output'
                                '.storage_meta}}'
                            ),
                        },
                        'target': 'task.input.job_spec.input.storage_meta',
                    }
                }
            ],
            'input': {
                'job_spec': {
                    **job_spec,
                    'pre_build_actions': job_spec.get(
                        'pre_build_actions', []).extend([
                            {
                                'description': 'download dir to parse',
                                'action': 'storage:download',
                                'params': {
                                    'storage_meta': {
                                        'template': (
                                            '{{ctx.job_spec.input.storage_meta}}'
                                        ),
                                    },
                                    'dest': {
                                        'template': (
                                            '{{ctx.job_dir}}/dir_to_parse')
                                    }
                                },
                            },
                            {
                                'description': 'set name of dir_to_parse',
                                'action': 'set_ctx_value',
                                'params': {
                                    'value': 'dir_to_parse',
                                    'target': 'data.input.dir_to_parse',
                                }
                            },
                        ])
                }
            }
        }
        self.assertEqual(
            self.flow_generator.generate_parse_job_task(flow=flow),
            expected_task
        )

class GenerateLoadJobTaskTestCase(BaseTestCase):
    def test_generates_expected_load_job_task(self):
        flow = MagicMock()
        job_spec = flow.data['flow_spec']['load_job_spec']
        expected_task = {
            'task_engine': JobTaskEngine.__name__,
            'pre_start_actions': [
                {
                    'action': 'set_ctx_value',
                    'description': 'wire output from job task to job input',
                    'params': {
                        'value': {
                            'template': (
                                '{{ctx.flow.tasks.confgen_run.output'
                                '.storage_meta}}'
                            ),
                        },
                        'target': 'task.input.job_spec.input.storage_meta',
                    }
                }
            ],
            'input': {
                'job_spec': {
                    **job_spec,
                    'pre_build_actions': job_spec.get(
                        'pre_build_actions', []).extend([
                            {
                                'action': 'storage:download',
                                'params': {
                                    'storage_meta': {
                                        'template': (
                                            '{{ctx.job_spec.input.storage_meta}}'
                                        ),
                                    },
                                    'dest': {
                                        'template': (
                                            '{{ctx.job_dir}}/dir_to_parse')
                                    }
                                },
                            },
                            {
                                'action': 'set_ctx_value',
                                'description': 'set name of dir_to_parse',
                                'params': {
                                    'value': 'dir_to_parse',
                                    'target': 'data.input.dir_to_parse',
                                }
                            },
                        ])
                }
            },
        }
        self.assertEqual(
            self.flow_generator.generate_load_job_task(flow=flow),
            expected_task
        )
