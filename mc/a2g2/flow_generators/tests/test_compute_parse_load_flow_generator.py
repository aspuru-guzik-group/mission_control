from collections import defaultdict
import unittest
from unittest.mock import patch, MagicMock

from .. import compute_parse_load_flow_generator

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_generator = \
                compute_parse_load_flow_generator.ComputeParseLoadFlowGenerator

class GetDependenciesTestCase(BaseTestCase):
    def test_returns_expected_dependencies(self):
        expected_dependencies = {
            'node_engines': set([self.flow_generator.job_node_engine])
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
            attr_name = 'generate_{job_name}_node'.format(job_name=job_name)
            flow_generator_patches[attr_name] = MagicMock(
                return_value=defaultdict(MagicMock))
        patchers = {'flow_generator': patch.multiple(self.flow_generator,
                                                     **flow_generator_patches)}
        return patchers

    def test_flow_has_expected_nodes(self):
        flow = self.flow_generator.generate_flow()
        for job_name in ['compute', 'parse', 'load']:
            node_name = '{job_name}_node'.format(job_name=job_name)
            self.assertEqual(
                flow.nodes[node_name],
                getattr(self.flow_generator,
                        'generate_%s' % node_name).return_value
            )
        self.assertEqual(flow.root_node_key, 'compute_node')

class GenerateComputeNodeTestCase(BaseTestCase):
    def test_generates_expected_compute_node(self):
        flow = MagicMock()
        expected_node = {
            'node_tasks': [
                {
                    'task_type': 'run_job',
                    'task_params': {
                        'job_spec': {
                            'job_params': flow.data['flow_spec'].get(
                                'compute_job_params', {}),
                            'job_tasks': [
                                {'type': 'job:execute'},
                                {
                                    'type': 'a2g2:storage:upload',
                                    'description': 'upload completed dir',
                                    'params': {
                                        'src': '{{job.completed_dir}}'
                                    }
                                },
                            ],
                        }
                    }
                }
            ]
        }
        self.assertEqual(
            self.flow_generator.generate_compute_node(flow=flow),
            expected_node
        )

class GenerateParseNode(BaseTestCase):
    def test_generates_expected_parse_node(self):
        flow = MagicMock()
        expected_node = {
            'node_tasks': [
                {
                    'task_type': 'mc:set_values',
                    'description': 'storage_meta plumbing',
                    'task_params': {
                        'value_specs': [
                            {
                                'target': 'node.data.storage_meta_json',
                                'src': (
                                    '{{flow.nodes.%s.output.storage_meta}}' % (
                                        'compute_node'
                                    )
                                )
                            }
                        ]
                    },
                },
                {
                    'task_type': 'run_job',
                    'task_params': {
                        'job_spec': {
                            'job_params': flow.data['flow_spec'].get(
                                'parse_job_params'),
                            'job_tasks': [
                                {
                                    'task_type': 'a2g2:storage:download',
                                    'task_params': {
                                        'storage_meta_json': '<UNSET>',
                                        'dest': '{{ctx.job_dir}}/dir_to_parse'
                                    }
                                },
                                {'task_type': 'job:execute'},
                                {
                                    'task_type': 'a2g2:storage:upload',
                                    'task_params': {
                                        'src': '{{job.completed_dir}}'
                                    }
                                },
                            ]
                        }
                    }
                }
            ]
        }
        self.assertEqual(
            self.flow_generator.generate_parse_node(flow=flow),
            expected_node
        )

class GenerateLoadNodeTestCase(BaseTestCase):
    def test_generates_expected_load_node(self):
        flow = MagicMock()
        expected_node = {
            'node_tasks': [
                {
                    'task_type': 'mc:set_values',
                    'description': 'storage_meta plumbing',
                    'task_params': {
                        'value_specs': [
                            {
                                'target': 'node.data.storage_meta_json',
                                'src': (
                                    '{{flow.nodes.%s.output.storage_meta}}' % (
                                        'parse_node'
                                    )
                                )
                            }
                        ]
                    },
                },
                {
                    'task_type': 'run_job',
                    'task_params': {
                        'job_spec': {
                            'job_params': flow.data['flow_spec'].get(
                                'load_job_params'),
                            'job_tasks': [
                                {
                                    'task_type': 'a2g2:storage:download',
                                    'task_params': {
                                        'storage_meta_json': '<UNSET>',
                                        'dest': '{{ctx.job_dir}}/dir_to_load'
                                    }
                                },
                                {'task_type': 'job:execute'},
                            ]
                        }
                    }
                }
            ]
        }
        self.assertEqual(
            self.flow_generator.generate_load_node(flow=flow),
            expected_node
        )
