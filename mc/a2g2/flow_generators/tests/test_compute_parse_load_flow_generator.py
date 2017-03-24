from collections import defaultdict
import unittest
from unittest.mock import patch, MagicMock
import yaml

from .. import compute_parse_load_flow_generator

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_generator = \
                compute_parse_load_flow_generator.ComputeParseLoadFlowGenerator

    def generate_flow_w_mock_job_params(self, params_key=None):
        flow = MagicMock()
        flow.data = {'flow_spec': {params_key: 'mock_job_params'}}
        return flow

    def dump_inline_yaml(self, obj=None):
        return yaml.dump(obj, default_style='"', default_flow_style=True)\
                .strip()

class GenerateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.expected_node_keys = ['compute', 'parse', 'load']
        self.patchers = self.setup_patchers()
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def setup_patchers(self):
        flow_generator_patches = {}
        for node_key in self.expected_node_keys:
            attr = 'generate_{node_key}_node'.format(node_key=node_key)
            flow_generator_patches[attr] = \
                    MagicMock(return_value=defaultdict(MagicMock))
        patchers = {'flow_generator': patch.multiple(self.flow_generator,
                                                     **flow_generator_patches)}
        return patchers

    def test_flow_has_expected_nodes(self):
        flow = self.flow_generator.generate_flow()
        for node_key in self.expected_node_keys:
            self.assertEqual(
                flow.nodes[node_key],
                getattr(self.flow_generator,
                        'generate_%s_node' % node_key).return_value
            )
        self.assertEqual(flow.root_node_key, 'compute')

class GenerateComputeNodeTestCase(BaseTestCase):
    def test_generates_expected_compute_node(self):
        flow = self.generate_flow_w_mock_job_params(
            params_key='compute_job_params')
        expected_node = yaml.load(
            '''
            node_tasks:
            - task_key: run_job
              task_params:
                job_spec:
                  job_params: %(job_params_yaml)s
                  job_tasks:
                  - task_key: execute
                    task_type: a2g2.tasks.jobs.execute
                  - task_key: put_artifact
                    task_params:
                      src: '{{ctx.tasks.execute.data.artifact}}'
                    task_type: a2g2.tasks.jobs.storage.put
                  - task_key: expose_artifact
                    task_params:
                      value_specs:
                      - dest: ctx.job.data.outputs.artifact
                        value: '{{ctx.tasks.put_artifact.data.artifact}}'
                    task_type: a2g2.tasks.set_values
              task_type: a2g2.tasks.nodes.run_job
            - task_key: expose_job_outputs
              task_params:
                value_specs:
                - dest: ctx.node.data.outputs
                  value: '{{ctx.tasks.run_job.data.outputs}}'
              task_type: a2g2.tasks.set_values
            ''' % {
                'job_params_yaml': self.dump_inline_yaml(
                    obj=flow.data['flow_spec'].get('compute_job_params', {}))
            }
        )
        self.assertEqual(
            self.flow_generator.generate_compute_node(flow=flow),
            expected_node
        )

class GenerateParseNode(BaseTestCase):
    def test_generates_expected_parse_node(self):
        flow = self.generate_flow_w_mock_job_params(
            params_key='parse_job_params')
        expected_node = yaml.load(
            '''
            node_tasks:
            - task_key: set_job_input_artifact
              task_type: a2g2.tasks.set_values
              task_params:
                value_specs:
                  - dest: ctx.node.tasks.run_job.task_params.job_spec
                      .job_tasks.0.task_params.artifact
                    value: '{{ctx.flow.nodes.compute.data.outputs.artifact}}'
            - task_key: run_job
              task_params:
                job_spec:
                  job_params: %(job_params_yaml)s
                  job_tasks:
                  - task_key: get_input_artifact
                    task_type: a2g2.tasks.jobs.storage.get
                    task_params:
                      artifact: TO_BE_WIRED_IN
                      dest: '{{ctx.job.dir}}/dir_to_parse'
                  - task_key: execute
                    task_type: a2g2.tasks.jobs.execute
                  - task_key: put_artifact
                    task_params:
                      src: '{{ctx.tasks.execute.data.artifact}}'
                    task_type: a2g2.tasks.jobs.storage.put
                  - task_key: expose_artifact
                    task_params:
                      value_specs:
                      - dest: ctx.job.data.outputs.artifact
                        value: '{{ctx.tasks.put_artifact.data.artifact}}'
                    task_type: a2g2.tasks.set_values
              task_type: a2g2.tasks.nodes.run_job
            - task_key: expose_job_outputs
              task_params:
                value_specs:
                - dest: ctx.node.data.outputs
                  value: '{{ctx.tasks.run_job.data.outputs}}'
              task_type: a2g2.tasks.set_values
            ''' % {
                'job_params_yaml': self.dump_inline_yaml(
                    obj=flow.data['flow_spec'].get('parse_job_params', {}))
            }
        )
        self.assertEqual(
            self.flow_generator.generate_parse_node(flow=flow),
            expected_node
        )

class GenerateLoadNodeTestCase(BaseTestCase):
    def test_generates_expected_load_node(self):
        flow = self.generate_flow_w_mock_job_params(
            params_key='load_job_params')
        expected_node = yaml.load(
            '''
            node_tasks:
            - task_key: set_job_input_artifact
              task_type: a2g2.tasks.set_values
              task_params:
                value_specs:
                  - dest: ctx.node.tasks.run_job.task_params.job_spec
                      .job_tasks.0.task_params.artifact
                    value: '{{ctx.flow.nodes.parse.data.outputs.artifact}}'
            - task_key: run_job
              task_params:
                job_spec:
                  job_params: %(job_params_yaml)s
                  job_tasks:
                  - task_key: get_input_artifact
                    task_type: a2g2.tasks.jobs.storage.get
                    task_params:
                      artifact: TO_BE_WIRED_IN
                      dest: '{{ctx.job.dir}}/dir_to_load'
                  - task_key: execute
                    task_type: a2g2.tasks.jobs.execute
            ''' % {
                'job_params_yaml': self.dump_inline_yaml(
                    obj=flow.data['flow_spec'].get('load_job_params', {}))
            }
        )
        self.assertEqual(
            self.flow_generator.generate_load_node(flow=flow),
            expected_node
        )
