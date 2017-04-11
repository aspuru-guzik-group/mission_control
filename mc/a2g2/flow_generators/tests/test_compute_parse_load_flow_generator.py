import unittest
import textwrap
import yaml

from .. import compute_parse_load_flow_generator

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow_generator = \
                compute_parse_load_flow_generator.ComputeParseLoadFlowGenerator

    def dump_inline_yaml(self, obj=None):
        return yaml.dump(obj, default_style='"', default_flow_style=True)\
                .strip()

class GenerateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.flow_spec = self.generate_flow_spec()
        self.flow = self.flow_generator.generate_flow(flow_spec=self.flow_spec)

    def generate_flow_spec(self):
        flow_spec = {}
        for job_key in ['compute', 'parse', 'load']:
            flow_spec['%s_job_spec' % job_key] = {
                'job_type': '%s_job_type' % job_key,
                'job_params': {('%s_param' % job_key): '%s_value' % job_key}
            }
        return flow_spec

    def test_generates_expected_flow(self):
        self.expected_flow = self.generate_expected_flow()
        self.assertEqual(self.flow.data, self.expected_flow.data)
        self.assertEqual(self.flow.nodes, self.expected_flow.nodes) 
        self.assertEqual(self.flow.root_node_key,
                         self.expected_flow.root_node_key) 
        self.assertEqual(self.flow.edges, self.expected_flow.edges)

    def generate_expected_flow(self):
        expected_flow = compute_parse_load_flow_generator.Flow()
        expected_flow.data['flow_spec'] = self.flow_spec
        expected_flow.nodes = self.generate_expected_flow_nodes()
        expected_flow.root_node_key = 'compute'
        expected_flow.edges = self.generate_expected_flow_edges()
        return expected_flow

    def generate_expected_flow_nodes(self):
        return yaml.load(textwrap.dedent(
            '''
            compute:
              key: compute
              node_tasks:
              - task_key: run_job
                task_params:
                  job_spec:
                    job_params: %(compute_job_params_yaml)s
                    job_type: %(compute_job_type)s
                task_type: a2g2.tasks.nodes.run_job
              - task_key: expose_job_outputs
                task_params:
                  value_specs:
                  - dest: ctx.node.data.serialized_artifact
                    value: '{{ctx.tasks.run_job.data.artifact|tojson}}'
                task_type: a2g2.tasks.set_values
              status: PENDING
            parse:
              key: parse
              node_tasks:
              - task_key: set_job_input_artifacts
                task_params:
                  value_specs:
                  - dest: ctx.tasks.run_job.task_params.job_spec.inputs.serialized_artifacts.input_dir
                    value: '{{ctx.flow.nodes.compute.data.serialized_artifact}}'
                  - dest: "ctx.tasks.run_job.task_params.job_spec.job_params.artifact"
                    value: '{{ctx.flow.nodes.compute.data.serialized_artifact}}'
                task_type: a2g2.tasks.set_values
              - task_key: run_job
                task_params:
                  job_spec:
                    inputs:
                      serialized_artifacts:
                        input_dir: WILL_BE_SET_FROM_PRECURSOR_TASK
                    job_params: %(parse_job_params_yaml)s
                    job_type: %(parse_job_type)s
                task_type: a2g2.tasks.nodes.run_job
              - task_key: expose_job_outputs
                task_params:
                  value_specs:
                  - dest: ctx.node.data.serialized_artifact
                    value: '{{ctx.tasks.run_job.data.artifact|tojson}}'
                task_type: a2g2.tasks.set_values
              status: PENDING
            load:
              key: load
              node_tasks:
              - task_key: set_job_input_artifacts
                task_params:
                  value_specs:
                  - dest: ctx.tasks.run_job.task_params.job_spec.inputs.serialized_artifacts.input_dir
                    value: '{{ctx.flow.nodes.parse.data.serialized_artifact}}'
                  - dest: ctx.tasks.run_job.task_params.job_spec.job_params.artifact
                    value: '{{ctx.flow.nodes.parse.data.serialized_artifact}}'
                task_type: a2g2.tasks.set_values
              - task_key: run_job
                task_params:
                  job_spec:
                    job_params: %(load_job_params_yaml)s
                    job_type: %(load_job_type)s
                task_type: a2g2.tasks.nodes.run_job
              - task_key: expose_job_outputs
                task_params:
                  value_specs:
                  - dest: ctx.node.data.serialized_artifact
                    value: '{{ctx.tasks.run_job.data.artifact|tojson}}'
                task_type: a2g2.tasks.set_values
              status: PENDING
            ''' % {
                'compute_job_type': self.flow_spec.get(
                    'compute_job_spec')['job_type'],
                'compute_job_params_yaml': self.dump_inline_yaml(
                    self.flow_spec.get('compute_job_spec')['job_params']),
                'parse_job_type': self.flow_spec.get(
                    'parse_job_spec')['job_type'],
                'parse_job_params_yaml': self.dump_inline_yaml({
                    **(self.flow_spec['parse_job_spec']['job_params']),
                    'dir_to_parse': 'inputs/dir_to_parse'
                }),
                'load_job_type': self.flow_spec.get(
                    'load_job_spec')['job_type'],
                'load_job_params_yaml': self.dump_inline_yaml(
                    self.flow_spec.get('load_job_spec')['job_params'])
            }
        ))

    def generate_expected_flow_edges(self):
        return {
            ('compute', 'parse'): {'src_key': 'compute', 'dest_key': 'parse'},
            ('parse', 'load'): {'src_key': 'parse', 'dest_key': 'load'}
        }
