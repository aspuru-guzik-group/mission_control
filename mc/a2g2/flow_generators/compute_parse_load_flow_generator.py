import yaml

from mc.flow_engines.flow import Flow

from . import base_flow_generator


class ComputeParseLoadFlowGenerator(base_flow_generator.BaseFlowGenerator):
    flow_type = 'ComputeParseLoadFlow'

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.data['flow_spec'] = flow_spec
        flow.add_node(
            as_root=True,
            key='compute',
            node=cls.generate_compute_node(flow=flow)
        )
        flow.add_node(
            key='parse',
            precursor_keys=['compute'],
            node=cls.generate_parse_node(flow=flow)
        )
        flow.add_node(
            key='load',
            precursor_keys=['parse'],
            node=cls.generate_load_node(flow=flow)
        )
        return flow

    @classmethod
    def generate_compute_node(cls, flow=None):
        compute_job_spec = flow.data['flow_spec']['compute_job_spec']
        node = yaml.load(
            '''
            node_tasks:
            - task_key: run_job
              task_params:
                job_spec:
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.tasks.nodes.run_job
            - task_key: expose_job_outputs
              task_params:
                value_specs:
                - dest: ctx.node.data.artifact
                  value: '{{ctx.tasks.run_job.data.artifact}}'
              task_type: a2g2.tasks.set_values
            ''' % {
                'job_type': compute_job_spec['job_type'],
                'job_params_yaml': cls.dump_inline_yaml(
                    compute_job_spec.get('compute_job_params', {}))
            }
        )
        return node

    @classmethod
    def dump_inline_yaml(cls, obj=None):
        return yaml.dump(obj, default_style='"', default_flow_style=True)\
                .strip()

    @classmethod
    def generate_parse_node(cls, flow=None):
        parse_job_spec = flow.data['flow_spec']['parse_job_spec']
        node = yaml.load(
            '''
            node_tasks:
            - task_key: set_job_input_artifacts
              task_type: a2g2.tasks.set_values
              task_params:
                value_specs:
                  - dest: ctx.node.tasks.run_job.task_params.job_spec
                      .inputs.artifacts.dir_to_parse
                    value: '{{ctx.flow.nodes.compute.data.artifact}}'
            - task_key: run_job
              task_params:
                job_spec:
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.tasks.nodes.run_job
            - task_key: expose_job_outputs
              task_params:
                value_specs:
                - dest: ctx.node.data.outputs
                  value: '{{ctx.tasks.run_job.data.outputs}}'
              task_type: a2g2.tasks.set_values
            ''' % {
                'job_type': parse_job_spec['job_type'],
                'job_params_yaml': cls.dump_inline_yaml({
                    **(parse_job_spec.get('job_params', {})),
                    'dir_to_parse': 'inputs/dir_to_parse'
                })
            }
        )
        return node

    @classmethod
    def generate_load_node(cls, flow=None):
        load_job_spec = flow.data['flow_spec']['load_job_spec']
        node = yaml.load(
            '''
            node_tasks:
            - task_key: set_job_input_artifact
              task_type: a2g2.tasks.set_values
              task_params:
                value_specs:
                  - dest: ctx.node.tasks.run_job.task_params.job_spec
                      .inputs.artifacts
                    value: '{{ctx.flow.nodes.parse.data.outputs.artifact}}'
            - task_key: run_job
              task_params:
                job_spec:
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.tasks.nodes.run_job
            - task_key: expose_job_outputs
              task_params:
                value_specs:
                - dest: ctx.node.data.outputs
                  value: '{{ctx.tasks.run_job.data.outputs}}'
              task_type: a2g2.tasks.set_values
            ''' % {
                'job_type': load_job_spec['job_type'],
                'job_params_yaml': cls.dump_inline_yaml(
                    load_job_spec.get('job_params', {}))
            }
        )
        return node
