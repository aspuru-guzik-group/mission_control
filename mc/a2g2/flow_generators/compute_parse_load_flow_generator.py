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
            - %(expose_job_outputs_task_yaml)s
            ''' % {
                'job_type': compute_job_spec['job_type'],
                'job_params_yaml': cls.dump_inline_yaml(
                    compute_job_spec.get('job_params', {})),
                'expose_job_outputs_task_yaml': cls.dump_inline_yaml(
                    cls.generate_expose_job_outputs_task()),
            }
        )
        return node

    @classmethod
    def generate_expose_job_outputs_task(cls):
        return yaml.load(
            '''
            task_key: expose_job_outputs
            task_params:
              value_specs:
              - dest: ctx.node.data.serialized_artifact
                value: '{{ctx.tasks.run_job.data.artifact|tojson}}'
            task_type: a2g2.tasks.set_values
            ''')

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
            - %(set_input_artifact_task_yaml)s
            - task_key: run_job
              task_params:
                job_spec:
                  inputs:
                    serialized_artifacts:
                      input_dir: WILL_BE_SET_FROM_PRECURSOR_TASK
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.tasks.nodes.run_job
            - %(expose_job_outputs_task_yaml)s
            ''' % {
                'set_input_artifact_task_yaml': cls.dump_inline_yaml(
                    cls.generate_set_input_artifact_task(src_node='compute')),
                'job_type': parse_job_spec['job_type'],
                'job_params_yaml': cls.dump_inline_yaml({
                    **(parse_job_spec.get('job_params', {})),
                    'dir_to_parse': 'inputs/dir_to_parse'
                }),
                'expose_job_outputs_task_yaml': cls.dump_inline_yaml(
                    cls.generate_expose_job_outputs_task()),
            }
        )
        return node

    @classmethod
    def generate_set_input_artifact_task(cls, src_node=None):
        return yaml.load(
            '''
            task_key: set_job_input_artifacts
            task_type: a2g2.tasks.set_values
            task_params:
              value_specs:
                - dest: "ctx.tasks.run_job.task_params.job_spec\\
                    .inputs.serialized_artifacts.input_dir"
                  value: "{{ctx.flow.nodes.%(src_node)s.data\\
                    .serialized_artifact}}"
            ''' % {'src_node': src_node}
        ) 

    @classmethod
    def generate_load_node(cls, flow=None):
        load_job_spec = flow.data['flow_spec']['load_job_spec']
        node = yaml.load(
            '''
            node_tasks:
            - %(set_input_artifact_task_yaml)s
            - task_key: run_job
              task_params:
                job_spec:
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.tasks.nodes.run_job
            - %(expose_job_outputs_task_yaml)s
            ''' % {
                'set_input_artifact_task_yaml': cls.dump_inline_yaml(
                    cls.generate_set_input_artifact_task(src_node='parse')),
                'job_type': load_job_spec['job_type'],
                'job_params_yaml': cls.dump_inline_yaml(
                    load_job_spec.get('job_params', {})),
                'expose_job_outputs_task_yaml': cls.dump_inline_yaml(
                    cls.generate_expose_job_outputs_task()),
            }
        )
        return node
