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
        node = yaml.load(
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
                'job_params_yaml': cls.dump_inline_yaml(
                    flow.data['flow_spec'].get('compute_job_params', {}))
            }
        )
        return node

    @classmethod
    def dump_inline_yaml(cls, obj=None):
        return yaml.dump(obj, default_style='"', default_flow_style=True)\
                .strip()

    @classmethod
    def generate_parse_node(cls, flow=None):
        node = yaml.load(
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
                'job_params_yaml': cls.dump_inline_yaml(
                    flow.data['flow_spec'].get('parse_job_params', {}))
            }
        )
        return node

    @classmethod
    def generate_load_node(cls, flow=None):
        node = yaml.load(
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
                'job_params_yaml': cls.dump_inline_yaml(
                    flow.data['flow_spec'].get('load_job_params', {}))
            }
        )
        return node
