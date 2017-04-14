import yaml
from mc.task_handlers.base_task_handler import BaseTaskHandler


class RunComputeParseLoadFlowTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        flow_task = self.generate_flow_task(task=task,
                                            task_context=task_context)
        task_context['node']['node_tasks'].append(flow_task)
        task['status'] = 'COMPLETED'

    def generate_flow_task(self, task=None, task_context=None):
        label = task['task_key'] + '__flow'
        flow_task = {
            'task_key': label,
            'task_type': 'a2g2.tasks.nodes.run_flow',
            'task_params': {
                'flow_spec': self.generate_flow_spec(task=task,
                                                     task_context=task_context,
                                                     label=label)
            }
        }
        return flow_task

    def generate_flow_spec(self, task=None, task_context=None, label=None):
        flow_spec = {}
        flow_spec['label'] = label
        flow_params = task['task_params'].get('flow_params', {})
        flow_spec['data'] = flow_params.get('data', {})
        flow_spec['cfg'] = flow_params.get('cfg', {})
        flow_spec['node_specs'] = self.generate_node_specs(
            flow_params=flow_params)
        return flow_spec

    def generate_node_specs(self, task=None, task_context=None,
                            flow_params=None):
        node_specs = []
        node_specs.append({
            'node': self.generate_compute_node(flow_params=flow_params),
            'precursor_keys': ['ROOT']
        })
        node_specs.append({
            'node': self.generate_parse_node(flow_params=flow_params),
            'precursor_keys': ['compute']
        })
        node_specs.append({
            'node': self.generate_load_node(flow_params=flow_params),
            'precursor_keys': ['parse']
        })
        return node_specs

    def generate_compute_node(self, flow_params=None):
        compute_job_spec = flow_params['compute_job_spec']
        node = yaml.load(
            '''
            node_key: compute
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
                'job_params_yaml': self.dump_inline_yaml(
                    compute_job_spec.get('job_params', {})),
                'expose_job_outputs_task_yaml': self.dump_inline_yaml(
                    self.generate_expose_job_outputs_task()),
            }
        )
        return node

    def generate_expose_job_outputs_task(self):
        return yaml.load(
            '''
            task_key: expose_job_outputs
            task_params:
              value_specs:
              - dest: ctx.node.data.artifact
                source: ctx.tasks.run_job.data.artifact
            task_type: a2g2.tasks.set_values
            ''')

    def dump_inline_yaml(self, obj=None):
        return yaml.dump(obj, default_style='"', default_flow_style=True)\
                .strip()

    def generate_parse_node(self, flow_params=None):
        parse_job_spec = flow_params['parse_job_spec']
        node = yaml.load(
            '''
            node_key: parse
            node_tasks:
            - %(set_input_artifact_task_yaml)s
            - task_key: run_job
              task_params:
                job_spec:
                  inputs:
                    artifacts:
                      input_dir: WILL_BE_SET_FROM_PRECURSOR_TASK
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.tasks.nodes.run_job
            - %(expose_job_outputs_task_yaml)s
            ''' % {
                'set_input_artifact_task_yaml': self.dump_inline_yaml(
                    self.generate_set_input_artifact_task(src_node='compute')),
                'job_type': parse_job_spec['job_type'],
                'job_params_yaml': self.dump_inline_yaml({
                    **(parse_job_spec.get('job_params', {})),
                    'dir_to_parse': 'inputs/dir_to_parse'
                }),
                'expose_job_outputs_task_yaml': self.dump_inline_yaml(
                    self.generate_expose_job_outputs_task()),
            }
        )
        return node

    def generate_set_input_artifact_task(self, src_node=None):
        return yaml.load(
            '''
            task_key: set_job_input_artifacts
            task_type: a2g2.tasks.set_values
            task_params:
              value_specs:
                - dest: "ctx.tasks.run_job.task_params.job_spec\\
                    .inputs.artifacts.input_dir"
                  source: "ctx.flow.nodes.%(src_node)s.data.artifact"
                - dest: "ctx.tasks.run_job.task_params.job_spec\\
                    .job_params.artifact"
                  source: "ctx.flow.nodes.%(src_node)s.data.artifact"
            ''' % {'src_node': src_node}
        ) 

    def generate_load_node(self, flow_params=None):
        load_job_spec = flow_params['load_job_spec']
        node = yaml.load(
            '''
            node_key: load
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
                'set_input_artifact_task_yaml': self.dump_inline_yaml(
                    self.generate_set_input_artifact_task(src_node='parse')),
                'job_type': load_job_spec['job_type'],
                'job_params_yaml': self.dump_inline_yaml(
                    load_job_spec.get('job_params', {})),
                'expose_job_outputs_task_yaml': self.dump_inline_yaml(
                    self.generate_expose_job_outputs_task()),
            }
        )
        return node
