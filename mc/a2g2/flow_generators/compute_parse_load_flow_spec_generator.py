import yaml


class ComputeParseLoadFlowSpecGenerator(object):
    def __init__(self, *args, flow_params=None, **kwargs):
        self.flow_params = flow_params

    def generate_flow_spec(self):
        flow_spec = {
            'label': self.flow_params.get('label'),
            'node_specs': self.generate_node_specs(),
        }
        return flow_spec

    def generate_node_specs(self):
        node_specs = []
        node_specs.append({'node': self.generate_compute_node(),
                           'precursor_keys': ['ROOT']})
        node_specs.append({'node': self.generate_parse_load_node(),
                           'precursor_keys': ['compute']})
        #node_specs.append({'node': self.generate_parse_node(),
                           #'precursor_keys': ['compute']})
        #node_specs.append({'node': self.generate_load_node(),
                           #'precursor_keys': ['parse']})
        return node_specs

    def generate_compute_node(self):
        compute_job_spec = self.flow_params['compute_job_spec']
        node = yaml.load(
            '''
            node_key: compute
            node_tasks:
            - task_key: run_job
              task_params:
                job_spec:
                  job_type: %(job_type)s
                  job_params: %(job_params_yaml)s
              task_type: a2g2.task_handlers.nodes.run_job_task_handler
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
            task_type: mc.task_handlers.set_values_task_handler
            ''')

    def dump_inline_yaml(self, obj=None):
        return yaml.dump(obj, default_style='"', default_flow_style=True)\
                .strip()

    def generate_parse_load_node(self):
        job_spec = self.flow_params['parse_load_job_spec']
        node = yaml.load(
            '''
            node_key: parse_load
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
              task_type: a2g2.task_handlers.nodes.run_job_task_handler
            - %(expose_job_outputs_task_yaml)s
            ''' % {
                'set_input_artifact_task_yaml': self.dump_inline_yaml(
                    self.generate_set_input_artifact_task(src_node='compute')),
                'job_type': job_spec['job_type'],
                'job_params_yaml': self.dump_inline_yaml(
                    job_spec.get('job_params', {})),
                'expose_job_outputs_task_yaml': self.dump_inline_yaml(
                    self.generate_expose_job_outputs_task()),
            }
        )
        return node


    def generate_parse_node(self):
        parse_job_spec = self.flow_params['parse_job_spec']
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
              task_type: a2g2.task_handlers.nodes.run_job_task_handler
            - %(expose_job_outputs_task_yaml)s
            ''' % {
                'set_input_artifact_task_yaml': self.dump_inline_yaml(
                    self.generate_set_input_artifact_task(src_node='compute')),
                'job_type': parse_job_spec['job_type'],
                'job_params_yaml': self.dump_inline_yaml(
                    parse_job_spec.get('job_params', {})
                ),
                'expose_job_outputs_task_yaml': self.dump_inline_yaml(
                    self.generate_expose_job_outputs_task()),
            }
        )
        return node

    def generate_set_input_artifact_task(self, src_node=None):
        return yaml.load(
            '''
            task_key: set_job_input_artifacts
            task_type: mc.task_handlers.set_values_task_handler
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

    def generate_load_node(self):
        load_job_spec = self.flow_params['load_job_spec']
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
              task_type: a2g2.task-handlers.nodes.run_job_task_handler
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

def generate_flow_spec(*args, flow_params=None, **kwargs):
    generator = ComputeParseLoadFlowSpecGenerator(*args,
                                                  flow_params=flow_params,
                                                  **kwargs)
    return generator.generate_flow_spec()
