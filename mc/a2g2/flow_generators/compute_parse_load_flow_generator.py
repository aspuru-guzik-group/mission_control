from mc.flow_engines.flow import Flow
from mc.a2g2.node_engines.job_node_engine import JobNodeEngine

from . import base_flow_generator


class ComputeParseLoadFlowGenerator(base_flow_generator.BaseFlowGenerator):
    flow_type = 'ComputeParseLoadFlow'
    job_node_engine = JobNodeEngine()
    job_node_engine_name = job_node_engine.__class__.__name__

    @classmethod
    def get_dependencies(cls):
        return {
            'node_engines': set([cls.job_node_engine]),
        }

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.data['flow_spec'] = flow_spec
        flow.add_node(
            as_root=True,
            key='compute_node',
            node=cls.generate_compute_node(flow=flow)
        )
        flow.add_node(
            key='parse_node',
            precursor_keys=['compute_node'],
            node=cls.generate_parse_node(flow=flow)
        )
        flow.add_node(
            key='load_node',
            precursor_keys=['parse_node'],
            node=cls.generate_load_node(flow=flow)
        )
        return flow

    @classmethod
    def generate_compute_node(cls, flow=None):
        node = {
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
        return node

    @classmethod
    def generate_parse_node(cls, flow=None):
        node = {
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
        return node

    @classmethod
    def generate_load_node(cls, flow=None):
        node = {
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
                                {'task_type': 'job:execute'}
                            ]
                        }
                    }
                }
            ]
        }
        return node
