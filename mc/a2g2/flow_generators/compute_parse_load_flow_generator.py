from mc.flow_engines.flow import Flow
from mc.a2g2.task_engines.job_task_engine import JobTaskEngine

from . import base_flow_generator


class ComputeParseLoadFlowGenerator(base_flow_generator.BaseFlowGenerator):
    flow_type = 'ComputeParseLoadFlow'
    job_task_engine = JobTaskEngine()
    job_task_engine_name = job_task_engine.__class__.__name__
    task_keys = {
        task_name: task_name for task_name in [
            'compute_job_task',
            'parse_job_task',
            'load_job_task'
        ]
    }

    @classmethod
    def get_dependencies(cls):
        return {
            'task_engines': set([cls.job_task_engine]),
        }

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.data['flow_spec'] = flow_spec
        flow.add_task(
            as_root=True,
            key=cls.task_keys['compute_job_task'],
            task=cls.generate_compute_job_task(flow=flow)
        )
        flow.add_task(
            key=cls.task_keys['parse_job_task'],
            precursor_keys=[cls.task_keys['compute_job_task']],
            task=cls.generate_parse_job_task(flow=flow)
        )
        flow.add_task(
            key=cls.task_keys['load_job_task'],
            precursor_keys=[cls.task_keys['parse_job_task']],
            task=cls.generate_load_job_task(flow=flow)
        )
        return flow

    @classmethod
    def generate_compute_job_task(cls, flow=None):
        job_spec = flow.data['flow_spec']['compute_job_spec']
        task = {
            'task_engine': cls.job_task_engine_name,
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
        return task

    @classmethod
    def generate_parse_job_task(cls, flow=None):
        job_spec = flow.data['flow_spec'].get('parse_job_spec', {})
        task = {
            'task_engine': cls.job_task_engine_name,
            'pre_start_actions': [
                {
                    'action': 'set_ctx_value',
                    'description': (
                        'wire output from compute job task to job input for'
                        ' this task'
                    ),
                    'params': {
                        'value': {
                            'template': (
                                '{{ctx.flow.tasks.{task_key}.output'
                                '.storage_meta}}'
                            ).format(task_key=cls.task_keys['compute_job_task'])
                        },
                        'target': 'task.input.job_spec.input.storage_meta',
                    }
                }
            ],
            'input': {
                'job_spec': {
                    **job_spec,
                    'pre_build_actions': job_spec.get(
                        'pre_build_actions', []) + [
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
                        ]
                }
            }
        }
        return task

    @classmethod
    def generate_load_job_task(cls, flow=None):
        job_spec = flow.data['flow_spec'].get('load_job_spec', {})
        task = {
            'task_engine': cls.job_task_engine_name,
            'pre_start_actions': [
                {
                    'action': 'set_ctx_value',
                    'description': (
                        'wire output from parse job task to job input for'
                        ' this task'
                    ),
                    'params': {
                        'value': {
                            'template': (
                                '{{ctx.flow.tasks.{task_key}.output'
                                '.storage_meta}}'
                            ).format(task_key=cls.task_keys['parse_job_task']),
                        },
                        'target': 'task.input.job_spec.input.storage_meta',
                    }
                }
            ],
            'input': {
                'job_spec': {
                    **job_spec,
                    'pre_build_actions': job_spec.get(
                        'pre_build_actions', []) + [
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
                        ]
                }
            },
        }
        return task
