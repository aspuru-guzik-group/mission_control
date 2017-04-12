class ReaxysFlowSpecGenerator():
    flow_type = 'reaxys'

    @classmethod
    def generate_flow_spec(cls):
        flow_spec = {'nodes': {}}
        flow_spec['nodes']['confgen'] = {
            'as_root': True,
            'node_tasks': [
                cls.generate_compute_parse_load_task(
                    flow_params={
                        'job_prefix': 'a2g2:job:confgen',
                        'computation_params': {},
                    }),
            ]
        }
        flow_spec['nodes']['confgen_demux'] = {
            'precursor_keys': ['confgen'],
            'node_tasks': [cls.generate_confgen_demux_task()]
        }

    @classmethod
    def generate_compute_parse_load_task(cls, flow_params=None):
        return {
            'task_key': 'run_flow',
            'task_params': {
                'flow_spec': {
                    'flow_type': 'a2g2:flow:compute_parse_load',
                    'flow_params': flow_params or {}
                }
            }
        }

    @classmethod
    def generate_confgen_demux_task(cls):
        return {
            'task_type': 'demux',
            'task_params': {
                'node_template': {
                    'node_tasks': [
                        {
                            'task_type': 'create_and_run_flow',
                            'task_params': {
                                'flow_spec': cls.generate_dft_flow_spec(),
                            }
                        }
                    ]
                },
                'template_mods': {
                    'dest': ('node_tpl.node_params.node_tasks.0'
                             '.task_params.flow_params.input'
                             '.molecule'),
                    'value': '{{demux.item}}'
                }
            }
        }

    @classmethod
    def generate_dft_flow_spec(cls):
        dft_flow_spec = {'nodes': {}}
        dft_flow_spec['nodes']['root'] = {'as_root': True}
        dft_flow_spec['nodes']['b3lyp_6_31gs_opt'] = {
            'precursor_keys': ['root'],
            'node_tasks': [
                cls.generate_compute_parse_load_task(
                    flow_params={
                        'job_prefix': 'a2g2:job:qcalc',
                        'computation_params': {
                            'method': 'b3lyp',
                            'basis': '6_31gs',
                        },
                    }
                ),
            ]
        }
        dft_flow_spec['nodes']['b3lyp_6_31gs_opt_t1'] = {
            'precursor_keys': ['root'],
            'node_tasks': [
                cls.generate_compute_parse_load_task(
                    flow_params={
                        'job_prefix': 'a2g2:job:qcalc',
                        'computation_params': {
                            'method': 'b3lyp',
                            'basis': '6_31gs',
                        },
                    }
                ),
            ]
        }
        dft_flow_spec['nodes']['b3lyp_6_31gs_tddft'] = {
            'precursor_keys': ['b3lyp_6_31gs_opt'],
            'node_tasks': [
                cls.generate_compute_parse_load_task(
                    flow_params={
                        'job_prefix': 'a2g2:job:qcalc',
                        'computation_params': {
                            'method': 'b3lyp',
                            'basis': '6_31gs',
                        },
                    }
                )
            ]
        }
        return dft_flow_spec
