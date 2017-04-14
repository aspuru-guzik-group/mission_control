from . import compute_parse_load_flow_spec_generator

class ReaxysFlowSpecGenerator(object):
    def __init__(self, *args, flow_params=None, **kwargs):
        self.flow_params = flow_params

    def generate_flow_spec(self):
        node_specs = []
        #node_specs.append({
            #'node': {
                #'node_key': 'confgen',
                #'data': {
                    #'conformer_chemthings': ['FAKE_%s' % i for i in range(3)]
                #}
            #},
            #'precursor_keys': ['ROOT']
        #})
        #node_specs.append({
            #'node': {
                #'node_key': 'confgen',
                #'node_tasks': [
                    #self.generate_compute_parse_load_task(
                        #flow_params={
                            #'compute_job_spec': {
                                #'job_type': 'a2g2.jobs.confgen.confgen',
                            #},
                            #'parse_job_spec': {
                                #'job_type': 'a2g2.jobs.confgen.parse',
                            #},
                            #'load_job_spec': {
                                #'job_type': 'a2g2.jobs.confgen.load',
                            #},
                        #}),
                #],
                #'data': { 
                    ## @TODO! REPLACE!
                    ##'conformer_chemthings': ['FAKE_%s' % i for i in range(3)]
                    #'conformer_chemthings': [1],
                #}
            #},
            #'precursor_keys': ['ROOT'],
        #})
        node_specs.append({
            'node': {
                'node_key': 'conformer_query',
                'node_tasks': [
                    {
                        'task_key': 'run_query_job',
                        'task_params': {
                            'job_spec': {
                                'job_type': 'a2g2.jobs.query',
                                'job_params': {},
                            }
                        },
                        'task_type': 'a2g2.tasks.nodes.run_job'
                    }
                ]
            },
            'precursor_keys': ['ROOT']
        })
        node_specs.append({
            'node': {'node_key': 'confgen_demux',
                     'node_tasks': self.generate_confgen_demux_tasks()},
            'precursor_keys': ['conformer_query']
        })
        flow_spec = {'label': 'reaxys_flow', 'node_specs': node_specs}
        return flow_spec

    def generate_compute_parse_load_task(self, flow_params=None):
        return {
            'task_key': 'run_compute_parse_load_flow',
            'task_type': 'a2g2.tasks.nodes.run_flow',
            'task_params': {
                'flow_spec': (
                    compute_parse_load_flow_spec_generator.generate_flow_spec(
                        flow_params=flow_params)
                ),
                'delete_flow_spec_after_creation': True,
            }
        }

    def generate_confgen_demux_tasks(self):
        confgen_demux_tasks = []
        wiring_task = {
            'task_key': 'wire_confgen_demux',
            'task_type': 'a2g2.tasks.set_values',
            'task_params': {
                'value_specs': [
                    {
                        'dest': 'ctx.tasks.confgen_demux.task_params.items',
                        'value': ('ctx.flow.nodes.confgen.data'
                                  '.conformer_chemthings')
                    }
                ]
            }
        }
        confgen_demux_tasks.append(wiring_task)
        demux_task = {
            'task_key': 'confgen_demux',
            'task_type': 'a2g2.tasks.nodes.demux',
            'task_params': {
                'items': 'TO BE WIRED',
                'demux_flow_params': {
                    'cfg': {
                        'fail_fast': False,
                    },
                    'data': {
                        'dft_flow_spec': self.generate_dft_flow_spec()
                    }
                },
                'node_spec_template_obj': {
                    'node': {
                        'node_tasks': [
                            {
                                'task_key': 'run_dft_flow',
                                'task_type': 'a2g2.tasks.nodes.run_flow',
                                'task_params': {
                                    '_task_params_type': 'interpolated_object',
                                    'object_template': {
                                        'delete_flow_spec_after_creation': True,
                                    },
                                    'interpolations': [
                                        {
                                            'value': ('ctx.task_context.flow'
                                                      '.data.dft_flow_spec'),
                                            'dest': 'ctx.object.flow_spec'
                                        }
                                    ]
                                }
                            }
                        ]
                    }
                },
                'substitutions': [
                    {
                        'dest': 'ctx.node_spec.node.data.conformer_chemthing',
                        'value': 'ctx.demux_item'
                    }
                ]
            }
        }
        confgen_demux_tasks.append(demux_task)
        return confgen_demux_tasks

    def generate_dft_flow_spec(self):
        node_specs = []
        node_specs.append({'node': {'node_key': 'root'}, 'as_root': True})
        node_specs.append({
            'node': {
                'node_key': 'mol3d_to_molecule',
                'node_tasks': [
                    'HERE!!!',
                    'make a custom task that translates a chemthing to a mol'
                ]
            },
            'precursor_keys': ['root'],
        })
        node_specs.append({
            'node': {
                'node_key': 'b3lyp_6_31gs_opt',
                'node_tasks': [
                    self.generate_compute_parse_load_task(
                        flow_params={
                            'label': 'b3lyp_6_31gs_opt_flow',
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.qchem',
                                'job_params': {
                                    'method': 'b3lyp',
                                    'basis': '6_31gs',
                                },
                            },
                            'parse_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.parse'
                            },
                            'load_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.load'
                            }
                        }
                    ),
                ]
            },
            'precursor_keys': ['mol3d_to_molecule'],
        })
        node_specs.append({
            'node': {
                'node_key': 'b3lyp_6_31gs_opt_t1',
                'node_tasks': [
                    self.generate_compute_parse_load_task(
                        flow_params={
                            'label': 'b3lyp_6_31gs_opt_t1_flow',
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.qchem',
                                'job_params': {
                                    'method': 'b3lyp',
                                    'basis': '6_31gs',
                                },
                            },
                            'parse_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.parse'
                            },
                            'load_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.load'
                            }
                        }
                    ),
                ]
            },
            'precursor_keys': ['root'],
        })
        node_specs.append({
            'node': {
                'node_key': 'b3lyp_6_31gs_tddft',
                'node_tasks': [
                    self.generate_compute_parse_load_task(
                        flow_params={
                            'label': 'b3lyp_6_31gs_tddft_flow',
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.qchem',
                                'job_params': {
                                    'method': 'b3lyp',
                                    'basis': '6_31gs',
                                },
                            },
                            'parse_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.parse'
                            },
                            'load_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.load'
                            }
                        }
                    )
                ]
            },
            'precursor_keys': ['b3lyp_6_31gs_opt'],
        })
        dft_flow_spec = {
            'label': 'dft_flow',
            'node_specs': node_specs
        }
        return dft_flow_spec

def generate_flow_spec(*args, flow_params=None, **kwargs):
    generator = ReaxysFlowSpecGenerator(*args, flow_params=flow_params,
                                        **kwargs)
    return generator.generate_flow_spec()
