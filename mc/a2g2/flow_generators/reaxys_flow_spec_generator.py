class ReaxysFlowSpecGenerator(object):
    def __init__(self, *args, flow_params=None, **kwargs):
        self.flow_params = flow_params

    def generate_flow_spec(self):
        node_specs = []
        #node_specs.append(self.generate_confgen_node_spec())
        node_specs.append(self.generate_conformer_query_node_spec())
        node_specs.append(self.generate_conformer_demux_node_spec())
        flow_spec = {'label': 'reaxys_flow', 'node_specs': node_specs}
        return flow_spec

    def generate_confgen_node_spec(self):
       return { 
            'node': {
                'node_key': 'confgen',
                'node_tasks': [
                    self.generate_compute_parse_load_task(
                        flow_params={
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.confgen.confgen',
                            },
                            'parse_job_spec': {
                                'job_type': 'a2g2.jobs.confgen.parse',
                            },
                            'load_job_spec': {
                                'job_type': 'a2g2.jobs.confgen.load',
                            },
                        }),
                ],
                'data': { 
                    # @TODO! REPLACE!
                    #'conformer_chemthings': ['FAKE_%s' % i for i in range(3)]
                    'conformer_chemthings': [1],
                }
            },
            'precursor_keys': ['ROOT'],
        }

    def generate_compute_parse_load_task(self, task_key=None, flow_params=None):
        return {
            'task_key': task_key or 'compute_parse_load_task',
            'task_type': ('a2g2.task_handlers.nodes'
                          '.run_compute_parse_load_flow_handler'),
            'task_params': {
                'flow_params': {
                    'delete_flow_spec_after_creation': True,
                    **flow_params,
                }
            }
        }

    def generate_conformer_query_node_spec(self):
        return {
            'node': {
                'node_key': 'conformer_query',
                'node_tasks': [
                    {
                        'task_key': 'run_conformer_query_job',
                        'task_params': {
                            'job_spec': {
                                'job_type': 'a2g2.jobs.a2g2_db.query',
                                'job_params': {},
                            }
                        },
                        'task_type': ('a2g2.task_handlers.nodes'
                                      '.run_job_task_handler')
                    },
                    {
                        'task_key': 'expose_chemthings',
                        'task_params': {
                            'value_specs': [
                                {
                                    'from_json': True,
                                    'source': ('ctx.tasks'
                                               '.run_conformer_query_job.data'
                                               '.stdout'),
                                    'dest': 'ctx.task.scratch.query_results',
                                },
                                {
                                    'source': ('ctx.task.scratch'
                                               '.query_results.hits'),
                                    'dest': 'ctx.node.data.conformer_chemthings'
                                }
                            ]
                        },
                        'task_type': 'mc.task_handlers.set_values_task_handler'
                    },
                ]
            },
            'precursor_keys': ['ROOT']
        }
    
    def generate_conformer_demux_node_spec(self):
        return {
            'node': {'node_key': 'confgen_demux',
                     'node_tasks': self.generate_confgen_demux_tasks()},
            'precursor_keys': ['conformer_query']
        }

    def generate_confgen_demux_tasks(self):
        confgen_demux_tasks = [
            {
                'task_key': 'wire_confgen_demux',
                'task_type': 'mc.task_handlers.set_values_task_handler',
                'task_params': {
                    'value_specs': [
                        {
                            'source': ('ctx.flow.nodes.conformer_query.data'
                                       '.conformer_chemthings'),
                            'dest': 'ctx.tasks.confgen_demux.task_params.items',
                        },
                    ]
                }
            },
            self.generate_demux_task()
        ]
        return confgen_demux_tasks

    def generate_demux_task(self):
        demux_task = {
            'task_key': 'confgen_demux',
            'task_type': 'a2g2.task_handlers.nodes.demux_task_handler',
            'task_params': {
                'items': 'TO BE WIRED',
                'demux_flow_params': {
                    'cfg': {
                        'fail_fast': False,
                    },
                    'data': {
                        'run_dft_flow_params': (
                            self.generate_run_dft_flow_params())
                    }
                },
                'node_spec_template': self.generate_demux_node_spec_template(),
                'substitutions': [
                    {
                        'source': 'ctx.demux_item',
                        'dest': 'ctx.node_spec.node.data.conformer_chemthing',
                    }
                ]
            }
        }
        return demux_task

    def generate_demux_node_spec_template(self):
        node_spec_template = {
            'node': {
                'data': {
                    'conformer_chemthing': 'TO BE WIRED',
                },
                'node_tasks': [
                    self.generate_molecule_to_chemthing_task(),
                    {
                        'task_key': 'wire_dft_flow_params_into_task',
                        'task_type': 'mc.task_handlers.set_values_task_handler',
                        'task_params': {
                            'value_specs': [
                                {
                                    'source': ('ctx.flow.data'
                                               '.run_dft_flow_params'),
                                    'dest': ('ctx.tasks.run_dft_flow'
                                             '.task_params')
                                },
                            ]
                        }
                    },
                    {
                        'task_key': 'wire_molecule_into_flow_spec',
                        'task_type': 'mc.task_handlers.set_values_task_handler',
                        'task_params': {
                            'value_specs': [
                                {
                                    'source': 'ctx.node.data.molecule',
                                    'dest': ('ctx.tasks.run_dft_flow'
                                             '.task_params.flow_spec'
                                             '.data.molecule')
                                },
                            ]
                        }
                    },
                    {
                        'task_key': 'run_dft_flow',
                        'task_type': ('a2g2.task_handlers.nodes'
                                      '.run_flow_task_handler'),
                        'task_params': 'TO BE WIRED',
                    },
                ]
            }
        }
        return node_spec_template

    def generate_molecule_to_chemthing_task(self):
        molecule_to_chemthing_task = {
            'task_key': 'generate_molecule_from_chemthing',
            'task_type': 'mc.task_handlers.set_values_task_handler',
            'task_params': {
                'value_specs': [
                    {
                        'from_value': True,
                        'value': {
                            'multiplicity': 0,
                            'charge': 0,
                        },
                        'dest': 'ctx.node.data.molecule'
                    },
                    {
                        'source': ('ctx.node.data'
                                   '.conformer_chemthing'
                                   '.props.a2g2:prop:atoms'
                                  ),
                        'dest': ('ctx.node.data.molecule'
                                 '.atoms')
                    },
                ]
            }
        }
        return molecule_to_chemthing_task

    def generate_run_dft_flow_params(self):
        run_dft_flow_params = {
            'flow_spec': self.generate_dft_flow_spec(),
            'delete_flow_spec_after_creation': True,
        }
        return run_dft_flow_params

    def generate_dft_flow_spec(self):
        dft_flow_spec = {
            'label': 'dft_flow',
            'node_specs': [
                self.generate_b3lyp_opt_node_spec(),
                #self.generate_b3lyp_opt_t1_node_spec(),
                #self.generate_b3lyp_tddft_node_spec(),
            ]
        }
        return dft_flow_spec

    def generate_b3lyp_opt_node_spec(self):
        node_key = 'b3lyp_6_31gs_opt'
        cpl_key = node_key + '_cpl'
        return {
            'node': {
                'node_key': node_key,
                'node_tasks': [
                    self.generate_wire_molecule_into_cpl_task(cpl_key=cpl_key),
                    self.generate_compute_parse_load_task(
                        task_key=cpl_key,
                        flow_params={
                            'label': 'b3lyp_6_31gs_opt_flow',
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.qchem',
                                'job_params': {
                                    'qchem_params': {
                                        'method': 'b3lyp',
                                        'basis': '6_31gs',
                                        'molecule': 'TO BE WIRED',
                                    }
                                },
                            },
                            'parse_load_job_spec': {
                                'job_type': 'a2g2.jobs.flow.task_list',
                                'job_params': {
                                    'tasks': [
                                        {
                                            'task_type': (
                                                'mc.a2g2.job_modules.qchem'
                                                '.task_handlers'
                                                '.parse_computation_meta'
                                                '_task_handler'
                                            ),
                                        },
                                        {
                                            'task_type': (
                                                'mc.a2g2.job_modules.qchem'
                                                '.task_handlers'
                                                '.parse_opt_coords_task_handler'
                                            ),
                                        },
                                        {
                                            'task_type': (
                                                'mc.task_handlers.'
                                                'log_task_handler'
                                            ),
                                        }
                                    ]
                                }
                            },
                        }
                    ),
                ]
            },
            'precursor_keys': ['ROOT'],
        }

    def generate_wire_molecule_into_cpl_task(self, cpl_key=None):
        return {
            'task_key': 'wire_molecule_into_compute_parse_load',
            'task_type': 'mc.task_handlers.set_values_task_handler',
            'task_params': {
                'value_specs': [{
                    'source': 'ctx.flow.data.molecule',
                    'dest': ('ctx.tasks.{cpl_key}'
                             '.task_params.flow_params'
                             '.compute_job_spec.job_params.qchem_params'
                             '.molecule').format(cpl_key=cpl_key)
                }]
            }
        }

    def generate_b3lyp_opt_t1_node_spec(self):
        node_key = 'b3lyp_6_31gs_opt_t1'
        cpl_key = node_key + '_cpl'
        return {
            'node': {
                'node_key': 'b3lyp_6_31gs_opt_t1',
                'node_tasks': [
                    self.generate_wire_molecule_into_cpl_task(cpl_key=cpl_key),
                    self.generate_compute_parse_load_task(
                        task_key=cpl_key,
                        flow_params={
                            'label': 'b3lyp_6_31gs_opt_t1_flow',
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.qchem',
                                'job_params': {
                                    'qchem_params': {
                                        'method': 'b3lyp',
                                        'basis': '6_31gs',
                                    }
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
            'precursor_keys': ['ROOT'],
        }


    def generate_b3lyp_tddft_node_spec(self):
        return {
            'node': {
                'node_key': 'b3lyp_6_31gs_tddft',
                'node_tasks': [
                    self.generate_compute_parse_load_task(
                        flow_params={
                            'label': 'b3lyp_6_31gs_tddft_flow',
                            'compute_job_spec': {
                                'job_type': 'a2g2.jobs.qchem.qchem',
                                'job_params': {
                                    'qchem_params': {
                                        'method': 'b3lyp',
                                        'basis': '6_31gs',
                                    },
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
        }

def generate_flow_spec(*args, flow_params=None, **kwargs):
    generator = ReaxysFlowSpecGenerator(*args, flow_params=flow_params,
                                        **kwargs)
    return generator.generate_flow_spec()
