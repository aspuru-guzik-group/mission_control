from .compute_parse_load_flow_generator import ComputeParseLoadFlowGenerator

class QChemFlowGenerator():
    flow_type = 'qchem'

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        job_type_prefix = 'a2g2.jobs.qchem'
        compute_parse_load_flow_spec = {
            'compute_job_spec': {
                'job_type': job_type_prefix + '.qchem',
                'job_params': {
                    'qchem_params': flow_spec['input']['qchem_params'],
                }
            },
            'parse_job_spec': {
                'job_type': job_type_prefix + '.parse',
                'job_params': {
                    'precursors': flow_spec['input'].get(
                        'precursors', {}),
                }
            },
            'load_job_spec': {
                'job_type': job_type_prefix + '.load',
            }
        }
        flow = ComputeParseLoadFlowGenerator.generate_flow(
            *args, flow_spec=compute_parse_load_flow_spec, **kwargs
        )
        return flow

    @classmethod
    def get_dependencies(cls):
        return {'flow_generator_classes': set([ComputeParseLoadFlowGenerator])}
