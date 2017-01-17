from flow_engines.flow import Flow

from ..tasks.flow import FlowTask
from ..flow_generators.run_and_load import RunAndLoadFlowGenerator


class ReaxysFlowGenerator(object):
    flow_type = 'reaxys'

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        flow = Flow()
        flow.add_task(key='confgen', 
                      as_root=True,
                      task=FlowTask(),
                      input={
                          'flow_spec': {
                              'flow_type': RunAndLoadFlowGenerator.flow_type,
                              'run_spec': {
                                  'job_type': 'confgen',
                                  'confgen': {
                                      'smiles': flow_spec['smiles'],
                                      'params': flow_spec['confgen_params'],
                                  }
                              },
                              'load_spec': {
                                  'job_type': 'confgen:load',
                              }
                          }
                      })
        return flow

    @classmethod
    def get_dependencies(cls):
        return {
            'task_classes': set([FlowTask]),
            'flow_generator_classes': set([RunAndLoadFlowGenerator])
        }
