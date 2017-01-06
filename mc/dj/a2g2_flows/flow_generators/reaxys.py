from flow_engines.flow import Flow

from ..tasks.a2g2_dao import A2G2_DAO_Task
from ..tasks.base import BaseTask
from ..tasks.flow import FlowTask
from ..flows.run_and_load import RunAndLoadFlowGenerator


class ReaxysFlowGenerator(object):
    flow_type = 'reaxys'
    label = 'Reaxys'
    description = 'Flow that runs reaxys pipeline.'

    @classmethod
    def generate_flow(cls, *arg, mol_key=None, **kwargs):
        flow = Flow()
        flow.add_task(key='fetch_mol', 
                      as_root=True,
                      task=A2G2_DAO_Task(),
                      static_input={
                          'query': {
                              'object': 'Mol', 'filters': [('key', mol_key)]
                          }
                      })
        flow.add_task(key='confgen_prep',
                      precursor='fetch_mol',
                      task=ConfgenPrepTask())
        flow.add_task(key='confgen',
                      precursor='confgen_prep',
                      task=FlowTask())

    @classmethod
    def get_dependencies(cls):
        return {
            'tasks': set([
                A2G2_DAO_Task,
                ConfgenPrepTask,
                FlowTask,
            ]),
            'flow_generators': set([RunAndLoadFlowGenerator])
        }

class ConfgenPrepTask(BaseTask):
    def tick(self, *args, **kwargs):
        self.output = {
            'smiles': self.input['result_set']['mols'][0]['smiles'],
            'confgen_spec': '!!!@TODO!!!'
        }
        self.status = 'COMPLETED'
