from flow_engines.flow import Flow

from ..nodes.a2g2_dao import A2G2_DAO_Node
from ..nodes.base import BaseNode
from ..nodes.flow_wrapper import FlowWrapper_Node
from .rpi import RPI_Flow_Module 


class Reaxys_Flow_Module(object):
    flow_type = 'reaxys'
    label = 'Reaxys'
    description = 'Flow that runs reaxys pipeline.'

    @classmethod
    def generate_flow(cls, *arg, mol_key=None, **kwargs):
        node_classes = cls.get_node_classes()
        flow = Flow()

        query_node = node_classes['reaxys:a2g2_dao'](
            id='fetch_mol',
            data={
                'input': {
                    'query': {
                        'object': 'Mol',
                        'filters': [('key', mol_key)]
                    }
                }
            })
        flow.add_node(node=query_node, as_root=True)

        confgen_setup_node = node_classes['reaxys:confgen_setup'](
            id='confgen_setup',
            data={
                'input': {
                    'src_node_id': 'fetch_mol',
                    'dest_node_id': 'confgen',
                }
            })
        flow.add_node(node=confgen_setup_node, precursor='fetch_mol')
        confgen_flow_node = node_classes['reaxys:flow'](
            id='confgen',
            flow=RPI_Flow_Module.generate_flow()
        )
        flow.add_node(confgen_flow_node, precursor='confgen_setup')

    @classmethod
    def get_node_classes(cls):
        node_classes = {
            'reaxys:a2g2_dao': A2G2_DAO_Node,
            'reaxys:reaxys_confgen_setup': Reaxys_Confgen_Setup_Node,
            'reaxys:flow': FlowWrapper_Node
        }
        return node_classes

class Reaxys_Confgen_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.flow.nodes[self.data['input']['src_node_id']]
        dest_node = self.flow.nodes[self.data['input']['dest_node_id']]
        dest_node.data['input'] = {
            'smiles': src_node.data['output']['result_set']['mols'][0]['smiles'],
            'confgen_spec': '!!!@TODO!!!'

        }
        self.status = 'COMPLETED'
