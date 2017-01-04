from workflow_engines.workflow import Workflow

from ..nodes.a2g2_dao import A2G2_DAO_Node
from ..nodes.base import BaseNode
from ..nodes.b3lyp import B3LYP_Node
from ..nodes.confgen import Confgen_Node, Confgen_Parse_Node

node_classes = {}

class Reaxys_Workflow_Module(object):
    workflow_type = 'reaxys'
    label = 'Reaxys'
    description = 'Workflow that runs reaxys pipeline.'

    @classmethod
    def generate_workflow(cls, *arg, mol_key=None, **kwargs):
        node_classes = cls.get_node_classes()
        workflow = Workflow()

        # Confgen Nodes
        workflow.add_node(
            node=node_classes['A2G2_DAO'](
                id='fetch_mol',
                data={
                    'input': {
                        'query': {
                            'object': 'Mol',
                            'filters': [('key', mol_key)]
                        }
                    }
                }
            ),
            as_root=True
        )
        workflow.add_node(
            node=node_classes['Reaxys_Confgen_Setup'](
                id='confgen_setup',
                data={
                    'input': {
                        'src_node_id': 'fetch_mol',
                        'dest_node_id': 'confgen',
                    }
                }
            ),
            precursor='fetch_mol'
        )
        workflow.add_node(
            node=node_classes['Confgen'](id='confgen'),
            precursor='confgen_setup'
        )
        workflow.add_node(
            node=node_classes['Reaxys_Confgen_Parse_Setup'](
                id='confgen_parse_setup',
                data={
                    'input': {
                        'src_node_id': 'confgen',
                        'dest_node_id': 'confgen_parse',
                    }
                }
            ),
            precursor='confgen'
        )
        workflow.add_node(
            node=node_classes['Confgen_Parse'](id='confgen_parse'),
            precursor='confgen_parse_setup'
        )
        workflow.add_node(
            node=node_classes['Reaxys_Confgen_Ingest_Setup'](
                id='confgen_ingest_setup',
                data={
                    'input': {
                        'src_node_id': 'confgen_parse',
                        'dest_node_id': 'confgen_ingest',
                        'mol_key': mol_key,
                    }
                }
            ),
            precursor='confgen_parse'
        )
        workflow.add_node(
            node=node_classes['A2G2_DAO'](id='confgen_ingest'),
            precursor='confgen_ingest_setup'
        )

        # B3LYP Nodes
        #workflow.add_node(
            #node=node_classes['Reaxys_B3LYP_Setup'](id='b3lyp_setup'),
            #precursor='confgen_ingest'
        #)
        #workflow.add_node(
            #node=node_classes['B3LYP'](id='b3lyp'),
            #precursor='b3lyp_setup'
        #)
        return workflow

    @classmethod
    def get_node_classes(cls):
        node_classes = {
            'A2G2_DAO': A2G2_DAO_Node,
            'Reaxys_Confgen_Setup': Reaxys_Confgen_Setup_Node,
            'Confgen': Confgen_Node,
            'Reaxys_Confgen_Parse_Setup': Reaxys_Confgen_Parse_Setup_Node,
            'Confgen_Parse': Confgen_Parse_Node,
            'Reaxys_Confgen_Ingest_Setup': Reaxys_Confgen_Ingest_Setup_Node,
            'Reaxys_B3LYP_Setup': Reaxys_B3LYP_Setup_Node,
            'B3LYP': B3LYP_Node,
        }
        return node_classes

class Reaxys_Confgen_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.workflow.nodes[self.data['input']['src_node_id']]
        dest_node = self.workflow.nodes[self.data['input']['dest_node_id']]
        dest_node.data['input'] = {
            'smiles': src_node.data['output']['result_set']['mols'][0]['smiles'],
            'confgen_spec': '!!!@TODO!!!'

        }
        self.status = 'COMPLETED'

class Reaxys_Confgen_Parse_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.workflow.nodes[self.data['input']['src_node_id']]
        dest_node = self.workflow.nodes[self.data['input']['dest_node_id']]
        dest_node.data['input'] = {
            'parse': {'dir': src_node.data['output']['dir']}
        }
        self.status = 'COMPLETED'

class Reaxys_Confgen_Ingest_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.workflow.nodes[self.data['input']['src_node_id']]
        dest_node = self.workflow.nodes[self.data['input']['dest_node_id']]
        dest_node.data['input'] = {
            'ingest': {
                'object': 'Geom',
                'item_kwargs': {'mol': self.data['input']['mol_key']},
                'dir': src_node.data['output']['dir'],
            }
        }
        self.status = 'COMPLETED'

class Reaxys_B3LYP_Setup_Node(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.workflow.nodes[self.data['input']['src']]
        dest_node = self.workflow.nodes[self.data['input']['dest']]
        src_node.data['input']['smiles'] = dest_node['output']['mol']['smiles']
        self.status = 'COMPLETED'
