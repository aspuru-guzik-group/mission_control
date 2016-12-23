from workflow_engines.workflow import Workflow

from ..nodes.a2g2_dao import A2G2_DAO_Node
from ..nodes.base import BaseNode

class WireConfgenNode(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.workflow.nodes[self.data['input']['src']]
        dest_node = self.workflow.nodes[self.data['input']['dest']]
        src_node.data['input']['smiles'] = dest_node['output']['mol']['smiles']
        self.status = 'COMPLETED'

class ReaxysWorkflowGenerator(object):
    def generate_workflow(self, mol_key=None):
        self.mol_key = mol_key
        self.workflow = Workflow()
        self.nodes = self.workflow.nodes # shortcut
        self.add_fetch_mol_node()
        self.add_wire_confgen_node()
        return self.workflow

    def add_fetch_mol_node(self):
        fetch_mol_query = {'object': 'Mol', 'filters': [('key', self.mol_key)]}
        fetch_mol_node = A2G2_DAO_Node(
            id='fetch_mol',
            status='PENDING',
            data={'input': {'query': fetch_mol_query}})
        self.workflow.add_node(fetch_mol_node)
        self.workflow.root_node = fetch_mol_node

    def add_wire_confgen_node(self):
        self.workflow.add_node(WireConfgenNode(
            id='wire_confgen',
            data={'input': {'src': 'fetch_mol', 'dest': 'confgen'}}))

def generate_reaxys_workflow(mol_key=None):
    generator = ReaxysWorkflowGenerator()
    return generator.generate_workflow(mol_key=mol_key)
