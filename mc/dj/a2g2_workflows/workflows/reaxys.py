from workflow_engines.workflow import Workflow

from ..nodes.a2g2_dao import A2G2_DAO_Node
from ..nodes.base import BaseNode

node_classes = {}

class Reaxys_WorkflowModule(object):
    workflow_type = 'reaxys'
    label = 'Reaxys'
    description = 'Workflow that runs reaxys pipeline.'

    @classmethod
    def get_node_classes(cls):
        return node_classes

    @classmethod
    def generate_workflow(cls, *arg, mol_key=None, **kwargs):
        node_classes = cls.get_node_classes()
        workflow = Workflow()
        workflow.add_node(
            node_classes['A2G2_DAO'](id='fetch_mol', data=None),
            root=True)
        workflow.add_node(
            node_classes['Reaxys_ConfgenSetupNode'](id='confgen_setup'),
            precursor='fetch_mol')
        workflow.add_node(
            node_classes['Confgen'](id='confgen'),
            precursor='confgen_setup')
        workflow.add_node(
            node_classes['Reaxys_B3LYPSetup'](id='b3lyp_setup'),
            precursor='confgen')
        workflow.add_node(
            node_classes['B3LYPNode'](id='b3lyp'),
            precursor='b3lyp_setup')

class ConfgenSetupNode(BaseNode):
    def tick(self, *args, **kwargs):
        src_node = self.workflow.nodes[self.data['input']['src']]
        dest_node = self.workflow.nodes[self.data['input']['dest']]
        src_node.data['input']['smiles'] = dest_node['output']['mol']['smiles']
        self.status = 'COMPLETED'
