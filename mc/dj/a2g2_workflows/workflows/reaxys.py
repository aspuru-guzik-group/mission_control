from workflow_engines.workflow import Workflow

from ..nodes.base import BaseNode


def generate_reaxys_workflow(mol_key=None):
    workflow = Workflow()

    root_node = BaseNode(id='root', data={'output': {'mol_key': mol_key}},
                         status='COMPLETED')
    workflow.add_node(root_node)
    workflow.root_node = root_node

    return workflow

