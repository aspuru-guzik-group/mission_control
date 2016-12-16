from collections import OrderedDict
from uuid import uuid4

from .workflow import Workflow


class BaseWorkflowRunner(object):
    def __init__(self):
        self.components_registry = OrderedDict()

    def register_component(self, test_fn=None, component=None, key=None):
        if key is None: key = str(uuid4())
        self.components_registry[key] = {'test_fn': test_fn,
                                         'component': component}

    def deserialize_workflow(self, serialized_workflow=None):
        workflow = Workflow()
        serialized_nodes = serialized_workflow.get('nodes', {})
        nodes = {self.deserialize_node(serialized_node)
                 for serialized_node in serialized_nodes}
        workflow.add_nodes(nodes=nodes)
        for serialized_node in serialized_nodes:
            for link_target_id in serialized_node.get('links_to', []):
                edge_spec = {'src_id': serialized_node['id'],
                             'dest_id': link_target_id}
                workflow.add_edge_from_spec(edge_spec=edge_spec)
        return workflow

    def deserialize_node(self, serialized_node):
        component = self.get_component_for_serialized_node(
            serialized_node=serialized_node)
        if not component:
            msg = "Could not find component for serialized_node '{}'".format( 
                serialized_node)
            raise Exception(msg)
        component_kwargs = self.get_component_kwargs_for_serialized_node(
            serialized_node=serialized_node)
        node = component(**component_kwargs)
        return node

    def get_component_for_serialized_node(self, serialized_node=None):
        for registry_entry in self.components_registry.values():
            if registry_entry['test_fn'](serialized_node):
                return registry_entry['component']
        return None

    def get_component_kwargs_for_serialized_node(self, serialized_node=None):
        return serialized_node
