from collections import OrderedDict
from uuid import uuid4

from .workflow import Workflow


class BaseWorkflowRunner(object):
    def __init__(self):
        self.node_class_registry = OrderedDict()

    def register_node_class(self, test_fn=None, node_class=None, key=None):
        if key is None: key = str(uuid4())
        self.node_class_registry[key] = {'test_fn': test_fn,
                                         'node_class': node_class}

    def deserialize_workflow(self, serialized_workflow=None):
        workflow = Workflow()
        nodes = [
            self.deserialize_node(serialized_node)
            for serialized_node in serialized_workflow.get('nodes', [])
        ]
        workflow.add_nodes(nodes=nodes)
        edges = [
            {'src': workflow.nodes[serialized_edge['src_id']],
             'dest': workflow.nodes[serialized_edge['dest_id']]}
            for serialized_edge in serialized_workflow.get('edges', [])
        ]
        workflow.add_edges(edges=edges)
        return workflow

    def deserialize_node(self, serialized_node):
        node_class = self.get_node_class_for_serialized_node(
            serialized_node=serialized_node)
        if not node_class:
            msg = "Could not find node_class for serialized_node '{}'".format( 
                serialized_node)
            raise Exception(msg)
        node_kwargs = self.get_node_kwargs_for_serialized_node(
            serialized_node=serialized_node)
        node = node_class(**node_kwargs)
        return node

    def get_node_class_for_serialized_node(self, serialized_node=None):
        for registry_entry in self.node_class_registry.values():
            if registry_entry['test_fn'](serialized_node):
                return registry_entry['node_class']
        return None

    def get_node_kwargs_for_serialized_node(self, serialized_node=None):
        return serialized_node

    def serialize_workflow(self, workflow=None):
        serialized_workflow = {
            'nodes': self.serialize_nodes(nodes=workflow.nodes.values()),
            'edges': self.serialize_edges(edges=workflow.edges.values()),
            'state': workflow.state,
            'status': workflow.status,
        }
        return serialized_workflow

    def serialize_nodes(self, nodes=None):
        serialized_nodes = [self.serialize_node(node=node) for node in nodes]
        return serialized_nodes

    def serialize_node(self, node=None):
        serialized_node = {
            'id': node.id,
            'state': node.state,
            'type': node.type,
            'status': node.status,
        }
        return serialized_node

    def serialize_edges(self, edges=None):
        serialized_edges = [self.serialize_edge(edge=edge) for edge in edges]
        return serialized_edges

    def serialize_edge(self, edge=None):
        serialized_edge = {
            'src_id': edge['src'].id,
            'dest_id': edge['dest'].id,
        }
        return serialized_edge
