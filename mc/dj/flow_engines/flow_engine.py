from collections import OrderedDict
from uuid import uuid4

from .flow import Flow


class FlowEngine(object):
    def __init__(self, ctx=None):
        self.node_class_registry = OrderedDict()
        self.ctx = ctx

    def register_node_class(self, node_class=None, key=None, test_fn=None):
        if key is None: key = str(uuid4())
        if not test_fn:
            def test_fn(serialized_node):
                return (serialized_node.get('node_type') == node_class.__name__)
        self.node_class_registry[key] = {'test_fn': test_fn,
                                         'node_class': node_class}

    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        nodes = [
            self.deserialize_node(serialized_node)
            for serialized_node in serialized_flow.get('nodes', [])
        ]
        flow.add_nodes(nodes=nodes)
        root_node_id = serialized_flow.get('root_node_id')
        if root_node_id is not None:
            flow.root_node = flow.nodes[root_node_id]
        edges = [
            {'src': flow.nodes[serialized_edge['src_id']],
             'dest': flow.nodes[serialized_edge['dest_id']]}
            for serialized_edge in serialized_flow.get('edges', [])
        ]
        flow.add_edges(edges=edges)
        flow.status = serialized_flow['status']
        return flow

    def deserialize_node(self, serialized_node):
        node_class = self.get_node_class_for_serialized_node(
            serialized_node=serialized_node)
        if not node_class:
            msg = "Could not find node_class for serialized_node '{}'".format( 
                serialized_node)
            raise Exception(msg)
        node_kwargs = self.get_node_kwargs_for_serialized_node(
            serialized_node=serialized_node)
        node = node_class(ctx=self.ctx, **node_kwargs)
        return node

    def get_node_class_for_serialized_node(self, serialized_node=None):
        for registry_entry in self.node_class_registry.values():
            if registry_entry['test_fn'](serialized_node):
                return registry_entry['node_class']
        return None

    def get_node_kwargs_for_serialized_node(self, serialized_node=None):
        return serialized_node

    def serialize_flow(self, flow=None):
        serialized_flow = {
            'nodes': self.serialize_nodes(nodes=flow.nodes.values()),
            'edges': self.serialize_edges(edges=flow.edges.values()),
            'state': flow.state,
            'status': flow.status,
        }
        if flow.root_node is not None:
            serialized_flow['root_node_id'] = flow.root_node.id
        return serialized_flow

    def serialize_nodes(self, nodes=None):
        serialized_nodes = [self.serialize_node(node=node) for node in nodes]
        return serialized_nodes

    def serialize_node(self, node=None):
        serialized_node = {
            'id': node.id,
            'state': node.state,
            'node_type': node.node_type,
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

    def tick_flow(self, flow=None, ctx=None):
        self.start_nearest_pending_nodes(flow=flow)
        self.tick_running_nodes(flow=flow, ctx=ctx)
        if not flow.has_incomplete_nodes():
            flow.status = 'COMPLETED'

    def start_nearest_pending_nodes(self, flow=None):
        for node in flow.get_nearest_pending_nodes():
            node.status = 'RUNNING'

    def tick_running_nodes(self, flow=None, ctx=None):
        ctx_with_jobs = {**ctx, 'jobs': flow.jobs}
        for node in flow.get_nodes_by_status(status='RUNNING'):
            self.tick_node(node=node, ctx=ctx_with_jobs)

    def tick_node(self, node=None, ctx=None):
        node.tick(engine=self, ctx=ctx)

