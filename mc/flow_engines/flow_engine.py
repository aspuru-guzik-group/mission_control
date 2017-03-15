import collections
import logging

from .flow import Flow


class FlowEngine(object):
    def __init__(self, action_processor=None, logger=None):
        self.action_processor = action_processor
        self.logger = logger or logging
        self.node_engine_registry = collections.OrderedDict()
        self.flow_generator_class_registry = collections.OrderedDict()

    def register_node_engine(self, node_engine=None, key=None):
        if not key:
            key = getattr(node_engine, 'node_engine',
                          node_engine.__class__.__name__)
        self.node_engine_registry[key] = node_engine

    def register_flow_generator_class(self, flow_generator_class=None):
        key = getattr(flow_generator_class, 'flow_type',
                      flow_generator_class.__name__)
        self.flow_generator_class_registry[key] = flow_generator_class
        self.register_flow_generator_class_dependencies(
            flow_generator_class=flow_generator_class)

    def register_flow_generator_class_dependencies(self,
                                                   flow_generator_class=None):
        if not hasattr(flow_generator_class, 'get_dependencies'): return
        dependencies = flow_generator_class.get_dependencies()
        for node_engine in dependencies.get('node_engines', []):
            self.register_node_engine(node_engine)
        for generator_class in dependencies.get('flow_generator_classes', []):
            self.register_flow_generator_class(
                flow_generator_class=generator_class)

    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        for attr in ['data', 'input', 'output', 'status', 'root_node_key']:
            setattr(flow, attr, serialized_flow.get(attr, None))
        for node in serialized_flow.get('nodes', []):
            flow.add_node(node=node)
        for edge in serialized_flow.get('edges', []):
            flow.add_edge(edge=edge)
        return flow

    def serialize_flow(self, flow=None):
        serialized_flow = {
            **{attr: getattr(flow, attr, None)
               for attr in ['data', 'input', 'output', 'status',
                            'root_node_key']},
            'nodes': list(flow.nodes.values()),
            'edges': list(flow.edges.values()),
        }
        return serialized_flow

    def tick_flow(self, flow=None, ctx=None):
        if flow.status == 'PENDING': self.start_flow(flow=flow)
        self.start_nearest_pending_nodes(flow=flow)
        self.tick_running_nodes(flow=flow, ctx=ctx)
        if not flow.has_incomplete_nodes():
            self.complete_flow(flow=flow)

    def start_flow(self, flow=None):
        flow.status = 'RUNNING'

    def start_nearest_pending_nodes(self, flow=None):
        for node in flow.get_nearest_pending_nodes():
            self.start_node(flow=flow, node=node)

    def start_node(self, flow=None, node=None):
        pre_start_actions = node.get('pre_start_actions', None)
        if pre_start_actions: self.process_actions(
            actions=pre_start_actions,
            ctx={'flow': flow, 'node': node}
        )
        node['status'] = 'RUNNING'

    def process_actions(self, actions=None, ctx=None):
        for action in actions:
            self.action_processor.process_action(action=action, ctx=ctx)

    def tick_running_nodes(self, flow=None, ctx=None):
        for node in flow.get_nodes_by_status(status='RUNNING'):
            self.tick_node(flow=flow, node=node, ctx=ctx)

    def tick_node(self, flow=None, node=None, ctx=None):
        ctx = ctx or {}
        try:
            node_engine = self.get_engine_for_node(node=node)
            node_engine.tick_node(node=node, ctx={**ctx, 'flow': flow})
            if node.get('status', None) == 'COMPLETED':
                self.complete_node(flow=flow, node=node)
        except Exception as e:
            msg = "tick failed for node with key '{key}': {e}".format(
                key=node['key'], e=e)
            self.logger.exception(msg)

    def complete_node(self, flow=None, node=None):
        actions = node.get('post_complete_actions', None)
        if actions: self.process_actions(
            actions=actions,
            ctx={'flow': flow, 'node': node}
        )

    def get_engine_for_node(self, node=None):
        node_engine = self.node_engine_registry.get(node['node_engine'], None)
        if not node_engine:
            raise Exception("Could not find node_engine for node '%s'" % node)
        return node_engine

    def complete_flow(self, flow=None):
        flow.output = self.get_flow_output(flow=flow)
        flow.status = 'COMPLETED'

    def get_flow_output(self, flow=None):
        _get_node_output = lambda t: t.get('output', None)
        tail_nodes = flow.get_tail_nodes()
        if len (tail_nodes) == 0: output = None
        elif len(tail_nodes) == 1: output = _get_node_output(tail_nodes[0])
        elif len(tail_nodes) > 1: output = {
            tail_node['key']: _get_node_output(tail_node)
            for tail_node in tail_nodes
        }
        return output

    def generate_flow(self, flow_spec=None):
        FlowGenerator = self.get_flow_generator_class_for_spec(flow_spec)
        if not FlowGenerator:
            msg = "Could not find FlowGenerator for flow_spec '{}'".format( 
                flow_spec)
            raise Exception(msg)
        flow_kwargs = self.get_flow_kwargs_for_flow_spec(flow_spec=flow_spec)
        flow = FlowGenerator.generate_flow(**flow_kwargs)
        return flow

    def get_flow_generator_class_for_spec(self, flow_spec=None):
        return self.flow_generator_class_registry.get(
            flow_spec['flow_type'], None)

    def get_flow_kwargs_for_flow_spec(self, flow_spec=None):
        return {'flow_spec': flow_spec}
