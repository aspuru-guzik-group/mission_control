import collections
import logging

from mc.task_runners.base_task_runner import BaseTaskRunner

from .flow import Flow


class FlowEngine(object):
    def __init__(self, task_handler=None, logger=None):
        self.logger = logger or logging
        self.task_handler = task_handler
        self.flow_generator_class_registry = collections.OrderedDict()

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

    def tick_flow(self, flow=None, flow_ctx=None):
        if flow.status == 'PENDING': self.start_flow(flow=flow)
        self.start_nearest_pending_nodes(flow=flow)
        self.tick_running_nodes(flow=flow, flow_ctx=flow_ctx)
        if not flow.has_incomplete_nodes():
            self.complete_flow(flow=flow)

    def start_flow(self, flow=None):
        flow.status = 'RUNNING'

    def start_nearest_pending_nodes(self, flow=None):
        for node in flow.get_nearest_pending_nodes():
            self.start_node(flow=flow, node=node)

    def start_node(self, flow=None, node=None):
        node['status'] = 'RUNNING'

    def tick_running_nodes(self, flow=None, flow_ctx=None):
        for node in flow.get_nodes_by_status(status='RUNNING'):
            if self.node_is_running(node=node):
                try:
                    self.tick_node(node=node, flow=flow, flow_ctx=flow_ctx)
                except Exception as exception:
                    self.fail_node(node=node,
                                   error=self.stringify_exception(exception))
            else:
                self.complete_node(node=node)

    def stringify_exception(self, exception=None):
        return '[%s] %s)' % (type(exception), exception)

    def node_is_running(self, node=None):
        return node['status'] == 'RUNNING'

    def tick_node(self, node=None, flow=None, flow_ctx=None):
        try:
            task_runner = self.get_task_runner_for_node(
                node=node, flow=flow, flow_ctx=flow_ctx)
            tasks_status = task_runner.tick_tasks()
            if tasks_status == 'COMPLETED':
                self.complete_node(flow=flow, node=node)
            else:
                node['status'] = tasks_status
        except Exception as error:
            self.fail_node(node=node, error=error)

    def get_task_runner_for_node(self, node=None, flow=None, flow_ctx=None):
        task_runner = BaseTaskRunner(
            get_tasks=(lambda: node.get('node_tasks', [])),
            get_task_context=(
                lambda :self.get_task_context(node=node, flow=flow,
                                              flow_ctx=flow_ctx)
            ),
            task_handler=self.task_handler
        )
        return task_runner

    def get_task_context(self, node=None, flow=None, flow_ctx=None):
        task_context = {'node': node, 'flow': flow, 'flow_ctx': flow_ctx}
        return task_context

    def fail_node(self, node=None, error=None):
        node['error'] = error
        node['status'] = 'FAILED'

    def complete_node(self, flow=None, node=None):
        node['status'] = 'COMPLETED'

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
