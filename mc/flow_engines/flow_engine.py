import logging

from mc.task_runners.base_task_runner import BaseTaskRunner

from .flow import Flow


class FlowEngine(object):
    simple_flow_serialization_attrs = ['data', 'status', 'root_node_key']

    def __init__(self, task_handler=None, logger=None):
        self.logger = logger or logging
        self.task_handler = task_handler

    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        for attr in self.simple_flow_serialization_attrs:
            setattr(flow, attr, serialized_flow.get(attr, None))
        for node in serialized_flow.get('nodes', []):
            flow.add_node(node=node)
        for edge in serialized_flow.get('edges', []):
            flow.add_edge(edge=edge)
        return flow

    def serialize_flow(self, flow=None):
        serialized_flow = {
            **{attr: getattr(flow, attr, None)
               for attr in self.simple_flow_serialization_attrs},
            'nodes': list(flow.nodes.values()),
            'edges': list(flow.edges.values()),
        }
        return serialized_flow

    def tick_flow(self, flow=None, flow_ctx=None):
        try:
            if flow.status == 'PENDING': self.start_flow(flow=flow)
            self.start_nearest_pending_nodes(flow=flow)
            self.tick_running_nodes(flow=flow, flow_ctx=flow_ctx)
            if not flow.has_incomplete_nodes(): self.complete_flow(flow=flow)
        except Exception as exception:
            self.fail_flow(flow=flow, error=self.stringify_exception(exception))

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
                self.tick_node(node=node, flow=flow, flow_ctx=flow_ctx)
            else:
                self.complete_node(node=node)

    def stringify_exception(self, exception=None):
        return '(%s: %s)' % (type(exception), exception)

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
        except Exception as exception:
            self.fail_node(node=node, error=self.stringify_exception(exception))

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
        msg = "Node with key '{node_key}' failed, error: {error}".format(
            node_key=node.get('node_key', '<unknown key>'),
            error=error
        )
        raise Exception(msg)

    def fail_flow(self, flow=None, error=None):
        flow.data['error'] = error
        flow.status = 'FAILED'

    def complete_node(self, flow=None, node=None):
        node['status'] = 'COMPLETED'

    def complete_flow(self, flow=None):
        flow.status = 'COMPLETED'
