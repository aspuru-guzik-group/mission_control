import json
import logging
import traceback

from mc.task_runners.base_task_runner import BaseTaskRunner

from .flow import Flow


class FlowEngine(object):
    simple_flow_serialization_attrs = ['data', 'label', 'status', 'cfg']

    class NodeError(Exception): pass

    def __init__(self, task_handler=None, logger=None, max_msg_len=1024):
        self.logger = logger or logging
        self.task_handler = task_handler
        self.max_msg_len = max_msg_len

    @classmethod
    def generate_flow(self, flow_spec=None):
        flow_kwargs = {k:v for k, v in flow_spec.items()
                       if k in self.simple_flow_serialization_attrs}
        flow = Flow(**flow_kwargs)
        for node_spec in flow_spec.get('node_specs', []):
            flow.add_node(**node_spec)
        return flow

    @classmethod
    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        flow_dict = json.loads(serialized_flow or '{}')
        for attr in self.simple_flow_serialization_attrs:
            setattr(flow, attr, flow_dict.get(attr, None))
        if flow.data is None: flow.data = {}
        for node_key, node in flow_dict.get('nodes', {}).items():
            if node_key == Flow.ROOT_NODE_KEY: continue
            flow.add_node(node=node)
        for edge in flow_dict.get('edges', []):
            flow.add_edge(edge=edge)
        return flow


    @classmethod
    def serialize_flow(self, flow=None):
        flow_dict = {
            **{attr: getattr(flow, attr, None)
               for attr in self.simple_flow_serialization_attrs},
            'nodes': {node_key: node for node_key, node in flow.nodes.items()
                      if node_key != Flow.ROOT_NODE_KEY},
            'edges': [edge for edge in flow.edges.values()],
        }
        return json.dumps(flow_dict)

    def run_flow(self, flow=None, flow_ctx=None):
        completed_statuses = set(['COMPLETED', 'FAILED'])
        while flow.status not in completed_statuses:
            self.tick_flow(flow=flow, flow_ctx=flow_ctx)

    def tick_flow(self, flow=None, flow_ctx=None):
        try:
            if flow.status == 'PENDING': self.start_flow(flow=flow)
            self.start_nearest_pending_nodes(flow=flow)
            self.tick_running_nodes(flow=flow, flow_ctx=flow_ctx)
            if not flow.has_incomplete_nodes(): self.complete_flow(flow=flow)
        except Exception as exception:
            fail_flow = True
            self.append_flow_error(flow=flow, error=traceback.format_exc())
            if isinstance(exception, self.NodeError):
                if not flow.cfg.get('fail_fast', True): fail_flow = False
            if fail_flow: self.fail_flow(flow=flow)

    def start_flow(self, flow=None):
        flow.status = 'RUNNING'

    def start_nearest_pending_nodes(self, flow=None):
        for node in flow.get_nearest_pending_nodes():
            try: self.start_node(flow=flow, node=node)
            except: self.fail_node(node=node, error=traceback.format_exc())

    def start_node(self, flow=None, node=None):
        node['status'] = 'RUNNING'

    def tick_running_nodes(self, flow=None, flow_ctx=None):
        for node in flow.get_nodes_by_status(status='RUNNING'):
            if self.node_is_running(node=node):
                self.tick_node(node=node, flow=flow, flow_ctx=flow_ctx)
            else:
                self.complete_node(node=node)

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
            self.fail_node(node=node, error=traceback.format_exc())

    def get_task_runner_for_node(self, node=None, flow=None, flow_ctx=None):
        task_runner = BaseTaskRunner(
            get_tasks=(lambda: node.get('node_tasks', [])),
            get_task_context=(
                lambda: self.get_task_context(node=node, flow=flow,
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
        raise self.NodeError(msg)

    def fail_flow(self, flow=None, error=None):
        if error: self.append_flow_error(flow=flow, error=error)
        flow.status = 'FAILED'

    def append_flow_error(self, flow=None, error=None):
        flow.data.setdefault('errors', [])
        flow.data['errors'].append(self.elide_text(error))

    def elide_text(self, text=None, max_len=None):
        max_len = max_len or self.max_msg_len
        if len(text) > max_len: text = text[0:max_len] + '...'
        return text

    def complete_node(self, flow=None, node=None):
        node['status'] = 'COMPLETED'

    def complete_flow(self, flow=None):
        if flow.data.get('errors'): self.fail_flow(flow=flow)
        else: flow.status = 'COMPLETED'
