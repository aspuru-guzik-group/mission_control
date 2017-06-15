import json
import logging
import traceback

from .flow import Flow


class FlowEngine(object):
    simple_flow_serialization_attrs = ['data', 'label', 'status', 'cfg']

    class FlowError(Exception):
        def __init__(self, *args, flow=None, **kwargs):
            self.flow = flow
            super().__init__(*args, **kwargs)

    class TaskError(Exception): pass

    def __init__(self, task_handler=None, logger=None, max_msg_len=None):
        self.logger = logger or logging
        self.task_handler = task_handler or self.get_default_task_handler()
        self.max_msg_len = max_msg_len

    def get_default_task_handler(self):
        from mc.task_handlers.module_dispatch_task_handler import \
                ModuleDispatchTaskHandler
        return ModuleDispatchTaskHandler()

    @classmethod
    def generate_flow(self, flow_spec=None):
        flow_kwargs = {k:v for k, v in flow_spec.items()
                       if k in self.simple_flow_serialization_attrs}
        flow = Flow(**flow_kwargs)
        for task in flow_spec.get('tasks', []): flow.add_task(task=task)
        return flow

    @classmethod
    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        flow_dict = json.loads(serialized_flow or '{}')
        for attr in self.simple_flow_serialization_attrs:
            setattr(flow, attr, flow_dict.get(attr, None))
        if flow.data is None: flow.data = {}
        for key, task in flow_dict.get('tasks', {}).items():
            if key == Flow.ROOT_TASK_KEY: continue
            flow.add_task(task=task)
        for edge in flow_dict.get('edges', []):
            flow.add_edge(edge=edge)
        return flow

    @classmethod
    def serialize_flow(self, flow=None):
        flow_dict = {
            **{attr: getattr(flow, attr, None)
               for attr in self.simple_flow_serialization_attrs},
            'tasks': {key: task for key, task in flow.tasks.items()
                      if key != Flow.ROOT_TASK_KEY},
            'edges': [edge for edge in flow.edges.values()],
        }
        return json.dumps(flow_dict)

    def run_flow(self, flow=None, task_context=None, check=False):
        completed_statuses = {'COMPLETED', 'FAILED'}
        while flow.status not in completed_statuses:
            self.tick_flow(flow=flow, task_context=task_context)
        if check and flow.status == 'FAILED':
            raise self.FlowError(flow=flow)

    def tick_flow(self, flow=None, task_context=None):
        try:
            if flow.status == 'PENDING': self.start_flow(flow=flow)
            self.start_nearest_pending_tasks(flow=flow)
            self.tick_running_tasks(flow=flow, task_context=task_context)
            if not flow.has_incomplete_tasks(): self.complete_flow(flow=flow)
        except Exception as exception:
            fail_flow = True
            self.append_flow_error(flow=flow, error=traceback.format_exc())
            if isinstance(exception, self.TaskError):
                if not flow.cfg.get('fail_fast', True): fail_flow = False
            if fail_flow: self.fail_flow(flow=flow)

    def tick_flow_until_has_no_pending(self, flow=None, task_context=None):
        self.tick_flow(flow=flow, task_context=task_context)
        while (flow.status in {'PENDING', 'RUNNING'}
               and len(flow.get_nearest_pending_tasks()) > 0):
            self.tick_flow(flow=flow, task_context=task_context)

    def start_flow(self, flow=None):
        flow.status = 'RUNNING'

    def start_nearest_pending_tasks(self, flow=None):
        for task in flow.get_nearest_pending_tasks():
            try: self.start_task(flow=flow, task=task)
            except: self.fail_task(task=task, error=traceback.format_exc())

    def start_task(self, flow=None, task=None):
        task['status'] = 'RUNNING'

    def tick_running_tasks(self, flow=None, task_context=None):
        for task in flow.get_tasks_by_status(status='RUNNING'):
            if self.task_is_running(task=task):
                self.tick_task(task=task, flow=flow, task_context=task_context)
            else: self.complete_task(task=task)

    def task_is_running(self, task=None):
        return task['status'] == 'RUNNING'

    def tick_task(self, task=None, flow=None, task_context=None):
        task_context = task_context or {}
        try:
            self.task_handler.tick_task(
                task=task,
                task_context={**task_context, 'task': task, 'flow': flow}
            )
            if task.get('status') == 'COMPLETED':
                self.complete_task(flow=flow, task=task)
        except Exception as exception:
            self.fail_task(task=task, error=traceback.format_exc())

    def fail_task(self, task=None, error=None):
        task['error'] = error
        task['status'] = 'FAILED'
        msg = "Task with key '{key}' failed, error: {error}".format(
            key=task.get('key', '<unknown key>'),
            error=error
        )
        raise self.TaskError(msg)

    def fail_flow(self, flow=None, error=None):
        if error: self.append_flow_error(flow=flow, error=error)
        flow.status = 'FAILED'

    def append_flow_error(self, flow=None, error=None):
        flow.data.setdefault('errors', [])
        flow.data['errors'].append(self.elide_text(error))

    def elide_text(self, text=None, max_len=None):
        max_len = max_len or self.max_msg_len
        if max_len is not None and len(text) > max_len:
            text = text[0:max_len] + '...'
        return text

    def complete_task(self, flow=None, task=None):
        task['status'] = 'COMPLETED'

    def complete_flow(self, flow=None):
        if flow.data.get('errors'): self.fail_flow(flow=flow)
        else: flow.status = 'COMPLETED'
