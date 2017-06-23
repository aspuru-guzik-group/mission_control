import logging
import traceback

from .flow import Flow


class FlowEngine(object):

    class FlowError(Exception):
        def __init__(self, *args, flow=None, **kwargs):
            self.flow = flow
            err = "\n".join([str(err) for err in flow.data.get('errors', [])])
            super().__init__(err, *args, **kwargs)

    class TaskError(Exception): pass

    def __init__(self, task_handler=None, logger=None, max_msg_len=None):
        self.logger = logger or logging
        self.task_handler = task_handler or self.get_default_task_handler()
        self.max_msg_len = max_msg_len

    def get_default_task_handler(self):
        from mc.task_handlers.mc_default_task_handler import \
                McDefaultTaskHandler
        return McDefaultTaskHandler()

    @classmethod
    def flow_spec_to_flow(cls, flow_spec=None, **kwargs):
        flow = Flow(**Flow.sanitize_flow_kwargs(flow_spec))
        for i, task in enumerate(flow_spec.get('tasks', [])):
            if 'precursors' not in task and 'sucessors' not in task:
                if i == 0: precursor = 'ROOT'
                else: precursor = flow_spec['tasks'][i - 1]['key']
                task['precursors'] = [precursor]
            flow.add_task(task=task)
        return flow

    @classmethod
    def flow_dict_to_flow(cls, flow_dict=None, **kwargs):
        return Flow.from_flow_dict(flow_dict=flow_dict)

    @classmethod
    def flow_to_flow_dict(cls, flow=None): return flow.to_flow_dict()

    def run_flow(self, flow=None, task_ctx=None, check=False):
        completed_statuses = {'COMPLETED', 'FAILED'}
        while flow.status not in completed_statuses:
            self.tick_flow(flow=flow, task_ctx=task_ctx)
        if check and flow.status == 'FAILED':
            raise self.FlowError(flow=flow)

    def tick_flow(self, flow=None, task_ctx=None):
        try:
            flow.data.setdefault('_tick_counter', 0)
            flow.data['_tick_counter'] += 1
            if flow.status == 'PENDING': self.start_flow(flow=flow)
            self.start_nearest_tickable_pending_tasks(flow=flow)
            self.tick_running_tasks(flow=flow, task_ctx=task_ctx)
            if not flow.has_incomplete_tasks(): self.complete_flow(flow=flow)
        except Exception as exception:
            fail_flow = True
            self.append_flow_error(flow=flow, error=traceback.format_exc())
            if isinstance(exception, self.TaskError):
                if not flow.cfg.get('fail_fast', True): fail_flow = False
            if fail_flow: self.fail_flow(flow=flow)

    def tick_flow_until_has_no_pending(self, flow=None, task_ctx=None):
        self.tick_flow(flow=flow, task_ctx=task_ctx)
        while (flow.status in {'PENDING', 'RUNNING'}
               and len(flow.get_nearest_tickable_pending_tasks()) > 0):
            self.tick_flow(flow=flow, task_ctx=task_ctx)

    def start_flow(self, flow=None):
        flow.status = 'RUNNING'

    def start_nearest_tickable_pending_tasks(self, flow=None):
        for task in flow.get_nearest_tickable_pending_tasks():
            try: self.start_task(flow=flow, task=task)
            except: self.fail_task(task=task, error=traceback.format_exc())

    def start_task(self, flow=None, task=None):
        task['status'] = 'RUNNING'

    def tick_running_tasks(self, flow=None, task_ctx=None):
        for task in flow.get_tasks_by_status(status='RUNNING'):
            if self.task_is_running(task=task):
                self.tick_task(task=task, flow=flow, task_ctx=task_ctx)
            else: self.complete_task(task=task)
            if flow.status == 'COMPLETED': break

    def task_is_running(self, task=None):
        return task['status'] == 'RUNNING'

    def tick_task(self, task=None, flow=None, task_ctx=None):
        task_ctx = task_ctx or {}
        try:
            if self.is_proxying_task(task=task):
                self.tick_proxying_task(proxying_task=task, flow=flow,
                                        task_ctx=task_ctx)
            else:
                self.task_handler.tick_task(task_ctx={**task_ctx, 'task': task,
                                                      'flow': flow})
            if task.get('status') == 'COMPLETED':
                self.complete_task(flow=flow, task=task)
        except Exception as exception:
            self.fail_task(task=task, error=traceback.format_exc())

    def is_proxying_task(self, task=None):
        return (
            ('proxied_task' in task)
            and
            (task['proxied_task'].get('status') not in {'COMPLETED', 'FAILED'})
        )

    def tick_proxying_task(self, proxying_task=None, flow=None, task_ctx=None):
        proxied_task = proxying_task['proxied_task']
        self.tick_task(task=proxied_task, flow=flow, task_ctx=task_ctx)

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
