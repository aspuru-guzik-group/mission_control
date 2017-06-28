import logging
import traceback

from mc.utils import debug_utils
from .flow import Flow


class FlowEngine(object):
    """
    A FlowEngine is what ticks a flow.
    (https://en.wikipedia.org/wiki/Visitor_pattern)
    """
    class FlowError(Exception):
        def __init__(self, *args, flow=None, **kwargs):
            self.flow = flow
            msg = "\n".join([str(err) for err in flow.data.get('errors', [])])
            super().__init__(msg, *args, **kwargs)

    class TaskError(Exception): pass

    def __init__(self, task_handler=None, logger=None, max_msg_len=None,
                 debug=False):
        """
        Args:
            logger [logging.Logger]: logger instance
            debug [bool]: whether to log debug statements. Default: False
            task_handler [task_handler]: task_handler instance.
                Default: get_default_task_handler return.
            max_msg_len [int]: maximum length of error messages. Default:
                no max length.
        """
        self.logger = logger or self.get_default_logger()
        self.debug = debug
        self.task_handler = task_handler or self.get_default_task_handler()
        self.max_msg_len = max_msg_len

    def debug_locals(self):
        if self.debug:
            debug_utils.debug_locals(frames_back=1, logger=self.logger)

    def get_default_logger(self):
        logger = logging.getLogger(name=str(self))
        logger.addHandler(logging.StreamHandler())
        return logger

    def get_default_task_handler(self):
        """
        Uses mc.task_handlers.mc_default_task_handler.McDefaultTaskHandler
        as the default task handler.
        """
        from mc.task_handlers.mc_default_task_handler import \
                McDefaultTaskHandler
        return McDefaultTaskHandler()

    @classmethod
    def flow_spec_to_flow(cls, flow_spec=None, **kwargs):
        return Flow.from_flow_spec(flow_spec=flow_spec)

    @classmethod
    def flow_dict_to_flow(cls, flow_dict=None, **kwargs):
        return Flow.from_flow_dict(flow_dict=flow_dict)

    @classmethod
    def flow_to_flow_dict(cls, flow=None): return flow.to_flow_dict()

    @classmethod
    def flow_spec_to_flow_dict(cls, flow_spec=None, **kwargs):
        flow = Flow.from_flow_spec(flow_spec=flow_spec)
        return cls.flow_to_flow_dict(flow=flow)

    def run_flow(self, flow=None, task_ctx=None, max_ticks=5e2):
        """
        Tick flow until its status is 'COMPLETED' or 'FAILED'.
        If it failed, raise a FlowError.
        
        Args:
            flow <flow>: flow to tick
            task_ctx [dict]: task_ctx to include when ticking flow
            max_ticks [int]: Fail if hit more than this number of ticks.
        """
        self.debug_locals()
        completed_statuses = {'COMPLETED', 'FAILED'}
        tick_counter = 0
        while flow.status not in completed_statuses:
            tick_counter += 1
            if self.debug: self.logger.debug("tick #{}".format(tick_counter))
            self.tick_flow(flow=flow, task_ctx=task_ctx)
            if tick_counter > max_ticks: raise Exception("Exceeed max ticks")
        if flow.status == 'FAILED': raise self.FlowError(flow=flow)
        if self.debug: self.logger.debug("complete")

    def tick_flow(self, flow=None, task_ctx=None):
        """Tick a flow.

        Args:
            flow <flow>: flow to tick
            task_ctx [dict]: task_ctx to include when ticking flow
        """
        self.debug_locals()
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
        """Tick a flow until it has no tickable pending tasks.

        Args:
            flow <flow>: flow to tick
            task_ctx [dict]: task_ctx to include when ticking flow
        """
        self.tick_flow(flow=flow, task_ctx=task_ctx)
        while (flow.status in {'PENDING', 'RUNNING'}
               and len(flow.get_nearest_tickable_pending_tasks()) > 0):
            self.tick_flow(flow=flow, task_ctx=task_ctx)

    def start_flow(self, flow=None):
        self.debug_locals()
        flow.status = 'RUNNING'

    def start_nearest_tickable_pending_tasks(self, flow=None):
        for task in flow.get_nearest_tickable_pending_tasks():
            try: self.start_task(flow=flow, task=task)
            except: self.fail_task(task=task, error=traceback.format_exc())

    def start_task(self, flow=None, task=None):
        self.debug_locals()
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
        """Tick a task.

        If task is a proxying task, we will tick the proxied task.

        Args:
            task <task>: task to tick
            flow <flow>: flow to tick
            task_ctx [dict]: task_ctx to include when ticking flow
        """
        self.debug_locals()
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
