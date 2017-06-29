import collections
import logging
import pprint
import textwrap

from .base_task_handler import BaseTaskHandler


class LogTaskHandler(BaseTaskHandler):
    """Writes a msg to a python logger.

    Can also dump task_ctx info, by setting dump_task_ctx, dump_flow in
    task_params.
    
    """
    def initial_tick(self):
        task_params = self.task.get('task_params', {})
        log_level_name = task_params.get('log_level', 'WARNING')
        log_level = getattr(logging, log_level_name)
        use_print = task_params.get('use_print')

        def _log(msg='', obj=None, label=None, pprint_kwargs=None, **kwargs):
            if obj is not None:
                obj_dump = self.dump_obj(obj=obj, pprint_kwargs=pprint_kwargs)
                msg += "\n" + textwrap.indent(obj_dump, '  ')
            if label is not None: msg = label + ':' + msg
            if use_print: print(msg)
            else: self.get_logger().log(log_level, msg)

        if 'msg' in task_params: _log(task_params['msg'])

        dump_task_params = task_params.get('dump_task')
        if dump_task_params:
            extra_log_params = {}
            if isinstance(dump_task_params, collections.abc.Mapping):
                extra_log_params = dump_task_params
            _log(label='TASK', obj=self.task, **extra_log_params)

        dump_task_ctx_params = task_params.get('dump_task_ctx')
        if dump_task_ctx_params: 
            extra_log_params = {}
            if isinstance(dump_task_ctx_params, collections.abc.Mapping):
                extra_log_params = dump_task_ctx_params
            _log(label='TASK_CONTEXT', obj=self.task_ctx, **extra_log_params)

        dump_flow_params = task_params.get('dump_flow')
        if dump_flow_params:
            flow = self.task_ctx.get('flow')
            if flow: flow_dict = flow.to_dict()
            else: flow_dict = None
            extra_log_params = {}
            if isinstance(dump_flow_params, collections.abc.Mapping):
                extra_log_params = dump_flow_params
            _log(obj=flow_dict, label='FLOW', **extra_log_params)

        self.task['status'] = 'COMPLETED'

    def dump_obj(self, obj=None, pprint_kwargs=None):
        return pprint.pformat(obj, **(pprint_kwargs or {}))

    def get_logger(self): return logging

TaskHandler = LogTaskHandler
