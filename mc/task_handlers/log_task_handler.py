import logging
import pprint
import textwrap

from .base_task_handler import BaseTaskHandler


class LogTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        task_params = task.get('task_params', {})
        log_level_name = task_params.get('log_level', 'INFO')
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
        if dump_task_params: _log(label='TASK', obj=task, **dump_task_params)

        dump_task_context_params = task_params.get('dump_task_context')
        if dump_task_context_params:
            _log(label='TASK_CONTEXT', obj=task_context,
                 **dump_task_context_params)

        dump_flow_params = task_params.get('dump_flow')
        if dump_flow_params:
            flow = task_context.get('flow')
            if flow: flow_dict = flow.to_dict()
            else: flow_dict = None
            _log(obj=flow_dict, label='FLOW', **dump_flow_params)

        task['status'] = 'COMPLETED'

    def dump_obj(self, obj=None, pprint_kwargs=None):
        return pprint.pformat(obj, **(pprint_kwargs or {}))

    def get_logger(self): return logging

TaskHandler = LogTaskHandler
