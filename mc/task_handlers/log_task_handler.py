import logging
import pprint
import textwrap

from .base_task_handler import BaseTaskHandler


class LogTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        task_params = task.get('task_params', {})
        log_level_name = task_params.get('log_level', 'INFO')
        log_level = getattr(logging, log_level_name)
        def _log(msg, obj=None):
            if obj is not None:
                msg += "\n" + textwrap.indent(self.format_obj(obj), '  ')
            logging.log(log_level, msg)
        if 'msg' in task_params: _log(task_params['msg'])
        if task_params.get('dump_task'): _log("task", obj=task)
        if task_params.get('dump_task_context'):
            _log("task_context", obj=task_context)
        task['status'] = 'COMPLETED'

    def format_obj(self, obj=None):
        return pprint.PrettyPrinter(indent=2).pformat(obj)
