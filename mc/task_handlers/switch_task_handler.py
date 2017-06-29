import operator

from .base_proxying_task_handler import BaseProxyingTaskHandler


class SwitchTaskHandler(BaseProxyingTaskHandler):
    """Task handler that performs a conditional switch.

    task_params:
        control_value: value to use when checking cases
        cases: list of case dicts, that look like:
            .. code-block:: python

            {'condition': condition_dict,
             'task': task to run if condition matches}

            where a condition_dict looks like:
            .. code-block:: python

            {'op': <operator from the python 'op' library  e.g 'eq'>,
             'arg': 2nd arg to operator}
        default_case: the case to use if no case matches.
    """

    def validate_task_params(self):
        assert self.task['task_params']['control_value']
        assert bool(
            (self.task['task_params'].get('cases', []) or
             self.task['task_params']['default_case'])
        ) is True

    class NoMatchingCaseError(Exception):
        def __init__(self, *args, **kwargs):
            error = "No case matched and no default case was provided."
            super().__init__(self, error, *args, **kwargs)

    def generate_proxied_task(self, *args, **kwargs):
        matching_case = self.get_matching_case()
        return matching_case['task']

    def get_matching_case(self):
        control_value = self.get_control_value()
        for case in self.get_cases():
            if self.evaluate_case(case=case, control_value=control_value):
                return case
        try: return self.get_default_case()
        except: raise self.NoMatchingCaseError()

    def get_control_value(self):
        return self.task['task_params']['control_value']

    def get_cases(self):
        return self.task['task_params'].get('cases', [])

    def evaluate_case(self, case=None, control_value=None):
        condition = case['condition']
        op = getattr(operator, condition['op'])
        return op(control_value, condition['arg'])

    def get_default_case(self):
        return self.task['task_params']['default_case']

    def on_proxied_task_finished(self):
        super().on_proxied_task_finished()
        self.task['data'] = self.proxied_task.get('data')

TaskHandler = SwitchTaskHandler
