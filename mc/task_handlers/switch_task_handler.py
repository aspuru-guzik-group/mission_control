import operator

from .base_proxying_task_handler import BaseProxyingTaskHandler


class SwitchTaskHandler(BaseProxyingTaskHandler):
    class NoMatchingCaseError(Exception):
        def __init__(self, *args, **kwargs):
            error = "No case matched and no default case was provided."
            super().__init__(self, error, *args, **kwargs)

    def get_proxied_task(self, *args, **kwargs):
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
        try: return self.task['task_params']['control_value']
        except: raise self.InvalidTaskParamsError()

    def get_cases(self):
        try: return self.task['task_params'].get('cases', {})
        except: raise self.InvalidTaskParamsError()

    def evaluate_case(self, case=None, control_value=None):
        op = getattr(operator, case['op'])
        return op(control_value, case['arg'])

    def get_default_case(self):
        return self.task['task_params']['default_case']

    def test_case(self, case=None, control_value=None):
        pass

    def on_proxied_task_finished(self):
        super().on_proxied_task_finished()
        self.task['data'] = self.proxied_task.get('data')

TaskHandler = SwitchTaskHandler
