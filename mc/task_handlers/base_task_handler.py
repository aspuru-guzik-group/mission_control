class BaseTaskHandler(object):
    def __init__(self, *args, **kwargs):
        pass

    def tick_task(self, task=None, job=None):
        try:
            self._ensure_task(task=task)
            self.increment_tick_counter()
            if task['data']['tick_counter'] == 1:
                task['status'] = 'RUNNING'
                self.initial_tick(task=task, job=job)
            else: self.intermediate_tick(task=task, job=job)
        except Exception as error:
            task['status'] = 'FAILED'
            task['data']['error'] = error

    def _ensure_task(self, task=None):
        task.setdefault('data', {})
        task['data'].setdefault('_tick_counter', 0)

    def increment_tick_counter(self, task=None):
        task['data']['_tick_counter'] += 1

    def initial_tick(self, task=None, job=None): raise NotImplementedError

    def intermediate_tick(self, task=None, job=None): raise NotImplementedError
