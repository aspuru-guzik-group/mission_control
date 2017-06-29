import abc
from .base_task_handler import BaseTaskHandler


class BaseProxyingTaskHandler(BaseTaskHandler):
    """
    Abstract base class for proxied tasks.

    Implementing classes must define :code:`generate_proxied_task`
    and may define :code:`on_proxied_task_finished`.

    """
    @property
    def proxied_task(self): return self.task.get('proxied_task')

    @proxied_task.setter
    def proxied_task(self, value): self.task['proxied_task'] = value

    def initial_tick(self, *args, **kwargs):
        self.proxied_task = self.generate_proxied_task()

    def intermediate_tick(self, *args, **kwargs):
        self.on_proxied_task_finished()

    @abc.abstractmethod
    def generate_proxied_task(self, *args, **kwargs):
        raise NotImplementedError()

    def on_proxied_task_finished(self, *args, **kwargs):
        self.propagate_proxied_task_status()

    def propagate_proxied_task_status(self):
        self.task['status'] = self.proxied_task.get('status')
