class TaskRegistry(object):
    def __init__(self):
        self.tasks = {}

    def register_task(self, task, key=None):
        key = key or task
        self.tasks[key] = task

    def unregister_task(self, task=None, key=None):
        key = key or task
        del self.tasks[key]

    def clear(self):
        self.tasks = {}

task_registry = TaskRegistry()

