class TaskRunner(object):
    def __init__(self, task_handlers=None):
        self.task_handlers = (task_handlers or {})

    def tick_task(self, task=None, job=None):
        task_handler = self.get_task_handler(task=task, job=job)
        task_handler.tick_task(task=task, job=job)

    def get_task_handler(self, task=None, job=None):
        try: return self.task_handlers[task['type']]
        except KeyError:
            raise Exception(
                "No handler found for task with type '{}'".format(task['type']))
