from mc.task_handlers.base_task_handler import BaseTaskHandler


class ExampleCountdownTaskHandler(BaseTaskHandler):
    def initial_tick(self, *args, task_ctx=None, **kwargs):
        self.task['data']['countdown'] = \
                self.task['task_params'].get('countdown_start', 3)
        print("Starting countdown for task with key '%s'" % self.task['key'])
        self.intermediate_tick(*args, task_ctx=task_ctx, **kwargs)
    
    def intermediate_tick(self, *args, task_ctx=None, **kwargs):
        countdown = self.task['data']['countdown']
        if countdown > 0:
            print("T minus %s" % countdown)
            self.task['data']['countdown'] -= 1
        else:
            print("Houston, we have lift off.")
            self.task['data']['my_result'] = 'Lift off'
            self.task['status'] = 'COMPLETED'

TaskHandler = ExampleCountdownTaskHandler
