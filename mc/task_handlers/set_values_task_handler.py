from uuid import uuid4

import jinja2

from .base_task_handler import BaseTaskHandler


class SetValuesTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, job=None):
        self.set_values(value_specs=task['params']['value_specs'],
                        job=job)
        task['status'] = 'COMPLETED'

    def set_values(self, value_specs=None, job=None):
        context = self.get_context(job=job)
        for value_spec in value_specs:
            self.set_context_value(value_spec=value_spec, context=context)

    def get_context(self, job=None):
        context = {
            'job': job,
            'keyed_tasks': self.get_keyed_tasks(job=job),
        }
        return context

    def get_keyed_tasks(self, job=None):
        keyed_tasks = {}
        for task in job.get('tasks', []):
            key = task.get('key', str(uuid4()))
            keyed_tasks[key] = task
        return keyed_tasks

    def set_context_value(self, value_spec=None, context=None):
        self.set_context_target_value(
            context=context,
            target=value_spec['target'],
            value=self.get_value_for_value_spec(
                value_spec=value_spec,
                context=context
            )
        )

    def get_value_for_value_spec(self, value_spec=None, context=None):
        if value_spec.get('is_raw_value', False): value = value_spec['value']
        else: value = self.render_template(template=value_spec['value'],
                                           context=context)
        return value

    def render_template(self, template=None, context=None):
        return jinja2.Template(template).render(**context)

    def set_context_target_value(self, context=None, target=None, value=None):
        self.set_value_from_dot_spec(obj=context, dot_spec=target, value=value)

    def get_attr_or_item(self, obj=None, key=None):
        if hasattr(obj, key): return getattr(obj, key)
        elif key in obj: return obj[key]
        else: return None

    def set_value_from_dot_spec(self, obj=None, value=None, dot_spec=None):
        path_elements = dot_spec.split('.')
        cursor = obj
        for path_element in path_elements[:-1]:
            next_cursor = self.get_attr_or_item(
                obj=cursor, key=path_element)
            if next_cursor is None:
                cursor[path_element] = {}
                next_cursor = cursor[path_element]
            cursor = next_cursor
        cursor[path_elements[-1]] = value

