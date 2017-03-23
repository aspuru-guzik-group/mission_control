import jinja2

from .base_task_handler import BaseTaskHandler


class SetValuesTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        self.set_values(value_specs=task['params']['value_specs'],
                        task_context=task_context)
        task['status'] = 'COMPLETED'

    def set_values(self, value_specs=None, task_context=None):
        for value_spec in value_specs:
            self.set_context_value(value_spec=value_spec, context=task_context)

    def set_context_value(self, value_spec=None, context=None):
        self.set_context_dest_value(
            context=context,
            dest=value_spec['dest'],
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
        return jinja2.Template(template).render(ctx=context)

    def set_context_dest_value(self, context=None, dest=None, value=None):
        self.set_value_from_dot_spec(obj={'ctx': context}, dot_spec=dest,
                                     value=value)

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

