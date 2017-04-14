import jinja2

class ContextValueSetter(object):
    def set_context_values(self, value_specs=None, context=None):
        for value_spec in (value_specs or []):
            self.set_context_value(value_spec=value_spec, context=context)

    def set_context_value(self, value_spec=None, context=None):
        try:
            self.set_context_dest_value(
                context=context,
                dest=value_spec['dest'],
                value=self.get_value_for_value_spec(
                    value_spec=value_spec,
                    context=context
                )
            )
        except Exception as exception:
            raise Exception("Unable to set_context_value:\n"
                            "  exception: {exception}\n"
                            "  value_spec: '{value_spec}'\n"
                            "  context: '{context}'".format(
                                exception=exception,
                                value_spec=value_spec,
                                context=self.elide_string(str(context))
                            ))

    def elide_string(self, string=None, max_len=100):
        elided_string = string
        if len(string) > max_len:
            elided_string = string[:max_len] + '...'
        return elided_string

    def get_value_for_value_spec(self, value_spec=None, context=None):
        if value_spec.get('as_raw_value'):
            value = value_spec['value']
        elif value_spec.get('as_template'):
            value = self.render_template(template=value_spec['value'],
                                         context=context)
        else: 
            value = self.get_value_from_dot_spec(obj={'ctx': context},
                                                 dot_spec=value_spec['value'])
        return value

    def render_template(self, template=None, context=None):
        return jinja2.Template(template, 
                               undefined=jinja2.StrictUndefined
                              ).render(ctx=context)

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
        pos = 0
        try:
            for path_element in path_elements[:-1]:
                next_cursor = self.get_attr_or_item(obj=cursor,
                                                    key=path_element)
                if next_cursor is None:
                    cursor[path_element] = {}
                    next_cursor = cursor[path_element]
                cursor = next_cursor
                pos += 1
            cursor[path_elements[-1]] = value
        except Exception as exception:
            raise Exception(("Unable to set value\n"
                             "  exception: {exception}\n"
                             "  error at path_element {pos} ({partial_path})"
                            ).format(exception=exception,
                                     pos=pos,
                                     partial_path=path_elements[:pos+1]))

    def get_value_from_dot_spec(self, obj=None, dot_spec=None):
        path_elements = dot_spec.split('.')
        cursor = obj
        pos = 0
        try:
            for path_element in path_elements:
                next_cursor = self.get_attr_or_item(obj=cursor, key=path_element)
                cursor = next_cursor
                pos += 1
            return cursor
        except Exception as exception:
            raise Exception(("Unable to get value\n"
                             "  exception: {exception}\n"
                             "  error at path_element {pos} ({partial_path})"
                            ).format(exception=exception,
                                     pos=pos,
                                     partial_path=path_elements[:pos+1]))

def set_context_values(value_specs=None, context=None):
    ContextValueSetter().set_context_values(value_specs=value_specs,
                                            context=context)