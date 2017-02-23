import collections
import jinja2


class ActionProcessor(object):
    def __init__(self, *args, handlers=None, **kwargs):
        self.handlers = {}
        self.setup_handlers(handlers=handlers)

    def setup_handlers(self, handlers=None):
        for key, handler in (handlers or {}).items():
            self.register_handler(key=key, handler=handler)
        if 'set_ctx_value' not in  self.handlers:
            self.register_handler(key='set_ctx_value',
                                  handler=self.set_ctx_value_handler)

    def register_handler(self, key=None, handler=None):
        self.handlers[key] = handler

    def set_ctx_value_handler(self, params=None, ctx=None): 
        self.set_ctx_value(ctx=ctx, **params)

    def set_ctx_value(self, ctx=None, target=None, value=None):
        value = ctx.transform_value(value=value)
        ctx.set(target=target, value=value)

    def process_action(self, action=None, ctx=None):
        handler = self.get_handler_for_action(action=action)
        wrapped_ctx = self._wrap_ctx_for_action(ctx=ctx)
        output = handler(params=action.get('params', None), ctx=wrapped_ctx)
        output_target = action.get('output_to_ctx_target', None)
        if output_target: self.set_ctx_value(ctx=wrapped_ctx,
                                             target=output_target,
                                             value=output)

    def _wrap_ctx_for_action(self, ctx=None):
        return ActionCtx(ctx_to_wrap=ctx)

    def get_handler_for_action(self, action=None):
        return self.handlers[action['action']]

class ActionCtx(object):
    def __init__(self, ctx_to_wrap=None):
        self.ctx = ctx_to_wrap

    def get(self, target=None):
        path_elements = target.split('.')
        cursor = self.ctx
        for path_element in path_elements:
            cursor = self.get_attr_or_item(obj=cursor, key=path_element)
        return cursor

    def get_attr_or_item(self, obj=None, key=None):
        if hasattr(obj, key): return getattr(obj, key)
        elif key in obj: return obj[key]
        else: return None

    def set(self, target=None, value=None):
        path_elements = target.split('.')
        cursor = self.ctx
        for path_element in path_elements[:-1]:
            next_cursor = self.get_attr_or_item(
                obj=cursor, key=path_element)
            if next_cursor is None:
                cursor[path_element] = {}
                next_cursor = cursor[path_element]
            cursor = next_cursor
        cursor[path_elements[-1]] = value

    def transform_value(self, value=None):
        if self.is_template_value(value):
            value = self.render_template(value['template'])
        return value

    def render_template(self, template=None):
        return jinja2.Template(template).render(ctx=self.ctx)

    def is_template_value(self, value=None):
        return isinstance(value, collections.Mapping) and ('template' in value)
