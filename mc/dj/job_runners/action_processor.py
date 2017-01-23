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
        ctx.set(target=target, value=value)

    def process_action(self, action=None, ctx=None):
        handler = self.handlers[action['type']]
        output = handler(params=action.get('params', None), ctx=ctx)
        output_target = action.get('output_to_ctx_target', None)
        if output_target: self.set_ctx_value(ctx=ctx, target=output_target,
                                             value=output)

