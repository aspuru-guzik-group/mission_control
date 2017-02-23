import jinja2


class CopyActionHandler(object):
    def __init__(self, transfer_client=None, render_tpl=None):
        self.transfer_client = transfer_client
        self.render_tpl = render_tpl or self._render_tpl

    def _render_tpl(self, tpl=None, ctx=None):
        return jinja2.Template(tpl).render(ctx=ctx)

    def __call__(self, params=None, ctx=None):
        formatted_params = self.format_params(params=params, ctx=ctx)
        return self.transfer_client.copy(src=formatted_params['src'],
                                         dest=formatted_params['dest'])

    def format_params(self, params=None, ctx=None):
        formatted_params = {key:self.render_tpl(tpl=params[key], ctx=ctx)
                            for key in ['src', 'dest']}
        return formatted_params
