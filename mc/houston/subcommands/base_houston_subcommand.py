import logging


class BaseHoustonSubcommand(object):
    def __init__(self, logger=None, args=None, kwargs=None, unparsed_args=None,
                 get_cfg=None, **kwargs_):
        self.logger = logger or logging.getLogger(__name__)
        self.args = args
        self.kwargs = kwargs
        self.unparsed_args = unparsed_args
        self.get_cfg = get_cfg

        self._cfg = ...

    @property
    def cfg(self):
        if self._cfg is ...: self._cfg = self.get_cfg()
        return self._cfg

    @classmethod
    def run(cls, args=None, kwargs=None, unparsed_args=None, get_cfg=None,
            **kwargs_):
        instance = cls(args=args, kwargs=kwargs, unparsed_args=unparsed_args,
                       get_cfg=get_cfg, **kwargs_)
        return instance._run()

    def _run(self): raise NotImplementedError


