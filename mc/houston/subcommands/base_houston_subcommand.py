import logging


class BaseHoustonSubcommand(object):
    def __init__(self, logger=None, args=None, kwargs=None, unparsed_args=None,
                 get_cfg=None, utils=..., **kwargs_):
        self.logger = logger or logging.getLogger(__name__)
        self.args = args
        self.kwargs = kwargs
        self.unparsed_args = unparsed_args
        self.get_cfg = get_cfg
        self.utils = utils

    @property
    def utils(self):
        if self._utils is ...: self._utils = self._get_utils()
        return self._utils

    @utils.setter
    def utils(self, new_value): self._utils = new_value

    def _get_utils(self):
        from . import _utils
        return _utils.HoustonSubcommandUtils(get_cfg=self.get_cfg)

    @classmethod
    def run(cls, args=None, kwargs=None, unparsed_args=None, get_cfg=None,
            **kwargs_):
        instance = cls(args=args, kwargs=kwargs, unparsed_args=unparsed_args,
                       get_cfg=get_cfg, **kwargs_)
        return instance._run()

    def _run(self): raise NotImplementedError


