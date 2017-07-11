import logging

from mc.utils.commands.command_parser import CommandParser


class BaseHoustonSubcommand(object):
    def __init__(self, logger=None, args=None, kwargs=None, unparsed_args=None,
                 load_cfg=None, utils=..., **kwargs_):
        self.logger = logger or logging.getLogger(__name__)
        self.args = args
        self.kwargs = kwargs
        self.unparsed_args = unparsed_args
        subcommand_kwargs, self.unparsed_args = \
                self.parse_subcommand_kwargs(subcommand_argv=unparsed_args)
        self.kwargs = {**self.kwargs, **(vars(subcommand_kwargs))}
        self.load_cfg = load_cfg
        self.utils = utils

    @property
    def utils(self):
        if self._utils is ...: self._utils = self._get_utils()
        return self._utils

    @utils.setter
    def utils(self, new_value): self._utils = new_value

    def _get_utils(self):
        from . import _utils
        return _utils.HoustonSubcommandUtils(load_cfg=self.load_cfg)

    def parse_subcommand_kwargs(self, subcommand_argv=None):
        parser = CommandParser()
        self.add_arguments(parser=parser)
        return parser.parse_known_args(subcommand_argv)

    def add_arguments(self, parser=None): pass

    @classmethod
    def run(cls, args=None, kwargs=None, unparsed_args=None, load_cfg=None,
            **kwargs_):
        instance = cls(args=args, kwargs=kwargs, unparsed_args=unparsed_args,
                       load_cfg=load_cfg, **kwargs_)
        return instance._run()

    def _run(self): raise NotImplementedError


