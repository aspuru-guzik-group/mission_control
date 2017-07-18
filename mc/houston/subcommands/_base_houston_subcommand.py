import logging

from mc.utils.commands.argument_parser import ArgumentParser


class BaseHoustonSubcommand(object):
    def __init__(self, logger=None, parsed_args=None, unparsed_args=None,
                 get_cfg=None, utils=..., **kwargs_):
        self.logger = logger or logging.getLogger(__name__)
        self.parsed_args = parsed_args
        self.unparsed_args = unparsed_args
        subcommand_parsed_args, self.unparsed_args = \
                self.parse_subcommand_args(subcommand_argv=unparsed_args)
        self.parsed_args = {**self.parsed_args,
                            **(vars(subcommand_parsed_args))}
        self.get_cfg = get_cfg
        self.utils = utils

    @property
    def utils(self):
        if self._utils is ...: self._utils = self._get_utils()
        return self._utils

    @utils.setter
    def utils(self, new_value): self._utils = new_value

    def _get_utils(self):
        from mc.houston import utils
        return utils.HoustonUtils(get_cfg=self.get_cfg)

    def parse_subcommand_args(self, subcommand_argv=None):
        parser = ArgumentParser()
        self.add_arguments(parser=parser)
        return parser.parse_known_args(subcommand_argv)

    def add_arguments(self, parser=None): pass

    @classmethod
    def run(cls, parsed_args=None, unparsed_args=None, get_cfg=None,
            logger=None, **kwargs_):
        instance = cls(parsed_args=parsed_args, unparsed_args=unparsed_args,
                       get_cfg=get_cfg, logger=logger, **kwargs_)
        return instance._run()

    def _run(self): raise NotImplementedError


