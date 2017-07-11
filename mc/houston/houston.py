#!/usr/bin/env python

import logging


from mc.utils.commands.subcommand_command import SubcommandCommand

from . import houston_cfg


class HoustonCommand(SubcommandCommand):
    class CfgError(Exception): pass

    def __init__(self, *args, logger=None, default_cfg_path=None,
                 subcommands_registry=..., **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logger or self._get_default_logger()
        self.default_cfg_path = default_cfg_path
        if subcommands_registry is ...:
            subcommands_registry = self._get_default_subcommands_registry()
        self.subcommands_registry = subcommands_registry
        self.subcommands = list(self.subcommands_registry.keys())
        
    def _get_default_logger(self):
        logger = logging.getLogger(__name__)
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)
        return logger

    def _get_default_subcommands_registry(self):
        from .subcommands._registry import SubcommandRegistry
        return SubcommandRegistry()

    def add_arguments(self, parser=None):
        parser.add_argument('-c', '--cfg_path', dest='cfg_path',
                            help="path to houston cfg file",
                            default=self.default_cfg_path)
        super().add_arguments(parser=parser)

    # Override
    def _get_subcommand_fn(self, subcommand=None):
        return self.subcommands_registry[subcommand]

    # Override
    def _call_subcommand_fn(self, subcommand_fn=None, args=None, kwargs=None,
                            unparsed_args=None):
        return subcommand_fn(
            logger=self.logger,
            args=args, kwargs=kwargs,
            unparsed_args=unparsed_args,
            load_cfg=self._get_load_cfg_fn(kwargs['cfg_path'])
        )

    def _get_load_cfg_fn(self, cfg_path=None):
        def load_cfg():
            try: return houston_cfg.HoustonCfg(cfg_path=cfg_path)
            except Exception as exc: raise self.CfgError() from exc
        return load_cfg

if __name__ == '__main__': HoustonCommand.run()
