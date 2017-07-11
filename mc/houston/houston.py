#!/usr/bin/env python

import logging


from mc.utils.commands.subcommand_command import SubcommandCommand

from . import houston_cfg


class HoustonCommand(SubcommandCommand):
    class CfgError(Exception): pass

    def __init__(self, *args, logger=None, subcommands_registry=..., **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logger or self._get_default_logger()
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
        from .subcommands.registry import SubcommandRegistry
        return SubcommandRegistry()

    def _get_cfg(self):
        if not hasattr(self, '_houston_cfg'):
            self._houston_cfg = self._load_houston_cfg()
        return self._houston_cfg

    def _load_houston_cfg(self):
        try: return houston_cfg.HoustonCfg()
        except Exception as exc: raise self.CfgError() from exc

    # Override
    def _get_subcommand_fn(self, subcommand=None):
        return self.subcommands_registry[subcommand]

    # Override
    def _call_subcommand_fn(self, subcommand_fn=None, args=None, kwargs=None,
                            unparsed_args=None):
        return subcommand_fn(args=args, kwargs=kwargs,
                             unparsed_args=unparsed_args,
                             get_cfg=self._get_cfg)

if __name__ == '__main__': HoustonCommand.run()
