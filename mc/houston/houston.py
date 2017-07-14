import logging

from . import houston_cfg


class Houston(object):
    class CfgError(Exception): pass

    def __init__(self, *args, logger=None, cfg=None, subcommands_registry=...,
                 **kwargs):
        super().__init__()
        self.logger = logger or self._get_default_logger()
        self.cfg = cfg
        if subcommands_registry is ...:
            subcommands_registry = self._get_default_subcommands_registry()
        self.subcommands_registry = subcommands_registry
        
    def _get_default_logger(self):
        logger = logging.getLogger(__name__)
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)
        return logger

    def _get_default_subcommands_registry(self):
        from .subcommands._registry import SubcommandRegistry
        return SubcommandRegistry()

    def call_command(self, command, *args, **kwargs):
        self._call_subcommand_fn(
            subcommand_fn=self._get_subcommand_fn(subcommand=command),
            args=None, kwargs=None,
            unparsed_args=[*args, self._kwargs_to_unparsed_args(**kwargs)]
        )

    def _kwargs_to_unparsed_args(self, **kwargs):
        return ['--{k}={v}'.format(k=k, v=v) for k, v in kwargs.items()]

    # Override for mc.utils.comamnds.subcommand_command
    def _get_subcommand_fn(self, subcommand=None):
        return self.subcommands_registry[subcommand]

    # Override for mc.utils.comamnds.subcommand_command
    def _call_subcommand_fn(self, subcommand_fn=None, args=None, kwargs=None,
                            unparsed_args=None):
        return subcommand_fn(
            logger=self.logger,
            args=args, kwargs=kwargs,
            unparsed_args=unparsed_args,
            load_cfg=self._get_load_cfg_fn(
                subcommand_fn=None, args=None, kwargs=None, unparsed_args=None)
        )

    def _get_load_cfg_fn(self, subcommand_fn=None, args=None, kwargs=None,
                         unparsed_args=None):
        def load_cfg(): 
            try:
                raw_cfg = self._get_raw_cfg(args=args, kwargs=kwargs,
                                            unparsed_args=None)
                return houston_cfg.HoustonCfg(cfg=raw_cfg)
            except Exception as exc: raise self.CfgError() from exc
        return load_cfg

    def _get_raw_cfg(self, *args, **kwargs): return self.cfg
