from mc.utils.config import Config
from mc.utils import logging_utils


class Houston(object):
    class ConfigError(Exception):
        pass

    def __init__(self, *args, logger=None, cfg=None, ensure_db=True,
                 ensure_job_dirs=True, **kwargs):
        super().__init__()
        self.logger = logger or logging_utils.get_default_logger(name=__name__)
        if cfg:
            self.set_cfg(cfg=cfg)
        if ensure_db:
            self.utils.db.ensure_tables()
        if ensure_job_dirs:
            self.utils.ensure_job_dirs()

    def set_cfg(self, cfg=None):
        self.cfg = Config(cfg=cfg)

    @property
    def utils(self):
        if not hasattr(self, '_utils'):
            from .utils import HoustonUtils
            self._utils = HoustonUtils(houston=self)
        return self._utils

    @utils.setter
    def utils(self, new_value): self._utils = new_value

    @property
    def subcommands(self):
        if not hasattr(self, '_subcommands'):
            self._subcommands = self.get_default_subcommands()
        return self._subcommands

    @subcommands.setter
    def subcommands(self, value): self._subcommands = value

    @classmethod
    def get_default_subcommands(cls):
        from .subcommands._registry import SubcommandRegistry
        return SubcommandRegistry()

    def call_command(self, command, *args, **kwargs):
        subcommand = self._get_subcommand(subcommand=command)
        return self._call_subcommand(
            subcommand,
            parsed_args={},
            unparsed_args=[*args, *self._kwargs_to_unparsed_args(**kwargs)]
        )

    def _kwargs_to_unparsed_args(self, **kwargs):
        args = []
        for k, v in kwargs.items():
            args.extend(self._kvp_to_args(k=k, v=v))
        return args

    def _kvp_to_args(self, k=None, v=None):
        if isinstance(v, bool):
            args = ['--{k}'.format(k=k)]
        elif isinstance(v, list):
            args = []
            for item in v:
                args.extend(self._kvp_to_args(k=k, v=item))
        else:
            args = ['--{k}={v}'.format(k=k, v=v)]
        return args

    def _get_subcommand(self, subcommand=None):
        return self.subcommands[subcommand]

    def _call_subcommand(self, subcommand, parsed_args=None,
                         unparsed_args=None):
        return subcommand.call(houston=self, parsed_args=parsed_args,
                               unparsed_args=unparsed_args)

    def run_command(self, command, *args, **kwargs):
        subcommand = self._get_subcommand(subcommand=command)
        return self._run_subcommand(subcommand, *args, **kwargs)

    def _run_subcommand(self, subcommand=None, *args, **kwargs):
        return subcommand.run(*args, houston=self, **kwargs)

    @property
    def db(self): return self.utils.db  # shortcut
