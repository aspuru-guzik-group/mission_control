import json
from pathlib import Path

from mc.utils.commands.argument_parser import ArgumentParser


class BaseSubcommand(object):
    class ENTRYPOINTS:
        RUN = 'RUN'
        CALL = 'CALL'

    def __init__(self, houston=None, logger=None, **kwargs):
        self.houston = houston
        self.logger = logger or self.houston.logger
        self.parsed_args = None
        self.unparsed_args = None

    @classmethod
    def run(cls, *args, houston=None, logger=None, **kwargs):
        self = cls(houston=houston, logger=logger)
        self.entrypoint = self.ENTRYPOINTS.RUN
        self.unparsed_args = args
        self.parsed_args = kwargs
        self._ensure_parsed_args()
        return self._run()

    def _ensure_parsed_args(self):
        for key, default_value in self._get_defaults().items():
            self.parsed_args.setdefault(key, default_value)

    def _get_defaults(self): return {}

    def _run(self): raise NotImplementedError

    @classmethod
    def call(cls, houston=None, logger=None, parsed_args=None,
             unparsed_args=None, **kwargs):
        self = cls(houston=houston, logger=logger, **kwargs)
        self.entrypoint = self.ENTRYPOINTS.CALL
        extra_parsed, self.unparsed_args = self.parse_argv(argv=unparsed_args)
        self.parsed_args = {**(parsed_args or {}), **(vars(extra_parsed))}
        self._ensure_parsed_args()
        return self._run()

    def parse_argv(self, argv=None):
        parser = ArgumentParser()
        self.add_arguments(parser=parser)
        return parser.parse_known_args(argv)

    def add_arguments(self, parser=None): pass

    @property
    def utils(self): return self.houston.utils

    @property
    def db(self): return self.utils.db

    @property
    def session(self): return self.db.session

    def json_path_arg(self, path):
        with Path(path).expanduser().open() as f:
            return json.load(f)
