import importlib
import os


class SubcommandRegistry(object):
    def __init__(self, setup=True):
        self._setup()

    def _setup(self):
        self.subcommands_dir = os.path.dirname(__file__)
        self._subcommands = {
            subcommand: ...
            for subcommand in self._get_subcommand_names()
        }

    def keys(self): return self._subcommands.keys()

    def __iter__(self): return iter(self.keys())

    def __getitem__(self, key):
        if self._subcommands[key] is ...:
            self._subcommands[key] = self._get_subcommand_fn(subcommand=key)
        return self._subcommands[key]

    def _get_subcommand_names(self):
        subcommands = [
            os.path.splitext(file_name)[0]
            for file_name in os.listdir(self.subcommands_dir)
            if not file_name.startswith('_') and file_name.endswith('.py')
        ]
        return subcommands

    def _get_subcommand_fn(self, subcommand=None):
        subcommand_module = self._load_subcommand_module(subcommand=subcommand)
        subcommand_cls = getattr(subcommand_module, 'Subcommand')
        return subcommand_cls

    def _load_subcommand_module(self, subcommand=None):
        return importlib.import_module('.' + subcommand, __package__)
