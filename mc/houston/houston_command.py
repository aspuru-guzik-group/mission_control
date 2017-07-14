#!/usr/bin/env python
from mc.utils import import_utils
from mc.utils.commands.subcommand_command import SubcommandCommand

from .houston import Houston


class HoustonCommand(Houston, SubcommandCommand):
    class CfgError(Exception): pass

    def __init__(self, *args, default_cfg_path=None, **kwargs):
        Houston.__init__(self, *args, **kwargs)
        self.default_cfg_path = default_cfg_path
        SubcommandCommand.__init__(self, *args, **kwargs)
        self.subcommands = list(self.subcommands_registry.keys())
        
    def add_arguments(self, parser=None):
        parser.add_argument('-c', '--cfg_path', dest='cfg_path',
                            help="path to houston cfg file",
                            default=self.default_cfg_path)
        super().add_arguments(parser=parser)

    def _get_raw_cfg(self):
        cfg_path = self.kwargs['cfg_path']
        cfg_module = import_utils.load_module_from_path(path=cfg_path)
        return cfg_module.__dict__

if __name__ == '__main__': HoustonCommand.run()
