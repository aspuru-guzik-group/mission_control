import argparse
import os

from mc.utils.commands.subcommand_command import SubcommandCommand
from mc.job_module_utils.job_module_command_dispatcher import (
    JobModuleCommandDispatcher)


class JobModuleCliCommand(SubcommandCommand):
    help = 'JobModule CLI'

    subcommands = ['run_jobdir']

    def run_jobdir(self, parsed_args=None, unparsed_args=None):
        run_jobdir_kwargs = self._parse_run_jobdir_args(
            unparsed_args=unparsed_args)
        command_dispatcher = self._get_command_dispatcher()
        command_dispatcher.run_jobdir(**run_jobdir_kwargs)

    def _get_command_dispatcher(self): return JobModuleCommandDispatcher()

    def _parse_run_jobdir_args(self, unparsed_args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument('--jobdir', default=os.environ['PWD'])
        return vars(parser.parse_args(args=unparsed_args))

if __name__ == '__main__': JobModuleCliCommand.run()
