import argparse

from mc.utils.commands.subcommand_command import SubcommandCommand
from mc.utils import import_utils
from .dispatcher import JobModuleCommandDispatcher


class JobModuleCliCommand(SubcommandCommand):
    help = 'JobModule CLI'

    subcommands = ['run_job_submission']

    def run_job_submission(self, args=None, kwargs=None, unparsed_args=None):
        run_job_submission_kwargs = self._parse_run_job_submission_args(
            unparsed_args=unparsed_args)
        if run_job_submission_kwargs.get('cfg_file_path'):
            submission_runner_cfg = import_utils.load_module_from_path(
                path=run_job_submission_kwargs['cfg_file_path'])
        else: submission_runner_cfg = ...
        JobModuleCommandDispatcher().run_job_submission(
            submission_dir=run_job_submission_kwargs['submission_dir'],
            submission_runner_cfg=submission_runner_cfg
        )

    def _parse_run_job_submission_args(self, unparsed_args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument('--submission_dir', required=True)
        parser.add_argument('--cfg_file_path')
        return vars(parser.parse_args(args=unparsed_args))

if __name__ == '__main__': JobModuleCliCommand.run()
