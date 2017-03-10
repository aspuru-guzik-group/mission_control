from unittest.mock import call, Mock

from mc.a2g2.utils.tests.test_base_command import BaseCommandBaseTestCase

from ..job_module_base_command import JobModuleBaseCommand


class JobModuleBaseCommandBaseTestCase(BaseCommandBaseTestCase):
    def setUp(self):
        self.job = self.generate_job()
        self.ctx_dir = 'some_ctx_dir'
        self.output_dir = 'some_output_dir'
        super().setUp()
        self.file_args = self.generate_file_args(target_dir=self.tmpdir)
        self.argv = self.generate_argv(arg_tuples=[
            *[(file_arg_name, file_arg['path']) 
              for file_arg_name, file_arg in self.file_args.items()
             ],
            ('ctx_dir', self.ctx_dir),
            ('output_dir', self.output_dir),
        ])
        self.command = self.generate_command()

    def generate_job(self): return {}

    def generate_command(self):
        class BasicJobModuleCommand(JobModuleBaseCommand):
            def execute_job(self): pass
        return BasicJobModuleCommand()

    def generate_file_args(self, target_dir=None):
        file_args = super().generate_file_args(target_dir=target_dir)
        extra_file_args = {
            'job': {
                'value': self.job,
            }
        }
        for file_arg_name, file_arg in extra_file_args.items():
            file_arg['path'] = self.generate_json_file_for_arg(
                target_dir=target_dir,
                arg_name=file_arg_name,
                arg_value=file_arg['value']
            )
            file_args[file_arg_name] = file_arg
        return file_args

    def test_calls_execute_job(self):
        self.command.execute_job = Mock()
        self.execute_command()
        self.assertEqual(self.command.execute_job.call_args,
                         self.get_expected_execute_job_call())

    def get_expected_execute_job_call(self):
        return call(job=self.job, cfg=self.cfg,
                    ctx_dir=self.ctx_dir, output_dir=self.output_dir)
