import textwrap
import unittest
from unittest.mock import call, MagicMock, patch

from .. import jobdir_builder


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.mod_mocks = {}
        self.job = MagicMock()
        self.cfg = MagicMock()
        self.output_dir = '/some/output/dir'
        self.decorate_job_spec_fn = MagicMock()
        self.builder = jobdir_builder.JobdirBuilder(
            job=self.job,
            cfg=self.cfg,
            output_dir=self.output_dir,
            decorate_job_spec_fn=self.decorate_job_spec_fn
        )

    def mockify_builder_attrs(self, attrs=None):
        for attr in attrs: setattr(self.builder, attr, MagicMock())

    def mockify_module_attrs(self, module=jobdir_builder, attrs=None):
        mocks = {}
        for attr in attrs:
            patcher = patch.object(module, attr)
            self.addCleanup(patcher.stop)
            mocks[attr] = patcher.start()
        self.mod_mocks.update(mocks)

class BuildJobdirTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_builder_attrs(attrs=['_write_entrypoint',
                                          '_generate_job_spec'])
        self.result = self.builder.build_jobdir()

    def test_writes_entrypoint(self):
        self.assertEqual(self.builder._write_entrypoint.call_args,
                         call())

    def test_generates_job_spec(self):
        self.assertEqual(self.builder._generate_job_spec.call_args,
                         call())

    def test_returns_job_spec(self):
        self.assertEqual(self.result,
                         self.builder._generate_job_spec.return_value)

class _WriteEntrypointTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_builder_attrs(attrs=['_generate_entrypoint_content'])
        self.mockify_module_attrs(attrs=['open', 'os'])
        self.builder._write_entrypoint()

    def test_writes_entrypoint_file(self):
        self.assertEqual(self.mod_mocks['open'].call_args,
                         call(self.builder.entrypoint_path, 'w'))
        expected_fh = self.mod_mocks['open'].return_value.__enter__.return_value
        self.assertEqual(
            expected_fh.write.call_args,
            call(self.builder._generate_entrypoint_content.return_value)
        )

    def test_makes_entrypoint_executable(self):
        self.assertEqual(self.mod_mocks['os'].chmod.call_args,
                         call(self.builder.entrypoint_path, 0o755))

class _GenerateEntrypointContentTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.result = self.builder._generate_entrypoint_content()

    def test_generates_expected_content(self):
        expected_content = textwrap.dedent(
            """
            #!/bin/bash
            ${RUN_JOBDIR_CMD_VARNAME}
            """
        ).lstrip().format(
            RUN_JOBDIR_CMD_VARNAME=self.builder.RUN_JOBDIR_CMD_VARNAME
        )
        self.assertEqual(self.result, expected_content)

class _GenerateJobSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.mockify_builder_attrs(attrs=['_generate_base_job_spec'])
        self.result = self.builder._generate_job_spec()

    def test_generates_base_job_spec(self):
        self.assertEqual(self.builder._generate_base_job_spec.call_args,
                         call())

    def test_decorates_base_job_spec(self):
        self.assertEqual(
            self.builder.decorate_job_spec_fn.call_args,
            call(job_spec=self.builder._generate_base_job_spec.return_value,
                 job=self.builder.job, cfg=self.builder.cfg)
        )

    def test_returns_decorated_job_spec(self):
        self.assertEqual(self.result,
                         self.builder.decorate_job_spec_fn.return_value)

class _GenerateBaseJobSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.result = self.builder._generate_base_job_spec()

    def test_returns_expected_job_spec(self):
        expected_job_spec = {
            'cfg_specs': {
                self.builder.RUN_JOBDIR_CMD_VARNAME: {'required': True}
            },
            'dir': self.builder.output_dir,
            'entrypoint': self.builder.ENTRYPOINT_NAME,
        }
        self.assertEqual(self.result, expected_job_spec)
