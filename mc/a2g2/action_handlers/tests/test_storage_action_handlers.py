import filecmp
import json
import os
import tempfile
import unittest
from unittest.mock import call, Mock

from .. import storage_action_handlers

class StorageActionHandlerTestCase(unittest.TestCase):
    def setUp(self):
        self.storage_client = Mock()
        self.ctx = self.generate_ctx()

    def generate_ctx(self):
        ctx = Mock()
        def transform_value(value): return value
        ctx.transform_value.side_effect = transform_value
        return ctx

class UploadActionHandlerTestCase(StorageActionHandlerTestCase):
    def setUp(self):
        super().setUp()
        self.storage_client.post_data.return_value = {'some': 'object'}
        self.src_dir = self.setup_src_dir()
        self.params = {'src': self.src_dir}

    def setup_src_dir(self):
        src_dir = tempfile.mkdtemp()
        for i in range(3):
            path = os.path.join(src_dir, 'file_%s' % i)
            with open(path, 'w') as f: f.write(str(i))
        return src_dir

    def do_upload(self):
        return storage_action_handlers.upload_action_handler(
            storage_client=self.storage_client,
            params=self.params, ctx=self.ctx)

    def test_posts_tgz_of_transformed_src(self):
        self.do_upload()
        self.assertEqual(self.ctx.transform_value.call_args,
                         call(self.params['src']))
        expected_data = storage_action_handlers.dir_to_tgz_bytes(
            dir_path=self.ctx.transform_value(self.params['src']))
        self.assertEqual(self.storage_client.post_data.call_args,
                         call(data=expected_data))

    def test_returns_serialized_result_of_post_data(self):
        result = self.do_upload()
        self.assertEqual(result,
                         json.dumps(self.storage_client.post_data.return_value))

class DownloadActionHandlerTestCase(StorageActionHandlerTestCase):
    def setUp(self):
        super().setUp()
        self.src_params = {'some': 'src params'}
        self.src_dir = self.setup_src_dir()
        self.src_tgz_bytes = storage_action_handlers.dir_to_tgz_bytes(
            dir_path=self.src_dir)
        self.storage_client.get_data.return_value = self.src_tgz_bytes
        self.dest_dir = tempfile.mkdtemp()
        self.dest = os.path.join(self.dest_dir, 'dest')
        self.params = {
            'storage_meta': json.dumps({'params': self.src_params}),
            'dest': self.dest
        }

    def setup_src_dir(self):
        src_dir = tempfile.mkdtemp()
        for i in range(3):
            path = os.path.join(src_dir, 'file_%s' % i)
            with open(path, 'w') as f: f.write(str(i))
        return src_dir

    def generate_ctx(self):
        ctx = Mock()
        def transform_value(value): return value
        ctx.transform_value.side_effect = transform_value
        return ctx

    def do_download(self):
        return storage_action_handlers.download_action_handler(
            storage_client=self.storage_client,
            params=self.params, ctx=self.ctx)

    def test_calls_get_data_with_deserialized_transformed_src_params(self):
        self.do_download()
        transformed_storage_meta = self.ctx.transform_value(
            self.params['storage_meta'])
        storage_meta = json.loads(transformed_storage_meta)
        self.assertEqual(self.storage_client.get_data.call_args,
                         call(storage_params=storage_meta['params']))

    def test_writes_tgz_to_transformed_dest(self):
        self.do_download()
        self.assert_dirs_equal(self.src_dir,
                               self.ctx.transform_value(self.params['dest']))

    def assert_dirs_equal(self, left, right):
        dircmp = filecmp.dircmp(left, right)
        self.assertEqual(dircmp.left_only, [])
        self.assertEqual(dircmp.right_only, [])

    def test_returns_transformed_dest(self):
        result = self.do_download()
        self.assertEqual(result, self.ctx.transform_value(self.params['dest']))
