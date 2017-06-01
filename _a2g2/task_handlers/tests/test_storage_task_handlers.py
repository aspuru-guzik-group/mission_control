from collections import defaultdict
import json
import os
import tempfile
import unittest
from unittest.mock import call, MagicMock

from mc.mc_utils import test_utils as mc_test_utils

from .. import storage_task_handlers

class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.task = self.generate_task()
        self.job = defaultdict(MagicMock)
        self.storage_client = MagicMock()
        self.task_handler = self.generate_task_handler()

    def generate_task(self):
        task = defaultdict(MagicMock)
        task['data'] = defaultdict(MagicMock)
        return task

    def generate_task_handler(self): pass

    def _do_initial_tick(self):
        self.task_handler.initial_tick(task=self.task, job=self.job)

class UploadTaskHandlerBaseTestCase(BaseTestCase):
    def generate_task_handler(self):
        return storage_task_handlers.UploadTaskHandler(
            storage_client=self.storage_client)

class UploadTaskHandlerInitialTickTestCase(UploadTaskHandlerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.storage_client.post_data.return_value = {'some': 'object'}
        self.src_dir = self.setup_src_dir()
        self.task['params'] = {'src': self.src_dir}

    def setup_src_dir(self):
        src_dir = tempfile.mkdtemp()
        for i in range(3):
            path = os.path.join(src_dir, 'file_%s' % i)
            with open(path, 'w') as f: f.write(str(i))
        return src_dir

    def test_posts_tgz_of_src(self):
        self._do_initial_tick()
        expected_data = storage_task_handlers.dir_to_tgz_bytes(
            dir_path=self.task['params']['src'])
        self.assertEqual(self.storage_client.post_data.call_args,
                         call(data=expected_data))

    def test_stores_serialized_result_of_post_data(self):
        self._do_initial_tick()
        self.assertEqual(self.task['data']['storage_meta_json'],
                         json.dumps(self.storage_client.post_data.return_value))

    def test_sets_task_status_to_completed(self):
        self._do_initial_tick()
        self.assertEqual(self.task['status'], 'COMPLETED')

class DownloadTaskHandlerBaseTestCase(BaseTestCase):
    def generate_task_handler(self):
        return storage_task_handlers.DownloadTaskHandler(
            storage_client=self.storage_client)

class DownloadTaskHandlerInitialTickTestCase(DownloadTaskHandlerBaseTestCase):
    def setUp(self):
        super().setUp()
        self.src_params = {'some': 'src params'}
        self.src_dir = self.setup_src_dir()
        self.src_tgz_bytes = storage_task_handlers.dir_to_tgz_bytes(
            dir_path=self.src_dir)
        self.storage_client.get_data.return_value = self.src_tgz_bytes
        self.dest_dir = tempfile.mkdtemp()
        self.dest = os.path.join(self.dest_dir, 'dest')
        self.task['params'] = {
            'storage_meta_json': json.dumps({'params': self.src_params}),
            'dest': self.dest
        }

    def setup_src_dir(self):
        src_dir = tempfile.mkdtemp()
        for i in range(3):
            path = os.path.join(src_dir, 'file_%s' % i)
            with open(path, 'w') as f: f.write(str(i))
        return src_dir

    def test_calls_get_data_with_deserialized_src_params(self):
        self._do_initial_tick()
        expected_storage_meta = json.loads(
            self.task['params']['storage_meta_json'])
        self.assertEqual(
            self.storage_client.get_data.call_args,
            call(storage_params=expected_storage_meta['params'])
        )

    def test_writes_tgz_to_dest(self):
        self._do_initial_tick()
        mc_test_utils.assert_dirs_equal(
            test=self, left=self.src_dir, right=self.task['params']['dest'])

    def test_stores_dest_as_path_to_download(self):
        self._do_initial_tick()
        self.assertEqual(self.task['data']['path_to_download'],
                        self.task['params']['dest'])

    def test_sets_task_status_to_completed(self):
        self._do_initial_tick()
        self.assertEqual(self.task['status'], 'COMPLETED')
