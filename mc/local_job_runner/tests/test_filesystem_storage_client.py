import json
import os
import tempfile
import unittest
from unittest.mock import call, DEFAULT, patch

from ..filesystem_storage_client import FileSystemStorageClient


class FileSystemStorageClientBaseTestCase(unittest.TestCase):
    def setUp(self):
        self.dirs = self._generate_dirs()
        self.client = FileSystemStorageClient(
            storage_root_path=self.dirs['storage_root'],
            state_file_path=os.path.join(self.dirs['serialization'],
                                         'state.json')
        )

    def _generate_dirs(self):
        tmp_dir = tempfile.mkdtemp(prefix='localfs.test')
        dirs = {}
        for subdir in ['storage_root', 'serialization']:
            subdir_path = os.path.join(tmp_dir, subdir)
            os.makedirs(subdir_path)
            dirs[subdir] = subdir_path
        return dirs

class PutTestCase(FileSystemStorageClientBaseTestCase):
    def setUp(self):
        super().setUp()
        self.src_path = 'src_path'
        self.patchers = { 
            'uuid': patch.multiple('uuid', uuid4=DEFAULT),
            'client': patch.multiple(self.client, start_copy=DEFAULT,
                                     serialize_state=DEFAULT),
        }
        self.mocks = {patcher_key: patcher.start()
                      for patcher_key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def test_generates_key(self):
        key = self.client.put(src_path=self.src_path)
        self.assertEqual(key, str(self.mocks['uuid']['uuid4'].return_value))

    def test_starts_copy(self):
        key = self.client.put(src_path=self.src_path)
        expected_tgt_path = os.path.join(self.client.storage_root_path, key)
        self.assertEqual(
            self.mocks['client']['start_copy'].call_args,
            call(src_path=self.src_path, tgt_path=expected_tgt_path))

    def test_stores_copy_proc_meta(self):
        key = self.client.put(src_path=self.src_path)
        self.assertEqual(
            self.client.state['processes'][key],
            self.mocks['client']['start_copy'].return_value)

    def test_serializes_state(self):
        self.client.put(src_path=self.src_path)
        self.assertEqual(self.mocks['client']['serialize_state'].call_count, 1)

class StartCopyTestCase(FileSystemStorageClientBaseTestCase):
    def setUp(self):
        super().setUp()
        self.src_path = 'src_path'
        self.tgt_path = 'tgt_path'
        self.patchers = { 
            'subprocess': patch.multiple('subprocess', Popen=DEFAULT),
        }
        self.mocks = {patcher_key: patcher.start()
                      for patcher_key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def test_starts_copy_proc(self):
        self.client.start_copy(src_path=self.src_path, tgt_path=self.tgt_path)
        expected_cmd = ['cp', '-a', self.src_path, self.tgt_path]
        self.assertEqual(
            self.mocks['subprocess']['Popen'].call_args,
            call(expected_cmd))

    def test_returns_proc_meta(self):
        proc_meta = self.client.start_copy(src_path=self.src_path,
                                           tgt_path=self.tgt_path)
        expected_proc_meta = {
            'pid': self.mocks['subprocess']['Popen'].return_value.pid}
        self.assertEqual(proc_meta, expected_proc_meta)

class StateSerializationTestCase(FileSystemStorageClientBaseTestCase):
    def test_serialize_state(self):
        self.client.state = {'pie': 'rhubarb'}
        self.client.serialize_state()
        with open(self.client.state_file_path) as f:
            serialized_state = f.read()
        expected_serialized_state = json.dumps(self.client.state)
        self.assertEqual(serialized_state, expected_serialized_state)

    def test_deserialize_state(self):
        serialized_state = {"pie": "rhubarb"}
        with open(self.client.state_file_path, 'w') as f:
            json.dump(serialized_state, f)
        deserialized_state = self.client.deserialize_state()
        expected_deserialized_state = {**serialized_state,
                                       **self.client.generate_initial_state()}
        self.assertEqual(deserialized_state, expected_deserialized_state)


class PollTestCase(FileSystemStorageClientBaseTestCase):
    def setUp(self):
        super().setUp()
        self.key = 'my_key'
        self.patchers = { 
            'client': patch.multiple(self.client, start_copy=DEFAULT,
                                     deserialize_state=DEFAULT,
                                     serialize_state=DEFAULT),
            'psutil': patch.multiple('psutil', Process=DEFAULT),
        }
        self.mocks = {patcher_key: patcher.start()
                      for patcher_key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def test_deserializes_state(self):
        self.client.poll(key=self.key)
        self.assertEqual(self.mocks['client']['deserialize_state'].call_count,
                         1)

    def test_get_process_info_for_known_process(self):
        pid = 'my_pid'
        self.mocks['client']['deserialize_state'].return_value = {
            'processes': {self.key: {'pid': pid}}
        }
        self.client.state['processes'][self.key] = {'pid': pid}
        poll_result = self.client.poll(key=self.key)
        self.assertEqual(self.mocks['psutil']['Process'].call_args,
                         call(pid=pid))
        mock_process_state = self.mocks['psutil']['Process']\
                .return_value.as_dict.return_value
        self.assertEqual(poll_result, mock_process_state)
        self.assertEqual(self.mocks['client']['serialize_state'].call_count, 1)

    def test_get_process_info_for_unknown_process(self):
        self.mocks['client']['deserialize_state'].return_value = {
            'processes': {}
        }
        poll_result = self.client.poll(key=self.key)
        self.assertEqual(poll_result, None)


if __name__ == '__main__':
    unittest.main()
