import json
import os
import subprocess
import uuid
from . import process_utils


class FileSystemStorageClient(object):
    def __init__(self, storage_root_path=None, state_file_path=None):
        self.storage_root_path = storage_root_path
        self.state_file_path = state_file_path
        self.state = self.generate_initial_state()

    def deserialize_state(self):
        deserialized_state = self.generate_initial_state()
        if os.path.exists(self.state_file_path):
            with open(self.state_file_path) as f:
                deserialized_state.update(json.load(f))
        return deserialized_state

    def generate_initial_state(self):
        return {'processes': {}}

    def serialize_state(self):
        with open(self.state_file_path, 'w') as f:
            json.dump(self.state, f)

    def put(self, src_path=None):
        key = str(uuid.uuid4())
        tgt_path = os.path.join(self.storage_root_path, key)
        proc_meta = self.start_copy(src_path=src_path, tgt_path=tgt_path)
        self.state['processes'][key] = proc_meta
        self.serialize_state()
        return key

    def start_copy(self, src_path=None, tgt_path=None):
        cmd = ['cp',  '-a', src_path, tgt_path]
        proc = subprocess.Popen(cmd)
        proc_meta = {'pid': proc.pid}
        return proc_meta

    def poll(self, key=None):
        self.state = self.deserialize_state()
        if key not in self.state['processes']:
            process_state = None
        else:
            process_state = process_utils.get_process_state(
                pid=self.state['processes'][key]['pid'])
            process_state['transferring'] = process_state['running']
            self.state['processes'][key] = process_state
            self.serialize_state()
        return process_state
