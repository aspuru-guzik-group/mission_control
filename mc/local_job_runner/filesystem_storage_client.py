import json
import os
import psutil
import subprocess
import uuid


class FileSystemStorageClient(object):
    def __init__(self, root_path=None, serialization_dir=None):
        self.root_path = root_path
        self.serialization_dir = serialization_dir
        self.serialized_state_path = os.path.join(self.serialization_dir,
                                                  'state.json')
        self.state = self.generate_initial_state()

    def deserialize_state(self):
        deserialized_state = self.generate_initial_state()
        if os.path.exists(self.serialized_state_path):
            with open(self.serialized_state_path) as f:
                deserialized_state.update(json.load(f))
        return deserialized_state

    def generate_initial_state(self):
        return {'processes': {}}

    def serialize_state(self):
        with open(self.serialized_state_path, 'w') as f:
            json.dump(self.state, f)

    def put(self, src_path=None):
        key = uuid.uuid4()
        tgt_path = os.path.join(self.root_path, key)
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
            copy_process_state = None
        else:
            copy_process_state = psutil.Process(
                pid=self.state['processes'][key]['pid']).as_dict()
            self.state['processes'][key] = copy_process_state
            self.serialize_state()
        return copy_process_state
