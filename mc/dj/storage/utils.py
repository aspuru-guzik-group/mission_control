import os
from uuid import uuid4


def get_storage_backend(params=None):
    backend = FileSystemBackend()
    return backend

class FileSystemBackend(object):
    def __init__(self, *args, base_dir=None, **kwargs):
        self.base_dir = base_dir

    def get_key(self): return str(uuid4())

    def post(self, data=None, params=None):
        key = self.get_key()
        with open(self.get_path_for_key(key), 'wb') as f: f.write(data)
        return {'key': key}

    def get_path_for_key(self, key=None):
        return os.path.join(self.base_dir, key)

    def get(self, params=None):
        key = params['key']
        with open(self.get_path_for_key(key), 'rb') as f:
            return {'data': f.read()}
