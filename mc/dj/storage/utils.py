import os
from uuid import uuid4


def get_storage_backend(params=None):
    backend = FileSystemBackend()
    return backend

class FileSystemBackend(object):
    def __init__(self, *args, base_dir=None, **kwargs):
        if not base_dir: base_dir = self.get_base_dir_from_settings()
        self.base_dir = base_dir

    def get_base_dir_from_settings(self):
        attr_name = 'STORAGE_FILESYSTEM_BACKEND_BASEDIR'
        from django.conf import settings
        try:
            return getattr(settings, attr_name)
        except:
            raise Exception(
                "Could not find value for '%s' in settings." % attr_name)

    def get_key(self): return str(uuid4())

    def post_data(self, data=None, params=None):
        key = self.get_key()
        with open(self.get_path_for_key(key), 'wb') as f: f.write(data)
        return {'key': key}

    def get_path_for_key(self, key=None):
        return os.path.join(self.base_dir, key)

    def get_data(self, params=None):
        key = params['key']
        with open(self.get_path_for_key(key), 'rb') as f: return f.read()
