import json
import io
import tarfile

from mc.task_handlers.base_task_handler import BaseTaskHandler

class StorageTaskHandler(BaseTaskHandler):
    def __init__(self, *args, storage_client=None, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.storage_client = storage_client

class UploadTaskHandler(StorageTaskHandler):
    def initial_tick(self, *args, task=None, job=None, **kwargs):
        dir_path = task['params']['src']
        tgz_bytes = dir_to_tgz_bytes(dir_path=dir_path)
        storage_meta = self.storage_client.post_data(data=tgz_bytes)
        task['data']['storage_meta_json'] = json.dumps(storage_meta)
        task['status'] = 'COMPLETED'

def dir_to_tgz_bytes(dir_path=None):
    mem_file = io.BytesIO()
    tgz = tarfile.open(mode='w:gz', fileobj=mem_file)
    tgz.add(dir_path, arcname='.')
    tgz.close()
    return mem_file.getvalue()

class DownloadTaskHandler(StorageTaskHandler):
    def initial_tick(self, *args, task=None, job=None, **kwargs):
        serialized_storage_meta = task['params']['storage_meta_json']
        storage_meta = json.loads(serialized_storage_meta)
        src_params = storage_meta['params']
        dest = task['params']['dest']
        tgz_bytes = self.storage_client.get_data(storage_params=src_params)
        tgz_bytes_to_dir(tgz_bytes=tgz_bytes, dir_path=dest)
        task['data']['path_to_download'] = dest
        task['status'] = 'COMPLETED'

def tgz_bytes_to_dir(tgz_bytes=None, dir_path=None):
    mem_file = io.BytesIO(tgz_bytes)
    tgz = tarfile.open(mode='r:gz', fileobj=mem_file)
    tgz.extractall(path=dir_path)
