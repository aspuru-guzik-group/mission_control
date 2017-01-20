import os
import tempfile


class CopyActionHandler(object):
    def __init__(self, transfer_client=None, mkdtemp=None):
        self.transfer_client = transfer_client
        self.mkdtemp = mkdtemp or tempfile.mkdtemp

    def __call__(self, params=None, ctx=None):
        intermediate_dir = self.get_intermediate_dir()
        self.transfer_client.copy(src=params['src'], dest=intermediate_dir)
        output = self.transfer_client.copy(src=intermediate_dir,
                                           dest=params['dest'])
        return output

    def get_intermediate_dir(self):
        tmp_dir = self.mkdtemp()
        intermediate_dir = os.path.join(tmp_dir, 'local_src')
        return intermediate_dir
