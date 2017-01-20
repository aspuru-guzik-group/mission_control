import os
import tempfile


class TransferClient(object):
    staging_dir_name = 'staging'

    def __init__(self, mkdtemp=None):
        self.mkdtemp = mkdtemp or tempfile.mkdtemp

    def copy(self, src=None, dest=None):
        tmpdir = self.mkdtemp()
        staging_dir = os.path.join(tmpdir, self.staging_dir_name)
        self.download(src=src, dest=staging_dir)
        return self.upload(src=staging_dir, dest=dest)

    def download(self, src=None, dest=None):
        subclient = self.get_subclient_for_transfer_spec(transfer_spec=src)
        return subclient.download(src=src, dest=dest)

    def upload(self, src=None, dest=None):
        subclient = self.get_subclient_for_transfer_spec(transfer_spec=dest)
        return subclient.upload(src=src, dest=dest)

    def get_subclient_for_transfer_spec(self, transfer_spec=None):
        raise NotImplementedError
