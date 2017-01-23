import os
import tempfile


class TransferClient(object):
    staging_dir_name = 'staging'

    def __init__(self, mkdtemp=None):
        self.mkdtemp = mkdtemp or tempfile.mkdtemp
        self.subclients_by_storage_key = {}

    def register_subclient(self, subclient=None, key=None):
        self.subclients_by_storage_key[key] = subclient

    def copy(self, src=None, dest=None):
        tmpdir = self.mkdtemp()
        staging_dir = os.path.join(tmpdir, self.staging_dir_name)
        self.download(src=src, dest=staging_dir)
        return self.upload(src=staging_dir, dest=dest)

    def download(self, src=None, dest=None):
        subclient = self.get_subclient_for_transfer_uri(transfer_uri=src)
        return subclient.download(src=src, dest=dest)

    def upload(self, src=None, dest=None):
        subclient = self.get_subclient_for_transfer_uri(transfer_uri=dest)
        return subclient.upload(src=src, dest=dest)

    def get_subclient_for_transfer_uri(self, transfer_uri=None):
        try:
            storage_key = transfer_uri.split(':')[0]
        except Exception as error:
            error_msg = "Invalid transfer_uri '%s': %s" % (transfer_uri, error)
            raise Exception(error_msg)
        subclient = self.subclients_by_storage_key.get(storage_key, None)
        if not subclient:
            error_msg = "Could not get subclient for storage_key '%s'" % (
                storage_key)
            raise Exception(error_msg)
        return subclient
