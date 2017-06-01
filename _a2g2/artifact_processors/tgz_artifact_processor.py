import base64
import io
import tarfile

from .base_artifact_processor import BaseArtifactProcessor

class TgzArtifactProcessor(BaseArtifactProcessor):
    def dir_to_artifact(self, _dir=None):
        return {
            'artifact_type': 'tgz:bytes',
            'artifact_params': {
                'encoded_bytes': base64.b64encode(
                    self.dir_to_tgz_bytes(_dir=_dir)).decode(),
                'encoding': 'base64'
            }
        }

    def dir_to_tgz_bytes(self, _dir=None):
        mem_file = io.BytesIO()
        tgz = tarfile.open(mode='w:gz', fileobj=mem_file)
        tgz.add(_dir, arcname='.')
        tgz.close()
        return mem_file.getvalue()

    def artifact_to_dir(self, artifact=None, dest=None):
        if artifact['artifact_type'] == 'tgz:bytes':
            encoded_bytes = artifact['artifact_params']['encoded_bytes']
            encoding = artifact['artifact_params'].get('encoding', 'base64')
            tgz_bytes = None
            if encoding == 'base64': tgz_bytes = base64.b64decode(encoded_bytes)
            self.tgz_bytes_to_dir(tgz_bytes=tgz_bytes, dest=dest)
        else: raise self.UnknownArtifactTypeError()

    def tgz_bytes_to_dir(self, tgz_bytes=None, dest=None, encoding='base64'):
        mem_file = io.BytesIO(tgz_bytes)
        tgz = tarfile.open(mode='r:gz', fileobj=mem_file)
        tgz.extractall(path=dest)
