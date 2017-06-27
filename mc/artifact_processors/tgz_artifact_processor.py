import base64
import io
import tarfile

from .base_artifact_processor import BaseArtifactProcessor

class TgzArtifactProcessor(BaseArtifactProcessor):
    """ArtifactProcessor that converts dir <=> inline compressed bytes."""

    ARTIFACT_TYPE = 'tgz:bytes'

    def dir_to_artifact(self, dir_=None, **kwargs):
        return {
            'artifact_type': self.ARTIFACT_TYPE,
            'artifact_params': {
                'encoded_bytes': base64.b64encode(
                    self.dir_to_tgz_bytes(dir_=dir_)).decode(),
                'encoding': 'base64'
            }
        }

    def dir_to_tgz_bytes(self, dir_=None):
        mem_file = io.BytesIO()
        tgz = tarfile.open(mode='w:gz', fileobj=mem_file)
        tgz.add(dir_, arcname='.')
        tgz.close()
        return mem_file.getvalue()

    def artifact_to_dir_(self, artifact=None, dest=None, **kwargs):
        if artifact['artifact_type'] == self.ARTIFACT_TYPE:
            encoded_bytes = artifact['artifact_params']['encoded_bytes']
            encoding = artifact['artifact_params'].get('encoding', 'base64')
            tgz_bytes = None
            if encoding == 'base64': tgz_bytes = base64.b64decode(encoded_bytes)
            self.tgz_bytes_to_dir_(tgz_bytes=tgz_bytes, dest=dest)
        else: raise self.UnknownArtifactTypeError()

    def tgz_bytes_to_dir_(self, tgz_bytes=None, dest=None, encoding='base64'):
        mem_file = io.BytesIO(tgz_bytes)
        tgz = tarfile.open(mode='r:gz', fileobj=mem_file)
        tgz.extractall(path=dest)
