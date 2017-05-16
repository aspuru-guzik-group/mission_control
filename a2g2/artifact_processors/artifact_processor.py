import base64
import io
import tarfile

class ArtifactProcessor(object):
    class UnknownArtifactTypeError(Exception): pass

    def dir_to_artifact(self, _dir=None): raise NotImplementedError

    def artifact_to_dir(self, artifact=None, dest=None):
        if artifact['artifact_type'] == 'tgz:bytes':
            encoded_bytes = artifact['artifact_params']['encoded_bytes']
            encoding = artifact['artifact_params'].get('encoding', 'base64')
            tgz_bytes = None
            if encoding == 'base64': tgz_bytes = base64.b64decode(encoded_bytes)
            self.tgz_bytes_to_dir(tgz_bytes=tgz_bytes, dest=dest)
        else:
            raise self.UnknownArtifactTypeError()

    def tgz_bytes_to_dir(self, tgz_bytes=None, dest=None, encoding='base64'):
        mem_file = io.BytesIO(tgz_bytes)
        tgz = tarfile.open(mode='r:gz', fileobj=mem_file)
        tgz.extractall(path=dest)
