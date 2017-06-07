import os

from .base_artifact_processor import BaseArtifactProcessor

class LocalPathArtifactProcessor(BaseArtifactProcessor):
    ARTIFACT_TYPE = 'local_path'

    def dir_to_artifact(self, dir_=None, **kwargs):
        return {
            'artifact_type': self.ARTIFACT_TYPE,
            'artifact_params': {'path': os.path.abspath(dir_)}
        }

    def artifact_to_dir(self, artifact=None, dest=None, **kwargs):
        if artifact['artifact_type'] == self.ARTIFACT_TYPE:
            os.symlink(artifact['params']['path'], dest)
        else: raise self.UnknownArtifactTypeError()
