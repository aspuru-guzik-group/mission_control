import os

from .base_artifact_processor import BaseArtifactProcessor

class OdysseyArtifactProcessor(BaseArtifactProcessor):
    def dir_to_artifact(self, _dir=None):
        return {
            'artifact_type': 'a2g2.artifacts.odyssey',
            'artifact_params': {'path': _dir}
        }

    def artifact_to_dir(self, artifact=None, dest=None):
        if artifact['artifact_type'] == 'a2g2.artifacts.odyssey':
            os.symlink(artifact['artifact_params']['path'], dest)
        else: raise self.UnknownArtifactTypeError()
