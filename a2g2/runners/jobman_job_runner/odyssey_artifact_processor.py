import os

from .artifact_processor import ArtifactProcessor

class OdysseyArtifactProcessor(object):
    def __init__(self):
        self.artifact_processor = ArtifactProcessor()

    def dir_to_artifact(self, _dir=None):
        return {
            'artifact_type': 'a2g2.artifacts.odyssey',
            'artifact_params': {'path': _dir}
        }

    def artifact_to_dir(self, artifact=None, dest=None):
        if artifact['artifact_type'] == 'a2g2.artifacts.odyssey':
            os.symlink(artifact['artifact_params']['path'], dest)
        else:
            self.artifact_processor.artifact_to_dir(
                artifact=artifact, dest=dest)
