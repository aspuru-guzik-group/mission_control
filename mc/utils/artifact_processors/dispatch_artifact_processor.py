from .base_artifact_processor import BaseArtifactProcessor

class DispatchArtifactProcessor(BaseArtifactProcessor):
    """ArtifactProcessor that dispatches to subprocessors, based on
    artifact type."""
    def __init__(self, processors=None):
        self.processors = processors or {}

    def dir_to_artifact(self, artifact_type=None, **kwargs):
        processor = self.get_processor_for_artifact_type(
            artifact_type=artifact_type)
        return processor.dir_to_artifact(**kwargs)

    def get_processor_for_artifact_type(self, artifact_type=None):
        try: return self.processors[artifact_type]
        except KeyError: raise self.UnknownArtifactTypeError()

    def artifact_to_dir(self, artifact=None, **kwargs):
        processor = self.get_processor_for_artifact_type(
            artifact_type=artifact['type'])
        return processor.artifact_to_dir(artifact=artifact, **kwargs)
