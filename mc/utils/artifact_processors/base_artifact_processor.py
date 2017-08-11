import abc


class BaseArtifactProcessor(abc.ABC):
    """Abstract base class for artifact processors."""

    ARTIFACT_TYPE = None

    class InvalidArtifactError(Exception): pass

    @abc.abstractmethod
    def dir_to_artifact(self, dir=None, **kwargs):
        raise NotImplementedError

    @abc.abstractmethod
    def artifact_to_dir(self, artifact=None, dest=None, **kwargs):
        raise NotImplementedError
