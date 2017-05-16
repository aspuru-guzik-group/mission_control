class BaseArtifactProcessor(object):
    class UnknownArtifactTypeError(Exception): pass

    def dir_to_artifact(self, _dir=None):
        raise NotImplementedError

    def artifact_to_dir(self, artifact=None, dest=None):
        raise NotImplementedError
