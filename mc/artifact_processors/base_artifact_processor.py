class BaseArtifactProcessor(object):

    ARTIFACT_TYPE = None

    class UnknownArtifactTypeError(Exception): pass

    def dir_to_artifact(self, dir_=None, **kwargs):
        raise NotImplementedError

    def artifact_todir_(self, artifact=None, dest=None, **kwargs):
        raise NotImplementedError
