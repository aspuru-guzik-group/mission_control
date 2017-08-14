class BaseArchiver(object):
    def __init__(self, *args, **kwargs):
        pass

    def ingest_dir(self, dir_=None): raise NotImplementedError

    def materialize_as_path(self, meta=None): raise NotImplementedError
