class BaseFlowGenerator(object):
    @classmethod
    def get_dependencies(cls):
        return {}

    @classmethod
    def generate_flow(cls, *args, flow_spec=None, **kwargs):
        raise NotImplementedError
