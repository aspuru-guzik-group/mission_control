from flow_engines.flow_engine import FlowEngine

class OdysseyPushRunner(object):
    def __init__(self, *args, run_setup=True, **kwargs):
        if run_setup: self.setup(**kwargs)

    def setup(self,
              flow_generator_registry=None, 
              flow_engine=None, 
              flow_client=None, 
              job_client=None,
              job_runner=None,
              flow_runner=None):
        self.flow_generator_registry = flow_generator_registry or \
                self.generate_flow_generator_registry()
        self.flow_engine = flow_engine or self.generate_flow_engine()
        self.flow_client = flow_client or self.generate_flow_client()
        self.job_client = job_client or self.generate_job_client()
        self.job_runner = job_runner or self.generate_job_runner()
        self.flow_runner = flow_runner or self.generate_flow_runner()

    def generate_flow_generator_registry(self):
        registry = set()
        from a2g2.flow_generators import reaxys
        registry.add(reaxys.ReaxysFlowGenerator)
        return registry

    def generate_flow_engine(self, flow_generator_classes=None):
        flow_engine = FlowEngine()
        for flow_generator_class in (flow_generator_classes or []):
            flow_engine.register_flow_generator_class(
                flow_generator_class=flow_generator_class)
        return flow_engine

    def generate_flow_client(self): pass
    def generate_job_client(self): pass
    def generate_job_runner(self): pass
    def generate_flow_runner(self): pass

    def create_flow(self, *args, flow=None, **kwargs):
        return self.flow_client.create_flow(flow=flow)
