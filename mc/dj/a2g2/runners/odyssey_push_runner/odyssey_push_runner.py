class OdysseyPushRunner(object):
    def __init__(self,
                 flow_engine=None, 
                 flow_client=None, 
                 job_client=None,
                 job_runner=None,
                 flow_runner=None):
        self.flow_engine = flow_engine or self.generate_flow_engine()
        self.flow_client = flow_client or self.generate_flow_client()
        self.job_client = job_client or self.generate_job_client()
        self.job_runner = job_runner or self.generate_job_runner()
        self.flow_runner = flow_runner or self.generate_flow_runner()

    def create_flow(self, *args, flow=None, **kwargs):
        return self.flow_client.create_flow(flow=flow)

    def generate_flow_engine(self): pass
    def generate_flow_client(self): pass
    def generate_job_client(self): pass
    def generate_job_runner(self): pass
    def generate_flow_runner(self): pass
