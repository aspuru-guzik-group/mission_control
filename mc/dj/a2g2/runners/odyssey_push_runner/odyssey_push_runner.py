class OdysseyPushRunner(object):
    def __init__(self, flow_client=None):
        self.flow_client = flow_client

    def create_flow(self, *args, flow=None, **kwargs):
        return self.flow_client.create_flow(flow=flow)
