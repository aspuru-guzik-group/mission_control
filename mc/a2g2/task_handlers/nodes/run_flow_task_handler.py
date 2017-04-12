from mc.task_handlers.base_task_handler import BaseTaskHandler


class RunFlowTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        flow_ctx = task_context['flow_ctx']
        create_flow_fn = flow_ctx['create_flow']
        flow_kwargs = {'flow_spec': task['task_params']['flow_spec']}
        flow = create_flow_fn(flow_kwargs=flow_kwargs)
        task['data']['flow_uuid'] = flow['uuid']

    def intermediate_tick(self, task=None, task_context=None):
        flow_ctx = task_context['flow_ctx']
        flow = self.get_flow(task=task, flow_ctx=flow_ctx)
        assert flow is not None
        if flow['status'] == 'COMPLETED':
            task['data']['artifact'] = flow['data'].get('artifact')
            task['status'] = 'COMPLETED'
        elif flow['status'] == 'FAILED':
            try: error = flow['data']['error']
            except KeyError: error = '<unknown>'
            raise Exception(error)

    def get_flow(self, task=None, flow_ctx=None):
        return flow_ctx['get_flow'](uuid=task['data']['flow_uuid'])
