import collections
import copy

from mc.mc_utils import context_value_setter
from mc.task_handlers.base_task_handler import BaseTaskHandler
from mc.flow_engines import flow_engine


class RunFlowTaskHandler(BaseTaskHandler):
    def initial_tick(self, task=None, task_context=None):
        # @TODO CLEAN THIS UP! TOTALLY A HACK.
        formatted_task_params = self.format_task_params(
            task['task_params'], task_context=task_context)
        flow_kwargs = self.generate_flow_kwargs(
            flow_spec=formatted_task_params['flow_spec'])
        created_flow = task_context['flow_ctx']['create_flow'](
            flow_kwargs=flow_kwargs)
        if formatted_task_params.get('delete_flow_spec_after_creation'):
            task['task_params'].pop('flow_spec', None)
        task['data']['flow_uuid'] = created_flow['uuid']

    def format_task_params(self, task_params=None, task_context=None):
        if isinstance(task_params, collections.Mapping) and \
           task_params.get('_task_params_type') == 'interpolated_object':
            obj = copy.deepcopy(task_params.get('object_template', {}))
            context_value_setter.set_context_values(
                value_specs=task_params.get('interpolations', []),
                context={'task_context': task_context, 'object': obj})
            formatted_task_params = obj
        else:
            formatted_task_params = {**task_params}
        return formatted_task_params

    def generate_flow_kwargs(self, flow_spec=None):
        flow_kwargs = {
            'label': flow_spec.get('label'),
            'serialization': flow_engine.FlowEngine.serialize_flow(
                flow=flow_engine.FlowEngine.generate_flow(flow_spec=flow_spec)
            )
        }
        return flow_kwargs

    def intermediate_tick(self, task=None, task_context=None):
        flow_ctx = task_context['flow_ctx']
        flow = self.get_flow(task=task, flow_ctx=flow_ctx)
        assert flow is not None
        if flow['status'] == 'COMPLETED':
            task['status'] = 'COMPLETED'
        elif flow['status'] == 'FAILED':
            try: error = flow['data']['error']
            except KeyError: error = '<unknown>'
            raise Exception(error)

    def get_flow(self, task=None, flow_ctx=None):
        return flow_ctx['get_flow'](uuid=task['data']['flow_uuid'])
