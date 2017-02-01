import collections
import logging

from .flow import Flow


class FlowEngine(object):
    def __init__(self, ctx=None, logger=None):
        self.ctx = ctx
        self.logger = logger or logging
        self.task_engine_registry = collections.OrderedDict()
        self.flow_generator_class_registry = collections.OrderedDict()

    def register_task_engine(self, task_engine=None):
        key = getattr(task_engine, 'task_engine', task_engine.__name__)
        self.task_engine_registry[key] = task_engine

    def register_flow_generator_class(self, flow_generator_class=None):
        key = getattr(flow_generator_class, 'flow_type',
                      flow_generator_class.__name__)
        self.flow_generator_class_registry[key] = flow_generator_class
        self.register_flow_generator_class_dependencies(
            flow_generator_class=flow_generator_class)

    def register_flow_generator_class_dependencies(self,
                                                   flow_generator_class=None):
        if not hasattr(flow_generator_class, 'get_dependencies'): return
        dependencies = flow_generator_class.get_dependencies()
        for task_engine in dependencies.get('task_enginees', []):
            self.register_task_engine(task_engine)
        for generator_class in dependencies.get('flow_generator_classes', []):
            self.register_flow_generator_class(
                flow_generator_class=generator_class)

    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        for attr in ['data', 'input', 'output', 'status', 'root_task_key']:
            setattr(flow, attr, serialized_flow.get(attr, None))
        for task in serialized_flow.get('tasks', []): flow.add_task(task=task)
        for edge in serialized_flow.get('edges', []): flow.add_edge(edge=edge)
        return flow

    def serialize_flow(self, flow=None):
        serialized_flow = {
            **{attr: getattr(flow, attr, None)
               for attr in ['data', 'input', 'output', 'status',
                            'root_task_key']},
            'tasks': flow.tasks.values(),
            'edges': flow.edges.values(),
        }
        return serialized_flow

    def tick_flow(self, flow=None, ctx=None):
        if flow.status == 'PENDING': self.start_flow(flow=flow)
        self.start_nearest_pending_tasks(flow=flow)
        self.tick_running_tasks(flow=flow, ctx=ctx)
        if not flow.has_incomplete_tasks():
            self.complete_flow(flow=flow)

    def start_flow(self, flow=None):
        flow.status = 'RUNNING'

    def start_nearest_pending_tasks(self, flow=None):
        for task in flow.get_nearest_pending_tasks():
            self.start_task(flow=flow, task=task)

    def start_task(self, flow=None, task=None):
        task['status'] = 'RUNNING'

    def tick_running_tasks(self, flow=None, ctx=None):
        for task in flow.get_tasks_by_status(status='RUNNING'):
            self.tick_task(task=task, ctx=ctx)

    def tick_task(self, task=None, ctx=None):
        try:
            task_engine = self.get_engine_for_task(task=task)
            task_engine.tick_task(task=task, ctx=ctx)
        except Exception as e:
            msg = "tick failed for task with key '{key}': {e}".format(
                key=task['key'], e=e)
            self.logger.exception(msg)

    def get_engine_for_task(self, task=None):
        task_engine =  self.task_engine_registry.get(task['task_engine'], None)
        if not task_engine:
            raise Exception("Could not find task_engine for task '%s'" % task)

    def complete_flow(self, flow=None):
        flow.output = self.get_flow_output(flow=flow)
        flow.status = 'COMPLETED'

    def get_flow_output(self, flow=None):
        _get_task_output = lambda t: t.get('output', None)
        tail_tasks = flow.get_tail_tasks()
        if len (tail_tasks) == 0: output = None
        elif len(tail_tasks) == 1: output = _get_task_output(tail_tasks[0])
        elif len(tail_tasks) > 1: output = {
            tail_task['key']: _get_task_output(tail_task)
            for tail_task in tail_tasks
        }
        return output

    def generate_flow(self, flow_spec=None):
        FlowGenerator = self.get_flow_generator_class_for_spec(flow_spec)
        if not FlowGenerator:
            msg = "Could not find FlowGenerator for flow_spec '{}'".format( 
                flow_spec)
            raise Exception(msg)
        flow_kwargs = self.get_flow_kwargs_for_flow_spec(flow_spec=flow_spec)
        flow = FlowGenerator.generate_flow(**flow_kwargs)
        return flow

    def get_flow_generator_class_for_spec(self, flow_spec=None):
        return self.flow_generator_class_registry.get(
            flow_spec['flow_type'], None)

    def get_flow_kwargs_for_flow_spec(self, flow_spec=None):
        return {'flow_spec': flow_spec}
