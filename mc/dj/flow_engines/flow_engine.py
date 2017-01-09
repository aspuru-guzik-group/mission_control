import collections
import logging

from .flow import Flow


class FlowEngine(object):
    def __init__(self, ctx=None, logger=None):
        self.ctx = ctx
        self.logger = logger or logging
        self.task_class_registry = collections.OrderedDict()
        self.flow_generator_class_registry = collections.OrderedDict()

    def register_task_class(self, task_class=None):
        key = getattr(task_class, 'task_type', task_class.__name__)
        self.task_class_registry[key] = task_class

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
        for task_class in dependencies.get('task_classes', []):
            self.register_task_class(task_class)
        for generator_class in dependencies.get('flow_generator_classes', []):
            self.register_flow_generator_class(
                flow_generator_class=generator_class)


    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        for attr in ['data', 'input', 'output', 'status']:
            setattr(flow, attr, serialized_flow.get(attr, None))
        tasks = [
            self.deserialize_task(serialized_task)
            for serialized_task in serialized_flow.get('tasks', [])
        ]
        flow.add_tasks(tasks=tasks)
        root_task_key = serialized_flow.get('root_task_key')
        if root_task_key is not None:
            flow.root_task = flow.tasks[root_task_key]
        edges = [
            {'src': flow.tasks[serialized_edge['src_key']],
             'dest': flow.tasks[serialized_edge['dest_key']]}
            for serialized_edge in serialized_flow.get('edges', [])
        ]
        flow.add_edges(edges=edges)
        return flow

    def deserialize_task(self, serialized_task):
        task_class = self.get_task_class_for_serialized_task(
            serialized_task=serialized_task)
        if not task_class:
            msg = "Could not find task_class for serialized_task '{}'".format( 
                serialized_task)
            raise Exception(msg)
        task_kwargs = self.get_task_kwargs_for_serialized_task(
            serialized_task=serialized_task)
        task = task_class(ctx=self.ctx, **task_kwargs)
        return task

    def get_task_class_for_serialized_task(self, serialized_task=None):
        return self.task_class_registry.get(serialized_task['task_type'], None)

    def get_task_kwargs_for_serialized_task(self, serialized_task=None):
        return serialized_task

    def serialize_flow(self, flow=None):
        serialized_flow = {
            **{attr: getattr(flow, attr, None)
               for attr in ['data', 'input', 'output', 'status']},
            'tasks': self.serialize_tasks(tasks=flow.tasks.values()),
            'edges': self.serialize_edges(edges=flow.edges.values()),
        }
        if flow.root_task is not None:
            serialized_flow['root_task_key'] = flow.root_task.key
        return serialized_flow

    def serialize_tasks(self, tasks=None):
        serialized_tasks = [self.serialize_task(task=task) for task in tasks]
        return serialized_tasks

    def serialize_task(self, task=None):
        serialized_task = {}
        for attr in ['key', 'data', 'input', 'output', 'status']:
            serialized_task[attr] = getattr(task, attr, None)
        serialized_task['task_type'] = getattr(task, 'task_type',
                                               task.__class__.__name__)
        return serialized_task

    def serialize_edges(self, edges=None):
        serialized_edges = [self.serialize_edge(edge=edge) for edge in edges]
        return serialized_edges

    def serialize_edge(self, edge=None):
        serialized_edge = {
            'src_key': edge['src'].key,
            'dest_key': edge['dest'].key,
        }
        return serialized_edge

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
        self.set_task_input(flow=flow, task=task)
        task.status = 'RUNNING'

    def set_task_input(self, flow=None, task=None):
        existing_input = getattr(task, 'input', None)
        if existing_input is None:
            _get_output = lambda t: getattr(t, 'output', None)
            precursors = flow.get_precursors(task=task)
            if len(precursors) == 0: task.input = None
            elif len(precursors) == 1: task.input = _get_output(precursors[0])
            else: task.input = [_get_output(precursor)
                                for precursor in precursors]

    def tick_running_tasks(self, flow=None, ctx=None):
        for task in flow.get_tasks_by_status(status='RUNNING'):
            self.tick_task(task=task, ctx=ctx)

    def tick_task(self, task=None, ctx=None):
        try:
            task.tick(engine=self, ctx=ctx)
        except Exception as e:
            msg = "tick failed for task with key '{key}': {e}".format(
                key=task.key, e=e)
            self.logger.exception(msg)

    def complete_flow(self, flow=None):
        flow.status = 'COMPLETED'
        flow.output = self.get_flow_output(flow=flow)

    def get_flow_output(self, flow=None):
        _get_task_output = lambda t: getattr(t, 'output', None)
        tail_tasks = flow.get_tail_tasks()
        if len (tail_tasks) == 0: output = None
        elif len(tail_tasks) == 1: output = _get_task_output(tail_tasks[0])
        elif len(tail_tasks) > 1: output = {
            tail_task.key: _get_task_output(tail_task)
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
