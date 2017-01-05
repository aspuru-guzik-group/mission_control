import collections
from uuid import uuid4

from .flow import Flow


class FlowEngine(object):
    def __init__(self, ctx=None):
        self.task_class_registry = collections.OrderedDict()
        self.ctx = ctx

    def register_task_class(self, task_class=None, key=None, test_fn=None):
        if key is None: key = str(uuid4())
        if not test_fn:
            def test_fn(serialized_task):
                return (serialized_task.get('task_type') == task_class.__name__)
        self.task_class_registry[key] = {'test_fn': test_fn,
                                         'task_class': task_class}

    def deserialize_flow(self, serialized_flow=None):
        flow = Flow()
        tasks = [
            self.deserialize_task(serialized_task)
            for serialized_task in serialized_flow.get('tasks', [])
        ]
        flow.add_tasks(tasks=tasks)
        root_task_id = serialized_flow.get('root_task_id')
        if root_task_id is not None:
            flow.root_task = flow.tasks[root_task_id]
        edges = [
            {'src': flow.tasks[serialized_edge['src_id']],
             'dest': flow.tasks[serialized_edge['dest_id']]}
            for serialized_edge in serialized_flow.get('edges', [])
        ]
        flow.add_edges(edges=edges)
        flow.status = serialized_flow['status']
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
        for registry_entry in self.task_class_registry.values():
            if registry_entry['test_fn'](serialized_task):
                return registry_entry['task_class']
        return None

    def get_task_kwargs_for_serialized_task(self, serialized_task=None):
        return serialized_task

    def serialize_flow(self, flow=None):
        serialized_flow = {
            'tasks': self.serialize_tasks(tasks=flow.tasks.values()),
            'edges': self.serialize_edges(edges=flow.edges.values()),
            'state': flow.state,
            'status': flow.status,
        }
        if flow.root_task is not None:
            serialized_flow['root_task_id'] = flow.root_task.id
        return serialized_flow

    def serialize_tasks(self, tasks=None):
        serialized_tasks = [self.serialize_task(task=task) for task in tasks]
        return serialized_tasks

    def serialize_task(self, task=None):
        serialized_task = {
            'id': task.id,
            'state': task.state,
            'task_type': task.task_type,
            'status': task.status,
        }
        return serialized_task

    def serialize_edges(self, edges=None):
        serialized_edges = [self.serialize_edge(edge=edge) for edge in edges]
        return serialized_edges

    def serialize_edge(self, edge=None):
        serialized_edge = {
            'src_id': edge['src'].id,
            'dest_id': edge['dest'].id,
        }
        return serialized_edge

    def tick_flow(self, flow=None, ctx=None):
        self.start_nearest_pending_tasks(flow=flow)
        self.tick_running_tasks(flow=flow, ctx=ctx)
        if not flow.has_incomplete_tasks():
            flow.status = 'COMPLETED'

    def start_nearest_pending_tasks(self, flow=None):
        for task in flow.get_nearest_pending_tasks():
            self.start_task(flow=flow, task=task)

    def start_task(self, flow=None, task=None):
        self.set_task_input(flow=flow, task=task)
        task.status = 'RUNNING'

    def set_task_input(self, flow=None, task=None):
        static_input = getattr(task, 'static_input', None)
        if static_input is not None:
            task.input = task.static_input
        else:
            _get_output = lambda t: getattr(t, 'output', None)
            precursors = flow.get_precursors(task=task)
            if len(precursors) == 0: task.input = None
            elif len(precursors) == 1: task.input = _get_output(precursors[0])
            else: task.input = [_get_output(precursor)
                                for precursor in precursors]

    def tick_running_tasks(self, flow=None, ctx=None):
        if not ctx: ctx = {}
        ctx_with_jobs = {**ctx, 'jobs': getattr(flow, 'jobs', {})}
        for task in flow.get_tasks_by_status(status='RUNNING'):
            self.tick_task(task=task, ctx=ctx_with_jobs)

    def tick_task(self, task=None, ctx=None):
        task.tick(engine=self, ctx=ctx)
