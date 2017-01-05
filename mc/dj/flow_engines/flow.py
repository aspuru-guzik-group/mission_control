import collections
from uuid import uuid4

class BaseTask(object):
    def __init__(self, id=None, flow=None, status=None, state=None,
                 **kwargs):
        self.id = id or str(uuid4())
        self.flow = flow
        self.status = status
        self.state = state

    @property
    def task_type(self): return self.__class__.__name__

class Flow(object):
    def __init__(self, jobs=None):
        self.jobs = jobs or {}
        self.tasks = {}
        self.edges = {}
        self.edges_by_task_id = collections.defaultdict(dict)
        self.state = None
        self.status = None
        self.root_task = None

    def add_tasks(self, tasks=None):
        for task in tasks: self.add_task(task=task)

    def add_task(self, task=None, as_root=False, precursor=None,
                 successor=None):
        assert hasattr(task, 'id')
        self.tasks[task.id] = task
        task.flow = self
        if as_root: self.root_task = task
        for precursor_task in self._ensure_iterable(precursor):
            self.add_edge(src=precursor_task, dest=task)
        for successor_task in self._ensure_iterable(successor):
            self.add_edge(src=task, dest=successor_task)
        return task

    def _ensure_iterable(self, obj=None):
        if obj is None: return []
        elif isinstance(obj, str): return [obj]
        elif isinstance(obj, collections.Iterable): return obj
        else: return [obj]

    def add_edges(self, edges=None):
        for edge in edges: self.add_edge(**edge)

    def add_edge(self, src=None, dest=None):
        if dest is self.root_task:
            raise Exception("Root task can not be an edge dest")
        src_id, dest_id = [self.ensure_task_id(task) for task in [src, dest]]
        edge_key = self.get_edge_key(src=src_id, dest=dest_id)
        edge = {'src': self.tasks[src_id], 'dest': self.tasks[dest_id]}
        self.edges[edge_key] = edge
        self.edges_by_task_id[src_id][edge_key] = edge
        self.edges_by_task_id[dest_id][edge_key] = edge

    def get_edge_key(self, src=None, dest=None):
        return (self.ensure_task_id(src), self.ensure_task_id(dest))

    def has_edge(self, src=None, dest=None):
        return self.get_edge_key(src=src, dest=dest) in self.edges

    def ensure_task_id(self, task_or_task_id):
        if isinstance(task_or_task_id, str) or isinstance(task_or_task_id, int):
            task_id = task_or_task_id
        else:
            task_id = task_or_task_id.id
        return task_id

    def get_precursors(self, task=None):
        task_edges = self.edges_by_task_id[task.id].values()
        precursors = [e['src'] for e in task_edges if e['dest'] is task]
        return precursors

    def get_successors(self, task=None):
        task_edges = self.edges_by_task_id[task.id].values()
        successors = [e['dest'] for e in task_edges if e['src'] is task]
        return successors

    def get_nearest_pending_tasks(self):
        nearest_pending_tasks = []
        cursors = [self.root_task]
        while len(cursors) > 0:
            next_cursors = []
            for cursor in cursors:
                if cursor.status == 'PENDING':
                    nearest_pending_tasks.append(cursor)
                elif cursor.status == 'COMPLETED':
                    successors = self.get_successors(task=cursor)
                    next_cursors.extend(successors)
            cursors = next_cursors
        return nearest_pending_tasks

    def filter_tasks(self, filters=None):
        result = []
        for task in self.tasks.values():
            passes_filters = True
            for _filter in filters:
                if not _filter(task):
                    passes_filters = False
                    break
            if passes_filters:
                result.append(task)
        return result

    def get_tasks_by_status(self, status=None):
        state_filter = lambda task: task.status == status
        return self.filter_tasks(filters=[state_filter])

    def has_incomplete_tasks(self):
        filter_fn = lambda task: task.status != 'COMPLETED'
        incomplete_tasks = self.filter_tasks(filters=[filter_fn])
        return len(incomplete_tasks) > 0

