import collections
from uuid import uuid4


class Flow(object):
    def __init__(self, *args, input=None, output=None, data=None, status=None,
                 **kwargs):
        self.input = input
        self.output = output
        self.data = data or {}
        self.status = status

        self.tasks = {}
        self.edges = {}
        self.edges_by_task_key = collections.defaultdict(
            lambda: collections.defaultdict(dict))
        self.root_task = None

    def add_tasks(self, tasks=None):
        for task in tasks: self.add_task(task=task)

    def add_task(self, task=None, key=None, input=None,
                 as_root=False, precursor=None, successor=None):
        key = self.determine_key_for_task(task=task, key=key)
        self.tasks[key] = task
        task.key = key
        task.flow = self
        if input is not None: task.input = input
        if as_root: self.root_task = task
        for precursor_task in self._ensure_iterable(precursor):
            self.add_edge(src=precursor_task, dest=task)
        for successor_task in self._ensure_iterable(successor):
            self.add_edge(src=task, dest=successor_task)
        return task

    def determine_key_for_task(self, task=None, key=None):
        existing_key = getattr(task, 'key', None)
        if existing_key:
            if key: raise Exception("Can not reassign task key.")
            else: key = existing_key
        else:
            if key is None: key = str(uuid4())
        return key

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
        src_key, dest_key = [self.ensure_task_key(task) for task in [src, dest]]
        edge_key = self.get_edge_key(src=src_key, dest=dest_key)
        edge = {'src': self.tasks[src_key], 'dest': self.tasks[dest_key]}
        self.edges[edge_key] = edge
        self.edges_by_task_key[src_key]['outgoing'][edge_key] = edge
        self.edges_by_task_key[dest_key]['incoming'][edge_key] = edge

    def get_edge_key(self, src=None, dest=None):
        return (self.ensure_task_key(src), self.ensure_task_key(dest))

    def has_edge(self, src=None, dest=None):
        return self.get_edge_key(src=src, dest=dest) in self.edges

    def ensure_task_key(self, task_or_task_key):
        if isinstance(task_or_task_key, str) or \
           isinstance(task_or_task_key, int):
            task_key = task_or_task_key
        else:
            task_key = task_or_task_key.key
        return task_key

    def get_precursors(self, task=None):
        task_edges = self.edges_by_task_key[task.key]
        precursors = [e['src'] for e in task_edges['incoming'].values()]
        return precursors

    def get_successors(self, task=None):
        task_edges = self.edges_by_task_key[task.key]
        successors = [e['dest'] for e in task_edges['outgoing'].values()]
        return successors

    def get_nearest_pending_tasks(self):
        if not self.root_task: return []
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
        status_filter = lambda task: task.status == status
        return self.filter_tasks(filters=[status_filter])

    def has_incomplete_tasks(self):
        filter_fn = lambda task: task.status != 'COMPLETED'
        incomplete_tasks = self.filter_tasks(filters=[filter_fn])
        return len(incomplete_tasks) > 0

    def get_tail_tasks(self):
        if len(self.tasks) == 1 and self.root_task.key in self.tasks:
            tail_tasks = [self.root_task]
        else:
            tail_tasks = [
                self.tasks[task_key]
                for task_key, task_edges in self.edges_by_task_key.items()
                if len(task_edges['outgoing']) == 0
            ]
        return tail_tasks
