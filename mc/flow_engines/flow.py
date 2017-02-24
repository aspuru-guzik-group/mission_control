import collections
from uuid import uuid4


class Flow(object):
    def __init__(self, *args, input=None, output=None, data=None, status=None,
                 **kwargs):
        self.input = input
        self.output = output
        self.data = data or {}
        self.status = status or 'PENDING'

        self.tasks = {}
        self.edges = {}
        self.edges_by_task_key = collections.defaultdict(
            lambda: collections.defaultdict(dict))
        self.root_task_key = None

    def add_task(self, task=None, key=None, as_root=False, precursor_keys=None,
                 successor_keys=None):
        self._prepare_task(task=task, key=key)
        self.tasks[task['key']] = task
        if as_root: self.root_task_key = task['key']
        for precursor_key in (precursor_keys or []):
            self.add_edge(edge={'src_key': precursor_key,
                                'dest_key': task['key']})
        for successor_key in (successor_keys or []):
            self.add_edge(edge={'src_key': task['key'],
                                'dest_key': successor_key})
        return task

    def _prepare_task(self, task=None, key=None):
        self._ensure_task_key(task=task, key=key)
        self._ensure_task_status(task=task)

    def _ensure_task_key(self, task=None, key=None):
        if key is not None: task['key'] = key
        elif task.get('key', None) is None:
            task['key'] = self.generate_task_key()

    def generate_task_key(self):
        return str(uuid4())

    def _ensure_task_status(self, task=None):
        if 'status' not in task: task['status'] = 'PENDING'

    def add_edge(self, edge=None):
        src_key, dest_key = (edge['src_key'], edge['dest_key'])
        if dest_key is self.root_task_key:
            raise Exception("Root task can not be an edge dest")
        edge_key = (src_key, dest_key)
        self.edges[edge_key] = edge
        self.edges_by_task_key[src_key]['outgoing'][edge_key] = edge
        self.edges_by_task_key[dest_key]['incoming'][edge_key] = edge

    def has_edge(self, src_key=None, dest_key=None):
        return (src_key, dest_key) in self.edges

    def get_precursors(self, task=None):
        task_edges = self.edges_by_task_key[task['key']]
        precursors = [self.tasks[edge['src_key']]
                      for edge in task_edges['incoming'].values()]
        return precursors

    def get_successors(self, task=None):
        task_edges = self.edges_by_task_key[task['key']]
        successors = [self.tasks[edge['dest_key']] 
                      for edge in task_edges['outgoing'].values()]
        return successors

    def get_nearest_pending_tasks(self):
        if not self.root_task_key: return []
        nearest_pending_tasks = []
        cursors = [self.tasks[self.root_task_key]]
        while len(cursors) > 0:
            next_cursors = []
            for cursor in cursors:
                if cursor['status'] == 'PENDING':
                    nearest_pending_tasks.append(cursor)
                elif cursor['status'] == 'COMPLETED':
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
        status_filter = lambda task: task['status'] == status
        return self.filter_tasks(filters=[status_filter])

    def has_incomplete_tasks(self):
        dead_statuses = ['COMPLETED', 'FAILED']
        filter_fn = lambda task: task['status'] not in dead_statuses
        incomplete_tasks = self.filter_tasks(filters=[filter_fn])
        return len(incomplete_tasks) > 0

    def get_tail_tasks(self):
        if len(self.tasks) == 1 and self.root_task_key in self.tasks:
            tail_tasks = [self.tasks[self.root_task_key]]
        else:
            tail_tasks = [
                self.tasks[task_key]
                for task_key, task_edges in self.edges_by_task_key.items()
                if len(task_edges['outgoing']) == 0
            ]
        return tail_tasks
