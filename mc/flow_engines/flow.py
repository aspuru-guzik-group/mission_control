import collections
from uuid import uuid4


class Flow(object):
    ROOT_TASK_KEY = 'ROOT'

    def __init__(self, *args, cfg=None, data=None, label=None, status=None,
                 **kwargs):
        self.cfg = cfg or {'fail_fast': True}
        self.data = data or {}
        self.label = label
        self.status = status or 'PENDING'

        self.tasks = {}
        self.edges = {}
        self.edges_by_key = collections.defaultdict(
            lambda: collections.defaultdict(dict))

        self.add_root_task()

    def to_dict(self):
        return {
            attr: getattr(self, attr, None)
            for attr in ['cfg', 'data', 'label', 'status', 'tasks', 'edges']
        }

    def add_root_task(self):
        self.add_task(task={'key': self.ROOT_TASK_KEY, 'status': 'COMPLETED'})

    def add_task(self, task=None):
        task.setdefault('key', self.generate_task_key())
        task.setdefault('status', 'PENDING')
        self.tasks[task['key']] = task
        for precursor_key in task.get('precursors', []):
            self.add_edge(edge={'src_key': precursor_key,
                                'dest_key': task['key']})
        for successor_key in task.get('successors', []):
            self.add_edge(edge={'src_key': task['key'],
                                'dest_key': successor_key})
        return task

    def generate_task_key(self): return str(uuid4())

    def add_edge(self, edge=None):
        src_key, dest_key = (edge['src_key'], edge['dest_key'])
        if dest_key is self.ROOT_TASK_KEY:
            raise Exception("Root task can not be an edge dest")
        edge_key = (src_key, dest_key)
        self.edges[edge_key] = edge
        self.edges_by_key[src_key]['outgoing'][edge_key] = edge
        self.edges_by_key[dest_key]['incoming'][edge_key] = edge

    def has_edge(self, src_key=None, dest_key=None):
        return (src_key, dest_key) in self.edges

    def get_precursors(self, task=None):
        task_edges = self.edges_by_key[task['key']]
        precursors = [self.tasks[edge['src_key']]
                      for edge in task_edges['incoming'].values()]
        return precursors

    def get_successors(self, task=None):
        task_edges = self.edges_by_key[task['key']]
        successors = [self.tasks[edge['dest_key']] 
                      for edge in task_edges['outgoing'].values()]
        return successors

    def get_nearest_pending_tasks(self):
        nearest_pending_tasks = []
        cursors = [self.tasks[self.ROOT_TASK_KEY]]
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

    def filter_tasks(self, filters=None, include_root_task=False):
        result = []
        for task in self.tasks.values():
            if not include_root_task and task['key'] == self.ROOT_TASK_KEY:
                continue
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
        if len(self.tasks) == 1 and self.ROOT_TASK_KEY in self.tasks:
            tail_tasks = [self.tasks[self.ROOT_TASK_KEY]]
        else:
            tail_tasks = [
                self.tasks[key]
                for key, task_edges in self.edges_by_key.items()
                if len(task_edges['outgoing']) == 0
            ]
        return tail_tasks
