import collections
from uuid import uuid4


class Flow(object):
    ROOT_TASK_KEY = 'ROOT'

    SIMPLE_ATTRS = {'key', 'cfg', 'data', 'label', 'status', 'depth'}

    @staticmethod
    def sanitize_flow_kwargs(flow_kwargs):
        flow_kwargs = flow_kwargs or {}
        special_keys = {'tasks', 'edges', 'graph'}
        return {k: v for k, v in flow_kwargs.items() if k not in special_keys}

    def __init__(self, key=None, cfg=None, data=None, label=None,
                 status=None, depth=0, add_root_task=True, **kwargs):
        self.key = key
        self.cfg = cfg or {'fail_fast': True}
        self.data = data or {}
        self.label = label
        self.status = status or 'PENDING'
        self.depth = depth or 0
        
        self.tasks = {}
        self._edges = {}
        self._edges_by_key = collections.defaultdict(
            lambda: collections.defaultdict(dict))

        if add_root_task: self.add_root_task()

    def add_root_task(self):
        self.add_task(task={'key': self.ROOT_TASK_KEY, 'status': 'COMPLETED'})

    def add_task(self, task=None):
        task.setdefault('key', self.generate_key())
        task.setdefault('status', 'PENDING')
        self.tasks[task['key']] = task
        for precursor_key in task.get('precursors', []):
            self.add_edge(edge={'src_key': precursor_key,
                                'dest_key': task['key']})
        for successor_key in task.get('successors', []):
            self.add_edge(edge={'src_key': task['key'],
                                'dest_key': successor_key})
        return task

    def generate_key(self): return str(uuid4())

    def add_edge(self, edge=None):
        src_key, dest_key = (edge['src_key'], edge['dest_key'])
        if dest_key is self.ROOT_TASK_KEY:
            raise Exception("Root task can not be an edge dest")
        edge_key = (src_key, dest_key)
        self._edges[edge_key] = edge
        self._edges_by_key[src_key]['outgoing'][edge_key] = edge
        self._edges_by_key[dest_key]['incoming'][edge_key] = edge

    @classmethod
    def from_flow_dict(cls, flow_dict=None, **kwargs):
        flow = cls(**cls.sanitize_flow_kwargs(flow_dict),
                   add_root_task=False)
        graph = flow_dict.get('graph', {})
        for task in graph.get('tasks', {}).values(): flow.add_task(task=task)
        for edge in graph.get('edges', []): flow.add_edge(edge=edge)
        return flow

    def to_flow_dict(self):
        flow_dict = {
            **{attr: getattr(self, attr, None) for attr in self.SIMPLE_ATTRS},
            'graph': {
                'tasks': {key: task for key, task in self.tasks.items()},
                'edges': [edge for edge in self._edges.values()],
            }
        }
        if 'key' in flow_dict and flow_dict['key'] is None: del flow_dict['key']
        return flow_dict

    def has_edge(self, src_key=None, dest_key=None):
        return (src_key, dest_key) in self._edges

    def get_precursors(self, task=None):
        task_edges = self._edges_by_key[task['key']]
        precursors = [self.tasks[edge['src_key']]
                      for edge in task_edges['incoming'].values()]
        return precursors

    def get_successors(self, task=None):
        task_edges = self._edges_by_key[task['key']]
        successors = [self.tasks[edge['dest_key']] 
                      for edge in task_edges['outgoing'].values()]
        return successors

    def get_nearest_tickable_pending_tasks(self):
        nearest_tickable_pending_tasks = []
        cursors = [self.tasks[self.ROOT_TASK_KEY]]
        while len(cursors) > 0:
            next_cursors = []
            for cursor in cursors:
                if cursor['status'] == 'PENDING':
                    nearest_tickable_pending_tasks.append(cursor)
                elif cursor['status'] == 'COMPLETED':
                    successors = self.get_successors(task=cursor)
                    next_cursors.extend(successors)
            cursors = next_cursors
        return nearest_tickable_pending_tasks

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
                for key, task_edges in self._edges_by_key.items()
                if len(task_edges['outgoing']) == 0
            ]
        return tail_tasks

    def get_tickable_tasks(self):
        return (self.get_running_tasks() +
                self.get_nearest_tickable_pending_tasks())

    def get_running_tasks(self):
        return self.get_tasks_by_status(status='RUNNING')
