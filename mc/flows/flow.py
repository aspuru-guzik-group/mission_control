import collections
from uuid import uuid4


class Flow(object):
    """A flow is modeled as a <directed acyclic graph>.
    (https://en.wikipedia.org/wiki/Directed_acyclic_graph)
    """
    ROOT_TASK_KEY = 'ROOT'

    SIMPLE_ATTRS = {'key', 'cfg', 'data', 'label', 'status', 'depth'}

    @staticmethod
    def sanitize_flow_kwargs(flow_kwargs):
        flow_kwargs = flow_kwargs or {}
        special_keys = {'tasks', 'edges', 'graph'}
        return {k: v for k, v in flow_kwargs.items() if k not in special_keys}

    def __init__(self, key=None, cfg=None, data=None, label=None,
                 status=None, depth=0, add_root_task=True, **kwargs):
        """
        Args:
            key (str): a key to identify the flow.
            cfg [dict]: a dict of flow cfg data. Defaults to:
                {'fail_fast': True}.
            data [dict]: dict of flow data.
            label [str]: flow label
            status [str]: flow status, defaults to 'PENDING'
            depth [int]: depth level, for prioritizing nested flows.
            add_root_task [bool]: whether to add root task. Default: True.
        """
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
        self._last_added_task = None
        if add_root_task: self.add_root_task()

    def add_root_task(self):
        self.add_task(task={'key': self.ROOT_TASK_KEY, 'status': 'COMPLETED'},
                      add_default_connections=False)

    def add_task(self, task=None, add_default_connections=True):
        """
        Args:
            task <dict>: a task dict
            add_default_connections [bool]: whether to add default precursor
                if precursors or successors are not specified. Default: True.
        Returns:
            task <dict>: altered task dict
        """
        task.setdefault('key', self.generate_key())
        task.setdefault('status', 'PENDING')
        self.tasks[task['key']] = task
        if add_default_connections:
            self.add_default_connections_to_task(task=task)
        for precursor_key in task.get('precursors', []):
            self.add_edge(edge={'src_key': precursor_key,
                                'dest_key': task['key']})
        for successor_key in task.get('successors', []):
            self.add_edge(edge={'src_key': task['key'],
                                'dest_key': successor_key})
        self._last_added_task = task
        return task

    def generate_key(self):
        """Generate keys as uuid.uuid4 strings"""
        return str(uuid4())

    def add_default_connections_to_task(self, task=None):
        """
        If task has no precursors or successors given, use last task added
        as precursor.
        """
        if ('precursors' not in task and 'successors' not in task
            and self._last_added_task is not None
            and task['key'] != self.ROOT_TASK_KEY
           ):
            task['precursors'] = [self._last_added_task['key']]

    def add_edge(self, edge=None):
        """
        Also adds entry to _edges_by_key.

        Args:
            edge <dict>: dict in this shape: {'src_key': <>, 'dest_key': <>}
        """
        src_key, dest_key = (edge['src_key'], edge['dest_key'])
        if dest_key is self.ROOT_TASK_KEY:
            raise Exception("Root task can not be an edge dest")
        edge_key = (src_key, dest_key)
        self._edges[edge_key] = edge
        self._edges_by_key[src_key]['outgoing'][edge_key] = edge
        self._edges_by_key[dest_key]['incoming'][edge_key] = edge

    @classmethod
    def from_flow_spec(cls, flow_spec=None):
        """Create flow from flow_spec."""
        flow = Flow(**Flow.sanitize_flow_kwargs(flow_spec))
        for i, task in enumerate(flow_spec.get('tasks', [])):
            flow.add_task(task=task)
        return flow

    @classmethod
    def from_flow_dict(cls, flow_dict=None, **kwargs):
        """Create flow from flow dict."""
        flow = cls(**cls.sanitize_flow_kwargs(flow_dict),
                   add_root_task=False)
        graph = flow_dict.get('graph', {})
        for task in graph.get('tasks', {}).values():
            flow.add_task(task=task, add_default_connections=False)
        for edge in graph.get('edges', []): flow.add_edge(edge=edge)
        if flow.ROOT_TASK_KEY not in flow.tasks: flow.add_root_task()
        return flow

    def to_flow_dict(self):
        """Create flow_dict from flow."""
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
        """
        Tickable pending tasks are tasks with status 'PENDING' which can be
        reached by starting from the flow root and advancing only through tasks
        with status 'COMPLETED'.
        """
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
        """
        Loop through all tasks, returning tasks which pass all of the filter
        functions provided.

        Args:
            filters <list(<fn>)>: list of filter functions
            include_root_task [bool]: whether to include the root task.
                Default: False
        """
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
        """
        Get tasks that are the edge of the graph.
        """
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
        """
        Tickable tasks are nearest tickable tasks + running tasks.
        """
        return (self.get_running_tasks() +
                self.get_nearest_tickable_pending_tasks())

    def get_running_tasks(self):
        return self.get_tasks_by_status(status='RUNNING')
