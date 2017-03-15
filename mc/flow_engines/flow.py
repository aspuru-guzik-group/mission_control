import collections
from uuid import uuid4


class Flow(object):
    def __init__(self, *args, input=None, output=None, data=None, status=None,
                 **kwargs):
        self.input = input
        self.output = output
        self.data = data or {}
        self.status = status or 'PENDING'

        self.nodes = {}
        self.edges = {}
        self.edges_by_node_key = collections.defaultdict(
            lambda: collections.defaultdict(dict))
        self.root_node_key = None

    def add_node(self, node=None, key=None, as_root=False,
                 precursor_keys=None, successor_keys=None):
        self._prepare_node(node=node, key=key)
        self.nodes[node['key']] = node
        if as_root: self.root_node_key = node['key']
        for precursor_key in (precursor_keys or []):
            self.add_edge(edge={'src_key': precursor_key,
                                'dest_key': node['key']})
        for successor_key in (successor_keys or []):
            self.add_edge(edge={'src_key': node['key'],
                                'dest_key': successor_key})
        return node

    def _prepare_node(self, node=None, key=None):
        self._ensure_node_key(node=node, key=key)
        self._ensure_node_status(node=node)

    def _ensure_node_key(self, node=None, key=None):
        if key is not None: node['key'] = key
        elif node.get('key', None) is None:
            node['key'] = self.generate_node_key()

    def generate_node_key(self):
        return str(uuid4())

    def _ensure_node_status(self, node=None):
        if 'status' not in node: node['status'] = 'PENDING'

    def add_edge(self, edge=None):
        src_key, dest_key = (edge['src_key'], edge['dest_key'])
        if dest_key is self.root_node_key:
            raise Exception("Root node can not be an edge dest")
        edge_key = (src_key, dest_key)
        self.edges[edge_key] = edge
        self.edges_by_node_key[src_key]['outgoing'][edge_key] = edge
        self.edges_by_node_key[dest_key]['incoming'][edge_key] = edge

    def has_edge(self, src_key=None, dest_key=None):
        return (src_key, dest_key) in self.edges

    def get_precursors(self, node=None):
        node_edges = self.edges_by_node_key[node['key']]
        precursors = [self.nodes[edge['src_key']]
                      for edge in node_edges['incoming'].values()]
        return precursors

    def get_successors(self, node=None):
        node_edges = self.edges_by_node_key[node['key']]
        successors = [self.nodes[edge['dest_key']] 
                      for edge in node_edges['outgoing'].values()]
        return successors

    def get_nearest_pending_nodes(self):
        if not self.root_node_key: return []
        nearest_pending_nodes = []
        cursors = [self.nodes[self.root_node_key]]
        while len(cursors) > 0:
            next_cursors = []
            for cursor in cursors:
                if cursor['status'] == 'PENDING':
                    nearest_pending_nodes.append(cursor)
                elif cursor['status'] == 'COMPLETED':
                    successors = self.get_successors(node=cursor)
                    next_cursors.extend(successors)
            cursors = next_cursors
        return nearest_pending_nodes

    def filter_nodes(self, filters=None):
        result = []
        for node in self.nodes.values():
            passes_filters = True
            for _filter in filters:
                if not _filter(node):
                    passes_filters = False
                    break
            if passes_filters:
                result.append(node)
        return result

    def get_nodes_by_status(self, status=None):
        status_filter = lambda node: node['status'] == status
        return self.filter_nodes(filters=[status_filter])

    def has_incomplete_nodes(self):
        dead_statuses = ['COMPLETED', 'FAILED']
        filter_fn = lambda node: node['status'] not in dead_statuses
        incomplete_nodes = self.filter_nodes(filters=[filter_fn])
        return len(incomplete_nodes) > 0

    def get_tail_nodes(self):
        if len(self.nodes) == 1 and self.root_node_key in self.nodes:
            tail_nodes = [self.nodes[self.root_node_key]]
        else:
            tail_nodes = [
                self.nodes[node_key]
                for node_key, node_edges in self.edges_by_node_key.items()
                if len(node_edges['outgoing']) == 0
            ]
        return tail_nodes
