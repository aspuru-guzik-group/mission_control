import collections
from uuid import uuid4

class BaseNode(object):
    def __init__(self, id=None, workflow=None, status=None, state=None,
                 **kwargs):
        self.id = id or str(uuid4())
        self.workflow = workflow
        self.status = status
        self.state = state

    @property
    def node_type(self): return self.__class__.__name__

class Workflow(object):
    def __init__(self, root_node=None):
        self.nodes = {}
        self.edges = {}
        self.edges_by_node_id = collections.defaultdict(dict)
        self.state = None
        self.status = None
        self.root_node = root_node

    def add_nodes(self, nodes=None):
        for node in nodes: self.add_node(node=node)

    def add_node(self, node=None, as_root=False, precursor=None,
                 successor=None):
        assert hasattr(node, 'id')
        self.nodes[node.id] = node
        node.workflow = self
        if as_root: self.root_node = node
        for precursor_node in self._ensure_iterable(precursor):
            self.add_edge(src=precursor_node, dest=node)
        for successor_node in self._ensure_iterable(successor):
            self.add_edge(src=node, dest=successor_node)
        return node

    def _ensure_iterable(self, obj=None):
        if obj is None: return []
        elif isinstance(obj, str): return [obj]
        elif isinstance(obj, collections.Iterable): return obj
        else: return [obj]

    def add_edges(self, edges=None):
        for edge in edges: self.add_edge(**edge)

    def add_edge(self, src=None, dest=None):
        if dest is self.root_node:
            raise Exception("Root node can not be an edge dest")
        src_id, dest_id = [self.ensure_node_id(node) for node in [src, dest]]
        edge_key = self.get_edge_key(src=src_id, dest=dest_id)
        edge = {'src': self.nodes[src_id], 'dest': self.nodes[dest_id]}
        self.edges[edge_key] = edge
        self.edges_by_node_id[src_id][edge_key] = edge
        self.edges_by_node_id[dest_id][edge_key] = edge

    def get_edge_key(self, src=None, dest=None):
        return (self.ensure_node_id(src), self.ensure_node_id(dest))

    def has_edge(self, src=None, dest=None):
        return self.get_edge_key(src=src, dest=dest) in self.edges

    def ensure_node_id(self, node_or_node_id):
        if isinstance(node_or_node_id, str) or isinstance(node_or_node_id, int):
            node_id = node_or_node_id
        else:
            node_id = node_or_node_id.id
        return node_id

    def add_child_nodes(self, parent_node=None, child_nodes=None):
        self.add_nodes(nodes=child_nodes)
        for child_node in child_nodes:
            self.add_edge(src=parent_node, dest=child_node)

    def get_child_nodes(self, parent_node=None):
        parent_node_edges = self.edges_by_node_id[parent_node.id].values()
        child_nodes = [edge['dest'] for edge in parent_node_edges
                       if edge['src'] is parent_node]
        return child_nodes

    def get_nearest_pending_nodes(self):
        nearest_pending_nodes = []
        cursors = [self.root_node]
        while len(cursors) > 0:
            next_cursors = []
            for cursor in cursors:
                if cursor.status == 'PENDING':
                    nearest_pending_nodes.append(cursor)
                elif cursor.status == 'COMPLETED':
                    child_nodes = self.get_child_nodes(parent_node=cursor)
                    next_cursors.extend(child_nodes)
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
        state_filter = lambda node: node.status == status
        return self.filter_nodes(filters=[state_filter])

    def has_incomplete_nodes(self):
        filter_fn = lambda node: node.status != 'COMPLETED'
        incomplete_nodes = self.filter_nodes(filters=[filter_fn])
        return len(incomplete_nodes) > 0

