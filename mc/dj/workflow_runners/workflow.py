from collections import defaultdict
from uuid import uuid4


class StaticNode(object):
    def __init__(self, id=None, status=None, state=None, **kwargs):
        self.id = id or str(uuid4())
        self.status = status
        self.state = state

class Workflow(object):
    def __init__(self, root_node=None):
        self.nodes = {}
        self.edges = {}
        self.edges_by_node_id = defaultdict(dict)
        self.state = None
        self.status = None
        self.root_node = root_node

    def add_nodes(self, nodes=None):
        for node in nodes: self.add_node(node=node)

    def add_node(self, node=None):
        assert hasattr(node, 'id')
        self.nodes[node.id] = node

    def add_edges(self, edges=None):
        for edge in edges: self.add_edge(edge=edge)

    def add_edge(self, edge=None):
        src, dest = [edge[k] for k in ('src','dest')]
        edge_key = (src.id, dest.id)
        self.edges[edge_key] = edge
        self.edges_by_node_id[src.id][edge_key] = edge
        self.edges_by_node_id[dest.id][edge_key] = edge

    def connect_nodes(self, src=None, dest=None):
        edge = {'src': src, 'dest': dest}
        self.add_edge(edge=edge)

    def add_child_nodes(self, parent_node=None, child_nodes=None):
        self.add_nodes(nodes=child_nodes)
        for child_node in child_nodes:
            self.connect_nodes(src=parent_node, dest=child_node)

    def get_child_nodes(self, parent_node=None):
        parent_node_edges = self.edges_by_node_id[parent_node.id].values()
        child_nodes = [edge['dest'] for edge in parent_node_edges
                       if edge['src'] is parent_node]
        return child_nodes

    def get_nearest_incomplete_nodes(self):
        nearest_incomplete_nodes = []
        cursors = [self.root_node]
        while len(cursors) > 0:
            next_cursors = []
            for cursor in cursors:
                if cursor.status == 'COMPLETED':
                    child_nodes = self.get_child_nodes(parent_node=cursor)
                    next_cursors.extend(child_nodes)
                else:
                    nearest_incomplete_nodes.append(cursor)
            cursors = next_cursors
        return nearest_incomplete_nodes
