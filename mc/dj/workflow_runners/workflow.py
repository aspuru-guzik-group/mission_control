

class Workflow(object):
    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.state = None
        self.status = None

    def add_nodes(self, nodes=None):
        for node in nodes: self.add_node(node=node)

    def add_node(self, node=None):
        assert hasattr(node, 'id')
        self.nodes[node.id] = node

    def add_edges(self, edges=None):
        for edge in edges: self.add_edge(edge=edge)

    def add_edge(self, edge=None):
        edge_key = (edge['src'].id, edge['dest'].id)
        self.edges[edge_key] = edge
