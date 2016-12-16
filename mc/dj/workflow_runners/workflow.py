

class Workflow(object):
    def __init__(self):
        self.nodes = {}
        self.edges = {}

    def add_nodes(self, nodes=None):
        for node in nodes: self.add_node(node=node)

    def add_node(self, node=None):
        assert hasattr(node, 'id')
        self.nodes[node.id] = node

    def add_edges_from_specs(self, edge_specs=None):
        for edge_spec in edge_specs:
            self.add_edge_from_spec(edge_spec=edge_spec)

    def add_edge_from_spec(self, edge_spec=None):
        src_id, dest_id = edge_spec['src_id'], edge_spec['dest_id']
        self.edges[(src_id, dest_id)] = {'src': self.nodes[src_id],
                                         'dest': self.nodes[dest_id]}
