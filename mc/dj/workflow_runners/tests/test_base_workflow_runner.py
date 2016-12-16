import unittest
#from unittest.mock import call, DEFAULT, MagicMock, patch

from ..base_workflow_runner import BaseWorkflowRunner


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.runner = BaseWorkflowRunner()

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start()
                 for key, patcher in patchers.items()}
        return mocks

class DeserializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_components()

    def setup_components(self):
        class BaseComponent(object):
            def __init__(self, *args, id=None, **kwargs):
                self.id = id

        self.components = {name: type(name, (BaseComponent,), {})
                           for name in ['Type1', 'Type2', 'Type3']}
        for name, component in self.components.items():
            def test_fn(serialized_node, *args, _name=name, **kwargs):
                return serialized_node['type'] == _name
            self.runner.register_component(test_fn=test_fn, component=component)

    def test_deserializes_workflow(self):
        serialized_workflow = {
            'nodes': [
                {'id': 'a', 'type': 'Type1', 'status': 'COMPLETED', 'data': 'a',
                 'links_to': ['b']},
                {'id': 'b', 'type': 'Type2', 'status': 'COMPLETED', 'data': 'b',
                 'links_to': ['c', 'd']},
                {'id': 'c', 'type': 'Type1', 'status': 'RUNNING', 'data': 'c'},
                {'id': 'd', 'type' : 'Type2', 'status': 'RUNNING', 'data': 'd'},
            ]
        }
        workflow = self.runner.deserialize_workflow(serialized_workflow)
        expected_nodes = {}
        for serialized_node in serialized_workflow['nodes']:
            expected_node = self.runner.deserialize_node(serialized_node)
            expected_nodes[expected_node.id] = expected_node

        def summarize_nodes(nodes):
            return {node.id: {'type': type(node), '__dict__': node.__dict__}
                    for node in nodes.values()}
        self.assertEqual(summarize_nodes(workflow.nodes),
                         summarize_nodes(expected_nodes))
        expected_edges = {}
        for serialized_node in serialized_workflow['nodes']:
            node_id = serialized_node['id']
            for link_target_id in serialized_node.get('links_to', []):
                expected_edges[(node_id, link_target_id)] = { 
                    'src': workflow.nodes[node_id],
                    'dest': workflow.nodes[link_target_id]
                }
        self.assertEqual(workflow.edges, expected_edges)

if __name__ == '__main__':
    unittest.main()
