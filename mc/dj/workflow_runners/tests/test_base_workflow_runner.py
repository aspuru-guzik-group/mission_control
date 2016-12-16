import unittest
#from unittest.mock import call, DEFAULT, MagicMock, patch

from ..base_workflow_runner import BaseWorkflowRunner
from ..workflow import Workflow


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
        class BaseNode(object):
            def __init__(self, *args, id=None, **kwargs):
                self.id = id

        self.node_classes = {name: type(name, (BaseNode,), {})
                             for name in ['Type1', 'Type2', 'Type3']}
        for name, node_class in self.node_classes.items():
            def test_fn(serialized_node, *args, _name=name, **kwargs):
                return serialized_node['type'] == _name
            self.runner.register_node_class(test_fn=test_fn,
                                            node_class=node_class)

    def test_deserializes_workflow(self):
        serialized_workflow = {
            'nodes': [
                {'id': 'a', 'type': 'Type1', 'status': 'COMPLETED'},
                {'id': 'b', 'type': 'Type2', 'status': 'COMPLETED'},
                {'id': 'c', 'type': 'Type1', 'status': 'RUNNING'},
                {'id': 'd', 'type' : 'Type2', 'status': 'RUNNING'},
            ],
            'edges': [
                {'src_id': 'a', 'dest_id': 'b'},
                {'src_id': 'b', 'dest_id': 'c'},
                {'src_id': 'b', 'dest_id': 'd'},
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
        for serialized_edge in serialized_workflow['edges']:
            src_id, dest_id = (serialized_edge['src_id'],
                               serialized_edge['dest_id'])
            expected_edges[(src_id, dest_id)] = {
                'src': workflow.nodes[src_id],
                'dest': workflow.nodes[dest_id]
            }
        self.assertEqual(workflow.edges, expected_edges)


class SerializationTestCase(BaseTestCase):
    def test_serializes_workflow(self):
        workflow = self.generate_workflow()
        serialization = self.runner.serialize_workflow(workflow)
        expected_serialization = {
            'nodes': [{'id': node.id, 'type': node.type, 'status': node.status,
                       'state': node.state}
                      for node in workflow.nodes.values()],
            'edges': [{'src_id': edge['src'].id, 'dest_id': edge['dest'].id}
                      for edge in workflow.edges.values()],
            'state': workflow.state,
            'status': workflow.status,
        }
        self.assertEqual(serialization, expected_serialization)

    def generate_workflow(self):
        workflow = Workflow()
        class BaseComponent(object):
            def __init__(self, id=None, type=None, status=None, state=None,
                         **kwargs):
                self.id = id
                self.type = type
                self.status = status
                self.state = state

        nodes = [BaseComponent(id=i, type='type_%s' % i, status='status_%s' % i,
                               state=i)
                 for i in range(3)]
        workflow.add_nodes(nodes=nodes)
        edges = [{'src': nodes[0], 'dest': nodes[1]},
                 {'src': nodes[1], 'dest': nodes[2]}]
        workflow.add_edges(edges=edges)
        workflow.state = 'workflow_state'
        workflow.status = 'workflow_status'
        return workflow

if __name__ == '__main__':
    unittest.main()
