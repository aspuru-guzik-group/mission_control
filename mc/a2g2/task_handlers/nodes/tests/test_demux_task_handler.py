from collections import defaultdict
import unittest
from unittest.mock import call, patch, MagicMock
from uuid import uuid4

from .. import demux_task_handler


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.node = defaultdict(MagicMock)
        self.task_params = defaultdict(MagicMock)
        self.task = defaultdict(MagicMock, **{'task_params': self.task_params})
        self.flow_ctx = {
            'create_flow': MagicMock(),
            'get_flow': MagicMock()
        }
        self.task_context = self.generate_task_context()
        self.task_handler = demux_task_handler.DemuxTaskHandler()

    def generate_flow(self, uuid=None, status='PENDING', **flow_state):
        flow = defaultdict(MagicMock)
        if not uuid: uuid = str(uuid4())
        flow.update({'uuid': uuid, 'status': status, **flow_state})
        return flow

    def generate_task_context(self, **kwargs):
        task_context = {
            'task': self.task,
            'node': self.node,
            'flow_ctx': self.flow_ctx,
            **kwargs
        }
        return task_context

class InitialTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.node['node_tasks'] = []
        self.task_handler.generate_flow_task = MagicMock()
        self.task_handler.initial_tick(task=self.task,
                                       task_context=self.task_context)

    def test_initial_tick_adds_flow_task(self):
        self.assertEqual(self.task_handler.generate_flow_task.call_args,
                         call(task=self.task, task_context=self.task_context))
        self.assertEqual(self.node['node_tasks'][-1],
                         self.task_handler.generate_flow_task.return_value)

    def test_completes_task(self):
        self.assertEqual(self.task['status'], 'COMPLETED')

class GenerateFlowTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task_handler.generate_flow_spec = MagicMock()
        self.flow_task = self.task_handler.generate_flow_task(
            task=self.task, task_context=self.task_context)

    def test_has_expected_task_type(self):
        self.assertEqual(self.flow_task['task_type'], 'a2g2:task:run_flow')

    def test_has_expected_task_params(self):
        self.assertEqual(self.task_handler.generate_flow_spec.call_args,
                         call(task=self.task, task_context=self.task_context))
        expected_task_params = {
            'flow_spec': self.task_handler.generate_flow_spec.return_value,
        }
        self.assertEqual(self.flow_task['task_params'], expected_task_params)

class GenerateFlowSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task['task_params']['items'] = [MagicMock() for i in range(3)]
        self.task_handler.generate_node_spec = MagicMock()
        self.flow_spec = self.task_handler.generate_flow_spec(
            task=self.task, task_context=self.task_context)

    def test_makes_expected_generate_node_spec_calls(self):
        expected_call_args_list = [
            call(item=item, index=i, task=self.task,
                 task_context=self.task_context)
            for i, item in enumerate(self.task['task_params'].get('items'))
        ]
        self.assertEqual(self.task_handler.generate_node_spec.call_args_list,
                         expected_call_args_list)

    def test_has_expected_node_specs(self):
        expected_node_specs = [
            self.task_handler.generate_node_spec.return_value
            for item in self.task['task_params'].get('items', [])
        ]
        self.assertEqual(self.flow_spec['node_specs'], expected_node_specs)


@patch.object(demux_task_handler, 'deepcopy')
@patch.object(demux_task_handler, 'set_context_values')
class GenerateNodeSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.substitutions = [MagicMock() for i in range(3)]
        self.task['task_params']['substitutions'] = self.substitutions
        self.item = MagicMock()
        self.index = 1

    def test_calls_set_context_values(self, mock_set_context_values,
                                      mock_deepcopy):
        self.task_handler.generate_node_spec(
            item=self.item, index=self.index, task=self.task,
            task_context=self.task_context)
        self.assertEqual(
            mock_deepcopy.call_args,
            call(self.task['task_params']['node_spec_template_obj']))
        expected_context = {
            'item': self.item,
            'index': self.index,
            'node_spec': mock_deepcopy.return_value,
        }
        self.assertEqual(
            mock_set_context_values.call_args,
            call(value_specs=self.task['task_params']['substitutions'],
                 context=expected_context))

    def test_returns_rendered_node_spec(self, mock_set_context_values,
                                        mock_deepcopy):
        node_spec = self.task_handler.generate_node_spec(
            item=self.item, index=self.index, task=self.task,
            task_context=self.task_context)
        self.assertEqual(node_spec, mock_deepcopy.return_value)

if __name__ == '__main__':
    unittest.main()
