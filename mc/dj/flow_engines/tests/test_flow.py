from collections import defaultdict
import unittest
from unittest.mock import call, DEFAULT, patch, Mock
from uuid import uuid4

from ..flow import Flow


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow = Flow()

    def tearDown(self):
        if hasattr(self, 'patchers'): self.stop_patchers(self.patchers)

    def stop_patchers(self, patchers):
        for patcher in patchers.values(): patcher.stop()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start() for key, patcher in patchers.items()}
        return mocks

    def generate_tasks(self, n=3):
        tasks = []
        for i in range(n):
            tasks.append(self.generate_task(key=i))
        return tasks

    def generate_task(self, key=None, **kwargs):
        if key is None: key = str(uuid4())
        task = Mock()
        task.configure_mock(key=key, **kwargs)
        return task

class AddTasksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'flow': patch.multiple(self.flow, add_task=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_adds_tasks(self):
        tasks = self.generate_tasks(n=3)
        self.flow.add_tasks(tasks=tasks)
        self.assertEqual(self.mocks['flow']['add_task'].call_args_list,
                         [call(task=task) for task in tasks])

class AddTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task = self.generate_task(key='task')

    def test_adds_task(self):
        self.assertFalse(self.task.key in self.flow.tasks)
        self.flow.add_task(task=self.task)
        self.assertEqual(self.flow.tasks[self.task.key], self.task)

    def test_returns_task(self):
        result = self.flow.add_task(task=self.task)
        self.assertEqual(result, self.task)

    def test_sets_flow_on_task(self):
        self.flow.add_task(task=self.task)
        self.assertEqual(self.task.flow, self.flow)

    def test_sets_key_on_task_wout_key(self):
        self.task.key = None
        key = 'some key'
        self.flow.add_task(task=self.task, key=key)
        self.assertEqual(self.task.key, key)

    def test_errors_on_key_reassignment(self):
        with self.assertRaises(Exception):
            self.flow.add_task(task=self.task, key='some key')

    def test_uses_task_key_if_no_key_provided(self):
        original_key = self.task.key
        self.flow.add_task(task=self.task)
        self.assertEqual(self.task.key, original_key)

    def test_sets_key_to_uuid_as_fallback(self):
        self.task.key = None
        self.flow.add_task(task=self.task)
        self.assertTrue(self.task.key is not None)

    def test_add_as_root(self):
        self.flow.add_task(task=self.task, as_root=True)
        self.assertEqual(self.task.flow, self.flow)
        self.assertEqual(self.flow.root_task, self.task)

    def test_add_with_single_precursor(self):
        precursor = self.flow.add_task(task=self.generate_task())
        self.flow.add_task(task=self.task, precursor=precursor)
        self.assertTrue(self.flow.has_edge(src=precursor, dest=self.task))

    def test_add_with_precursor_key(self):
        precursor = self.generate_task(key='precursor_key')
        self.flow.add_task(task=precursor)
        self.flow.add_task(task=self.task, precursor='precursor_key')
        self.assertTrue(self.flow.has_edge(src=precursor, dest=self.task))

    def test_add_with_many_precursors(self):
        precursors = [self.flow.add_task(task=self.generate_task(key=i))
                      for i in range(1)]
        self.flow.add_task(task=self.task, precursor=precursors)
        for precursor in precursors:
            self.assertTrue(
                self.flow.has_edge(src=precursor, dest=self.task))

    def test_add_with_single_successor(self):
        successor = self.flow.add_task(task=self.generate_task())
        self.flow.add_task(task=self.task, successor=successor)
        self.assertTrue(self.flow.has_edge(src=self.task, dest=successor))

    def test_add_with_successsor_id(self):
        successor = self.generate_task(key='successor_id')
        self.flow.add_task(task=successor)
        self.flow.add_task(task=self.task, successor='successor_id')
        self.assertTrue(self.flow.has_edge(src=self.task, dest=successor))

    def test_add_with_many_successors(self):
        successors = [self.flow.add_task(task=self.generate_task(key=i))
                      for i in range(3)]
        self.flow.add_task(task=self.task, successor=successors)
        for successor in successors:
            self.assertTrue(
                self.flow.has_edge(src=self.task, dest=successor))

class AddEdgesTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'flow': patch.multiple(self.flow, add_edge=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)

    def test_adds_edges(self):
        edges = [{'src': i, 'dest': i} for i in range(3)]
        self.flow.add_edges(edges=edges)
        self.assertEqual(self.mocks['flow']['add_edge'].call_args_list,
                         [call(**edge) for edge in edges])

class AddEdgeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.tasks = self.generate_tasks(n=2)
        self.flow.add_tasks(tasks=self.tasks)
        self.edge = {'src': self.tasks[0], 'dest': self.tasks[1]}
        self.expected_edge_key = (self.edge['src'].key, self.edge['dest'].key)

    def test_add_edge(self):
        self.assertFalse(self.expected_edge_key in self.flow.edges)
        self.flow.add_edge(**self.edge)
        self.assertEqual(self.flow.edges[self.expected_edge_key], self.edge)

    def test_updates_edges_by_task_key(self):
        self.assertTrue(self.tasks[0].key not in self.flow.edges_by_task_key)
        self.assertTrue(self.tasks[1].key not in self.flow.edges_by_task_key)
        self.flow.add_edge(**self.edge)
        for task in self.tasks:
            self.assertEqual(
                self.flow.edges_by_task_key[task.key][self.expected_edge_key],
                self.edge)

class GetNearestPendingTasksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow()

    def setup_flow(self):
        self.flow.add_task(
            self.generate_task(key='ROOT', status='COMPLETED'), as_root=True)
        self.linear_branches = {
            '1': [
                self.generate_task(key='1.1', status='COMPLETED'),
                self.generate_task(key='1.2', status='PENDING'),
                self.generate_task(key='1.3', status='PENDING')
            ],
            '2': [
                self.generate_task(key='2.1', status='COMPLETED'),
                self.generate_task(key='2.2', status='COMPLETED'),
                self.generate_task(key='2.3', status='PENDING'),
                self.generate_task(key='2.4', status='PENDING'),
            ],
        }
        for branch_tasks in self.linear_branches.values():
            self.add_linear_branch_to_flow(flow=self.flow,
                                               branch_tasks=branch_tasks)

    def add_linear_branch_to_flow(self, flow=None, branch_tasks=None):
        tail = flow.root_task
        for task in branch_tasks:
            tail = flow.add_task(task=task, precursor=[tail])

    def test_get_nearest_pending_tasks(self):
        expected_nearest_pending_tasks = [self.linear_branches['1'][1],
                                          self.linear_branches['2'][2]]
        self.assertEqual(set(self.flow.get_nearest_pending_tasks()),
                         set(expected_nearest_pending_tasks))

class GetSuccessorsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow(flow=self.flow)

    def setup_flow(self, flow=None):
        self.flow.add_task(self.generate_task(key='ROOT'), as_root=True)
        task_1 = self.generate_task(key='1')
        self.flow.add_task(task_1, precursor=flow.root_task)
        task_1_1 = self.generate_task(key='1_1')
        self.flow.add_task(task_1_1, precursor=task_1)
        task_1_2 = self.generate_task(key='1_2')
        self.flow.add_task(task_1_2, precursor=task_1)

    def test_gets_successors(self):
        successors = self.flow.get_successors(task=self.flow.tasks['1'])
        expected_successors = [self.flow.tasks['1_1'], self.flow.tasks['1_2']]
        self.assertEqual(set(successors), set(expected_successors))

class FilterTasksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()

    def test_filters(self):
        tasks = defaultdict(list)
        for i in range(4):
            if i % 2: answer = 'yes'
            else: answer = 'no'
            task = self.generate_task(answer=answer)
            tasks[answer].append(task)
            self.flow.add_task(task)
        _filter = lambda task: task.answer == 'yes'
        result = self.flow.filter_tasks(filters=[_filter])
        expected_result = tasks['yes']
        self.assertEqual(set(result), set(expected_result))

class GetTasksByStatusTestCase(BaseTestCase):
    def test_gets_tasks_by_status(self):
        tasks = defaultdict(list)
        statuses = ['status_%s' % i for i in range(3)]
        for status in statuses:
            for i in range(3):
                task = self.generate_task(status=status)
                tasks[status].append(task)
                self.flow.add_task(task)
        for status in statuses:
            result = self.flow.get_tasks_by_status(status=status)
            expected_result = tasks[status]
            self.assertEqual(set(result), set(expected_result))

class HasIncompleteTasksTestCase(BaseTestCase):
    def test_has_incomplete(self):
        self.flow.add_task(self.generate_task(status='PENDING'))
        self.assertTrue(self.flow.has_incomplete_tasks())

    def test_does_not_have_incomplete(self):
        self.flow.add_task(self.generate_task(status='COMPLETED'))
        self.assertFalse(self.flow.has_incomplete_tasks())


if __name__ == '__main__':
    unittest.main()
