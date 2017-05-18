from collections import defaultdict
import unittest
from unittest.mock import Mock
from uuid import uuid4

from ..flow import Flow


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.flow = Flow()

    def generate_tasks(self, n=3):
        tasks = []
        for i in range(n):
            tasks.append(self.generate_task(key=i))
        return tasks

    def generate_task(self, key=None, **kwargs):
        if key is None: key = str(uuid4())
        task = {'key': key, **kwargs}
        return task

    def _task_list_to_dict(self, task_list):
        return {task['key']: task for task in task_list}

class AddTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.task = self.generate_task(key='task')

    def test_adds_task(self):
        self.assertFalse(self.task['key'] in self.flow.tasks)
        self.flow.add_task(task=self.task)
        self.assertEqual(self.flow.tasks[self.task['key']], self.task)

    def test_returns_task(self):
        result = self.flow.add_task(task=self.task)
        self.assertEqual(result, self.task)

    def test_generates_key_if_unset(self):
        self.flow.generate_task_key = Mock()
        del self.task['key']
        self.flow.add_task(task=self.task)
        self.assertEqual(self.task['key'],
                         self.flow.generate_task_key.return_value)

    def test_adds_edges_for_precursors(self):
        precursors = [self.flow.add_task(task=self.generate_task(key=i))
                      for i in range(1)]
        self.flow.add_task(
            task=self.task,
            precursors=[precursor['key'] for precursor in precursors])
        for precursor in precursors:
            self.assertTrue(self.flow.has_edge(src_key=precursor['key'],
                                               dest_key=self.task['key']))

    def test_adds_edges_for_successor(self):
        successors = [self.flow.add_task(task=self.generate_task(key=i))
                      for i in range(3)]
        self.flow.add_task(
            task=self.task,
            successors=[successor['key'] for successor in successors])
        for successor in successors:
            self.assertTrue(self.flow.has_edge(src_key=self.task['key'],
                                               dest_key=successor['key']))

class AddEdgeTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.tasks = [self.flow.add_task(task=task)
                      for task in self.generate_tasks(n=2)]
        self.edge = {'src_key': self.tasks[0]['key'],
                     'dest_key': self.tasks[1]['key']}
        self.expected_edge_key = (self.edge['src_key'], self.edge['dest_key'])

    def test_add_edge(self):
        self.assertFalse(self.expected_edge_key in self.flow.edges)
        self.flow.add_edge(edge=self.edge)
        self.assertEqual(self.flow.edges[self.expected_edge_key], self.edge)

    def test_updates_edges_by_key(self):
        self.assertTrue(self.tasks[0]['key'] not in self.flow.edges_by_key)
        self.assertTrue(self.tasks[1]['key'] not in self.flow.edges_by_key)
        self.flow.add_edge(edge=self.edge)
        self.assertEqual(
            self.flow.edges_by_key[self.edge['src_key']]['outgoing'],
            {self.expected_edge_key: self.edge})
        self.assertEqual(
            self.flow.edges_by_key[self.edge['dest_key']]['incoming'],
            {self.expected_edge_key: self.edge})

class GetNearestPendingTasksTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow()

    def setup_flow(self):
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
        tail = flow.tasks[flow.ROOT_TASK_KEY]
        for task in branch_tasks:
            tail = flow.add_task(task=task, precursors=[tail['key']])

    def test_get_nearest_pending_tasks(self):
        expected_nearest_pending_tasks = [self.linear_branches['1'][1],
                                          self.linear_branches['2'][2]]
        self.assertEqual(self.flow.get_nearest_pending_tasks(),
                         expected_nearest_pending_tasks)

class GetSuccessorsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.setup_flow(flow=self.flow)

    def setup_flow(self, flow=None):
        task_1 = self.generate_task(key='1')
        self.flow.add_task(task=task_1, precursors=[flow.ROOT_TASK_KEY])
        task_1_1 = self.generate_task(key='1_1')
        self.flow.add_task(task=task_1_1, precursors=[task_1['key']])
        task_1_2 = self.generate_task(key='1_2')
        self.flow.add_task(task=task_1_2, precursors=[task_1['key']])

    def test_gets_successors(self):
        successors = self.flow.get_successors(task=self.flow.tasks['1'])
        expected_successors = [self.flow.tasks['1_1'], self.flow.tasks['1_2']]
        self.assertEqual(self._task_list_to_dict(successors),
                         self._task_list_to_dict(expected_successors))

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
            self.flow.add_task(task=task)
        _filter = lambda task: task['answer'] == 'yes'
        result = self.flow.filter_tasks(filters=[_filter])
        expected_result = tasks['yes']
        self.assertEqual(self._task_list_to_dict(result),
                         self._task_list_to_dict(expected_result))

class GetTasksByStatusTestCase(BaseTestCase):
    def test_gets_tasks_by_status(self):
        tasks = defaultdict(list)
        statuses = ['status_%s' % i for i in range(3)]
        for status in statuses:
            for i in range(3):
                task = self.generate_task(status=status)
                tasks[status].append(task)
                self.flow.add_task(task=task)
        for status in statuses:
            result = self.flow.get_tasks_by_status(status=status)
            expected_result = tasks[status]
            self.assertEqual(self._task_list_to_dict(result),
                             self._task_list_to_dict(expected_result))

class HasIncompleteTasksTestCase(BaseTestCase):
    def test_has_incomplete(self):
        self.flow.add_task(task=self.generate_task(status='PENDING'))
        self.assertTrue(self.flow.has_incomplete_tasks())

    def test_does_not_have_incomplete(self):
        self.flow.add_task(task=self.generate_task(status='COMPLETED'))
        self.assertFalse(self.flow.has_incomplete_tasks())

class GetTailTasksTestCase(BaseTestCase):
    def test_gets_tail_tasks(self):
        tail_tasks = [
            self.flow.add_task(task=self.generate_task(),
                               precursors=[self.flow.ROOT_TASK_KEY])
            for i in range(3)]
        self.assertTrue(self.flow.get_tail_tasks(), tail_tasks)

    def test_handles_flow_with_only_root_task(self):
        self.assertTrue(self.flow.get_tail_tasks(),
                        self.flow.tasks[self.flow.ROOT_TASK_KEY])

if __name__ == '__main__':
    unittest.main()
