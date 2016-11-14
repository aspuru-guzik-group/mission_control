from django.test import TestCase
from django.test.utils import override_settings
from fireworks.core.firework import FireTaskBase

from .urls import urlpatterns
assert urlpatterns
from .utils import task_registry


@override_settings(ROOT_URLCONF=__name__)
class CreateTaskFormTestCase(TestCase):
    def setUp(self):
        self.url_for_create = '/create/'
    def tearDown(self):
        task_registry.clear()

    def test_form_exists(self):
        response = self.client.get(self.url_for_create)
        self.assertEqual(response.status_code, 200)

    def test_lists_registered_tasks(self):
        tasks = [type('Task', (FireTaskBase,), {'_fw_name': 'task_%s' % i})
                 for i in range(3)]
        for task in tasks:
            task_registry.register_task(task)
        response = self.client.get(self.url_for_create)
        for task in tasks:
            self.assertTrue(task._fw_name in response.content.decode())

    def test_create_task(self):
        class TestTask(FireTaskBase):
            _fw_name = 'test_task'
        task_registry.register_task(TestTask, key='test_task')
        self.client.post(self.url_for_create, {'task': 'test_task'})
        self.assertTrue(isinstance(add_wf_call_arg, TestTask))
