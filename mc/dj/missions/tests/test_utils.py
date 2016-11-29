from django.test import TestCase
from unittest.mock import call, MagicMock, patch, DEFAULT

from ..models import Workflow
from .. import utils

class StartWorkflowTestCase(TestCase):
    def setUp(self):
        self.patchers = {
            'utils': patch.multiple(utils, get_workflow_runner=DEFAULT)
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def test_creates_workflow(self):
        runner_key = 'some runner'
        workflow = utils.start_workflow(runner_key=runner_key)
        self.assertEqual(list(Workflow.objects.all()), [workflow])
        self.assertEqual(workflow.runner_key, runner_key)

    def test_calls_runner_tick(self):
        workflow = utils.start_workflow(runner_key='junk')
        mock_runner = self.mocks['utils']['get_workflow_runner'].return_value
        self.assertEqual(mock_runner.tick.call_args, call(workflow=workflow))

class GetWorkflowRunnerTestCase(TestCase):
    def test_calls_loader(self):
        runner_loader = MagicMock()
        runner_key = 'some runner'
        utils.get_workflow_runner(runner_key=runner_key,
                                  runner_loader=runner_loader)
        self.assertEqual(runner_loader.load_runner.call_args, 
                         call(key=runner_key))
