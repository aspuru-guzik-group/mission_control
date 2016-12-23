import json
from django.conf.urls import url, include
from django.test import TestCase, override_settings

from missions.models import Workflow as WorkflowModel
from workflow_engines.workflow import Workflow, BaseNode
from workflow_runners.base_workflow_runner import BaseWorkflowRunner
from workflow_engines.workflow_engine import WorkflowEngine
from workflow_client.workflow_client import MissionControlWorkflowClient


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('missions.urls')),
]

@override_settings(ROOT_URLCONF=__name__)
class WorkflowRunnerE2ETestCase(TestCase):
    def setUp(self):
        self.configure_request_client()
        self.keyed_node_classes = self.generate_keyed_node_classes()
        self.workflow_engine = self.generate_workflow_engine(
            keyed_node_classes=self.keyed_node_classes)
        self.workflow_client = MissionControlWorkflowClient(
            base_url='/%s' % BASE_PATH,
            request_client=self.client
        )
        self.workflow_runner = self.generate_workflow_runner(
            workflow_client=self.workflow_client,
            workflow_engine=self.workflow_engine)
        self.populate_workflows()

    def configure_request_client(self):
        orig_patch = self.client.patch
        def json_patch(path, data=None, **kwargs):
            return orig_patch(path, json.dumps(data),
                              content_type='application/json', **kwargs)
        self.client.patch = json_patch

    def generate_keyed_node_classes(self):
        class BaseTestNode(BaseNode):
            def tick(self, *args, **kwargs):
                print("tick", self.__class__.__name__)
                if not hasattr(self.__class__, 'ticks'):
                    self.__class__.ticks = 0
                self.__class__.ticks += 1

        node_classes = []
        class Node_0(BaseTestNode):
            def tick(self, *args, **kwargs):
                super().tick()
                next_node = Node_1(status='PENDING')
                self.workflow.add_node(next_node, precursors=[self])
                self.status = 'COMPLETED'
        node_classes.append(Node_0)

        class Node_1(BaseTestNode):
            allow_completion = False

            def tick(self, *args, **kwargs):
                super().tick()
                if self.__class__.allow_completion:
                    next_node = Node_2(status='PENDING')
                    self.workflow.add_node(next_node, precursors=[self])
                    self.status = 'COMPLETED'
        node_classes.append(Node_1)

        class Node_2(BaseTestNode):
            allow_completion = False

            def tick(self, *args, **kwargs):
                super().tick()
                if self.__class__.allow_completion:
                    self.status = 'COMPLETED'
        node_classes.append(Node_2)

        keyed_node_classes = {node_class.__name__: node_class
                              for node_class in node_classes}
        return keyed_node_classes

    def generate_workflow_engine(self, keyed_node_classes=None):
        workflow_engine = WorkflowEngine()
        for node_class in keyed_node_classes.values():
            workflow_engine.register_node_class(node_class=node_class)
        return workflow_engine

    def generate_workflow_runner(self, workflow_engine=None,
                                 workflow_client=None):
        workflow_runner = BaseWorkflowRunner(workflow_engine=workflow_engine,
                                             workflow_client=workflow_client)
        return workflow_runner

    def populate_workflows(self):
        workflow = Workflow()
        root_node = self.keyed_node_classes['Node_0'](
            workflow=workflow, status='PENDING')
        workflow.add_node(root_node)
        workflow.root_node = root_node
        serialization = self.workflow_engine.serialize_workflow(workflow)
        workflow_model = WorkflowModel.objects.create(
            serialization=json.dumps(serialization))
        self.workflow_uuid = workflow_model.uuid

    def test_run_workflows(self):
        states = []
        states.append(self.get_state(prev_state=None))

        def tick_and_capture_state(label=None):
            print("tick_and_capture", label)
            self.workflow_runner.tick()
            states.append(self.get_state(prev_state=states[-1]))

        tick_and_capture_state(label='a')
        self.assertEqual(
            states[-1]['node_classes']['Node_0']['ticks'], 1)

        tick_and_capture_state(label='b')
        self.assertEqual(
            states[-1]['node_classes']['Node_0']['ticks'], 1)
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 1)

        tick_and_capture_state(label='c')
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 2)
        self.assertFalse(
            'ticks' in states[-1]['node_classes']['Node_2'])

        self.keyed_node_classes['Node_1'].allow_completion = True
        tick_and_capture_state(label='d')
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 3)

        tick_and_capture_state(label='e')
        self.assertEqual(
            states[-1]['node_classes']['Node_1']['ticks'], 3)
        self.assertEqual(
            states[-1]['node_classes']['Node_2']['ticks'], 1)

        self.keyed_node_classes['Node_2'].allow_completion = True
        tick_and_capture_state(label='f')
        self.assertEqual(
            states[-1]['node_classes']['Node_2']['ticks'], 2)
        self.assertEqual(
            states[-1]['workflows'][self.workflow_uuid]['status'], 'COMPLETED')

        tick_and_capture_state(label='g')
        self.assertEqual(
            states[-2]['workflows'][self.workflow_uuid]['modified'],
            states[-1]['workflows'][self.workflow_uuid]['modified'])

    def get_state(self, prev_state=None):
        state = {
            'workflows': {
                workflow_model.uuid: workflow_model.__dict__
                for workflow_model in WorkflowModel.objects.all()
            },
            'node_classes': {
                key: node_class.__dict__
                for key, node_class in self.keyed_node_classes.items()
            }
        }
        return state
