from django.conf.urls import url, include
from django.test import TestCase, override_settings

from missions.models import Workflow as WorkflowModel
from workflow_runners.workflow import Workflow, BaseNode
from workflow_runners.base_workflow_runner import BaseWorkflowRunner


BASE_PATH = 'test_api'
urlpatterns = [
    url(r'^%s' % BASE_PATH, include('jobs.urls')),
]

class TestWorkflowRootNode(BaseNode):
    def tick(self):
        next_node = Node_1(workflow=self.workflow)
        self.workflow.add_node(next_node)
        self.workflow.connect_nodes(src=self, dest=next_node)
        self.status = 'COMPLETED'

class Node_1(BaseNode):
    def tick(self):
        next_node = Node_2(workflow=self.workflow)
        self.workflow.add_node(next_node)
        self.workflow.connect_nodes(src=self, dest=next_node)
        self.status = 'COMPLETED'

class Node_2(BaseNode):
    def tick(self):
        self.status = 'COMPLETED'

@override_settings(ROOT_URLCONF=__name__)
class WorkflowRunner_TestCase(TestCase):
    def setUp(self):
        self.workflow_runner = self.generate_workflow_runner()
        self.populate_workflows()

    def generate_workflow_runner(self):
        workflow_runner = BaseWorkflowRunner()
        for node_class in [TestWorkflowRootNode, Node_1, Node_2]:
            workflow_runner.register_node_class(node_class=node_class)
        return workflow_runner

    def populate_workflows(self):
        workflow = Workflow()
        workflow.root_node = TestWorkflowRootNode(workflow=workflow)
        serialization = self.workflow_runner.serialize_workflow(workflow)
        self.workflow_model = WorkflowModel.objects.create(
            serialization=serialization)
        self.claimed_workflow_model = WorkflowModel.objects.create(claimed=True)

    def test_run_workflows(self):
        states = []
        states.append(self.get_state(prev_state=None))
        self.workflow_runner.tick()
        states.append(self.get_state(prev_state=states[0]))
        self.fail()

    def get_state(self, prev_state=None):
        state = {
            'workflows': {workflow_model.id: workflow_model
                          for workflow_model in WorkflowModel.objects.all()},
        }
        return state
