from unittest.mock import call, DEFAULT, MagicMock, patch
from django.test import TestCase

from ..models import WorkflowSpec, Workflow, WorkflowJob
from ..spiff_workflow_engine import SpiffWorkflowEngine


class BaseTestCase(TestCase):
    def setUp(self):
        self.engine = SpiffWorkflowEngine()

    def start_patchers(self, patchers):
        mocks = {key: patcher.start() for key, patcher in patchers.items()}
        return mocks

    def tearDown(self):
        if hasattr(self, 'patchers'):
            for patcher in self.patchers.values(): patcher.stop()

class GenerateWorkflowFromSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine,
                                     deserialize_spiff_spec=DEFAULT,
                                     serialize_spiff_workflow=DEFAULT),
            'SpiffWorkflow': patch.multiple('SpiffWorkflow', Workflow=DEFAULT),
        }
        self.mocks = self.start_patchers(self.patchers)
        self.spec = WorkflowSpec(serialization='some serialization')

    def test_generates_workflow(self):
        workflow_model = self.engine.generate_workflow_from_spec(spec=self.spec)
        mock_deserialize_spiff_spec = self.mocks['engine'].get(
            'deserialize_spiff_spec')
        self.assertEqual(
            mock_deserialize_spiff_spec.call_args,
            call(serialization=self.spec.serialization))
        mock_spiff_Workflow = self.mocks['SpiffWorkflow']['Workflow']
        expected_spiff_spec = mock_deserialize_spiff_spec.return_value
        self.assertEqual(mock_spiff_Workflow.call_args,
                         call(expected_spiff_spec))
        expected_spiff_workflow = mock_spiff_Workflow.return_value
        mock_serialize_spiff_workflow = self.mocks['engine'].get(
            'serialize_spiff_workflow')
        self.assertEqual(mock_serialize_spiff_workflow.call_args,
                         call(spiff_workflow=expected_spiff_workflow))
        expected_serialization = mock_serialize_spiff_workflow.return_value
        expected_workflow_attrs = {
            'serialization': expected_serialization,
        }
        workflow_attrs = {attr: getattr(workflow_model, attr)
                          for attr in expected_workflow_attrs}
        self.assertEqual(workflow_attrs, expected_workflow_attrs)

class DeserializeSpiffSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine, get_serializer=DEFAULT),
            'SpiffWorkflow.specs': patch.multiple('SpiffWorkflow.specs',
                                                 WorkflowSpec=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mock_serialization = MagicMock()

    def test_deserializes_spiff_spec(self):
        spiff_spec = self.engine.deserialize_spiff_spec(
            serialization=self.mock_serialization)
        mock_spiff_WorkflowSpec = self.mocks['SpiffWorkflow.specs'].get(
            'WorkflowSpec')
        mock_get_serializer = self.mocks['engine']['get_serializer']
        self.assertEqual(
            mock_spiff_WorkflowSpec.deserialize.call_args,
            call(self.mock_serialization, mock_get_serializer.return_value))
        self.assertEqual(spiff_spec,
                         mock_spiff_WorkflowSpec.deserialize.return_value)

class SerializeSpiffWorkflowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine, get_serializer=DEFAULT),
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mock_spiff_workflow = MagicMock()

    def test_serializes_spiff_workflow(self):
        serialization = self.engine.serialize_spiff_workflow(
            spiff_workflow=self.mock_spiff_workflow)
        self.assertEqual(
            self.mock_spiff_workflow.serialize.call_args,
            call(self.mocks['engine']['get_serializer'].return_value))
        self.assertEqual(serialization,
                         self.mock_spiff_workflow.serialize.return_value)

class RunWorkflowTickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine, get_serializer=DEFAULT,
                                     serialize_spiff_workflow=DEFAULT,
                                     deserialize_spiff_workflow=DEFAULT,
                                     start_spiff_task=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mock_workflow = MagicMock()
        self.task_groups = self.generate_task_groups()

    def generate_task_groups(self):
        task_groups = {'all': [], 'running': [], 'pending': []}
        for i in range(6):
            task = MagicMock()
            running = ((i % 2) == 0)
            task.get_data.return_value = running
            task_groups['all'].append(task)
            if running:
                task_groups['running'].append(task)
            else:
                task_groups['pending'].append(task)
        return task_groups

    def test_deserializes_workflow(self):
        initial_serialization = self.mock_workflow.serialization
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        self.assertEqual(
            self.mocks['engine']['deserialize_spiff_workflow'].call_args,
            call(serialization=initial_serialization))

    def test_starts_pending_tasks(self):
        mock_deserialize = self.mocks['engine']['deserialize_spiff_workflow']
        deserialized_workflow = mock_deserialize.return_value
        deserialized_workflow.get_tasks.return_value = self.task_groups['all']
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        call_args_list = self.mocks['engine']['start_spiff_task'].call_args_list
        expected_call_args_list = [
            call(spiff_task=task, workflow=self.mock_workflow)
            for task in self.task_groups['pending']]
        self.assertEqual(call_args_list, expected_call_args_list)

    def test_serializes_workflow(self):
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        mock_deserialize = self.mocks['engine']['deserialize_spiff_workflow']
        expected_spiff_workflow = mock_deserialize.return_value
        self.assertEqual(
            self.mocks['engine']['serialize_spiff_workflow'].call_args,
            call(spiff_workflow=expected_spiff_workflow))

    def test_saves_workflow(self):
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        self.assertEqual(self.mock_workflow.save.call_count, 1)

class DeserializeSpiffWorkflow(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine, get_serializer=DEFAULT),
            'SpiffWorkflow': patch.multiple('SpiffWorkflow', Workflow=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mock_serialization = MagicMock()

    def test_deserializes_spiff_workflow(self):
        spiff_workflow = self.engine.deserialize_spiff_workflow(
            serialization=self.mock_serialization)
        mock_spiff_Workflow = self.mocks['SpiffWorkflow']['Workflow']
        mock_get_serializer = self.mocks['engine']['get_serializer']
        self.assertEqual(
            mock_spiff_Workflow.deserialize.call_args,
            call(self.mock_serialization, mock_get_serializer.return_value))
        self.assertEqual(
            spiff_workflow,
            mock_spiff_Workflow.deserialize.return_value)

class StartSpiffTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_spec = {'junk': 'data'}
        self.spiff_task = MagicMock()
        self.spiff_task.get_data.return_value = self.job_spec
        self.workflow = Workflow.objects.create()

    def test_creates_job_for_task(self):
        self.engine.start_spiff_task(spiff_task=self.spiff_task,
                                     workflow=self.workflow)
        self.assertEqual(WorkflowJob.objects.count(), 1)
        workflow_job = WorkflowJob.objects.all()[0]
        self.assertEqual(workflow_job.workflow, self.workflow)
        self.assertEqual(workflow_job.job.spec,
                         self.spiff_task.get_data('job_spec'))

    def test_marks_task_as_running(self):
        self.engine.start_spiff_task(spiff_task=self.spiff_task,
                                     workflow=self.workflow)
        self.assertEqual(self.spiff_task.set_data.call_args,
                         call(_running=True))
