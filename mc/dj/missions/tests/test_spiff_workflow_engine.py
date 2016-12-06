import uuid
from unittest.mock import call, DEFAULT, MagicMock, patch

from django.test import TestCase

from ..models import WorkflowSpec, Workflow, WorkflowJob, WorkflowStatuses
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
            'engine': patch.multiple(
                self.engine,
                deserialize_spiff_spec=DEFAULT,
                serialize_spiff_workflow=DEFAULT,
                complete_spiff_workflow_start_task=DEFAULT
            ),
            'SpiffWorkflow': patch.multiple('SpiffWorkflow', Workflow=DEFAULT),
        }
        self.mocks = self.start_patchers(self.patchers)
        self.spec = WorkflowSpec(serialization='some serialization')
        self.expected_spiff_workflow = self.mocks['SpiffWorkflow']['Workflow']\
                .return_value

    def test_creates_spiff_workflow_from_spec(self):
        self.engine.generate_workflow_from_spec(spec=self.spec)
        mock_deserialize_spiff_spec = self.mocks['engine'].get(
            'deserialize_spiff_spec')
        self.assertEqual(
            mock_deserialize_spiff_spec.call_args,
            call(serialization=self.spec.serialization))
        mock_spiff_Workflow = self.mocks['SpiffWorkflow']['Workflow']
        expected_spiff_spec = mock_deserialize_spiff_spec.return_value
        self.assertEqual(mock_spiff_Workflow.call_args,
                         call(expected_spiff_spec))

    def test_completes_start_task_on_spiff_workflow(self):
        self.engine.generate_workflow_from_spec(spec=self.spec)
        mock_complete_start_task_fn = self.mocks['engine'].get(
            'complete_spiff_workflow_start_task')
        self.assertEqual(mock_complete_start_task_fn.call_args,
                         call(spiff_workflow=self.expected_spiff_workflow))

    def test_created_workflow_has_expected_serialization(self):
        workflow_model = self.engine.generate_workflow_from_spec(spec=self.spec)
        mock_serialize_spiff_workflow = self.mocks['engine'].get(
            'serialize_spiff_workflow')
        self.assertEqual(mock_serialize_spiff_workflow.call_args,
                         call(spiff_workflow=self.expected_spiff_workflow))
        expected_serialization = mock_serialize_spiff_workflow.return_value
        self.assertEqual(workflow_model.serialization, expected_serialization)

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
            call(mock_get_serializer.return_value, self.mock_serialization))
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
        mock_deserialize = self.mocks['engine']['deserialize_spiff_workflow']
        self.deserialized_workflow = mock_deserialize.return_value
        self.deserialized_workflow.is_completed.return_value = False
        self.task_groups = self.generate_task_groups()

    def generate_task_groups(self):
        task_groups = {'all': [], 'running': [], 'pending': []}
        for i in range(6):
            task = MagicMock()
            running = ((i % 2) == 0)
            task._get_internal_data.return_value = running
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
        self.deserialized_workflow.get_tasks.return_value = \
                self.task_groups['all']
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        call_args_list = self.mocks['engine']['start_spiff_task'].call_args_list
        expected_call_args_list = [
            call(spiff_task=task, workflow=self.mock_workflow)
            for task in self.task_groups['pending']]
        self.assertEqual(call_args_list, expected_call_args_list)

    def test_serializes_workflow(self):
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        self.assertEqual(
            self.mocks['engine']['serialize_spiff_workflow'].call_args,
            call(spiff_workflow=self.deserialized_workflow))

    def test_marks_completed_workflow_as_completed(self):
        self.deserialized_workflow.is_completed.return_value = True
        self.engine.run_workflow_tick(workflow=self.mock_workflow)
        self.assertEqual(self.mock_workflow.status,
                         WorkflowStatuses.Completed.name)

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
            call(mock_get_serializer.return_value, self.mock_serialization))
        self.assertEqual(
            spiff_workflow,
            mock_spiff_Workflow.deserialize.return_value)

class StartSpiffTaskTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job_spec = {'junk': 'data'}
        self.spiff_task = MagicMock()
        self.spiff_task.get_spec_data.return_value = self.job_spec
        self.workflow = Workflow.objects.create()

    def test_creates_job_for_task(self):
        self.engine.start_spiff_task(spiff_task=self.spiff_task,
                                     workflow=self.workflow)
        self.assertEqual(WorkflowJob.objects.count(), 1)
        workflow_job = WorkflowJob.objects.all()[0]
        self.assertEqual(workflow_job.workflow, self.workflow)
        self.assertEqual(workflow_job.job.spec,
                         self.spiff_task.get_spec_data('job_spec'))
        self.assertEqual(workflow_job.meta,
                         {'task_uuid': str(self.spiff_task.id)})

    def test_marks_task_as_running(self):
        self.engine.start_spiff_task(spiff_task=self.spiff_task,
                                     workflow=self.workflow)
        self.assertEqual(self.spiff_task._set_internal_data.call_args,
                         call(_running=True))

class NormalizeSpecSerializationTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine, serialize_spiff_spec=DEFAULT,
                                     deserialize_spiff_spec=DEFAULT),
        }
        self.mocks = self.start_patchers(self.patchers)
        self.serialization = 'some serialization'

    def test_returns_serialization_of_deserialization(self):
        normalized_serialization = self.engine.normalize_spec_serialization(
            serialization=self.serialization)
        mock_deserialize_fn = self.mocks['engine']['deserialize_spiff_spec']
        mock_serialize_fn = self.mocks['engine']['serialize_spiff_spec']
        self.assertEqual(mock_serialize_fn.call_args,
                         call(spiff_spec=mock_deserialize_fn.return_value))
        expected_normalized_serialization = mock_serialize_fn.return_value
        self.assertEqual(normalized_serialization,
                         expected_normalized_serialization)

class SerializeSpiffSpecTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine, get_serializer=DEFAULT),
        }
        self.mocks = self.start_patchers(self.patchers)
        self.mock_spiff_spec = MagicMock()

    def test_serializes_spiff_spec(self):
        serialization = self.engine.serialize_spiff_spec(
            spiff_spec=self.mock_spiff_spec)
        self.assertEqual(
            self.mock_spiff_spec.serialize.call_args,
            call(self.mocks['engine']['get_serializer'].return_value))
        self.assertEqual(serialization,
                         self.mock_spiff_spec.serialize.return_value)

class ProcessFinishedWorkflowJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'engine': patch.multiple(self.engine,
                                     deserialize_spiff_workflow=DEFAULT,
                                     serialize_spiff_workflow=DEFAULT)
        }
        self.mocks = self.start_patchers(self.patchers)
        self.workflow = MagicMock()
        self.workflow_jobs = self.generate_workflow_jobs(workflow=self.workflow)
        self.spiff_workflow = self.mocks['engine'].get(
            'deserialize_spiff_workflow').return_value

    def generate_workflow_jobs(self, workflow):
        workflow_jobs = []
        for i in range(3):
            workflow_job = MagicMock()
            workflow_job.workflow = workflow
            workflow_job.meta = {'task_uuid': str(uuid.uuid4())}
            workflow_jobs.append(workflow_job)
        return workflow_jobs

    def test_deserializes_workflow(self):
        initial_serialization = self.workflow.serialization
        self.engine.process_finished_workflow_jobs(
            workflow=self.workflow, workflow_jobs=self.workflow_jobs)
        self.assertEqual(
            self.mocks['engine']['deserialize_spiff_workflow'].call_args,
            call(serialization=initial_serialization))

    def test_completes_corresponding_tasks(self):
        self.engine.process_finished_workflow_jobs(
            workflow=self.workflow, workflow_jobs=self.workflow_jobs)
        call_args_list = self.spiff_workflow.complete_task_from_id\
                .call_args_list
        expected_call_args_list = [
            call(uuid.UUID(workflow_job.meta['task_uuid']))
            for workflow_job in self.workflow_jobs]
        self.assertEqual(call_args_list, expected_call_args_list)

    def test_reserializes_workflow(self):
        self.engine.process_finished_workflow_jobs(
            workflow=self.workflow, workflow_jobs=self.workflow_jobs)
        self.assertEqual(
            self.mocks['engine']['serialize_spiff_workflow'].call_args,
            call(spiff_workflow=self.spiff_workflow))

    def test_saves_worflow(self):
        self.engine.process_finished_workflow_jobs(
            workflow=self.workflow, workflow_jobs=self.workflow_jobs)
        self.assertEqual(self.workflow.save.call_count, 1)
