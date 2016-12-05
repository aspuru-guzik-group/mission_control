import SpiffWorkflow
import SpiffWorkflow.specs
from SpiffWorkflow.storage import JSONSerializer, XmlSerializer

from jobs.models import Job
from .models import Workflow, WorkflowJob

TASK_RUNNING_KEY = '_running'


class SpiffWorkflowEngine(object):
    def generate_workflow_from_spec(self, spec=None):
        spiff_spec = self.deserialize_spiff_spec(
            serialization=spec.serialization)
        spiff_workflow = SpiffWorkflow.Workflow(spiff_spec)
        workflow = Workflow(
            serialization=self.serialize_spiff_workflow(
                spiff_workflow=spiff_workflow)
        )
        return workflow

    def deserialize_spiff_spec(self, serialization=None, format=None):
        return SpiffWorkflow.specs.WorkflowSpec.deserialize(
            serialization, self.get_serializer(format=format))

    def serialize_spiff_spec(self, spiff_spec=None, format=None):
        return spiff_spec.serialize(self.get_serializer(format=format))

    def get_serializer(self, format='json'):
        if not format: format = 'json'
        if format == 'json':
            return JSONSerializer()
        elif format == 'xml':
            return XmlSerializer()

    def serialize_spiff_workflow(self, spiff_workflow=None):
        return spiff_workflow.serialize(self.get_serializer())

    def run_workflow_tick(self, workflow=None):
        spiff_workflow = self.deserialize_spiff_workflow(
            serialization=workflow.serialization)
        for spiff_task in self.get_pending_tasks(spiff_workflow=spiff_workflow):
            self.start_spiff_task(spiff_task=spiff_task, workflow=workflow)
        workflow.serialization = self.serialize_spiff_workflow(
            spiff_workflow=spiff_workflow)
        workflow.save()

    def get_pending_tasks(self, spiff_workflow=None):
        ready_tasks = spiff_workflow.get_tasks(SpiffWorkflow.Task.READY)
        pending_tasks = [task for task in ready_tasks
                         if not task.get_data(TASK_RUNNING_KEY)]
        return pending_tasks

    def deserialize_spiff_workflow(self, serialization=None):
        return SpiffWorkflow.Workflow.deserialize(
            serialization, self.get_serializer())

    def start_spiff_task(self, spiff_task=None, workflow=None):
        WorkflowJob.objects.create(
            job=Job.objects.create(spec=spiff_task.get_data('job_spec')),
            workflow=workflow)
        spiff_task.set_data(**{TASK_RUNNING_KEY: True})


    def normalize_spec_serialization(self, serialization=None,
                                     input_format=None):
        spiff_spec = self.deserialize_spiff_spec(serialization=serialization,
                                                 format=input_format)
        normalized_serialization = self.serialize_spiff_spec(
            spiff_spec=spiff_spec)
        return normalized_serialization
