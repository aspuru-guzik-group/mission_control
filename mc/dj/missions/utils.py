import collections

from jobs.models import JobStatuses
from .models import WorkflowStatuses, WorkflowJob, WorkflowSpec
from .spiff_workflow_engine import SpiffWorkflowEngine


def start_workflow(spec_key=None, mission=None):
    workflow_spec = WorkflowSpec.objects.get(key=spec_key)
    workflow = generate_workflow_from_spec(spec=workflow_spec)
    workflow.spec = workflow_spec
    workflow.status = WorkflowStatuses.Running.name
    workflow.mission = mission
    workflow.save()
    run_workflow_tick(workflow=workflow)
    return workflow

def generate_workflow_from_spec(spec=None):
    engine = get_workflow_engine()
    workflow = engine.generate_workflow_from_spec(spec=spec)
    return workflow

def get_workflow_engine():
    return SpiffWorkflowEngine()

def run_workflow_tick(workflow=None):
    engine = get_workflow_engine()
    engine.run_workflow_tick(workflow=workflow)

def poll_workflow_jobs():
    finished_job_statuses = [JobStatuses.Completed.name]
    newly_finished_workflow_jobs_qset = WorkflowJob.objects.filter(
        finished=False, job__status__in=finished_job_statuses)
    # Capture jobs in list because update operation will change the queryset.
    newly_finished_workflow_jobs = list(newly_finished_workflow_jobs_qset.all())
    newly_finished_workflow_jobs_qset.update(finished=True)
    process_finished_workflow_jobs(workflow_jobs=newly_finished_workflow_jobs)

def process_finished_workflow_jobs(workflow_jobs=None):
    workflow_jobs_by_workflow = group_workflow_jobs_by_workflow(workflow_jobs)
    for workflow, workflow_jobs_group in workflow_jobs_by_workflow.items():
        engine = get_workflow_engine()
        engine.process_finished_workflow_jobs(workflow=workflow,
                                              workflow_jobs=workflow_jobs_group)
        run_workflow_tick(workflow=workflow)

def group_workflow_jobs_by_workflow(workflow_jobs=None):
    workflow_jobs_by_workflow = collections.defaultdict(list)
    for workflow_job in workflow_jobs:
        workflow_jobs_by_workflow[workflow_job.workflow].append(workflow_job)
    return workflow_jobs_by_workflow

def create_spec(key=None, serialization=None, serialization_format='xml'):
    engine = get_workflow_engine()
    normalized_serialization = engine.normalize_spec_serialization(
        serialization=serialization, input_format=serialization_format)
    spec = WorkflowSpec.objects.create(key=key,
                                       serialization=normalized_serialization)
    return spec

