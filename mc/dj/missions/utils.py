from jobs.models import Job, JobStatuses

from .models import Workflow, WorkflowStatuses, WorkflowJob, WorkflowSpec
from .spiff_workflow_engine import SpiffWorkflowEngine


def start_workflow(spec_key=None, mission=None):
    workflow_spec = WorkflowSpec.objects.get(key=spec_key)
    workflow = Workflow(spec=workflow_spec, mission=mission)
    workflow.serialization = serialize_workflow(workflow=workflow)
    workflow.status = WorkflowStatuses.Running.name
    workflow.save()
    run_workflow_tick(workflow=workflow)
    return workflow

def serialize_workflow(workflow=None):
    engine = get_workflow_engine()
    return engine.serialize_workflow(workflow=workflow)

def get_workflow_engine():
    return SpiffWorkflowEngine()

def run_workflow_tick(workflow=None):
    engine = get_workflow_engine()
    engine.run_workflow_tick(workflow=workflow)

def poll_workflow_jobs():
    finished_job_statuses = [JobStatuses.Completed.name]
    workflow_jobs_to_mark_as_finished_qset = WorkflowJob.objects.filter(
        finished=False, job__status__in=finished_job_statuses)
    # Need to capture jobs, because update operation will change the queryset.
    workflow_jobs_to_mark_as_finished = list(
        workflow_jobs_to_mark_as_finished_qset.all())
    workflow_jobs_to_mark_as_finished_qset.update(finished=True)
    for workflow_job in workflow_jobs_to_mark_as_finished:
        run_workflow_tick(workflow=workflow_job.workflow)
