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
    workflow_jobs_to_mark_as_finished_qset = WorkflowJob.objects.filter(
        finished=False, job__status__in=finished_job_statuses)
    # Need to capture jobs, because update operation will change the queryset.
    workflow_jobs_to_mark_as_finished = list(
        workflow_jobs_to_mark_as_finished_qset.all())
    workflow_jobs_to_mark_as_finished_qset.update(finished=True)
    for workflow_job in workflow_jobs_to_mark_as_finished:
        run_workflow_tick(workflow=workflow_job.workflow)
