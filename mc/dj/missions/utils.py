import glob

from jobs.models import Job, JobStatuses
from .models import Workflow, WorkflowJob, WorkflowRunner


def start_workflow(runner_key=None, mission=None):
    workflow = Workflow.objects.create(runner_key=runner_key, mission=mission)
    run_workflow_tick(workflow=workflow)
    return workflow

def run_workflow_tick(workflow=None):
    runner = get_workflow_runner(runner_key=workflow.runner_key)
    runner.run_workflow_tick(workflow=workflow)

def get_workflow_runner(runner_key=None, runner_loader=None):
    if not runner_loader:
        runner_loader = HardcodedWorkflowRunnerLoader()
    return runner_loader.load_runner(key=runner_key)

def poll_workflow_jobs():
    finished_job_statuses = [JobStatuses.Completed.name]
    workflow_jobs_to_mark_as_finished_qset = WorkflowJob.objects.filter(
        finished=False, job__status__in=finished_job_statuses)
    # Need to capture jobs, because update operation will change the queryset.
    workflow_jobs_to_mark_as_finished = list(
        workflow_jobs_to_mark_as_finished_qset.all())
    workflow_jobs_to_mark_as_finished_qset.update(finished=True)
    run_workflow_job_workflow_ticks(
        workflow_jobs=workflow_jobs_to_mark_as_finished)

def run_workflow_job_workflow_ticks(workflow_jobs=None):
    for workflow_job in workflow_jobs:
        run_workflow_tick(workflow=workflow_job.workflow)

# @TODO: change this to be dynamic!
class HardcodedWorkflowRunnerLoader(object):
    def load_runner(self, key=None):
        return StubPopulateNamesWorkflowRunner()

class StubPopulateNamesWorkflowRunner(object):
    def run_workflow_tick(self, workflow=None):
        if not workflow.jobs.exists():
            self.generate_GenerateNames_job(workflow=workflow)
        else:
            last_finished_job = workflow.last_finished_job
            if not last_finished_job: return
            if workflow.last_finished_job.type == 'GenerateNames':
                self.generate_IngestNames_job(workflow=workflow)
            elif workflow.last_finished_job.type == 'IngestNames':
                self.finalize(workflow=workflow)

    def generate_GenerateNames_job(self, workflow=None):
        WorkflowJob.objects.create(
            workflow=workflow,
            job=Job.objects.create(type='GenerateNames')
        )

    def generate_IngestNames_job(self, workflow=None):
        WorkflowJob.objects.create(
            workflow=workflow,
            job=Job.objects.create(type='IngestNames')
        )

    def finalize(self, workflow=None):
        workflow.finished = True
        workflow.save()

def sync_workflow_runners(runner_dirs=None):
    result = {'created': {}, 'synced': {}, 'missing': {}, 'errors': {}}
    if not runner_dirs:
        runner_dirs = get_default_runner_dirs()
    paths_in_dirs = set(find_runner_files(dirs=runner_dirs))
    paths_in_db = set([runner.path for runner in WorkflowRunner.objects.all()])

    paths_in_dirs_only = paths_in_dirs.difference(paths_in_db)
    for path in paths_in_dirs_only:
        try:
            result['created'][path] = add_workflow_runner(runner_path=path)
        except Exception as error:
            result['errors'][path] = error

    paths_in_db_only = paths_in_db.difference(paths_in_dirs)
    for path in paths_in_db_only:
        error = MissingRunnerFileError(path)
        WorkflowRunner.objects.filter(path=path).update(error=error)
        result['errors'][path] = error

    paths_in_dirs_and_db = paths_in_dirs.intersection(paths_in_db)
    for path in paths_in_dirs_and_db:
        try: result['synced'][path] = sync_workflow_runner(runner_path=path)
        except Exception as error:
            WorkflowRunner.objects.filter(path=path).update(error=error)
            result['errors'][path] = error

    return result

class MissingRunnerFileError(Exception): pass

def find_runner_files(dirs=None):
    runner_file_paths = []
    for _dir in dirs:
        runner_file_paths.extend(list(glob.glob('%s/**/*_wf_runner.py' % _dir,
                                                recursive=True)))
    return runner_file_paths

def get_default_runner_dirs():
    raise NotImplementedError()

def add_workflow_runner(runner_path=None):
    raise NotImplementedError()

def sync_workflow_runner(runner_path=None):
    raise NotImplementedError()


