from jobs.models import Job

from .models import Workflow


def start_workflow(runner_key=None, mission=None):
    workflow = Workflow.objects.create(runner_key=runner_key, mission=mission)
    runner = get_workflow_runner(runner_key=workflow.runner_key)
    runner.tick(workflow=workflow)
    return workflow

def get_workflow_runner(runner_key=None, runner_loader=None):
    if not runner_loader:
        runner_loader = HardcodedWorkflowRunnerLoader()
    return runner_loader.load_runner(key=runner_key)

# @TODO: change this to be dynamic!
class HardcodedWorkflowRunnerLoader(object):
    def load_runner(self, key=None):
        return StubPopulateNamesWorkflowRunner()

class StubPopulateNamesWorkflowRunner(object):
    def tick(self, workflow=None):
        if not workflow.jobs.exists():
            self.generate_GenerateNames_job(workflow=workflow)
        else:
            if workflow.last_finished_job.type == 'GenerateNames':
                self.generate_IngestNames_job(workflow=workflow)
            elif workflow.last_finished_job.type == 'IngestNames':
                self.finalize(workflow=workflow)

    def generate_GenerateNames_job(self, workflow=None):
        Job.objects.create(workflow=workflow, type='GenerateNames')

    def generate_IngestNames_job(self, workflow=None):
        Job.objects.create(workflow=workflow, type='IngestNames')

    def finalize(self, workflow=None):
        workflow.finished = True
        workflow.save()
