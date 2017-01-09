from flow_engines.flow import Flow

from .base import BaseFlowGenerator
from ..tasks.base import BaseTask
from ..tasks.job import JobTask


class RunAndLoadFlowGenerator(BaseFlowGenerator):
    flow_type = 'run_and_load'
    label = 'Run-and-Load Flow'
    description = """Flow that runs 2 jobs:
        (1) a generic job, and
        (2) a 'load' job"""

    @classmethod
    def generate_flow(cls, *arg, flow_spec=None, **kwargs):
        flow_kwargs = {'status': 'PENDING', **flow_spec.get('flow_kwargs', {})}
        flow = Flow(**flow_kwargs) 
        flow.data['flow_spec'] = flow_spec
        flow.add_task(
            key='run',
            as_root=True,
            task=JobTask(),
            input={'job_spec': flow_spec['run_spec']}
        )
        flow.add_task(
            key='load_prep',
            precursor='run',
            task=LoadPrepTask(flow_spec=flow_spec)
        )
        flow.add_task(
            key='load',
            precursor='load_prep',
            task=JobTask()
        )
        return flow

    @classmethod
    def get_dependencies(cls, *args, **kwargs):
        return {
            'task_classes': set([JobTask, LoadPrepTask])
        }

class LoadPrepTask(BaseTask):
    def tick(self, *args, **kwargs):
        load_spec = self.flow.data['flow_spec'].get('load_spec', {})
        self.output = {'job_spec': {**load_spec, 'dir': self.input['dir']}}
        self.status = 'COMPLETED'

