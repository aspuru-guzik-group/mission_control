import collections
import json
import unittest
from unittest.mock import DEFAULT, MagicMock

from django.test import override_settings

from missions.models import Flow as FlowModel, Job as JobModel
from ..odyssey_push_runner import OdysseyPushRunner
from . import e2e_utils
from mc.flow_engines.flow import Flow
from mc.a2g2.node_engines.job_node_engine import JobNodeEngine
from mc.a2g2.node_engines.flow_node_engine import FlowNodeEngine


urlpatterns = e2e_utils.urlpatterns
assert urlpatterns

@override_settings(ROOT_URLCONF=__name__)
class OdysseyPushRunnerE2ETestCase(e2e_utils.BaseTestCase):
    def setUp(self):
        super().setUp()
        self.execution_client = self.generate_execution_client()
        self.task_runner = self.generate_task_runner()
        self.job_submission_factory = MagicMock()
        self.flow_generator_classes = self.generate_flow_generator_classes()
        self.runner = self.generate_runner()
        self.flow_spec = self.generate_flow_spec()
        self.parent_flow_model = self.generate_flow_model()

    def generate_execution_client(self):
        execution_client = MagicMock()
        execution_client.start_execution.return_value = \
                collections.defaultdict(MagicMock)
        return execution_client

    def generate_task_runner(self):
        task_runner = MagicMock()
        return task_runner

    def generate_flow_generator_classes(self):
        flow_generator_classes = {}
        class FlowGeneratorA:
            flow_type = 'FlowA'

            @classmethod
            def generate_flow(cls, *args, flow_spec=None, **kwargs):
                flow = Flow()
                flow.add_node(
                    key='a.1',
                    node={
                        'node_engine': JobNodeEngine.__name__,
                        'input': {
                            'job_spec': {
                                'type': 'Job_A',
                            }
                        },
                        'status': 'PENDING',
                    },
                    as_root=True,
                )
                return flow

            @classmethod
            def get_dependencies(cls):
                return {
                    'node_engines': set([JobNodeEngine()])
                }
        flow_generator_classes['A'] = FlowGeneratorA

        class FlowGeneratorB:
            flow_type = 'FlowB'

            @classmethod
            def generate_flow(cls, *args, flow_spec=None, **kwargs):
                flow = Flow()
                flow.add_node(
                    key='b.1',
                    as_root=True,
                    node={
                        'node_engine': FlowNodeEngine.__name__,
                        'input': {
                            'flow_spec': {
                                'flow_type': FlowGeneratorA.flow_type
                            }
                        },
                        'status': 'PENDING'
                    }
                )
                return flow

            @classmethod
            def get_dependencies(cls):
                return {
                    'node_engines': set([FlowNodeEngine()]),
                    'flow_generator_classes': set([FlowGeneratorA])
                }
        flow_generator_classes['B'] = FlowGeneratorB
        return flow_generator_classes

    def generate_runner(self):
        mc_server_url = '/%s/' % e2e_utils.BASE_PATH
        runner = OdysseyPushRunner(
            task_runner=self.task_runner,
            request_client=self.client,
            run_setup=True,
            mc_server_url=mc_server_url,
            flow_generator_classes=self.flow_generator_classes.values(),
            job_submission_factory=self.job_submission_factory,
            job_runner_kwargs={'execution_client': self.execution_client}
        )
        return runner

    def generate_flow_spec(self):
        flow_spec = {'flow_type': self.flow_generator_classes['B'].flow_type}
        return flow_spec

    def generate_flow_model(self):
        flow_model = FlowModel.objects.create(
            spec=json.dumps(self.flow_spec)
        )
        return flow_model

    def test_runs_flows(self):
        self.runner.tick()
        self.assert_parent_flow_model_attr('status', 'RUNNING')
        self.child_flow_model = self.get_child_flow_model()
        self.assertTrue(self.child_flow_model is not None)

        self.runner.tick()
        self.assert_parent_flow_model_attr('status', 'RUNNING')
        self.job_model = self.get_job_model()
        self.assertTrue(self.job_model is not None)

        self.runner.tick()
        self.assert_job_model_attr('status', 'RUNNING')
        self.assert_child_flow_model_attr('status', 'RUNNING')
        self.assert_parent_flow_model_attr('status', 'RUNNING')

        self.mock_job_execution()
        self.runner.tick()

        self.runner.tick()

        self.assert_job_model_attr('status', 'COMPLETED')

        self.runner.tick()
        self.assert_child_flow_model_attr('status', 'COMPLETED')
        self.runner.tick()
        self.assert_parent_flow_model_attr('status', 'COMPLETED')

    def get_parent_flow_model(self):
        return FlowModel.objects.get(uuid=self.parent_flow_model.uuid)

    def assert_parent_flow_model_attr(self, attr, expected):
        parent_flow_model = self.get_parent_flow_model()
        self.assertEqual(getattr(parent_flow_model, attr), expected)

    def get_child_flow_model(self):
        return FlowModel.objects.exclude(uuid=self.parent_flow_model.uuid)[0]

    def assert_child_flow_model_attr(self, attr, expected):
        child_flow_model = self.get_child_flow_model()
        self.assertEqual(getattr(child_flow_model, attr), expected)

    def get_job_model(self):
        job_models = JobModel.objects.all()
        job_model = [job_model for job_model in job_models
                     if job_model.job_spec['type'] == 'Job_A'][0]
        return job_model

    def assert_job_model_attr(self, attr, expected):
        job_model = self.get_job_model()
        self.assertEqual(getattr(job_model, attr), expected)

    def mock_job_execution(self):
        def mock_tick_task(task=None, job=None, **kwargs):
            if task['type'] == 'job:execute':
                task['status'] = 'COMPLETED'
                job['execution'] = {'result': {'result': 'COMPLETED'}}
            return DEFAULT
        self.runner.task_runner.tick_task.side_effect = mock_tick_task

if __name__ == '__main__':
    unittest.main()
