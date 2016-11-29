from django.test import TestCase
from unittest.mock import call, patch, DEFAULT

from jobs.models import Job
from ..models import Workflow
from .. import utils


class BaseTestCase(TestCase):
    def setUp(self):
        self.runner = utils.StubPopulateNamesWorkflowRunner()

class TickTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.patchers = {
            'runner': patch.multiple(self.runner,
                                     generate_GenerateNames_job=DEFAULT,
                                     generate_IngestNames_job=DEFAULT,
                                     finalize=DEFAULT) 
        }
        self.mocks = {key: patcher.start()
                      for key, patcher in self.patchers.items()}

    def tearDown(self):
        for patcher in self.patchers.values(): patcher.stop()

    def test_GenerateNames(self):
        workflow = Workflow.objects.create()
        self.runner.tick(workflow=workflow)
        self.assertEqual(
            self.mocks['runner']['generate_GenerateNames_job'].call_args,
            call(workflow=workflow))

    def test_IngestNames(self):
        workflow = Workflow.objects.create()
        Job.objects.create(type='GenerateNames', workflow=workflow,
                           finished=True)
        self.runner.tick(workflow=workflow)
        self.assertEqual(
            self.mocks['runner']['generate_IngestNames_job'].call_args,
            call(workflow=workflow))

    def test_finalize(self):
        workflow = Workflow.objects.create()
        Job.objects.create(type='IngestNames', workflow=workflow,
                           finished=True)
        self.runner.tick(workflow=workflow)
        self.assertEqual(
            self.mocks['runner']['finalize'].call_args,
            call(workflow=workflow))

class GenerateGenerateNamesJobTestCase(BaseTestCase):
    def test_generates_GenerateNames_job(self, workflow=None):
        workflow = Workflow.objects.create()
        self.runner.generate_GenerateNames_job(workflow=workflow)
        job = Job.objects.first()
        self.assertEqual(job.type, 'GenerateNames')
        self.assertEqual(job.workflow, workflow)

class GenerateIngestNamesTestCase(BaseTestCase):
    def test_generates_IngestNames_job(self, workflow=None):
        workflow = Workflow.objects.create()
        self.runner.generate_IngestNames_job(workflow=workflow)
        job = Job.objects.first()
        self.assertEqual(job.type, 'IngestNames')
        self.assertEqual(job.workflow, workflow)

class FinalizeTestCase(BaseTestCase):
    def test_sets_status(self):
        workflow = Workflow.objects.create()
        self.runner.finalize(workflow=workflow)
        self.assertEqual(workflow.finished, True)
