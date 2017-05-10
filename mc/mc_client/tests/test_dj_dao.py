import unittest
from unittest.mock import call, MagicMock

from .. import dj_dao


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.models = MagicMock()
        self.serializers = MagicMock()
        self.queue_utils = MagicMock()
        self.dao = dj_dao.DjDao(models=self.models,
                                serializers=self.serializers,
                                queue_utils=self.queue_utils)

    def generate_mocks(self, n=3):
        return [MagicMock() for i in range(n)]

class CreateFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.serialize_flow_model = MagicMock()
        self.flow = MagicMock()

    def _create_flow(self):
        return self.dao.create_flow(flow=self.flow)

    def test_dispatches_to_flow_model_create(self):
        self._create_flow()
        self.assertEqual(self.models['Flow'].objects.create.call_args,
                         call(**self.flow))

    def test_returns_serialized_result(self):
        result = self._create_flow()
        self.assertEqual(
            self.dao.serialize_flow_model.call_args,
            call(flow_model=self.models['Flow'].objects.create.return_value))
        self.assertEqual(result, self.dao.serialize_flow_model.return_value)

class GetFlowsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.flow_models = self.generate_mocks()
        self.dao.get_flow_models = MagicMock(return_value=self.flow_models)
        self.dao.serialize_flow_model = MagicMock()

    def _get_flows(self):
        return self.dao.get_flows(query=self.query)

    def test_gets_flow_models(self):
        self._get_flows()
        self.assertEqual(self.dao.get_flow_models.call_args,
                         call(query=self.query))

    def test_returns_serialized_flow_models(self):
        result = self._get_flows()
        self.assertEqual(result,
                         [self.dao.serialize_flow_model.return_value
                          for flow_model in self.flow_models])

class GetFlowModelsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.flow_query_to_queryset = MagicMock()

    def _get_flow_models(self):
        return self.dao.get_flow_models(query=self.query)

    def test_generates_queryset(self):
        self._get_flow_models()
        self.assertEqual(self.dao.flow_query_to_queryset.call_args,
                         call(query=self.query))

    def test_returns_queryset_results(self):
        result = self._get_flow_models()
        self.assertEqual(
            result,
            self.dao.flow_query_to_queryset.return_value.all.return_value)

class PatchFlowTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.key = MagicMock()
        self.patches = MagicMock()
        self.dao.serialize_flow_model = MagicMock()

    def _patch_flow(self):
        return self.dao.patch_flow(key=self.key, patches=self.patches)

    def test_dispatches_to_flow_model_update(self):
        self._patch_flow()
        self.assertEqual(self.models['Flow'].objects.get.call_args,
                         call(uuid=self.key))
        self.assertEqual(
            self.models['Flow'].objects.get.return_value.update.call_args,
            call(**self.patches))

    def test_returns_serialized_result(self):
        result = self._patch_flow()
        self.assertEqual(
            self.dao.serialize_flow_model.call_args,
            call(flow_model=self.models['Flow'].objects.get.return_value)
        )
        self.assertEqual(result, self.dao.serialize_flow_model.return_value)

class CreateJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.serialize_job_model = MagicMock()
        self.job = MagicMock()

    def _create_job(self):
        return self.dao.create_job(job=self.job)

    def test_dispatches_to_job_model_create(self):
        self._create_job()
        self.assertEqual(self.models['Job'].objects.create.call_args,
                         call(**self.job))

    def test_returns_serialized_result(self):
        result = self._create_job()
        self.assertEqual(
            self.dao.serialize_job_model.call_args,
            call(job_model=self.models['Job'].objects.create.return_value))
        self.assertEqual(result, self.dao.serialize_job_model.return_value)

class GetJobsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.job_models = self.generate_mocks()
        self.dao.get_job_models = MagicMock(return_value=self.job_models)
        self.dao.serialize_job_model = MagicMock()

    def _get_jobs(self):
        return self.dao.get_jobs(query=self.query)

    def test_gets_job_models(self):
        self._get_jobs()
        self.assertEqual(self.dao.get_job_models.call_args,
                         call(query=self.query))

    def test_returns_serialized_job_models(self):
        result = self._get_jobs()
        self.assertEqual(result,
                         [self.dao.serialize_job_model.return_value
                          for job_model in self.job_models])

class GetJobModelsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.query = MagicMock()
        self.dao.job_query_to_queryset = MagicMock()

    def _get_job_models(self):
        return self.dao.get_job_models(query=self.query)

    def test_generates_queryset(self):
        self._get_job_models()
        self.assertEqual(self.dao.job_query_to_queryset.call_args,
                         call(query=self.query))

    def test_returns_queryset_results(self):
        result = self._get_job_models()
        self.assertEqual(
            result,
            self.dao.job_query_to_queryset.return_value.all.return_value)

class PatchJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.key = MagicMock()
        self.patches = MagicMock()
        self.dao.serialize_job_model = MagicMock()

    def _patch_job(self):
        return self.dao.patch_job(key=self.key, patches=self.patches)

    def test_dispatches_to_job_model_update(self):
        self._patch_job()
        self.assertEqual(self.models['Job'].objects.get.call_args,
                         call(uuid=self.key))
        self.assertEqual(
            self.models['Job'].objects.get.return_value.update.call_args,
            call(**self.patches))

    def test_returns_serialized_result(self):
        result = self._patch_job()
        self.assertEqual(
            self.dao.serialize_job_model.call_args,
            call(job_model=self.models['Job'].objects.get.return_value)
        )
        self.assertEqual(result, self.dao.serialize_job_model.return_value)

class FlushTestCase(BaseTestCase):
    def test_dispatches_to_model_deletes(self):
        self.dao.flush_mc_db()
        self.assertEqual(self.models['Flow'].objects.delete.call_args, call())
        self.assertEqual(self.models['Job'].objects.delete.call_args, call())
        self.assertEqual(self.models['Queue'].objects.delete.call_args, call())

class CreateQueueTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.dao.serialize_queue_model = MagicMock()
        self.queue = MagicMock()

    def _create_queue(self):
        return self.dao.create_queue(queue=self.queue)

    def test_dispatches_to_queue_model_create(self):
        self._create_queue()
        self.assertEqual(self.models['Queue'].objects.create.call_args,
                         call(**self.queue))

    def test_returns_serialized_result(self):
        result = self._create_queue()
        self.assertEqual(
            self.dao.serialize_queue_model.call_args,
            call(queue_model=self.models['Queue'].objects.create.return_value))
        self.assertEqual(result, self.dao.serialize_queue_model.return_value)

class ClaimQueueItemsTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.queue_key = 'some_key'
        self.queue = MagicMock()
        self.models['Queue'].objects.get.return_value = self.queue

    def _claim_queue_items(self):
        return self.dao.claim_queue_items(queue_key=self.queue_key)

    def test_gets_queue(self):
        self._claim_queue_items()
        self.assertEqual(self.models['Queue'].objects.get.call_args,
                         call(uuid=self.queue_key))

    def test_claims_items(self):
        self._claim_queue_items()
        self.assertEqual(self.queue_utils.claim_queue_items.call_args,
                         call(queue=self.queue, models=self.models))

    def test_returns_serialized_items(self):
        result = self._claim_queue_items()
        self.assertEqual(
            self.queue_utils.serialize_queue_items.call_args,
            call(queue=self.queue,
                 items=self.queue_utils.claim_queue_items.return_value,
                 serializers=self.serializers)
        )
        self.assertEqual(result,
                         self.queue_utils.serialize_queue_items.return_value)

if __name__ == '__main__':
    unittest.main()


